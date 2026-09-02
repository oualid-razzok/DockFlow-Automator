"""Docking result analysis: interactions, RMSD, clustering, efficiency.

The contact model works directly on PDBQT atom types (the receptor from
:class:`~dockflow_core.preparator.ReceptorPreparator` keeps polar hydrogens
with ``HD`` types, so hydrogen bonds can be detected geometrically):

* hydrogen bond : ``HD`` on one side within 3.5 A of an acceptor
  (``NA/OA/SA/OS/N/O/S``) on the other side.
* hydrophobic   : ``C/A`` on both sides within 4.5 A.
* ionic         : ASP/GLU carboxylate O or ARG/LYS protonated N within 4.5 A
  of an oppositely charged ligand atom.
* metal         : receptor metal (Zn/Mg/...) within 3.0 A of a ligand
  N/O/S.

Heavy math (pairwise distances, Kabsch RMSD) is transparently offloaded to
the C++ accelerator module when it is installed.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .models import Contact, DockingResult
from .pdbio import Atom, parse_pdbqt, split_pdbqt_models
from .utils import get_logger

logger = get_logger("analyzer")

__all__ = [
    "analyze_interactions",
    "classify_pair",
    "contact_summary",
    "kabsch_rmsd",
    "direct_rmsd",
    "cluster_poses",
    "ligand_efficiency",
    "analyze_docking_result",
    "write_contacts_csv",
    "ResidueContactRow",
]

ACCEPTOR_TYPES = {"NA", "OA", "SA", "OS", "N", "O", "S"}
DONOR_TYPES = {"HD"}
HYDROPHOBIC_TYPES = {"C", "A"}
METALS = {"Mg", "Ca", "Mn", "Fe", "Zn", "NI", "CU", "CO"}

_NEG_RESNAME_O = {"ASP", "GLU"}   # carboxylate oxygens OD1/OD2/OE1/OE2
_POS_RESNAME_N = {"ARG", "LYS"}   # NH1/NH2/NZ
_LIG_CATION_TYPES = {"N", "NA", "NS"}
_LIG_ANION_TYPES = {"OA", "OS", "O", "SA", "S"}

HBOND_CUTOFF = 3.5
HYDROPHOBIC_CUTOFF = 4.5
IONIC_CUTOFF = 4.5
METAL_CUTOFF = 3.0
CONTACT_CUTOFF = 5.0


# ---------------------------------------------------------------------------
# Accelerated kernels
# ---------------------------------------------------------------------------
def _pairwise_min_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """For every row of ``a``, the distance to the closest row of ``b``."""
    try:
        import dockflow_bindings as _dfb  # type: ignore

        return _dfb.pairwise_min_dist(a, b)
    except Exception:  # noqa: BLE001 - pure numpy fallback
        diff = a[:, None, :] - b[None, :, :]
        return np.sqrt((diff * diff).sum(axis=2)).min(axis=1)


def _contacts_within(a: np.ndarray, b: np.ndarray, cutoff: float):
    """Yield (i, j, distance) for all pairs closer than ``cutoff``."""
    try:
        import dockflow_bindings as _dfb  # type: ignore

        return _dfb.min_contacts(a, b, cutoff)
    except Exception:  # noqa: BLE001
        diff = a[:, None, :] - b[None, :, :]
        dist = np.sqrt((diff * diff).sum(axis=2))
        ii, jj = np.where(dist <= cutoff)
        return [(int(i), int(j), float(dist[i, j]))
                for i, j in zip(ii, jj, strict=True)]


# ---------------------------------------------------------------------------
# RMSD
# ---------------------------------------------------------------------------
def direct_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """Plain RMSD without superposition."""
    a = np.asarray(coords_a, dtype=float).reshape(-1, 3)
    b = np.asarray(coords_b, dtype=float).reshape(-1, 3)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch for RMSD: {a.shape} vs {b.shape}")
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def kabsch_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float:
    """RMSD after optimal superposition (Kabsch algorithm)."""
    a = np.asarray(coords_a, dtype=float).reshape(-1, 3)
    b = np.asarray(coords_b, dtype=float).reshape(-1, 3)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch for RMSD: {a.shape} vs {b.shape}")
    try:
        import dockflow_bindings as _dfb  # type: ignore

        return float(_dfb.kabsch_rmsd(a, b))
    except Exception:  # noqa: BLE001
        pass
    pa = a - a.mean(axis=0)
    pb = b - b.mean(axis=0)
    covariance = pa.T @ pb
    u, _, vt = np.linalg.svd(covariance)
    d = np.sign(np.linalg.det(u @ vt))
    rotation = u @ np.diag([1.0, 1.0, d]) @ vt
    aligned = pa @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - pb) ** 2, axis=1))))


# ---------------------------------------------------------------------------
# Contact detection
# ---------------------------------------------------------------------------
def classify_pair(
    ligand_type: str,
    receptor_type: str,
    distance: float,
    receptor_resname: str,
    receptor_atom_name: str,
) -> str | None:
    """Classify one ligand-receptor atom pair; ``None`` when uninteresting."""
    lig = ligand_type.strip().upper()
    rec = receptor_type.strip().upper()
    resname = receptor_resname.strip().upper()
    rec_name = receptor_atom_name.strip().upper()
    # Hydrogen bonds.
    if lig in DONOR_TYPES and rec in ACCEPTOR_TYPES and distance <= HBOND_CUTOFF:
        return "hbond"
    if rec in DONOR_TYPES and lig in ACCEPTOR_TYPES and distance <= HBOND_CUTOFF:
        return "hbond"
    # Metal coordination (receptor metal to ligand heteroatom).
    if rec in {m.upper() for m in METALS} and lig in _LIG_ANION_TYPES | _LIG_CATION_TYPES \
            and distance <= METAL_CUTOFF:
        return "metal"
    # Ionic contacts via charged residues.
    if resname in _NEG_RESNAME_O and rec_name.startswith(("OD", "OE")) \
            and lig in _LIG_CATION_TYPES and distance <= IONIC_CUTOFF:
        return "ionic"
    if resname in _POS_RESNAME_N and rec_name.startswith(("NH", "NZ")) \
            and lig in _LIG_ANION_TYPES and distance <= IONIC_CUTOFF:
        return "ionic"
    # Hydrophobic contacts.
    if lig in HYDROPHOBIC_TYPES and rec in HYDROPHOBIC_TYPES \
            and distance <= HYDROPHOBIC_CUTOFF:
        return "hydrophobic"
    return None


def analyze_interactions(
    ligand_atoms: Sequence[Atom],
    receptor_atoms: Sequence[Atom],
    cutoff: float = CONTACT_CUTOFF,
) -> list[Contact]:
    """Detect all ligand-receptor contacts within ``cutoff``."""
    if not ligand_atoms or not receptor_atoms:
        return []
    lig_coords = np.array([[a.x, a.y, a.z] for a in ligand_atoms], dtype=float)
    rec_coords = np.array([[a.x, a.y, a.z] for a in receptor_atoms], dtype=float)
    pairs = _contacts_within(lig_coords, rec_coords, cutoff)
    contacts: list[Contact] = []
    for i, j, distance in pairs:
        lig_atom = ligand_atoms[i]
        rec_atom = receptor_atoms[j]
        kind = classify_pair(
            lig_atom.atom_type, rec_atom.atom_type, distance,
            rec_atom.resname, rec_atom.name,
        )
        if kind is None:
            continue
        contacts.append(
            Contact(
                ligand_atom_index=i,
                ligand_atom_name=lig_atom.name.strip(),
                ligand_atom_type=lig_atom.atom_type,
                receptor_atom_name=rec_atom.name.strip(),
                receptor_resname=rec_atom.resname.strip(),
                receptor_chain=rec_atom.chain,
                receptor_resseq=rec_atom.resseq,
                receptor_atom_type=rec_atom.atom_type,
                distance=round(distance, 2),
                kind=kind,
            )
        )
    contacts.sort(key=lambda c: c.distance)
    return contacts



@dataclass
class ResidueContactRow:
    """Aggregated contacts of one receptor residue."""

    chain: str
    resname: str
    resseq: int
    hbonds: int = 0
    hydrophobic: int = 0
    ionic: int = 0
    metal: int = 0
    closest: float = 999.0

    @property
    def total(self) -> int:
        return self.hbonds + self.hydrophobic + self.ionic + self.metal


def contact_summary(contacts: Sequence[Contact]) -> list[ResidueContactRow]:
    """Aggregate contacts per receptor residue, most-contacted first."""
    rows: dict[tuple[str, str, int], ResidueContactRow] = {}
    for contact in contacts:
        key = (contact.receptor_chain, contact.receptor_resname, contact.receptor_resseq)
        row = rows.get(key)
        if row is None:
            row = ResidueContactRow(
                chain=contact.receptor_chain,
                resname=contact.receptor_resname,
                resseq=contact.receptor_resseq,
            )
            rows[key] = row
        if contact.kind == "hbond":
            row.hbonds += 1
        elif contact.kind == "hydrophobic":
            row.hydrophobic += 1
        elif contact.kind == "ionic":
            row.ionic += 1
        elif contact.kind == "metal":
            row.metal += 1
        row.closest = min(row.closest, contact.distance)
    return sorted(rows.values(), key=lambda r: (-r.total, r.closest))


# ---------------------------------------------------------------------------
# Pose clustering & efficiency
# ---------------------------------------------------------------------------
def cluster_poses(
    pose_coordinates: Sequence[np.ndarray],
    cutoff: float = 2.0,
) -> list[list[int]]:
    """Greedy RMSD clustering (best poses are expected first)."""
    coords = [np.asarray(c, dtype=float).reshape(-1, 3) for c in pose_coordinates]
    if not coords:
        return []
    if len({c.shape for c in coords}) > 1:
        raise ValueError("all poses must have the same number of atoms to cluster")
    clusters: list[list[int]] = []
    representatives: list[np.ndarray] = []
    for index, c in enumerate(coords):
        assigned = False
        for cluster_index, rep in enumerate(representatives):
            rmsd = kabsch_rmsd(c, rep)
            if rmsd <= cutoff:
                clusters[cluster_index].append(index)
                assigned = True
                break
        if not assigned:
            clusters.append([index])
            representatives.append(c.copy())
    return clusters


def ligand_efficiency(affinity: float, num_heavy_atoms: int) -> float | None:
    """Affinity per heavy atom (kcal/mol per atom)."""
    if num_heavy_atoms and num_heavy_atoms > 0:
        return affinity / num_heavy_atoms
    return None


# ---------------------------------------------------------------------------
# Higher level analysis
# ---------------------------------------------------------------------------
@dataclass
class PoseAnalysis:
    """Analysis of a single docked pose."""

    pose_index: int = 0
    affinity: float = 0.0
    num_contacts: int = 0
    num_hbonds: int = 0
    num_hydrophobic: int = 0
    num_ionic: int = 0
    num_metal: int = 0
    contacts: list[Contact] = field(default_factory=list)
    residue_rows: list[ResidueContactRow] = field(default_factory=list)
    ligand_efficiency: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "pose_index": self.pose_index,
            "affinity": self.affinity,
            "num_contacts": self.num_contacts,
            "num_hbonds": self.num_hbonds,
            "num_hydrophobic": self.num_hydrophobic,
            "num_ionic": self.num_ionic,
            "num_metal": self.num_metal,
            "ligand_efficiency": self.ligand_efficiency,
            "contacts": [c.to_dict() for c in self.contacts],
            "residues": [
                {
                    "chain": r.chain, "resname": r.resname, "resseq": r.resseq,
                    "hbonds": r.hbonds, "hydrophobic": r.hydrophobic,
                    "ionic": r.ionic, "metal": r.metal, "closest": r.closest,
                }
                for r in self.residue_rows
            ],
        }


def analyze_docking_result(
    result: DockingResult,
    receptor_pdbqt: str | Path,
    top_poses: int = 3,
    cutoff: float = CONTACT_CUTOFF,
    ligand_heavy_atoms: int | None = None,
) -> list[PoseAnalysis]:
    """Analyse the top poses of one docking result against a receptor PDBQT."""
    if result.out_path is None or not Path(result.out_path).is_file():
        return []
    receptor_atoms = parse_pdbqt(receptor_pdbqt).atoms
    pose_files = split_pdbqt_models(result.out_path, Path(result.out_path).parent,
                                    result.ligand_name or "pose")
    analyses: list[PoseAnalysis] = []
    for pose_file in pose_files[:top_poses]:
        pose_data = parse_pdbqt(pose_file)
        pose_atoms = pose_data.atoms
        vina_result = pose_data.models[0].vina_result if pose_data.models else None
        affinity = vina_result.affinity if vina_result else 0.0
        contacts = analyze_interactions(pose_atoms, receptor_atoms, cutoff=cutoff)
        heavy = ligand_heavy_atoms
        if heavy is None:
            heavy = sum(1 for a in pose_atoms if a.element.upper() != "H")
        analysis = PoseAnalysis(
            pose_index=len(analyses) + 1,
            affinity=affinity,
            num_contacts=len(contacts),
            num_hbonds=sum(1 for c in contacts if c.kind == "hbond"),
            num_hydrophobic=sum(1 for c in contacts if c.kind == "hydrophobic"),
            num_ionic=sum(1 for c in contacts if c.kind == "ionic"),
            num_metal=sum(1 for c in contacts if c.kind == "metal"),
            contacts=contacts,
            residue_rows=contact_summary(contacts),
            ligand_efficiency=ligand_efficiency(affinity, heavy),
        )
        analyses.append(analysis)
    return analyses


def write_contacts_csv(contacts: Sequence[Contact], path: str | Path) -> Path:
    """Write one row per contact for spreadsheet inspection."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(target, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["kind", "distance", "ligand_atom", "ligand_type",
             "receptor_chain", "receptor_resname", "receptor_resseq",
             "receptor_atom", "receptor_type"]
        )
        for contact in contacts:
            writer.writerow(
                [contact.kind, f"{contact.distance:.2f}", contact.ligand_atom_name,
                 contact.ligand_atom_type, contact.receptor_chain,
                 contact.receptor_resname, contact.receptor_resseq,
                 contact.receptor_atom_name, contact.receptor_atom_type]
            )
    return target


def write_analysis_json(analyses: Sequence[PoseAnalysis], path: str | Path) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = {"poses": [a.to_dict() for a in analyses]}
    target.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return target
