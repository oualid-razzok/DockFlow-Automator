"""Docking result analysis: contacts, RMSD, clustering, efficiency.

The contact model is **geometric / contact-based**, not a rigorous
hydrogen-bond or ionic-bond analysis: it works directly on PDBQT atom
types (the receptor from :class:`~dockflow_core.preparator.ReceptorPreparator`
keeps polar hydrogens with ``HD`` types) and classifies a pair purely by
distance.  No angle criterion is applied.

* geometric H-bond contact (donor-acceptor): ``HD`` on one side within
  3.5 A of an acceptor (``NA/OA/SA/OS/N/O/S``) on the other side.
* hydrophobic contact : ``C/A`` on both sides within 4.5 A.
* ionic contact     : ASP/GLU carboxylate O or ARG/LYS protonated N within
  4.5 A of an oppositely charged ligand atom.
* metal contact     : receptor metal (Zn/Mg/...) within 3.0 A of a ligand
  N/O/S.

See ``docs/interaction_criteria.md`` for the full criteria table and the
limitations of a distance-only model.

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
from .pdbio import Atom, parse_pdb, parse_pdbqt, split_pdbqt_models
from .utils import get_logger

logger = get_logger("analyzer")

__all__ = [
    "chemical_element",
    "analyze_interactions",
    "classify_pair",
    "contact_summary",
    "CONTACT_KIND_LABELS",
    "kabsch_rmsd",
    "direct_rmsd",
    "symmetry_tolerant_rmsd",
    "extract_reference_atoms",
    "crystal_rmsd",
    "cluster_poses",
    "pose_cluster_summary",
    "docking_score_efficiency",
    "ligand_efficiency",
    "analyze_docking_result",
    "write_contacts_csv",
    "ResidueContactRow",
]

ACCEPTOR_TYPES = {"NA", "OA", "SA", "OS", "N", "O", "S"}
DONOR_TYPES = {"HD"}
HYDROPHOBIC_TYPES = {"C", "A"}
METALS = {"Mg", "Ca", "Mn", "Fe", "Zn", "NI", "CU", "CO"}

# Display labels: every user-facing surface (report.md, CLI prints) must
# qualify contacts as geometric, never as verified bonds (audit item 13).
CONTACT_KIND_LABELS = {
    "hbond": "geometric H-bond contact (D-A <= 3.5 A)",
    "hydrophobic": "hydrophobic contact (<= 4.5 A)",
    "ionic": "ionic contact (<= 4.5 A)",
    "metal": "metal contact (<= 3.0 A)",
    "close": "close contact",
}

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
# Symmetry-tolerant RMSD (redocking validation)
# ---------------------------------------------------------------------------
# AutoDock PDBQT atom types -> chemical elements.  Meeko emits extra
# *virtual* atoms for guanidinium-like groups (types G0-G3, atom name "G")
# which have no counterpart in the crystal pose and must be excluded from
# RMSD computations; aromatic carbons are typed "A", heteroatom acceptors
# "NA/OA/SA/OS", guanidinium carbons "CG0-CG3".
_AD4_TYPE_TO_ELEMENT = {
    "A": "C", "C": "C", "CG0": "C", "CG1": "C", "CG2": "C", "CG3": "C",
    "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S", "SE": "Se",
    "P": "P", "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "H": "H", "HD": "H", "HS": "H",
    "SI": "Si", "B": "B",
    "MG": "Mg", "CA": "Ca", "MN": "Mn", "FE": "Fe", "ZN": "Zn",
    "NI": "Ni", "CU": "Cu", "CO": "Co",
}
_REAL_ELEMENTS = {
    "H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE", "NA", "MG",
    "AL", "SI", "P", "S", "CL", "AR", "K", "CA", "SC", "TI", "V", "CR",
    "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR",
    "KR", "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD",
    "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE", "CS", "BA", "LA",
    "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER",
    "TM", "YB", "LU", "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU",
    "HG", "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH",
    "PA", "U", "NP", "PU",
}


def chemical_element(atom: Atom) -> str | None:
    """Chemical element of a PDB/PDBQT atom; ``None`` for virtual atoms.

    PDBQT ligands carry AutoDock types in the element column (``A`` for
    aromatic carbon, ``G0`` for meeko's guanidinium virtual atoms, ...),
    so the type mapping is authoritative when present.
    """
    for candidate in (atom.atom_type, atom.element):
        token = (candidate or "").strip().upper()
        if not token:
            continue
        if token in _AD4_TYPE_TO_ELEMENT:
            return _AD4_TYPE_TO_ELEMENT[token]
        if token in _REAL_ELEMENTS:
            return token
    return None  # virtual (G0-G3, W, XX) or unrecognised


def _greedy_same_element_assignment(
    a: np.ndarray, b: np.ndarray
) -> list[tuple[int, int]]:
    """Greedy mutually-closest pairing between two same-element atom groups."""
    if len(a) != len(b):
        raise ValueError("greedy assignment needs equal group sizes")
    diff = a[:, None, :] - b[None, :, :]
    dist = np.sqrt((diff * diff).sum(axis=2))
    pairs: list[tuple[int, int]] = []
    used_a: set[int] = set()
    used_b: set[int] = set()
    order = sorted(
        (dist[i, j], i, j) for i in range(len(a)) for j in range(len(b))
    )
    for _, i, j in order:
        if i in used_a or j in used_b:
            continue
        used_a.add(i)
        used_b.add(j)
        pairs.append((int(i), int(j)))
    return pairs


def symmetry_tolerant_rmsd(
    pose_atoms: Sequence[Atom],
    reference_atoms: Sequence[Atom],
    max_refinements: int = 4,
    return_correspondence: bool = False,
) -> float | list[tuple[int, int]]:
    """Symmetry-aware heavy-atom Kabsch RMSD between a pose and a reference.

    Standard Kabsch RMSD assumes a fixed 1:1 atom correspondence, which
    over-penalises symmetric ligands (a 180-degree flip of a carboxylate or
    a ring nitrogen permutation is the same chemistry).  This function:

    1. matches atoms by element (groups must have equal counts),
    2. superposes the current correspondence with Kabsch,
    3. re-matches same-element atoms by mutual nearest distance under the
       current superposition (greedy), and
    4. iterates (3)-(4) until the correspondence is stable.

    The result is the RMSD of the best correspondence found, i.e. the
    smallest RMSD over chemically equivalent atom mappings - the same
    convention used by redocking benchmarks (``compute_rmsd``-style tools).
    With ``return_correspondence=True`` the final (pose_idx, ref_idx)
    correspondence over the heavy-atom lists is returned instead (used to
    seed the graph-automorphism RMSD, peer item 17).
    """
    pose_heavy = [a for a in pose_atoms
                  if chemical_element(a) not in (None, "H")]
    ref_heavy = [a for a in reference_atoms
                 if chemical_element(a) not in (None, "H")]
    pose_elements = sorted(chemical_element(a) for a in pose_heavy)
    ref_elements = sorted(chemical_element(a) for a in ref_heavy)
    if pose_elements != ref_elements:
        raise ValueError(
            "cannot compute crystal RMSD: pose and reference ligands have "
            f"different atom sets ({pose_elements} vs {ref_elements})"
        )
    pose_xyz = np.array([[a.x, a.y, a.z] for a in pose_heavy], dtype=float)
    ref_xyz = np.array([[a.x, a.y, a.z] for a in ref_heavy], dtype=float)
    pose_elem = np.array([chemical_element(a) for a in pose_heavy])
    ref_elem = np.array([chemical_element(a) for a in ref_heavy])

    # Initial correspondence: element by element, input order.
    correspondence = [
        (i, j)
        for element in dict.fromkeys(pose_elem)
        for i, j in zip(
            (k for k in range(len(pose_heavy)) if pose_elem[k] == element),
            (k for k in range(len(ref_heavy)) if ref_elem[k] == element),
            strict=True,
        )
    ]
    best_rmsd = np.inf
    for _ in range(max(1, max_refinements)):
        order_a = np.array([i for i, _ in correspondence], dtype=int)
        order_b = np.array([j for _, j in correspondence], dtype=int)
        # Superpose the paired subsets, apply to ALL pose atoms.
        pa = pose_xyz[order_a]
        pb = ref_xyz[order_b]
        centroid_a = pa.mean(axis=0)
        centroid_b = pb.mean(axis=0)
        covariance = (pa - centroid_a).T @ (pb - centroid_b)
        u, _, vt = np.linalg.svd(covariance)
        d = np.sign(np.linalg.det(u @ vt))
        rotation = u @ np.diag([1.0, 1.0, d]) @ vt
        transformed = (pose_xyz - centroid_a) @ rotation + centroid_b
        # Re-match by element under the current superposition.
        new_correspondence: list[tuple[int, int]] = []
        for element in dict.fromkeys(pose_elem):
            pose_idx = np.where(pose_elem == element)[0]
            ref_idx = np.where(ref_elem == element)[0]
            assignment = _greedy_same_element_assignment(
                transformed[pose_idx], ref_xyz[ref_idx]
            )
            new_correspondence += [
                (int(pose_idx[i]), int(ref_idx[j])) for i, j in assignment
            ]
        deltas = transformed[[i for i, _ in new_correspondence]] - \
            ref_xyz[[j for _, j in new_correspondence]]
        rmsd = float(np.sqrt(np.mean(np.sum(deltas * deltas, axis=1))))
        best_rmsd = min(best_rmsd, rmsd)
        if new_correspondence == correspondence:
            break
        correspondence = new_correspondence
    if return_correspondence:
        return sorted(correspondence)
    return best_rmsd


def extract_reference_atoms(
    structure_path: str | Path, resname: str, chain: str | None = None
) -> list[Atom]:
    """Heavy atoms of ONE ligand residue instance extracted from a structure.

    The co-crystallized pose of the reference ligand (e.g. ``XK2`` inside
    ``1HVR.pdb``) is the ground truth for redocking validation.  Alternate
    conformations (altLoc) are resolved to the highest-occupancy variant so
    duplicated atoms do not corrupt the atom count.  When the ligand is
    present in several copies (homodimers, multi-copy entries), ONE
    instance is used (the largest) - the reference must describe a single
    binding site, not the union of two.
    """

    atoms = parse_pdb(structure_path)
    selected = [
        a
        for a in atoms
        if a.resname.strip().upper() == resname.strip().upper()
        and (chain is None or a.chain == chain)
        and chemical_element(a) != "H"
    ]
    if not selected:
        raise ValueError(
            f"reference ligand {resname!r} not found in {structure_path}"
        )
    # Resolve alternate conformations / duplicate atom names to the
    # highest-occupancy variant so duplicated atoms do not corrupt the
    # heavy-atom count used for the RMSD correspondence.
    best: dict[tuple[str, int, str, str], Atom] = {}
    for atom in selected:
        key = (atom.chain, atom.resseq, atom.icode, atom.name.strip())
        if key not in best or atom.occupancy > best[key].occupancy:
            best[key] = atom
    selected = list(best.values())
    # One residue instance only (largest) when several copies exist.
    instances: dict[tuple[str, int], list[Atom]] = {}
    for atom in selected:
        instances.setdefault((atom.chain.strip(), int(atom.resseq)), []).append(atom)
    if len(instances) > 1:
        selected = max(instances.values(), key=len)
    return selected


def crystal_rmsd(
    pose_atoms: Sequence[Atom], reference_atoms: Sequence[Atom]
) -> float | None:
    """Symmetry-aware heavy-atom RMSD of a pose vs the crystal pose.

    Uses RDKit ``GetBestRMS`` (graph-automorphism aware, peer item 17) when
    RDKit is importable and both atom sets can be converted into a
    proximity-bonded molecule; otherwise falls back to the element-greedy
    Kabsch heuristic.  See :func:`crystal_rmsd_with_method` for the variant
    that also reports which method was used.

    Returns ``None`` (with a logged debug note) when the pose and reference
    are not comparable, so callers can treat redocking validation as
    optional rather than fatal.
    """
    return crystal_rmsd_with_method(pose_atoms, reference_atoms)[0]


def _rdkit_mol_from_heavy_atoms(atoms: Sequence[Atom]):
    """Build an RDKit molecule from heavy ``Atom`` records.

    The atoms are written to an in-memory PDB block (real element symbols)
    and parsed with ``proximityBonding=True`` so connectivity is inferred
    from the 3D coordinates.  Returns ``None`` when RDKit is unavailable,
    the block cannot be parsed, or the result is bondless.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return None
    lines = []
    for serial, atom in enumerate(atoms, start=1):
        element = chemical_element(atom)
        if element is None:
            return None
        name = f" {element}  " if len(element) == 1 else f"{element}  "
        lines.append(
            f"ATOM  {serial:5d} {name} LIG A   1    "
            f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}  "
        )
    block = "\n".join(lines) + "\nEND\n"
    mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False,
                               proximityBonding=True)
    if mol is None or mol.GetNumAtoms() != len(atoms) or mol.GetNumBonds() == 0:
        return None
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, catchErrors=True)
        Chem.FastFindRings(mol)
    except Exception:  # noqa: BLE001
        return None
    return mol


def rdkit_symmetry_rmsd(
    pose_atoms: Sequence[Atom], reference_atoms: Sequence[Atom]
) -> float | None:
    """Graph-automorphism-aware RMSD via RDKit ``GetBestRMS`` (item 17).

    The reference's proximity-bonded molecular graph is transplanted onto
    the pose: pose coordinates are assigned to reference atoms through the
    *geometry-matched* element correspondence (the greedy Kabsch matching
    from :func:`symmetry_tolerant_rmsd`, so the initial assignment is
    chemically sensible rather than input-ordered), and ``GetBestRMS``
    then minimises the RMSD over all automorphisms of the shared graph
    (ring flips, carboxylate turns, terminal permutations).  This is the
    standard symmetry-correct ligand RMSD.  Returns ``None`` when the
    RDKit path is not usable.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolAlign
    except ImportError:
        return None
    pose_heavy = [a for a in pose_atoms
                  if chemical_element(a) not in (None, "H")]
    ref_heavy = [a for a in reference_atoms
                 if chemical_element(a) not in (None, "H")]
    ref_elements = sorted(chemical_element(a) for a in ref_heavy)
    pose_elements = sorted(chemical_element(a) for a in pose_heavy)
    if pose_elements != ref_elements or not ref_heavy:
        return None
    ref_mol = _rdkit_mol_from_heavy_atoms(ref_heavy)
    if ref_mol is None:
        return None
    # Seed assignment with the geometry-matched element correspondence
    # (the greedy Kabsch matching): the pairing decides which pose
    # coordinate lands on which reference atom.  GetBestRMS re-permutes
    # via graph automorphisms, but a sensible seed keeps even asymmetric
    # ligands correct (input-order seeding scrambles them).
    try:
        correspondence = symmetry_tolerant_rmsd(
            pose_atoms, reference_atoms, return_correspondence=True)
    except ValueError:
        return None
    conf = Chem.Conformer(ref_mol.GetNumAtoms())
    for pose_index, ref_index in correspondence:
        pose_atom = pose_heavy[pose_index]
        conf.SetAtomPosition(ref_index,
                             (float(pose_atom.x), float(pose_atom.y),
                              float(pose_atom.z)))
    pose_mol = Chem.Mol(ref_mol)
    pose_mol.RemoveAllConformers()
    pose_mol.AddConformer(conf, assignId=True)
    return float(rdMolAlign.GetBestRMS(pose_mol, ref_mol))


def crystal_rmsd_with_method(
    pose_atoms: Sequence[Atom], reference_atoms: Sequence[Atom]
) -> tuple[float | None, str | None]:
    """Symmetry-aware RMSD plus the method actually used (peer item 17).

    Returns ``(rmsd, method)`` where method is
    ``"rdkit-getbestrms-graph-automorphism"`` (preferred: exact minimisation
    over molecular graph symmetries) or ``"element-greedy-kabsch"`` (the
    dependency-free fallback).  ``(None, None)`` when not computable.
    """
    try:
        rmsd = rdkit_symmetry_rmsd(pose_atoms, reference_atoms)
        if rmsd is not None:
            return rmsd, "rdkit-getbestrms-graph-automorphism"
    except Exception:  # noqa: BLE001
        logger.debug("rdkit GetBestRMS path failed", exc_info=True)
    try:
        return (symmetry_tolerant_rmsd(pose_atoms, reference_atoms),
                "element-greedy-kabsch")
    except ValueError as exc:
        logger.debug("crystal RMSD not computable: %s", exc)
        return None, None


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
    """Classify one ligand-receptor atom pair; ``None`` when uninteresting.

    The returned kind is a *geometric contact* classification (distance +
    atom type only, no angle criterion) - see the module docstring.
    """
    lig = ligand_type.strip().upper()
    rec = receptor_type.strip().upper()
    resname = receptor_resname.strip().upper()
    rec_name = receptor_atom_name.strip().upper()
    # Geometric hydrogen-bond contacts (donor-acceptor distance criterion).
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
    geom_hbond: int = 0   # geometric H-bond contacts (distance criterion only)
    hydrophobic: int = 0
    ionic: int = 0
    metal: int = 0
    closest: float = 999.0

    @property
    def total(self) -> int:
        return self.geom_hbond + self.hydrophobic + self.ionic + self.metal


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
            row.geom_hbond += 1
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


def pose_cluster_summary(
    pose_coordinates: Sequence[np.ndarray],
    affinities: Sequence[float],
    cutoff: float = 2.0,
) -> list[dict[str, Any]]:
    """Cluster summary rows for a pose ensemble (audit item 20).

    Shows whether N poses are one converged cluster or several distinct
    binding modes: cluster id, member poses (1-based), mean and best
    affinity, intra-cluster spread (max pairwise Kabsch RMSD) and the
    representative pose (best affinity, since poses arrive ranked).
    """
    coords = [np.asarray(c, dtype=float).reshape(-1, 3) for c in pose_coordinates]
    if not coords:
        return []
    clusters = cluster_poses(coords, cutoff=cutoff)
    rows: list[dict[str, Any]] = []
    for cluster_id, members in enumerate(clusters, start=1):
        # guard against affinities lists shorter than the pose files
        # (mocked engines, truncated results)
        member_affinities = [float(affinities[i]) for i in members
                             if i < len(affinities)] or [0.0]
        spread = 0.0
        for position, i in enumerate(members):
            for j in members[position + 1:]:
                spread = max(spread, kabsch_rmsd(coords[i], coords[j]))
        rows.append(
            {
                "cluster": cluster_id,
                "poses": [i + 1 for i in members],
                "representative_pose": members[0] + 1,
                "mean_affinity": round(
                    sum(member_affinities) / len(member_affinities), 3
                ),
                "best_affinity": round(min(member_affinities), 3),
                "intra_cluster_rmsd_spread": round(float(spread), 3),
            }
        )
    return rows


def docking_score_efficiency(affinity: float, num_heavy_atoms: int) -> float | None:
    """Docking score per heavy atom: ``affinity / num_heavy_atoms``.

    A *proxy* derived from the Vina docking score, **not** the
    experimentally derived ligand efficiency (Delta G / heavy-atom count
    from binding assays).  Useful for comparing ligands of very different
    sizes under the same scoring function; it carries no thermodynamic
    meaning.  Reports ``None`` for empty molecules.
    """
    if num_heavy_atoms and num_heavy_atoms > 0:
        return affinity / num_heavy_atoms
    return None


# Deprecated alias kept for one release (renamed in 0.2.0, audit item 14).
ligand_efficiency = docking_score_efficiency


# ---------------------------------------------------------------------------
# Higher level analysis
# ---------------------------------------------------------------------------
@dataclass
class PoseAnalysis:
    """Analysis of a single docked pose."""

    pose_index: int = 0
    affinity: float = 0.0
    num_contacts: int = 0
    num_geom_hbond: int = 0    # geometric H-bond contacts (distance only)
    num_hydrophobic: int = 0
    num_ionic: int = 0
    num_metal: int = 0
    contacts: list[Contact] = field(default_factory=list)
    residue_rows: list[ResidueContactRow] = field(default_factory=list)
    # Proxy: docking score / heavy-atom count, NOT experimental ligand
    # efficiency (renamed in 0.2.0; audit item 14).
    docking_score_efficiency: float | None = None
    # Symmetry-aware heavy-atom RMSD to the co-crystallized reference
    # pose (None when no reference was available; audit item 18).
    crystal_rmsd: float | None = None
    # Which RMSD method produced crystal_rmsd (peer item 17):
    # "rdkit-getbestrms-graph-automorphism" or "element-greedy-kabsch".
    crystal_rmsd_method: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "pose_index": self.pose_index,
            "affinity": self.affinity,
            "num_contacts": self.num_contacts,
            "num_geom_hbond_contacts": self.num_geom_hbond,
            "num_hydrophobic_contacts": self.num_hydrophobic,
            "num_ionic_contacts": self.num_ionic,
            "num_metal_contacts": self.num_metal,
            "docking_score_efficiency": self.docking_score_efficiency,
            "crystal_rmsd": self.crystal_rmsd,
            "crystal_rmsd_method": self.crystal_rmsd_method,
            "contacts": [c.to_dict() for c in self.contacts],
            "residues": [
                {
                    "chain": r.chain, "resname": r.resname, "resseq": r.resseq,
                    "geom_hbond_contacts": r.geom_hbond,
                    "hydrophobic_contacts": r.hydrophobic,
                    "ionic_contacts": r.ionic, "metal_contacts": r.metal,
                    "closest": r.closest,
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
    reference_atoms: Sequence[Atom] | None = None,
) -> list[PoseAnalysis]:
    """Analyse the top poses of one docking result against a receptor PDBQT.

    When ``reference_atoms`` (the co-crystallized ligand, see
    :func:`extract_reference_atoms`) is given, a symmetry-tolerant heavy-atom
    ``crystal_rmsd`` is computed for **every** pose and written back onto
    ``result.poses[i].crystal_rmsd`` so ``summary.csv`` can report it.
    """
    if result.out_path is None or not Path(result.out_path).is_file():
        return []
    receptor_atoms = parse_pdbqt(receptor_pdbqt).atoms
    pose_files = split_pdbqt_models(result.out_path, Path(result.out_path).parent,
                                    result.ligand_name or "pose")
    analyses: list[PoseAnalysis] = []
    rmsd_method: str | None = None
    for index, pose_file in enumerate(pose_files):
        pose_data = parse_pdbqt(pose_file)
        pose_atoms = pose_data.atoms
        # Redocking validation for every pose, not just the top ones, so
        # summary.csv can carry a crystal_rmsd column.  Symmetry-aware
        # RMSD (peer item 17): RDKit GetBestRMS when possible, else the
        # element-greedy Kabsch fallback; the method is recorded.
        if reference_atoms is not None and index < len(result.poses):
            value, method = crystal_rmsd_with_method(pose_atoms, reference_atoms)
            result.poses[index].crystal_rmsd = value
            if method:
                rmsd_method = method
        if index >= top_poses:
            continue
        vina_result = pose_data.models[0].vina_result if pose_data.models else None
        affinity = vina_result.affinity if vina_result else 0.0
        contacts = analyze_interactions(pose_atoms, receptor_atoms, cutoff=cutoff)
        heavy = ligand_heavy_atoms
        if heavy is None:
            heavy = sum(1 for a in pose_atoms
                        if chemical_element(a) not in (None, "H"))
        analysis = PoseAnalysis(
            pose_index=len(analyses) + 1,
            affinity=affinity,
            num_contacts=len(contacts),
            num_geom_hbond=sum(1 for c in contacts if c.kind == "hbond"),
            num_hydrophobic=sum(1 for c in contacts if c.kind == "hydrophobic"),
            num_ionic=sum(1 for c in contacts if c.kind == "ionic"),
            num_metal=sum(1 for c in contacts if c.kind == "metal"),
            contacts=contacts,
            residue_rows=contact_summary(contacts),
            docking_score_efficiency=docking_score_efficiency(affinity, heavy),
            crystal_rmsd=result.poses[index].crystal_rmsd
            if index < len(result.poses) else None,
            crystal_rmsd_method=rmsd_method,
        )
        analyses.append(analysis)
    return analyses


def pose_coordinates(result: DockingResult) -> list[np.ndarray]:
    """Heavy-atom coordinates of every pose in a docking result."""
    if result.out_path is None or not Path(result.out_path).is_file():
        return []
    pose_files = split_pdbqt_models(result.out_path, Path(result.out_path).parent,
                                    result.ligand_name or "pose")
    coords: list[np.ndarray] = []
    for pose_file in pose_files:
        atoms = [a for a in parse_pdbqt(pose_file).atoms
                 if chemical_element(a) not in (None, "H")]
        if atoms:
            coords.append(np.array([[a.x, a.y, a.z] for a in atoms], dtype=float))
    return coords


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


# ---------------------------------------------------------------------------
# Interaction fingerprints (peer item 18)
# ---------------------------------------------------------------------------
def pose_ifp_keys(contacts: Sequence[Contact]) -> set[tuple[str, str, int, str]]:
    """The IFP bit keys of one pose: ``(chain, resname, resseq, kind)``.

    One bit per receptor residue x geometric contact type - the standard
    PLIF (protein-ligand interaction fingerprint) layout.  The keys are
    contact-derived, so an IFP says "which residues the pose touches and
    how", not "how strongly" (contact counts are not energies).
    """
    return {
        (contact.receptor_chain, contact.receptor_resname,
         int(contact.receptor_resseq), contact.kind)
        for contact in contacts
    }


def ifp_vocabulary(contact_sets: Sequence[set[tuple[str, str, int, str]]]
                   ) -> list[tuple[str, str, int, str]]:
    """Sorted union vocabulary over several poses (stable bit order)."""
    union: set[tuple[str, str, int, str]] = set()
    for keys in contact_sets:
        union |= keys
    return sorted(union)


def ifp_bits(keys: set[tuple[str, str, int, str]],
             vocabulary: Sequence[tuple[str, str, int, str]]) -> list[int]:
    """Bit list for one pose against a vocabulary (1 = contact present)."""
    wanted = set(keys)
    return [1 if bit in wanted else 0 for bit in vocabulary]


def tanimoto_ifp(bits_a: Sequence[int], bits_b: Sequence[int]) -> float:
    """Tanimoto coefficient of two IFP bit vectors (dependency-free)."""
    if len(bits_a) != len(bits_b):
        raise ValueError("IFP vectors must share one vocabulary")
    ones_a = sum(bits_a)
    ones_b = sum(bits_b)
    both = sum(1 for a, b in zip(bits_a, bits_b, strict=True) if a and b)
    if ones_a == 0 and ones_b == 0:
        return 1.0  # two contact-free poses are identical (both empty)
    return both / (ones_a + ones_b - both)


def ifp_similarity_matrix(bit_vectors: Sequence[Sequence[int]]
                           ) -> list[list[float]]:
    """Pairwise Tanimoto matrix over IFP bit vectors."""
    size = len(bit_vectors)
    matrix = [[1.0] * size for _ in range(size)]
    for i in range(size):
        for j in range(i + 1, size):
            similarity = tanimoto_ifp(bit_vectors[i], bit_vectors[j])
            matrix[i][j] = similarity
            matrix[j][i] = similarity
    return matrix


def cluster_ifp(bit_vectors: Sequence[Sequence[int]],
                threshold: float = 0.7) -> list[list[int]]:
    """Single-linkage clustering of poses by IFP Tanimoto similarity.

    Heuristic clustering by *interaction* similarity (not geometry): two
    poses belong to one cluster when their fingerprints overlap by at
    least ``threshold``.  Returns a list of clusters as lists of pose
    indices (0-based, input order).
    """
    size = len(bit_vectors)
    if size == 0:
        return []
    matrix = ifp_similarity_matrix(bit_vectors)
    assigned: list[int | None] = [None] * size
    clusters: list[list[int]] = []
    for i in range(size):
        if assigned[i] is not None:
            continue
        members = [i]
        assigned[i] = len(clusters)
        frontier = [i]
        while frontier:
            current = frontier.pop()
            for j in range(size):
                if assigned[j] is None and matrix[current][j] >= threshold:
                    assigned[j] = len(clusters)
                    members.append(j)
                    frontier.append(j)
        clusters.append(sorted(members))
    return clusters


def write_ifp_outputs(ligand_name: str,
                      analyses: Sequence[PoseAnalysis],
                      analysis_dir: str | Path) -> list[Path]:
    """Write the per-ligand IFP artifacts (peer item 18b).

    * ``<ligand>_ifp.csv``  - one row per pose, one column per bit
      (vocabulary = union over that ligand's analysed poses),
    * ``<ligand>_pose<N>_ifp.dat`` - RDKit ``ExplicitBitVect.to_binary``
      per pose over the same vocabulary (skipped without RDKit).

    Returns the list of files written (CSV always; .dat when RDKit is
    importable).
    """
    analysis_dir = Path(analysis_dir)
    analysis_dir.mkdir(parents=True, exist_ok=True)
    contact_sets = [pose_ifp_keys(analysis.contacts) for analysis in analyses]
    vocabulary = ifp_vocabulary(contact_sets)
    header = ["pose"] + [
        f"{chain or '-'}:{resname}{resseq}:{kind}"
        for chain, resname, resseq, kind in vocabulary
    ]
    csv_path = analysis_dir / f"{ligand_name}_ifp.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for analysis, keys in zip(analyses, contact_sets, strict=True):
            writer.writerow([analysis.pose_index]
                            + ifp_bits(keys, vocabulary))
    written = [csv_path]
    try:
        from rdkit.DataStructs import ExplicitBitVect

        for analysis, keys in zip(analyses, contact_sets, strict=True):
            bits = ifp_bits(keys, vocabulary)
            vector = ExplicitBitVect(len(bits))
            for position, bit in enumerate(bits):
                if bit:
                    vector.SetBit(position)
            path = analysis_dir / f"{ligand_name}_pose{analysis.pose_index}_ifp.dat"
            path.write_bytes(vector.ToBinary())
            written.append(path)
    except ImportError:
        logger.debug("RDKit unavailable: IFP .dat export skipped")
    return written
