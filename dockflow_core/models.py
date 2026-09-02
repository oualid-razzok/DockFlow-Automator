"""Data models shared between the pipeline, engine, GUI and reports."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class ProteinRecord:
    """A downloaded / loaded macromolecular target."""

    identifier: str = ""
    source: str = "unknown"  # pdb | uniprot-alphafold | file
    path: Path | None = None
    title: str = ""
    resolution: float | None = None
    uniprot_ids: list[str] = field(default_factory=list)
    ligand_codes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["path"] = str(self.path) if self.path else None
        return data


@dataclass
class LigandRecord:
    """One ligand in a docking run, tracked from download to poses."""

    identifier: str = ""
    source: str = "smiles"  # smiles | file | pubchem | zinc | pdb_ligand
    value: str = ""  # the SMILES / id / path
    path: Path | None = None  # downloaded / intermediate 3D structure
    pdbqt_path: Path | None = None
    status: str = "pending"  # pending | downloaded | prepared | docked | error
    error: str | None = None
    num_rotatable_bonds: int | None = None

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["path"] = str(self.path) if self.path else None
        data["pdbqt_path"] = str(self.pdbqt_path) if self.pdbqt_path else None
        return data


@dataclass
class PoseRecord:
    """One docking pose of one ligand."""

    model: int = 1
    affinity: float = 0.0  # kcal/mol
    rmsd_lb: float = 0.0  # lower-bound RMSD to best mode
    rmsd_ub: float = 0.0  # upper-bound RMSD to best mode
    remarks: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class DockingResult:
    """Result of docking a single ligand with AutoDock Vina."""

    ligand: LigandRecord | None = None
    ligand_name: str = ""
    poses: list[PoseRecord] = field(default_factory=list)
    out_path: Path | None = None  # multi-pose PDBQT
    log_text: str = ""
    log_path: Path | None = None
    runtime: float = 0.0
    backend: str = "unknown"
    error: str | None = None

    def __post_init__(self) -> None:
        if not self.ligand_name and self.ligand is not None:
            self.ligand_name = self.ligand.identifier

    @property
    def ok(self) -> bool:
        return self.error is None and bool(self.poses)

    @property
    def best_affinity(self) -> float | None:
        return self.poses[0].affinity if self.poses else None

    def to_dict(self) -> dict[str, Any]:
        return {
            "ligand": self.ligand.to_dict() if self.ligand else None,
            "ligand_name": self.ligand_name,
            "poses": [p.to_dict() for p in self.poses],
            "out_path": str(self.out_path) if self.out_path else None,
            "log_path": str(self.log_path) if self.log_path else None,
            "runtime": round(self.runtime, 3),
            "backend": self.backend,
            "error": self.error,
            "best_affinity": self.best_affinity,
        }


@dataclass
class Contact:
    """A ligand-receptor contact detected during analysis."""

    ligand_atom_index: int = 0
    ligand_atom_name: str = ""
    ligand_atom_type: str = ""
    receptor_atom_name: str = ""
    receptor_resname: str = ""
    receptor_chain: str = ""
    receptor_resseq: int = 0
    receptor_atom_type: str = ""
    distance: float = 0.0
    kind: str = "close"  # hbond | hydrophobic | ionic | metal | close

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
