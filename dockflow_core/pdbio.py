"""Pure-Python PDB and PDBQT reading / writing.

This module is deliberately dependency-free (only the standard library) so
that every other component - preparation, grid box, docking, analysis and
visualization - can rely on it in any environment, including CI and Docker.

Column conventions follow the PDB v3.3 specification and the AutoDock PDBQT
extension (partial charge in columns 71-76, AD4 atom type in columns 78-79).

Column map (1-based, inclusive):

    1 -  6  record type       "ATOM  " / "HETATM"
    7 - 11  atom serial
    12       space
   13 - 16  atom name
   17       altLoc
   18 - 20  residue name
   21       space
   22       chain id
   23 - 26  residue sequence number
   27       insertion code
   28 - 30  spaces
   31 - 38  x
   39 - 46  y
   47 - 54  z
   55 - 60  occupancy
   61 - 66  b-factor
   67 - 70  spaces (PDBQT)
   71 - 76  partial charge (PDBQT, %6.3f)
   77       space (PDBQT)
   78 - 79  AD4 atom type (PDBQT)
   77 - 78  element symbol (PDB, right-justified)
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import TextIO

__all__ = [
    "Atom",
    "PDBQTModel",
    "PDBQTFile",
    "VinaResultRecord",
    "AMINO_ACIDS",
    "NUCLEIC_ACIDS",
    "WATER_RESIDUES",
    "parse_pdb",
    "parse_pdb_models",
    "parse_pdbqt",
    "split_pdbqt_models",
    "write_pdb",
    "write_pdbqt",
    "format_pdb_line",
    "format_pdbqt_line",
    "parse_pdbqt_results",
    "element_from_atom_name",
]

# ---------------------------------------------------------------------------
# Residue dictionaries
# ---------------------------------------------------------------------------
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL", "HYP",  # common modified residues treated as polymer
}
NUCLEIC_ACIDS = {
    "A", "C", "G", "T", "U", "DA", "DC", "DG", "DT", "DU",
    "RA", "RC", "RG", "RU",
}
WATER_RESIDUES = {"HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SOL", "SPC"}

_TWO_LETTER_ELEMENTS = {
    "CL", "BR", "NA", "MG", "SI", "LI", "K", "CA", "MN", "FE", "ZN",
    "SE", "CU", "NI", "CO", "CD", "HG", "AL", "AU", "AG", "PT", "PB",
    "SN", "BA", "SR", "CS", "SB", "BI", "TE", "I",
}

_PDB_ID_RE = re.compile(r"^\d[A-Za-z0-9]{3}$")


def is_valid_pdb_id(pdb_id: str) -> bool:
    return bool(_PDB_ID_RE.match(pdb_id.strip()))


def is_polymer(resname: str) -> bool:
    return resname in AMINO_ACIDS or resname in NUCLEIC_ACIDS


# ---------------------------------------------------------------------------
# Atom record
# ---------------------------------------------------------------------------
@dataclass
class Atom:
    """A single atom from a PDB / PDBQT file (one model)."""

    serial: int = 0
    name: str = "UNK"
    altloc: str = ""
    resname: str = "UNK"
    chain: str = ""
    resseq: int = 0
    icode: str = ""
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    occupancy: float = 1.0
    bfactor: float = 0.0
    element: str = ""
    charge: float = 0.0  # PDBQT partial charge
    atom_type: str = ""  # PDBQT AD4 type
    record_type: str = "ATOM"
    model: int = 1

    # -- conveniences -------------------------------------------------------
    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def is_water(self) -> bool:
        return self.resname.strip() in WATER_RESIDUES

    @property
    def is_polymer(self) -> bool:
        return is_polymer(self.resname.strip())

    @property
    def is_metal(self) -> bool:
        return self.element.upper() in _METALS

    @property
    def is_hydrogen(self) -> bool:
        return self.element.upper() == "H"

    def residue_key(self) -> tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)

    def identity_key(self) -> tuple[str, str, int, str, str]:
        return (self.chain, self.resname, self.resseq, self.name.strip(), self.altloc)

    def clone(self) -> Atom:
        return replace(self)


_METALS = {
    "LI", "BE", "NA", "MG", "AL", "K", "CA", "SC", "TI", "V", "CR", "MN",
    "FE", "CO", "NI", "CU", "ZN", "GA", "RB", "SR", "Y", "ZR", "NB", "MO",
    "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "CS", "BA", "LA", "CE",
    "PR", "ND", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",
    "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI",
}


def element_from_atom_name(name: str, resname: str = "") -> str:
    """Guess the element from an atom name when column 77-78 is absent.

    Follows the usual PDB heuristics: strip digits, prefer two-letter
    symbols, but treat backbone "CA" of amino/nucleic acids as carbon.
    """
    cleaned = name.strip().lstrip("0123456789")
    if not cleaned:
        return "H"
    upper = cleaned.upper()
    # Backbone / side-chain atoms of polymers are single-letter elements.
    if resname and is_polymer(resname):
        if upper in {"CA", "CB", "CG", "CD", "CE", "CF", "CH", "CI", "CJ", "CK", "CL", "CM",
                     "CN", "CQ", "CR", "CS", "CT", "CU", "CV", "CX", "CY", "CZ", "C1", "C2",
                     "C3", "C4", "C5", "C6", "C7", "C8", "C9"}:
            return "C"
        if upper.startswith("N"):
            return "N"
        if upper.startswith("O"):
            return "O"
        if upper.startswith("S"):
            return "S"
    if len(upper) >= 2 and upper[:2] in _TWO_LETTER_ELEMENTS:
        return upper[:2].capitalize()
    return upper[0]


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------
def _safe_float(text: str, default: float = 0.0) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return default


def _safe_int(text: str, default: int = 0) -> int:
    try:
        return int(text)
    except (TypeError, ValueError):
        return default


def _parse_atom_line(line: str, model: int) -> Atom:
    line = line.ljust(54)
    name = line[12:16].strip()
    resname = line[17:20].strip()
    element_col = line[76:78].strip() if len(line) >= 78 else ""
    element = element_col if element_col else element_from_atom_name(name, resname)
    charge = 0.0
    atom_type = ""
    if len(line) >= 79:
        atom_type = line[77:79].strip()
        charge = _safe_float(line[70:76])
    elif len(line) >= 76:  # charge without type (some flavours)
        charge = _safe_float(line[70:76])
    if not atom_type:
        atom_type = element.upper()
    return Atom(
        serial=_safe_int(line[6:11]),
        name=name,
        altloc=line[16:17].strip(),
        resname=resname or "UNK",
        chain=line[21:22].strip(),
        resseq=_safe_int(line[22:26]),
        icode=line[26:27].strip(),
        x=_safe_float(line[30:38]),
        y=_safe_float(line[38:46]),
        z=_safe_float(line[46:54]),
        occupancy=_safe_float(line[54:60], 1.0) or 1.0,
        bfactor=_safe_float(line[60:66]),
        element=element.capitalize() if element else "C",
        charge=charge,
        atom_type=atom_type,
        record_type=line[0:6].strip() or "ATOM",
        model=model,
    )


def _iter_structure_lines(source: str | Path | TextIO | Iterable[str]) -> Iterator[str]:
    """Yield structure lines from a file path, raw text, file object or iterable.

    Long multi-line strings (raw structure text) are never mistaken for
    paths, and ``stat`` errors on exotic strings are swallowed.
    """
    if isinstance(source, Path):
        with open(source, encoding="utf-8", errors="replace") as fh:
            yield from fh
        return
    if isinstance(source, str):
        if "\n" not in source and len(source) < 1024:
            candidate = Path(source)
            try:
                if candidate.is_file():
                    with open(candidate, encoding="utf-8", errors="replace") as fh:
                        yield from fh
                    return
            except OSError:
                pass
        yield from source.splitlines()
        return
    if hasattr(source, "read"):
        yield from source.read().splitlines()  # type: ignore[union-attr]
        return
    yield from (str(line) for line in source)  # type: ignore[arg-type]


def parse_pdb_models(source: str | Path | TextIO | Iterable[str]) -> list[list[Atom]]:
    """Parse a PDB file into models (a file without MODEL records yields one)."""
    models: list[list[Atom]] = []
    current: list[Atom] = []
    model_index = 1
    in_model_block = False
    for line in _iter_structure_lines(source):
        record = line[0:6]
        if record == "MODEL ":
            if current and not in_model_block:
                models.append(current)
                current = []
            in_model_block = True
            model_index = _safe_int(line[10:14], model_index)
            continue
        if record == "ENDMDL":
            models.append(current)
            current = []
            in_model_block = False
            model_index += 1
            continue
        if record in ("ATOM  ", "HETATM"):
            current.append(_parse_atom_line(line, model_index))
        elif record == "END   ":
            break
    if current and (not models or in_model_block):
        models.append(current)
    if not models:
        models = [[]] if not current else models
    return [m for m in models if m is not None] or [[]]


def parse_pdb(source: str | Path | TextIO | Iterable[str], model: int = 1) -> list[Atom]:
    """Parse a PDB file and return the atoms of the requested model."""
    models = parse_pdb_models(source)
    if not models:
        return []
    if model == 1:
        return models[0]
    for atoms in models:
        if atoms and atoms[0].model == model:
            return atoms
    return models[0]


# ---------------------------------------------------------------------------
# PDBQT parsing
# ---------------------------------------------------------------------------
@dataclass
class VinaResultRecord:
    """``REMARK VINA RESULT:    -9.423      0.000      0.000`` of one pose."""

    model: int
    affinity: float
    rmsd_lb: float
    rmsd_ub: float


@dataclass
class PDBQTModel:
    """One MODEL/ENDMDL block of a PDBQT file."""

    index: int
    atoms: list[Atom] = field(default_factory=list)
    remarks: list[str] = field(default_factory=list)

    @property
    def vina_result(self) -> VinaResultRecord | None:
        for remark in self.remarks:
            if "VINA RESULT" in remark.upper():
                numbers = re.findall(r"-?\d+\.?\d*", remark.split(":", 1)[-1])
                if len(numbers) >= 3:
                    return VinaResultRecord(
                        model=self.index,
                        affinity=float(numbers[0]),
                        rmsd_lb=float(numbers[1]),
                        rmsd_ub=float(numbers[2]),
                    )
        return None


@dataclass
class PDBQTFile:
    """A parsed PDBQT file: models + global remarks + torsion info."""

    models: list[PDBQTModel] = field(default_factory=list)
    remarks: list[str] = field(default_factory=list)
    torsdof: int = 0

    @property
    def atoms(self) -> list[Atom]:
        return self.models[0].atoms if self.models else []

    def vina_results(self) -> list[VinaResultRecord]:
        results = [m.vina_result for m in self.models]
        return [r for r in results if r is not None]


def parse_pdbqt(source: str | Path | TextIO | Iterable[str]) -> PDBQTFile:
    """Parse a (possibly multi-pose) PDBQT file."""
    data = PDBQTFile()
    current = PDBQTModel(index=1)
    model_index = 1
    for line in _iter_structure_lines(source):
        stripped = line.strip()
        record = line[0:6]
        if record == "MODEL ":
            if current.atoms or current.remarks:
                data.models.append(current)
            model_index = _safe_int(line[10:14], model_index)
            current = PDBQTModel(index=model_index)
            continue
        if record == "ENDMDL":
            data.models.append(current)
            current = PDBQTModel(index=model_index + 1)
            model_index += 1
            continue
        if record in ("ATOM  ", "HETATM"):
            current.atoms.append(_parse_atom_line(line, model_index))
        elif stripped.startswith("REMARK"):
            current.remarks.append(stripped)
        elif stripped.startswith("TORSDOF"):
            parts = stripped.split()
            data.torsdof = _safe_int(parts[1]) if len(parts) > 1 else 0
        elif stripped in ("ROOT", "ENDROOT", "END", "TER"):
            continue
        # BRANCH/BRANCH lines are ignored for coordinate work.
    if current.atoms or current.remarks:
        if not (data.models and data.models[-1] is current):
            data.models.append(current)
    if not data.models:
        data.models = [PDBQTModel(index=1)]
    if data.models and not data.models[0].atoms and not data.models[0].remarks:
        # single-model file without MODEL records: atoms landed in one model
        pass
    return data


def parse_pdbqt_results(source: str | Path) -> list:
    """Convenience: extract Vina pose scores from a docked PDBQT file.

    Returns a list of :class:`dockflow_core.models.PoseRecord`.
    """
    from .models import PoseRecord

    data = parse_pdbqt(source)
    poses: list[PoseRecord] = []
    for model in data.models:
        result = model.vina_result
        if result is None:
            continue
        poses.append(
            PoseRecord(
                model=model.index,
                affinity=result.affinity,
                rmsd_lb=result.rmsd_lb,
                rmsd_ub=result.rmsd_ub,
                remarks=list(model.remarks),
            )
        )
    return poses


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------
def _atom_name_field(name: str) -> str:
    """Format the 4-character atom-name field (columns 13-16).

    Short names that start with a letter are right-shifted by one column
    (" CA "), full 4-character names and names starting with a digit fill
    the field from the left.
    """
    name = name.strip()
    if not name:
        return "    "
    if len(name) < 4 and name[0].isalpha():
        return f" {name[:3]:<3}"
    return f"{name[:4]:<4}"


def format_pdb_line(atom: Atom, serial: int | None = None) -> str:
    """Format a standard PDB ATOM/HETATM line (78 columns, no charge).

    Fixed-width columns follow the PDB v3.3 spec: altLoc (17), chain (22)
    and iCode (27) always occupy exactly one column, even when blank.
    """
    record = atom.record_type.strip() or "ATOM"
    if record not in ("ATOM", "HETATM"):
        record = "ATOM"
    serial = serial if serial is not None else atom.serial
    altloc = (atom.altloc or " ")[:1]
    chain = (atom.chain or " ")[:1]
    icode = (atom.icode or " ")[:1]
    head = (
        f"{record:<6}"
        f"{serial:5d}"
        f" "
        f"{_atom_name_field(atom.name)}"
        f"{altloc}"
        f"{atom.resname.strip()[:3]:>3}"
        f" {chain}"
        f"{atom.resseq % 10000:4d}"
        f"{icode}"
        f"   "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"{atom.occupancy:6.2f}{atom.bfactor:6.2f}"
    )
    assert len(head) == 66, f"PDB head line must be 66 columns, got {len(head)}"
    element = (atom.element or "C").strip().upper()
    if len(element) == 1:
        element = f" {element}"
    return f"{head}          {element:>2}"


def format_pdbqt_line(atom: Atom, serial: int | None = None) -> str:
    """Format a PDBQT ATOM/HETATM line.

    Partial charge occupies columns 71-76 (%6.3f) and the AD4 atom type
    columns 78-79, matching the output of MGLTools prepare_receptor4.py /
    prepare_ligand4.py and of AutoDock Vina itself.
    """
    line = format_pdb_line(atom, serial)
    atom_type = (atom.atom_type or atom.element or "C").strip()[:2]
    charge = atom.charge
    if charge > 99.999:
        charge = 99.999
    if charge < -99.999:
        charge = -99.999
    body = line[:66]
    return f"{body:<66}    {charge:6.3f} {atom_type:<2}"


def write_pdb(
    atoms: Sequence[Atom],
    path: str | Path,
    header: Iterable[str] = (),
    include_models: bool = False,
) -> Path:
    """Write atoms to a PDB file (TER/END appended automatically)."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = [f"REMARK   {h}" if not h.startswith("REMARK") else h for h in header]
    serial = 0
    last_chain: str | None = None
    last_residue: tuple | None = None
    models = {}
    if include_models:
        for atom in atoms:
            models.setdefault(atom.model, []).append(atom)
    groups: list[Iterable[Atom]] = (
        [models[m] for m in sorted(models)] if include_models else [atoms]
    )
    for group in groups:
        if include_models and len(groups) > 1:
            index = next(iter(group)).model
            lines.append(f"MODEL     {index:4d}")
        for atom in group:
            serial += 1
            lines.append(format_pdb_line(atom, serial))
            residue = atom.residue_key()
            chain_changed = last_chain is not None and atom.chain != last_chain
            residue_changed = last_residue is not None and residue != last_residue
            if (chain_changed or residue_changed) and atom.is_polymer:
                lines.append(_ter_line(serial, last_residue, last_chain))
            last_chain, last_residue = atom.chain, residue
        if include_models and len(groups) > 1:
            lines.append("ENDMDL")
    if last_residue is not None:
        lines.append(_ter_line(serial + 1, last_residue, last_chain))
    lines.append("END")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def _ter_line(serial: int, residue: tuple | None, chain: str | None) -> str:
    resname, resseq, icode = (residue if residue else ("UNK", 0, ""))
    chain = chain or ""
    return (
        f"TER   {serial:5d}      {resname:>3} "
        f"{chain[:1]}{resseq % 10000:4d}{icode[:1]}"
    )


def write_pdbqt(
    atoms: Sequence[Atom],
    path: str | Path,
    remarks: Iterable[str] = (),
    header: Iterable[str] = (),
    models: Sequence[Sequence[Atom]] | None = None,
    torsdof: int | None = None,
    root_block: bool = False,
) -> Path:
    """Write atoms to a PDBQT file.

    Args:
        atoms: atoms of the (single) structure; ignored when ``models`` given.
        path: output file path.
        remarks: REMARK lines written before the atoms.
        header: extra header lines written verbatim at the top.
        models: optional list of atom lists (one per MODEL block) for
            multi-pose output (ROOT/BRANCH blocks are not reconstructed).
        torsdof: optional TORSDOF value appended at the end.
        root_block: emit ROOT/ENDROOT around the atom block.  AutoDock Vina
            requires this for **ligand** PDBQT files but rejects it in rigid
            **receptor** files, so it is opt-in (split pose files use it).
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = list(header)
    lines += [f"REMARK  {r}" if not str(r).startswith("REMARK") else str(r) for r in remarks]
    groups: list[Sequence[Atom]] = list(models) if models is not None else [atoms]
    serial = 0
    for group_index, group in enumerate(groups):
        if len(groups) > 1:
            lines.append(f"MODEL     {group_index + 1:4d}")
        if root_block:
            lines.append("ROOT")
        for atom in group:
            serial += 1
            lines.append(format_pdbqt_line(atom, serial))
        if root_block:
            lines.append("ENDROOT")
        if len(groups) > 1:
            lines.append("ENDMDL")
    if torsdof is not None:
        lines.append(f"TORSDOF {int(torsdof)}")
    lines.append("END")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def split_pdbqt_models(source: str | Path, out_dir: str | Path, basename: str) -> list[Path]:
    """Split a multi-pose PDBQT into one file per pose.

    Each output keeps its own REMARK VINA RESULT (and other remarks), which
    is what visualization tools and per-pose analysis want.
    """
    data = parse_pdbqt(source)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for model in data.models:
        if not model.atoms:
            continue
        path = out / f"{basename}_pose{model.index:02d}.pdbqt"
        write_pdbqt(model.atoms, path, remarks=model.remarks, torsdof=data.torsdof,
                    root_block=True)
        written.append(path)
    return written


# ---------------------------------------------------------------------------
# Bulk helpers used across the code base
# ---------------------------------------------------------------------------
def filter_atoms(
    atoms: Sequence[Atom],
    chains: Sequence[str] | None = None,
    drop_waters: bool = True,
    keep_hetero: bool = False,
    keep_resnames: Sequence[str] | None = None,
    records: Sequence[str] | None = None,
) -> list[Atom]:
    """Filter an atom list by chain / water / hetero rules (order preserved)."""
    keep_res = {r.strip().upper() for r in (keep_resnames or ())}
    chain_set = {c.strip() for c in chains} if chains else None
    record_set = {r.strip().upper() for r in records} if records else None
    result: list[Atom] = []
    for atom in atoms:
        if chain_set is not None and atom.chain not in chain_set:
            continue
        if drop_waters and atom.is_water:
            continue
        if not keep_hetero and not atom.is_polymer:
            if atom.resname.strip().upper() not in keep_res:
                continue
        if record_set is not None and atom.record_type.upper() not in record_set:
            continue
        result.append(atom)
    return result


def het_resnames(atoms: Sequence[Atom], exclude_water: bool = True) -> list[str]:
    """Unique HETATM residue names in order of appearance (for ligand pickers)."""
    seen: set[str] = set()
    order: list[str] = []
    for atom in atoms:
        name = atom.resname.strip()
        if not name or not name.isalnum():
            continue
        if atom.is_polymer:
            continue
        if exclude_water and atom.is_water:
            continue
        if name not in seen:
            seen.add(name)
            order.append(name)
    return order


def renumber_serials(atoms: Sequence[Atom]) -> list[Atom]:
    """Return copies of atoms with serial numbers 1..N."""
    return [replace(atom, serial=i + 1) for i, atom in enumerate(atoms)]

