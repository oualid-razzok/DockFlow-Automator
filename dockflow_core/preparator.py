"""Receptor and ligand preparation for AutoDock Vina.

This module re-implements the parameterization logic of the legacy MGLTools
scripts (``prepare_receptor4.py`` and ``prepare_ligand4.py``) on top of
modern, maintained toolkits:

===========================  ==============================================
MGLTools flag                DockFlow equivalent
===========================  ==============================================
``prepare_receptor4.py``
  ``-A hydrogens``           ``ReceptorPrepOptions.add_hydrogens=True``
  ``-U nphs`` (default)      ``ReceptorPrepOptions.merge_nonpolar_h=True``
  ``-U lps``                 lone pairs removed (always; Vina ignores them)
  ``-U altloc``              ``ReceptorPrepOptions.altloc="A"|"best"``
  ``-C`` (keep chains)       ``ReceptorPrepOptions.chains=[...]``
  ``-w`` (keep water)        ``ReceptorPrepOptions.keep_water=True``
  ``-e`` (no charges)        ``ReceptorPrepOptions.charge_model="zero"``
  ``-r``/``-o``              ``prepare(input, output_dir)`` paths
``prepare_ligand4.py``
  ``-l`` / ``-o``            input / output paths of ``LigandPreparator``
  ``-A checkhydrogens``      hydrogens always added by RDKit/Meeko
  ``-U nphs``                Meeko default (non-polar H merged)
  ``-B bonds``               Meeko torsion detection (O-O bonds etc.)
  ``-Z`` (ph)                ``LigandPrepOptions.protonate=True`` (dimorphite-dl)
  ``-F flexible_amide``      Meeko handles flexibility internally
===========================  ==============================================

Hydrogen addition and Gasteiger charging are performed by one of several
interchangeable "hydration engines", selected automatically:

1. ``openbabel``   - OpenBabel Python bindings (openbabel-wheel / conda-forge)
2. ``rdkit``       - RDKit (PDB reader + AddHs + Gasteiger)
3. ``openbabel-cli``- the ``obabel`` executable via a MOL2 round-trip
4. ``none``        - dependency-free fallback: no new hydrogens, zero
                     charges, distance-based bonds (for CI/testing only)

With ``engine="auto"`` every engine failure (a broken/incompatible optional
dependency such as an ``AttributeError`` inside the OpenBabel bindings, an
``ImportError``, or a parse failure) is caught, reported as a WARNING naming
the failed engine and the fallback chosen, and the next engine in the chain
is tried instead - the run must degrade gracefully, never crash.
"""

from __future__ import annotations

import tempfile
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .pdbio import (
    Atom,
    is_polymer,
    parse_pdb,
    renumber_serials,
    write_pdb,
    write_pdbqt,
)
from .utils import (
    DockFlowError,
    ensure_dir,
    get_logger,
    is_importable,
    run_command,
    which,
)

logger = get_logger("preparator")

__all__ = [
    "PreparationError",
    "ReceptorPrepOptions",
    "ReceptorPrepResult",
    "ReceptorPreparator",
    "LigandPrepOptions",
    "LigandPrepResult",
    "LigandPreparator",
    "COVALENT_RADII",
    "ELEMENT_SYMBOLS",
    "element_symbol",
]

class PreparationError(DockFlowError):
    """Structure preparation failed."""


# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
@dataclass
class ReceptorPrepOptions:
    """Options mirroring ``prepare_receptor4.py``."""

    chains: list[str] | None = None          # keep only these chains (None = all)
    keep_water: bool = False                 # ``-w``
    keep_hetero: bool = False                # keep every HETATM (cofactors, metals)
    keep_resnames: list[str] = field(default_factory=list)  # always keep these HETATM resnames
    remove_resnames: list[str] = field(default_factory=list)
    remove_residues: list[tuple[str, int]] = field(default_factory=list)  # (chain, resseq)
    altloc: str = "best"                     # "best" | "A" | "B" ... | "" (keep all)
    add_hydrogens: bool = True               # ``-A hydrogens``
    merge_nonpolar_h: bool = True            # ``-U nphs``
    charge_model: str = "gasteiger"          # "gasteiger" | "zero"  (``-e``)
    engine: str = "auto"                     # auto|openbabel|openbabel-cli|rdkit|none
    keep_unknown_types: bool = True          # keep atoms with non-AD4 elements (warn)


@dataclass
class ReceptorPrepResult:
    """Outcome of receptor preparation."""

    pdbqt_path: Path | None = None
    pdb_path: Path | None = None
    engine: str = "none"
    atoms_in: int = 0
    atoms_out: int = 0
    waters_removed: int = 0
    hetero_removed: int = 0
    removed_resnames: dict[str, int] = field(default_factory=dict)
    hydrogens_added: int = 0
    warnings: list[str] = field(default_factory=list)
    # Structure quality checks (work order items 6, 7, 10):
    missing_residues: list[dict] = field(default_factory=list)
    disulfides: list[dict] = field(default_factory=list)
    reduced_cysteines: list[dict] = field(default_factory=list)
    # Flexible-receptor split (peer item 8): rigid + flex PDBQT pair.
    rigid_pdbqt_path: Path | None = None
    flex_pdbqt_path: Path | None = None
    flexible_residues: list[dict] = field(default_factory=list)
    metal_coordination: list[dict] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.pdbqt_path is not None and self.pdbqt_path.is_file()


@dataclass
class LigandPrepOptions:
    """Options mirroring ``prepare_ligand4.py`` (via Meeko)."""

    embed_3d: bool = True                    # ETKDG embedding when no conformer exists
    minimize: bool = True                    # MMFF94 (UFF fallback) minimisation
    minimize_steps: int = 500
    keep_largest_fragment: bool = True
    remove_salts: bool = True
    # pH-aware protonation (peer item 11): True = major microspecies at
    # ``ph`` via dimorphite-dl; "enumerate" = dock every populated state.
    protonate: bool | str = False
    ph: float = 7.4                          # target pH +/- 0.5
    # Tautomer enumeration (peer item 12)
    enumerate_tautomers: bool = False
    max_tautomers: int = 16
    # Stereochemistry (peer item 13)
    enumerate_stereoisomers: bool = False
    max_stereoisomers: int = 8
    charge_model: str = "gasteiger"
    random_seed: int = 42


@dataclass
class LigandPrepResult:
    """Outcome of ligand preparation (one ligand)."""

    identifier: str = ""
    pdbqt_path: Path | None = None
    sdf_path: Path | None = None
    smiles: str = ""
    num_rotatable_bonds: int | None = None
    num_heavy_atoms: int | None = None
    num_atoms: int | None = None
    charge_model: str = "gasteiger"
    engine: str = "meeko"
    # Provenance of the scientific decisions (audit item 26):
    input_had_explicit_hydrogens: bool | None = None
    protonation_source: str = "input protonation state (no pH adjustment)"
    num_protonation_variants: int = 1
    warnings: list[str] = field(default_factory=list)
    error: str | None = None
    # State-chemistry provenance (work order items 11-14):
    variants: list[dict] = field(default_factory=list)
    protonation: dict | None = None
    tautomers: dict | None = None
    stereocenters: dict | None = None
    salt_stripping: dict | None = None

    @property
    def ok(self) -> bool:
        return self.error is None and self.pdbqt_path is not None


# ---------------------------------------------------------------------------
# Element table (Z -> symbol)
# ---------------------------------------------------------------------------
# Local periodic table so receptor preparation never depends on a toolkit
# helper that optional wheels may not ship (modern openbabel wheels, e.g.
# 3.1.1.23, removed the OBElementTable SWIG wrapper entirely).
ELEMENT_SYMBOLS: dict[int, str] = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
    9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
    16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti",
    23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu",
    30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr",
    37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc",
    44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La",
    58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd",
    65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu",
    72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt",
    79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At",
    86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U",
    93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es",
    100: "Fm", 101: "Md", 102: "No", 103: "Lr",
}


def element_symbol(atomic_num: int) -> str:
    """Element symbol for an atomic number from the local table ("X" if unknown)."""
    return ELEMENT_SYMBOLS.get(int(atomic_num), "X")


# ---------------------------------------------------------------------------
# Engine atom graph (toolkit-independent)
# ---------------------------------------------------------------------------
COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07,
    "F": 0.57, "CL": 1.02, "BR": 1.20, "I": 1.39, "NA": 1.66, "MG": 1.41,
    "K": 2.03, "CA": 1.76, "MN": 1.39, "FE": 1.32, "ZN": 1.22, "CU": 1.32,
    "NI": 1.24, "CO": 1.26, "SE": 1.20, "CD": 1.44, "HG": 1.44, "B": 0.84,
    "SI": 1.11, "LI": 1.28, "AL": 1.21, "AU": 1.36, "AG": 1.45, "PT": 1.75,
}


@dataclass
class EngineAtom:
    """Toolkit-independent atom with connectivity, charge and flags."""

    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    element: str = "C"
    charge: float = 0.0
    name: str = "UNK"
    resname: str = "UNK"
    chain: str = ""
    resseq: int = 0
    icode: str = ""
    record_type: str = "ATOM"
    occupancy: float = 1.0
    bfactor: float = 0.0
    is_aromatic: bool = False
    is_donor: bool = False
    is_acceptor: bool = False
    neighbors: list[int] = field(default_factory=list)
    atom_type: str = ""

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def is_hydrogen(self) -> bool:
        return self.element.upper() == "H"

    def heavy_neighbors(self, graph: Sequence[EngineAtom]) -> list[int]:
        return [j for j in self.neighbors if not graph[j].is_hydrogen]


def bond_by_distance(a: EngineAtom, b: EngineAtom) -> bool:
    """Distance-based bond heuristic (covalent radii + tolerance)."""
    if a.element == b.element and a.is_hydrogen and b.is_hydrogen:
        return False
    ra = COVALENT_RADII.get(a.element.upper(), 0.8)
    rb = COVALENT_RADII.get(b.element.upper(), 0.8)
    dx, dy, dz = a.x - b.x, a.y - b.y, a.z - b.z
    distance = (dx * dx + dy * dy + dz * dz) ** 0.5
    if a.is_hydrogen or b.is_hydrogen:
        return distance < (ra + rb + 0.35)
    return distance < (ra + rb + 0.45)


def build_graph(atoms: Sequence[Atom]) -> list[EngineAtom]:
    """Convert parsed PDB atoms into an EngineAtom graph bonded by distance."""
    graph = [
        EngineAtom(
            x=a.x, y=a.y, z=a.z,
            element=a.element.upper(),
            charge=0.0,
            name=a.name.strip(),
            resname=a.resname.strip(),
            chain=a.chain,
            resseq=a.resseq,
            icode=a.icode,
            record_type=a.record_type,
            occupancy=a.occupancy,
            bfactor=a.bfactor,
        )
        for a in atoms
    ]
    n = len(graph)
    for i in range(n):
        for j in range(i + 1, n):
            if bond_by_distance(graph[i], graph[j]):
                graph[i].neighbors.append(j)
                graph[j].neighbors.append(i)
    return graph


# ---------------------------------------------------------------------------
# Hydration engines
# ---------------------------------------------------------------------------
class BaseEngine:
    """Interface for hydrogen-addition / charging engines."""

    name = "base"

    @staticmethod
    def available() -> bool:  # pragma: no cover - overridden
        return False

    def process(self, atoms: Sequence[Atom], charge_model: str) -> list[EngineAtom]:
        raise NotImplementedError


class PassThroughEngine(BaseEngine):
    """Dependency-free fallback: keeps existing atoms, zero charges.

    Bonds are perceived from geometry so that AD4 typing and non-polar
    hydrogen merging still work.  Use only when no chemistry toolkit is
    installed - charges will be 0.0 and Vina scoring will be degraded.
    """

    name = "none"

    @staticmethod
    def available() -> bool:
        return True

    def process(self, atoms: Sequence[Atom], charge_model: str) -> list[EngineAtom]:
        graph = build_graph(atoms)
        if charge_model == "zero":
            return graph
        for atom in graph:  # no charge engine available
            atom.charge = 0.0
        return graph


class RDKitEngine(BaseEngine):
    """RDKit: PDB reading, AddHs, Gasteiger charges, aromatic perception."""

    name = "rdkit"

    @staticmethod
    def available() -> bool:
        return is_importable("rdkit")

    def process(self, atoms: Sequence[Atom], charge_model: str) -> list[EngineAtom]:
        from rdkit import Chem
        from rdkit.Chem import rdPartialCharges

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_pdb = Path(tmpdir) / "input.pdb"
            write_pdb(atoms, tmp_pdb)
            mol = Chem.MolFromPDBFile(
                str(tmp_pdb), removeHs=False, sanitize=True, proximityBonding=True
            )
            if mol is None:
                mol = Chem.MolFromPDBFile(
                    str(tmp_pdb), removeHs=False, sanitize=False, proximityBonding=True
                )
                if mol is None:
                    raise PreparationError("RDKit could not parse the filtered receptor PDB")
                try:
                    mol.UpdatePropertyCache(strict=False)
                    Chem.SanitizeMol(
                        mol,
                        Chem.SanitizeFlags.SANITIZE_ALL
                        ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                        ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
                    )
                except Exception as exc:  # noqa: BLE001 - RDKit raises many types
                    raise PreparationError(f"RDKit sanitisation failed: {exc}") from exc
            mol = Chem.AddHs(mol, addCoords=True)
            if charge_model == "gasteiger":
                try:
                    rdPartialCharges.ComputeGasteigerCharges(mol)
                except Exception as exc:  # noqa: BLE001
                    logger.warning("Gasteiger charges failed: %s", exc)
            return self._graph_from_rdkit(mol, charge_model)

    @staticmethod
    def _graph_from_rdkit(mol, charge_model: str) -> list[EngineAtom]:

        conf = mol.GetConformer()
        graph: list[EngineAtom] = []
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)
            mono = atom.GetMonomerInfo()
            name, resname, chain, resseq, icode, record = (
                f"{atom.GetSymbol()}{idx + 1}", "UNK", "", 0, "", "HETATM"
            )
            if mono is not None and hasattr(mono, "GetResidueName"):
                name = mono.GetName().strip() or name
                resname = mono.GetResidueName().strip() or "UNK"
                chain = (mono.GetChainId() or "").strip()
                resseq = mono.GetResidueNumber() or 0
                icode = (mono.GetInsertionCode() or "").strip()
                record = "ATOM" if is_polymer(resname) else "HETATM"
            charge = 0.0
            if charge_model == "gasteiger" and atom.HasProp("_GasteigerCharge"):
                try:
                    charge = float(atom.GetProp("_GasteigerCharge"))
                except ValueError:
                    charge = 0.0
            element = atom.GetSymbol().upper()
            neighbors_h = any(
                neighbor.GetAtomicNum() == 1 for neighbor in atom.GetNeighbors()
            )
            symbol = atom.GetSymbol()
            acceptor = symbol in ("O", "S") and atom.GetFormalCharge() >= 0
            if symbol == "N":
                acceptor = not neighbors_h
            donor = symbol in ("N", "O") and neighbors_h
            graph.append(
                EngineAtom(
                    x=float(pos.x), y=float(pos.y), z=float(pos.z),
                    element=element, charge=charge, name=name, resname=resname,
                    chain=chain, resseq=int(resseq), icode=icode, record_type=record,
                    is_aromatic=bool(atom.GetIsAromatic()),
                    is_donor=bool(donor), is_acceptor=bool(acceptor),
                    neighbors=[n.GetIdx() for n in atom.GetNeighbors()],
                )
            )
        # Added hydrogens may lack monomer info: inherit from heavy neighbour.
        for atom in graph:
            if atom.resname != "UNK" or not atom.is_hydrogen:
                continue
            for j in atom.neighbors:
                if not graph[j].is_hydrogen and graph[j].resname != "UNK":
                    target = graph[j]
                    atom.resname, atom.chain = target.resname, target.chain
                    atom.resseq, atom.icode = target.resseq, target.icode
                    atom.record_type = target.record_type
                    break
        return graph


class OpenBabelPythonEngine(BaseEngine):
    """OpenBabel Python bindings (openbabel-wheel or conda-forge)."""

    name = "openbabel"

    @staticmethod
    def available() -> bool:
        return is_importable("openbabel")

    def process(self, atoms: Sequence[Atom], charge_model: str) -> list[EngineAtom]:
        from openbabel import openbabel as ob

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_pdb = Path(tmpdir) / "input.pdb"
            write_pdb(atoms, tmp_pdb)
            conv = ob.OBConversion()
            conv.SetInAndOutFormats("pdb", "pdb")
            mol = ob.OBMol()
            if not conv.ReadFile(mol, str(tmp_pdb)):
                raise PreparationError("OpenBabel could not read the filtered receptor PDB")
            mol.AddHydrogens(False, False, 7.4)
            if charge_model == "gasteiger":
                charge_model_obj = ob.OBChargeModel.FindType("gasteiger")
                if charge_model_obj is not None:
                    charge_model_obj.ComputeCharges(mol)
            return self._graph_from_obmol(mol, charge_model)

    @staticmethod
    def _graph_from_obmol(mol, charge_model: str) -> list[EngineAtom]:
        from openbabel import openbabel as ob

        graph: list[EngineAtom] = []
        idx_map: dict[int, int] = {}
        for atom in ob.OBMolAtomIter(mol):
            index = atom.GetIdx() - 1  # OB is 1-based
            idx_map[atom.GetIdx()] = index
            # Local Z -> symbol table: OBElementTable was removed from modern
            # openbabel wheels (3.1.1.23+); module-level GetSymbol is not
            # available in every distribution either.
            element = element_symbol(atom.GetAtomicNum()).upper()
            res = atom.GetResidue()
            name = "UNK"
            resname, chain, resseq, icode = "UNK", "", 0, ""
            if res is not None:
                try:
                    name = res.GetAtomID(atom).strip() or f"{element}{index + 1}"
                except Exception:  # noqa: BLE001
                    name = f"{element}{index + 1}"
                resname = (res.GetName() or "UNK").strip()
                # OBResidue.GetChainID() was dropped from the SWIG wrapper in
                # modern wheels (3.1.1.23); GetChain() returns the same string
                # id and exists in both old and new binding generations.
                try:
                    chain = (res.GetChain() or "").strip()
                except AttributeError:
                    chain = (res.GetChainID() or "").strip()
                resseq = int(res.GetNum() or 0)
                get_icode = getattr(res, "GetInsertionCode", None)
                if callable(get_icode):
                    # openbabel 3.1.1.23 returns the C terminator '\x00'
                    # for empty insertion codes - strip it or every NUL
                    # leaks into written PDBQT files.
                    icode = (get_icode() or "").replace("\x00", "").strip()
            charge = 0.0
            if charge_model == "gasteiger":
                try:
                    charge = float(atom.GetPartialCharge())
                except (TypeError, ValueError):
                    charge = 0.0
            is_donor = bool(getattr(atom, "IsHbondDonor", lambda: False)())
            is_acceptor = bool(getattr(atom, "IsHbondAcceptor", lambda: False)())
            graph.append(
                EngineAtom(
                    x=atom.GetX(), y=atom.GetY(), z=atom.GetZ(),
                    element=element, charge=charge, name=name, resname=resname,
                    chain=chain, resseq=resseq, icode=icode,
                    record_type="ATOM" if is_polymer(resname) else "HETATM",
                    is_aromatic=bool(atom.IsAromatic()),
                    is_donor=is_donor, is_acceptor=is_acceptor,
                    neighbors=[],
                )
            )
        for atom in ob.OBMolAtomIter(mol):
            i = idx_map[atom.GetIdx()]
            for neighbor in ob.OBAtomAtomIter(atom):
                j = idx_map.get(neighbor.GetIdx())
                if j is not None and j != i:
                    graph[i].neighbors.append(j)
        return graph


def _parse_mol2(text: str) -> tuple[list[EngineAtom], list[tuple[int, int, str]]]:
    """Parse a (subset of the) MOL2 format produced by ``obabel``."""
    graph: list[EngineAtom] = []
    bonds: list[tuple[int, int, str]] = []
    section = ""
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith("@<TRIPOS>"):
            section = line
            continue
        if not line or section.endswith("MOLECULE") or section.endswith("INFO"):
            continue
        if section.endswith("ATOM"):
            parts = line.split()
            if len(parts) < 8:
                continue
            atom_id, name, x, y, z, sybyl, subst_id, subst_name = parts[:8]
            charge = float(parts[8]) if len(parts) > 8 else 0.0
            resname = "UNK"
            resnum = 0
            sub = subst_name
            digits = sub.lstrip("0123456789")
            if digits:
                resname = digits[:3].upper()
            try:
                resnum = int(sub[: len(sub) - len(digits)] or subst_id)
            except ValueError:
                resnum = int(subst_id)
            element = sybyl.split(".")[0].upper()
            if element not in COVALENT_RADII and len(element) > 2:
                element = element[:2]
            graph.append(
                EngineAtom(
                    x=float(x), y=float(y), z=float(z),
                    element=element, charge=charge, name=name, resname=resname,
                    chain="", resseq=resnum,
                    record_type="ATOM" if is_polymer(resname) else "HETATM",
                    is_aromatic=sybyl.endswith(".ar"),
                )
            )
        elif section.endswith("BOND"):
            parts = line.split()
            if len(parts) >= 4:
                bonds.append((int(parts[1]) - 1, int(parts[2]) - 1, parts[3]))
    for i, j, _kind in bonds:
        if 0 <= i < len(graph) and 0 <= j < len(graph):
            graph[i].neighbors.append(j)
            graph[j].neighbors.append(i)
    return graph, bonds


class OpenBabelCLIEngine(BaseEngine):
    """``obabel`` executable round-trip (PDB -> MOL2 with charges).

    Chain identifiers are not carried through MOL2; use the Python bindings
    or RDKit engines when chains matter.
    """

    name = "openbabel-cli"

    @staticmethod
    def available() -> bool:
        return which("obabel") is not None

    def process(self, atoms: Sequence[Atom], charge_model: str) -> list[EngineAtom]:
        executable = which("obabel")
        if executable is None:
            raise PreparationError("obabel executable not found on PATH")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_pdb = Path(tmpdir) / "input.pdb"
            out_mol2 = Path(tmpdir) / "out.mol2"
            write_pdb(atoms, tmp_pdb)
            args = [executable, str(tmp_pdb), "-O", str(out_mol2), "-h"]
            if charge_model == "gasteiger":
                args += ["--partialcharge", "gasteiger"]
            result = run_command(args, timeout=120)
            if not result.ok or not out_mol2.is_file():
                raise PreparationError(
                    f"obabel failed (rc={result.returncode}): {result.stdout[-400:]}"
                )
            graph, _bonds = _parse_mol2(out_mol2.read_text(encoding="utf-8", errors="replace"))
            if not graph:
                raise PreparationError("obabel produced an empty MOL2 file")
            return graph


_ENGINES: dict[str, type[BaseEngine]] = {
    "openbabel": OpenBabelPythonEngine,
    "rdkit": RDKitEngine,
    "openbabel-cli": OpenBabelCLIEngine,
    "none": PassThroughEngine,
}

_AUTO_ENGINE_ORDER = ("openbabel", "rdkit", "openbabel-cli", "none")


def _engine_chain(preferred: str = "auto") -> list[BaseEngine]:
    """Ordered engines to attempt for one preparation request.

    ``auto`` returns every available engine in priority order, always ending
    with the dependency-free ``none`` fallback.  An explicit engine returns
    exactly that engine so a specific request fails loudly instead of
    silently substituting a different chemistry.
    """
    if preferred != "auto":
        return [select_engine(preferred)]
    chain: list[BaseEngine] = []
    for name in _AUTO_ENGINE_ORDER[:-1]:
        engine_cls = _ENGINES[name]
        if engine_cls.available():
            chain.append(engine_cls())
    chain.append(PassThroughEngine())  # always available last resort
    return chain


# ---------------------------------------------------------------------------
# Engine comparability (peer item 4)
# ---------------------------------------------------------------------------
# Engines listed in ``comparable_to`` produce scientifically equivalent
# receptors for the same input (same toolkit, different binding); engines
# in ``incomparable_to`` produce hydrogens/charges that differ enough that
# docking scores are NOT comparable across them (README section "Known
# scientific differences per preparation engine").
_ENGINE_COMPARABILITY: dict[str, tuple[list[str], list[str]]] = {
    "openbabel": (["openbabel", "openbabel-cli"], ["rdkit", "none"]),
    "openbabel-cli": (["openbabel", "openbabel-cli"], ["rdkit", "none"]),
    "rdkit": (["rdkit"], ["openbabel", "openbabel-cli", "none"]),
    "none": ([], ["openbabel", "openbabel-cli", "rdkit"]),
}


def engine_comparability(engine_name: str) -> dict:
    """Which engines this run's results are (in)comparable with (item 4).

    ``engine_name`` is the *resolved* engine recorded in the manifest
    (``"rdkit"``, ``"openbabel"``, ``"none (hydrogens disabled)"`` ...);
    the base token decides the comparability sets.
    """
    base = (engine_name or "none").split()[0].strip().lower()
    comparable, incomparable = _ENGINE_COMPARABILITY.get(base, ([base], []))
    return {
        "resolved_engine": base,
        "comparable_to": comparable,
        "incomparable_to": incomparable,
    }


def select_engine(preferred: str = "auto") -> BaseEngine:
    """Pick a hydration engine by name or automatically by availability."""
    if preferred != "auto":
        engine_cls = _ENGINES.get(preferred)
        if engine_cls is None:
            raise PreparationError(
                f"unknown preparation engine {preferred!r} "
                f"(valid: {', '.join(_ENGINES)})"
            )
        if not engine_cls.available():
            raise PreparationError(
                f"engine {preferred!r} is not available in this environment"
            )
        return engine_cls()
    for name in ("openbabel", "rdkit", "openbabel-cli", "none"):
        engine_cls = _ENGINES[name]
        if engine_cls.available():
            return engine_cls()
    return PassThroughEngine()  # pragma: no cover - unreachable


def resolved_engine_name(preferred: str = "auto") -> str:
    """The engine ``auto`` (or an explicit request) will actually use.

    Used by the degraded-mode gates (final pass, item 8): the CLI
    ``--allow-degraded`` flag and the GUI modal both need to know BEFORE
    the run whether the resolved engine is the dependency-free ``none``
    fallback.  Raises :class:`PreparationError` for unknown/unavailable
    explicit engines - the caller decides whether that is fatal.
    """
    return select_engine(preferred).name


# ---------------------------------------------------------------------------
# AD4 atom typing (mirrors MGLTools' fuse/smarts typing, pragmatic subset)
# ---------------------------------------------------------------------------
_VINA_ATOM_TYPES = {
    "H", "HD", "C", "A", "N", "NA", "NS", "O", "OA", "OS",
    "F", "Mg", "P", "S", "SA", "Cl", "Ca", "Mn", "Fe", "Zn", "Br", "I",
}


def assign_ad4_types(graph: Sequence[EngineAtom], keep_unknown: bool = True) -> list[str]:
    """Assign AutoDock/Vina atom types to every EngineAtom.

    Returns the list of warnings for atoms with types outside the Vina set.
    """
    warnings: list[str] = []
    unknown: set[str] = set()
    for i, atom in enumerate(graph):
        element = atom.element.upper()
        if element == "H":
            heavy = next((graph[j] for j in atom.neighbors if not graph[j].is_hydrogen), None)
            heavy_element = heavy.element.upper() if heavy else "C"
            if heavy_element in ("N", "O", "S"):
                # S-attached hydrogens are typed HD (not the AD4-only HS
                # type): Vina 1.2.x PDBQT parsers reject HS, and OpenBabel
                # itself types S-H as HD.  Vina has no separate S-H donor
                # term, so the H is treated as a generic polar hydrogen.
                atom.atom_type = "HD"
            else:
                atom.atom_type = "H"
        elif element == "C":
            atom.atom_type = "A" if atom.is_aromatic else "C"
        elif element == "N":
            has_h = any(graph[j].is_hydrogen for j in atom.neighbors)
            if has_h or (not atom.is_acceptor and not is_n_acceptor(graph, i)):
                atom.atom_type = "N"
            else:
                atom.atom_type = "NA"
        elif element == "O":
            has_phosphorus_neighbor = any(
                graph[j].element.upper() == "P" for j in atom.neighbors
            )
            atom.atom_type = "OS" if has_phosphorus_neighbor else "OA"
        elif element == "S":
            atom.atom_type = "SA" if atom.is_acceptor and not any(
                graph[j].is_hydrogen for j in atom.neighbors
            ) else "S"
        elif element in ("F", "CL", "BR", "I"):
            atom.atom_type = {"F": "F", "CL": "Cl", "BR": "Br", "I": "I"}[element]
        elif element == "P":
            atom.atom_type = "P"
        else:
            candidate = element.capitalize()
            if candidate in ("Mg", "Ca", "Mn", "Fe", "Zn", "Ni", "Cu", "Co"):
                candidate = candidate if candidate in _VINA_ATOM_TYPES else candidate
            atom.atom_type = candidate
        if atom.atom_type not in _VINA_ATOM_TYPES:
            unknown.add(atom.atom_type)
            if not keep_unknown:
                atom.atom_type = ""
    for utype in sorted(unknown):
        warnings.append(
            f"atom type {utype!r} is not in the Vina type set; Vina may reject the "
            "receptor unless you remove those atoms (see keep_resnames/keep_hetero)"
        )
    return warnings


def is_n_acceptor(graph: Sequence[EngineAtom], index: int) -> bool:
    """Pragmatic nitrogen acceptor check (amide N is a non-acceptor)."""
    atom = graph[index]
    for j in atom.neighbors:
        neighbor = graph[j]
        if neighbor.element.upper() != "C":
            continue
        for k in neighbor.neighbors:
            if k == index:
                continue
            partner = graph[k]
            if partner.element.upper() in ("O", "N", "S") and not partner.is_hydrogen:
                return False  # amide-like C(=X)-N
    return True


def merge_nonpolar_hydrogens(graph: list[EngineAtom]) -> tuple[list[EngineAtom], int]:
    """Merge non-polar hydrogens into their heavy atoms (``-U nphs``).

    Charges of the merged hydrogens are added onto the heavy atom, exactly
    like ``prepare_receptor4.py`` does before writing the PDBQT.
    """
    merged_count = 0
    keep: list[bool] = [True] * len(graph)
    for i, atom in enumerate(graph):
        if not atom.is_hydrogen:
            continue
        heavy = next((j for j in atom.neighbors if not graph[j].is_hydrogen), None)
        if heavy is None:
            continue
        if graph[heavy].element.upper() in ("C",):
            graph[heavy].charge += atom.charge
            keep[i] = False
            merged_count += 1
    result: list[EngineAtom] = []
    remap: dict[int, int] = {}
    for i, atom in enumerate(graph):
        if keep[i]:
            remap[i] = len(result)
            result.append(atom)
    for atom in result:
        atom.neighbors = [remap[j] for j in atom.neighbors if j in remap]
    return result, merged_count


# ---------------------------------------------------------------------------
# Receptor preparation
# ---------------------------------------------------------------------------
def removed_species_warning(removed: dict[str, int]) -> str | None:
    """WARNING text for removed metals/cofactors, or None if only waters were removed."""
    notable = sorted(
        (name, count) for name, count in removed.items()
        if name.upper() not in {"HOH", "DOD", "WAT", "H2O", "TIP", "TIP3", "SOL"}
    )
    if not notable:
        return None
    species = ", ".join(f"{count} {name}" for name, count in notable)
    return (
        f"removed hetero species: {species}. These may be structurally or "
        "catalytically required (metals/cofactors); keep them with "
        "receptor.keep_resnames or keep_hetero if they belong in the binding site"
    )


class ReceptorPreparator:
    """Prepare a receptor PDB -> PDBQT following prepare_receptor4.py logic."""

    def __init__(self, options: ReceptorPrepOptions | None = None) -> None:
        self.options = options or ReceptorPrepOptions()

    # -- public API ---------------------------------------------------------
    def prepare(
        self,
        input_path: str | Path,
        output_dir: str | Path,
        basename: str = "receptor",
    ) -> ReceptorPrepResult:
        """Run the full preparation pipeline.

        Steps: parse -> filter (chains/altloc/waters/hetero) -> add
        hydrogens & Gasteiger charges -> merge non-polar hydrogens ->
        AD4 typing -> write PDBQT + clean PDB.
        """
        options = self.options
        input_path = Path(input_path)
        if not input_path.is_file():
            raise PreparationError(f"receptor input not found: {input_path}")
        output_dir = ensure_dir(output_dir)
        result = ReceptorPrepResult(engine="none")

        atoms = parse_pdb(input_path)
        result.atoms_in = len(atoms)
        if not atoms:
            raise PreparationError(f"no atoms parsed from {input_path}")

        atoms, stats = self._filter_atoms(atoms)
        result.waters_removed = stats["waters"]
        result.hetero_removed = stats["hetero"]
        result.warnings.extend(stats["warnings"])
        result.removed_resnames = dict(stats["removed_resnames"])

        # Structure quality checks (work order items 6, 7, 10) - run on the
        # filtered atom set (the actual docking receptor) so kept metals are
        # the ones whose coordination is recorded.
        from .structure_qc import (
            detect_disulfides,
            detect_metal_coordination,
            detect_missing_residues,
            detect_reduced_cysteines,
            missing_residue_warning,
            reduced_cys_warning,
        )
        result.missing_residues = detect_missing_residues(atoms)
        result.disulfides = detect_disulfides(atoms)
        result.reduced_cysteines = detect_reduced_cysteines(
            atoms, result.disulfides)
        result.metal_coordination = detect_metal_coordination(atoms)
        for message in (missing_residue_warning(result.missing_residues),
                        reduced_cys_warning(result.reduced_cysteines)):
            if message:
                result.warnings.append(message)
                logger.warning("%s", message)

        engines = _engine_chain(options.engine) if options.add_hydrogens \
            else [PassThroughEngine()]
        graph: list[EngineAtom] | None = None
        engine: BaseEngine | None = None
        last_error: Exception | None = None
        for position, candidate in enumerate(engines):
            try:
                graph = candidate.process(atoms, options.charge_model)
            except Exception as exc:  # noqa: BLE001 - optional deps must never crash
                # Broken optional dependency (AttributeError/ImportError inside
                # the bindings, parse failure, ...): fall through to the next
                # engine and record both the failure and the fallback.
                last_error = exc
                fallback = engines[position + 1].name if position + 1 < len(engines) else "none"
                message = (
                    f"preparation engine '{candidate.name}' failed "
                    f"({type(exc).__name__}: {exc}); falling back to '{fallback}'"
                )
                logger.warning("%s", message)
                result.warnings.append(message)
                continue
            engine = candidate
            break
        if engine is None or graph is None:
            raise PreparationError(
                f"receptor preparation failed on every engine "
                f"(last error: {type(last_error).__name__}: {last_error})"
            )
        result.engine = engine.name if options.add_hydrogens else "none (hydrogens disabled)"
        if engine.name == "none" and options.add_hydrogens:
            result.warnings.append(
                "Receptor prepared with the `none` engine: charges are 0.0 "
                "and hydrogens are absent. Results are geometry-only and "
                "not comparable to a charged, protonated receptor "
                "(install openbabel-wheel or rdkit)."
            )

        hydrogens_before = sum(1 for a in atoms if a.element.upper() == "H")
        hydrogens_after = sum(1 for g in graph if g.is_hydrogen)
        result.hydrogens_added = max(0, hydrogens_after - hydrogens_before)

        if options.merge_nonpolar_h:
            graph, merged = merge_nonpolar_hydrogens(graph)
            logger.debug("merged %d non-polar hydrogens", merged)

        type_warnings = assign_ad4_types(graph, keep_unknown=options.keep_unknown_types)
        result.warnings.extend(type_warnings)

        if not options.keep_unknown_types:
            graph = [g for g in graph if g.atom_type]

        atoms_out = self._to_atoms(graph)
        result.atoms_out = len(atoms_out)
        if not atoms_out:
            raise PreparationError("receptor is empty after preparation")

        # Disulfide preservation check (peer item 7b): every detected S-S
        # bridge must survive preparation with both SG atoms intact.
        from .structure_qc import broken_disulfide_warning
        broken = broken_disulfide_warning(atoms_out, result.disulfides)
        if broken:
            result.warnings.append(broken)
            logger.warning("%s", broken)

        result.pdbqt_path = output_dir / f"{basename}.pdbqt"
        write_pdbqt(
            atoms_out,
            result.pdbqt_path,
            remarks=[f"prepared by DockFlow-Automator (engine={result.engine})"],
        )
        result.pdb_path = output_dir / f"{basename}_clean.pdb"
        write_pdb(atoms_out, result.pdb_path)
        logger.info(
            "receptor prepared: %d -> %d atoms (engine=%s, +%dH) -> %s",
            result.atoms_in,
            result.atoms_out,
            result.engine,
            result.hydrogens_added,
            result.pdbqt_path.name,
        )
        return result

    # -- steps --------------------------------------------------------------
    def _filter_atoms(self, atoms: Sequence[Atom]) -> tuple[list[Atom], dict[str, Any]]:
        options = self.options
        warnings: list[str] = []
        removed_resnames: dict[str, int] = {}
        stats = {"waters": 0, "hetero": 0, "warnings": warnings,
                 "removed_resnames": removed_resnames}
        selected: list[Atom] = []
        chain_set = {c.strip() for c in options.chains} if options.chains else None
        remove_res = {(c, int(r)) for c, r in options.remove_residues}
        remove_names = {n.strip().upper() for n in options.remove_resnames}
        keep_res = {n.strip().upper() for n in options.keep_resnames}

        atoms = self._resolve_altlocs(atoms, options.altloc)

        def _drop(atom: Atom, bucket: str) -> None:
            stats[bucket] += 1
            resname = atom.resname.strip().upper()
            removed_resnames[resname] = removed_resnames.get(resname, 0) + 1

        for atom in atoms:
            resname = atom.resname.strip().upper()
            if chain_set is not None and atom.chain not in chain_set:
                continue
            if atom.is_water:
                if options.keep_water:
                    selected.append(atom)
                else:
                    _drop(atom, "waters")
                continue
            if (atom.chain, atom.resseq) in remove_res or resname in remove_names:
                _drop(atom, "hetero")
                continue
            if not atom.is_polymer and not options.keep_hetero and resname not in keep_res:
                _drop(atom, "hetero")
                continue
            selected.append(atom)

        if stats["waters"]:
            logger.debug("removed %d water molecules", stats["waters"])
        if stats["hetero"]:
            logger.debug("removed %d hetero atoms (ligands/cofactors)", stats["hetero"])
        removal_warning = removed_species_warning(removed_resnames)
        if removal_warning is not None:
            # Only warn about non-water species: waters are a routine removal,
            # metals/cofactors are a scientific decision the user must see.
            logger.warning("%s", removal_warning)
            warnings.append(removal_warning)
        if not selected:
            warnings.append("all atoms were filtered out; check chains/remove options")
        return selected, stats

    @staticmethod
    def _resolve_altlocs(atoms: Sequence[Atom], policy: str) -> list[Atom]:
        """Keep one atom per altLoc group (``-U altloc``).

        ``best`` keeps the alternative with the highest average occupancy
        (blank altLoc wins over any letter), an explicit letter keeps that
        alternative, and ``""`` keeps everything.
        """
        if policy == "":
            return list(atoms)
        groups: dict[tuple, list[Atom]] = {}
        order: list[tuple] = []
        for atom in atoms:
            key = (atom.chain, atom.resname, atom.resseq, atom.icode, atom.name.strip())
            if key not in groups:
                groups[key] = []
                order.append(key)
            groups[key].append(atom)
        result: list[Atom] = []
        for key in order:
            candidates = groups[key]
            blanks = [a for a in candidates if not a.altloc]
            if len(candidates) == 1:
                result.append(candidates[0])
                continue
            if policy == "best":
                if blanks:
                    result.append(blanks[0])
                else:
                    result.append(max(candidates, key=lambda a: a.occupancy))
            else:
                wanted = [a for a in candidates if a.altloc == policy]
                if wanted:
                    result.append(wanted[0])
                elif blanks:
                    result.append(blanks[0])
                else:
                    result.append(candidates[0])
        return result

    @staticmethod
    def _to_atoms(graph: Sequence[EngineAtom]) -> list[Atom]:
        """Convert the typed engine graph into final Atom records."""
        atoms = [
            Atom(
                serial=i + 1,
                name=engine_atom.name,
                altloc="",
                resname=engine_atom.resname,
                chain=engine_atom.chain,
                resseq=engine_atom.resseq,
                icode=engine_atom.icode,
                x=engine_atom.x,
                y=engine_atom.y,
                z=engine_atom.z,
                occupancy=engine_atom.occupancy,
                bfactor=engine_atom.bfactor,
                element=engine_atom.element.capitalize(),
                charge=round(engine_atom.charge, 3),
                atom_type=engine_atom.atom_type,
                record_type=engine_atom.record_type,
            )
            for i, engine_atom in enumerate(graph)
        ]
        return renumber_serials(atoms)


# ---------------------------------------------------------------------------
# Ligand preparation (prepare_ligand4.py equivalent on Meeko + RDKit)
# ---------------------------------------------------------------------------
class LigandPreparator:
    """Prepare ligands (SMILES / SDF / MOL2 / PDB) into PDBQT via Meeko."""

    def __init__(self, options: LigandPrepOptions | None = None) -> None:
        self.options = options or LigandPrepOptions()

    # -- public API ---------------------------------------------------------
    def prepare(
        self,
        source: str | Path,
        output_dir: str | Path,
        identifier: str | None = None,
    ) -> LigandPrepResult:
        """Prepare one ligand from a file path or a SMILES string."""
        output_dir = ensure_dir(output_dir)
        name = identifier or _derive_identifier(source)
        result = LigandPrepResult(identifier=name, engine="meeko")
        try:
            mols = self._load_molecules(source, result)
            if not mols:
                raise PreparationError(f"no valid molecules in {source!r}")
            mol = mols[0]
            # Hydrogen provenance (audit item 26): was the input already
            # protonated, or did preparation add all hydrogens (making the
            # protonation state RDKit/Meeko's choice)?
            try:
                result.input_had_explicit_hydrogens = any(
                    a.GetAtomicNum() == 1 for a in mol.GetAtoms()
                )
            except Exception:  # noqa: BLE001 - provenance must never fail a run
                result.input_had_explicit_hydrogens = None
            # Salt / fragment provenance BEFORE stripping (peer item 14):
            # what was removed is recorded, and suspiciously large
            # counter-ions raise a warning.
            from .ligand_states import salt_report
            result.salt_stripping = salt_report(mol)
            if result.salt_stripping:
                warning = result.salt_stripping.get("large_fragment_warning")
                if warning:
                    result.warnings.append(warning)
            # State chemistry (peer items 11-13): protonation -> tautomers
            # -> stereoisomers, each level optional, cartesian product capped.
            states = self._state_variants(mol, result)
            result.num_protonation_variants = sum(
                1 for state in states if "prot" in state.tag.split("_")) or 1
            outputs: list[Path] = []
            for index, state in enumerate(states):
                suffix = "" if len(states) == 1 else f"_{state.tag}"
                variant_name = f"{name}{suffix}"
                pdbqt_path, sdf_path = self._prepare_variant(
                    state.mol, output_dir, variant_name, result)
                outputs.append(pdbqt_path)
                if sdf_path is not None and result.sdf_path is None:
                    result.sdf_path = sdf_path
                result.variants.append({
                    "name": variant_name,
                    "kind": state.kind,
                    "tag": state.tag,
                    "index": index + 1,
                    "smiles": state.smiles,
                    "num_heavy_atoms": result.num_heavy_atoms,
                    "num_rotatable_bonds": result.num_rotatable_bonds,
                    "pdbqt": str(pdbqt_path),
                })
            result.pdbqt_path = outputs[0]
            if len(outputs) > 1:
                result.warnings.append(
                    f"{len(outputs)} ligand states prepared "
                    f"({' '.join(p.stem for p in outputs)})"
                )
        except PreparationError as exc:
            result.error = str(exc)
            result.status = "error" if hasattr(result, "status") else None
            raise
        except Exception as exc:  # noqa: BLE001 - wrap third-party errors
            result.error = f"{type(exc).__name__}: {exc}"
            raise PreparationError(result.error) from exc
        return result

    def prepare_library(
        self,
        sdf_path: str | Path,
        output_dir: str | Path,
        max_records: int | None = None,
    ) -> list[LigandPrepResult]:
        """Prepare every record of a multi-record SDF (batch virtual screening)."""
        from rdkit import Chem

        output_dir = ensure_dir(output_dir)
        results: list[LigandPrepResult] = []
        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
        for index, mol in enumerate(supplier):
            if max_records is not None and index >= max_records:
                break
            if mol is None:
                results.append(
                    LigandPrepResult(identifier=f"record_{index + 1}",
                                     error="invalid SDF record")
                )
                continue
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"ligand_{index + 1}"
            name = (name or f"ligand_{index + 1}").strip().replace(" ", "_")[:60] or \
                f"ligand_{index + 1}"
            result = LigandPrepResult(identifier=name, engine="meeko")
            try:
                pdbqt_path, sdf_path = self._prepare_variant(mol, output_dir, name, result)
                result.pdbqt_path = pdbqt_path
                result.sdf_path = sdf_path
            except PreparationError as exc:
                result.error = str(exc)
            results.append(result)
        return results

    # -- internals ------------------------------------------------------------
    def _load_molecules(self, source: str | Path, result: LigandPrepResult) -> list[Any]:
        """Return RDKit molecules for a SMILES string or a structure file."""
        from rdkit import Chem

        path = Path(source)
        if isinstance(source, str) and not path.exists():
            text = source.strip()
            if not text:
                raise PreparationError("empty ligand source")
            mol = Chem.MolFromSmiles(text)
            if mol is None:
                raise PreparationError(f"invalid SMILES: {text!r}")
            result.smiles = text
            return [mol]
        if not path.is_file():
            raise PreparationError(f"ligand file not found: {path}")
        suffix = path.suffix.lower()
        if suffix in (".sdf", ".sd"):
            mols = [m for m in Chem.SDMolSupplier(str(path), removeHs=False,
                                                 sanitize=True) if m is not None]
            if mols and mols[0].HasProp("_Name"):
                result.identifier = result.identifier or mols[0].GetProp("_Name")
        elif suffix == ".mol2":
            mol = Chem.MolFromMol2File(str(path), removeHs=False, sanitize=True)
            mols = [mol] if mol is not None else []
        elif suffix in (".pdb", ".ent"):
            mol = Chem.MolFromPDBFile(str(path), removeHs=False, sanitize=True,
                                      proximityBonding=True)
            mols = [mol] if mol is not None else []
        elif suffix == ".smi":
            text = path.read_text(encoding="utf-8").splitlines()
            line = next((ln for ln in text if ln.strip()), "")
            smiles = line.split()[0] if line.split() else ""
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise PreparationError(f"invalid SMILES in {path}")
            result.smiles = smiles
            mols = [mol]
        else:
            raise PreparationError(f"unsupported ligand format: {suffix}")
        if not mols:
            raise PreparationError(f"no valid molecules parsed from {path}")
        return mols

    def _state_variants(self, mol: Any, result: LigandPrepResult) -> list[Any]:
        """Protonation -> tautomer -> stereo chain (peer items 11-13).

        Each level is optional; the cartesian product is capped by the
        per-level maxima.  The chain order matches the chemistry: choose
        the protonation state first, then tautomer, then stereoisomer.
        All state reports are written onto ``result`` for the manifest.
        """
        from .ligand_states import (
            StateVariant,
            protonation_states,
            stereocenter_report,
            stereoisomer_states,
            tautomer_states,
        )
        options = self.options

        # -- level 1: protonation --------------------------------------
        mode = None
        if options.protonate:
            mode = "enumerate" if options.protonate == "enumerate" else "major"
        if mode:
            mols, report = protonation_states(
                mol, mode=mode, ph=options.ph)
            result.protonation = report
            result.protonation_source = report["source"]
            if report.get("note"):
                result.warnings.append(f"protonation: {report['note']}")
            level1 = [
                StateVariant(mol=m, tag=f"prot{i}" if len(mols) > 1 else "",
                             kind="protonation", index=i,
                             smiles=report["states"][i - 1]["smiles"],
                             detail={"net_charge":
                                     report["states"][i - 1].get("net_charge")})
                for i, m in enumerate(mols, start=1)
            ]
        else:
            result.protonation = {
                "mode": "off",
                "source": "input protonation state (no pH adjustment)",
            }
            result.protonation_source = "input protonation state (no pH adjustment)"
            level1 = [StateVariant(mol=mol, kind="input", index=1,
                                   smiles=result.smiles or _quick_smiles(mol))]

        # -- level 2: tautomers -----------------------------------------
        if options.enumerate_tautomers:
            expanded: list = []
            for state in level1:
                mols, report = tautomer_states(
                    state.mol, max_tautomers=options.max_tautomers)
                if result.tautomers is None:
                    result.tautomers = report
                for j, t in enumerate(mols, start=1):
                    tag = "_".join(part for part in (state.tag, f"taut{j}")
                                    if part) if len(mols) > 1 else state.tag
                    expanded.append(StateVariant(
                        mol=t, tag=tag,
                        kind="tautomer" if len(mols) > 1 else state.kind,
                        index=j, smiles=report["tautomers"][j - 1]["smiles"],
                        detail={"tautomer": j}))
            level1 = expanded
        else:
            result.tautomers = {
                "source": "single input tautomer (no enumeration)",
                "tautomers": [],
            }

        # -- level 3: stereoisomers --------------------------------------
        report = stereocenter_report(mol)
        result.stereocenters = report
        if report["undefined"] and not options.enumerate_stereoisomers:
            count = report["possible_stereoisomers"]
            result.warnings.append(
                f"Ligand {result.identifier} has "
                f"{report['undefined']} undefined stereocenters; docking "
                f"will use the deposited geometry but the result is one of "
                f"{count} possible stereoisomers. Consider "
                f"enumerate_stereoisomers: true."
            )
        if options.enumerate_stereoisomers and report["undefined"]:
            expanded = []
            for state in level1:
                mols, stereo_report = stereoisomer_states(
                    state.mol, max_isomers=options.max_stereoisomers)
                enumerated = stereo_report.get("enumerated") or []
                if not enumerated:
                    continue
                if result.stereocenters is not None:
                    result.stereocenters = {
                        **result.stereocenters, "enumeration": stereo_report}
                for k, iso in enumerate(mols, start=1):
                    tag = "_".join(part for part in (state.tag, f"stereo{k}")
                                    if part) if len(mols) > 1 else state.tag
                    expanded.append(StateVariant(
                        mol=iso, tag=tag, kind="stereoisomer", index=k,
                        smiles=enumerated[k - 1]["smiles"],
                        detail={"stereoisomer": k}))
            level1 = expanded
        return level1

    def _prepare_variant(
        self,
        mol: Any,
        output_dir: Path,
        name: str,
        result: LigandPrepResult,
    ) -> tuple[Path, Path | None]:
        """Embed, minimise and write the PDBQT (+SDF) for one molecule."""
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors

        options = self.options
        mol = self._standardize(mol, result)
        if mol.GetNumConformers() == 0 and options.embed_3d:
            mol = self._embed_3d(mol, result)
        elif mol.GetNumConformers() == 0:
            raise PreparationError(
                f"ligand {name!r} has no 3D coordinates and embed_3d is disabled"
            )
        if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
            mol = Chem.AddHs(mol, addCoords=True)
        if options.minimize:
            self._minimize(mol, result)

        result.num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        result.num_heavy_atoms = mol.GetNumHeavyAtoms()
        result.num_atoms = mol.GetNumAtoms()
        if not result.smiles:
            result.smiles = Chem.MolToSmiles(mol)

        sdf_path = output_dir / f"{name}.sdf"
        writer = Chem.SDWriter(str(sdf_path))
        writer.SetKekulize(True)
        writer.write(mol)
        writer.close()

        pdbqt_text = self._meeko_pdbqt(mol, result)
        pdbqt_path = output_dir / f"{name}.pdbqt"
        pdbqt_path.write_text(pdbqt_text, encoding="utf-8")
        return pdbqt_path, sdf_path

    def _standardize(self, mol: Any, result: LigandPrepResult) -> Any:
        """Sanitize, deduplicate, keep the largest fragment, add hydrogens."""
        from rdkit import Chem

        options = self.options
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:  # noqa: BLE001
            raise PreparationError(f"RDKit sanitisation failed: {exc}") from exc
        if options.remove_salts or options.keep_largest_fragment:
            fragments = list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False))
            if len(fragments) > 1:
                # Salt stripping provenance was captured before this call
                # (prepare() records result.salt_stripping); the warning
                # text for the removed count stays actionable.
                removed = len(fragments) - 1
                if not any("salt/solvent fragment" in w for w in result.warnings):
                    result.warnings.append(
                        f"removed {removed} salt/solvent fragment(s); see "
                        f"manifest ligands[*].prep.salt_stripping for what "
                        f"was kept and dropped"
                    )
                mol = max(fragments, key=lambda m: m.GetNumAtoms())
                try:
                    Chem.SanitizeMol(mol)
                except Exception:  # noqa: BLE001
                    pass
        mol = Chem.AddHs(mol)
        return mol

    def _embed_3d(self, mol: Any, result: LigandPrepResult) -> Any:
        """ETKDGv3 embedding with MMFF/UFF minimisation to pick one conformer."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        params = AllChem.ETKDGv3()
        params.randomSeed = self.options.random_seed
        params.useRandomCoords = False
        conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=4, params=params))
        if not conf_ids:
            params.useRandomCoords = True
            conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=4, params=params))
        if not conf_ids:
            raise PreparationError("ETKDG 3D embedding failed for this ligand")
        energies = self._minimize(mol, result, return_energies=True)
        best = conf_ids[int(min(range(len(conf_ids)),
                                key=lambda i: energies[i]))] if energies else conf_ids[0]
        # Keep only the best conformer.
        kept = Chem.Mol(mol)
        kept.RemoveAllConformers()
        kept.AddConformer(mol.GetConformer(best), assignId=True)
        return kept

    def _minimize(self, mol: Any, result: LigandPrepResult,
                  return_energies: bool = False) -> list[float] | None:
        """MMFF94 (UFF fallback) optimisation of every conformer."""
        from rdkit.Chem import AllChem

        try:
            use_mmff = AllChem.MMFFHasAllMoleculeParams(mol)
        except Exception:  # noqa: BLE001
            use_mmff = False
        if use_mmff:
            results = AllChem.MMFFOptimizeMoleculeConfs(
                mol, maxIters=self.options.minimize_steps, mmffVariant="MMFF94"
            )
        else:
            if return_energies:
                result.warnings.append("MMFF parameters unavailable; using UFF")
            results = AllChem.UFFOptimizeMoleculeConfs(
                mol, maxIters=self.options.minimize_steps
            )
        energies = [float(entry[1]) if entry else 0.0 for entry in (results or [])]
        if return_energies:
            return energies
        return None

    @staticmethod
    def _meeko_pdbqt(mol: Any, result: LigandPrepResult) -> str:
        """Run Meeko's MoleculePreparation + PDBQTWriterLegacy.

        This is the modern replacement for ``prepare_ligand4.py``: Gasteiger
        charges, non-polar hydrogen merging and torsion tree construction.
        """
        try:
            from meeko import MoleculePreparation, PDBQTWriterLegacy
        except ImportError as exc:
            raise PreparationError(
                "meeko is required for ligand preparation "
                "(pip install 'dockflow-automator[prep]')"
            ) from exc
        try:
            preparator = MoleculePreparation()
        except TypeError:
            preparator = MoleculePreparation  # very old API fallback
        try:
            setups = preparator(mol)
        except TypeError:
            preparator = MoleculePreparation()
            setups = preparator(mol)
        if not setups:
            raise PreparationError("Meeko produced no ligand setup")
        pdbqt_text, is_ok, error_msg = PDBQTWriterLegacy.write_string(setups[0])
        if not is_ok:
            raise PreparationError(f"Meeko failed to write PDBQT: {error_msg}")
        return pdbqt_text


def _quick_smiles(mol: Any) -> str:
    """Best-effort SMILES for provenance records (empty on failure)."""
    try:
        from rdkit import Chem

        return Chem.MolToSmiles(mol)
    except Exception:  # noqa: BLE001
        return ""


def _derive_identifier(source: str | Path) -> str:
    """Build a filesystem-safe ligand identifier from a SMILES or path."""
    text = str(source)
    path = Path(text)
    if path.exists() or (len(text) < 4 and text.isalnum()):
        return path.stem or "ligand"
    import hashlib

    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()[:10]
    return f"ligand_{digest}"



