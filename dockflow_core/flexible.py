"""Flexible receptor residues for Vina side-chain docking (peer item 8).

Vina supports two-file flexible docking: a rigid receptor PDBQT plus a
flex PDBQT holding selected side chains as torsion trees wrapped in
``BEGIN_RES``/``END_RES`` blocks.  This module generates both files from
a prepared receptor PDBQT and a residue selection.

Design decisions (ADR-0003):

* The rotatable side-chain bonds are taken from a per-residue-type chi
  table (chi1..chi4 as used by standard side-chain flexibility).  This
  is deterministic and matches what ``prepare_flexreceptor4.py`` did
  for the common cases, without a ring/aromaticity perception pass over
  the receptor graph.
* The flexible part of a residue is its side chain from CB (or SG for
  CYS) onward; the backbone (N, CA, C, O, OXT) always stays rigid.
* GLY/ALA/PRO side chains have no rotatable bonds - requesting them
  produces a warning, not an error.
* ``auto`` selection derives the residue list from the geometric
  contacts between the co-crystallized reference ligand and the
  receptor (residues with any contact become flexible) - a documented
  heuristic, overridable with an explicit list.

Score comparability: flexible-receptor docking scores are NOT
comparable with rigid-receptor docking scores of the same ligand (the
search space changed); every output that reports a flex-res run says so.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

from .pdbio import Atom, parse_pdbqt, write_pdbqt
from .utils import DockFlowError, get_logger

if TYPE_CHECKING:  # pragma: no cover - typing only
    pass

logger = get_logger("flexible")

# Rotatable side-chain bonds per residue type: (atom_a, atom_b) pairs,
# the standard chi angles used by side-chain flexibility.
CHI_BONDS: dict[str, list[tuple[str, str]]] = {
    "ARG": [("CB", "CG"), ("CG", "CD"), ("CD", "NE")],
    "ASN": [("CB", "CG")],
    "ASP": [("CB", "CG")],
    "CYS": [("CB", "SG")],
    "GLN": [("CB", "CG"), ("CG", "CD"), ("CD", "OE1")],
    "GLU": [("CB", "CG"), ("CG", "CD"), ("CD", "OE1")],
    "HIS": [("CB", "CG")],
    "ILE": [("CB", "CG1"), ("CB", "CG2")],
    "LEU": [("CB", "CG")],
    "LYS": [("CB", "CG"), ("CG", "CD"), ("CD", "CE")],
    "MET": [("CB", "CG"), ("CG", "SD"), ("SD", "CE")],
    "PHE": [("CB", "CG")],
    "SER": [("CB", "OG")],
    "THR": [("CB", "OG1"), ("CB", "CG2")],
    "TRP": [("CB", "CG")],
    "TYR": [("CB", "CG")],
    "VAL": [("CB", "CG1"), ("CB", "CG2")],
}

BACKBONE = {"N", "CA", "C", "O", "OXT", "H", "HA", "H2", "H3", "1H", "2H", "3H"}
# Side chains without rotatable bonds (GLY/ALA: nothing to flex; PRO: ring).
NO_FLEX = {"GLY", "ALA", "PRO"}


def parse_flexible_residues(
    selection: Sequence[dict] | str | None,
) -> list[dict[str, str | int]]:
    """Normalise the ``receptor.flexible_residues`` YAML value.

    Accepts a list of ``{chain, resseq}`` dicts or the string ``"auto"``
    (resolved later from reference-ligand contacts).
    """
    if selection is None:
        return []
    if isinstance(selection, str):
        if selection.strip().lower() == "auto":
            return [{"auto": True}]
        raise DockFlowError(
            f"receptor.flexible_residues must be a list of {{chain, resseq}} "
            f"or 'auto', got {selection!r}"
        )
    if not isinstance(selection, (list, tuple)):
        raise DockFlowError(
            "receptor.flexible_residues must be a list of {chain, resseq} "
            "mappings or 'auto'"
        )
    residues: list[dict[str, str | int]] = []
    for entry in selection:
        if not isinstance(entry, dict) or "resseq" not in entry:
            raise DockFlowError(
                f"flexible residue entries must be {{chain, resseq}} dicts, "
                f"got {entry!r}"
            )
        residues.append({
            "chain": str(entry.get("chain", "")).strip(),
            "resseq": int(entry["resseq"]),
        })
    return residues


def auto_flexible_residues(
    receptor_pdbqt: str | Path,
    reference_ligand_pdb: str | Path,
    contact_cutoff: float = 4.5,
) -> list[dict[str, str | int]]:
    """Residues contacting the reference ligand (heuristic 'auto')."""
    import numpy as np

    from .pdbio import parse_pdb

    receptor_atoms = [a for a in parse_pdbqt(receptor_pdbqt).atoms
                      if not a.is_hydrogen]
    ligand_atoms = [a for a in parse_pdb(reference_ligand_pdb)
                    if not a.is_hydrogen]
    if not ligand_atoms:
        return []
    ligand_xyz = np.array([[a.x, a.y, a.z] for a in ligand_atoms], dtype=float)
    residues: set[tuple[str, int]] = set()
    for atom in receptor_atoms:
        if (atom.resname or "").strip().upper() not in CHI_BONDS:
            continue  # only side chains with rotatable bonds
        distance = float(np.linalg.norm(
            ligand_xyz - np.array([atom.x, atom.y, atom.z]), axis=1).min())
        if distance <= contact_cutoff:
            residues.add((atom.chain.strip(), int(atom.resseq)))
    return [{"chain": chain, "resseq": resseq}
            for chain, resseq in sorted(residues)]


def _residue_key(atom: Atom) -> tuple[str, str, int]:
    return (atom.chain.strip(), (atom.resname or "").strip().upper(),
            int(atom.resseq))


def build_flexible_pdbqts(
    receptor_pdbqt: str | Path,
    residues: Sequence[dict],
    output_dir: str | Path,
    basename: str = "receptor",
) -> tuple[Path, Path, list[dict]]:
    """Split a prepared receptor into rigid + flexible PDBQT files.

    Returns ``(rigid_path, flex_path, report)`` where report lists the
    residues actually made flexible and any skipped ones (no rotatable
    side chain / not found) with reasons - never silently dropped.
    """
    receptor_pdbqt = Path(receptor_pdbqt)
    atoms = parse_pdbqt(receptor_pdbqt).atoms
    if not atoms:
        raise DockFlowError(f"no atoms parsed from {receptor_pdbqt}")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    by_residue: dict[tuple[str, str, int], list[Atom]] = {}
    for atom in atoms:
        by_residue.setdefault(_residue_key(atom), []).append(atom)

    flex_atoms: list[Atom] = []
    rigid_atoms: list[Atom] = []
    report: list[dict] = []

    wanted = {(str(r.get("chain", "")).strip(), int(r["resseq"]))
              for r in residues}
    resolved: set[tuple[str, int]] = set()
    for key, residue_atoms in by_residue.items():
        chain, resname, resseq = key
        if (chain, resseq) not in wanted:
            rigid_atoms.extend(residue_atoms)
            continue
        resolved.add((chain, resseq))
        if resname in NO_FLEX:
            report.append({"chain": chain, "resseq": resseq, "resname": resname,
                           "status": "skipped",
                           "reason": f"{resname} has no rotatable side-chain "
                                     f"bonds"})
            rigid_atoms.extend(residue_atoms)
            continue
        chi_bonds = CHI_BONDS.get(resname)
        if not chi_bonds:
            report.append({"chain": chain, "resseq": resseq, "resname": resname,
                           "status": "skipped",
                           "reason": "no chi-angle table for this residue type"})
            rigid_atoms.extend(residue_atoms)
            continue
        side, backbone = _split_side_chain(residue_atoms)
        if not side:
            report.append({"chain": chain, "resseq": resseq, "resname": resname,
                           "status": "skipped", "reason": "no side-chain atoms"})
            rigid_atoms.extend(residue_atoms)
            continue
        flex_atoms.extend(side)
        rigid_atoms.extend(backbone)
        report.append({"chain": chain, "resseq": resseq, "resname": resname,
                       "status": "flexible",
                       "rotatable_bonds": sum(
                           1 for a, b in chi_bonds
                           if a in {(x.name or "").strip().upper()
                                    for x in residue_atoms}
                           and b in {(x.name or "").strip().upper()
                                     for x in residue_atoms}),
                       "atoms_flexible": len(side)})

    for chain, resseq in sorted(wanted - resolved):
        report.append({"chain": chain, "resseq": resseq, "resname": "?",
                       "status": "skipped", "reason": "residue not found"})

    if not flex_atoms:
        raise DockFlowError(
            "no flexible side chains could be generated from "
            f"{receptor_pdbqt} (requested: {sorted(wanted)}); see the report "
            "entries for why each residue was skipped"
        )

    rigid_path = output_dir / f"{basename}_rigid.pdbqt"
    write_pdbqt(rigid_atoms, rigid_path,
                remarks=["rigid part (DockFlow flexible-receptor split)"])
    flex_path = output_dir / f"{basename}_flex.pdbqt"
    _write_flex_file(flex_atoms, path=flex_path)
    return rigid_path, flex_path, report


def _write_flex_file(atoms: Sequence[Atom], path: Path) -> None:
    """Write the flexible-residue PDBQT (BEGIN_RES/ROOT/BRANCH blocks).

    Per-residue emission lives in :func:`_write_flex_residue`; this
    function owns the shared serial counter (``emit``) and the output
    lines.  No TORSDOF is written - Vina's flex parser rejects the tag.
    """
    from .pdbio import format_pdbqt_line

    lines: list[str] = []
    state = {"serial": 0}

    def emit(atom: Atom) -> int:
        state["serial"] += 1
        lines.append(format_pdbqt_line(atom, state["serial"]))
        return state["serial"]

    def next_serial() -> int:
        return state["serial"] + 1

    by_residue: dict[tuple[str, str, int], list[Atom]] = {}
    for atom in atoms:
        by_residue.setdefault(_residue_key(atom), []).append(atom)
    for key in sorted(by_residue):
        chain, resname, resseq = key
        _write_flex_residue(by_residue[key], resname, chain, resseq,
                            emit=emit, next_serial=next_serial, lines=lines)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_flex_residue(residue_atoms: Sequence[Atom], resname: str,
                        chain: str, resseq: int, emit, next_serial,
                        lines: list[str]) -> int:
    """Emit BEGIN_RES..END_RES for one residue; return its branch count.

    The side chain is partitioned into rigid segments by removing the chi
    bonds from the residue graph; every chi bond becomes a nested BRANCH
    (chains like LYS CB-CG-CD-CE produce branch-in-branch trees exactly
    like ``prepare_flexreceptor4.py`` output).  BRANCH lines reference
    written atom serials; atom lines use :func:`format_pdbqt_line`.
    """
    from .pdbio import format_pdbqt_line  # noqa: F401 - used via emit()

    serials: dict[int, int] = {}
    names = {(a.name or "").strip().upper(): a for a in residue_atoms
             if not a.is_hydrogen}
    lines.append(f"BEGIN_RES {resname} {chain} {resseq}")

    adjacency = _residue_adjacency(residue_atoms)
    chi_pairs = [(names[a], names[b]) for a, b in CHI_BONDS.get(resname, [])
                 if a in names and b in names]
    chi_edges = {(id(a), id(b)) for a, b in chi_pairs} | \
                {(id(b), id(a)) for a, b in chi_pairs}

    def flood(start: Atom) -> list[Atom]:
        """Atoms reachable from ``start`` without crossing chi bonds."""
        seen = {id(start)}
        stack = [start]
        tree = [start]
        while stack:
            current = stack.pop()
            for neighbor in adjacency.get(id(current), []):
                if id(neighbor) in seen:
                    continue
                if (id(current), id(neighbor)) in chi_edges:
                    continue
                seen.add(id(neighbor))
                tree.append(neighbor)
                stack.append(neighbor)
        return tree

    root_name = "CB" if "CB" in names else next(iter(names))
    root_segment = flood(names[root_name])
    segment_of: dict[int, int] = {id(a): 0 for a in root_segment}
    segments: dict[int, list[Atom]] = {0: root_segment}
    bond_segments: list[tuple[Atom, Atom, int]] = []
    for index, (atom_a, atom_b) in enumerate(chi_pairs, start=1):
        segment = [a for a in flood(atom_b) if id(a) not in segment_of]
        heavy_ids = {id(a) for a in segment}
        segment += [a for a in residue_atoms if a.is_hydrogen and any(
            id(n) in heavy_ids for n in adjacency.get(id(a), [])
            if n is not a) and id(a) not in segment_of]
        segments[index] = segment
        for atom in segment:
            segment_of[id(atom)] = index
        bond_segments.append((atom_a, atom_b, index))

    lines.append("ROOT")
    for atom in root_segment:
        serials[id(atom)] = emit(atom)
    lines.append("ENDROOT")

    def emit_branches(parent_segment: int) -> int:
        count = 0
        for atom_a, atom_b, index in bond_segments:
            if segment_of.get(id(atom_a)) != parent_segment:
                continue
            segment = segments[index]
            if not segment or id(atom_b) not in segment_of:
                continue
            serial_a = serials.get(id(atom_a), 0)
            serial_b = next_serial()
            lines.append(_branch_line(serial_a, serial_b))
            for atom in segment:
                serials[id(atom)] = emit(atom)
            count += 1 + emit_branches(index)
            lines.append(_end_branch_line(serial_a, serial_b))
        return count

    branches = emit_branches(0)
    for atom in residue_atoms:  # safety net: anything not covered
        if id(atom) not in serials:
            serials[id(atom)] = emit(atom)
    lines.append(f"END_RES {resname} {chain} {resseq}")
    return branches


def _split_side_chain(residue_atoms: Sequence[Atom]
                      ) -> tuple[list[Atom], list[Atom]]:
    """Side-chain atoms (CB onward + their H) vs backbone atoms."""
    names = {(a.name or "").strip().upper() for a in residue_atoms}
    if "CB" not in names:
        return [], list(residue_atoms)
    adjacency = _residue_adjacency(residue_atoms)
    side_heavy = {"CB"}
    changed = True
    while changed:
        changed = False
        for atom in residue_atoms:
            name = (atom.name or "").strip().upper()
            if name in side_heavy or name in BACKBONE or atom.is_hydrogen:
                continue
            neighbors = adjacency.get(id(atom), [])
            if any((n.name or "").strip().upper() in side_heavy
                   for n in neighbors):
                side_heavy.add(name)
                changed = True
    side: list[Atom] = []
    backbone: list[Atom] = []
    for atom in residue_atoms:
        name = (atom.name or "").strip().upper()
        if atom.is_hydrogen:
            neighbors = adjacency.get(id(atom), [])
            parent = next((n for n in neighbors if not n.is_hydrogen), None)
            parent_name = ((parent.name or "").strip().upper()
                           if parent else "CA")
            (side if parent_name in side_heavy else backbone).append(atom)
        elif name in side_heavy:
            side.append(atom)
        else:
            backbone.append(atom)
    return side, backbone


def _residue_adjacency(residue_atoms: Sequence[Atom]) -> dict[int, list[Atom]]:
    """Covalent-radius bond graph within one residue."""
    import numpy as np

    covalent = {"C": 0.77, "N": 0.75, "O": 0.73, "S": 1.02, "H": 0.37}
    xyz = np.array([[a.x, a.y, a.z] for a in residue_atoms], dtype=float)
    adjacency: dict[int, list[Atom]] = {id(a): [] for a in residue_atoms}
    for i, a in enumerate(residue_atoms):
        for j in range(i + 1, len(residue_atoms)):
            b = residue_atoms[j]
            radius = (covalent.get((a.element or "C").upper(), 0.8)
                      + covalent.get((b.element or "C").upper(), 0.8))
            distance = float(np.linalg.norm(xyz[i] - xyz[j]))
            if distance <= radius + 0.45:
                adjacency[id(a)].append(b)
                adjacency[id(b)].append(a)
    return adjacency


def _attached_to(hydrogen: Atom, parent: Atom) -> bool:
    import numpy as np

    return float(np.linalg.norm(
        np.array([hydrogen.x, hydrogen.y, hydrogen.z])
        - np.array([parent.x, parent.y, parent.z]))) <= 1.3


def _branch_subtree(first: Atom, exclude: Atom,
                    residue_atoms: Sequence[Atom]) -> list[Atom]:
    """Atoms reachable from ``first`` without crossing back to ``exclude``."""
    adjacency = _residue_adjacency(residue_atoms)
    seen = {id(first)}
    frontier = [first]
    tree = [first]
    while frontier:
        current = frontier.pop()
        for neighbor in adjacency.get(id(current), []):
            if neighbor is exclude or id(neighbor) in seen:
                continue
            seen.add(id(neighbor))
            tree.append(neighbor)
            frontier.append(neighbor)
    # subtree heavy atoms + hydrogens attached to them (heavy parents stay out)
    heavy = {id(a) for a in tree}
    result = [a for a in residue_atoms
              if id(a) in heavy or (a.is_hydrogen and any(
                  id(n) in heavy for n in adjacency.get(id(a), [])
                  if n is not a))]
    return result


def _branch_line(serial_a: int, serial_b: int) -> str:
    return f"BRANCH {serial_a:5d} {serial_b:5d}"


def _end_branch_line(serial_a: int, serial_b: int) -> str:
    return f"ENDBRANCH {serial_a:5d} {serial_b:5d}"
