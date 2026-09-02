"""Grid box (search space) definition and computation.

A grid box is fully described by a center and a size (both in Angstrom).
DockFlow can derive it from:

* the bounding box of a reference (co-crystallized) ligand plus padding,
* a set of active-site residues plus padding,
* explicit center/size values (e.g. from literature or a previous run).

Boxes can be exported to / imported from AutoDock Vina config files so the
exact search space is reproducible on the command line::

    center_x = 12.3
    center_y = -4.5
    center_z = 21.7
    size_x = 22
    size_y = 24
    size_z = 20
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .pdbio import Atom, parse_pdb, parse_pdbqt
from .utils import get_logger

logger = get_logger("gridbox")


@dataclass
class GridBox:
    """Rectangular search space: center + size (Angstrom)."""

    center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    size: tuple[float, float, float] = (20.0, 20.0, 20.0)
    source: str = "explicit"
    padding: float = 0.0

    # -- geometry -----------------------------------------------------------
    @property
    def min_corner(self) -> np.ndarray:
        return np.asarray(self.center) - np.asarray(self.size) / 2.0

    @property
    def max_corner(self) -> np.ndarray:
        return np.asarray(self.center) + np.asarray(self.size) / 2.0

    @property
    def volume(self) -> float:
        sx, sy, sz = self.size
        return float(sx) * float(sy) * float(sz)

    def expanded(self, margin: float) -> GridBox:
        size = tuple(max(1.0, s + 2.0 * margin) for s in self.size)
        return GridBox(center=self.center, size=size, source=self.source, padding=self.padding)

    def contains(self, point: Sequence[float], slack: float = 0.0) -> bool:
        p = np.asarray(point, dtype=float)
        return bool(
            np.all(p >= self.min_corner - slack) and np.all(p <= self.max_corner + slack)
        )

    def corner_points(self) -> np.ndarray:
        """The 8 box corners as an (8, 3) array."""
        lo, hi = self.min_corner, self.max_corner
        corners = [
            (x, y, z)
            for x in (lo[0], hi[0])
            for y in (lo[1], hi[1])
            for z in (lo[2], hi[2])
        ]
        return np.asarray(corners, dtype=float)

    # -- serialization ------------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        return {
            "center": list(self.center),
            "size": list(self.size),
            "source": self.source,
            "padding": self.padding,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> GridBox:
        return cls(
            center=tuple(data.get("center", (0.0, 0.0, 0.0))),  # type: ignore[arg-type]
            size=tuple(data.get("size", (20.0, 20.0, 20.0))),  # type: ignore[arg-type]
            source=data.get("source", "explicit"),
            padding=float(data.get("padding", 0.0)),
        )

    def __str__(self) -> str:
        c = ", ".join(f"{v:8.3f}" for v in self.center)
        s = ", ".join(f"{v:8.3f}" for v in self.size)
        return f"GridBox(center=({c}), size=({s}), volume={self.volume:.0f} A^3, {self.source})"


# ---------------------------------------------------------------------------
# Box computation
# ---------------------------------------------------------------------------
def box_from_coordinates(
    coordinates: Sequence[Sequence[float]],
    padding: float = 4.0,
    source: str = "coordinates",
    min_size: float = 10.0,
) -> GridBox:
    """Bounding box of a point cloud plus per-axis padding."""
    arr = np.asarray(coordinates, dtype=float).reshape(-1, 3)
    if arr.size == 0:
        raise ValueError("box_from_coordinates: empty coordinate set")
    low = arr.min(axis=0)
    high = arr.max(axis=0)
    size = np.maximum(high - low + 2.0 * padding, min_size)
    center = (high + low) / 2.0
    return GridBox(
        center=tuple(float(v) for v in center),  # type: ignore[arg-type]
        size=tuple(float(v) for v in size),  # type: ignore[arg-type]
        source=source,
        padding=float(padding),
    )


def box_from_atoms(atoms: Sequence[Atom], padding: float = 4.0, source: str = "atoms") -> GridBox:
    coords = [(a.x, a.y, a.z) for a in atoms]
    return box_from_coordinates(coords, padding=padding, source=source)


def _load_atoms(path: str | Path) -> list[Atom]:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".pdbqt", ".qt"}:
        return parse_pdbqt(path).atoms
    return parse_pdb(path)


def box_from_reference_ligand(
    ligand_path: str | Path, padding: float = 4.0
) -> GridBox:
    """Grid box from a reference ligand file (PDB or PDBQT)."""
    atoms = _load_atoms(ligand_path)
    if not atoms:
        raise ValueError(f"no atoms found in reference ligand {ligand_path}")
    return box_from_atoms(atoms, padding=padding, source=f"ligand:{Path(ligand_path).name}")


def box_from_pocket(
    structure_path: str | Path,
    resname: str,
    chain: str | None = None,
    resseq: int | None = None,
    padding: float = 4.0,
) -> GridBox:
    """Grid box from a co-crystallized ligand found inside a structure.

    Args:
        structure_path: PDB (or PDBQT) file of the target.
        resname: 3-letter residue name of the ligand, e.g. ``MK1``.
        chain: optional chain restriction.
        resseq: optional residue-number restriction (disambiguation).
        padding: padding added around the ligand bounding box.
    """
    atoms = _load_atoms(structure_path)
    selected = [
        a
        for a in atoms
        if a.resname.strip().upper() == resname.strip().upper()
        and (chain is None or a.chain == chain)
        and (resseq is None or a.resseq == resseq)
        and not a.is_hydrogen
    ]
    if not selected:
        raise ValueError(
            f"ligand residue {resname!r} not found in {structure_path} "
            f"(chain={chain!r}, resseq={resseq!r})"
        )
    return box_from_atoms(selected, padding=padding, source=f"pocket:{resname}")


def box_from_residues(
    structure_path: str | Path | Sequence[Atom],
    chain: str | None,
    residues: Sequence[int],
    padding: float = 5.0,
) -> GridBox:
    """Grid box enclosing all atoms of the given active-site residues.

    Args:
        structure_path: PDB/PDBQT path or an already-parsed atom list.
        chain: chain identifier (None = any chain).
        residues: residue sequence numbers to include.
        padding: padding around the selected atoms.
    """
    if isinstance(structure_path, (str, Path)):
        atoms = _load_atoms(structure_path)
    else:
        atoms = list(structure_path)
    wanted = {int(r) for r in residues}
    selected = [
        a
        for a in atoms
        if a.resseq in wanted and (chain is None or a.chain == chain) and not a.is_hydrogen
    ]
    if not selected:
        raise ValueError(f"no atoms for residues {sorted(wanted)} chain={chain!r}")
    return box_from_atoms(selected, padding=padding, source="residues")


def box_from_structure(
    structure_path: str | Path, padding: float = 4.0, atoms_limit: int = 20000
) -> GridBox:
    """Grid box around the whole structure (last-resort default)."""
    atoms = [a for a in _load_atoms(structure_path) if not a.is_hydrogen][:atoms_limit]
    if not atoms:
        raise ValueError(f"no atoms in {structure_path}")
    return box_from_atoms(atoms, padding=padding, source="structure")


def box_union(boxes: Sequence[GridBox]) -> GridBox:
    """Smallest box containing all input boxes."""
    if not boxes:
        raise ValueError("box_union: no boxes")
    lows = np.min([b.min_corner for b in boxes], axis=0)
    highs = np.max([b.max_corner for b in boxes], axis=0)
    center = (lows + highs) / 2.0
    size = highs - lows
    return GridBox(center=tuple(center), size=tuple(size), source="union")  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Vina config files
# ---------------------------------------------------------------------------
_BOX_KEYS = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
_NUMBER_RE = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def box_to_vina_config(
    box: GridBox,
    path: str | Path,
    extra: dict[str, Any] | None = None,
    receptor: str | Path | None = None,
    ligand: str | Path | None = None,
    out: str | Path | None = None,
) -> Path:
    """Write a Vina ``--config`` file describing this box.

    Args:
        box: the search space.
        path: destination file.
        extra: additional ``key = value`` pairs (exhaustiveness, cpu, ...).
        receptor/ligand/out: optional file references written into the config.
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = ["# DockFlow-Automator grid box config"]
    cx, cy, cz = (float(v) for v in box.center)
    sx, sy, sz = (float(v) for v in box.size)
    lines += [
        f"center_x = {cx:.3f}",
        f"center_y = {cy:.3f}",
        f"center_z = {cz:.3f}",
        f"size_x = {sx:.3f}",
        f"size_y = {sy:.3f}",
        f"size_z = {sz:.3f}",
    ]
    if receptor:
        lines.append(f"receptor = {_quote(str(receptor))}")
    if ligand:
        lines.append(f"ligand = {_quote(str(ligand))}")
    if out:
        lines.append(f"out = {_quote(str(out))}")
    for key, value in (extra or {}).items():
        if value is not None:
            lines.append(f"{key} = {value}")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def _quote(value: str) -> str:
    return f'"{value}"' if (" " in value or "=" in value) else value


def read_vina_config(path: str | Path) -> dict[str, str]:
    """Read a Vina config file into a plain ``key -> value`` dict."""
    result: dict[str, str] = {}
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        result[key.strip()] = value.strip().strip('"').strip("'")
    return result


def box_from_vina_config(path: str | Path) -> GridBox:
    """Parse a Vina config file into a :class:`GridBox`."""
    data = read_vina_config(path)
    missing = [k for k in _BOX_KEYS if k not in data]
    if missing:
        raise ValueError(f"vina config {path} missing keys: {', '.join(missing)}")
    center = tuple(float(data[f"center_{a}"]) for a in "xyz")  # type: ignore[arg-type]
    size = tuple(float(data[f"size_{a}"]) for a in "xyz")  # type: ignore[arg-type]
    return GridBox(center=center, size=size, source=f"config:{Path(path).name}")


def distance_to_box(box: GridBox, point: Sequence[float]) -> float:
    """Distance from a point to the box surface (0 if inside)."""
    p = np.asarray(point, dtype=float)
    lo, hi = box.min_corner, box.max_corner
    deltas = np.maximum(np.maximum(lo - p, p - hi), 0.0)
    return float(math.sqrt(float(np.sum(deltas**2))))
