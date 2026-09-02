"""Grid box computation and Vina config tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dockflow_core.gridbox import (
    GridBox,
    box_from_coordinates,
    box_from_pocket,
    box_from_reference_ligand,
    box_from_residues,
    box_from_structure,
    box_from_vina_config,
    box_to_vina_config,
    box_union,
    distance_to_box,
    read_vina_config,
)
from dockflow_core.pdbio import parse_pdb


def test_box_from_coordinates():
    box = box_from_coordinates([[0, 0, 0], [2, 2, 2]], padding=4.0)
    assert box.center == (1.0, 1.0, 1.0)
    assert box.size == (10.0, 10.0, 10.0)
    assert box.volume == 1000.0


def test_box_minimum_size():
    box = box_from_coordinates([[1, 1, 1]], padding=1.0, min_size=12.0)
    assert all(s == 12.0 for s in box.size)


def test_box_from_coordinates_requires_points():
    with pytest.raises(ValueError):
        box_from_coordinates([])


def test_box_contains_and_distance():
    box = GridBox(center=(0, 0, 0), size=(10, 10, 10))
    assert box.contains((0, 0, 0))
    assert box.contains((4.9, 4.9, 4.9))
    assert not box.contains((5.1, 0, 0))
    assert distance_to_box(box, (0, 0, 0)) == 0.0
    assert abs(distance_to_box(box, (10, 0, 0)) - 5.0) < 1e-9


def test_box_from_residues(receptor_pdb_path: Path):
    box = box_from_residues(receptor_pdb_path, chain="A", residues=[1, 2], padding=2.0)
    # residues 1-2 span x in [0, 4.02]
    assert -2.0 <= box.min_corner[0] <= 0.0
    assert box.max_corner[0] >= 4.02
    assert box.source == "residues"


def test_box_from_residues_atoms_input(receptor_pdb_path: Path):
    atoms = parse_pdb(receptor_pdb_path)
    box = box_from_residues(atoms, chain=None, residues=[99], padding=3.0)
    assert box.center == (9.0, 9.0, 9.0)


def test_box_from_pocket(receptor_pdb_path: Path):
    box = box_from_pocket(receptor_pdb_path, "BEN", padding=4.0)
    # benzene carbons span roughly x 11.2-13.7
    assert box.min_corner[0] <= 11.0 and box.max_corner[0] >= 14.0
    assert box.source == "pocket:BEN"
    with pytest.raises(ValueError):
        box_from_pocket(receptor_pdb_path, "XXX")


def test_box_from_reference_ligand_file(ligand_pdbqt_path: Path):
    box = box_from_reference_ligand(ligand_pdbqt_path, padding=2.0)
    assert box.min_corner[0] <= 11.1 - 2.0


def test_box_from_structure(receptor_pdb_path: Path):
    box = box_from_structure(receptor_pdb_path, padding=1.0)
    assert box.min_corner[0] <= 0.0 and box.max_corner[0] >= 18.0


def test_vina_config_roundtrip(tmp_path: Path):
    box = GridBox(center=(1.5, -2.5, 3.25), size=(22, 24, 20.5))
    path = box_to_vina_config(box, tmp_path / "box.txt",
                              extra={"exhaustiveness": 16, "cpu": 4})
    data = read_vina_config(path)
    assert data["center_x"] == "1.500"
    assert data["size_z"] == "20.500"
    assert data["exhaustiveness"] == "16"
    again = box_from_vina_config(path)
    assert again.center == box.center
    assert again.size == box.size


def test_vina_config_missing_keys(tmp_path: Path):
    path = tmp_path / "bad.txt"
    path.write_text("center_x = 0\n", encoding="utf-8")
    with pytest.raises(ValueError):
        box_from_vina_config(path)


def test_box_union():
    a = GridBox(center=(0, 0, 0), size=(10, 10, 10))
    b = GridBox(center=(20, 0, 0), size=(10, 10, 10))
    union = box_union([a, b])
    assert union.center == (10.0, 0.0, 0.0)
    assert union.size == (30.0, 10.0, 10.0)


def test_box_serialization():
    box = GridBox(center=(1, 2, 3), size=(4, 5, 6), source="test")
    data = box.to_dict()
    again = GridBox.from_dict(data)
    assert again.center == box.center and again.size == box.size
    assert again.source == "test"


def test_corner_points():
    box = GridBox(center=(0, 0, 0), size=(2, 2, 2))
    corners = box.corner_points()
    assert corners.shape == (8, 3)
    assert set(np.unique(corners[:, 0])) == {-1.0, 1.0}
