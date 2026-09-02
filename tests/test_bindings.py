"""C++ accelerator binding tests (skipped when the wheel is not built)."""

from __future__ import annotations

import numpy as np
import pytest

dfb = pytest.importorskip("dockflow_bindings",
                          reason="C++ bindings not built in this environment")


def test_version():
    assert dfb.__version__


def test_parse_pdbqt_atoms(ligand_pdbqt_text: str):
    atoms = dfb.parse_pdbqt_atoms(ligand_pdbqt_text)
    assert len(atoms) == 5
    assert atoms[0]["name"] == "C1"
    assert atoms[1]["atom_type"] == "OA"
    assert abs(atoms[1]["charge"] - (-0.641)) < 1e-9
    assert atoms[3]["resname"] == "LIG"


def test_parse_pdbqt_matches_python_parser(docked_pdbqt_text: str):
    from dockflow_core.pdbio import parse_pdbqt

    cpp_atoms = dfb.parse_pdbqt_atoms(docked_pdbqt_text)
    py_atoms = [a for model in parse_pdbqt(docked_pdbqt_text).models
                for a in model.atoms]
    assert len(cpp_atoms) == len(py_atoms)
    for cpp, py in zip(cpp_atoms, py_atoms, strict=True):
        assert cpp["name"] == py.name.strip()
        assert cpp["resname"] == py.resname
        assert abs(cpp["x"] - py.x) < 1e-6
        assert abs(cpp["charge"] - py.charge) < 1e-6
        assert cpp["atom_type"] == py.atom_type


def test_grid_box_matches_python():
    from dockflow_core.gridbox import box_from_coordinates

    rng = np.random.default_rng(5)
    coords = rng.uniform(-10, 10, size=(100, 3))
    cpp = dfb.grid_box(coords, padding=4.0)
    py = box_from_coordinates(coords, padding=4.0)
    assert np.allclose(cpp["center"], py.center)
    assert np.allclose(cpp["size"], py.size)


def test_grid_box_validation():
    with pytest.raises(ValueError):
        dfb.grid_box(np.zeros((0, 3)), padding=1.0)
    with pytest.raises(ValueError):
        dfb.grid_box(np.zeros((5, 4)), padding=1.0)


def test_pairwise_min_dist():
    a = np.array([[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]])
    b = np.array([[1.0, 0.0, 0.0], [50.0, 50.0, 50.0]])
    result = dfb.pairwise_min_dist(a, b)
    assert abs(result[0] - 1.0) < 1e-12
    # [10,10,10] to [1,0,0] is sqrt(9^2+10^2+10^2) = sqrt(281)
    assert abs(result[1] - np.sqrt(281.0)) < 1e-9


def test_min_contacts_sorted():
    rng = np.random.default_rng(9)
    a = rng.normal(size=(20, 3))
    b = rng.normal(size=(30, 3))
    contacts = dfb.min_contacts(a, b, cutoff=5.0)
    assert all(c[2] <= 5.0 for c in contacts)
    # sorted by ligand index then distance
    indices = [c[0] for c in contacts]
    assert indices == sorted(indices)


def test_direct_rmsd():
    a = np.zeros((6, 3))
    b = np.ones((6, 3))
    assert abs(dfb.direct_rmsd(a, b) - np.sqrt(3)) < 1e-12


def test_kabsch_rmsd_rigid_transform():
    rng = np.random.default_rng(21)
    points = rng.normal(scale=2.0, size=(35, 3))
    theta, phi = 0.9, 2.1
    r1 = np.array([[np.cos(theta), -np.sin(theta), 0],
                   [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    r2 = np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)],
                   [0, np.sin(phi), np.cos(phi)]])
    rotation = r1 @ r2
    moved = points @ rotation.T + rng.normal(scale=3.0, size=3)
    assert dfb.kabsch_rmsd(points, moved) < 1e-6


def test_kabsch_matches_numpy():
    from dockflow_core.analyzer import kabsch_rmsd as py_kabsch

    rng = np.random.default_rng(33)
    for _ in range(5):
        a = rng.normal(size=(25, 3))
        b = rng.normal(size=(25, 3))
        assert abs(dfb.kabsch_rmsd(a, b) - py_kabsch(a, b)) < 1e-6


def test_ligand_efficiency():
    assert dfb.ligand_efficiency(-9.0, 30) == pytest.approx(-0.3)
    with pytest.raises(ValueError):
        dfb.ligand_efficiency(-9.0, 0)


def test_box_corners():
    corners = dfb.box_corners(np.array([0.0, 0.0, 0.0]), np.array([2.0, 2.0, 2.0]))
    assert len(corners) == 8
    flat = np.array(corners)
    assert set(np.unique(flat[:, 0])) == {-1.0, 1.0}
