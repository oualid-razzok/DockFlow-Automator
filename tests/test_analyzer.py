"""Analyzer tests: contacts, RMSD, clustering, reports."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from dockflow_core.analyzer import (
    ResidueContactRow,
    analyze_docking_result,
    analyze_interactions,
    classify_pair,
    cluster_poses,
    contact_summary,
    direct_rmsd,
    kabsch_rmsd,
    ligand_efficiency,
    write_analysis_json,
    write_contacts_csv,
)
from dockflow_core.models import DockingResult
from dockflow_core.pdbio import Atom
from dockflow_core.preparator import ReceptorPreparator, ReceptorPrepOptions


def _atom(name, resname, chain, resseq, x, y, z, element, atom_type):
    return Atom(name=name, resname=resname, chain=chain, resseq=resseq,
                x=x, y=y, z=z, element=element, atom_type=atom_type)


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------
def test_classify_pair_rules():
    assert classify_pair("HD", "OA", 3.0, "VAL", "O") == "hbond"
    assert classify_pair("NA", "HD", 3.2, "SER", "HG") == "hbond"
    assert classify_pair("HD", "OA", 4.0, "VAL", "O") is None  # too far
    assert classify_pair("A", "C", 4.0, "LEU", "CD1") == "hydrophobic"
    assert classify_pair("A", "OA", 4.0, "LEU", "O") is None
    assert classify_pair("N", "OA", 4.0, "ASP", "OD1") == "ionic"
    assert classify_pair("OA", "N", 4.4, "LYS", "NZ") == "ionic"
    assert classify_pair("OA", "Zn", 2.5, "ZN", "ZN") == "metal"
    assert classify_pair("C", "C", 3.0, "GLY", "CA") == "hydrophobic"


def test_analyze_interactions_hbond():
    ligand = [_atom("O1", "LIG", "A", 1, 0.0, 0.0, 0.0, "O", "OA")]
    receptor = [
        _atom("H", "SER", "A", 42, 2.0, 0.0, 0.0, "H", "HD"),
        _atom("CA", "GLY", "A", 9, 9.0, 9.0, 9.0, "C", "C"),
    ]
    contacts = analyze_interactions(ligand, receptor, cutoff=5.0)
    assert len(contacts) == 1
    contact = contacts[0]
    assert contact.kind == "hbond"
    assert contact.receptor_resname == "SER" and contact.receptor_resseq == 42
    assert contact.ligand_atom_type == "OA"
    assert contact.distance == 2.0


def test_analyze_interactions_metal_and_ionic():
    ligand = [
        _atom("O1", "LIG", "A", 1, 0.0, 0.0, 0.0, "O", "OA"),
        _atom("N1", "LIG", "A", 1, 10.0, 0.0, 0.0, "N", "N"),
    ]
    receptor = [
        _atom("ZN", "ZN", "A", 99, 2.0, 0.0, 0.0, "Zn", "Zn"),
        _atom("OD1", "ASP", "A", 25, 12.0, 0.0, 0.0, "O", "OA"),
    ]
    contacts = analyze_interactions(ligand, receptor, cutoff=6.0)
    kinds = {c.kind for c in contacts}
    assert "metal" in kinds and "ionic" in kinds


def test_analyze_interactions_empty():
    assert analyze_interactions([], []) == []


# ---------------------------------------------------------------------------
# RMSD
# ---------------------------------------------------------------------------
def test_direct_rmsd():
    a = np.zeros((4, 3))
    b = np.ones((4, 3))
    assert abs(direct_rmsd(a, b) - np.sqrt(3)) < 1e-12


def test_direct_rmsd_shape_mismatch():
    with pytest.raises(ValueError):
        direct_rmsd(np.zeros((3, 3)), np.zeros((4, 3)))


def test_kabsch_rmsd_rigid_transform():
    rng = np.random.default_rng(7)
    points = rng.normal(scale=2.5, size=(40, 3))
    theta = 1.234
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1],
    ])
    moved = points @ rotation.T + np.array([4.0, -2.0, 9.0])
    assert kabsch_rmsd(points, moved) < 1e-8


def test_kabsch_rmsd_reflection():
    rng = np.random.default_rng(7)
    points = rng.normal(size=(30, 3))
    mirrored = points * np.array([1.0, -1.0, 1.0])
    # a mirror is not a rotation; rmsd after proper rotation is nonzero
    assert kabsch_rmsd(points, mirrored) > 1e-3


def test_kabsch_matches_bindings_if_present():
    pytest.importorskip("dockflow_bindings")
    rng = np.random.default_rng(11)
    points = rng.normal(size=(25, 3))
    theta = 0.3
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1],
    ])
    moved = points @ rotation.T + np.array([1.0, 2.0, 3.0])
    import dockflow_bindings as dfb

    assert abs(dfb.kabsch_rmsd(points, moved) - kabsch_rmsd(points, moved)) < 1e-6


# ---------------------------------------------------------------------------
# Clustering & efficiency
# ---------------------------------------------------------------------------
def test_cluster_poses_separates():
    rng = np.random.default_rng(3)
    base = rng.normal(size=(10, 3))
    # Kabsch removes rigid transformations, so poses must differ in internal
    # geometry.  Move single atoms far away (as a docking pose with one arm
    # flipped would) to create genuinely distinct poses.
    pose_one_moved = base.copy()
    pose_one_moved[2] += np.array([12.0, 0.0, 0.0])
    pose_three_moved = base.copy()
    pose_three_moved[[1, 5, 8]] += np.array([0.0, 9.0, 9.0])
    poses = [
        base,                            # cluster 1
        base + rng.normal(scale=0.05, size=base.shape),  # cluster 1 (noise)
        pose_one_moved,                  # distinct pose
        pose_three_moved,                # another distinct pose
    ]
    clusters = cluster_poses(poses, cutoff=2.0)
    assert len(clusters) == 3
    assert clusters[0][:2] == [0, 1]


def test_cluster_poses_empty():
    assert cluster_poses([]) == []


def test_ligand_efficiency():
    assert ligand_efficiency(-9.0, 30) == pytest.approx(-0.3)
    assert ligand_efficiency(-9.0, 0) is None


# ---------------------------------------------------------------------------
# Contact summary + writers
# ---------------------------------------------------------------------------
def test_contact_summary():
    from dockflow_core.models import Contact

    contacts = [
        Contact(receptor_resname="ASP", receptor_chain="A", receptor_resseq=25,
                kind="hbond", distance=2.9),
        Contact(receptor_resname="ASP", receptor_chain="A", receptor_resseq=25,
                kind="ionic", distance=3.6),
        Contact(receptor_resname="LEU", receptor_chain="A", receptor_resseq=10,
                kind="hydrophobic", distance=4.2),
    ]
    rows = contact_summary(contacts)
    assert rows[0].resname == "ASP" and rows[0].total == 2
    assert rows[0].hbonds == 1 and rows[0].ionic == 1
    assert rows[1].resname == "LEU"


def test_write_contacts_csv(tmp_path: Path):
    from dockflow_core.models import Contact

    contacts = [
        Contact(ligand_atom_index=0, ligand_atom_name="O1",
                ligand_atom_type="OA", receptor_atom_name="OD1",
                receptor_resname="ASP", receptor_chain="A", receptor_resseq=25,
                receptor_atom_type="OA", distance=2.9, kind="hbond"),
    ]
    path = write_contacts_csv(contacts, tmp_path / "contacts.csv")
    content = path.read_text(encoding="utf-8")
    assert "hbond" in content and "ASP" in content and "2.90" in content


def test_write_analysis_json(tmp_path: Path):
    from dockflow_core.analyzer import PoseAnalysis

    analyses = [PoseAnalysis(pose_index=1, affinity=-9.0, num_contacts=2)]
    path = write_analysis_json(analyses, tmp_path / "analysis.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["poses"][0]["affinity"] == -9.0


# ---------------------------------------------------------------------------
# End-to-end analysis against a prepared receptor
# ---------------------------------------------------------------------------
def test_analyze_docking_result(tmp_path: Path, receptor_pdb_path: Path,
                                docked_pdbqt_path: Path):
    prep = ReceptorPreparator(
        ReceptorPrepOptions(engine="none", charge_model="zero",
                            keep_resnames=["BEN"])
    ).prepare(receptor_pdb_path, tmp_path)
    result = DockingResult(
        ligand_name="lig",
        poses=[
            __import__("dockflow_core.models", fromlist=["PoseRecord"]).PoseRecord(
                model=1, affinity=-9.423),
        ],
        out_path=docked_pdbqt_path,
    )
    analyses = analyze_docking_result(result, prep.pdbqt_path, top_poses=3)
    assert len(analyses) == 3
    first = analyses[0]
    assert first.affinity == -9.423
    assert first.num_contacts >= 1
    assert first.residue_rows  # interactions with the mini receptor


def test_residue_contact_row_totals():
    row = ResidueContactRow(chain="A", resname="GLU", resseq=12,
                            hbonds=2, hydrophobic=3, ionic=1, metal=1,
                            closest=2.8)
    assert row.total == 7
