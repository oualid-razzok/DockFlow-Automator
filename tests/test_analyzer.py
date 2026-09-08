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
    docking_score_efficiency,
    kabsch_rmsd,
    ligand_efficiency,
    write_analysis_json,
    write_contacts_csv,
)
from dockflow_core.models import DockingResult
from dockflow_core.pdbio import Atom
from dockflow_core.preparator import ReceptorPreparator, ReceptorPrepOptions

# The two crystal_rmsd tests below assert the real RDKit GetBestRMS path
# itself and therefore need RDKit installed (the ``prep`` extra).  Without
# it the analyzer falls back to element-greedy Kabsch BY DESIGN - that
# fallback is covered by test_crystal_rmsd_falls_back_when_rdkit_unavailable
# - so the rdkit-path assertions SKIP instead of failing in minimal
# environments (e.g. the CI gui job, which installs only [test,gui,viz]).
try:
    import rdkit  # noqa: F401  (imported only to detect availability)
except ImportError:
    _RDKIT_ABSENT = True
else:
    _RDKIT_ABSENT = False

requires_rdkit = pytest.mark.skipif(
    _RDKIT_ABSENT, reason="rdkit not installed (prep extra)")


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
    # renamed to docking_score_efficiency in 0.2.0 (audit item 14); the old
    # name stays as a deprecated alias for one release
    assert docking_score_efficiency(-9.0, 30) == pytest.approx(-0.3)
    assert docking_score_efficiency(-9.0, 0) is None
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
    assert rows[0].geom_hbond == 1 and rows[0].ionic == 1
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
                            geom_hbond=2, hydrophobic=3, ionic=1, metal=1,
                            closest=2.8)
    assert row.total == 7


# ---------------------------------------------------------------------------
# Analytic RMSD cases with known values (peer item 15)
# ---------------------------------------------------------------------------
def _atoms_from_coords(coords, element="C"):
    return [
        Atom(name=f"{element}{i + 1}", resname="LIG", chain="A", resseq=1,
             x=x, y=y, z=z, element=element)
        for i, (x, y, z) in enumerate(coords)
    ]


def test_kabsch_rmsd_identity_is_zero():
    coords = [(0.0, 0.0, 0.0), (1.5, 0.2, 0.0), (0.3, 1.4, 0.7), (2.1, 1.1, -0.4)]
    a = np.array(coords)
    assert kabsch_rmsd(a, a.copy()) == pytest.approx(0.0, abs=1e-12)


def test_kabsch_rmsd_translation_invariant():
    coords = [(0.0, 0.0, 0.0), (1.5, 0.2, 0.0), (0.3, 1.4, 0.7)]
    a = np.array(coords)
    b = a + np.array([10.0, -20.0, 30.0])
    assert kabsch_rmsd(a, b) == pytest.approx(0.0, abs=1e-9)


def test_kabsch_rmsd_rotation_invariant():
    # rigid rotation of the whole set: optimal superposition removes it
    coords = [(0.0, 0.0, 0.0), (1.5, 0.2, 0.0), (0.3, 1.4, 0.7), (2.1, 1.1, -0.4)]
    a = np.array(coords)
    theta = 0.9  # radians
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta), np.cos(theta), 0.0],
        [0.0, 0.0, 1.0],
    ])
    b = a @ rotation.T + np.array([5.0, 5.0, -3.0])
    # 1e-6: the C++ accelerated path computes in single-precision territory
    assert kabsch_rmsd(a, b) == pytest.approx(0.0, abs=1e-6)


def test_kabsch_rmsd_known_translation_without_superposition():
    # direct RMSD of a uniform 1 A displacement over 4 atoms = 1 A
    coords = np.array([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0),
                       (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)])
    shifted = coords + np.array([1.0, 0.0, 0.0])
    assert direct_rmsd(coords, shifted) == pytest.approx(1.0)


def test_kabsch_rmsd_scales_with_separation():
    # stretching a 4-point set: optimal superposition cannot absorb a
    # centroid-preserving stretch, and a bigger stretch must give a bigger
    # RMSD (monotonicity under a rigid superposition)
    a = np.array([(0.0, 0.0, 0.0), (1.0, 0.2, 0.0), (2.0, 0.1, 0.3), (3.0, 0.4, 0.1)])
    centroid = a.mean(axis=0)
    small = centroid + (a - centroid) * 1.05
    large = centroid + (a - centroid) * 1.5
    r_small = kabsch_rmsd(a, small)
    r_large = kabsch_rmsd(a, large)
    assert 0.0 < r_small < r_large
    # upper bound: RMSD cannot exceed the direct (unsuperposed) value
    assert r_large <= direct_rmsd(a, large) + 1e-9


def test_symmetry_tolerant_rmsd_identity_and_translation():
    coords = [(0.0, 0.0, 0.0), (1.5, 0.2, 0.0), (0.3, 1.4, 0.7)]
    from dockflow_core.analyzer import symmetry_tolerant_rmsd

    a = _atoms_from_coords(coords)
    assert symmetry_tolerant_rmsd(a, a) == pytest.approx(0.0, abs=1e-9)
    moved = _atoms_from_coords([(x + 4.0, y - 2.0, z + 1.0) for x, y, z in coords])
    assert symmetry_tolerant_rmsd(a, moved) == pytest.approx(0.0, abs=1e-9)


def test_symmetry_tolerant_rmsd_symmetric_pair_180_rotation():
    """180-degree flip of a symmetric two-atom group must cost 0 RMSD.

    Plain Kabsch with fixed correspondence would report 2.0 A (each atom
    swapped position); the symmetry-aware matcher must find the equivalent
    permutation and report 0.  This is the carboxylate-flip case.
    """
    from dockflow_core.analyzer import symmetry_tolerant_rmsd

    reference = _atoms_from_coords(
        [(-1.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 0.0, 3.0)],
        element="O",
    )
    reference += _atoms_from_coords([(0.0, 0.0, 4.0)], element="C")
    # same geometry, oxygens swapped (a 180-degree rotation of the O-O pair
    # around the C axis is achieved by permuting the two equivalent atoms)
    pose = _atoms_from_coords(
        [(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 0.0, 3.0)],
        element="O",
    )
    pose += _atoms_from_coords([(0.0, 0.0, 4.0)], element="C")
    assert symmetry_tolerant_rmsd(pose, reference) == pytest.approx(0.0, abs=1e-6)


def test_symmetry_tolerant_rmsd_displacement_bounds():
    """One atom displaced by 1 A over 4: bounded and monotone in the shift.

    The exact value is not analytic (the optimal superposition absorbs part
    of the displacement), but it must be positive, at most the
    no-superposition value (0.5 A) and larger for a bigger displacement.
    """
    from dockflow_core.analyzer import symmetry_tolerant_rmsd

    coords = [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (3.0, 0.0, 0.0), (4.5, 0.0, 0.0)]
    reference = _atoms_from_coords(coords)

    def _displaced(shift: float) -> list:
        return _atoms_from_coords(
            [(x, y + (shift if i == 0 else 0.0), z) for i, (x, y, z) in enumerate(coords)]
        )

    small = symmetry_tolerant_rmsd(_displaced(1.0), reference)
    large = symmetry_tolerant_rmsd(_displaced(2.0), reference)
    direct_small = (1.0 ** 2 / 4) ** 0.5
    assert 0.0 < small <= direct_small + 1e-9
    assert small < large


def test_symmetry_tolerant_rmsd_rejects_different_atom_sets():
    from dockflow_core.analyzer import symmetry_tolerant_rmsd

    with pytest.raises(ValueError, match="different atom sets"):
        symmetry_tolerant_rmsd(_atoms_from_coords([(0, 0, 0), (1, 0, 0)], "C"),
                               _atoms_from_coords([(0, 0, 0)], "C"))


def test_crystal_rmsd_returns_none_on_mismatch():
    from dockflow_core.analyzer import crystal_rmsd

    assert crystal_rmsd(_atoms_from_coords([(0, 0, 0)], "C"),
                        _atoms_from_coords([(0, 0, 0)], "N")) is None


def test_extract_reference_atoms(receptor_pdb_path: Path):
    from dockflow_core.analyzer import extract_reference_atoms

    atoms = extract_reference_atoms(receptor_pdb_path, "XK2" if False else "BEN")
    assert len(atoms) == 7
    assert all(a.resname.strip().upper() == "BEN" for a in atoms)
    with pytest.raises(ValueError, match="not found"):
        extract_reference_atoms(receptor_pdb_path, "ZZZ")


def test_analyze_docking_result_with_crystal_reference(
    docked_pdbqt_path: Path, receptor_pdb_path: Path, tmp_path: Path
):
    """Redocking validation end to end: crystal_rmsd lands on every pose."""
    from dockflow_core.analyzer import analyze_docking_result
    from dockflow_core.models import DockingResult, PoseRecord

    prep = ReceptorPreparator(
        ReceptorPrepOptions(engine="none", charge_model="zero", keep_resnames=["BEN"])
    ).prepare(receptor_pdb_path, tmp_path)
    # reference = the crystal pose of the same C+O fragment as the fixture
    # ligand (pose 1 coordinates == the crystal pose -> RMSD 0)
    reference = _atoms_from_coords([(12.5, 9.1, 8.1)], element="C")
    reference += _atoms_from_coords([(13.1, 9.8, 8.7)], element="O")
    result = DockingResult(
        ligand_name="lig",
        poses=[
            PoseRecord(model=1, affinity=-9.423),
            PoseRecord(model=2, affinity=-8.711),
            PoseRecord(model=3, affinity=-7.905),
        ],
        out_path=docked_pdbqt_path,
    )
    analyses = analyze_docking_result(
        result, prep.pdbqt_path, top_poses=3, reference_atoms=reference
    )
    assert result.poses[0].crystal_rmsd == pytest.approx(0.0, abs=1e-9)
    assert all(pose.crystal_rmsd is not None for pose in result.poses)
    # the fixture poses 2 and 3 are rigid copies (same C-O bond geometry),
    # so superposition-invariant RMSD is ~0 for them too - which is exactly
    # what a pose-recovery metric must report for pure rigid displacement
    assert result.poses[1].crystal_rmsd == pytest.approx(0.0, abs=1e-6)
    assert result.poses[2].crystal_rmsd == pytest.approx(0.0, abs=1e-6)
    assert analyses[0].crystal_rmsd == pytest.approx(result.poses[0].crystal_rmsd)


def test_analyze_docking_result_reference_mismatch_is_not_fatal(
    docked_pdbqt_path: Path, receptor_pdb_path: Path, tmp_path: Path
):
    """A non-comparable reference (decoys) leaves crystal_rmsd None."""
    from dockflow_core.analyzer import analyze_docking_result, extract_reference_atoms
    from dockflow_core.models import DockingResult, PoseRecord

    prep = ReceptorPreparator(
        ReceptorPrepOptions(engine="none", charge_model="zero", keep_resnames=["BEN"])
    ).prepare(receptor_pdb_path, tmp_path)
    benzene = extract_reference_atoms(receptor_pdb_path, "BEN")  # C6O vs C,O
    result = DockingResult(
        ligand_name="lig",
        poses=[PoseRecord(model=1, affinity=-9.423)],
        out_path=docked_pdbqt_path,
    )
    analyses = analyze_docking_result(
        result, prep.pdbqt_path, top_poses=3, reference_atoms=benzene
    )
    assert analyses  # analysis itself still works
    assert all(pose.crystal_rmsd is None for pose in result.poses)


def test_pose_cluster_summary_known_structure():
    """Two conformations + one distinct pose = three clusters with known members.

    ``cluster_poses`` measures Kabsch (superposition) RMSD, so rigid-body
    placement does NOT separate poses - a pose translated 5 A away is still
    the same binding conformation once optimally superposed.  Poses must
    therefore differ in *internal geometry* to land in separate clusters:
    here two poses share conformation A (tiny thermal noise), two share
    conformation B (an arm flipped), and one is a third distinct conformation.
    """
    from dockflow_core.analyzer import pose_cluster_summary

    rng = np.random.default_rng(5)
    base = rng.normal(size=(12, 3))
    arm_flipped = base.copy()
    arm_flipped[[2, 5, 8]] += np.array([0.0, 9.0, 9.0])  # different conformation
    distinct = base.copy()
    distinct[[1, 4, 7, 10]] += np.array([9.0, 0.0, -9.0])  # yet another one
    poses = [
        base,                                                   # conformation A
        base + rng.normal(scale=0.05, size=base.shape),          # A + noise
        arm_flipped,                                            # conformation B
        arm_flipped + rng.normal(scale=0.05, size=base.shape),   # B + noise
        distinct,                                               # conformation C
    ]
    affinities = [-9.0, -8.9, -8.0, -7.9, -6.0]
    clusters = pose_cluster_summary(poses, affinities, cutoff=2.0)
    assert len(clusters) == 3
    assert clusters[0]["poses"] == [1, 2]
    assert clusters[1]["poses"] == [3, 4]
    assert clusters[2]["poses"] == [5]
    assert clusters[0]["representative_pose"] == 1
    assert clusters[0]["mean_affinity"] == pytest.approx(-8.95)
    assert clusters[2]["intra_cluster_rmsd_spread"] == pytest.approx(0.0)
    # noise-only spread inside a cluster stays far below the 2 A cutoff
    assert 0.0 < clusters[0]["intra_cluster_rmsd_spread"] < 0.5


def test_kabsch_rmsd_matches_bindings_both_paths():
    """NumPy and C++ Kabsch agree on an analytic rotation (item 15/35).

    Only runs where the C++ accelerator is built (dev machines, the
    bindings-equivalence CI job); the plain test matrix skips it because
    it never compiles the optional extension module.
    """
    pytest.importorskip(
        "dockflow_bindings", reason="C++ accelerator not built in this job"
    )
    coords = np.array([
        [0.0, 0.0, 0.0], [1.5, 0.2, 0.0], [0.3, 1.4, 0.7],
        [2.1, 1.1, -0.4], [-0.7, 0.5, 1.2],
    ])
    theta = 0.7
    rotation = np.array([
        [np.cos(theta), -np.sin(theta), 0.0],
        [np.sin(theta), np.cos(theta), 0.0],
        [0.0, 0.0, 1.0],
    ])
    b = coords @ rotation.T + np.array([3.0, -2.0, 1.0])
    # the pure-numpy path (bindings temporarily hidden) and the accelerated
    # path must both report superposition-invariance
    assert kabsch_rmsd(coords, b) == pytest.approx(0.0, abs=1e-9)
    import builtins

    real_import = builtins.__import__

    def _no_bindings(name, *args, **kwargs):
        if name == "dockflow_bindings":
            raise ImportError("hidden for test")
        return real_import(name, *args, **kwargs)

    import dockflow_bindings  # noqa: F401  (present in dev/CI builds)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(builtins, "__import__", _no_bindings)
        # force the pure-numpy fallback path
        from dockflow_core import analyzer as _analyzer

        assert _analyzer.kabsch_rmsd(coords, b) == pytest.approx(0.0, abs=1e-9)


# -- Symmetry-aware RMSD via RDKit GetBestRMS (peer item 17) ---------------
@requires_rdkit
def test_crystal_rmsd_with_method_prefers_rdkit_path():
    """With RDKit importable the graph-automorphism method must be used."""
    import math

    from dockflow_core.analyzer import crystal_rmsd_with_method

    radius = 1.39 / (2 * math.sin(math.pi / 6))
    ring = _atoms_from_coords(
        [(radius * math.cos(i * math.pi / 3), radius * math.sin(i * math.pi / 3),
          0.0) for i in range(6)]
    )
    shifted = _atoms_from_coords(
        [(radius * math.cos((i + 1) * math.pi / 3),
          radius * math.sin((i + 1) * math.pi / 3), 0.0) for i in range(6)]
    )
    value, method = crystal_rmsd_with_method(shifted, ring)
    assert method == "rdkit-getbestrms-graph-automorphism"
    # a 60-degree ring relabel is the same chemistry: RMSD ~ 0
    assert value == pytest.approx(0.0, abs=1e-3)


def test_crystal_rmsd_falls_back_when_rdkit_unavailable(monkeypatch):
    """Without RDKit the element-greedy Kabsch heuristic is used and named."""
    from dockflow_core import analyzer

    monkeypatch.setattr(analyzer, "rdkit_symmetry_rmsd",
                        lambda pose, ref: None)
    coords = [(0.0, 0.0, 0.0), (1.5, 0.2, 0.0), (0.3, 1.4, 0.7)]
    a = _atoms_from_coords(coords)
    value, method = analyzer.crystal_rmsd_with_method(a, a)
    assert method == "element-greedy-kabsch"
    assert value == pytest.approx(0.0, abs=1e-9)


@requires_rdkit
def test_crystal_rmsd_getbestrms_benzoate_symmetry():
    """Benzoate-like ligand: swapped carboxylate oxygens cost ~0 RMSD.

    The two terminal oxygens of a carboxylate are graph-equivalent, so a
    180-degree flip of the COO group is the same chemistry (peer item 17:
    "a benzoic acid flipped 180 degrees").  Fixed-correspondence Kabsch
    would report >1 A; the automorphism-aware RMSD must report ~0.
    """
    import math

    from dockflow_core.analyzer import crystal_rmsd_with_method, kabsch_rmsd

    radius = 1.39 / (2 * math.sin(math.pi / 6))
    ring = [(radius * math.cos(i * math.pi / 3),
             radius * math.sin(i * math.pi / 3), 0.0) for i in range(6)]
    carboxy = [(0.0, 2.80, 0.0)]
    oxygens = [(1.2, 3.5, 0.0), (-1.2, 3.5, 0.0)]
    reference = (_atoms_from_coords(ring)
                 + _atoms_from_coords(carboxy)
                 + _atoms_from_coords(oxygens, element="O"))
    # pose: identical geometry, terminal oxygens swapped
    pose = (_atoms_from_coords(ring)
            + _atoms_from_coords(carboxy)
            + _atoms_from_coords(oxygens[::-1], element="O"))
    # fixed-correspondence Kabsch would penalise the swap
    fixed = kabsch_rmsd(
        [(a.x, a.y, a.z) for a in pose], [(a.x, a.y, a.z) for a in reference])
    assert fixed > 0.5
    value, method = crystal_rmsd_with_method(pose, reference)
    assert method == "rdkit-getbestrms-graph-automorphism"
    assert value == pytest.approx(0.0, abs=1e-3)


def test_crystal_rmsd_with_method_incomparable_returns_none():
    from dockflow_core.analyzer import crystal_rmsd_with_method

    a = _atoms_from_coords([(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)], element="C")
    b = _atoms_from_coords([(0.0, 0.0, 0.0)], element="N")
    assert crystal_rmsd_with_method(a, b) == (None, None)


def test_extract_reference_atoms_single_instance(tmp_path):
    """Two copies of the ligand: reference is ONE instance (largest)."""
    from dockflow_core.analyzer import extract_reference_atoms

    lines = []
    serial = 1
    for copy, chain in enumerate("AB"):
        for i in range(4):
            lines.append(
                f"ATOM  {serial:5d}  C{i+1}  XXX {chain}  {100 + copy:3d}    "
                f"{copy * 20.0 + i * 1.5:8.3f}    0.000    0.000  1.00  0.00           C"
            )
            serial += 1
    path = tmp_path / "multi.pdb"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    reference = extract_reference_atoms(path, "XXX")
    assert len(reference) == 4  # one instance only


# ---------------------------------------------------------------------------
# Symmetry-aware RMSD: analytically known cases (final pass, item 5)
# ---------------------------------------------------------------------------
def _rdkit_heavy_atoms(smiles: str, seed: int = 42) -> list[Atom]:
    """Deterministic heavy-Atom coordinates of an RDKit molecule."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"RDKit could not parse {smiles!r}"
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    assert AllChem.EmbedMolecule(mol, params) == 0
    conf = mol.GetConformer()
    atoms: list[Atom] = []
    for index, atom in enumerate(mol.GetAtoms()):
        if atom.GetAtomicNum() == 1:
            continue
        pos = conf.GetAtomPosition(index)
        symbol = atom.GetSymbol()
        atoms.append(Atom(name=f"{symbol}{index}", resname="LIG", chain="A",
                          resseq=1, x=pos.x, y=pos.y, z=pos.z,
                          element=symbol))
    return atoms


def _rotate(atoms: list[Atom], rotation: np.ndarray,
            translation=(0.0, 0.0, 0.0)) -> list[Atom]:
    rotated = []
    for a in atoms:
        xyz = np.asarray([a.x, a.y, a.z]) @ rotation + np.asarray(translation)
        rotated.append(
            Atom(name=a.name, resname=a.resname, chain=a.chain,
                 resseq=a.resseq, x=float(xyz[0]), y=float(xyz[1]),
                 z=float(xyz[2]), element=a.element))
    return rotated


def _molecule_axes(atoms: list[Atom]) -> tuple[np.ndarray, np.ndarray]:
    """(long in-plane axis, plane normal) of a planar molecule (SVD)."""
    coords = np.array([[a.x, a.y, a.z] for a in atoms])
    centered = coords - coords.mean(axis=0)
    _u, _s, vt = np.linalg.svd(centered, full_matrices=False)
    return vt[0], vt[2]  # largest variance = long axis, smallest = normal


def _rotation_about(axis: np.ndarray, degrees: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    theta = np.deg2rad(degrees)
    k = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return (np.eye(3) + np.sin(theta) * k
            + (1 - np.cos(theta)) * (k @ k))


def test_symmetry_rmsd_known_cases():
    """Analytic symmetry cases with hand-constructed molecules (item 5).

    Each case has a KNOWN answer, so the symmetry-aware RMSD
    (``crystal_rmsd_with_method`` -> RDKit GetBestRMS, graph automorphisms)
    is validated against analytic expectations, not just internal
    consistency:

    (a) benzene rigidly rotated 60 deg about its ring normal: a C6
        automorphism - symmetry-aware ~ 0; a naive UNALIGNED fixed-order
        RMSD (no superposition, the careless-script baseline) reports the
        full ~1.4 A ring shift;
    (b) naphthalene flipped 180 deg about its long in-plane axis: a C2
        automorphism - symmetry-aware ~ 0, naive unaligned > 0;
    (c) aspirin with the two carboxylic-acid oxygens relabeled (the atom
        order tools emit for symmetric groups): same chemistry -
        symmetry-aware ~ 0, fixed-correspondence Kabsch (superposed,
        identity order) > 0;
    (d) a chiral ligand with NO automorphisms (2-butanol), rigidly
        displaced: symmetry-aware == naive (~ 0 both, as any superposed
        RMSD must be for a rigid motion);
    (e) the same chiral ligand with one atom genuinely displaced so the
        RMSD is ~1.5 A: both methods must agree (no symmetry to exploit).
    """
    pytest.importorskip("rdkit")
    from dockflow_core.analyzer import crystal_rmsd_with_method, kabsch_rmsd

    def _naive(pose: list[Atom], ref: list[Atom]) -> float:
        return kabsch_rmsd([(a.x, a.y, a.z) for a in pose],
                           [(a.x, a.y, a.z) for a in ref])

    def _naive_unaligned(pose: list[Atom], ref: list[Atom]) -> float:
        return float(np.sqrt(np.mean([
            (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
            for a, b in zip(pose, ref, strict=True)])))

    def _symmetric(pose: list[Atom], ref: list[Atom]) -> float:
        value, method = crystal_rmsd_with_method(pose, ref)
        assert method == "rdkit-getbestrms-graph-automorphism"
        assert value is not None
        return value

    # (a) benzene rotated 60 deg about the ring normal ---------------------
    benzene = _rdkit_heavy_atoms("c1ccccc1")
    _long, normal = _molecule_axes(benzene)
    rotated = _rotate(benzene, _rotation_about(normal, 60.0),
                      (0.5, -0.4, 0.3))
    assert _symmetric(rotated, benzene) == pytest.approx(0.0, abs=1e-3)
    assert _naive_unaligned(rotated, benzene) > 0.5
    # (superposed fixed-order Kabsch is ~0 for ANY rigid motion - the
    # automorphism matters because real pose/reference atom orders and
    # partial-symmetry arrangements are NOT rigid motions of each other)

    # (b) naphthalene flipped 180 deg about the long in-plane axis ---------
    naphthalene = _rdkit_heavy_atoms("c1ccc2ccccc2c1")
    long_axis, _normal = _molecule_axes(naphthalene)
    flipped = _rotate(naphthalene, _rotation_about(long_axis, 180.0),
                      (-0.2, 0.6, 0.1))
    assert _symmetric(flipped, naphthalene) == pytest.approx(0.0, abs=1e-3)
    assert _naive_unaligned(flipped, naphthalene) > 0.3

    # (c) aspirin with the two acid oxygens swapped -----------------------
    aspirin = _rdkit_heavy_atoms("CC(=O)Oc1ccccc1C(=O)O")
    swapped = list(aspirin)
    for carbon in [a for a in aspirin if a.element == "C"]:
        oxygens = [(i, a) for i, a in enumerate(aspirin)
                   if a.element == "O" and _shares_bond(a, carbon)]
        if len(oxygens) != 2:
            continue
        # the COOH group (NOT the acetyl ester): both oxygens bonded
        # only to this carbon - swapping them is a graph automorphism
        others = [a for a in aspirin if a is not carbon]
        if any(_shares_bond(o, other) for _i, o in oxygens
               for other in others if other not in (o,)):
            continue
        (i, first), (j, second) = oxygens
        swapped[i], swapped[j] = second, first
        break
    else:
        pytest.fail("aspirin test: carboxylic acid oxygens not found")
    assert _symmetric(swapped, aspirin) == pytest.approx(0.0, abs=1e-3)
    assert _naive(swapped, aspirin) > 0.3

    # (d) chiral ligand, no automorphisms, rigid displacement -------------
    butanol = _rdkit_heavy_atoms("CC[C@H](C)O")
    displaced = _rotate(butanol, _rotation_about(np.array([0.2, 1.0, -0.3]),
                                                 37.0), (1.0, 2.0, -1.0))
    naive_d = _naive(displaced, butanol)
    sym_d = _symmetric(displaced, butanol)
    assert naive_d == pytest.approx(0.0, abs=1e-6)
    assert sym_d == pytest.approx(naive_d, abs=1e-3)

    # (e) same ligand, one atom displaced so RMSD ~ 1.5 A -----------------
    n_heavy = len(butanol)
    delta = 1.5 * np.sqrt(n_heavy)
    shifted = list(butanol)
    oxygen = next(i for i, a in enumerate(butanol) if a.element == "O")
    last = butanol[oxygen]
    shifted[oxygen] = Atom(name=last.name, resname="LIG", chain="A",
                           resseq=1, x=last.x, y=last.y + delta, z=last.z,
                           element="O")
    naive_e = _naive(shifted, butanol)
    sym_e = _symmetric(shifted, butanol)
    assert 1.0 < sym_e < 2.0
    assert abs(sym_e - naive_e) < 0.1  # both methods agree


def _shares_bond(a: Atom, b: Atom, cutoff: float = 2.0) -> bool:
    """Crude distance-based bonding for the aspirin test setup."""
    return (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2 < cutoff ** 2

