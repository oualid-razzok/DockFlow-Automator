"""Tests for interaction fingerprints (peer item 18) and consensus
ranking (peer item 19)."""

import json
from pathlib import Path

import pytest

from dockflow_core.analyzer import (
    PoseAnalysis,
    cluster_ifp,
    ifp_bits,
    ifp_vocabulary,
    pose_ifp_keys,
    tanimoto_ifp,
    write_ifp_outputs,
)
from dockflow_core.models import Contact
from dockflow_core.ranking import consensus_ranking, ifp_similarity_report, poses_from_interactions
from dockflow_core.utils import DockFlowError


def _contact(chain, resname, resseq, kind):
    return Contact(receptor_chain=chain, receptor_resname=resname,
                   receptor_resseq=resseq, kind=kind)


# -- item 18: IFP bitvectors --------------------------------------------------
def test_pose_ifp_keys_one_bit_per_residue_and_kind():
    contacts = [
        _contact("A", "ASP", 25, "hbond"),
        _contact("A", "ASP", 25, "hbond"),   # duplicate: same bit
        _contact("A", "ASP", 25, "ionic"),   # different kind: new bit
        _contact("B", "LYS", 90, "hydrophobic"),
    ]
    keys = pose_ifp_keys(contacts)
    assert keys == {("A", "ASP", 25, "hbond"), ("A", "ASP", 25, "ionic"),
                    ("B", "LYS", 90, "hydrophobic")}


def test_ifp_bits_and_vocabulary_are_stable():
    pose_a = {("A", "ASP", 25, "hbond"), ("B", "LYS", 90, "hydrophobic")}
    pose_b = {("A", "ASP", 25, "ionic")}
    vocabulary = ifp_vocabulary([pose_a, pose_b])
    assert len(vocabulary) == 3
    bits_a = ifp_bits(pose_a, vocabulary)
    assert bits_a == [1, 0, 1] or sum(bits_a) == 2
    bits_b = ifp_bits(pose_b, vocabulary)
    assert sum(bits_b) == 1


def test_tanimoto_ifp_identical_disjoint_and_overlap():
    vec_a = [1, 0, 1, 1]
    assert tanimoto_ifp(vec_a, vec_a) == pytest.approx(1.0)
    assert tanimoto_ifp(vec_a, [0, 1, 0, 0]) == pytest.approx(0.0)
    # 2 shared bits, 1+2 unique -> 2/3
    assert tanimoto_ifp(vec_a, [1, 0, 1, 0]) == pytest.approx(2.0 / 3.0)
    assert tanimoto_ifp([0, 0], [0, 0]) == 1.0  # both empty: identical


def test_cluster_ifp_single_linkage():
    # poses 1,2 share 2 bits; 3 is disjoint
    vectors = [[1, 1, 0], [1, 1, 0], [0, 0, 1]]
    clusters = cluster_ifp(vectors, threshold=0.5)
    assert sorted(clusters) == [[0, 1], [2]]
    # threshold 1.0: only exact duplicates join
    assert cluster_ifp(vectors, threshold=1.0) == [[0, 1], [2]]
    assert cluster_ifp([[1, 0], [0, 1]], threshold=0.5) == [[0], [1]]


def test_write_ifp_outputs_csv_and_dat(tmp_path):
    analyses = [
        PoseAnalysis(pose_index=1, contacts=[
            _contact("A", "ASP", 25, "hbond")]),
        PoseAnalysis(pose_index=2, contacts=[
            _contact("A", "ASP", 25, "hbond"),
            _contact("B", "LYS", 90, "hydrophobic")]),
    ]
    written = write_ifp_outputs("lig1", analyses, tmp_path)
    csv_path = tmp_path / "lig1_ifp.csv"
    assert csv_path in written
    lines = csv_path.read_text().splitlines()
    assert lines[0].startswith("pose,A:ASP25:hbond")
    assert lines[1] == "1,1,0"
    assert lines[2] == "2,1,1"
    dat_files = [p for p in written if p.suffix == ".dat"]
    if dat_files:  # RDKit available
        assert len(dat_files) == 2
        from rdkit.DataStructs import ExplicitBitVect

        vector = ExplicitBitVect(bytes(dat_files[0].read_bytes()))
        assert vector.GetNumBits() == 2
        assert vector.GetNumOnBits() == 1


# -- item 18c: ifp similarity report ------------------------------------------
def _write_interactions(run_dir: Path, payload: dict) -> None:
    (run_dir / "analysis").mkdir(parents=True, exist_ok=True)
    (run_dir / "analysis" / "interactions.json").write_text(
        json.dumps(payload), encoding="utf-8")


def _pose_payload(poses_spec):
    """{ligand: [(pose_index, affinity, contacts), ...]} -> interactions.json payload."""
    payload = {}
    for ligand, poses in poses_spec.items():
        entries = []
        for pose_index, affinity, contacts in poses:
            entries.append({
                "pose_index": pose_index,
                "affinity": affinity,
                "contacts": [
                    {"receptor_chain": c[0], "receptor_resname": c[1],
                     "receptor_resseq": c[2], "kind": c[3]}
                    for c in contacts
                ],
            })
        payload[ligand] = entries
    return payload


def test_ifp_similarity_report_clusters_and_csv(tmp_path):
    _write_interactions(tmp_path, _pose_payload({
        "ligA": [
            (1, -9.0, [("A", "ASP", 25, "hbond"), ("A", "ASP", 29, "ionic")]),
            (2, -8.5, [("A", "ASP", 25, "hbond")]),
            (3, -7.0, [("B", "LYS", 90, "hydrophobic")]),
        ],
    }))
    report = ifp_similarity_report(tmp_path, threshold=0.5)
    entry = report["ligA"]
    assert entry["n_poses"] == 3
    clusters = [c["poses"] for c in entry["clusters"]]
    assert [1, 2] in clusters and [3] in clusters
    csv_path = tmp_path / "analysis" / "ifp_similarity.csv"
    assert csv_path.is_file()
    # tanimoto(1,2) = 1 shared / (2 + 1 - 1) = 0.5
    body = csv_path.read_text().splitlines()
    assert "0.5" in " ".join(body)


def test_ifp_similarity_report_missing_json(tmp_path):
    with pytest.raises(DockFlowError):
        ifp_similarity_report(tmp_path)


# -- item 19: consensus ranking ------------------------------------------------
def _pose_pdbqt(ligand: str, n_models: int) -> Path:
    text = []
    for model in range(1, n_models + 1):
        text.append(f"MODEL {model}")
        text.append(
            "ATOM      1  C   LIG A   1      "
            f"{0.0:8.3f}{0.0:8.3f}{model * 0.5:8.3f}"
            "  1.00  0.00     0.000 C")
        text.append(
            "ATOM      2  N   LIG A   1      "
            f"{1.5:8.3f}{0.0:8.3f}{model * 0.5:8.3f}"
            "  1.00  0.00     0.000 N")
        text.append("ENDMDL")
    path = tmp_global / "docking" / f"{ligand}_out.pdbqt"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(text) + "\n", encoding="utf-8")
    return path


tmp_global = None


def test_consensus_ranking_orders_poses(tmp_path):
    global tmp_global
    tmp_global = tmp_path
    _pose_pdbqt("ligA", 3)
    _pose_pdbqt("ligB", 2)
    _write_interactions(tmp_path, _pose_payload({
        "ligA": [
            (1, -10.0, [("A", "ASP", 25, "hbond"), ("A", "ASP", 29, "ionic")]),
            (2, -9.0, [("A", "ASP", 25, "hbond"), ("A", "ASP", 29, "ionic")]),
            (3, -6.0, [("B", "LYS", 90, "hydrophobic")]),
        ],
        "ligB": [
            (1, -8.0, [("A", "ASP", 25, "hbond")]),
            (2, -7.5, [("B", "LYS", 90, "hydrophobic")]),
        ],
    }))
    ranked = consensus_ranking(tmp_path)
    assert len(ranked) == 5
    assert ranked[0]["rank"] == 1
    # the strongest affinity with the most shared interactions leads
    assert ranked[0]["ligand"] == "ligA"
    fields = {"ligand", "pose_index", "affinity_kcal_mol", "score_z",
              "ifp_representativeness", "geometry_closeness", "consensus",
              "rank"}
    assert fields <= set(ranked[0].keys())
    csv_path = tmp_path / "docking" / "ranked.csv"
    assert csv_path.is_file()
    # consensus is a weighted sum of the components
    for row in ranked:
        expected = (0.5 * row["score_z"] + 0.3 * row["ifp_representativeness"]
                    + 0.2 * row["geometry_closeness"])
        assert row["consensus"] == pytest.approx(expected, abs=0.005)


def test_consensus_ranking_weights_are_configurable(tmp_path):
    global tmp_global
    tmp_global = tmp_path
    _pose_pdbqt("ligC", 2)
    _write_interactions(tmp_path, _pose_payload({
        "ligC": [
            (1, -9.0, [("A", "ASP", 25, "hbond")]),
            (2, -8.0, [("A", "ASP", 25, "hbond")]),
        ],
    }))
    ranked = consensus_ranking(tmp_path, score_weight=1.0,
                               ifp_weight=0.0, geometry_weight=0.0)
    # score-only ranking: best affinity first
    assert ranked[0]["affinity_kcal_mol"] <= ranked[1]["affinity_kcal_mol"]
    for row in ranked:
        assert row["consensus"] == pytest.approx(row["score_z"], abs=0.005)


def test_poses_from_interactions_flattens_payload():
    poses = poses_from_interactions({
        "lig": [{"pose_index": 1, "affinity": -9.0,
                 "contacts": [{"receptor_chain": "A", "receptor_resname": "ASP",
                               "receptor_resseq": 25, "kind": "hbond"}]}]
    })
    assert poses[0]["ligand"] == "lig"
    assert ("A", "ASP", 25, "hbond") in poses[0]["ifp_keys"]
