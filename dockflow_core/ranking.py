"""Consensus pose/ligand ranking (peer item 19) and IFP analysis CLI
support (peer item 18c).

Ranking philosophy: the Vina score alone is a weak ranker; combining it
with interaction-fingerprint representativeness (how typical the pose's
contacts are) and geometric cluster compactness (RMSD to the pose
cluster centroid) gives a consensus that is less sensitive to a single
scoring function's error.  All three components are heuristic; the
weights are configurable and every component is reported next to the
consensus so no number is taken on faith.

Formula (per pose, higher = better ranked)::

    consensus = w_score * z(affinity)
              + w_ifp   * mean_tanimoto_to_other_poses
              + w_geom  * (1 - rmsd_to_cluster_centroid / max_rmsd)

with ``z(affinity)`` the z-score over all ranked poses, the IFP term
the mean Tanimoto similarity of the pose's fingerprint to every other
analysed pose (interaction centrality), and the geometry term the
normalised closeness of the pose to its Kabsch pose-cluster centroid.
Default weights: score 0.5, IFP 0.3, geometry 0.2.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .analyzer import cluster_ifp, ifp_bits, ifp_vocabulary, tanimoto_ifp
from .utils import DockFlowError, get_logger

logger = get_logger("ranking")


def load_interactions(run_dir: str | Path) -> dict[str, list[dict]]:
    """analysis/interactions.json -> {ligand: [pose dicts]}."""
    path = Path(run_dir) / "analysis" / "interactions.json"
    if not path.is_file():
        raise DockFlowError(
            f"no interactions.json in {path}; run the analysis stage first"
        )
    return json.loads(path.read_text(encoding="utf-8"))


def poses_from_interactions(payload: dict[str, list[dict]]
                            ) -> list[dict[str, Any]]:
    """Flatten to [{ligand, pose_index, affinity, contacts, ifp_keys}]."""
    poses: list[dict[str, Any]] = []
    for ligand, entries in payload.items():
        for entry in entries:
            contacts = [
                {
                    "receptor_chain": c.get("receptor_chain", c.get("chain", "")),
                    "receptor_resname": c.get("receptor_resname", c.get("resname", "")),
                    "receptor_resseq": c.get("receptor_resseq", c.get("resseq", 0)),
                    "kind": c.get("kind", ""),
                }
                for c in entry.get("contacts", [])
            ]
            keys = {
                (contact["receptor_chain"], contact["receptor_resname"],
                 int(contact["receptor_resseq"] or 0), contact["kind"])
                for contact in contacts
            }
            poses.append({
                "ligand": ligand,
                "pose_index": entry.get("pose_index", 0),
                "affinity": entry.get("affinity", 0.0),
                "ifp_keys": keys,
            })
    return poses


def ifp_similarity_report(run_dir: str | Path, threshold: float = 0.7
                          ) -> dict[str, Any]:
    """Per-ligand IFP clustering + pose-pose Tanimoto table (item 18c).

    Returns a report dict (ligand -> clusters + similarity rows) and
    writes ``analysis/ifp_similarity.csv`` next to interactions.json.
    """
    payload = load_interactions(run_dir)
    report: dict[str, Any] = {}
    rows: list[list[Any]] = []
    for ligand, entries in payload.items():
        contact_sets = []
        for entry in entries:
            contacts = entry.get("contacts", [])
            keys = {
                (c.get("receptor_chain", c.get("chain", "")),
                 c.get("receptor_resname", c.get("resname", "")),
                 int(c.get("receptor_resseq", c.get("resseq", 0)) or 0),
                 c.get("kind", ""))
                for c in contacts
            }
            contact_sets.append(keys)
        vocabulary = ifp_vocabulary(contact_sets)
        bit_vectors = [ifp_bits(keys, vocabulary) for keys in contact_sets]
        clusters = cluster_ifp(bit_vectors, threshold=threshold)
        similarity_rows = []
        for i, vector_a in enumerate(bit_vectors):
            for j, vector_b in enumerate(bit_vectors):
                if j > i:
                    similarity_rows.append(
                        [i + 1, j + 1, round(tanimoto_ifp(vector_a, vector_b), 3)])
            for j, vector_b in enumerate(bit_vectors):
                if j >= i:
                    continue
                rows.append([ligand, i + 1, j + 1,
                             round(tanimoto_ifp(vector_a, vector_b), 3)])
        report[ligand] = {
            "n_poses": len(bit_vectors),
            "n_bits": len(vocabulary),
            "clusters": [
                {"cluster": index + 1, "poses": [p + 1 for p in members]}
                for index, members in enumerate(clusters)
            ],
            "pairwise_tanimoto": similarity_rows,
            "note": "heuristic clustering by interaction similarity "
                    f"(single linkage, Tanimoto >= {threshold}); no "
                    "thermodynamic meaning",
        }
    out_path = Path(run_dir) / "analysis" / "ifp_similarity.csv"
    if rows:
        with open(out_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["ligand", "pose_a", "pose_b", "tanimoto_ifp"])
            writer.writerows(rows)
    return report


def consensus_ranking(run_dir: str | Path,
                      score_weight: float = 0.5,
                      ifp_weight: float = 0.3,
                      geometry_weight: float = 0.2,
                      max_reference_rmsd: float = 4.0,
                      ) -> list[dict[str, Any]]:
    """Consensus ranking of every analysed pose (peer item 19).

    Components (all reported, all heuristic):
      * ``score_z``: z-score of the pose affinity over all analysed poses,
      * ``ifp_representativeness``: mean IFP Tanimoto to the other poses
        of the same ligand (interaction centrality),
      * ``geometry_closeness``: ``1 - rmsd_to_centroid / max`` where the
        centroid is that of the pose's Kabsch cluster (geometrically
        central poses score higher).

    Poses are ranked by ``consensus``; the caller aggregates per ligand
    (best pose wins).  Writes ``docking/ranked.csv``.
    """
    run_dir = Path(run_dir)
    poses = poses_from_interactions(load_interactions(run_dir))
    if not poses:
        raise DockFlowError("no analysed poses in interactions.json")

    # -- score component: z-score over all poses --------------------------
    affinities = [float(pose["affinity"]) for pose in poses]
    mean = sum(affinities) / len(affinities)
    variance = sum((a - mean) ** 2 for a in affinities) / len(affinities)
    std = variance ** 0.5 or 1.0

    # -- IFP component: per-ligand mean Tanimoto to other poses -----------
    for ligand in {pose["ligand"] for pose in poses}:
        ligand_poses = [p for p in poses if p["ligand"] == ligand]
        contact_sets = [p["ifp_keys"] for p in ligand_poses]
        vocabulary = ifp_vocabulary(contact_sets)
        bit_vectors = [ifp_bits(keys, vocabulary) for keys in contact_sets]
        for index, pose in enumerate(ligand_poses):
            if len(bit_vectors) > 1:
                similarities = [
                    tanimoto_ifp(bit_vectors[index], other)
                    for other_index, other in enumerate(bit_vectors)
                    if other_index != index
                ]
                pose["ifp_representativeness"] = sum(similarities) / len(similarities)
            else:
                pose["ifp_representativeness"] = 1.0

    # -- geometry component: RMSD to the pose-cluster centroid ------------
    _attach_geometry_closeness(run_dir, poses, max_reference_rmsd)

    ranked: list[dict[str, Any]] = []
    for pose in poses:
        score_z = (float(pose["affinity"]) - mean) / std
        # affinities are negative; more negative = better -> flip the z
        score_z = -score_z
        consensus = (score_weight * score_z
                     + ifp_weight * float(pose.get("ifp_representativeness", 0.0))
                     + geometry_weight * float(pose.get("geometry_closeness", 0.0)))
        ranked.append({
            "ligand": pose["ligand"],
            "pose_index": pose["pose_index"],
            "affinity_kcal_mol": float(pose["affinity"]),
            "score_z": round(score_z, 3),
            "ifp_representativeness": round(
                float(pose.get("ifp_representativeness", 0.0)), 3),
            "geometry_closeness": round(
                float(pose.get("geometry_closeness", 0.0)), 3),
            "consensus": round(consensus, 4),
        })
    ranked.sort(key=lambda row: (-row["consensus"], row["ligand"]))
    for position, row in enumerate(ranked, start=1):
        row["rank"] = position
    out_path = run_dir / "docking" / "ranked.csv"
    with open(out_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ranked[0].keys()))
        writer.writeheader()
        writer.writerows(ranked)
    return ranked


def _attach_geometry_closeness(run_dir: Path, poses: list[dict],
                               max_reference_rmsd: float) -> None:
    """Per ligand: RMSD of each pose to its Kabsch cluster centroid."""
    from .analyzer import kabsch_rmsd, pose_coordinates
    from .models import DockingResult, LigandRecord

    docking_dir = run_dir / "docking"
    for ligand in {pose["ligand"] for pose in poses}:
        out_path = docking_dir / f"{ligand}_out.pdbqt"
        if not out_path.is_file():
            for pose in poses:
                if pose["ligand"] == ligand:
                    pose["geometry_closeness"] = 0.0
            continue
        record = LigandRecord(identifier=ligand)
        result = DockingResult(ligand=record, ligand_name=ligand,
                               out_path=out_path)
        coords = pose_coordinates(result)
        ligand_poses = [p for p in poses if p["ligand"] == ligand]
        if len(coords) < len(ligand_poses):
            for pose in ligand_poses:
                pose["geometry_closeness"] = 0.0
            continue
        # cluster by Kabsch RMSD at 2.0 A (same heuristic as the report)

        n = len(ligand_poses)
        distances = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if coords[i].shape != coords[j].shape:
                    distance = max_reference_rmsd
                else:
                    distance = kabsch_rmsd(coords[i], coords[j])
                distances[i][j] = distances[j][i] = distance
        assigned: list[int | None] = [None] * n
        centroids: dict[int, list[int]] = {}
        for i in range(n):
            if assigned[i] is not None:
                continue
            members = [i]
            assigned[i] = len(centroids)
            frontier = [i]
            while frontier:
                current = frontier.pop()
                for j in range(n):
                    if assigned[j] is None and distances[current][j] <= 2.0:
                        assigned[j] = len(centroids)
                        members.append(j)
                        frontier.append(j)
            centroids[len(centroids)] = members
        for pose_index, pose in enumerate(ligand_poses):
            members = centroids[assigned[pose_index]]
            if len(members) > 1:
                spread = [
                    distances[pose_index][other] for other in members
                    if other != pose_index
                ]
                rmsd = sum(spread) / len(spread)
            else:
                rmsd = 0.0
            pose["geometry_closeness"] = round(
                max(0.0, 1.0 - rmsd / max_reference_rmsd), 3)
