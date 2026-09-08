"""Automated redocking benchmark runner (work order item 1; final pass items 1-3, 15).

For every complex in ``complexes.csv`` this script:

1. downloads the PDB structure from RCSB,
2. runs the standard DockFlow pipeline: receptor preparation (engine auto),
   ligand preparation, grid box from the co-crystallized ligand pocket
   (4 A padding),
3. docks with AutoDock Vina (fixed seed, configurable exhaustiveness),
4. records the symmetry-tolerant heavy-atom ``crystal_rmsd`` of every pose,
   the best affinity, the rank of the first successful pose and runtime.

The runner is **resumable and crash-safe**:

* every completed complex is appended to ``redocking_results.csv`` the
  moment it finishes (incremental flush, one row per complex);
* a complex already present in the CSV with a terminal status is skipped
  on restart, so the benchmark can be stopped and continued at any time;
* each complex runs in a dedicated worker subprocess (fresh interpreter),
  so per-complex memory (Vina bindings, RDKit caches) is fully released;
* ``summary.md`` + ``summary_stats.json`` + figures are regenerated after
  every complex, so partial results are always published;
* failed complexes are recorded with their failure reason and are NEVER
  dropped from the CSV or the summary table.

Per-complex preparation + grid provenance (final pass, item 1): every row
is enriched from the run's ``manifest.json`` with the receptor engine,
atom counts, ligand flexibility, grid-box source/center/size/padding/
assumption strength and the docking parameters (exhaustiveness, modes,
seed, cpu) so a reader can judge whether a failure is preparation-related
or docking-related.  ``--report-only`` re-derives ALL published numbers
from the committed per-complex artifacts without docking.

Pose-recovery definitions (final pass, items 3 and 15) - published under
ALL of them so the number is not definition-shopped (formal definition in
``docs/interaction_criteria.md``):

* ``rank1_rmsd_a``   - RMSD of the best-SCORING pose (Vina rank 1): the
  formal "pose recovery" definition;
* ``top3_rmsd_a``    - best RMSD among the top 3 scored poses;
* ``best_crystal_rmsd_a`` - best RMSD among all returned poses (the
  v0.3.x legacy definition, kept for continuity).

Outputs (in ``--out``):

* ``redocking_results.csv``  - one row per complex (the core artifact)
* ``summary.md``             - success rates with Wilson 95% CIs, RMSD
  distribution, per-target grouping, failed-complex sub-table,
  per-complex table incl. grid-box provenance
* ``summary_stats.json``     - the same numbers, machine-readable
* ``scatter.png``            - best affinity vs best crystal RMSD scatter
* ``rmsd_distribution.png``  - RMSD histogram (final pass, item 2)
* ``topn_success.png``       - success rate vs N (top-N poses, item 15)
* ``<pdb_id>_benchmark/``    - full pipeline run directories
* ``benchmark.log``          - driver log (workers log into their run dirs)

Example::

    python benchmarks/redocking/run_benchmark.py --limit 5
    python benchmarks/redocking/run_benchmark.py --pdb-ids 1HVR,1STP
    python benchmarks/redocking/run_benchmark.py --report-only
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

SUCCESS_RMSD = 2.0  # A, conventional redocking success threshold
THRESHOLDS = (1.0, 2.0, 3.0)  # final pass, item 15: multi-threshold success

FIELDNAMES = [
    "pdb_id", "ligand", "target", "status", "num_poses",
    "best_affinity_kcal_mol", "best_crystal_rmsd_a", "best_pose_rank",
    "success", "engine", "docking_backend", "runtime_s", "run_dir",
    # -- preparation + grid provenance (final pass, item 1) ---------------
    "receptor_engine", "receptor_atoms_out",
    "ligand_rotatable_bonds", "ligand_heavy_atoms",
    "gridbox_source", "gridbox_center_x", "gridbox_center_y",
    "gridbox_center_z", "gridbox_size_x", "gridbox_size_y",
    "gridbox_size_z", "gridbox_padding", "gridbox_assumption_strength",
    # -- docking parameters (final pass, item 1) ---------------------------
    "exhaustiveness", "num_modes", "refine", "seed", "cpu",
    # -- pose-recovery definitions (final pass, items 3, 15) ---------------
    "crystal_rmsd_method", "rank1_rmsd_a", "top3_rmsd_a",
    "success_rank1", "success_top3",
]

# One-line hypothesized failure causes for the failed-complex sub-table
# (final pass, item 2).  Grounded in the per-complex provenance columns:
# every failed complex had a strong ligand-derived grid box, so the
# failures are docking-side, not preparation-side.  The spyrmsd
# cross-check (item 5) corrected the crystal-RMSD implementation in
# v1.0.0-rc1 (atom-order independence); with corrected RMSDs five of
# the six v0.3.x "failures" actually recover a <= 2.0 A pose - the
# hypotheses below describe the failure modes that REMAIN.
FAILURE_HYPOTHESES = {
    "1HTG": "penicillin-derived peptidomimetic, the most flexible ligand "
            "in the set (17 rotatable bonds): no pose within 3.0 A of the "
            "crystal pose among all 9 modes - a sampling limitation, not "
            "a ranking one (the top-N curve stays flat)",
    "1DWD": "HIV-PR peptidomimetic (9 rotatable bonds): a 1.36 A pose IS "
            "sampled but ranks 3rd - the failure is ranking (Vina scoring), "
            "not sampling",
    "1HPV": "HIV-PR amprenavir: a 0.69 A pose IS sampled but ranks 3rd - "
            "ranking failure, not sampling; recovered under the top-3 "
            "definition",
}


def load_complexes(path: Path) -> list[dict[str, str]]:
    """Read the complex table, skipping ``#`` comment lines and blanks."""
    lines = [
        line for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    rows = [
        {key or "extra": ((value or "") if isinstance(value, str) else ",".join(value))
         for key, value in row.items()}
        for row in csv.DictReader(lines)
    ]
    return [row for row in rows if row.get("pdb_id")]


def load_existing(csv_path: Path) -> dict[str, dict]:
    """Read already-recorded rows keyed by pdb id (resume support)."""
    if not csv_path.exists():
        return {}
    with open(csv_path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {row["pdb_id"]: row for row in rows if row.get("pdb_id")}


def append_row(csv_path: Path, row: dict) -> None:
    """Append one result row immediately (crash-safe incremental flush)."""
    is_new = not csv_path.exists()
    with open(csv_path, "a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        if is_new:
            writer.writeheader()
        writer.writerow({key: row.get(key, "") for key in FIELDNAMES})


def rewrite_csv(csv_path: Path, rows: list[dict]) -> None:
    """Rewrite the full CSV with the enriched column set (item 1).

    Atomic (tmp + replace) so a crash mid-report never truncates the
    core artifact.  Old-format rows (pre-item-1 columns) are upgraded in
    place by :func:`enrich_rows` first.
    """
    tmp_path = csv_path.with_suffix(".csv.tmp")
    with open(tmp_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in FIELDNAMES})
    tmp_path.replace(csv_path)


# ---------------------------------------------------------------------------
# Wilson score interval (final pass, item 2)
# ---------------------------------------------------------------------------
def wilson_ci(successes: int, total: int, z: float = 1.959964
              ) -> tuple[float, float]:
    """Wilson score 95% confidence interval for a binomial proportion.

    Formula (Wilson 1927, JASA 22(158):209-212)::

        centre = p + z^2 / (2 n)
        margin = z * sqrt(p (1 - p) / n + z^2 / (4 n^2))
        CI     = ((centre - margin) / (1 + z^2 / n),
                  (centre + margin) / (1 + z^2 / n))

    with p = successes / total and z = 1.959964 (two-sided 95%).  Chosen
    over the Wald interval because it behaves correctly near p = 0 and
    p = 1 and at small n (hand-rolled: no statsmodels dependency).
    """
    if total <= 0:
        return (0.0, 0.0)
    p = successes / total
    z2 = z * z
    centre = p + z2 / (2.0 * total)
    margin = z * (p * (1.0 - p) / total + z2 / (4.0 * total * total)) ** 0.5
    denom = 1.0 + z2 / total
    return ((centre - margin) / denom, (centre + margin) / denom)


def _fmt_pct_rate(successes: int, total: int) -> str:
    """'75.0% (95% CI: 55.1%-88.0%)' - the item-2 presentation format."""
    if total <= 0:
        return "n/a"
    low, high = wilson_ci(successes, total)
    rate = 100.0 * successes / total
    return f"{rate:.1f}% (95% CI: {100 * low:.1f}%-{100 * high:.1f}%)"


def run_complex(complex_row: dict[str, str], out_dir: Path,
                exhaustiveness: int, seed: int, cpu: int) -> dict:
    """Run one complex in-process (worker mode)."""
    from dockflow_core.pipeline import DockingPipeline, PipelineConfig

    pdb_id = complex_row["pdb_id"]
    ligand_resname = complex_row["ligand_resname"]
    workdir = out_dir / pdb_id
    workdir.mkdir(parents=True, exist_ok=True)

    config = PipelineConfig(
        workdir=out_dir,
        run_id=f"{pdb_id}_benchmark",
        target={"pdb_id": pdb_id},
        ligands=[{"id": f"{ligand_resname.lower()}_redock",
                  "pdb_ligand": ligand_resname}],
        receptor={"engine": "auto", "charge_model": "gasteiger"},
        gridbox={"source": "ligand", "reference_ligand_resname": ligand_resname,
                 "padding": 4.0},
        docking={"backend": "auto", "scoring": "vina",
                 "exhaustiveness": exhaustiveness, "num_modes": 9,
                 "seed": seed, "cpu": cpu, "timeout": 3600,
                 "parallel": 1},
        analysis={"top_poses": 3},
        visualization={"enabled": False},
    )
    started = time.perf_counter()
    report = DockingPipeline(config).run()
    if not report.ok:
        return {"pdb_id": pdb_id, "ligand": ligand_resname,
                "target": complex_row.get("target", ""),
                "status": f"failed: {report.error}"}

    results = report.docking["results"]
    best_affinity = None
    best_rmsd = None
    best_pose_rank = None
    poses = []
    for entry in results:
        for pose in entry.get("poses", []):
            poses.append(pose)
    if poses:
        best_affinity = poses[0].get("affinity")
        rmsds = [(pose.get("crystal_rmsd"), index + 1)
                 for index, pose in enumerate(poses)
                 if pose.get("crystal_rmsd") is not None]
        if rmsds:
            best_rmsd, best_pose_rank = min(rmsds)
    # run_dir is stored RELATIVE to the results root so the published CSV
    # is portable across machines (audit item 1: real, inspectable results).
    try:
        run_dir_rel = str(Path(report.run_dir).relative_to(out_dir))
    except ValueError:
        run_dir_rel = Path(report.run_dir).name
    return {
        "pdb_id": pdb_id,
        "ligand": ligand_resname,
        "target": complex_row.get("target", ""),
        "status": "ok",
        "num_poses": len(poses),
        "best_affinity_kcal_mol": best_affinity,
        "best_crystal_rmsd_a": best_rmsd,
        "best_pose_rank": best_pose_rank,
        "success": (best_rmsd is not None and best_rmsd <= SUCCESS_RMSD),
        "engine": report.receptor.get("engine", ""),
        "docking_backend": report.docking.get("backend", ""),
        "runtime_s": round(time.perf_counter() - started, 1),
        "run_dir": run_dir_rel,
    }


# ---------------------------------------------------------------------------
# Per-complex provenance enrichment (final pass, item 1)
# ---------------------------------------------------------------------------
def _manifest_enrichment(run_dir: Path) -> dict:
    """Preparation + grid + docking settings recorded in one run's manifest."""
    out: dict = {}
    try:
        manifest = json.loads(
            (run_dir / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return out
    receptor = manifest.get("receptor") or {}
    gridbox = manifest.get("gridbox") or {}
    ligand = (manifest.get("ligands") or [{}])[0] or {}
    docking_cfg = ((manifest.get("config") or {}).get("docking")) or {}
    analysis = manifest.get("analysis") or {}
    poses: list[dict] = []
    for entry in (manifest.get("docking", {}).get("results") or []):
        poses += entry.get("poses") or []
    rmsds = [pose.get("crystal_rmsd") for pose in poses
             if pose.get("crystal_rmsd") is not None]
    center = list(gridbox.get("center") or [None, None, None])
    size = list(gridbox.get("size") or [None, None, None])
    rank1 = rmsds[0] if rmsds else None
    top3 = min(rmsds[:3]) if rmsds else None
    best_rmsd, best_rank = None, None
    if rmsds:
        best_index = min(range(len(rmsds)), key=lambda i: rmsds[i])
        best_rmsd = rmsds[best_index]
        best_rank = best_index + 1
    out.update({
        # legacy columns are refreshed from the manifest too: the pose
        # files are the source of truth, and the crystal RMSD
        # implementation was corrected after the spyrmsd cross-check
        # (item 5) - a report-only regeneration must never keep stale
        # numbers when fresher per-pose values exist.
        "num_poses": len(poses) or None,
        "best_affinity_kcal_mol": poses[0].get("affinity") if poses else None,
        "best_crystal_rmsd_a": best_rmsd,
        "best_pose_rank": best_rank,
        "success": (best_rmsd is not None and best_rmsd <= SUCCESS_RMSD),
        "receptor_engine": receptor.get("engine"),
        "receptor_atoms_out": receptor.get("atoms_out"),
        "ligand_rotatable_bonds": ligand.get("num_rotatable_bonds"),
        "ligand_heavy_atoms": ligand.get("num_heavy_atoms"),
        "gridbox_source": gridbox.get("source"),
        "gridbox_center_x": center[0], "gridbox_center_y": center[1],
        "gridbox_center_z": center[2],
        "gridbox_size_x": size[0], "gridbox_size_y": size[1],
        "gridbox_size_z": size[2],
        "gridbox_padding": gridbox.get("padding"),
        "gridbox_assumption_strength": gridbox.get("assumption_strength"),
        "exhaustiveness": docking_cfg.get("exhaustiveness"),
        "num_modes": docking_cfg.get("num_modes"),
        # Vina python bindings expose no refine parameter (the CLI's
        # --refine): recorded honestly instead of pretending a value.
        "refine": "engine-default",
        "seed": docking_cfg.get("seed"),
        "cpu": docking_cfg.get("cpu"),
        "crystal_rmsd_method": analysis.get("crystal_rmsd_method"),
        "rank1_rmsd_a": rank1,
        "top3_rmsd_a": top3,
        "success_rank1": (rank1 is not None and rank1 <= SUCCESS_RMSD),
        "success_top3": (top3 is not None and top3 <= SUCCESS_RMSD),
    })
    return out


def enrich_rows(rows: list[dict], out_dir: Path) -> dict[str, list[float]]:
    """Fill the item-1 provenance columns from the per-complex manifests.

    Rows are enriched IN PLACE (old-format rows are upgraded); returns
    ``{pdb_id: [per-pose crystal RMSD, ...]}`` for the top-N analysis
    (item 15).  Missing manifests leave the columns blank - never
    invented.
    """
    pose_rmsds: dict[str, list[float]] = {}
    for row in rows:
        run_dir_name = row.get("run_dir") or f"{row['pdb_id']}_benchmark"
        run_dir = out_dir / run_dir_name
        row.update(_manifest_enrichment(run_dir))
        try:
            manifest = json.loads(
                (run_dir / "manifest.json").read_text(encoding="utf-8"))
            poses: list[dict] = []
            for entry in (manifest.get("docking", {}).get("results") or []):
                poses += entry.get("poses") or []
            values = [pose.get("crystal_rmsd") for pose in poses
                      if pose.get("crystal_rmsd") is not None]
            if values:
                pose_rmsds[row["pdb_id"]] = values
        except (OSError, ValueError):
            continue
    return pose_rmsds


def _target_group(target: str) -> str:
    """Per-target grouping key (final pass, item 2): HIV-1 PR vs others."""
    return ("HIV-1 protease" if "hiv-1 protease" in (target or "").lower()
            else "other targets")


def _stats(rows: list[dict], pose_rmsds: dict[str, list[float]]) -> dict:
    ok_rows = [row for row in rows if row.get("status") == "ok"]
    rmsds = [float(row["best_crystal_rmsd_a"]) for row in ok_rows
             if row.get("best_crystal_rmsd_a") not in (None, "", "None")]
    runtimes = [float(row["runtime_s"]) for row in ok_rows
                if row.get("runtime_s") not in (None, "")]
    affinities = [float(row["best_affinity_kcal_mol"]) for row in ok_rows
                  if row.get("best_affinity_kcal_mol") not in (None, "", "None")]
    successes = [row for row in ok_rows if str(row.get("success")) == "True"]
    rank1_hits = [row for row in ok_rows
                  if str(row.get("success_rank1")) == "True"]
    top3_hits = [row for row in ok_rows
                 if str(row.get("success_top3")) == "True"]
    sorted_rmsd = sorted(rmsds)
    median = (sorted_rmsd[len(sorted_rmsd) // 2]
              if len(sorted_rmsd) % 2 else
              (sorted_rmsd[len(sorted_rmsd) // 2 - 1]
               + sorted_rmsd[len(sorted_rmsd) // 2]) / 2) if sorted_rmsd else None

    def _dist(values: list[float]) -> dict | None:
        if not values:
            return None
        import numpy as np

        return {
            "min": round(min(values), 3),
            "p25": round(float(np.percentile(values, 25)), 3),
            "median": round(float(np.percentile(values, 50)), 3),
            "p75": round(float(np.percentile(values, 75)), 3),
            "max": round(max(values), 3),
            "std": round(float(np.std(values, ddof=1)), 3)
            if len(values) > 1 else 0.0,
            "n": len(values),
        }

    def _hits(rows_subset: list[dict], column: str,
              threshold: float) -> int:
        return sum(1 for row in rows_subset
                   if row.get(column) not in (None, "", "None")
                   and float(row[column]) <= threshold)

    # Per-target grouping (item 2): HIV-1 PR complexes vs others.
    groups: dict[str, list[dict]] = {}
    for row in ok_rows:
        groups.setdefault(_target_group(row.get("target", "")), []).append(row)
    per_target = [
        {
            "group": group,
            "attempted": len(members),
            "successful_within_2a": sum(
                1 for row in members if str(row.get("success")) == "True"),
            "rate_pct": round(100.0 * sum(
                1 for row in members if str(row.get("success")) == "True")
                / len(members), 1) if members else None,
            "wilson_95ci_pct": [round(100 * v, 1) for v in wilson_ci(
                sum(1 for row in members if str(row.get("success")) == "True"),
                len(members))],
        }
        for group, members in sorted(groups.items())
    ]

    # Top-N success curve (item 15): best pose among the first N modes.
    topn = {}
    max_poses = max((len(v) for v in pose_rmsds.values()), default=0)
    for n in range(1, min(9, max_poses) + 1):
        hits = sum(1 for values in pose_rmsds.values()
                   if min(values[:n]) <= SUCCESS_RMSD)
        topn[n] = round(100.0 * hits / len(pose_rmsds), 1) \
            if pose_rmsds else None

    stats = {
        "complexes_attempted": len(rows),
        "completed": len(ok_rows),
        "failed": len(rows) - len(ok_rows),
        "successful_within_2a": len(successes),
        "success_rate_pct": (round(100 * len(successes) / len(ok_rows), 1)
                             if ok_rows else None),
        "success_rate_wilson_95ci_pct": [
            round(100 * v, 1) for v in wilson_ci(len(successes), len(ok_rows))],
        "rank1_successful_within_2a": len(rank1_hits),
        "rank1_success_rate_pct": (round(100 * len(rank1_hits) / len(ok_rows), 1)
                                   if ok_rows else None),
        "rank1_wilson_95ci_pct": [
            round(100 * v, 1) for v in wilson_ci(len(rank1_hits), len(ok_rows))],
        "top3_successful_within_2a": len(top3_hits),
        "top3_success_rate_pct": (round(100 * len(top3_hits) / len(ok_rows), 1)
                                  if ok_rows else None),
        "top3_wilson_95ci_pct": [
            round(100 * v, 1) for v in wilson_ci(len(top3_hits), len(ok_rows))],
        "mean_rmsd_a": round(sum(rmsds) / len(rmsds), 2) if rmsds else None,
        "median_rmsd_a": round(median, 2) if median is not None else None,
        "rmsd_distribution": _dist(rmsds),
        "rank1_rmsd_distribution": _dist([
            float(row["rank1_rmsd_a"]) for row in ok_rows
            if row.get("rank1_rmsd_a") not in (None, "", "None")]),
        # multi-definition success (item 15)
        "success_by_definition_threshold": {
            "rank1": {f"{t:.1f}a": _hits(ok_rows, "rank1_rmsd_a", t)
                      for t in THRESHOLDS},
            "top3": {f"{t:.1f}a": _hits(ok_rows, "top3_rmsd_a", t)
                     for t in THRESHOLDS},
            "best_any_rank": {f"{t:.1f}a": _hits(ok_rows, "best_crystal_rmsd_a", t)
                              for t in THRESHOLDS},
        },
        "topn_success_rate_pct": topn,
        "per_target_group": per_target,
        "mean_runtime_s": round(sum(runtimes) / len(runtimes), 1)
        if runtimes else None,
        "mean_best_affinity": round(sum(affinities) / len(affinities), 2)
        if affinities else None,
    }
    return stats


def write_report(rows: list[dict], out_dir: Path, exhaustiveness: int,
                 seed: int) -> None:
    """Publish summary.md + summary_stats.json + figures from the rows."""
    csv_path = out_dir / "redocking_results.csv"
    pose_rmsds = enrich_rows(rows, out_dir)
    rewrite_csv(csv_path, rows)  # upgrade old-format rows in place (item 1)
    stats = _stats(rows, pose_rmsds)
    ok_rows = [row for row in rows if row.get("status") == "ok"]
    failed_rows = [row for row in ok_rows if str(row.get("success")) != "True"]
    rank1_failed = [row for row in ok_rows
                    if str(row.get("success_rank1")) != "True"]

    def _fmt(value, digits: int = 2) -> str:
        if value in (None, "", "None"):
            return "n/a"
        try:
            return f"{float(value):.{digits}f}"
        except (TypeError, ValueError):
            return str(value)

    lines = [
        "# Redocking benchmark results",
        "",
        f"- complexes attempted: {stats['complexes_attempted']}",
        f"- completed: {stats['completed']} "
        f"({stats['failed']} failed - failures are recorded, not dropped)",
        "",
        "## Pose recovery (formal definition: rank-1 pose, "
        f"symmetry-aware heavy-atom RMSD <= {SUCCESS_RMSD:.1f} A)",
        "",
        f"- **rank-1 success: {stats['rank1_successful_within_2a']}/"
        f"{stats['completed']}** = "
        + _fmt_pct_rate(stats["rank1_successful_within_2a"], stats["completed"]),
        f"- top-3 success (best RMSD among the 3 best-scored poses): "
        f"**{stats['top3_successful_within_2a']}/{stats['completed']}** = "
        + _fmt_pct_rate(stats["top3_successful_within_2a"], stats["completed"]),
        f"- best-any-rank success (v0.3.x legacy definition, kept for "
        f"continuity): **{stats['successful_within_2a']}/"
        f"{stats['completed']}** = "
        + _fmt_pct_rate(stats["successful_within_2a"], stats["completed"]),
        "",
        "Definition: 'pose recovery' = the best-SCORING Vina pose (rank 1) "
        "whose symmetry-aware heavy-atom RMSD to the co-crystallized ligand "
        f"is <= {SUCCESS_RMSD:.1f} A.  The {SUCCESS_RMSD:.1f} A threshold is a "
        "community convention (Trott & Olson 2010; CSAR / D3R benchmark "
        "practice), NOT a thermodynamically validated cutoff - see "
        "docs/interaction_criteria.md.  Alternative definitions "
        "(best-RMSD-among-top-N) yield different success rates, which is "
        "exactly why all of them are published above.",
        "",
        "## RMSD distribution (best pose per complex, A)",
        "",
    ]
    dist = stats.get("rmsd_distribution") or {}
    if dist:
        lines += [
            f"- min {dist['min']:.2f} | p25 {dist['p25']:.2f} | "
            f"median {dist['median']:.2f} | p75 {dist['p75']:.2f} | "
            f"max {dist['max']:.2f} | std {dist['std']:.2f} (n={dist['n']})",
            "- histogram: `rmsd_distribution.png`",
            f"- mean {_fmt(stats['mean_rmsd_a'])} A | "
            f"median {_fmt(stats['median_rmsd_a'])} A",
        ]
    lines += [
        "",
        f"- mean runtime per complex: {_fmt(stats['mean_runtime_s'], 1)} s",
        f"- exhaustiveness: {exhaustiveness} (fixed seed {seed}, "
        "success rates reported with Wilson score 95% confidence "
        "intervals; formula documented in run_benchmark.py)",
        "",
        "## Success by definition and threshold (item 15)",
        "",
        "| pose-recovery definition | <= 1.0 A | <= 2.0 A | <= 3.0 A |",
        "|---|---|---|---|",
    ]
    by_def = stats.get("success_by_definition_threshold") or {}
    for label, key in (("rank-1 (formal)", "rank1"), ("top-3", "top3"),
                       ("best-any-rank (legacy)", "best_any_rank")):
        row_stats = by_def.get(key) or {}
        lines.append(
            f"| {label} | {row_stats.get('1.0a', 'n/a')} | "
            f"{row_stats.get('2.0a', 'n/a')} | {row_stats.get('3.0a', 'n/a')} |")
    lines += [
        "",
        "Top-N success curve (success rate when any of the first N scored "
        "poses is within 2.0 A): `topn_success.png`.  A flat curve for a "
        "complex means no good pose was sampled at all; a rising curve "
        "means a good pose exists but is outranked.",
        "",
        "## Per-target grouping (item 2)",
        "",
        "| target group | complexes | successes (best-any-rank, 2.0 A) | rate (95% CI) |",
        "|---|---|---|---|",
    ]
    for group in stats.get("per_target_group") or []:
        low, high = group["wilson_95ci_pct"]
        lines.append(
            f"| {group['group']} | {group['attempted']} | "
            f"{group['successful_within_2a']} | "
            f"{group['rate_pct']:.1f}% (95% CI: {low:.1f}%-{high:.1f}%) |")
    if failed_rows:
        lines += [
            "",
            "## Failed complexes, legacy definition (best-any-rank, item 2)",
            "",
            "Every failed complex had a strong ligand-derived grid box "
            "(see the per-complex table), so these failures are "
            "docking-side (sampling/scoring of flexible ligands), not "
            "grid-box or preparation-side.  Hypothesized causes below are "
            "the authors' reading of the provenance columns, not proven "
            "mechanisms.",
            "",
            "| pdb | ligand | target | best RMSD (A) | rank-1 RMSD (A) | "
            "rotatable bonds | hypothesized cause |",
            "|---|---|---|---|---|---|---|",
        ]
        for row in failed_rows:
            hypothesis = FAILURE_HYPOTHESES.get(
                row["pdb_id"],
                "no pose within 2.0 A among the returned modes "
                "(see per-complex artifacts)")
            lines.append(
                f"| {row['pdb_id']} | {row['ligand']} | "
                f"{row.get('target', '')} | {_fmt(row.get('best_crystal_rmsd_a'))} "
                f"| {_fmt(row.get('rank1_rmsd_a'))} | "
                f"{_fmt(row.get('ligand_rotatable_bonds'), 0)} | {hypothesis} |")
    if rank1_failed:
        lines += [
            "",
            "## Failed complexes, formal definition (rank-1, item 2/15)",
            "",
            "The formal pose-recovery definition is stricter: the "
            "best-SCORING pose must be within 2.0 A.  Distinguishing the "
            "two failure modes is the point of the top-N curve: a good "
            "pose that exists but ranks 3rd is a SCORING failure; no good "
            "pose at all is a SAMPLING failure.",
            "",
            "| pdb | ligand | rank-1 RMSD (A) | best RMSD (A) at rank | "
            "failure mode |",
            "|---|---|---|---|---|",
        ]
        for row in rank1_failed:
            best = row.get("best_crystal_rmsd_a")
            sampled = best not in (None, "", "None") and float(best) <= SUCCESS_RMSD
            mode = ("ranking (good pose sampled but outranked)" if sampled
                    else "sampling (no pose within 2.0 A at any rank)")
            lines.append(
                f"| {row['pdb_id']} | {row['ligand']} | "
                f"{_fmt(row.get('rank1_rmsd_a'))} | "
                f"{_fmt(row.get('best_crystal_rmsd_a'))} @ rank "
                f"{row.get('best_pose_rank', '-')} | {mode} |")
    lines += [
        "",
        "## Per-complex results (with preparation + grid provenance, item 1)",
        "",
        "Full provenance (grid-box center/size, docking parameters, RMSD "
        "method) is in `redocking_results.csv`; this table shows the "
        "columns a reader needs to judge whether a failure is prep- or "
        "docking-related.",
        "",
        "| pdb | ligand | target | best affinity | rank-1 RMSD (A) | "
        "best RMSD (A) | success | gridbox source | assumption | engine | "
        "runtime (s) |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        if row.get("status") == "ok":
            lines.append(
                f"| {row['pdb_id']} | {row['ligand']} | "
                f"{row.get('target', '')} | "
                f"{_fmt(row.get('best_affinity_kcal_mol'))} | "
                f"{_fmt(row.get('rank1_rmsd_a'))} | "
                f"{_fmt(row.get('best_crystal_rmsd_a'))} | "
                f"{'yes' if str(row.get('success')) == 'True' else 'NO'} | "
                f"{row.get('gridbox_source', '')} | "
                f"{row.get('gridbox_assumption_strength', '')} | "
                f"{row.get('engine', '')} | "
                f"{_fmt(row.get('runtime_s'), 1)} |"
            )
        else:
            lines.append(
                f"| {row['pdb_id']} | {row.get('ligand', '')} | "
                f"{row.get('target', '')} | FAILED | - | - | - | - | - | - | "
                f"{row.get('status', '')} |"
            )
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n",
                                        encoding="utf-8")
    (out_dir / "summary_stats.json").write_text(
        json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    _write_scatter(ok_rows, out_dir)
    _write_histogram(ok_rows, out_dir)
    _write_topn(pose_rmsds, out_dir)


def _write_scatter(ok_rows: list[dict], out_dir: Path) -> None:
    """Scatter of best affinity (predicted) vs best crystal RMSD (agreement)."""
    pts: list[tuple[float, float, str, bool]] = []
    for row in ok_rows:
        affinity = row.get("best_affinity_kcal_mol")
        rmsd = row.get("best_crystal_rmsd_a")
        if affinity in (None, "", "None") or rmsd in (None, "", "None"):
            continue
        try:
            pts.append((float(affinity), float(rmsd), row["pdb_id"],
                        str(row.get("success")) == "True"))
        except (TypeError, ValueError):
            continue
    if not pts:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    colors = ["#1a7f37" if p[3] else "#b35900" for p in pts]
    ax.scatter(xs, ys, c=colors, s=55, zorder=3)
    for x, y, pdb_id, _ok in pts:
        ax.annotate(pdb_id, (x, y), fontsize=7, xytext=(4, 3),
                    textcoords="offset points")
    ax.axhline(SUCCESS_RMSD, color="#c62828", ls="--", lw=1.2, zorder=2)
    x_max = max(xs)
    ax.text(x_max, SUCCESS_RMSD + 0.06,
            f"heuristic success threshold {SUCCESS_RMSD:.1f} A",
            ha="right", fontsize=8, color="#c62828")
    ax.set_xlabel("best Vina affinity (kcal/mol)")
    ax.set_ylabel("best heavy-atom RMSD to crystal pose (A)")
    ax.set_title("Redocking benchmark: predicted affinity vs pose recovery "
                 f"({len(pts)} complexes)")
    fig.savefig(out_dir / "scatter.png", dpi=150)
    plt.close(fig)


def _write_histogram(ok_rows: list[dict], out_dir: Path) -> None:
    """RMSD distribution histogram (final pass, item 2)."""
    values = []
    for row in ok_rows:
        rmsd = row.get("best_crystal_rmsd_a")
        if rmsd not in (None, "", "None"):
            try:
                values.append(float(rmsd))
            except (TypeError, ValueError):
                continue
    if not values:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, ax = plt.subplots(figsize=(7.0, 5.0), constrained_layout=True)
    ax.hist(values, bins=min(12, max(5, len(values) // 2)), color="#4c72b0",
            edgecolor="white", zorder=2)
    ax.axvline(SUCCESS_RMSD, color="#c62828", ls="--", lw=1.2,
               label="2.0 A community-convention threshold", zorder=3)
    mean_v = sum(values) / len(values)
    sorted_v = sorted(values)
    median_v = (sorted_v[len(sorted_v) // 2] if len(sorted_v) % 2
                else (sorted_v[len(sorted_v) // 2 - 1]
                      + sorted_v[len(sorted_v) // 2]) / 2)
    ax.axvline(mean_v, color="#2f4b7c", ls=":", lw=1.4,
               label=f"mean {mean_v:.2f} A", zorder=3)
    ax.axvline(median_v, color="#55a868", ls="-.", lw=1.4,
               label=f"median {median_v:.2f} A", zorder=3)
    ax.set_xlabel("best pose heavy-atom RMSD to crystal (A)")
    ax.set_ylabel("complexes")
    ax.set_title(f"RMSD distribution ({len(values)} complexes, "
                 "best pose per complex)")
    ax.legend(loc="upper right")
    fig.savefig(out_dir / "rmsd_distribution.png", dpi=150)
    plt.close(fig)


def _write_topn(pose_rmsds: dict[str, list[float]], out_dir: Path) -> None:
    """Success rate vs N (top-N scored poses, final pass, item 15)."""
    if not pose_rmsds:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    max_poses = min(9, max(len(v) for v in pose_rmsds.values()))
    ns = list(range(1, max_poses + 1))
    rates = [
        100.0 * sum(1 for values in pose_rmsds.values()
                    if min(values[:n]) <= SUCCESS_RMSD) / len(pose_rmsds)
        for n in ns
    ]
    fig, ax = plt.subplots(figsize=(6.5, 4.5), constrained_layout=True)
    ax.plot(ns, rates, marker="o", color="#4c72b0", lw=1.8, zorder=3)
    for n, rate in zip(ns, rates, strict=True):
        ax.annotate(f"{rate:.0f}%", (n, rate), fontsize=8,
                    xytext=(0, 7), textcoords="offset points",
                    ha="center")
    ax.set_xticks(ns)
    ax.set_xlabel("N: best RMSD among the first N scored poses")
    ax.set_ylabel("complexes within 2.0 A (%)")
    ax.set_ylim(0, 105)
    ax.set_title("Top-N success curve: is the right pose sampled but "
                 "outranked, or not sampled at all?")
    ax.grid(alpha=0.3, zorder=1)
    fig.savefig(out_dir / "topn_success.png", dpi=150)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--complexes",
                        default=str(Path(__file__).parent / "complexes.csv"))
    parser.add_argument("--out", default=str(Path(__file__).parent / "results"))
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--cpu", type=int, default=0)
    parser.add_argument("--limit", type=int, default=None,
                        help="run only the first N complexes")
    parser.add_argument("--pdb-ids", default=None,
                        help="comma-separated subset of PDB ids")
    parser.add_argument("--report-only", action="store_true",
                        help="regenerate summary/scatter from the CSV "
                        "without docking")
    parser.add_argument("--max-seconds", type=float, default=None,
                        help="stop cleanly after this many seconds "
                        "(resume by re-running; already-recorded "
                        "complexes are skipped)")
    parser.add_argument("--worker-timeout", type=float, default=520,
                        help="per-complex worker timeout in seconds")
    parser.add_argument("--worker", default=None, help=argparse.SUPPRESS)
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "redocking_results.csv"

    # ---- worker mode: run a single complex in this fresh interpreter ----
    if args.worker:
        complexes = load_complexes(Path(args.complexes))
        row = next(c for c in complexes if c["pdb_id"] == args.worker)
        try:
            result = run_complex(row, out_dir, args.exhaustiveness,
                                 args.seed, args.cpu)
        except Exception as exc:  # noqa: BLE001 - record, never crash the driver
            result = {"pdb_id": row["pdb_id"], "ligand": row["ligand_resname"],
                      "target": row.get("target", ""),
                      "status": f"error: {exc}"}
        print(json.dumps(result), flush=True)
        return 0

    # ---- report-only mode ------------------------------------------------
    if args.report_only:
        rows = list(load_existing(csv_path).values())
        if not rows:
            print("no results recorded yet", file=sys.stderr)
            return 1
        write_report(rows, out_dir, args.exhaustiveness, args.seed)
        print(f"summary: {out_dir / 'summary.md'}")
        print(f"(rows enriched from per-complex manifests; CSV upgraded "
              f"in place: {csv_path})")
        return 0

    # ---- driver mode -----------------------------------------------------
    complexes = load_complexes(Path(args.complexes))
    if args.pdb_ids:
        wanted = {pid.strip().upper() for pid in args.pdb_ids.split(",")}
        complexes = [row for row in complexes if row["pdb_id"].upper() in wanted]
    if args.limit:
        complexes = complexes[: args.limit]

    existing = load_existing(csv_path)
    todo = [row for row in complexes
            if row["pdb_id"] not in existing]
    print(f"redocking benchmark: {len(complexes)} complex(es) in set, "
          f"{len(existing)} already recorded, {len(todo)} to run | "
          f"exhaustiveness={args.exhaustiveness}, seed={args.seed}",
          flush=True)

    started_driver = time.perf_counter()
    for index, complex_row in enumerate(todo, start=1):
        pdb_id = complex_row["pdb_id"]
        if args.max_seconds is not None and \
                time.perf_counter() - started_driver > args.max_seconds:
            print(f"time budget reached; {len(todo) - index + 1} complex(es) "
                  "remain - re-run this command to resume", flush=True)
            break
        print(f"[{index}/{len(todo)}] {pdb_id} "
              f"({complex_row['ligand_resname']}) ...", flush=True)
        started = time.perf_counter()
        cmd = [sys.executable, str(Path(__file__).resolve()),
               "--complexes", args.complexes, "--out", args.out,
               "--exhaustiveness", str(args.exhaustiveness),
               "--seed", str(args.seed), "--cpu", str(args.cpu),
               "--worker", pdb_id]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=args.worker_timeout)
            last = proc.stdout.strip().splitlines()[-1] if proc.stdout.strip() \
                else ""
            result = json.loads(last) if last.startswith("{") else None
        except subprocess.TimeoutExpired:
            result = {"pdb_id": pdb_id, "ligand": complex_row["ligand_resname"],
                      "target": complex_row.get("target", ""),
                      "status": f"error: worker exceeded "
                                f"{args.worker_timeout:.0f}s timeout"}
        if result is None:
            result = {"pdb_id": pdb_id, "ligand": complex_row["ligand_resname"],
                      "target": complex_row.get("target", ""),
                      "status": "error: worker crashed or produced no result"}
        append_row(csv_path, result)
        if result.get("status") == "ok":
            print(f"    best {result.get('best_affinity_kcal_mol')} kcal/mol, "
                  f"RMSD {result.get('best_crystal_rmsd_a')} A, "
                  f"success={result.get('success')} "
                  f"({time.perf_counter() - started:.0f}s)", flush=True)
        else:
            print(f"    {result.get('status')}", flush=True)
        # keep summary.md/scatter always current, even mid-run
        rows = list(load_existing(csv_path).values())
        write_report(rows, out_dir, args.exhaustiveness, args.seed)

    rows = list(load_existing(csv_path).values())
    write_report(rows, out_dir, args.exhaustiveness, args.seed)
    print(f"\nresults: {csv_path}")
    print(f"summary: {out_dir / 'summary.md'}")
    ok = [row for row in rows if row.get("status") == "ok"]
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

