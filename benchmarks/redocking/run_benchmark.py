"""Automated redocking benchmark runner (work order item 1).

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
* ``summary.md`` + ``summary_stats.json`` + ``scatter.png`` are regenerated
  after every complex, so partial results are always published;
* failed complexes are recorded with their failure reason and are NEVER
  dropped from the CSV or the summary table.

Outputs (in ``--out``):

* ``redocking_results.csv``  - one row per complex (the core artifact)
* ``summary.md``             - success rate, mean/median RMSD, mean runtime,
  per-complex table (including failures with reasons)
* ``summary_stats.json``     - the same numbers, machine-readable
* ``scatter.png``            - best affinity vs best crystal RMSD scatter
* ``<pdb_id>/``              - full pipeline run directories for inspection
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

FIELDNAMES = [
    "pdb_id", "ligand", "target", "status", "num_poses",
    "best_affinity_kcal_mol", "best_crystal_rmsd_a", "best_pose_rank",
    "success", "engine", "docking_backend", "runtime_s", "run_dir",
]


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


def _stats(rows: list[dict]) -> dict:
    ok_rows = [row for row in rows if row.get("status") == "ok"]
    rmsds = [float(row["best_crystal_rmsd_a"]) for row in ok_rows
             if row.get("best_crystal_rmsd_a") not in (None, "", "None")]
    runtimes = [float(row["runtime_s"]) for row in ok_rows
                if row.get("runtime_s") not in (None, "")]
    affinities = [float(row["best_affinity_kcal_mol"]) for row in ok_rows
                  if row.get("best_affinity_kcal_mol") not in (None, "", "None")]
    successes = [row for row in ok_rows if str(row.get("success")) == "True"]
    sorted_rmsd = sorted(rmsds)
    median = (sorted_rmsd[len(sorted_rmsd) // 2]
              if len(sorted_rmsd) % 2 else
              (sorted_rmsd[len(sorted_rmsd) // 2 - 1]
               + sorted_rmsd[len(sorted_rmsd) // 2]) / 2) if sorted_rmsd else None
    return {
        "complexes_attempted": len(rows),
        "completed": len(ok_rows),
        "failed": len(rows) - len(ok_rows),
        "successful_within_2a": len(successes),
        "success_rate_pct": (round(100 * len(successes) / len(ok_rows), 1)
                             if ok_rows else None),
        "mean_rmsd_a": round(sum(rmsds) / len(rmsds), 2) if rmsds else None,
        "median_rmsd_a": round(median, 2) if median is not None else None,
        "mean_runtime_s": round(sum(runtimes) / len(runtimes), 1)
        if runtimes else None,
        "mean_best_affinity": round(sum(affinities) / len(affinities), 2)
        if affinities else None,
    }


def write_report(rows: list[dict], out_dir: Path, exhaustiveness: int,
                 seed: int) -> None:
    """Publish summary.md + summary_stats.json + scatter.png from the rows."""
    stats = _stats(rows)
    ok_rows = [row for row in rows if row.get("status") == "ok"]

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
        f"- successful (best pose RMSD <= {SUCCESS_RMSD:.1f} A, heuristic "
        f"community convention): **{stats['successful_within_2a']}/"
        f"{stats['completed']}**"
        + (f" ({stats['success_rate_pct']:.0f}%)" if ok_rows else ""),
        f"- mean best-pose RMSD: {_fmt(stats['mean_rmsd_a'])} A | "
        f"median: {_fmt(stats['median_rmsd_a'])} A",
        f"- mean runtime per complex: {_fmt(stats['mean_runtime_s'], 1)} s",
        f"- exhaustiveness: {exhaustiveness} (fixed seed {seed})",
        "",
        "## Per-complex results",
        "",
        "| pdb | ligand | target | best affinity | best RMSD (A) | "
        "success | runtime (s) | engine |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        if row.get("status") == "ok":
            lines.append(
                f"| {row['pdb_id']} | {row['ligand']} | "
                f"{row.get('target', '')} | "
                f"{_fmt(row.get('best_affinity_kcal_mol'))} | "
                f"{_fmt(row.get('best_crystal_rmsd_a'))} | "
                f"{'yes' if str(row.get('success')) == 'True' else 'NO'} | "
                f"{_fmt(row.get('runtime_s'), 1)} | "
                f"{row.get('engine', '')} |"
            )
        else:
            lines.append(
                f"| {row['pdb_id']} | {row.get('ligand', '')} | "
                f"{row.get('target', '')} | FAILED | - | - | - | "
                f"{row.get('status', '')} |"
            )
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n",
                                        encoding="utf-8")
    (out_dir / "summary_stats.json").write_text(
        json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    _write_scatter(ok_rows, out_dir)


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

    # ---- report-only mode ----
    if args.report_only:
        rows = list(load_existing(csv_path).values())
        if not rows:
            print("no results recorded yet", file=sys.stderr)
            return 1
        write_report(rows, out_dir, args.exhaustiveness, args.seed)
        print(f"summary: {out_dir / 'summary.md'}")
        return 0

    # ---- driver mode ----
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
