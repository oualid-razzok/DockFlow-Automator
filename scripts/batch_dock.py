#!/usr/bin/env python3
"""Batch docking for ligand libraries (standalone utility).

Accepts a directory of prepared PDBQT ligands, a multi-record SDF, or a CSV
with a ``smiles`` column, docks every ligand against one receptor with
AutoDock Vina (CLI backend, parallelised across processes), and writes a
ranked CSV summary.

Examples::

    # directory of prepared ligands
    python scripts/batch_dock.py \\
        --receptor prepared/receptor.pdbqt \\
        --ligands "prepared/ligands/*.pdbqt" \\
        --center 12.3,-4.5,21.7 --size 22,24,20 \\
        --exhaustiveness 16 --parallel 4 --out-dir runs/batch

    # CSV library (column 'smiles', optional 'id')
    python scripts/batch_dock.py --receptor r.pdbqt --csv library.csv \\
        --center 0,0,0 --size 24,24,24 --out-dir runs/batch
"""

from __future__ import annotations

import argparse
import glob
import multiprocessing as mp
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from dockflow_core.docker_engine import VinaConfig, VinaEngine  # noqa: E402
from dockflow_core.utils import setup_logging  # noqa: E402


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch-dock a ligand library against one receptor.")
    parser.add_argument("--receptor", required=True, help="receptor PDBQT")
    parser.add_argument("--ligands", default=None,
                        help="glob of ligand PDBQT files (quote it!)")
    parser.add_argument("--csv", default=None,
                        help="CSV with a 'smiles' column (prepared with Meeko)")
    parser.add_argument("--sdf", default=None, help="multi-record SDF library")
    parser.add_argument("--center", required=True, help="grid center x,y,z")
    parser.add_argument("--size", default="24,24,24", help="grid size x,y,z")
    parser.add_argument("--out-dir", default="batch_docking")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--cpu", type=int, default=1,
                        help="threads per docking job")
    parser.add_argument("--parallel", type=int, default=None,
                        help="concurrent jobs (default: cpu count // 2)")
    parser.add_argument("--vina", default="vina", help="vina executable")
    parser.add_argument("--smina", default=None, help="use smina executable")
    parser.add_argument("--prepare-csv", action="store_true",
                        help="prepare SMILES through Meeko before docking")
    return parser.parse_args(argv)


def _triple(text: str) -> tuple[float, float, float]:
    parts = [p.strip() for p in text.split(",")]
    if len(parts) != 3:
        raise SystemExit(f"expected x,y,z but got {text!r}")
    return tuple(float(p) for p in parts)  # type: ignore[return-value]


def collect_ligands(args: argparse.Namespace, out_dir: Path) -> list[Path]:
    files: list[Path] = []
    if args.ligands:
        matches = sorted(glob.glob(args.ligands))
        if not matches:
            raise SystemExit(f"no files match {args.ligands}")
        files.extend(Path(m) for m in matches)
    if args.sdf:
        files.extend(_prepare_sdf(Path(args.sdf), out_dir))
    if args.csv:
        files.extend(_prepare_csv(Path(args.csv), out_dir, args.prepare_csv))
    if not files:
        raise SystemExit("no ligands: pass --ligands, --sdf or --csv")
    return files


def _prepare_sdf(sdf_path: Path, out_dir: Path) -> list[Path]:
    from dockflow_core.preparator import LigandPreparator

    results = LigandPreparator().prepare_library(sdf_path, out_dir / "prepared")
    return [r.pdbqt_path for r in results if r.ok and r.pdbqt_path]


def _prepare_csv(csv_path: Path, out_dir: Path, prepare: bool) -> list[Path]:
    import csv

    rows: list[dict] = []
    with open(csv_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if row.get("smiles"):
                rows.append(row)
    if not rows:
        raise SystemExit(f"no 'smiles' column in {csv_path}")
    if not prepare:
        raise SystemExit(
            "CSV input needs --prepare-csv to convert SMILES into PDBQT first"
        )
    smi_path = out_dir / "library.smi"
    with open(smi_path, "w", encoding="utf-8") as fh:
        for index, row in enumerate(rows):
            name = row.get("id") or row.get("name") or f"ligand_{index + 1}"
            fh.write(f"{row['smiles']} {name}\n")
    return _prepare_sdf(smi_path, out_dir)


def _dock_one(job: dict) -> dict:
    """Worker function executed inside a child process (CLI backend)."""
    config = VinaConfig(
        center=tuple(job["center"]),
        size=tuple(job["size"]),
        exhaustiveness=job["exhaustiveness"],
        num_modes=job["num_modes"],
        seed=job["seed"],
        cpu=job["cpu"],
        timeout=job["timeout"],
    )
    engine = VinaEngine(
        config, backend=job["backend"], workdir=job["out_dir"],
        vina_exec=job["vina_exec"], smina_exec=job["smina_exec"],
    )
    result = engine.dock(job["receptor"], job["ligand"])
    return {
        "ligand": result.ligand_name,
        "best": result.best_affinity,
        "poses": len(result.poses),
        "runtime": result.runtime,
        "error": result.error,
        "out": str(result.out_path) if result.out_path else None,
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    setup_logging("INFO")
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ligand_files = collect_ligands(args, out_dir)
    print(f"receptor : {args.receptor}")
    print(f"ligands  : {len(ligand_files)}")
    print(f"grid     : center {_triple(args.center)} size {_triple(args.size)}")

    backend = "smina" if args.smina else "cli"
    vina_exec = args.smina or args.vina
    jobs = [
        {
            "center": _triple(args.center),
            "size": _triple(args.size),
            "exhaustiveness": args.exhaustiveness,
            "num_modes": args.num_modes,
            "seed": args.seed,
            "cpu": max(1, args.cpu),
            "timeout": 3600.0,
            "backend": backend,
            "vina_exec": vina_exec,
            "smina_exec": args.smina,
            "out_dir": str(out_dir),
            "receptor": args.receptor,
            "ligand": str(path),
        }
        for path in ligand_files
    ]
    workers = args.parallel or max(1, (mp.cpu_count() or 2) // 2)
    workers = min(workers, len(jobs))
    print(f"workers  : {workers} (backend={backend})")

    started = time.perf_counter()
    results: list[dict] = []
    if workers == 1:
        for job in jobs:
            result = _dock_one(job)
            results.append(result)
            _report(result)
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(workers) as pool:
            for result in pool.imap_unordered(_dock_one, jobs):
                results.append(result)
                _report(result)

    ranked = sorted(
        [r for r in results if r["best"] is not None],
        key=lambda r: r["best"],
    ) + [r for r in results if r["best"] is None]
    summary = out_dir / "batch_summary.csv"
    _write_csv(ranked, summary)
    print(f"\nfinished in {time.perf_counter() - started:.1f}s")
    print(f"summary  : {summary}")
    print("\ntop 10 ligands:")
    for index, result in enumerate(ranked[:10], start=1):
        if result["best"] is not None:
            print(f"  {index:3d}. {result['ligand']:<30s} {result['best']:8.2f} kcal/mol")
    return 0 if any(r["best"] is not None for r in results) else 1


def _report(result: dict) -> None:
    if result["error"]:
        print(f"  FAILED {result['ligand']}: {result['error']}")
    elif result["best"] is not None:
        print(f"  {result['ligand']}: {result['best']:.2f} kcal/mol "
              f"({result['poses']} poses, {result['runtime']:.0f}s)")


def _write_csv(results: list[dict], path: Path) -> None:
    import csv

    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["rank", "ligand", "best_affinity", "num_poses",
                         "runtime_s", "error", "out_pdbqt"])
        for index, result in enumerate(results, start=1):
            writer.writerow([
                index, result["ligand"],
                f"{result['best']:.3f}" if result["best"] is not None else "",
                result["poses"], f"{result['runtime']:.1f}",
                result["error"] or "", result["out"] or "",
            ])


if __name__ == "__main__":
    sys.exit(main())
