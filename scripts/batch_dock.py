#!/usr/bin/env python3
"""Batch docking for ligand libraries (standalone utility).

Accepts a directory of prepared PDBQT ligands, a multi-record SDF, or a CSV
with a ``smiles`` column, docks every ligand against one receptor with
AutoDock Vina (CLI backend, parallelised across processes), and writes a
ranked CSV summary.

Large-library support (work-order item 23):

* **chunked streaming** - SDF/SMILES libraries are prepared and docked
  in chunks (``--chunk-size``, default 500), so memory stays constant
  for 100k+ record libraries (records are never all in RAM);
* **checkpoint resume** - every finished ligand is appended to
  ``batch_results.csv`` the moment its job returns; re-running the same
  command skips ligands already recorded, so an interrupted screen
  continues where it stopped (``--no-resume`` to start over);
* the ranked ``batch_summary.csv`` is rebuilt from the incremental
  results file at the end.

Examples::

    # directory of prepared ligands
    python scripts/batch_dock.py \\
        --receptor prepared/receptor.pdbqt \\
        --ligands "prepared/ligands/*.pdbqt" \\
        --center 12.3,-4.5,21.7 --size 22,24,20 \\
        --exhaustiveness 16 --parallel 4 --out-dir runs/batch

    # 100k-record SDF library, chunked
    python scripts/batch_dock.py --receptor r.pdbqt \\
        --sdf library.sdf --center 0,0,0 --size 24,24,24 \\
        --chunk-size 500 --out-dir runs/screen

    # CSV library (column 'smiles', optional 'id')
    python scripts/batch_dock.py --receptor r.pdbqt --csv library.csv \\
        --center 0,0,0 --size 24,24,24 --out-dir runs/batch --prepare-csv
"""

from __future__ import annotations

import argparse
import csv
import glob
import multiprocessing as mp
import sys
import time
from collections.abc import Iterator
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from dockflow_core.docker_engine import VinaConfig, VinaEngine  # noqa: E402
from dockflow_core.utils import setup_logging  # noqa: E402

RESULT_FIELDS = ["ligand", "best_affinity", "num_poses", "runtime_s",
                 "error", "out_pdbqt"]


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
    parser.add_argument("--backend", default="cli", choices=["cli", "python"],
                        help="docking backend: vina CLI subprocesses or the "
                             "python bindings (per worker process)")
    parser.add_argument("--prepare-csv", action="store_true",
                        help="prepare SMILES through Meeko before docking")
    parser.add_argument("--chunk-size", type=int, default=500,
                        help="records prepared + docked per chunk "
                             "(SDF/SMILES libraries; constant memory)")
    parser.add_argument("--max-records", type=int, default=None,
                        help="stop after N records (smoke tests / sampling)")
    parser.add_argument("--no-resume", action="store_true",
                        help="ignore previously recorded results and "
                             "re-dock everything")
    parser.add_argument("--skip-first", type=int, default=0,
                        help="skip the first N ligands of the source list "
                             "(slurm array task slicing)")
    parser.add_argument("--take", type=int, default=None,
                        help="process at most N ligands after the skip "
                             "(slurm array task slicing)")
    return parser.parse_args(argv)


def _triple(text: str) -> tuple[float, float, float]:
    parts = [p.strip() for p in text.split(",")]
    if len(parts) != 3:
        raise SystemExit(f"expected x,y,z but got {text!r}")
    return tuple(float(p) for p in parts)  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# chunked library preparation (item 23: constant memory for 100k+ records)
# ---------------------------------------------------------------------------
def iter_ligand_batches(args: argparse.Namespace, out_dir: Path,
                        chunk_size: int) -> Iterator[list[Path]]:
    """Yield ligand PDBQT paths in bounded chunks.

    PDBQT globs stream from disk; SDF/SMILES sources are prepared chunk
    by chunk so at most ``chunk_size`` molecules are in memory at once.
    """
    if args.ligands:
        matches = sorted(glob.glob(args.ligands))
        if not matches:
            raise SystemExit(f"no files match {args.ligands}")
        batch: list[Path] = []
        for match in matches:
            batch.append(Path(match))
            if len(batch) >= chunk_size:
                yield batch
                batch = []
        if batch:
            yield batch
        return

    if args.sdf:
        yield from _iter_prepared_chunks(Path(args.sdf), out_dir, chunk_size,
                                          args.max_records)
        return

    if args.csv:
        rows: list[dict] = []
        with open(args.csv, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                if row.get("smiles"):
                    rows.append(row)
        if not rows:
            raise SystemExit(f"no 'smiles' column in {args.csv}")
        if not args.prepare_csv:
            raise SystemExit(
                "CSV input needs --prepare-csv to convert SMILES into PDBQT first"
            )
        smi_path = out_dir / "library.smi"
        with open(smi_path, "w", encoding="utf-8") as fh:
            for index, row in enumerate(rows):
                name = row.get("id") or row.get("name") or f"ligand_{index + 1}"
                fh.write(f"{row['smiles']} {name}\n")
        yield from _iter_prepared_chunks(smi_path, out_dir, chunk_size,
                                          args.max_records)
        return

    raise SystemExit("no ligands: pass --ligands, --sdf or --csv")


def _iter_prepared_chunks(source: Path, out_dir: Path, chunk_size: int,
                          max_records: int | None) -> Iterator[list[Path]]:
    """Prepare an SDF/SMI library chunk by chunk (bounded memory)."""
    from rdkit import Chem

    preparator = _ligand_preparator()
    prepared_dir = out_dir / "prepared"
    supplier = Chem.ForwardSDMolSupplier(str(source), removeHs=False,
                                         sanitize=True) \
        if source.suffix.lower() == ".sdf" else \
        Chem.SmilesMolSupplier(str(source), titleLine=True)
    batch: list = []
    seen = 0
    for index, mol in enumerate(supplier):
        if max_records is not None and seen >= max_records:
            break
        if mol is None:
            print(f"  skipping invalid record {index + 1}")
            continue
        seen += 1
        batch.append(mol)
        if len(batch) >= chunk_size:
            yield _prepare_batch(preparator, batch, prepared_dir)
            batch = []
    if batch:
        yield _prepare_batch(preparator, batch, prepared_dir)


def _prepare_batch(preparator, mols: list, prepared_dir: Path) -> list[Path]:
    """Prepare one in-memory chunk of molecules to PDBQT (item 23)."""
    from dockflow_core.preparator import LigandPrepResult
    from dockflow_core.utils import ensure_dir

    prepared_dir = ensure_dir(prepared_dir)
    paths: list[Path] = []
    for index, mol in enumerate(mols):
        name = (mol.GetProp("_Name") if mol.HasProp("_Name") else "") or ""
        name = name.strip().replace(" ", "_")[:60] or f"ligand_{index + 1}"
        result = LigandPrepResult(identifier=name, engine="meeko")
        try:
            pdbqt_path, _ = preparator._prepare_variant(
                mol, prepared_dir, name, result)
            result.pdbqt_path = pdbqt_path
        except Exception as exc:  # noqa: BLE001 - one bad record never stops a screen
            result.error = str(exc)
        if result.pdbqt_path is not None:
            paths.append(Path(result.pdbqt_path))
        else:
            print(f"  prep failed for {name}: {result.error}")
    return paths


def _ligand_preparator():
    from dockflow_core.preparator import LigandPreparator

    return LigandPreparator()


# ---------------------------------------------------------------------------
# checkpoint / results bookkeeping
# ---------------------------------------------------------------------------
def load_recorded(results_csv: Path) -> dict[str, dict]:
    """Rows already recorded by a previous (interrupted) run, keyed by name."""
    if not results_csv.is_file():
        return {}
    with open(results_csv, newline="", encoding="utf-8") as fh:
        return {row["ligand"]: row for row in csv.DictReader(fh)
                if row.get("ligand")}


def append_rows(results_csv: Path, rows: list[dict]) -> None:
    """Append finished ligands immediately (crash-safe incremental flush)."""
    is_new = not results_csv.exists()
    with open(results_csv, "a", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=RESULT_FIELDS)
        if is_new:
            writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in RESULT_FIELDS})


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
    results_csv = out_dir / "batch_results.csv"
    recorded = {} if args.no_resume else load_recorded(results_csv)
    if recorded:
        print(f"resume  : {len(recorded)} ligand(s) already recorded in "
              f"{results_csv.name} (use --no-resume to re-dock)")

    print(f"receptor : {args.receptor}")
    print(f"grid     : center {_triple(args.center)} size {_triple(args.size)}")

    backend = args.smina or args.backend
    vina_exec = args.smina or args.vina

    def _job(path: Path) -> dict:
        return {
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

    workers = args.parallel or max(1, (mp.cpu_count() or 2) // 2)
    print(f"workers  : {workers} (backend={backend}, "
          f"chunk size {args.chunk_size})")

    started = time.perf_counter()
    total_docked = total_skipped = 0
    slice_seen = 0  # ligands consumed for --skip-first/--take slicing
    for chunk_index, batch in enumerate(iter_ligand_batches(args, out_dir,
                                                            args.chunk_size)):
        # slurm-style slicing: drop the first --skip-first ligands and
        # stop after --take (only meaningful for ordered PDBQT globs;
        # SDF chunking + resume makes each task's work idempotent)
        if args.skip_first or args.take is not None:
            slice_seen += len(batch)
            if slice_seen <= args.skip_first:
                continue
            if args.skip_first:
                batch = batch[len(batch) - (slice_seen - args.skip_first):]
            if args.take is not None:
                remaining = args.take - (slice_seen - len(batch) - args.skip_first)
                if remaining <= 0:
                    break
                batch = batch[:remaining]
        pending = [path for path in batch
                   if path.stem.removesuffix("_out") not in recorded]
        total_skipped += len(batch) - len(pending)
        if not pending:
            continue
        jobs = [_job(path) for path in pending]
        chunk_rows: list[dict] = []
        if workers == 1:
            for job in jobs:
                chunk_rows.append(_dock_one(job))
        else:
            ctx = mp.get_context("spawn")
            with ctx.Pool(min(workers, len(jobs))) as pool:
                chunk_rows = list(pool.imap_unordered(_dock_one, jobs))
        append_rows(results_csv, chunk_rows)
        total_docked += len(chunk_rows)
        for result in chunk_rows:
            _report(result)
        print(f"  [chunk {chunk_index + 1}] {total_docked} ligand(s) docked, "
              f"{total_skipped} skipped, "
              f"{time.perf_counter() - started:.0f}s elapsed", flush=True)

    # final ranked summary rebuilt from the incremental results file
    all_rows = load_recorded(results_csv)
    ranked = sorted(
        [r for r in all_rows.values() if r.get("best_affinity")],
        key=lambda r: float(r["best_affinity"]),
    ) + [r for r in all_rows.values() if not r.get("best_affinity")]
    summary = out_dir / "batch_summary.csv"
    _write_csv(ranked, summary)
    print(f"\nfinished in {time.perf_counter() - started:.1f}s")
    print(f"results  : {results_csv} ({len(all_rows)} rows)")
    print(f"summary  : {summary}")
    print("\ntop 10 ligands:")
    for index, result in enumerate(ranked[:10], start=1):
        if result.get("best_affinity"):
            print(f"  {index:3d}. {result['ligand']:<30s} "
                  f"{float(result['best_affinity']):8.2f} kcal/mol")
    return 0 if any(r.get("best_affinity") for r in all_rows.values()) else 1


def _report(result: dict) -> None:
    if result["error"]:
        print(f"  FAILED {result['ligand']}: {result['error']}")
    elif result["best"] is not None:
        print(f"  {result['ligand']}: {result['best']:.2f} kcal/mol "
              f"({result['poses']} poses, {result['runtime']:.0f}s)")


def _write_csv(results: list[dict], path: Path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["rank", "ligand", "best_affinity", "num_poses",
                         "runtime_s", "error", "out_pdbqt"])
        for index, result in enumerate(results, start=1):
            writer.writerow([
                index, result["ligand"],
                f"{float(result['best_affinity']):.3f}"
                if result.get("best_affinity") else "",
                result.get("num_poses", ""),
                f"{float(result['runtime_s']):.1f}"
                if result.get("runtime_s") else "",
                result.get("error") or "", result.get("out_pdbqt") or "",
            ])


if __name__ == "__main__":
    sys.exit(main())
