"""Regenerate the published redocking benchmark and ASSERT it matches (item 14).

Every published number must be regenerable from the committed scripts +
datasets.  This script makes that claim checkable:

1. **dataset pinning**: verifies the SHA-256 of every committed raw PDB
   against ``complexes.csv`` (column ``pdb_sha256``), so silent RCSB
   drift is detected;
2. **docking regeneration**: re-runs the full 24-complex benchmark into
   a scratch directory (``results-regen/``) with the published settings
   (exhaustiveness 8, seed 2026) and the thread count PINNED to the
   value the published results were produced with (``--cpu 2``): Vina's
   search is thread-count dependent (the v0.3.1 CI lesson - same seed,
   different thread count, different trajectory), so pinning the thread
   count makes the regeneration machine-independent;
3. **assertion**: compares the regenerated CSV against the committed one
   on every numeric column (best affinity, best/Rank-1/top-3 RMSD) with
   tolerance 1e-6 by default, and on the status/success flags exactly;
4. **cross-check** (optional, ``--with-crosscheck``): re-validates the
   symmetry-aware RMSD implementation against spyrmsd on the
   regenerated poses.

Usage::

    python benchmarks/redocking/regenerate.py
    python benchmarks/redocking/regenerate.py --tolerance 1e-3
    python benchmarks/redocking/regenerate.py --skip-docking
        # only verify dataset pins + report-only regeneration

Exit code 0 = the published numbers are regenerable; non-zero = they are
not (the failure is printed, never hidden).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
RESULTS = HERE / "results"
PUBLISHED_CPU = 2  # the published results were produced with 2 threads

NUMERIC_COLUMNS = [
    "best_affinity_kcal_mol", "best_crystal_rmsd_a", "rank1_rmsd_a",
    "top3_rmsd_a",
]
EXACT_COLUMNS = ["status", "success", "success_rank1", "success_top3",
                 "engine", "num_poses"]


def load_complexes() -> list[dict]:
    lines = [line for line in (HERE / "complexes.csv").read_text(
        encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")]
    return [row for row in csv.DictReader(lines) if row.get("pdb_id")]


def verify_dataset_pins(complexes: list[dict]) -> list[str]:
    """SHA-256 of every committed raw PDB must match complexes.csv."""
    problems: list[str] = []
    for row in complexes:
        pdb_id = row["pdb_id"]
        raw = RESULTS / f"{pdb_id}_benchmark" / "raw" / f"{pdb_id.lower()}.pdb"
        if not raw.is_file():
            problems.append(f"{pdb_id}: raw PDB not committed")
            continue
        digest = hashlib.sha256(raw.read_bytes()).hexdigest()
        if digest != row.get("pdb_sha256"):
            problems.append(
                f"{pdb_id}: raw PDB SHA-256 drift (committed file hashes to "
                f"{digest[:12]}, complexes.csv pins "
                f"{str(row.get('pdb_sha256'))[:12]})")
    return problems


def load_rows(path: Path) -> dict[str, dict]:
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {row["pdb_id"]: row for row in rows if row.get("pdb_id")}


def compare_csvs(regenerated: Path, published: Path,
                 tolerance: float) -> list[str]:
    regen = load_rows(regenerated)
    committed = load_rows(published)
    problems: list[str] = []
    if set(regen) != set(committed):
        problems.append(
            f"complex sets differ: regen-only {sorted(set(regen) - set(committed))}"
            f", committed-only {sorted(set(committed) - set(regen))}")
    for pdb_id in sorted(set(regen) & set(committed)):
        fresh, want = regen[pdb_id], committed[pdb_id]
        for column in EXACT_COLUMNS:
            if str(fresh.get(column, "")) != str(want.get(column, "")):
                problems.append(
                    f"{pdb_id}.{column}: regen={fresh.get(column)!r} vs "
                    f"published={want.get(column)!r}")
        for column in NUMERIC_COLUMNS:
            fresh_value = fresh.get(column)
            want_value = want.get(column)
            if (fresh_value or "") == "" or (want_value or "") == "":
                if (fresh_value or "") != (want_value or ""):
                    problems.append(
                        f"{pdb_id}.{column}: missing on one side "
                        f"(regen={fresh_value!r}, published={want_value!r})")
                continue
            delta = abs(float(fresh_value) - float(want_value))
            if delta > tolerance:
                problems.append(
                    f"{pdb_id}.{column}: |delta| {delta:.3e} > {tolerance:g} "
                    f"(regen={fresh_value}, published={want_value})")
    return problems


def run_benchmark(scratch: Path, cpu: int) -> None:
    cmd = [sys.executable, str(HERE / "run_benchmark.py"),
           "--out", str(scratch), "--cpu", str(cpu),
           "--exhaustiveness", "8", "--seed", "2026"]
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tolerance", type=float, default=1e-6,
                        help="max |delta| for numeric columns (default 1e-6)")
    parser.add_argument("--cpu", type=int, default=PUBLISHED_CPU,
                        help=f"thread count to pin (default {PUBLISHED_CPU}: "
                             "the published results' thread count)")
    parser.add_argument("--skip-docking", action="store_true",
                        help="verify dataset pins only (no re-docking)")
    parser.add_argument("--with-crosscheck", action="store_true",
                        help="also run the spyrmsd cross-check on the "
                             "regenerated poses")
    args = parser.parse_args()

    complexes = load_complexes()
    print(f"1) verifying dataset pins for {len(complexes)} complexes ...")
    problems = verify_dataset_pins(complexes)
    for problem in problems:
        print(f"   PIN FAILURE: {problem}", file=sys.stderr)
    if problems:
        return 1
    print("   all raw PDB SHA-256 digests match complexes.csv")

    if args.skip_docking:
        print("2) --skip-docking: done")
        return 0

    scratch = HERE / "results-regen"
    print(f"2) regenerating the benchmark into {scratch} "
          f"(cpu pinned to {args.cpu}) ...")
    run_benchmark(scratch, args.cpu)

    print("3) comparing regenerated vs published CSV ...")
    problems = compare_csvs(scratch / "redocking_results.csv",
                            RESULTS / "redocking_results.csv",
                            args.tolerance)
    for problem in problems:
        print(f"   MISMATCH: {problem}", file=sys.stderr)
    if problems:
        print(f"REGENERATION FAILED: {len(problems)} mismatch(es) beyond "
              f"tolerance {args.tolerance:g}", file=sys.stderr)
        return 1

    stats = json.loads((scratch / "summary_stats.json").read_text())
    published_stats = json.loads(
        (RESULTS / "summary_stats.json").read_text())
    for key in ("completed", "rank1_successful_within_2a",
                "successful_within_2a"):
        if stats.get(key) != published_stats.get(key):
            print(f"   STATS MISMATCH: {key} regen={stats.get(key)} vs "
                  f"published={published_stats.get(key)}", file=sys.stderr)
            return 1

    print("REGENERATION PASSED: every published number reproduced "
          f"within {args.tolerance:g}")
    if args.with_crosscheck:
        cmd = [sys.executable, str(HERE / "rmsd_crosscheck.py"),
               "--results", str(scratch), "--out", str(scratch)]
        print(" ".join(cmd), flush=True)
        return subprocess.run(cmd, check=False).returncode
    return 0


if __name__ == "__main__":
    sys.exit(main())
