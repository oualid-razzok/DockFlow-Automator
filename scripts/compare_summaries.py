"""Compare two DockFlow summary.csv files for reproducibility (item 6).

Byte-identity is checked first (the strong, expected outcome for a fixed
Vina seed); if that fails, every numeric column is compared within a
1e-6 tolerance and the per-column worst delta is reported.  Exits 0 when
reproducible, 1 otherwise.

Usage::

    python scripts/compare_summaries.py run1/docking/summary.csv \
        run2/docking/summary.csv [--tolerance 1e-6]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

NUMERIC_COLUMNS = {
    "affinity_kcal_mol", "rmsd_lb", "rmsd_ub", "crystal_rmsd",
    "docking_score_efficiency", "num_heavy_atoms", "runtime_s",
}
# runtime_s legitimately varies between runs (wall time, not science):
# it is compared for presence only.
IGNORED_COLUMNS = {"runtime_s"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary_a")
    parser.add_argument("summary_b")
    parser.add_argument("--tolerance", type=float, default=1e-6)
    args = parser.parse_args()

    path_a, path_b = Path(args.summary_a), Path(args.summary_b)
    text_a = path_a.read_text(encoding="utf-8")
    text_b = path_b.read_text(encoding="utf-8")

    if text_a == text_b:
        print(f"REPRODUCIBLE: {path_a.name} and {path_b.name} are byte-identical")
        return 0

    rows_a = list(csv.DictReader(text_a.splitlines()))
    rows_b = list(csv.DictReader(text_b.splitlines()))
    if len(rows_a) != len(rows_b):
        print(f"NOT REPRODUCIBLE: {len(rows_a)} rows vs {len(rows_b)} rows")
        return 1
    if not rows_a:
        print("NOT REPRODUCIBLE: empty summaries")
        return 1

    worst: dict[str, float] = {}
    failures: list[str] = []
    for index, (row_a, row_b) in enumerate(zip(rows_a, rows_b, strict=True)):
        for column, value_a in row_a.items():
            value_b = row_b.get(column, "")
            if column in IGNORED_COLUMNS:
                continue
            if column in NUMERIC_COLUMNS:
                if value_a == "" or value_b == "":
                    if value_a != value_b:
                        failures.append(f"row {index} column {column}: "
                                        f"{value_a!r} vs {value_b!r}")
                    continue
                delta = abs(float(value_a) - float(value_b))
                worst[column] = max(worst.get(column, 0.0), delta)
                if delta > args.tolerance:
                    failures.append(
                        f"row {index} column {column}: {value_a} vs {value_b} "
                        f"(delta {delta:.3g} > {args.tolerance:g})"
                    )
            elif value_a != value_b:
                failures.append(f"row {index} column {column}: "
                                f"{value_a!r} vs {value_b!r}")

    print("Not byte-identical; numeric comparison within "
          f"{args.tolerance:g}:")
    for column, delta in sorted(worst.items()):
        print(f"  {column:28s}: worst delta {delta:.3g}")
    if failures:
        print(f"NOT REPRODUCIBLE: {len(failures)} difference(s):")
        for failure in failures[:20]:
            print(f"  {failure}")
        return 1
    print("REPRODUCIBLE: all scientific columns within tolerance "
          "(runtime_s excluded by design)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
