"""Re-run the ANALYSIS stage of every benchmark complex (final pass, item 5).

Why this exists: the spyrmsd cross-check (``rmsd_crosscheck.py``) found
that ``rdkit_symmetry_rmsd`` seeded its atom correspondence geometrically
and could lock in a topologically wrong assignment when the pose and the
crystal reference list atoms in different orders - over-penalising
symmetric ligands by up to ~0.4 A.  The analyzer now builds the pose
molecule with its own graph, making the RMSD atom-order independent.

Docking itself is untouched (bit-identical, seed 2026); this script only
re-runs the analysis + report stages over the committed pose files so
every published crystal RMSD, pose-level artifact and per-complex report
reflects the corrected implementation::

    python benchmarks/redocking/reanalyze.py            # all 24 complexes
    python benchmarks/redocking/reanalyze.py --pdb-ids 1HVR,1STP

After it finishes, regenerate the top-level summary with::

    python benchmarks/redocking/run_benchmark.py --report-only

and re-validate with::

    python benchmarks/redocking/rmsd_crosscheck.py
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dockflow_core.pipeline import DockingPipeline, PipelineConfig  # noqa: E402


def load_complexes(path: Path) -> list[dict[str, str]]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    return [row for row in csv.DictReader(lines) if row.get("pdb_id")]


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--complexes",
                        default=str(Path(__file__).parent / "complexes.csv"))
    parser.add_argument("--out", default=str(Path(__file__).parent / "results"))
    parser.add_argument("--pdb-ids", default=None)
    args = parser.parse_args()

    out_dir = Path(args.out)
    complexes = load_complexes(Path(args.complexes))
    if args.pdb_ids:
        wanted = {pid.strip().upper() for pid in args.pdb_ids.split(",")}
        complexes = [row for row in complexes
                     if row["pdb_id"].upper() in wanted]

    failures: list[str] = []
    for complex_row in complexes:
        pdb_id = complex_row["pdb_id"]
        run_id = f"{pdb_id}_benchmark"
        if not (out_dir / run_id / "docking").is_dir():
            print(f"{pdb_id}: no run directory - skipped", file=sys.stderr)
            continue
        config = PipelineConfig(
            workdir=out_dir,
            run_id=run_id,
            target={"pdb_id": pdb_id},
            ligands=[{"id": f"{complex_row['ligand_resname'].lower()}_redock",
                      "pdb_ligand": complex_row["ligand_resname"]}],
            receptor={"engine": "auto", "charge_model": "gasteiger"},
            gridbox={"source": "ligand",
                     "reference_ligand_resname": complex_row["ligand_resname"],
                     "padding": 4.0},
            docking={"backend": "auto", "scoring": "vina",
                     "exhaustiveness": 8, "num_modes": 9, "seed": 2026,
                     "cpu": 0, "timeout": 3600, "parallel": 1},
            analysis={"top_poses": 3},
            visualization={"enabled": False},
        )
        report = DockingPipeline(config).run(stages=["analyze"])
        poses = report.docking.get("results", [{}])[0].get("poses", []) \
            if report.ok else []
        best = min((p.get("crystal_rmsd") for p in poses
                    if p.get("crystal_rmsd") is not None), default=None)
        status = "ok" if report.ok else f"FAILED: {report.error}"
        print(f"{pdb_id}: {status} | {len(poses)} poses | "
              f"best RMSD {best if best is None else round(best, 3)} A")
        if not report.ok:
            failures.append(pdb_id)

    if failures:
        print(f"\n{len(failures)} complex(es) failed: {', '.join(failures)}",
              file=sys.stderr)
        return 1
    print("\nall complexes re-analysed; next: run_benchmark.py --report-only")
    return 0


if __name__ == "__main__":
    sys.exit(main())
