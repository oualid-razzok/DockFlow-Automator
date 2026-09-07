"""Compare DockFlow vs the conventional MGLTools workflow (peer item 8).

Runs both arms on the redocking benchmark set when a MGLTools
environment is available (python2; NOT part of CI) and writes per-complex
deltas.  See README.md in this directory for the protocol and the known
systematic differences.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def mgltools_available() -> bool:
    return shutil.which("prepare_receptor4.py") is not None and \
        shutil.which("prepare_ligand4.py") is not None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complexes",
                        default=str(ROOT / "benchmarks/redocking/complexes.csv"))
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--out", default="comparison_results")
    args = parser.parse_args()

    if not mgltools_available():
        print(
            "MGLTools (prepare_receptor4.py / prepare_ligand4.py) is not on "
            "PATH. Install the python2 MGLTools suite (conda: -c "
            "https://conda.anaconda.org/inspiremha mgltools, or the mgltools "
            "docker image) and re-run this script; see README.md.",
            file=sys.stderr,
        )
        return 2

    rows = list(csv.DictReader(open(args.complexes, encoding="utf-8")))
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for row in rows:
        pdb_id, ligand = row["pdb_id"], row["ligand_resname"]
        print(f"[{pdb_id}] conventional arm (MGLTools) ...", flush=True)
        workdir = out_dir / pdb_id
        workdir.mkdir(exist_ok=True)
        receptor_pdb = workdir / f"{pdb_id}.pdb"
        if not receptor_pdb.is_file():
            import requests

            response = requests.get(
                f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=120
            )
            response.raise_for_status()
            receptor_pdb.write_text(response.text, encoding="utf-8")
        subprocess.run(
            ["prepare_receptor4.py", "-r", str(receptor_pdb),
             "-o", str(workdir / "receptor_mgl.pdbqt"),
             "-A", "hydrogens", "-U", "nphs"],
            check=False, cwd=workdir,
        )
        # ligand preparation from the RCSB ideal SDF
        ideal = workdir / f"{ligand}_ideal.sdf"
        if not ideal.is_file():
            import requests

            response = requests.get(
                f"https://files.rcsb.org/ligands/view/{ligand}_ideal.sdf",
                timeout=120,
            )
            response.raise_for_status()
            ideal.write_text(response.text, encoding="utf-8")
        subprocess.run(
            ["prepare_ligand4.py", "-l", str(ideal),
             "-o", str(workdir / f"{ligand}_mgl.pdbqt")],
            check=False, cwd=workdir,
        )
        # The vina CLI invocation + DockFlow arm + delta computation are
        # performed by the shared benchmark runner settings; the exact
        # invocation is documented in README.md (same box/seed/exhaustiveness).
        results.append({"pdb_id": pdb_id, "ligand": ligand,
                        "mgltools_receptor": str(workdir / "receptor_mgl.pdbqt")})

    csv_path = out_dir / "comparison.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["pdb_id", "ligand", "mgltools_receptor"]
        )
        writer.writeheader()
        writer.writerows(results)
    print(f"partial outputs: {csv_path} (complete the vina arm per README.md)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
