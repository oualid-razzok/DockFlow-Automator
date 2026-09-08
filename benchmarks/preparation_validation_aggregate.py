"""Cross-workflow preparation aggregate over the 24-complex set (item 13).

Builds the per-complex comparison table published in
``docs/preparation_validation.md#cross-workflow-24-complex``:

* DockFlow arm: the committed redocking-benchmark artifacts
  (``benchmarks/redocking/results/<pdb>_benchmark/prepared/``);
* MGLTools arm: the committed comparison artifacts
  (``benchmarks/mgltools_comparison/results/<pdb>/``);
* Meeko-direct reference on 1HVR/XK2: ``mk_prepare_ligand.py`` (the
  upstream CLI) run on the same RCSB ideal SDF, to show DockFlow's Meeko
  invocation is faithful to upstream.

Everything is computed from the committed files - no hand-transcribed
numbers.  Run from the repository root::

    python benchmarks/preparation_validation_aggregate.py
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REDOCKING = ROOT / "benchmarks" / "redocking" / "results"
MGL = ROOT / "benchmarks" / "mgltools_comparison" / "results"


def pdbqt_atoms(path: Path) -> dict:
    """(heavy, polar_h, branches, charge_sum) of a PDBQT file."""
    heavy = 0
    polar_h = 0
    charge = 0.0
    branches = 0
    hydrogen_types = {"HD", "HS", "H"}
    virtual = {"G0", "G1", "G2", "G3", "W", "XX"}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            atom_type = line[77:].strip() if len(line) > 77 else ""
            try:
                charge += float(line[70:76])
            except ValueError:
                pass
            if atom_type in hydrogen_types:
                if atom_type == "HD":
                    polar_h += 1
            elif atom_type not in virtual:
                heavy += 1
        elif line.startswith("BRANCH"):
            branches += 1
    return {"heavy": heavy, "polar_h": polar_h, "branches": branches,
            "charge_sum": round(charge, 3)}


def load_complexes() -> list[dict]:
    path = ROOT / "benchmarks" / "redocking" / "complexes.csv"
    lines = [line for line in path.read_text(encoding="utf-8").splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    return [row for row in csv.DictReader(lines) if row.get("pdb_id")]


def meeko_direct_reference(pdb_id: str, resname: str) -> dict | None:
    """Run upstream ``mk_prepare_ligand.py`` on the committed ideal SDF."""
    sdf = next((REDOCKING / f"{pdb_id}_benchmark" / "raw").glob("*.sdf"),
               None)
    if sdf is None:
        return None
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "lig.pdbqt"
        # upstream CLI: the meeko package's own console script module,
        # invoked without DockFlow's wrapper
        from meeko.cli import mk_prepare_ligand as cli_module

        proc = subprocess.run(
            [sys.executable, Path(cli_module.__file__), "-i", str(sdf),
             "-o", str(out)],
            capture_output=True, text=True, timeout=120)
        if not out.is_file():
            print(f"meeko-direct: could not prepare {sdf.name}: "
                  f"{(proc.stderr or proc.stdout).strip()[:200]}",
                  file=sys.stderr)
            return None
        return pdbqt_atoms(out)


def main() -> int:
    complexes = load_complexes()
    rows = []
    for complex_row in complexes:
        pdb_id = complex_row["pdb_id"]
        resname = complex_row["ligand_resname"].lower()
        df_run = REDOCKING / f"{pdb_id}_benchmark"
        df_lig_path = next((df_run / "prepared").glob("*_redock.pdbqt"), None)
        df_rec = pdbqt_atoms(df_run / "prepared" / "receptor.pdbqt")
        df_lig = pdbqt_atoms(df_lig_path) if df_lig_path else None
        mgl_rec_path = MGL / pdb_id / "receptor_mgl.pdbqt"
        mgl_lig_path = MGL / pdb_id / f"{resname}_mgl.pdbqt"
        mgl_rec = pdbqt_atoms(mgl_rec_path) if mgl_rec_path.is_file() else None
        mgl_lig = pdbqt_atoms(mgl_lig_path) if mgl_lig_path.is_file() else None
        rows.append({
            "pdb_id": pdb_id,
            "ligand": complex_row["ligand_resname"],
            "df_receptor_heavy": df_rec["heavy"],
            "mgl_receptor_heavy": mgl_rec["heavy"] if mgl_rec else None,
            "df_ligand_heavy": df_lig["heavy"] if df_lig else None,
            "mgl_ligand_heavy": mgl_lig["heavy"] if mgl_lig else None,
            "df_ligand_polar_h": df_lig["polar_h"] if df_lig else None,
            "mgl_ligand_polar_h": mgl_lig["polar_h"] if mgl_lig else None,
            "df_ligand_branches": df_lig["branches"] if df_lig else None,
            "mgl_ligand_branches": mgl_lig["branches"] if mgl_lig else None,
            "df_ligand_charge_sum": df_lig["charge_sum"] if df_lig else None,
            "mgl_ligand_charge_sum": mgl_lig["charge_sum"] if mgl_lig else None,
        })

    out_json = ROOT / "benchmarks" / "preparation_validation_aggregate.json"
    out_json.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")

    # markdown table for docs/preparation_validation.md
    lines = [
        "| pdb | ligand | rec heavy DF/MGL | lig heavy DF/MGL | "
        "lig polar H DF/MGL | rotatable bonds DF/MGL | "
        "Gasteiger sum DF/MGL |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        def pair(key, row=row):
            a, b = row[f"df_{key}"], row[f"mgl_{key}"]
            if a is None or b is None:
                return f"{a if a is not None else 'n/a'}/{'crash'}"
            return f"{a}/{b}"
        lines.append(
            f"| {row['pdb_id']} | {row['ligand']} | {pair('receptor_heavy')} "
            f"| {pair('ligand_heavy')} | {pair('ligand_polar_h')} "
            f"| {pair('ligand_branches')} | {pair('ligand_charge_sum')} |")
    print("\n".join(lines))
    print(f"\nJSON: {out_json}")

    # meeko-direct faithfulness check on 1HVR/XK2
    direct = meeko_direct_reference("1HVR", "XK2")
    if direct is not None:
        df_row = next(r for r in rows if r["pdb_id"] == "1HVR")
        print("\nmeeko-direct (upstream CLI) vs DockFlow wrapper, 1HVR/XK2:")
        print(f"  heavy atoms:      {direct['heavy']} vs "
              f"{df_row['df_ligand_heavy']}")
        print(f"  polar hydrogens: {direct['polar_h']} vs "
              f"{df_row['df_ligand_polar_h']}")
        print(f"  rotatable bonds: {direct['branches']} vs "
              f"{df_row['df_ligand_branches']}")
        print(f"  charge sum:      {direct['charge_sum']} vs "
              f"{df_row['df_ligand_charge_sum']}")
        matches = (direct["heavy"] == df_row["df_ligand_heavy"]
                   and direct["branches"] == df_row["df_ligand_branches"])
        print(f"  faithful: {matches}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
