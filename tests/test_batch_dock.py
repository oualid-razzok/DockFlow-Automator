"""Chunked streaming + checkpoint resume for batch_dock (item 23).

The heavy part (docking) is covered by the scientific tests; these
tests exercise the library-streaming and bookkeeping mechanics that
make 100k+ record libraries safe.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(REPO_ROOT))

from batch_dock import (  # noqa: E402
    append_rows,
    iter_ligand_batches,
    load_recorded,
    parse_args,
)

needs_rdkit = pytest.importorskip("rdkit")


def _write_sdf(path: Path, count: int) -> Path:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    with Chem.SDWriter(str(path)) as writer:
        for index in range(count):
            mol = Chem.MolFromSmiles("c1ccccc1")  # benzene, minimal
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=index + 1)
            mol.SetProp("_Name", f"benzene_{index + 1:03d}")
            writer.write(mol)
    return path


def test_sdf_streams_in_chunks(tmp_path):
    """5 records with chunk size 2 -> batches [2, 2, 1] (bounded memory)."""
    sdf = _write_sdf(tmp_path / "library.sdf", 5)
    args = parse_args([
        "--receptor", "r.pdbqt", "--sdf", str(sdf),
        "--center", "0,0,0", "--out-dir", str(tmp_path / "out"),
        "--chunk-size", "2",
    ])
    sizes = [len(batch) for batch in
             iter_ligand_batches(args, tmp_path / "out", 2)]
    assert sizes == [2, 2, 1]
    prepared = list((tmp_path / "out" / "prepared").glob("*.pdbqt"))
    assert len(prepared) == 5  # every record became a PDBQT


def test_max_records_limits_stream(tmp_path):
    sdf = _write_sdf(tmp_path / "library.sdf", 4)
    args = parse_args([
        "--receptor", "r.pdbqt", "--sdf", str(sdf),
        "--center", "0,0,0", "--out-dir", str(tmp_path / "out"),
        "--chunk-size", "10", "--max-records", "3",
    ])
    sizes = [len(batch) for batch in
             iter_ligand_batches(args, tmp_path / "out", 10)]
    assert sizes == [3]


def test_checkpoint_roundtrip(tmp_path):
    """append_rows -> crash -> load_recorded sees exactly those rows."""
    results_csv = tmp_path / "batch_results.csv"
    rows = [
        {"ligand": "lig1", "best_affinity": -9.5, "num_poses": 9,
         "runtime_s": 12.0, "error": None, "out_pdbqt": "out/lig1_out.pdbqt"},
        {"ligand": "lig2", "best_affinity": None, "num_poses": 0,
         "runtime_s": 1.0, "error": "no poses", "out_pdbqt": None},
    ]
    append_rows(results_csv, rows)
    recorded = load_recorded(results_csv)
    assert set(recorded) == {"lig1", "lig2"}
    assert recorded["lig1"]["best_affinity"] == "-9.5"
    # a second append (simulating a resumed run) keeps both generations
    append_rows(results_csv, [
        {"ligand": "lig3", "best_affinity": -8.1, "num_poses": 3,
         "runtime_s": 5.0, "error": None, "out_pdbqt": "out/lig3_out.pdbqt"}])
    assert set(load_recorded(results_csv)) == {"lig1", "lig2", "lig3"}


def test_resume_skips_recorded_ligands(tmp_path):
    """The pending-filter inside main() skips already-recorded stems."""
    from batch_dock import load_recorded

    sdf = _write_sdf(tmp_path / "library.sdf", 3)
    out_dir = tmp_path / "out"
    args = parse_args([
        "--receptor", "r.pdbqt", "--sdf", str(sdf),
        "--center", "0,0,0", "--out-dir", str(out_dir),
        "--chunk-size", "10",
    ])
    batch = next(iter_ligand_batches(args, out_dir, 10))
    assert len(batch) == 3
    # pretend two of the three finished in a previous, interrupted run
    results_csv = out_dir / "batch_results.csv"
    results_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(results_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["ligand", "best_affinity", "num_poses",
                            "runtime_s", "error", "out_pdbqt"])
        writer.writeheader()
        for path in batch[:2]:
            writer.writerow({
                "ligand": path.stem, "best_affinity": "-7.0",
                "num_poses": "1", "runtime_s": "1.0", "error": "",
                "out_pdbqt": f"{path.stem}_out.pdbqt"})
    recorded = load_recorded(results_csv)
    pending = [p for p in batch
               if p.stem.removesuffix("_out") not in recorded]
    assert len(pending) == 1  # only the unfinished ligand re-docks
