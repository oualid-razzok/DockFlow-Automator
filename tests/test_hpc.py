"""HPC template generation tests (work-order item 22 / ADR-0007).

The sbatch template is TEXT: these tests lint it (placeholders resolved,
no stray format braces, correct batch_dock flags) without a scheduler.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from generate_sbatch import build_template, chunk_count, parse_args  # noqa: E402


def _args(**overrides) -> object:
    base = [
        "--receptor", "runs/r/prepared/receptor.pdbqt",
        "--ligand-glob", "screening/prepared/*.pdbqt",
        "--center=-8.7,15.5,27.9",  # '=' form: values starting with '-'
    ]
    for key, value in overrides.items():
        flag = f"--{key.replace('_', '-')}"
        if value is True:  # store_true flags take no value
            base.append(flag)
        else:
            base += [flag, str(value)]
    return parse_args(base)


def test_template_resolves_all_placeholders():
    text = build_template(_args())
    assert "{" not in text and "}" not in text  # every placeholder resolved
    assert "#SBATCH --array=0-9" in text
    assert "scripts/batch_dock.py" in text
    assert "--ligands" in text
    assert "--skip-first $((" in text  # array-task slicing wired in
    assert "--take" in text


def test_template_sdf_source_mode():
    text = build_template(_args(sdf_flags=True, out="x"))
    assert '--sdf "$LIGAND_SOURCE"' in text
    assert "--ligands" not in text


def test_template_overrides():
    text = build_template(_args(chunk_size=250, array="0-99%8",
                                exhaustiveness=16, seed=7, cpus=4,
                                mem="8G", partition="gpu-short",
                                out_dir="my_screen"))
    assert "CHUNK=250" in text
    assert "#SBATCH --array=0-99%8" in text
    assert "--exhaustiveness 16" in text
    assert "--seed 7" in text
    assert "#SBATCH --cpus-per-task=4" in text
    assert "#SBATCH --mem=8G" in text
    assert "#SBATCH --partition=gpu-short" in text
    assert "--out-dir my_screen" in text


def test_chunk_count_math():
    assert chunk_count(100_000, 500) == 200
    assert chunk_count(500, 500) == 1
    assert chunk_count(501, 500) == 2
    assert chunk_count(0, 500) == 1  # degenerate but safe


def test_generated_file_is_valid_bash_text(tmp_path):
    import subprocess

    out = tmp_path / "job.sbatch"
    args = _args(out=str(out), dockflow_home=str(REPO_ROOT))
    out.write_text(build_template(args), encoding="utf-8")
    # bash -n: syntax check only (no scheduler needed)
    proc = subprocess.run(["bash", "-n", str(out)],
                          capture_output=True, text=True, timeout=30)
    assert proc.returncode == 0, proc.stderr
