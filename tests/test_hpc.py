"""HPC template generation tests (work-order item 22 / ADR-0007).

The sbatch template is TEXT: these tests lint it (placeholders resolved,
no stray format braces, correct batch_dock flags) without a scheduler.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

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


def test_template_is_lf_only():
    """The sbatch targets Linux clusters: no CR bytes anywhere (a Windows
    text-mode write would inject CRLF and break bash/slurm - CI found this
    on windows-latest)."""
    text = build_template(_args())
    assert "\r" not in text, "sbatch template must not contain CR bytes"
    assert text.startswith("#!/bin/bash\n"), "missing bash shebang"


def _bash_exe() -> str | None:
    """A real POSIX bash for ``bash -n`` syntax checks, or None.

    On Windows the ``bash`` found on PATH is usually the WSL launcher stub
    (``C:\\Windows\\System32\\bash.exe``) which, with no distribution
    installed, prints a UTF-16 error and exits 1 - that stub made this test
    fail on windows-latest CI.  Git-for-Windows ships a real bash; GH
    runners have it installed.  Anything else -> skip the syntax check.
    """
    import shutil

    if os.name != "nt":
        return shutil.which("bash")
    for candidate in (
        r"C:\Program Files\Git\bin\bash.exe",
        r"C:\Program Files (x86)\Git\bin\bash.exe",
    ):
        if os.path.isfile(candidate):
            return candidate
    return None


def test_generated_file_is_valid_bash_text(tmp_path):
    import subprocess

    out = tmp_path / "job.sbatch"
    args = _args(out=str(out), dockflow_home=str(REPO_ROOT))
    text = build_template(args)
    # LF endings on every platform: sbatch targets Linux clusters.
    out.write_bytes(text.encode("utf-8"))

    bash = _bash_exe()
    if bash is None:
        pytest.skip("no POSIX bash available (Windows without Git-for-Windows)")
    # Git-bash opens mixed-separator paths; POSIX needs no translation.
    bash_path = str(out).replace("\\", "/") if os.name == "nt" else str(out)
    # bash -n: syntax check only (no scheduler needed)
    proc = subprocess.run([bash, "-n", bash_path],
                          capture_output=True, text=True, timeout=30)
    assert proc.returncode == 0, proc.stderr
