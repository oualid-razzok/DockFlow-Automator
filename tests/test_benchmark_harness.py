"""Benchmark-harness tests: worker-failure forensics (never swallow a crash).

The v1.0.0 nightly regression this guards against: every redocking worker
died hard on CI runners and the driver reported only the opaque
``error: worker crashed or produced no result`` because it *discarded* the
worker's captured output.  These tests lock in:

* signal-aware exit-code naming (a SIGSEGV reads as ``killed by SIGSEGV``);
* output-tail extraction (the crash output sits at the END of the stream);
* the ``worker_failures.log`` artifact written next to the results;
* the worker-mode guarantee that a worker ALWAYS prints exactly one JSON
  line, even when its result cannot be serialized or its input table is
  unreadable (verified end-to-end through a real subprocess).
"""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_BENCHMARK = REPO_ROOT / "benchmarks" / "redocking" / "run_benchmark.py"


def _load_run_benchmark():
    spec = importlib.util.spec_from_file_location("run_benchmark", RUN_BENCHMARK)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def harness():
    return _load_run_benchmark()


def test_exit_code_text_names_the_signal(harness):
    assert harness._exit_code_text(1) == "exit code 1"
    assert harness._exit_code_text(0) == "exit code 0"
    assert harness._exit_code_text(-11) == "killed by SIGSEGV (-11)"
    assert harness._exit_code_text(-9) == "killed by SIGKILL (-9)"


def test_tail_keeps_the_last_non_empty_lines(harness):
    text = "\n".join(f"line {i}" for i in range(10)) + "\n\n\n"
    tail = harness._tail(text, 3)
    assert tail.splitlines() == ["line 7", "line 8", "line 9"]
    assert harness._tail("", 5) == ""


def test_report_worker_failure_writes_log_and_prints_detail(
        harness, tmp_path, capsys):
    proc = SimpleNamespace(
        returncode=-11,
        stderr="INFO start\nINFO maps\nSegmentation fault",
        stdout="Computing Vina grid ... done.",
    )
    harness._report_worker_failure(tmp_path, "1HVR", proc, note=None)

    printed = capsys.readouterr().out
    assert "1HVR" in printed
    assert "killed by SIGSEGV (-11)" in printed
    assert "Segmentation fault" in printed
    assert "Computing Vina grid" in printed

    log_file = tmp_path / "worker_failures.log"
    assert log_file.is_file()
    content = log_file.read_text(encoding="utf-8")
    assert "worker failure: 1HVR" in content
    assert "killed by SIGSEGV (-11)" in content


def test_worker_mode_always_prints_one_json_line(tmp_path):
    """End-to-end: a worker with an unreadable complexes table must still
    answer with a JSON status line (and exit 0), never a bare traceback."""
    proc = subprocess.run(
        [sys.executable, str(RUN_BENCHMARK),
         "--complexes", str(tmp_path / "no-such-table.csv"),
         "--out", str(tmp_path / "out"),
         "--exhaustiveness", "8", "--seed", "2026", "--cpu", "2",
         "--worker", "1HVR"],
        capture_output=True, text=True, timeout=120,
        cwd=REPO_ROOT,
    )
    assert proc.returncode == 0, proc.stderr
    lines = [line for line in proc.stdout.strip().splitlines() if line.strip()]
    assert lines, "worker printed nothing"
    payload = json.loads(lines[-1])
    assert payload["pdb_id"] == "1HVR"
    assert payload["status"].startswith("error:")
