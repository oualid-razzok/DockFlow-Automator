"""Provenance tests: stage audit trail and environment fingerprint (items 4-5)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dockflow_core.provenance import (
    StageRecorder,
    scientific_environment,
    write_environment_sidecar,
)


def test_stage_recorder_records_done_stage():
    recorder = StageRecorder()
    with recorder.stage("download") as record:
        record.detail = "1HVR"
    stages = recorder.to_list()
    assert len(stages) == 1
    entry = stages[0]
    assert entry["name"] == "download"
    assert entry["status"] == "done"
    assert entry["started_at"] and entry["ended_at"]
    assert entry["duration_s"] is not None and entry["duration_s"] >= 0.0
    assert entry["detail"] == "1HVR"


def test_stage_recorder_records_failure_and_reraises():
    recorder = StageRecorder()
    with pytest.raises(ValueError, match="boom"), recorder.stage("docking"):
        raise ValueError("boom")
    entry = recorder.to_list()[0]
    assert entry["status"] == "failed"
    assert "ValueError: boom" in entry["detail"]
    assert entry["duration_s"] is not None


def test_stage_recorder_allows_skipped_status():
    recorder = StageRecorder()
    with recorder.stage("visualization") as record:
        record.status = "skipped"
        record.detail = "disabled in config"
    entry = recorder.to_list()[0]
    assert entry["status"] == "skipped"


def test_stage_recorder_engine_backend_fields():
    recorder = StageRecorder()
    with recorder.stage("prepare_receptor") as record:
        record.engine = "openbabel"
    with recorder.stage("docking") as record:
        record.backend = "python"
    stages = {entry["name"]: entry for entry in recorder.to_list()}
    assert stages["prepare_receptor"]["engine"] == "openbabel"
    assert stages["docking"]["backend"] == "python"


def test_stage_recorder_iso_timestamps():
    recorder = StageRecorder()
    with recorder.stage("gridbox"):
        pass
    entry = recorder.to_list()[0]
    # ISO-8601 UTC with Z suffix
    assert entry["started_at"].endswith("Z") and "T" in entry["started_at"]
    assert entry["ended_at"] >= entry["started_at"]


def test_scientific_environment_fingerprint_shape():
    environment = scientific_environment()
    for key in ("dockflow", "python", "platform", "numpy", "vina", "meeko",
                "rdkit", "gemmi", "openbabel_wheel", "pandas",
                "dockflow_bindings", "git_sha", "captured_at"):
        assert key in environment, f"missing fingerprint key: {key}"
    # absent optional packages report None, never raise
    assert environment["dockflow"]
    assert environment["python"]
    assert environment["captured_at"].endswith("Z")


def test_scientific_environment_never_raises():
    """The fingerprint must be collectable in any environment."""
    import subprocess

    code = (
        "from dockflow_core.provenance import scientific_environment;"
        "print(scientific_environment()['dockflow'])"
    )
    result = subprocess.run(
        ["python", "-c", code], capture_output=True, text=True, timeout=60, check=False
    )
    assert result.returncode == 0, result.stderr


def test_environment_sidecar_written(tmp_path: Path):
    environment = scientific_environment()
    path = write_environment_sidecar(tmp_path, environment)
    assert path.name == "environment.json"
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert loaded["dockflow"] == environment["dockflow"]
    assert list(loaded) == sorted(loaded)  # stable, diff-friendly ordering
