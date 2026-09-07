"""Run provenance: stage audit trail and scientific environment fingerprint.

Two concerns live here:

* **Stage audit trail** (audit item 4): every pipeline stage records its
  wall time, ISO-8601 start/stop timestamps, exit status and the resolved
  engine/backend, and the records land in ``manifest["stages"]``.
* **Environment fingerprint** (audit item 5): the full scientific-stack
  version vector (vina, meeko, rdkit, openbabel-wheel, gemmi, numpy,
  python, platform, DockFlow itself and the C++ bindings) recorded under
  ``manifest["environment"]`` and mirrored to a sidecar ``environment.json``
  next to the manifest so two runs can be diffed directly.

  Without this fingerprint "same DockFlow version" does not imply "same
  scientific result": a Vina or Meeko minor bump can shift scores.
"""

from __future__ import annotations

import platform
import subprocess
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .utils import get_logger, module_version, openbabel_version

logger = get_logger("provenance")

__all__ = [
    "StageRecord",
    "StageRecorder",
    "scientific_environment",
    "write_environment_sidecar",
]


def _utc_now() -> str:
    """Current time as an ISO-8601 UTC string (second resolution)."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _bindings_state() -> dict[str, Any]:
    """Importability + version of the optional C++ accelerator."""
    try:
        import dockflow_bindings  # type: ignore

        return {
            "importable": True,
            "version": str(getattr(dockflow_bindings, "__version__", "built")),
        }
    except Exception:  # noqa: BLE001 - optional accelerator
        return {"importable": False, "version": None}


def _git_sha() -> str | None:
    """HEAD commit of the DockFlow checkout when running from a git repo."""
    import shutil

    if shutil.which("git") is None:
        return None
    try:
        result = subprocess.run(  # noqa: S603 - fixed argv, read-only probe
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=10, check=False,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return None


def scientific_environment() -> dict[str, Any]:
    """Full scientific-stack version fingerprint for the current process.

    Deliberately reports ``None`` (never raises) for absent optional
    packages: the fingerprint describes reality, and a missing toolkit is a
    scientific fact worth recording rather than an error.
    """
    return {
        "dockflow": __version__,
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "platform": f"{platform.system()} {platform.release()} ({platform.machine()})",
        "numpy": module_version("numpy"),
        "vina": module_version("vina"),
        "meeko": module_version("meeko"),
        "rdkit": module_version("rdkit"),
        "gemmi": module_version("gemmi"),
        "openbabel_wheel": openbabel_version(),
        "pandas": module_version("pandas"),
        "dockflow_bindings": _bindings_state(),
        "git_sha": _git_sha(),
        "captured_at": _utc_now(),
    }


def write_environment_sidecar(run_dir: str | Path, environment: dict[str, Any]) -> Path:
    """Write ``environment.json`` next to the manifest for easy diffing."""
    import json

    path = Path(run_dir) / "environment.json"
    path.write_text(json.dumps(environment, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Stage recorder
# ---------------------------------------------------------------------------
@dataclass
class StageRecord:
    """Audit record for one pipeline stage."""

    name: str
    status: str = "pending"  # pending | running | done | failed | skipped
    started_at: str | None = None
    ended_at: str | None = None
    duration_s: float | None = None
    engine: str | None = None   # resolved preparation engine (prep stages)
    backend: str | None = None  # resolved docking backend (dock stage)
    detail: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "status": self.status,
            "started_at": self.started_at,
            "ended_at": self.ended_at,
            "duration_s": self.duration_s,
            "engine": self.engine,
            "backend": self.backend,
            "detail": self.detail,
        }


class StageRecorder:
    """Collects per-stage timings and outcomes for the manifest.

    Usage::

        recorder = StageRecorder()
        with recorder.stage("docking") as stage:
            ...
            stage.backend = engine.backend
        manifest["stages"] = recorder.to_list()

    The context manager timestamps entry/exit, computes wall time and flips
    the status to ``failed`` (preserving the exception) or ``done``
    automatically; the body may override the status (e.g. ``skipped``).
    """

    def __init__(self) -> None:
        self._records: dict[str, StageRecord] = {}
        self._t0: dict[str, float] = {}

    @contextmanager
    def stage(self, name: str) -> Iterator[StageRecord]:
        record = StageRecord(name=name, status="running", started_at=_utc_now())
        self._records[name] = record
        self._t0[name] = time.perf_counter()
        try:
            yield record
        except Exception as exc:  # noqa: BLE001 - record, then re-raise
            record.status = "failed"
            record.detail = record.detail or f"{type(exc).__name__}: {exc}"
            self._close(record, name)
            raise
        else:
            if record.status == "running":
                record.status = "done"
            self._close(record, name)

    def _close(self, record: StageRecord, name: str) -> None:
        record.ended_at = _utc_now()
        t0 = self._t0.get(name)
        record.duration_s = round(time.perf_counter() - t0, 3) if t0 is not None else None

    def mark_skipped(self, name: str, detail: str = "") -> None:
        """Flag a recorded stage as skipped (e.g. visualization disabled)."""
        record = self._records.get(name)
        if record is not None:
            record.status = "skipped"
            if detail:
                record.detail = detail

    def to_list(self) -> list[dict[str, Any]]:
        """All stage records, in first-run order, manifest-ready."""
        return [record.to_dict() for record in self._records.values()]
