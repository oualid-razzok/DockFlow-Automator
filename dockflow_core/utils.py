"""Shared utilities: exceptions, logging, subprocess helpers, file helpers."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

__all__ = [
    "DockFlowError",
    "ExternalToolError",
    "setup_logging",
    "get_logger",
    "run_command",
    "which",
    "ensure_dir",
    "sha256_file",
    "module_version",
    "openbabel_version",
    "is_importable",
    "VersionReport",
    "CommandResult",
]

log = logging.getLogger("dockflow")


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------
class DockFlowError(Exception):
    """Base class for every error raised by DockFlow-Automator."""


class ExternalToolError(DockFlowError):
    """An external tool (vina, obabel, pymol, ...) is missing or failed."""


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
_LOG_FORMAT = "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s"


class JsonLogFormatter(logging.Formatter):
    """One JSON event per line for log aggregation on HPC/batch systems.

    Emits ``timestamp, level, logger, message`` plus optional ``stage``,
    ``run_id`` and ``ligand`` fields (attached via ``extra=``) and the
    exception traceback when present.
    """

    def format(self, record: logging.LogRecord) -> str:
        payload: dict[str, object] = {
            "timestamp": datetime.fromtimestamp(
                record.created, tz=timezone.utc
            ).isoformat(timespec="milliseconds"),
            "level": record.levelname,
            "logger": record.name,
            "message": record.getMessage(),
        }
        for key in ("stage", "run_id", "ligand"):
            value = getattr(record, key, None)
            if value:
                payload[key] = value
        if record.exc_info:
            payload["exception"] = self.formatException(record.exc_info)
        return json.dumps(payload, ensure_ascii=False)


def setup_logging(
    level: str | int = "INFO",
    logfile: str | os.PathLike | None = None,
    log_format: str = "text",
) -> None:
    """Configure the root ``dockflow`` logger.

    Idempotent: re-invocations (e.g. the pipeline adding a log file after
    the CLI configured the console) replace the handler set so exactly one
    console handler and one optional file handler exist, and keep the most
    verbose level requested so far.

    ``log_format="json"`` switches every handler to
    :class:`JsonLogFormatter` (one JSON event per line - audit item 31).
    """
    logger = logging.getLogger("dockflow")
    requested = logging.getLevelName(level) if isinstance(level, str) else level
    current = logger.level
    if not isinstance(current, int) or current == 0 or requested < current:
        logger.setLevel(requested)
    formatter: logging.Formatter
    if log_format == "json":
        formatter = JsonLogFormatter()
    else:
        formatter = logging.Formatter(_LOG_FORMAT)
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(formatter)
    handlers: list[logging.Handler] = [console]
    if logfile:
        Path(logfile).parent.mkdir(parents=True, exist_ok=True)
        filehandler = logging.FileHandler(logfile, encoding="utf-8")
        filehandler.setFormatter(formatter)
        handlers.append(filehandler)
    for existing in list(logger.handlers):
        logger.removeHandler(existing)
    for handler in handlers:
        logger.addHandler(handler)


def get_logger(name: str = "dockflow") -> logging.Logger:
    """Return a namespaced child logger (``dockflow.<name>``)."""
    return logging.getLogger(f"dockflow.{name}").getChild("")


# ---------------------------------------------------------------------------
# Subprocess execution
# ---------------------------------------------------------------------------
@dataclass
class CommandResult:
    """Result of an external command execution."""

    command: list[str]
    returncode: int
    stdout: str = ""
    stderr: str = ""
    runtime: float = 0.0

    @property
    def ok(self) -> bool:
        return self.returncode == 0

    def __str__(self) -> str:  # pragma: no cover - debugging helper
        return f"$ {' '.join(self.command)} -> rc={self.returncode}"


def run_command(
    command: Sequence[str],
    timeout: float | None = None,
    cwd: str | os.PathLike | None = None,
    env: dict[str, str] | None = None,
    on_output: Callable[[str], None] | None = None,
) -> CommandResult:
    """Run an external command, streaming combined stdout/stderr line by line.

    Args:
        command: argv list, e.g. ``["vina", "--version"]``.
        timeout: kill the process after this many seconds.
        cwd: working directory for the child process.
        env: extra environment variables merged into ``os.environ``.
        on_output: callback invoked with each output line (live progress).

    Returns:
        CommandResult with combined output.

    Raises:
        ExternalToolError: if the executable cannot be started or times out.
    """
    merged_env = dict(os.environ)
    if env:
        merged_env.update(env)
    started = time.perf_counter()
    try:
        proc = subprocess.Popen(  # noqa: S603 - caller-controlled argv
            list(command),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=str(cwd) if cwd else None,
            env=merged_env,
        )
    except OSError as exc:  # FileNotFoundError, PermissionError and Windows'
        # "not a valid Win32 application" (WinError 193) are all OSError
        # subclasses - a shebang script or wrong-architecture binary must not
        # crash the caller, it must become a reportable tool error.
        raise ExternalToolError(f"cannot execute {command[0]} ({exc})") from exc

    lines: list[str] = []
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.rstrip("\n")
            lines.append(line)
            if on_output:
                try:
                    on_output(line)
                except Exception:  # pragma: no cover - callback must not kill us
                    log.debug("output callback failed", exc_info=True)
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
        raise ExternalToolError(
            f"command timed out after {timeout}s: {' '.join(map(str, command))}"
        ) from None
    runtime = time.perf_counter() - started
    result = CommandResult(
        command=[str(c) for c in command],
        returncode=proc.returncode,
        stdout="\n".join(lines),
        runtime=runtime,
    )
    if result.stdout:
        log.debug("\n".join(f"    {ln}" for ln in result.stdout.splitlines()[:40]))
    return result


def which(executable: str | os.PathLike) -> str | None:
    """shutil.which with extra explicit-path support."""
    exe = str(executable)
    if os.path.isfile(exe) and (os.access(exe, os.X_OK) or _windows_executable(exe)):
        return exe
    return shutil.which(exe)


_WINDOWS_EXECUTABLE_SUFFIXES = {".exe", ".bat", ".cmd", ".com"}


def _windows_executable(path: str) -> bool:
    """Windows defence-in-depth: PATHEXT-derived suffixes are always runnable.

    ``os.access(path, os.X_OK)`` normally covers these, but a customised
    ``PATHEXT`` on some CI runners or corporate machines can make it lie;
    an explicit suffix check keeps ``vina.exe`` / ``vina.bat`` detection
    reliable.  Always False on POSIX.
    """
    if sys.platform != "win32":
        return False
    return os.path.splitext(path)[1].lower() in _WINDOWS_EXECUTABLE_SUFFIXES


# ---------------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------------
def ensure_dir(path: str | os.PathLike) -> Path:
    """Create a directory (and parents) and return it as Path."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def sha256_file(path: str | os.PathLike, chunk: int = 1 << 16) -> str:
    """Compute the SHA-256 hex digest of a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def timestamped_run_id(prefix: str = "run") -> str:
    """Build a sortable run id such as ``run_20260902-153001``."""
    return time.strftime(f"{prefix}_%Y%m%d-%H%M%S")


def atomic_write_text(path: str | os.PathLike, text: str) -> Path:
    """Write text atomically (temp file + rename)."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(target.suffix + ".part")
    tmp.write_text(text, encoding="utf-8")
    tmp.replace(target)
    return target


def shorten(path: str | os.PathLike, keep: int = 2) -> str:
    """Shorten a long path for display, keeping the last ``keep`` components."""
    parts = Path(path).parts
    return str(Path(*parts[-keep:])) if len(parts) > keep else str(path)


def first_existing(candidates: Iterable[str | os.PathLike]) -> Path | None:
    """Return the first existing path among candidates."""
    for candidate in candidates:
        p = Path(candidate)
        if p.exists():
            return p
    return None


# ---------------------------------------------------------------------------
# Misc
# ---------------------------------------------------------------------------
@dataclass
class VersionReport:
    """Environment / tool version report produced by ``dockflow info``."""

    entries: dict[str, str] = field(default_factory=dict)

    def add(self, name: str, version: str | None, extra: str = "") -> None:
        value = version if version else "not installed"
        if extra:
            value = f"{value} ({extra})"
        self.entries[name] = value

    def as_text(self) -> str:
        width = max(len(name) for name in self.entries) if self.entries else 10
        lines = ["DockFlow-Automator environment report", "-" * 48]
        lines += [f"{name.ljust(width)} : {version}" for name, version in self.entries.items()]
        return "\n".join(lines)


def module_version(module_name: str) -> str | None:
    """Return the installed version of an importable module, if any."""
    try:
        from importlib import metadata

        return metadata.version(module_name)
    except Exception:  # pragma: no cover - not installed / metadata issues
        return None


def openbabel_version() -> str | None:
    """Version of the OpenBabel Python bindings, ``None`` when unusable.

    Detection is based on the actual import (the PyPI distribution is named
    ``openbabel-wheel`` while the import name is ``openbabel``, so a plain
    metadata lookup reports ``not installed`` even when the bindings work).
    A binding that imports but lacks the core ``OBMol`` class is treated as
    unusable rather than reporting a misleading version.
    """
    try:
        from openbabel import openbabel as ob  # noqa: F401
    except Exception:  # noqa: BLE001 - ImportError or wheel-specific breakage
        return None
    if not hasattr(ob, "OBMol"):
        return None
    for distribution in ("openbabel-wheel", "openbabel"):
        version = module_version(distribution)
        if version:
            return version
    return "importable (version unknown)"


def is_importable(module_name: str) -> bool:
    """True if ``module_name`` can be imported without importing it now."""
    import importlib.util

    try:
        return importlib.util.find_spec(module_name) is not None
    except (ImportError, ValueError, ModuleNotFoundError):
        return False
