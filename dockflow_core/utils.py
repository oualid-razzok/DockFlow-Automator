"""Shared utilities: exceptions, logging, subprocess helpers, file helpers."""

from __future__ import annotations

import hashlib
import logging
import os
import shutil
import subprocess
import sys
import time
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

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


def setup_logging(level: str | int = "INFO", logfile: str | os.PathLike | None = None) -> None:
    """Configure the root ``dockflow`` logger (idempotent)."""
    logger = logging.getLogger("dockflow")
    logger.setLevel(logging.getLevelName(level) if isinstance(level, str) else level)
    if logger.handlers:  # already configured
        return
    formatter = logging.Formatter(_LOG_FORMAT)
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(formatter)
    logger.addHandler(console)
    if logfile:
        Path(logfile).parent.mkdir(parents=True, exist_ok=True)
        filehandler = logging.FileHandler(logfile, encoding="utf-8")
        filehandler.setFormatter(formatter)
        logger.addHandler(filehandler)


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


def is_importable(module_name: str) -> bool:
    """True if ``module_name`` can be imported without importing it now."""
    import importlib.util

    try:
        return importlib.util.find_spec(module_name) is not None
    except (ImportError, ValueError, ModuleNotFoundError):
        return False
