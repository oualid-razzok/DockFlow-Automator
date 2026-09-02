"""Global application configuration (work directories, tool paths, CPU count).

Configuration resolution order (highest priority first):

1. Explicit values passed to :class:`AppConfig`.
2. Environment variables: ``DOCKFLOW_HOME``, ``DOCKFLOW_WORKDIR``,
   ``DOCKFLOW_VINA``, ``DOCKFLOW_SMINA``, ``DOCKFLOW_PYMOL``,
   ``DOCKFLOW_OBABEL``, ``DOCKFLOW_CPU``.
3. The YAML file ``$DOCKFLOW_HOME/config.yaml`` (created on first save).
4. Built-in defaults.
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .utils import get_logger, which

logger = get_logger("config")

CONFIG_FILENAME = "config.yaml"


def default_home() -> Path:
    """Directory holding config, caches and default runs."""
    env = os.environ.get("DOCKFLOW_HOME")
    if env:
        return Path(env).expanduser()
    return Path.home() / ".dockflow"


def default_workdir() -> Path:
    """Default directory for docking runs."""
    env = os.environ.get("DOCKFLOW_WORKDIR")
    if env:
        return Path(env).expanduser()
    return default_home() / "runs"


@dataclass
class AppConfig:
    """User-level application configuration."""

    home: Path = field(default_factory=default_home)
    workdir: Path = field(default_factory=default_workdir)
    cache_dir: Path = field(default=None)  # type: ignore[assignment]
    vina_exec: str | None = None
    smina_exec: str | None = None
    pymol_exec: str | None = None
    obabel_exec: str | None = None
    cpu: int = 0  # 0 -> use all cores
    parallel: int = 1  # ligands docked concurrently (CLI backend)
    log_level: str = "INFO"

    def __post_init__(self) -> None:
        self.home = Path(self.home)
        self.workdir = Path(self.workdir)
        if self.cache_dir is None:
            self.cache_dir = self.home / "cache"
        self.cache_dir = Path(self.cache_dir)
        # Environment overrides for tool executables.
        self.vina_exec = os.environ.get("DOCKFLOW_VINA", self.vina_exec)
        self.smina_exec = os.environ.get("DOCKFLOW_SMINA", self.smina_exec)
        self.pymol_exec = os.environ.get("DOCKFLOW_PYMOL", self.pymol_exec)
        self.obabel_exec = os.environ.get("DOCKFLOW_OBABEL", self.obabel_exec)
        env_cpu = os.environ.get("DOCKFLOW_CPU")
        if env_cpu:
            try:
                self.cpu = int(env_cpu)
            except ValueError:
                logger.warning("ignoring non-integer DOCKFLOW_CPU=%r", env_cpu)

    # -- persistence --------------------------------------------------------
    @property
    def config_path(self) -> Path:
        return self.home / CONFIG_FILENAME

    def save(self) -> Path:
        """Persist the configuration to ``$DOCKFLOW_HOME/config.yaml``."""
        self.home.mkdir(parents=True, exist_ok=True)
        data = {k: str(v) if isinstance(v, Path) else v for k, v in asdict(self).items()}
        with open(self.config_path, "w", encoding="utf-8") as fh:
            yaml.safe_dump(data, fh, sort_keys=False)
        return self.config_path

    @classmethod
    def load(cls, path: str | Path | None = None) -> AppConfig:
        """Load configuration from a YAML file (missing file -> defaults)."""
        target = Path(path) if path else default_home() / CONFIG_FILENAME
        if not target.is_file():
            return cls()
        try:
            with open(target, encoding="utf-8") as fh:
                data: dict[str, Any] = yaml.safe_load(fh) or {}
        except yaml.YAMLError as exc:
            logger.warning("could not parse %s (%s); using defaults", target, exc)
            return cls()
        valid = {f for f in cls.__dataclass_fields__}  # noqa: C416 - readability
        data = {k: v for k, v in data.items() if k in valid and v is not None}
        cfg = cls(**data)
        logger.debug("loaded config from %s", target)
        return cfg

    # -- tool resolution ----------------------------------------------------
    def resolve_tool(self, name: str, configured: str | None) -> str | None:
        """Resolve a tool executable: explicit setting -> PATH lookup."""
        if configured:
            found = which(configured)
            if found:
                return found
            logger.warning("configured %s executable not found: %s", name, configured)
        return which(name)


_cached_config: AppConfig | None = None


def get_config(reload: bool = False) -> AppConfig:
    """Return the process-wide application configuration (singleton)."""
    global _cached_config
    if _cached_config is None or reload:
        _cached_config = AppConfig.load()
    return _cached_config
