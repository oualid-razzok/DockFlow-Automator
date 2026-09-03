"""DockFlow-Automator core package.

Unified, automated molecular docking:

    target/ligand download -> preparation -> grid box -> docking -> 3D view

The package is intentionally dependency-light at import time.  Heavy
scientific dependencies (meeko, rdkit, openbabel, vina, pymol) are imported
lazily inside the functions that need them, and every module degrades
gracefully with informative errors when an optional tool is missing.

Quick start::

    from dockflow_core import get_config
    from dockflow_core.pipeline import DockingPipeline, PipelineConfig

    pipeline = DockingPipeline(PipelineConfig.from_yaml("run.yaml"))
    report = pipeline.run()
"""

from .config import AppConfig, get_config
from .utils import DockFlowError


def _resolve_version() -> str:
    """Prefer the installed distribution metadata; fall back to a literal.

    This keeps ``__version__`` in sync with ``pyproject.toml`` when the
    package is installed (regular or editable) while still working when the
    repository is used straight from a source checkout without installation.
    """
    try:
        from importlib import metadata

        return metadata.version("dockflow-automator")
    except Exception:  # pragma: no cover - source checkout without install
        return "0.1.1"


__version__ = _resolve_version()

__all__ = [
    "__version__",
    "AppConfig",
    "get_config",
    "DockFlowError",
]
