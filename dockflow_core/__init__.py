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

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "AppConfig",
    "get_config",
    "DockFlowError",
]
