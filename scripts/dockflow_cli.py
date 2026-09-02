#!/usr/bin/env python3
"""Standalone DockFlow CLI launcher.

Works with or without installation: when the package is not installed,
the repository root is prepended to ``sys.path`` so ``dockflow_core`` and
``dockflow_gui`` resolve from source.

Usage:
    python scripts/dockflow_cli.py <command> [options]
    python scripts/dockflow_cli.py --help
"""

from __future__ import annotations

import sys
from pathlib import Path


def _ensure_importable() -> None:
    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))


def main() -> int:
    _ensure_importable()
    from dockflow_core.cli import main as cli_main

    return cli_main(sys.argv[1:])


if __name__ == "__main__":
    sys.exit(main())
