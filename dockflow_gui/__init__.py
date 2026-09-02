"""DockFlow-Automator desktop GUI (PyQt6).

This package intentionally does not import PyQt6 at import time so that
``import dockflow_gui`` never fails on headless machines.  Use
:func:`dockflow_gui.app.main` to launch the application.
"""

__version__ = "0.1.0"

__all__ = ["main"]


def main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin wrapper
    """Launch the DockFlow desktop application."""
    from .app import main as _main

    return _main(argv)
