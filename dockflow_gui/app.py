"""QApplication bootstrap for the DockFlow desktop application."""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    """Create the application and show the main window."""
    parser = argparse.ArgumentParser(prog="dockflow-gui",
                                     description="DockFlow-Automator desktop app")
    parser.add_argument("--workdir", default=None, help="default working directory")
    parser.add_argument("--debug", action="store_true", help="verbose logging")
    parser.add_argument("--dark", action="store_true", help="dark statusbar styling")
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    from dockflow_core.utils import setup_logging

    setup_logging("DEBUG" if args.debug else "INFO")

    from PyQt6.QtWidgets import QApplication

    from .main_window import MainWindow
    from .resources import APP_STYLESHEET, app_icon

    app = QApplication(argv if argv is not None else sys.argv)
    app.setApplicationName("DockFlow-Automator")
    app.setOrganizationName("DockFlow")
    app.setWindowIcon(app_icon())
    app.setStyleSheet(APP_STYLESHEET)
    if args.workdir:
        import os

        os.environ["DOCKFLOW_WORKDIR"] = args.workdir
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
