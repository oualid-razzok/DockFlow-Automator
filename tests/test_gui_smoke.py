"""GUI smoke tests (offscreen platform, marked 'gui')."""

from __future__ import annotations

import os

import pytest

pytest.importorskip("PyQt6.QtWidgets",
                    reason="PyQt6 or its system GL libraries unavailable")
pytestmark = pytest.mark.gui


@pytest.fixture(scope="module")
def qapp():
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PyQt6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    yield app


def test_widgets_instantiate(qapp):
    from dockflow_gui.widgets import (
        ContactTableWidget,
        DockParamsWidget,
        GridBoxPreview,
        GridBoxWidget,
        LigandTableWidget,
        LogPanel,
        ResultsTableWidget,
        StepBar,
    )

    assert StepBar().states() == ["pending"] * 6
    params = DockParamsWidget()
    assert params.config_dict()["exhaustiveness"] == 8
    assert params.config_dict()["seed"] is None  # "random"
    box_widget = GridBoxWidget()
    box = box_widget.box()
    assert box.size == (22.0, 22.0, 22.0)
    preview = GridBoxPreview()
    preview.set_data([[0, 0, 0], [1, 1, 1]], [], box)
    LigandTableWidget()
    ResultsTableWidget()
    ContactTableWidget()
    LogPanel().append_log("hello", level="info")


def test_step_bar_states(qapp):
    from dockflow_gui.widgets import StepBar

    bar = StepBar()
    bar.set_state(0, "done")
    bar.set_state(1, "active")
    assert bar.states()[0] == "done"
    assert bar.states()[1] == "active"


def test_ligand_table_content(qapp):

    from dockflow_core.models import LigandRecord
    from dockflow_gui.widgets import LigandTableWidget

    table = LigandTableWidget()
    table.set_ligands([
        LigandRecord(identifier="lig1", source="smiles", value="CCO",
                     status="prepared"),
        LigandRecord(identifier="lig2", source="pubchem", value="aspirin",
                     status="error", error="boom"),
    ])
    assert table.rowCount() == 2
    assert table.item(0, 0).text() == "lig1"
    assert table.item(1, 3).text() == "error"


def test_main_window_smoke(qapp):
    from dockflow_gui.main_window import MainWindow

    window = MainWindow()
    assert window.pages.count() == 6
    assert window.windowTitle() == "DockFlow-Automator"
    # simulated workflow state transitions
    window.stepbar.set_state(0, "done")
    window._log("smoke test", "info")
    assert "smoke test" in window.log_panel.toPlainText()
    window.close()


def test_worker_thread_roundtrip(qapp):
    from dockflow_gui.threads import Worker

    results = []
    worker = Worker(lambda: 42)
    worker.signals.result.connect(results.append)
    worker.start()
    worker.wait(5000)
    assert results == [42]


def test_worker_thread_error(qapp):
    from dockflow_gui.threads import Worker

    errors = []
    worker = Worker(lambda: 1 / 0)
    worker.signals.error.connect(lambda message, tb: errors.append(message))
    worker.start()
    worker.wait(5000)
    assert errors and "ZeroDivisionError" in errors[0]
