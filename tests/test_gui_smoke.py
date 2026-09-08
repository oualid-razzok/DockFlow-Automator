"""GUI smoke tests (offscreen platform, marked 'gui')."""

from __future__ import annotations

import os
import time

import pytest

try:
    import PyQt6.QtWidgets  # noqa: F401
except Exception as _exc:  # ImportError: package missing OR system GL libs broken
    pytest.skip(f"PyQt6 or its system GL libraries unavailable ({_exc})",
                allow_module_level=True)
pytestmark = pytest.mark.gui


def _wait_for_worker(qapp, worker, timeout_s: float = 5.0) -> None:
    """Wait for a QThread while keeping the event loop alive.

    Cross-thread Qt signals are queued into the receiving thread's event
    loop; a plain ``worker.wait()`` blocks that loop and the connected
    callbacks would never fire.  ``processEvents()`` between polls drains
    the queue, exactly like ``app.exec()`` does in the running app.
    """
    deadline = time.monotonic() + timeout_s
    while not worker.isFinished() and time.monotonic() < deadline:
        qapp.processEvents()
        time.sleep(0.005)
    qapp.processEvents()  # drain pending queued signal deliveries


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


def test_main_window_smoke(qapp, tmp_path):
    from dockflow_gui.main_window import MainWindow
    from dockflow_gui.tutorial import mark_tutorial_seen

    settings = _fresh_settings(tmp_path)
    mark_tutorial_seen(settings)  # no modal tour inside a smoke test
    window = MainWindow(settings=settings)
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
    _wait_for_worker(qapp, worker)
    assert results == [42]


def test_worker_thread_error(qapp):
    from dockflow_gui.threads import Worker

    errors = []
    worker = Worker(lambda: 1 / 0)
    worker.signals.error.connect(lambda message, tb: errors.append(message))
    worker.start()
    _wait_for_worker(qapp, worker)
    assert errors and "ZeroDivisionError" in errors[0]


# ---------------------------------------------------------------- item 21
def _fresh_settings(tmp_path):
    from PyQt6.QtCore import QSettings

    return QSettings(str(tmp_path / "settings.ini"),
                     QSettings.Format.IniFormat)


def test_tutorial_shown_on_fresh_settings(qapp, tmp_path, monkeypatch):
    """Fresh QSettings -> the tutorial decision fires and marks seen."""
    from dockflow_gui import tutorial

    settings = _fresh_settings(tmp_path)
    assert tutorial.should_show_tutorial(settings) is True

    shown = []

    class _FakeDialog:
        def __init__(self, parent):
            pass

        def exec(self):  # non-blocking stand-in for the modal loop
            shown.append(True)
            return 0

    monkeypatch.setattr(tutorial, "TutorialDialog", _FakeDialog)
    displayed = tutorial.maybe_show_tutorial(None, settings)
    assert displayed is True
    assert shown == [True]
    # dismissed -> not shown again on the SAME settings
    assert tutorial.should_show_tutorial(settings) is False
    assert tutorial.maybe_show_tutorial(None, settings) is False
    assert shown == [True]  # second call did not open a dialog


def test_tutorial_dialog_content(qapp):
    from dockflow_gui.tutorial import TUTORIAL_STEPS, TutorialDialog

    dialog = TutorialDialog()
    assert dialog._stack.count() == 4
    titles = [step["title"] for step in TUTORIAL_STEPS]
    assert titles[0].startswith("Welcome")
    # every page points at a concrete GUI location and mentions guidance
    for step in TUTORIAL_STEPS:
        assert step["where"].strip()
        assert len(step["body"]) > 80
    texts = " ".join(dialog.page_texts())
    assert "manifest" in texts          # scientific decisions are exposed
    assert "USER_GUIDE" in texts        # doc pointers
    assert "heuristic" in texts.lower()  # honest caveats
    # navigation works offscreen
    assert dialog._btn_back.isEnabled() is False
    dialog._go_next()
    assert dialog._stack.currentIndex() == 1
    assert dialog._btn_back.isEnabled() is True


def test_main_window_tutorial_integration(qapp, tmp_path, monkeypatch):
    """MainWindow consults the injected settings for first-run detection."""
    from dockflow_gui import tutorial
    from dockflow_gui.main_window import MainWindow

    settings = _fresh_settings(tmp_path)
    shown = []

    class _FakeDialog:
        def __init__(self, parent):
            pass

        def exec(self):
            shown.append(True)
            return 0

    monkeypatch.setattr(tutorial, "TutorialDialog", _FakeDialog)
    window = MainWindow(settings=settings)
    assert shown == [True]
    assert tutorial.should_show_tutorial(window.settings) is False
    # a second window on the same settings does not re-show the tutorial
    window2 = MainWindow(settings=window.settings)
    assert shown == [True]
    window.close()
    window2.close()


# ---------------------------------------------------------------- item 30
def test_gui_exposes_engine_comparability(qapp, tmp_path):
    """The prepare page shows the manifest's comparability sets (item 30)."""
    from dockflow_gui.main_window import MainWindow
    from dockflow_gui.tutorial import mark_tutorial_seen

    settings = _fresh_settings(tmp_path)
    mark_tutorial_seen(settings)
    window = MainWindow(settings=settings)
    text = window.rec_comparability.text()
    assert "comparab" in text.lower()
    # switching to the 'none' engine warns about incomparability
    window.rec_engine.setCurrentText("none")
    text = window.rec_comparability.text()
    assert "ZERO charges" in text
    # switching to openbabel names its comparability set
    window.rec_engine.setCurrentText("openbabel")
    text = window.rec_comparability.text()
    assert "comparable:" in text
    window.close()


def test_gui_gridbox_assumption_label(qapp, tmp_path):
    """Manual box edits re-annotate the pocket-assumption provenance."""
    from dockflow_gui.main_window import MainWindow
    from dockflow_gui.tutorial import mark_tutorial_seen

    settings = _fresh_settings(tmp_path)
    mark_tutorial_seen(settings)
    window = MainWindow(settings=settings)
    assert window.gridbox_assumption.text()  # non-empty from construction
    # programmatic annotation (what a co-crystal box derivation does)
    window._set_gridbox_assumption("pocket:XK2(A301)")
    assert "strong assumption" in window.gridbox_assumption.text()
    assert "manifest.gridbox.assumption_strength" in window.gridbox_assumption.text()
    # a manual edit hands the assumption back to the user
    from dockflow_core.gridbox import GridBox

    window._on_box_changed(GridBox(center=(0, 0, 0), size=(20, 20, 20)))
    assert "explicit" in window.gridbox_assumption.text()
    window.close()
