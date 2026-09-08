"""First-run tutorial for the DockFlow GUI (work-order item 21).

Presentation ONLY: this module contains the welcome dialog's text and
layout and the QSettings bookkeeping that decides whether to show it.
It deliberately holds no docking/pipeline logic - the guiding principle
of the project is "automate mechanical decisions, expose scientific
decisions", and the tutorial's job is to point at where those decisions
live (which is documentation, not logic).

First-run detection is version-keyed: bumping ``TUTORIAL_VERSION`` (new
tutorial content in a release) re-shows the dialog once for existing
users, while ``Help > First-run tutorial`` always shows it.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QLabel,
    QPushButton,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

TUTORIAL_VERSION = 1
SEEN_KEY = "tutorial/seen_version"

#: The four tutorial steps, in order.  ``where`` names the GUI step-bar
#: entry (the left-hand navigation) so the text can reference concrete
#: screen coordinates instead of vague prose.  All content is text-only.
TUTORIAL_STEPS: list[dict[str, str]] = [
    {
        "title": "Welcome to DockFlow",
        "where": "the left step bar",
        "body": (
            "DockFlow automates the mechanical parts of molecular "
            "docking: downloading, preparing receptors and ligands, "
            "defining the search space, running AutoDock Vina and "
            "analysing the results.\n\n"
            "It does NOT make scientific decisions for you: protonation "
            "states, grid-box assumptions and success thresholds are "
            "recorded in every run's manifest and surfaced as warnings, "
            "so you can override them.\n\n"
            "This short tour shows where things live.  Full reference: "
            "USER_GUIDE.md."
        ),
    },
    {
        "title": "1. Pick a target",
        "where": "step bar: 'Target'",
        "body": (
            "Enter a PDB ID (e.g. 1HVR) and DockFlow fetches the "
            "structure and its co-crystal ligands from RCSB.\n\n"
            "The 'Ligands' step is where you add what to dock: the "
            "co-crystal ligand for redocking validation, PubChem names, "
            "or SMILES.\n\n"
            "Preparation options (engine, protonation, tautomers, salts) "
            "are visible on the 'Prepare' page - every choice ends up in "
            "manifest.json."
        ),
    },
    {
        "title": "2. Run the docking",
        "where": "step bar: 'Grid box' then 'Dock'",
        "body": (
            "The grid box can come from the co-crystal ligand (strongest "
            "assumption), active-site residues, or your explicit center "
            "and size - the page shows the provenance and warns when the "
            "assumption is weak.\n\n"
            "The 'Dock' page sets exhaustiveness, seed and backend.  A "
            "fixed seed makes runs reproducible; the progress bar and "
            "the activity log at the bottom show live status.\n\n"
            "Runs are checkpointed: an interrupted run resumes where it "
            "stopped."
        ),
    },
    {
        "title": "3. Interpret the results",
        "where": "step bar: 'Results'",
        "body": (
            "The Results page ranks poses by score and, for redocking, "
            "by symmetry-corrected RMSD to the crystal pose (2.0 A is "
            "the conventional - heuristic - success mark).\n\n"
            "Open the run folder for report.md, manifest.json and "
            "visualization/interactive.html - a self-contained 3D viewer "
            "(needs internet for the 3Dmol.js CDN).\n\n"
            "Remember: Vina scores are heuristic ranking functions, not "
            "binding free energies.  Interaction analysis is geometric "
            "contact counting - see docs/interaction_criteria.md."
        ),
    },
]


def should_show_tutorial(settings) -> bool:
    """True when the current tutorial version has not been seen yet.

    ``DOCKFLOW_NO_TUTORIAL=1`` disables it entirely (CI / headless
    automation), because a modal dialog without an event loop hangs.
    """
    import os

    if os.environ.get("DOCKFLOW_NO_TUTORIAL"):
        return False
    seen = settings.value(SEEN_KEY, 0, type=int)
    return int(seen) < TUTORIAL_VERSION


def mark_tutorial_seen(settings, version: int = TUTORIAL_VERSION) -> None:
    """Record that the tutorial was shown (dismissed counts as seen)."""
    settings.setValue(SEEN_KEY, int(version))
    settings.sync()


class TutorialDialog(QDialog):
    """A 4-page welcome tour; navigation via Next/Back/Finish buttons."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Welcome to DockFlow")
        self.setModal(True)
        self.setMinimumSize(560, 420)

        self._stack = QStackedWidget()
        for step in TUTORIAL_STEPS:
            self._stack.addWidget(self._build_page(step))
        self._page_indicator = QLabel()
        self._page_indicator.setAlignment(Qt.AlignmentFlag.AlignCenter)

        self._btn_back = QPushButton("Back")
        self._btn_next = QPushButton("Next")
        self._btn_back.clicked.connect(self._go_back)
        self._btn_next.clicked.connect(self._go_next)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.button(QDialogButtonBox.StandardButton.Close).clicked.connect(
            self.reject)

        nav = QVBoxLayout(self)
        nav.addWidget(self._stack, 1)
        nav.addWidget(self._page_indicator)
        row = QDialogButtonBox()
        row.addButton(self._btn_back,
                      QDialogButtonBox.ButtonRole.ActionRole)
        row.addButton(self._btn_next,
                      QDialogButtonBox.ButtonRole.ActionRole)
        row.addButton(buttons.button(
            QDialogButtonBox.StandardButton.Close),
            QDialogButtonBox.ButtonRole.RejectRole)
        nav.addWidget(row)
        self._update_nav()

    # ------------------------------------------------------------------ pages
    def _build_page(self, step: dict[str, str]) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(28, 24, 28, 24)
        title = QLabel(step["title"])
        title.setStyleSheet("font-size: 18px; font-weight: 700;")
        title.setWordWrap(True)
        where = QLabel(f"Where: {step['where']}")
        where.setStyleSheet("color: #666; font-style: italic;")
        body = QLabel(step["body"])
        body.setWordWrap(True)
        body.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(title)
        layout.addWidget(where)
        layout.addWidget(body, 1)
        return page

    def _update_nav(self) -> None:
        index = self._stack.currentIndex()
        total = self._stack.count()
        self._btn_back.setEnabled(index > 0)
        self._btn_next.setText("Finish" if index == total - 1 else "Next")
        self._page_indicator.setText(f"{index + 1} / {total}")

    def _go_back(self) -> None:
        self._stack.setCurrentIndex(max(0, self._stack.currentIndex() - 1))
        self._update_nav()

    def _go_next(self) -> None:
        if self._stack.currentIndex() == self._stack.count() - 1:
            self.accept()
            return
        self._stack.setCurrentIndex(self._stack.currentIndex() + 1)
        self._update_nav()

    # ------------------------------------------------------------------ API
    def page_texts(self) -> list[str]:
        """All page texts (smoke tests assert content presence)."""
        return [step["body"] for step in TUTORIAL_STEPS]


def maybe_show_tutorial(parent: QWidget, settings) -> bool:
    """Show the first-run tutorial once; True when it was displayed.

    Both finishing and dismissing mark the tutorial as seen - the user
    has seen the invitation either way.  Re-showing is always possible
    via Help > First-run tutorial.
    """
    if not should_show_tutorial(settings):
        return False
    dialog = TutorialDialog(parent)
    dialog.exec()
    mark_tutorial_seen(settings)
    return True
