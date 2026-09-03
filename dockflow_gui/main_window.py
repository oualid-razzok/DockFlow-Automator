"""Main window of the DockFlow-Automator desktop application.

A six-step workflow (Target -> Ligands -> Prepare -> Grid box -> Docking ->
Results) with a persistent log dock, non-blocking workers and an
interactive grid-box preview.  Heavy science happens in
:mod:`dockflow_core`; this module only orchestrates it.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from PyQt6.QtCore import QSettings, Qt
from PyQt6.QtGui import QColor, QKeySequence
from PyQt6.QtWidgets import (
    QComboBox,
    QDialog,
    QDockWidget,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSplitter,
    QStackedWidget,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)

from dockflow_core.config import get_config
from dockflow_core.docker_engine import detect_backends
from dockflow_core.models import LigandRecord, ProteinRecord
from dockflow_core.pdbio import parse_pdb
from dockflow_core.utils import get_logger

from .resources import ACCENT, DANGER, SUCCESS
from .threads import FunctionWorker, PipelineWorker
from .widgets import (
    ContactTableWidget,
    DockParamsWidget,
    GridBoxPreview,
    GridBoxWidget,
    LigandTableWidget,
    LogPanel,
    ResultsTableWidget,
    StepBar,
    make_caption,
)

logger = get_logger("gui")

_PREP_ENGINE_LABEL = {
    "auto": "auto (best available)",
    "openbabel": "OpenBabel python bindings",
    "openbabel-cli": "OpenBabel CLI (obabel)",
    "rdkit": "RDKit",
    "none": "none (no charges; testing only)",
}


@dataclass
class GuiState:
    """Everything the wizard accumulates between steps."""

    workdir: Path = field(default_factory=lambda: Path("dockflow_runs").resolve())
    target: ProteinRecord | None = None
    target_atoms: list = field(default_factory=list)
    ligands: list[LigandRecord] = field(default_factory=list)
    receptor_pdbqt: Path | None = None
    receptor_pdb: Path | None = None
    receptor_engine: str = ""
    grid_box: Any = None
    results: list = field(default_factory=list)
    selected_ligand: int | None = None
    prepared_dir: Path | None = None


class MainWindow(QMainWindow):
    """The DockFlow application window."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("DockFlow-Automator")
        self.resize(1240, 860)
        self.config = get_config()
        self.state = GuiState(workdir=self.config.workdir.resolve()
                              if self.config.workdir.exists() or True else Path("."))
        self._workers: list = []
        self.settings = QSettings("DockFlow", "DockFlow-Automator")

        self._build_ui()
        self._build_menus()
        self._restore_geometry()
        self._refresh_engine_labels()
        logger.info("GUI ready (workdir=%s)", self.state.workdir)

    # ------------------------------------------------------------------ layout
    def _build_ui(self) -> None:
        self.stepbar = StepBar()
        self.stepbar.stepActivated.connect(self._goto_step)

        self.pages = QStackedWidget()
        self.page_target = self._build_target_page()
        self.page_ligands = self._build_ligands_page()
        self.page_prepare = self._build_prepare_page()
        self.page_gridbox = self._build_gridbox_page()
        self.page_dock = self._build_dock_page()
        self.page_results = self._build_results_page()
        for page in (
            self.page_target, self.page_ligands, self.page_prepare,
            self.page_gridbox, self.page_dock, self.page_results,
        ):
            self.pages.addWidget(page)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self.stepbar)
        splitter.addWidget(self.pages)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([210, 1030])
        self.setCentralWidget(splitter)

        # log dock
        self.log_panel = LogPanel()
        self.log_dock = QDockWidget("Activity log", self)
        self.log_dock.setWidget(self.log_panel)
        self.log_dock.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetMovable
            | QDockWidget.DockWidgetFeature.DockWidgetClosable
        )
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.log_dock)

        # status bar
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.progress.setMaximumWidth(260)
        self.progress.setVisible(False)
        self.status_label = QLabel("ready")
        self.statusBar().addWidget(self.status_label, 1)
        self.statusBar().addPermanentWidget(self.progress)
        self.stepbar.set_state(0, "active")

    def _wrap_page(self, title: str, caption: str, body: QWidget) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(18, 14, 18, 14)
        header = QLabel(title)
        header.setStyleSheet(f"font-size: 17px; font-weight: 700; color: {ACCENT};")
        layout.addWidget(header)
        layout.addWidget(make_caption(caption))
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setWidget(body)
        layout.addWidget(scroll, 1)
        return page

    # ------------------------------------------------------------------ pages
    def _build_target_page(self) -> QWidget:
        body = QWidget()
        form = QFormLayout(body)
        self.target_kind = QComboBox()
        self.target_kind.addItems(["PDB id", "UniProt accession", "Local file"])
        self.target_edit = QComboBox()
        self.target_edit.setEditable(True)
        self.target_edit.setPlaceholderText("e.g. 1HVR")
        self.target_edit.addItem("1HVR")
        self.target_edit.addItem("1AKO")
        self.workdir_edit = QComboBox()
        self.workdir_edit.setEditable(True)
        self.workdir_edit.addItem(str(self.state.workdir))
        browse_target = QPushButton("Browse…")
        browse_target.setProperty("secondary", True)
        browse_workdir = QPushButton("Browse…")
        browse_workdir.setProperty("secondary", True)
        self.fetch_target_btn = QPushButton("Fetch target")
        self.target_info = QLabel("")
        self.target_info.setProperty("muted", True)
        self.target_info.setWordWrap(True)
        target_row = QHBoxLayout()
        target_row.addWidget(self.target_kind)
        target_row.addWidget(self.target_edit, 1)
        target_row.addWidget(browse_target)
        workdir_row = QHBoxLayout()
        workdir_row.addWidget(self.workdir_edit, 1)
        workdir_row.addWidget(browse_workdir)
        form.addRow("Target source", target_row)
        form.addRow("Working directory", workdir_row)
        form.addRow(self.fetch_target_btn)
        form.addRow(self.target_info)
        browse_target.clicked.connect(self._browse_target_file)
        browse_workdir.clicked.connect(self._browse_workdir)
        self.fetch_target_btn.clicked.connect(self._fetch_target)
        self.target_kind.currentTextChanged.connect(self._target_kind_changed)
        return self._wrap_page(
            "1 - Target",
            "Fetch a crystal structure from the RCSB PDB, resolve a UniProt "
            "accession to its best experimental structure (or AlphaFold model), "
            "or load a local PDB file.",
            body,
        )

    def _build_ligands_page(self) -> QWidget:
        body = QWidget()
        layout = QVBoxLayout(body)
        layout.setContentsMargins(0, 0, 0, 0)
        self.ligand_table = LigandTableWidget()
        layout.addWidget(LigandTableWidget.with_buttons(self.ligand_table))
        self.ligand_note = QLabel("")
        self.ligand_note.setProperty("muted", True)
        self.ligand_note.setWordWrap(True)
        layout.addWidget(self.ligand_note)
        self.ligand_table.addSmilesRequested.connect(self._add_ligand_smiles)
        self.ligand_table.addFilesRequested.connect(self._add_ligand_files)
        self.ligand_table.addPubchemRequested.connect(self._add_ligand_pubchem)
        self.ligand_table.addZincRequested.connect(self._add_ligand_zinc)
        self.ligand_table.addCocystalRequested.connect(self._add_ligand_cocystal)
        self.ligand_table.removeRequested.connect(self._remove_ligands)
        return self._wrap_page(
            "2 - Ligands",
            "Collect ligands from SMILES, local files (SDF/MOL2/PDB/PDBQT), "
            "PubChem, ZINC22, or straight from the co-crystallized ligand of "
            "the loaded target.",
            body,
        )

    def _build_prepare_page(self) -> QWidget:
        body = QWidget()
        layout = QVBoxLayout(body)
        layout.setContentsMargins(0, 0, 0, 0)
        # receptor options
        rec_box = QWidget()
        rec_form = QFormLayout(rec_box)
        self.rec_chains = QComboBox()
        self.rec_chains.setEditable(True)
        self.rec_chains.setPlaceholderText("empty = all chains")
        self.rec_altloc = QComboBox()
        self.rec_altloc.addItems(["best", "A", "B", "keep all"])
        self.rec_engine = QComboBox()
        self.rec_engine.addItems(list(_PREP_ENGINE_LABEL.keys()))
        self.rec_engine.setCurrentIndex(0)
        self.rec_keep_water = QPushButton("Keep waters")
        self.rec_keep_water.setCheckable(True)
        self.rec_keep_hetero = QPushButton("Keep all HETATM")
        self.rec_keep_hetero.setCheckable(True)
        self.rec_keep_metals = QComboBox()
        self.rec_keep_metals.setEditable(True)
        self.rec_keep_metals.setPlaceholderText("e.g. ZN, MG (kept even without HETATM)")
        rec_form.addRow("chains", self.rec_chains)
        rec_form.addRow("altLoc policy", self.rec_altloc)
        rec_form.addRow("preparation engine", self.rec_engine)
        toggles = QHBoxLayout()
        toggles.addWidget(self.rec_keep_water)
        toggles.addWidget(self.rec_keep_hetero)
        rec_form.addRow(toggles)
        rec_form.addRow("always keep residues", self.rec_keep_metals)
        self.prep_receptor_btn = QPushButton("Prepare receptor (PDBQT)")
        rec_form.addRow(self.prep_receptor_btn)
        self.receptor_status = QLabel("")
        self.receptor_status.setProperty("muted", True)
        self.receptor_status.setWordWrap(True)
        rec_form.addRow(self.receptor_status)
        # ligand options
        lig_box = QWidget()
        lig_form = QFormLayout(lig_box)
        self.lig_minimize = QPushButton("MMFF minimisation")
        self.lig_minimize.setCheckable(True)
        self.lig_minimize.setChecked(True)
        self.lig_protonate = QPushButton("dimorphite-dl protonation")
        self.lig_protonate.setCheckable(False)
        self.lig_protonate.setToolTip("requires the optional dimorphite-dl package")
        lig_form.addRow(self.lig_minimize)
        lig_form.addRow(self.lig_protonate)
        self.prep_ligands_btn = QPushButton("Prepare all ligands (Meeko)")
        lig_form.addRow(self.prep_ligands_btn)
        self.ligand_prep_status = QLabel("")
        self.ligand_prep_status.setProperty("muted", True)
        self.ligand_prep_status.setWordWrap(True)
        lig_form.addRow(self.ligand_prep_status)
        for widget, title, caption in (
            (rec_box, "Receptor (prepare_receptor4.py logic)",
             "Chain filtering, alternate locations, water/heteroatom removal, "
             "hydrogen addition, Gasteiger charges, non-polar H merging, AD4 types."),
            (lig_box, "Ligands (prepare_ligand4.py logic)",
             "Meeko + RDKit: sanitisation, largest fragment, 3D embedding "
             "(ETKDG), MMFF minimisation, torsion trees."),
        ):
            from PyQt6.QtWidgets import QGroupBox

            group = QGroupBox(title)
            inner = QVBoxLayout(group)
            inner.addWidget(widget)
            inner.addWidget(make_caption(caption))
            layout.addWidget(group)
        self.prep_receptor_btn.clicked.connect(self._prepare_receptor)
        self.prep_ligands_btn.clicked.connect(self._prepare_ligands)
        return self._wrap_page(
            "3 - Preparation",
            "Turn the raw target and ligands into AutoDock-ready PDBQT files.",
            body,
        )

    def _build_gridbox_page(self) -> QWidget:
        body = QWidget()
        layout = QHBoxLayout(body)
        layout.setContentsMargins(0, 0, 0, 0)
        self.gridbox_widget = GridBoxWidget()
        self.gridbox_preview = GridBoxPreview()
        layout.addWidget(self.gridbox_widget)
        layout.addWidget(self.gridbox_preview, 1)
        self.gridbox_widget.boxChanged.connect(self._on_box_changed)
        self.gridbox_widget.computeFromLigandRequested.connect(
            self._compute_box_from_cocystal)
        self.gridbox_widget.computeFromResiduesRequested.connect(
            self._compute_box_from_residues)
        return self._wrap_page(
            "4 - Grid box",
            "Define the Vina search space.  Derive it automatically from the "
            "co-crystallized ligand or active-site residues, or set it manually; "
            "the preview updates live (drag to rotate).",
            body,
        )

    def _build_dock_page(self) -> QWidget:
        body = QWidget()
        layout = QVBoxLayout(body)
        layout.setContentsMargins(0, 0, 0, 0)
        from PyQt6.QtWidgets import QGroupBox

        params_group = QGroupBox("Run")
        params_layout = QHBoxLayout(params_group)
        self.dock_params = DockParamsWidget()
        self.dock_run_btn = QPushButton("Start docking")
        self.dock_cancel_btn = QPushButton("Cancel")
        self.dock_cancel_btn.setProperty("danger", True)
        self.dock_cancel_btn.setEnabled(False)
        side = QVBoxLayout()
        side.addWidget(self.dock_run_btn)
        side.addWidget(self.dock_cancel_btn)
        side.addStretch(1)
        self.engine_label = QLabel("")
        self.engine_label.setProperty("muted", True)
        self.engine_label.setWordWrap(True)
        side.addWidget(self.engine_label)
        params_layout.addWidget(self.dock_params, 1)
        params_layout.addLayout(side)
        layout.addWidget(params_group)
        self.dock_summary = QTableWidget(0, 4)
        self.dock_summary.setHorizontalHeaderLabels(
            ["ligand", "best affinity", "poses", "runtime (s)"]
        )
        self.dock_summary.verticalHeader().setVisible(False)
        self.dock_summary.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        layout.addWidget(self.dock_summary, 1)
        self.dock_run_btn.clicked.connect(self._run_docking)
        self.dock_cancel_btn.clicked.connect(self._cancel_docking)
        return self._wrap_page(
            "5 - Docking",
            "Run AutoDock Vina on every prepared ligand inside the grid box. "
            "The table fills in as ligands finish.",
            body,
        )

    def _build_results_page(self) -> QWidget:
        body = QWidget()
        layout = QVBoxLayout(body)
        layout.setContentsMargins(0, 0, 0, 0)
        top = QHBoxLayout()
        self.results_ligand = QComboBox()
        top.addWidget(QLabel("ligand:"))
        top.addWidget(self.results_ligand, 1)
        self.analyze_btn = QPushButton("Analyse interactions")
        self.render_btn = QPushButton("Render PNG")
        self.open_pymol_btn = QPushButton("Open in PyMOL")
        self.export_btn = QPushButton("Export CSV…")
        self.open_folder_btn = QPushButton("Open run folder")
        for button in (self.analyze_btn, self.render_btn, self.open_pymol_btn,
                       self.export_btn, self.open_folder_btn):
            button.setProperty("secondary", True)
            top.addWidget(button)
        layout.addLayout(top)
        middle = QHBoxLayout()
        self.results_table = ResultsTableWidget()
        self.contact_table = ContactTableWidget()
        middle.addWidget(self.results_table, 1)
        middle.addWidget(self.contact_table, 1)
        layout.addLayout(middle, 1)
        self.preview_label = QLabel("no render yet")
        self.preview_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.preview_label.setMinimumHeight(280)
        self.preview_label.setStyleSheet(
            f"background: #ffffff; border: 1px dashed {DANGER}; color: #9aa0a6;"
        )
        layout.addWidget(self.preview_label, 2)
        self.results_ligand.currentIndexChanged.connect(self._on_results_ligand_changed)
        self.analyze_btn.clicked.connect(self._analyze_current)
        self.render_btn.clicked.connect(self._render_current)
        self.open_pymol_btn.clicked.connect(self._open_pymol_current)
        self.export_btn.clicked.connect(self._export_csv)
        self.open_folder_btn.clicked.connect(self._open_run_folder)
        return self._wrap_page(
            "6 - Results",
            "Inspect poses and affinities, analyse receptor interactions, "
            "render ray-traced images or open the scene in PyMOL.",
            body,
        )

    # ------------------------------------------------------------------ menus
    def _build_menus(self) -> None:
        from PyQt6.QtGui import QAction

        file_menu = self.menuBar().addMenu("&File")
        new_run = QAction("&New run…", self)
        new_run.setShortcut(QKeySequence("Ctrl+N"))
        new_run.triggered.connect(self._new_run)
        file_menu.addAction(new_run)
        settings_action = QAction("&Settings…", self)
        settings_action.triggered.connect(self._open_settings)
        file_menu.addAction(settings_action)
        file_menu.addSeparator()
        quit_action = QAction("E&xit", self)
        quit_action.setShortcut(QKeySequence("Ctrl+Q"))
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

        tools_menu = self.menuBar().addMenu("&Tools")
        run_pipeline = QAction("Run automated &pipeline (YAML)…", self)
        run_pipeline.triggered.connect(self._run_pipeline_yaml)
        tools_menu.addAction(run_pipeline)
        env_report = QAction("&Environment report", self)
        env_report.triggered.connect(self._show_env_report)
        tools_menu.addAction(env_report)

        help_menu = self.menuBar().addMenu("&Help")
        about = QAction("&About DockFlow", self)
        about.triggered.connect(self._show_about)
        help_menu.addAction(about)

    # ------------------------------------------------------------------ misc UI
    def _goto_step(self, index: int) -> None:
        self.pages.setCurrentIndex(index)

    def _set_status(self, message: str) -> None:
        self.status_label.setText(message)

    def _set_busy(self, busy: bool, button=None) -> None:
        if button is not None:
            button.setEnabled(not busy)
        self.fetch_target_btn.setEnabled(not busy or button is not self.fetch_target_btn)

    def _log(self, message: str, level: str = "info") -> None:
        self.log_panel.append_log(message, level=level)
        logger.info("%s", message)

    def _restore_geometry(self) -> None:
        geometry = self.settings.value("window/geometry")
        if geometry is not None:
            self.restoreGeometry(geometry)
        state = self.settings.value("window/state")
        if state is not None:
            self.restoreState(state)

    def closeEvent(self, event) -> None:  # noqa: N802 - Qt naming
        self.settings.setValue("window/geometry", self.saveGeometry())
        self.settings.setValue("window/state", self.saveState())
        if any(w.isRunning() for w in self._workers):
            answer = QMessageBox.question(
                self, "Workers running",
                "Background jobs are still running. Quit anyway?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            if answer != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
        for worker in self._workers:
            worker.quit()
        event.accept()

    # ------------------------------------------------------------------ workers
    def _start_worker(self, fn, on_result=None, on_progress=None, on_log=None,
                      busy_button=None, on_finished=None, worker_class=None,
                      cancel_event=None):
        worker_cls = worker_class or FunctionWorker
        kwargs = {}
        if cancel_event is not None:
            kwargs["cancel_event"] = cancel_event
        worker = worker_cls(fn, **kwargs)
        worker.signals.log.connect(
            lambda line, lvl="info": self._log(line, lvl))
        if on_progress:
            worker.signals.progress.connect(on_progress)
        if on_result:
            worker.signals.result.connect(on_result)
        worker.signals.error.connect(self._on_worker_error)
        if busy_button is not None:
            busy_button.setEnabled(False)
            worker.signals.finished.connect(lambda: busy_button.setEnabled(True))
        if on_finished:
            worker.signals.finished.connect(on_finished)
        worker.signals.finished.connect(
            lambda w=worker: self._workers.remove(w) if w in self._workers else None)
        self._workers.append(worker)
        worker.start()
        return worker

    def _on_worker_error(self, message: str, traceback_text: str) -> None:
        self._log(message, level="error")
        self._set_status(f"error: {message}")
        QMessageBox.critical(self, "DockFlow error", message)
        logger.debug("%s", traceback_text)

    def _show_progress(self, percent: int, message: str) -> None:
        self.progress.setVisible(True)
        self.progress.setValue(percent)
        self._set_status(message)
        if percent >= 100:
            self.progress.setVisible(False)

    # ------------------------------------------------------------------ step 1
    def _target_kind_changed(self, text: str) -> None:
        if text == "PDB id":
            self.target_edit.setPlaceholderText("e.g. 1HVR")
        elif text == "UniProt accession":
            self.target_edit.setPlaceholderText("e.g. P29978")
        else:
            self.target_edit.setPlaceholderText("path to a .pdb file")
        self.fetch_target_btn.setText("Load local file" if text == "Local file"
                                      else "Fetch target")

    def _browse_target_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose target structure", str(self.state.workdir),
            "PDB files (*.pdb *.ent);;PDBQT files (*.pdbqt);;All files (*)")
        if path:
            self.target_kind.setCurrentText("Local file")
            self.target_edit.setCurrentText(path)

    def _browse_workdir(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Choose working directory", str(self.state.workdir))
        if path:
            self.workdir_edit.setCurrentText(path)
            self.state.workdir = Path(path)

    def _fetch_target(self) -> None:
        kind = self.target_kind.currentText()
        value = self.target_edit.currentText().strip()
        if not value:
            QMessageBox.warning(self, "Missing input", "Enter a target id or file path.")
            return
        self.state.workdir = Path(self.workdir_edit.currentText().strip())
        raw_dir = self.state.workdir / "raw"
        if kind == "Local file":
            self._load_local_target(value)
            return

        def job():
            from dockflow_core.downloader import PDBDownloader

            downloader = PDBDownloader(cache_dir=self.config.cache_dir)
            if kind == "UniProt":
                return downloader.fetch_alphafold_model(value, raw_dir)
            return downloader.fetch_structure(value, raw_dir)

        self._log(f"fetching target {value!r}…")
        self._start_worker(job, on_result=self._on_target_loaded,
                           busy_button=self.fetch_target_btn)

    def _load_local_target(self, path: str) -> None:
        target_path = Path(path)
        if not target_path.is_file():
            QMessageBox.warning(self, "Not found", f"File not found:\n{path}")
            return

        def job():
            from dockflow_core.models import ProteinRecord

            atoms = parse_pdb(target_path)
            if not atoms:
                raise ValueError(f"no atoms parsed from {target_path}")
            codes: list[str] = []
            seen: set[str] = set()
            for atom in atoms:
                name = atom.resname.strip()
                if atom.is_polymer or atom.is_water or not name.isalnum():
                    continue
                if name not in seen:
                    seen.add(name)
                    codes.append(name)
            return ProteinRecord(identifier=target_path.stem, source="file",
                                 path=target_path, ligand_codes=codes)

        self._start_worker(job, on_result=self._on_target_loaded,
                           busy_button=self.fetch_target_btn)

    def _on_target_loaded(self, record: ProteinRecord) -> None:
        self.state.target = record
        self.state.target_atoms = parse_pdb(record.path) if record.path else []
        title = record.title or "(no title)"
        resolution = f", {record.resolution:.2f} A" if record.resolution else ""
        self.target_info.setText(
            f"target {record.identifier} ({record.source}) - {title}{resolution}\n"
            f"{len(self.state.target_atoms)} atoms"
            + (f", ligands: {', '.join(record.ligand_codes)}"
               if record.ligand_codes else "")
        )
        self._log(f"target loaded: {record.identifier} "
                  f"({len(self.state.target_atoms)} atoms)", "ok")
        self.stepbar.set_state(0, "done")
        self.stepbar.set_state(1, "active")
        self._goto_step(1)
        if record.ligand_codes:
            self.ligand_note.setText(
                "Co-crystallized ligands available: "
                + ", ".join(record.ligand_codes)
                + "  (use 'Co-crystal ligand…' to dock them)"
            )

    # ------------------------------------------------------------------ step 2
    def _add_ligand_smiles(self) -> None:
        from PyQt6.QtWidgets import QDialogButtonBox

        dialog = QDialog(self)
        dialog.setWindowTitle("Add ligand from SMILES")
        form = QFormLayout(dialog)
        smiles = QLineEdit()
        smiles.setPlaceholderText("CC(=O)Oc1ccccc1C(=O)O  (aspirin)")
        name = QLineEdit()
        name.setPlaceholderText("optional identifier")
        form.addRow("SMILES", smiles)
        form.addRow("name", name)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        form.addRow(buttons)
        if dialog.exec() != QDialog.DialogCode.Accepted or not smiles.text().strip():
            return
        identifier = name.text().strip() or f"smiles_{len(self.state.ligands) + 1}"
        self.state.ligands.append(
            LigandRecord(identifier=identifier, source="smiles",
                         value=smiles.text().strip(), status="pending")
        )
        self._refresh_ligand_table()

    def _add_ligand_files(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Choose ligand files", str(self.state.workdir),
            "Ligands (*.sdf *.sd *.mol2 *.pdb *.pdbqt *.smi);;All files (*)")
        for path in paths:
            self.state.ligands.append(
                LigandRecord(identifier=Path(path).stem, source="file",
                             value=path, path=Path(path),
                             status="prepared" if path.endswith(".pdbqt") else "pending")
            )
        if paths:
            self._refresh_ligand_table()

    def _add_ligand_pubchem(self) -> None:
        text, ok = QInputDialog.getText(
            self, "PubChem ligand", "Name, CID or SMILES:")
        if not ok or not text.strip():
            return

        def job():
            from dockflow_core.downloader import LigandDownloader

            return LigandDownloader().fetch_pubchem(text.strip(),
                                                    self.state.workdir / "raw")

        self._log(f"downloading {text.strip()} from PubChem…")
        self._start_worker(job, on_result=self._on_ligand_downloaded)

    def _add_ligand_zinc(self) -> None:
        text, ok = QInputDialog.getText(
            self, "ZINC ligand", "ZINC22 identifier (e.g. ZINC000000000001):")
        if not ok or not text.strip():
            return

        def job():
            from dockflow_core.downloader import LigandDownloader

            return LigandDownloader().fetch_zinc(text.strip(),
                                                 self.state.workdir / "raw")

        self._log(f"downloading {text.strip()} from ZINC22…")
        self._start_worker(job, on_result=self._on_ligand_downloaded)

    def _add_ligand_cocystal(self) -> None:
        record = self.state.target
        if record is None or not record.ligand_codes:
            QMessageBox.information(
                self, "No co-crystal ligand",
                "Load a target with a co-crystallized ligand first.")
            return
        code, ok = QInputDialog.getItem(
            self, "Co-crystal ligand", "Residue code:", record.ligand_codes, 0, False)
        if not ok:
            return

        def job():
            from dockflow_core.downloader import PDBDownloader

            path = PDBDownloader().fetch_ligand(code, self.state.workdir / "raw")
            return LigandRecord(identifier=code.lower(), source="pdb_ligand",
                                value=code, path=path, status="downloaded")

        self._log(f"downloading ligand {code} from RCSB…")
        self._start_worker(job, on_result=self._on_ligand_downloaded)

    def _on_ligand_downloaded(self, ligand: LigandRecord) -> None:
        self.state.ligands.append(ligand)
        self._log(f"ligand {ligand.identifier} downloaded -> {ligand.path}", "ok")
        self._refresh_ligand_table()

    def _remove_ligands(self) -> None:
        rows = self.ligand_table.selected_rows()
        for row in sorted(rows, reverse=True):
            if 0 <= row < len(self.state.ligands):
                self.state.ligands.pop(row)
        self._refresh_ligand_table()

    def _refresh_ligand_table(self) -> None:
        self.ligand_table.set_ligands(self.state.ligands)
        if self.state.ligands:
            self.stepbar.set_state(1, "done" if any(
                ligand.status in ("prepared", "docked") for ligand in self.state.ligands
            ) else "active")

    # ------------------------------------------------------------------ step 3
    def _receptor_options(self):
        from dockflow_core.preparator import ReceptorPrepOptions

        chains_text = self.rec_chains.currentText().strip()
        altloc = self.rec_altloc.currentText()
        keep_res = [
            token.strip().upper()
            for token in self.rec_keep_metals.currentText().split(",")
            if token.strip()
        ]
        return ReceptorPrepOptions(
            chains=[c.strip() for c in chains_text.split(",")] if chains_text else None,
            keep_water=self.rec_keep_water.isChecked(),
            keep_hetero=self.rec_keep_hetero.isChecked(),
            keep_resnames=keep_res,
            altloc="" if altloc == "keep all" else altloc,
            engine=self.rec_engine.currentText(),
        )

    def _prepare_receptor(self) -> None:
        target = self.state.target
        if target is None or target.path is None:
            QMessageBox.warning(self, "No target", "Load a target structure first (step 1).")
            return
        options = self._receptor_options()
        out_dir = self.state.workdir / "prepared"
        input_path = target.path

        def job(on_log=None):
            from dockflow_core.preparator import ReceptorPreparator

            result = ReceptorPreparator(options).prepare(input_path, out_dir)
            if on_log:
                for warning in result.warnings:
                    on_log(warning)
            return result

        self._log("preparing receptor…")
        self._start_worker(job, on_result=self._on_receptor_prepared,
                           busy_button=self.prep_receptor_btn)

    def _on_receptor_prepared(self, result) -> None:
        self.state.receptor_pdbqt = result.pdbqt_path
        self.state.receptor_pdb = result.pdb_path
        self.state.receptor_engine = result.engine
        self.state.prepared_dir = result.pdbqt_path.parent if result.pdbqt_path else None
        self.receptor_status.setText(
            f"PDBQT: {result.pdbqt_path}\nengine: {result.engine} | "
            f"{result.atoms_in} -> {result.atoms_out} atoms | +{result.hydrogens_added} H"
        )
        self._log(f"receptor ready: {result.pdbqt_path.name} "
                  f"(engine {result.engine}, {result.atoms_out} atoms)", "ok")
        self.stepbar.set_state(2, "done")
        self.stepbar.set_state(3, "active")
        self._update_gridbox_preview()

    def _prepare_ligands(self) -> None:
        if not self.state.ligands:
            QMessageBox.warning(self, "No ligands", "Add at least one ligand (step 2).")
            return
        pending = [
            (index, ligand)
            for index, ligand in enumerate(self.state.ligands)
            if ligand.status in ("pending", "downloaded") and ligand.source != "pdbqt"
        ]
        if not pending:
            self._log("all ligands are already prepared", "warn")
            return
        out_dir = self.state.workdir / "prepared"

        def job(on_log=None):
            from dockflow_core.preparator import LigandPreparator, LigandPrepOptions

            options = LigandPrepOptions(minimize=self.lig_minimize.isChecked(),
                                         protonate=False)
            preparator = LigandPreparator(options)
            outcomes = []
            for index, ligand in pending:
                try:
                    source = ligand.value if ligand.source == "smiles" else ligand.path
                    prep = preparator.prepare(source, out_dir, ligand.identifier)
                    outcomes.append((index, prep, None))
                    if on_log:
                        on_log(f"{ligand.identifier}: {prep.num_atoms} atoms, "
                               f"{prep.num_rotatable_bonds} rotatable bonds")
                except Exception as exc:  # noqa: BLE001
                    outcomes.append((index, None, str(exc)))
                    if on_log:
                        on_log(f"{ligand.identifier}: FAILED {exc}")
            return outcomes

        self._log(f"preparing {len(pending)} ligand(s) with Meeko…")
        self._start_worker(job, on_result=self._on_ligands_prepared,
                           busy_button=self.prep_ligands_btn)

    def _on_ligands_prepared(self, outcomes) -> None:
        ok = 0
        for index, prep, error in outcomes:
            ligand = self.state.ligands[index]
            if error is None and prep is not None:
                ligand.pdbqt_path = prep.pdbqt_path
                ligand.status = "prepared"
                ligand.num_rotatable_bonds = prep.num_rotatable_bonds
                ok += 1
            else:
                ligand.status = "error"
                ligand.error = error
        self.ligand_prep_status.setText(
            f"prepared {ok}/{len(outcomes)} ligands -> "
            f"{self.state.workdir / 'prepared'}"
        )
        self._refresh_ligand_table()
        self._log(f"ligand preparation finished: {ok}/{len(outcomes)} ok",
                  "ok" if ok else "error")
        if ok:
            self.stepbar.set_state(2, "done")

    # ------------------------------------------------------------------ step 4
    def _on_box_changed(self, box) -> None:
        self.state.grid_box = box
        self._update_gridbox_preview()

    def _update_gridbox_preview(self) -> None:
        atoms = []
        colors: list[QColor] = []
        if self.state.receptor_pdb and self.state.receptor_pdb.is_file():
            receptor_atoms = [
                a for a in parse_pdb(self.state.receptor_pdb) if not a.is_hydrogen
            ]
            stride = max(1, len(receptor_atoms) // 1500)
            atoms.extend(receptor_atoms[::stride])
            colors.extend([QColor("#9aa0a6")] * len(receptor_atoms[::stride]))
        for ligand in self.state.ligands:
            if ligand.pdbqt_path and ligand.pdbqt_path.is_file():
                from dockflow_core.pdbio import parse_pdbqt

                lig_atoms = parse_pdbqt(ligand.pdbqt_path).atoms
                atoms.extend(lig_atoms)
                colors.extend([QColor("#f59e0b")] * len(lig_atoms))
        if atoms:
            points = np.array([[a.x, a.y, a.z] for a in atoms])
        else:
            points = np.zeros((0, 3))
        self.gridbox_preview.set_data(points, colors, self.state.grid_box)

    def _compute_box_from_cocystal(self) -> None:
        target = self.state.target
        if target is None or target.path is None or not target.ligand_codes:
            QMessageBox.information(
                self, "No co-crystal ligand",
                "Load a target containing a co-crystallized ligand first.")
            return
        from dockflow_core.gridbox import box_from_pocket

        try:
            box = box_from_pocket(target.path, target.ligand_codes[0],
                                  padding=self.gridbox_widget.padding.value())
            self.gridbox_widget.set_box(box)
            self._log(f"grid box from ligand {target.ligand_codes[0]}: "
                      f"center {tuple(round(v, 1) for v in box.center)}", "ok")
        except (ValueError, Exception) as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Grid box failed", str(exc))

    def _compute_box_from_residues(self) -> None:
        if not self.state.target_atoms:
            QMessageBox.information(self, "No target", "Load a target structure first.")
            return
        text, ok = QInputDialog.getText(
            self, "Active-site residues",
            "chain (empty = any) and residue numbers.\nExample:  A : 32,48,84,90")
        if not ok:
            return
        chain = ""
        residues_text = text
        if ":" in text:
            chain, residues_text = (part.strip() for part in text.split(":", 1))
        residues = [r for r in (t.strip() for t in residues_text.split(",")) if r]
        if not residues:
            return
        from dockflow_core.gridbox import box_from_residues

        try:
            box = box_from_residues(self.state.target_atoms, chain or None,
                                    [int(r) for r in residues],
                                    padding=self.gridbox_widget.padding.value())
            self.gridbox_widget.set_box(box)
            self._log(f"grid box from residues {residues}: "
                      f"center {tuple(round(v, 1) for v in box.center)}", "ok")
        except ValueError as exc:
            QMessageBox.warning(self, "Grid box failed", str(exc))

    # ------------------------------------------------------------------ step 5
    def _refresh_engine_labels(self) -> None:
        try:
            backends = detect_backends(self.config.vina_exec, self.config.smina_exec)
        except Exception:  # noqa: BLE001
            backends = []
        lines = ["Vina backends detected:"]
        for backend in backends:
            mark = "✓" if backend.available else "✗"
            lines.append(f"  {mark} {backend.name}: {backend.version or backend.detail}")
        self.engine_label.setText("\n".join(lines))
        self._log(" / ".join(
            f"{b.name}={'ok' if b.available else 'missing'}" for b in backends))

    def _run_docking(self) -> None:
        if not self.state.receptor_pdbqt or not self.state.receptor_pdbqt.is_file():
            QMessageBox.warning(self, "Receptor missing", "Prepare the receptor first (step 3).")
            return
        prepared = [ligand for ligand in self.state.ligands if ligand.pdbqt_path]
        if not prepared:
            QMessageBox.warning(self, "Ligands missing", "Prepare ligands first (step 3).")
            return
        if self.state.grid_box is None:
            QMessageBox.warning(self, "Grid box missing", "Define the grid box first (step 4).")
            return
        params = self.dock_params.config_dict()
        receptor = self.state.receptor_pdbqt
        out_dir = self.state.workdir / "docking"
        box = self.state.grid_box
        ligands_snapshot = list(prepared)

        import threading

        cancel_event = threading.Event()

        def job(on_progress=None, cancel_event=None):
            from dockflow_core.docker_engine import VinaConfig, VinaEngine

            config = VinaConfig.from_gridbox(
                box,
                scoring=params["scoring"],
                exhaustiveness=params["exhaustiveness"],
                num_modes=params["num_modes"],
                refine=params["refine"],
                seed=params["seed"],
                cpu=params["cpu"],
                energy_range=params["energy_range"],
            )
            engine = VinaEngine(config, backend="auto", workdir=out_dir,
                                vina_exec=self.config.vina_exec,
                                smina_exec=self.config.smina_exec)
            return engine.dock_batch(
                receptor,
                [ligand.pdbqt_path for ligand in ligands_snapshot],
                out_dir=out_dir,
                ligand_records=ligands_snapshot,
                progress=on_progress,
                stop_event=cancel_event,
            )

        self.dock_summary.setRowCount(0)
        self.dock_run_btn.setEnabled(False)
        self.dock_cancel_btn.setEnabled(True)
        self._log(f"docking {len(ligands_snapshot)} ligand(s) "
                  f"(exhaustiveness={params['exhaustiveness']})…")
        self._dock_worker = self._start_worker(
            job,
            on_result=self._on_docking_finished,
            on_progress=self._show_progress,
            busy_button=self.dock_run_btn,
            cancel_event=cancel_event,
        )
        self._dock_cancel_event = cancel_event
        self.stepbar.set_state(4, "active")

    def _cancel_docking(self) -> None:
        event = getattr(self, "_dock_cancel_event", None)
        if event is not None:
            event.set()
        worker = getattr(self, "_dock_worker", None)
        if worker is not None:
            worker.cancel()
        self._log("cancellation requested…", "warn")
        self.dock_cancel_btn.setEnabled(False)

    def _on_docking_finished(self, results) -> None:
        self.state.results = results
        self.dock_cancel_btn.setEnabled(False)
        self.dock_run_btn.setEnabled(True)
        self.dock_summary.setRowCount(len(results))
        from PyQt6.QtWidgets import QTableWidgetItem as Item

        for row, result in enumerate(results):
            affinity = result.best_affinity
            values = (
                result.ligand_name,
                f"{affinity:.2f}" if affinity is not None else "n/a",
                str(len(result.poses)),
                f"{result.runtime:.1f}",
            )
            for column, text in enumerate(values):
                item = Item(text)
                if column == 1 and affinity is not None:
                    item.setForeground(
                        QColor(SUCCESS) if affinity <= -7 else
                        QColor("#d97706") if affinity <= -5 else QColor("#6b7280")
                    )
                self.dock_summary.setItem(row, column, item)
        self.dock_summary.resizeColumnsToContents()
        ok_count = sum(1 for r in results if r.ok)
        self._log(f"docking finished: {ok_count}/{len(results)} ligands docked",
                  "ok" if ok_count else "error")
        if ok_count:
            self.stepbar.set_state(4, "done")
            self.stepbar.set_state(5, "active")
            self._populate_results_page()
            self._goto_step(5)

    # ------------------------------------------------------------------ step 6
    def _populate_results_page(self) -> None:
        ok_results = [r for r in self.state.results if r.ok]
        self.results_ligand.blockSignals(True)
        self.results_ligand.clear()
        for result in sorted(ok_results, key=lambda r: r.best_affinity or 0):
            self.results_ligand.addItem(
                f"{result.ligand_name}  ({result.best_affinity:.2f} kcal/mol)"
            )
        self.results_ligand.blockSignals(False)
        if ok_results:
            self._on_results_ligand_changed(0)

    def _current_result(self):
        index = self.results_ligand.currentIndex()
        ok_results = [r for r in self.state.results if r.ok]
        ok_results.sort(key=lambda r: r.best_affinity or 0)
        if 0 <= index < len(ok_results):
            return ok_results[index]
        return None

    def _on_results_ligand_changed(self, index: int) -> None:
        result = self._current_result()
        if result is None:
            return
        heavy = None
        if result.out_path and result.out_path.is_file():
            from dockflow_core.pdbio import parse_pdbqt

            atoms = parse_pdbqt(result.out_path).atoms
            heavy = sum(1 for a in atoms if a.element.upper() != "H")
        self.results_table.set_poses(result, ligand_heavy_atoms=heavy)
        self.contact_table.set_contacts([])
        self.preview_label.setText("no render yet - press 'Render PNG'")

    def _analyze_current(self) -> None:
        result = self._current_result()
        if result is None or result.out_path is None or not self.state.receptor_pdbqt:
            QMessageBox.information(self, "Nothing to analyse",
                                    "Dock a ligand first.")
            return
        receptor = self.state.receptor_pdbqt
        out_path = result.out_path
        ligand_name = result.ligand_name

        def job():
            from dockflow_core.analyzer import analyze_docking_result

            analyses = analyze_docking_result(
                self._make_result(out_path, ligand_name),
                receptor, top_poses=3,
            )
            return analyses

        self._start_worker(job, on_result=self._on_analysis_done)

    @staticmethod
    def _make_result(out_path, ligand_name):
        from dockflow_core.models import DockingResult
        from dockflow_core.pdbio import parse_pdbqt_results

        return DockingResult(ligand_name=ligand_name,
                             poses=parse_pdbqt_results(out_path),
                             out_path=out_path)

    def _on_analysis_done(self, analyses) -> None:
        if not analyses:
            QMessageBox.information(self, "No poses", "No poses available to analyse.")
            return
        best = analyses[0]
        self.contact_table.set_contacts(best.contacts)
        self._log(f"analysis: pose {best.pose_index} has {best.num_contacts} contacts "
                  f"({best.num_hbonds} hbonds, {best.num_hydrophobic} hydrophobic)",
                  "ok")

    def _render_current(self) -> None:
        result = self._current_result()
        if result is None or result.out_path is None:
            QMessageBox.information(self, "Nothing to render", "Dock a ligand first.")
            return
        receptor = self.state.receptor_pdbqt or self.state.receptor_pdb
        if receptor is None:
            QMessageBox.warning(self, "No receptor", "Prepare a receptor first.")
            return
        out_dir = self.state.workdir / "visualization"
        box = self.state.grid_box
        poses_path = result.out_path
        ligand_name = result.ligand_name
        affinities = [p.affinity for p in result.poses]

        def job():
            from dockflow_core.visualizer import render_best_poses

            return render_best_poses(receptor, poses_path, box, out_dir, ligand_name,
                                     affinities=affinities, engine="auto", top=5)

        self._log(f"rendering poses of {ligand_name}…")
        self._start_worker(job, on_result=self._on_render_done)

    def _on_render_done(self, images) -> None:
        if not images:
            self._log("no images produced", "warn")
            return
        from PyQt6.QtGui import QPixmap

        image = images[0]
        pixmap = QPixmap(str(image))
        if not pixmap.isNull():
            self.preview_label.setPixmap(
                pixmap.scaled(
                    self.preview_label.width(), self.preview_label.height(),
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation,
                )
            )
        self._log(f"rendered {len(images)} image(s) -> {image.parent}", "ok")

    def _open_pymol_current(self) -> None:
        result = self._current_result()
        if result is None or result.out_path is None:
            QMessageBox.information(self, "Nothing to show", "Dock a ligand first.")
            return
        try:
            from dockflow_core.visualizer import (
                open_session,
                render_complex,
            )

            viz_dir = self.state.workdir / "visualization"
            session = viz_dir / f"{result.ligand_name}_session.pse"
            render_complex(
                self.state.receptor_pdbqt, [result.out_path],
                self.state.grid_box, viz_dir / "tmp_preview.png",
                session_path=session, engine="pymol",
            )
            if session.is_file():
                open_session(session, self.config.pymol_exec)
                self._log(f"opened PyMOL session {session}", "ok")
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "PyMOL unavailable", str(exc))

    def _export_csv(self) -> None:
        if not self.state.results:
            QMessageBox.information(self, "Nothing to export", "Run docking first.")
            return
        from dockflow_core.docker_engine import write_summary_csv

        default = str(self.state.workdir / "docking" / "summary.csv")
        path, _ = QFileDialog.getSaveFileName(self, "Export results", default,
                                              "CSV files (*.csv)")
        if not path:
            return
        write_summary_csv(self.state.results, path)
        self._log(f"results exported to {path}", "ok")

    def _open_run_folder(self) -> None:
        folder = self.state.workdir
        if not folder.exists():
            QMessageBox.information(self, "Not found", f"{folder} does not exist yet.")
            return
        import subprocess
        import sys

        if sys.platform.startswith("win"):
            os.startfile(str(folder))  # noqa: S606
        elif sys.platform == "darwin":
            subprocess.Popen(["open", str(folder)])  # noqa: S603
        else:
            subprocess.Popen(["xdg-open", str(folder)])  # noqa: S603

    # ------------------------------------------------------------------ menus etc
    def _new_run(self) -> None:
        answer = QMessageBox.question(
            self, "New run", "Reset the wizard state (files on disk are kept)?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
        if answer != QMessageBox.StandardButton.Yes:
            return
        self.state = GuiState(workdir=self.state.workdir)
        for index in range(6):
            self.stepbar.set_state(index, "pending")
        self.stepbar.set_state(0, "active")
        self.target_info.setText("")
        self.ligand_table.setRowCount(0)
        self.receptor_status.setText("")
        self.ligand_prep_status.setText("")
        self.dock_summary.setRowCount(0)
        self.results_ligand.clear()
        self.results_table.setRowCount(0)
        self.contact_table.setRowCount(0)
        self._goto_step(0)
        self._log("state reset")

    def _open_settings(self) -> None:
        dialog = QDialog(self)
        dialog.setWindowTitle("Settings")
        form = QFormLayout(dialog)
        vina = QComboBox()
        vina.setEditable(True)
        vina.setCurrentText(self.config.vina_exec or "vina")
        smina = QComboBox()
        smina.setEditable(True)
        smina.setCurrentText(self.config.smina_exec or "smina")
        pymol = QComboBox()
        pymol.setEditable(True)
        pymol.setCurrentText(self.config.pymol_exec or "pymol")
        cpu = QComboBox()
        cpu.setEditable(True)
        cpu.setCurrentText(str(self.config.cpu))
        form.addRow("vina executable", vina)
        form.addRow("smina executable", smina)
        form.addRow("pymol executable", pymol)
        form.addRow("cpu threads (0 = all)", cpu)
        from PyQt6.QtWidgets import QDialogButtonBox

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        form.addRow(buttons)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            self.config.vina_exec = vina.currentText() or None
            self.config.smina_exec = smina.currentText() or None
            self.config.pymol_exec = pymol.currentText() or None
            try:
                self.config.cpu = int(cpu.currentText() or 0)
            except ValueError:
                self.config.cpu = 0
            self.config.save()
            self._log("settings saved", "ok")
            self._refresh_engine_labels()

    def _run_pipeline_yaml(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Pipeline YAML", str(self.state.workdir), "YAML files (*.yaml *.yml)")
        if not path:
            return
        try:
            from dockflow_core.pipeline import DockingPipeline, PipelineConfig

            config = PipelineConfig.from_yaml(path)
            config.workdir = self.state.workdir
            self._pipeline = DockingPipeline(config)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Invalid config", str(exc))
            return
        worker = PipelineWorker(self._pipeline)
        worker.signals.step.connect(
            lambda step, status: self._log(f"pipeline {step}: {status}"))
        worker.signals.progress.connect(self._show_progress)
        worker.signals.result.connect(self._on_pipeline_done)
        worker.signals.error.connect(self._on_worker_error)
        self._workers.append(worker)
        worker.start()
        self._log(f"automated pipeline started ({path})", "ok")

    def _on_pipeline_done(self, report) -> None:
        if report.ok:
            self._log(f"pipeline finished -> {report.run_dir}", "ok")
            QMessageBox.information(
                self, "Pipeline finished",
                f"Run {report.run_id} completed.\nReport: {report.run_dir / 'report.md'}")
        else:
            QMessageBox.warning(self, "Pipeline failed", str(report.error))

    def _show_env_report(self) -> None:
        import contextlib
        import io

        from dockflow_core.cli import cmd_info

        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            cmd_info(None)
        QMessageBox.information(self, "Environment report", buffer.getvalue())

    def _show_about(self) -> None:
        from dockflow_core import __version__

        QMessageBox.about(
            self, "About DockFlow-Automator",
            f"<h3 style='color:{ACCENT}'>DockFlow-Automator {__version__}</h3>"
            "<p>Unified, automated molecular docking: download, prepare, "
            "grid box, dock (AutoDock Vina), analyze and visualize.</p>"
            "<p>Built on Meeko, RDKit, OpenBabel, AutoDock Vina and "
            "open-source PyMOL.</p>",
        )



