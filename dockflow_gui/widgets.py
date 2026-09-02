"""Reusable widgets: step bar, ligand table, grid-box editor + 3D preview,
docking parameters, results table and log panel."""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QColor
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from dockflow_core.gridbox import GridBox
from dockflow_core.models import LigandRecord

from .resources import ACCENT, BORDER, DANGER, MUTED, SUCCESS, WARNING

# ---------------------------------------------------------------------------
# Step indicator
# ---------------------------------------------------------------------------
_STEPS = ["Target", "Ligands", "Prepare", "Grid box", "Docking", "Results"]


class StepBar(QWidget):
    """Vertical step indicator with pending / active / done / error states."""

    stepSelected = pyqtSignal(int)
    stepActivated = pyqtSignal(int)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._states = ["pending"] * len(_STEPS)
        self.setMinimumWidth(190)
        self.setMinimumHeight(320)

    def set_state(self, index: int, state: str) -> None:
        if 0 <= index < len(self._states):
            self._states[index] = state
            self.update()

    def states(self) -> list[str]:
        return list(self._states)

    def mousePressEvent(self, event) -> None:
        position = event.position().toPoint()
        row_height = self.height() / len(_STEPS)
        index = int(position.y() // row_height)
        if 0 <= index < len(_STEPS):
            self.stepSelected.emit(index)
            self.stepActivated.emit(index)

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt naming
        from PyQt6.QtGui import QFont, QPainter, QPen

        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        row_height = self.height() / len(_STEPS)
        center_x = 28
        title_font = QFont(self.font())
        title_font.setWeight(QFont.Weight.DemiBold)
        title_font.setPointSizeF(self.font().pointSizeF() * 1.15)
        for index, name in enumerate(_STEPS):
            y = row_height * (index + 0.5)
            state = self._states[index]
            color = {
                "pending": QColor(BORDER),
                "active": QColor(ACCENT),
                "done": QColor(SUCCESS),
                "error": QColor(DANGER),
                "skipped": QColor(MUTED),
            }.get(state, QColor(BORDER))
            # connecting line
            if index < len(_STEPS) - 1:
                pen = QPen(QColor(BORDER))
                pen.setWidth(2)
                painter.setPen(pen)
                painter.drawLine(int(center_x), int(y + 12), int(center_x),
                                 int(row_height * (index + 1.5) - 12))
            # circle
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(color)
            painter.drawEllipse(int(center_x - 11), int(y - 11), 22, 22)
            # label inside circle
            painter.setPen(QColor("white"))
            label_font = QFont(self.font())
            label_font.setBold(True)
            painter.setFont(label_font)
            painter.drawText(int(center_x - 11), int(y - 11), 22, 22,
                             Qt.AlignmentFlag.AlignCenter,
                             "✓" if state == "done" else str(index + 1))
            # title
            painter.setPen(QColor("#1f2933" if state != "pending" else MUTED))
            painter.setFont(title_font)
            painter.drawText(int(center_x + 20), int(y - 12), 140, 24,
                             Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                             name)
        painter.end()


# ---------------------------------------------------------------------------
# Ligand table
# ---------------------------------------------------------------------------
class LigandTableWidget(QTableWidget):
    """Editable ligand list with add / remove helpers."""

    addSmilesRequested = pyqtSignal()
    addFilesRequested = pyqtSignal()
    addPubchemRequested = pyqtSignal()
    addZincRequested = pyqtSignal()
    addCocystalRequested = pyqtSignal()
    removeRequested = pyqtSignal()

    def __init__(self, parent=None) -> None:
        super().__init__(0, 4, parent)
        self.setHorizontalHeaderLabels(["identifier", "source", "value", "status"])
        self.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setMinimumHeight(220)

    # -- layout with buttons ------------------------------------------------
    @staticmethod
    def with_buttons(table: LigandTableWidget) -> QWidget:
        wrap = QWidget()
        layout = QVBoxLayout(wrap)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(table)
        buttons = QHBoxLayout()
        add_smiles = QPushButton("Add SMILES")
        add_files = QPushButton("Add files…")
        add_pubchem = QPushButton("PubChem…")
        add_zinc = QPushButton("ZINC…")
        add_cocystal = QPushButton("Co-crystal ligand…")
        remove = QPushButton("Remove")
        remove.setProperty("secondary", True)
        for button, signal in (
            (add_smiles, table.addSmilesRequested),
            (add_files, table.addFilesRequested),
            (add_pubchem, table.addPubchemRequested),
            (add_zinc, table.addZincRequested),
            (add_cocystal, table.addCocystalRequested),
            (remove, table.removeRequested),
        ):
            buttons.addWidget(button)
            button.clicked.connect(signal.emit)
        layout.addLayout(buttons)
        return wrap

    # -- content -------------------------------------------------------------
    def set_ligands(self, ligands: Sequence[LigandRecord]) -> None:
        self.setRowCount(len(ligands))
        for row, ligand in enumerate(ligands):
            value = ligand.value if ligand.source != "file" else (ligand.path or "")
            status = ligand.status
            for column, text in enumerate(
                (ligand.identifier, ligand.source, str(value), status)
            ):
                item = QTableWidgetItem(text)
                if status == "error":
                    item.setForeground(QColor(DANGER))
                elif status in ("prepared", "docked"):
                    item.setForeground(QColor(SUCCESS))
                elif status in ("downloaded", "pending"):
                    item.setForeground(QColor(WARNING))
                self.setItem(row, column, item)
        self.resizeColumnsToContents()
        self.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)

    def selected_rows(self) -> list[int]:
        return sorted({index.row() for index in self.selectedIndexes()})


# ---------------------------------------------------------------------------
# Grid box editor + interactive preview
# ---------------------------------------------------------------------------
class GridBoxWidget(QWidget):
    """Center/size spin boxes + padding + auto-compute buttons."""

    boxChanged = pyqtSignal(object)  # GridBox
    computeFromLigandRequested = pyqtSignal()
    computeFromResiduesRequested = pyqtSignal()

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        group = QGroupBox("Search space (Vina grid box)")
        form = QFormLayout(group)
        self.spins: dict[str, QDoubleSpinBox] = {}
        for axis in "xyz":
            center = _spin(0.0, -1000.0, 1000.0, 0.1)
            size = _spin(22.0, 5.0, 120.0, 0.5)
            self.spins[f"center_{axis}"] = center
            self.spins[f"size_{axis}"] = size
            form.addRow(f"center {axis.upper()} (A)", center)
            form.addRow(f"size {axis.upper()} (A)", size)
        self.padding = _spin(4.0, 0.0, 20.0, 0.5)
        form.addRow("padding (A)", self.padding)
        buttons = QHBoxLayout()
        from_ligand = QPushButton("From co-crystal ligand")
        from_residues = QPushButton("From residues…")
        apply = QPushButton("Apply box")
        apply.setProperty("secondary", True)
        buttons.addWidget(from_ligand)
        buttons.addWidget(from_residues)
        buttons.addWidget(apply)
        form.addRow(buttons)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(group)
        from_ligand.clicked.connect(self.computeFromLigandRequested.emit)
        from_residues.clicked.connect(self.computeFromResiduesRequested.emit)
        apply.clicked.connect(self._emit_box)
        for spin in self.spins.values():
            spin.valueChanged.connect(self._emit_box)

    def _emit_box(self) -> None:
        self.boxChanged.emit(self.box())

    def box(self) -> GridBox:
        return GridBox(
            center=tuple(self.spins[f"center_{a}"].value() for a in "xyz"),  # type: ignore[arg-type]
            size=tuple(self.spins[f"size_{a}"].value() for a in "xyz"),  # type: ignore[arg-type]
            source="gui",
            padding=self.padding.value(),
        )

    def set_box(self, box: GridBox) -> None:
        for index, axis in enumerate("xyz"):
            self.spins[f"center_{axis}"].setValue(float(box.center[index]))
            self.spins[f"size_{axis}"].setValue(float(box.size[index]))
        self._emit_box()


class GridBoxPreview(QWidget):
    """Lightweight interactive 3D preview of atoms + grid box.

    Rotation: drag with the left mouse button.  Zoom: mouse wheel.
    Pure QPainter - no OpenGL dependency.
    """

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setMinimumSize(320, 300)
        self._points: np.ndarray = np.zeros((0, 3))
        self._colors: list[QColor] = []
        self._box: GridBox | None = None
        self._yaw = 0.6
        self._pitch = 0.45
        self._zoom = 1.0
        self._hint = "drag to rotate - wheel to zoom"

    # -- public --------------------------------------------------------------
    def set_data(self, points: np.ndarray, colors: Sequence[QColor],
                 box: GridBox | None) -> None:
        self._points = np.asarray(points, dtype=float).reshape(-1, 3)
        self._colors = list(colors)[: len(self._points)]
        self._box = box
        self.update()

    # -- interaction ---------------------------------------------------------
    def mousePressEvent(self, event) -> None:
        self._last = event.position()

    def mouseMoveEvent(self, event) -> None:
        position = event.position()
        dx = position.x() - self._last.x()
        dy = position.y() - self._last.y()
        self._last = position
        self._yaw += dx * 0.01
        self._pitch = max(-1.4, min(1.4, self._pitch + dy * 0.01))
        self.update()

    def wheelEvent(self, event) -> None:  # noqa: N802 - Qt naming
        delta = event.angleDelta().y()
        self._zoom = max(0.25, min(8.0, self._zoom * (1.1 if delta > 0 else 0.9)))
        self.update()

    # -- rendering -----------------------------------------------------------
    def _rotation(self) -> np.ndarray:
        cy, sy = math.cos(self._yaw), math.sin(self._yaw)
        cp, sp = math.cos(self._pitch), math.sin(self._pitch)
        yaw = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
        pitch = np.array([[1, 0, 0], [0, cp, -sp], [0, sp, cp]])
        return pitch @ yaw

    def paintEvent(self, event) -> None:  # noqa: N802 - Qt naming
        from PyQt6.QtGui import QFont, QPainter, QPen

        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#ffffff"))
        width, height = self.width(), self.height()
        rotation = self._rotation()
        # fit scale: use box or data extent
        reference = self._box.corner_points() if self._box else self._points
        if reference.size == 0:
            painter.setPen(QColor(MUTED))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "no data")
            painter.end()
            return
        projected = reference @ rotation.T
        span = max(projected[:, 0].ptp() or 1.0, projected[:, 1].ptp() or 1.0)
        scale = min(width, height) * 0.38 * self._zoom / (span / 2)
        center = (width / 2, height / 2)
        painter.setPen(QPen(QColor(MUTED)))
        hint_font = QFont(self.font())
        hint_font.setPointSizeF(8)
        painter.setFont(hint_font)
        painter.drawText(8, height - 10, self._hint)
        # grid box wireframe
        if self._box is not None:
            corners = self._box.corner_points() @ rotation.T
            px = corners[:, 0] * scale + center[0]
            py = -corners[:, 1] * scale + center[1]
            edges = [
                (0, 1), (1, 2), (2, 3), (3, 0),
                (4, 5), (5, 6), (6, 7), (7, 4),
                (0, 4), (1, 5), (2, 6), (3, 7),
            ]
            pen = QPen(QColor(ACCENT))
            pen.setWidthF(1.4)
            pen.setStyle(Qt.PenStyle.DashLine)
            painter.setPen(pen)
            for a, b in edges:
                painter.drawLine(
                    int(px[a]), int(py[a]), int(px[b]), int(py[b])
                )
        # atoms
        if self._points.size:
            proj = self._points @ rotation.T
            px = proj[:, 0] * scale + center[0]
            py = -proj[:, 1] * scale + center[1]
            depth = proj[:, 2]
            z_lo, z_hi = float(depth.min()), float(depth.max()) or 1.0
            order = np.argsort(depth)
            painter.setPen(Qt.PenStyle.NoPen)
            for index in order:
                color = self._colors[int(index)] if index < len(self._colors) \
                    else QColor(MUTED)
                t = (depth[index] - z_lo) / (z_hi - z_lo) if z_hi > z_lo else 0.5
                radius = 1.6 + 2.6 * (1.0 - t)
                painter.setBrush(color)
                painter.drawEllipse(
                    int(px[index] - radius), int(py[index] - radius),
                    int(2 * radius), int(2 * radius),
                )
        painter.end()


# ---------------------------------------------------------------------------
# Docking parameters
# ---------------------------------------------------------------------------
class DockParamsWidget(QWidget):
    """AutoDock Vina parameter editors with sensible bounds."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        group = QGroupBox("Vina parameters")
        form = QFormLayout(group)
        self.exhaustiveness = QSpinBox()
        self.exhaustiveness.setRange(1, 256)
        self.exhaustiveness.setValue(8)
        self.num_modes = QSpinBox()
        self.num_modes.setRange(1, 100)
        self.num_modes.setValue(9)
        self.refine = QSpinBox()
        self.refine.setRange(0, 50)
        self.refine.setValue(5)
        self.cpu = QSpinBox()
        self.cpu.setRange(0, 256)
        self.cpu.setValue(0)
        self.cpu.setToolTip("0 = use all cores")
        self.seed = QLineEdit("random")
        self.seed.setToolTip("integer seed or 'random'")
        self.scoring = QComboBox()
        self.scoring.addItems(["vina", "vinardo", "ad4"])
        self.energy_range = _spin(3.0, 0.5, 20.0, 0.5)
        form.addRow("exhaustiveness", self.exhaustiveness)
        form.addRow("num_modes", self.num_modes)
        form.addRow("refine (inner iters)", self.refine)
        form.addRow("cpu threads", self.cpu)
        form.addRow("seed", self.seed)
        form.addRow("scoring function", self.scoring)
        form.addRow("energy_range (kcal/mol)", self.energy_range)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(group)

    def config_dict(self) -> dict:
        seed_text = self.seed.text().strip()
        seed = None if seed_text.lower() in ("", "random") else int(seed_text)
        return {
            "exhaustiveness": self.exhaustiveness.value(),
            "num_modes": self.num_modes.value(),
            "refine": self.refine.value(),
            "cpu": self.cpu.value(),
            "seed": seed,
            "scoring": self.scoring.currentText(),
            "energy_range": self.energy_range.value(),
        }


#
# ---------------------------------------------------------------------------
# Results table
# ---------------------------------------------------------------------------
class ResultsTableWidget(QTableWidget):
    """Pose table for one docking result (model, affinity, RMSDs)."""

    poseSelected = pyqtSignal(int)

    def __init__(self, parent=None) -> None:
        super().__init__(0, 5, parent)
        self.setHorizontalHeaderLabels(
            ["pose", "affinity (kcal/mol)", "rmsd l.b.", "rmsd u.b.", "LE (kcal/heavy)"]
        )
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.verticalHeader().setVisible(False)

    def set_poses(self, result, ligand_heavy_atoms: int | None = None) -> None:
        from PyQt6.QtGui import QFont

        self.setRowCount(len(result.poses))
        font = QFont(self.font())
        for row, pose in enumerate(result.poses):
            le = ""
            if ligand_heavy_atoms:
                le = f"{pose.affinity / ligand_heavy_atoms:.2f}"
            values = (
                str(pose.model),
                f"{pose.affinity:.2f}",
                f"{pose.rmsd_lb:.2f}",
                f"{pose.rmsd_ub:.2f}",
                le,
            )
            for column, text in enumerate(values):
                item = QTableWidgetItem(text)
                if column == 1:
                    item.setFont(font)
                    color = (
                        QColor(SUCCESS) if pose.affinity <= -7 else
                        QColor(WARNING) if pose.affinity <= -5 else QColor(MUTED)
                    )
                    item.setForeground(color)
                self.setItem(row, column, item)
        self.resizeColumnsToContents()
        self.horizontalHeader().setStretchLastSection(True)

    def selected_pose(self) -> int | None:
        rows = self.selectedIndexes()
        return rows[0].row() if rows else None


# ---------------------------------------------------------------------------
# Contact table
# ---------------------------------------------------------------------------
class ContactTableWidget(QTableWidget):
    """Interaction table (kind, residue, atoms, distance)."""

    def __init__(self, parent=None) -> None:
        super().__init__(0, 5, parent)
        self.setHorizontalHeaderLabels(
            ["kind", "receptor residue", "receptor atom", "ligand atom", "distance (A)"]
        )
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.verticalHeader().setVisible(False)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)

    def set_contacts(self, contacts: Sequence) -> None:
        self.setRowCount(len(contacts))
        kinds = {
            "hbond": (SUCCESS, "hydrogen bond"),
            "hydrophobic": (MUTED, "hydrophobic"),
            "ionic": (WARNING, "ionic"),
            "metal": (DANGER, "metal coordination"),
            "close": (MUTED, "close contact"),
        }
        for row, contact in enumerate(contacts):
            color, label = kinds.get(contact.kind, (MUTED, contact.kind))
            residue = f"{contact.receptor_resname}{contact.receptor_resseq}" \
                      f" {contact.receptor_chain}".strip()
            values = (
                label,
                residue,
                contact.receptor_atom_name,
                contact.ligand_atom_name,
                f"{contact.distance:.2f}",
            )
            for column, text in enumerate(values):
                item = QTableWidgetItem(text)
                if column == 0:
                    item.setForeground(QColor(color))
                self.setItem(row, column, item)
        self.resizeColumnsToContents()
        self.horizontalHeader().setStretchLastSection(True)


# ---------------------------------------------------------------------------
# Log panel
# ---------------------------------------------------------------------------
class LogPanel(QPlainTextEdit):
    """Append-only log panel with level prefixes."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setReadOnly(True)
        self.setMaximumBlockCount(4000)

    def append_log(self, message: str, level: str = "info") -> None:
        import datetime

        stamp = datetime.datetime.now().strftime("%H:%M:%S")
        prefix = {"info": "·", "ok": "✓", "warn": "!", "error": "✗"}.get(level, "·")
        self.appendPlainText(f"{stamp} {prefix} {message}")


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _spin(value: float, minimum: float, maximum: float, step: float) -> QDoubleSpinBox:
    spin = QDoubleSpinBox()
    spin.setRange(minimum, maximum)
    spin.setValue(value)
    spin.setSingleStep(step)
    spin.setDecimals(2)
    return spin


def make_caption(text: str) -> QLabel:
    label = QLabel(text)
    label.setProperty("muted", True)
    label.setWordWrap(True)
    return label

