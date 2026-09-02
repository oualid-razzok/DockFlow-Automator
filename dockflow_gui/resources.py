"""Stylesheet, colors and icon drawing (no binary assets)."""

from __future__ import annotations

ACCENT = "#2563eb"
ACCENT_HOVER = "#1d4ed8"
DANGER = "#dc2626"
SUCCESS = "#16a34a"
WARNING = "#d97706"
BG = "#f5f6f8"
CARD = "#ffffff"
BORDER = "#d4d7dc"
TEXT = "#1f2933"
MUTED = "#6b7280"

APP_STYLESHEET = f"""
QMainWindow, QDialog {{
    background: {BG};
}}
QGroupBox {{
    background: {CARD};
    border: 1px solid {BORDER};
    border-radius: 8px;
    margin-top: 12px;
    padding: 12px 8px 8px 8px;
    font-weight: 600;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    left: 12px;
    padding: 0 4px;
    color: {TEXT};
}}
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox, QTextEdit, QPlainTextEdit {{
    background: {CARD};
    border: 1px solid {BORDER};
    border-radius: 6px;
    padding: 5px 8px;
    color: {TEXT};
}}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
    border: 2px solid {ACCENT};
    padding: 4px 7px;
}}
QPushButton {{
    background: {ACCENT};
    color: white;
    border: none;
    border-radius: 6px;
    padding: 7px 16px;
    font-weight: 600;
}}
QPushButton:hover {{ background: {ACCENT_HOVER}; }}
QPushButton:disabled {{ background: #b6c2d6; }}
QPushButton[danger="true"] {{ background: {DANGER}; }}
QPushButton[secondary="true"] {{
    background: {CARD};
    color: {TEXT};
    border: 1px solid {BORDER};
}}
QPushButton[secondary="true"]:hover {{ border-color: {ACCENT}; color: {ACCENT}; }}
QTableWidget {{
    background: {CARD};
    border: 1px solid {BORDER};
    border-radius: 6px;
    gridline-color: {BORDER};
    selection-background-color: #dbeafe;
    selection-color: {TEXT};
}}
QHeaderView::section {{
    background: #e9edf2;
    border: none;
    border-bottom: 1px solid {BORDER};
    padding: 6px;
    font-weight: 600;
}}
QProgressBar {{
    background: {CARD};
    border: 1px solid {BORDER};
    border-radius: 6px;
    text-align: center;
    height: 18px;
}}
QProgressBar::chunk {{
    background: {ACCENT};
    border-radius: 5px;
}}
QDockWidget {{
    color: {TEXT};
    titlebar-close-icon: none;
    titlebar-normal-icon: none;
}}
QPlainTextEdit, QTextEdit {{
    font-family: 'DejaVu Sans Mono', 'Consolas', monospace;
    font-size: 12px;
}}
QLabel[muted="true"] {{ color: {MUTED}; }}
QStatusBar {{
    background: {CARD};
    border-top: 1px solid {BORDER};
}}
QTabWidget::pane {{ border: 1px solid {BORDER}; border-radius: 6px; }}
"""


def app_icon():
    """Draw a simple molecular icon (two bonded circles) as a QIcon."""
    from PyQt6.QtGui import QBrush, QColor, QIcon, QPainter, QPen, QPixmap

    pixmap = QPixmap(64, 64)
    pixmap.fill(QColor(0, 0, 0, 0))
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    pen = QPen(QColor(ACCENT))
    pen.setWidth(3)
    painter.setPen(pen)
    painter.setBrush(QBrush(QColor(ACCENT)))
    painter.drawEllipse(8, 18, 24, 24)
    painter.setBrush(QBrush(QColor("#f59e0b")))
    painter.drawEllipse(34, 24, 20, 20)
    painter.drawLine(32, 30, 34, 34)
    painter.end()
    return QIcon(pixmap)
