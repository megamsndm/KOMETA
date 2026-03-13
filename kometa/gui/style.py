# kometa/gui/style.py
# Design tokens and Qt stylesheet

"""
=============================================================================
KOMETA — KMC Organic Molecular Exciton Transport Analysis
=============================================================================
Version 1.0
Developer: Dr. Mohsen Doust Mohammadi

Run:    python kometa_gui.py
Build:  pyinstaller --onefile --windowed --name KOMETA kometa_gui.py
=============================================================================
"""

import sys, os, json, traceback, importlib.util, subprocess, platform

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget, QDialog,
    QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QCheckBox,
    QFileDialog, QTextEdit, QProgressBar, QSplitter,
    QGroupBox, QScrollArea, QMessageBox, QFrame, QAction,
    QStatusBar, QDoubleSpinBox, QSpinBox, QSizePolicy
)
from PyQt5.QtCore  import Qt, QThread, pyqtSignal
from PyQt5.QtGui   import QFont, QColor, QPalette, QPixmap, QFontDatabase

import matplotlib
matplotlib.use('Qt5Agg')

# =============================================================================
# PALETTE
# =============================================================================
C_BG      = "#07111F"
C_PANEL   = "#0D1E33"
C_CARD    = "#112440"
C_BORDER  = "#1E3A5A"
C_TEAL    = "#00D4EE"
C_TEAL2   = "#0099BB"
C_GOLD    = "#F5A623"
C_GOLD2   = "#C07800"
C_GREEN   = "#00E87A"
C_GREEN2  = "#00A050"
C_RED     = "#FF4C4C"
C_WHITE   = "#ECF3FF"
C_SILVER  = "#7A95B8"
C_MUTED   = "#2E4A6A"

# =============================================================================
# STYLESHEET
# =============================================================================
QSS = f"""
* {{
    font-family: 'Palatino Linotype', 'Book Antiqua', Palatino, serif;
    font-size: 14pt;
    color: {C_WHITE};
}}
QMainWindow, QDialog, QWidget {{
    background: {C_BG};
}}
/* ─── Menu bar ─────────────────────────────── */
QMenuBar {{
    background: {C_PANEL};
    border-bottom: 2px solid {C_TEAL};
    padding: 3px 8px;
    font-size: 14pt;
}}
QMenuBar::item {{
    padding: 6px 16px;
    border-radius: 4px;
    background: transparent;
}}
QMenuBar::item:selected {{
    background: {C_BORDER};
    color: {C_TEAL};
}}
QMenu {{
    background: {C_PANEL};
    border: 1px solid {C_TEAL};
    padding: 4px;
    font-size: 13pt;
}}
QMenu::item {{
    padding: 8px 32px 8px 18px;
    border-radius: 3px;
}}
QMenu::item:selected {{
    background: {C_TEAL};
    color: {C_BG};
}}
QMenu::separator {{
    height: 1px;
    background: {C_MUTED};
    margin: 4px 0;
}}
/* ─── Tabs ─────────────────────────────────── */
QTabWidget::pane {{
    border: 1px solid {C_BORDER};
    background: {C_PANEL};
    border-radius: 6px;
    top: -1px;
}}
QTabBar::tab {{
    background: {C_CARD};
    color: {C_SILVER};
    padding: 10px 20px;
    border-top-left-radius: 5px;
    border-top-right-radius: 5px;
    font-size: 13pt;
    min-width: 110px;
    margin-right: 2px;
}}
QTabBar::tab:selected {{
    background: {C_TEAL};
    color: {C_BG};
    font-weight: bold;
}}
QTabBar::tab:hover:!selected {{
    background: {C_MUTED};
    color: {C_WHITE};
}}
/* ─── Group boxes ───────────────────────────── */
QGroupBox {{
    border: 1px solid {C_TEAL2};
    border-radius: 6px;
    margin-top: 16px;
    padding-top: 12px;
    font-size: 13pt;
    font-weight: bold;
    color: {C_TEAL};
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    left: 14px;
    padding: 0 6px;
}}
/* ─── Inputs ─────────────────────────────────── */
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{
    background: {C_CARD};
    border: 1px solid {C_BORDER};
    border-radius: 4px;
    padding: 6px 10px;
    color: {C_WHITE};
    font-size: 13pt;
    min-height: 30px;
    selection-background-color: {C_TEAL};
    selection-color: {C_BG};
}}
QLineEdit:focus, QComboBox:focus,
QSpinBox:focus, QDoubleSpinBox:focus {{
    border: 1.5px solid {C_TEAL};
}}
QComboBox::drop-down {{ border: none; width: 26px; }}
QComboBox QAbstractItemView {{
    background: {C_CARD};
    color: {C_WHITE};
    selection-background-color: {C_TEAL};
    selection-color: {C_BG};
    font-size: 13pt;
    border: 1px solid {C_TEAL};
}}
QSpinBox::up-button, QSpinBox::down-button,
QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {{
    background: {C_MUTED};
    border: none;
    width: 20px;
    border-radius: 2px;
}}
/* ─── Buttons ────────────────────────────────── */
QPushButton {{
    background: {C_TEAL};
    color: {C_BG};
    border: none;
    border-radius: 5px;
    padding: 8px 22px;
    font-weight: bold;
    font-size: 13pt;
    min-height: 32px;
}}
QPushButton:hover  {{ background: #22E8FF; }}
QPushButton:pressed {{ background: {C_TEAL2}; }}
QPushButton:disabled {{ background: {C_MUTED}; color: {C_SILVER}; }}
QPushButton#run_btn {{
    background: qlineargradient(x1:0,y1:0,x2:1,y2:0,
        stop:0 {C_GREEN2}, stop:1 {C_TEAL2});
    color: {C_BG};
    font-size: 14pt;
    font-weight: bold;
    padding: 10px 28px;
    min-height: 40px;
    max-height: 44px;
    border-radius: 6px;
    letter-spacing: 1px;
}}
QPushButton#run_btn:hover {{
    background: qlineargradient(x1:0,y1:0,x2:1,y2:0,
        stop:0 {C_GREEN}, stop:1 #00CCDD);
}}
QPushButton#stop_btn {{
    background: {C_RED};
    color: white;
    font-size: 13pt;
    font-weight: bold;
    padding: 8px 20px;
    min-height: 40px;
    max-height: 44px;
    border-radius: 6px;
}}
QPushButton#stop_btn:hover {{ background: #FF7777; }}
QPushButton#browse_btn {{
    background: {C_CARD};
    color: {C_TEAL};
    border: 1px solid {C_TEAL2};
    font-size: 11pt;
    padding: 4px 10px;
    min-height: 24px;
    max-height: 30px;
}}
QPushButton#browse_btn:hover {{ background: {C_MUTED}; }}
/* ─── Console ────────────────────────────────── */
QTextEdit {{
    background: #040C16;
    color: #66FFB2;
    font-family: 'Consolas', 'Courier New', monospace;
    font-size: 12pt;
    border: 1px solid {C_BORDER};
    border-radius: 4px;
    line-height: 1.5;
}}
/* ─── Progress bar ───────────────────────────── */
QProgressBar {{
    background: {C_CARD};
    border: 1px solid {C_BORDER};
    border-radius: 5px;
    text-align: center;
    color: {C_WHITE};
    font-size: 12pt;
    min-height: 24px;
}}
QProgressBar::chunk {{
    background: qlineargradient(x1:0,y1:0,x2:1,y2:0,
        stop:0 {C_TEAL2}, stop:1 {C_TEAL});
    border-radius: 4px;
}}
/* ─── Checkboxes ─────────────────────────────── */
QCheckBox {{
    spacing: 10px;
    font-size: 13pt;
    color: {C_WHITE};
}}
QCheckBox::indicator {{
    width: 20px; height: 20px;
    border: 1.5px solid {C_BORDER};
    border-radius: 4px;
    background: {C_CARD};
}}
QCheckBox::indicator:checked {{
    background: {C_TEAL};
    border: 1.5px solid {C_TEAL};
}}
/* ─── Scroll bars ────────────────────────────── */
QScrollArea {{ border: none; background: transparent; }}
QScrollBar:vertical {{
    background: {C_PANEL};
    width: 10px;
    border-radius: 5px;
    margin: 0;
}}
QScrollBar::handle:vertical {{
    background: {C_MUTED};
    border-radius: 5px;
    min-height: 28px;
}}
QScrollBar::add-line:vertical,
QScrollBar::sub-line:vertical {{ height: 0px; }}
/* ─── Status bar ─────────────────────────────── */
QStatusBar {{
    background: {C_PANEL};
    color: {C_SILVER};
    border-top: 1px solid {C_BORDER};
    font-size: 11pt;
    padding: 2px 10px;
}}
/* ─── Splitter ───────────────────────────────── */
QSplitter::handle:horizontal {{
    background: {C_BORDER};
    width: 2px;
}}
/* ─── Labels ─────────────────────────────────── */
QLabel#unit {{
    color: {C_GOLD};
    font-size: 12pt;
    font-style: italic;
    min-width: 70px;
}}
QLabel#info {{
    color: {C_GOLD};
    font-size: 12pt;
    font-style: italic;
    padding: 8px 14px;
    background: #141F0A;
    border-left: 3px solid {C_GOLD};
    border-radius: 3px;
}}
QLabel#section {{
    color: {C_TEAL};
    font-size: 16pt;
    font-weight: bold;
}}
"""

# =============================================================================
# RICH TEXT LABEL  — renders HTML superscripts/subscripts
# =============================================================================
def rich(text):
    """Return a QLabel that renders HTML (for sub/superscripts)."""
    l = QLabel(); l.setText(text); l.setTextFormat(Qt.RichText)
    l.setFont(QFont("Palatino Linotype", 13))
    return l

def unit_lbl(text):
    l = rich(text); l.setObjectName("unit"); return l

def info_lbl(text):
    l = QLabel(); l.setText(text); l.setTextFormat(Qt.RichText)
    l.setObjectName("info"); l.setWordWrap(True); return l

def scrolled(w):
    sc = QScrollArea(); sc.setWidgetResizable(True); sc.setWidget(w); return sc

def make_le(default="", tip="", w=200):
    e = QLineEdit(str(default)); e.setToolTip(tip); e.setMinimumWidth(w); return e

def make_spin(val=0, lo=0, hi=10**8, tip="", w=180):
    s = QSpinBox(); s.setRange(lo, hi); s.setValue(int(val))
    s.setToolTip(tip); s.setMinimumWidth(w); return s

def make_dspin(val=0.0, lo=-1e9, hi=1e9, dec=4, step=0.001, tip="", w=180):
    s = QDoubleSpinBox(); s.setRange(lo, hi); s.setDecimals(dec)
    s.setSingleStep(step); s.setValue(float(val))
    s.setToolTip(tip); s.setMinimumWidth(w); return s

def make_group(title):
    box = QGroupBox(title)
    g = QGridLayout()
    g.setColumnStretch(3, 1)
    g.setHorizontalSpacing(16)
    g.setVerticalSpacing(10)
    box.setLayout(g)
    return box, g

def add_row(g, label_html, widget, unit_html="", tip=""):
    """Add a row using rich HTML label + widget + unit."""
    r = g.rowCount()
    lbl = rich(label_html); lbl.setToolTip(tip)
    g.addWidget(lbl,    r, 0)
    g.addWidget(widget, r, 1)
    if unit_html:
        g.addWidget(unit_lbl(unit_html), r, 2)

def parse_list(text):
    try:    return [float(x.strip()) for x in text.split(',') if x.strip()]
    except: return [0.0]

def l2s(v):
    if isinstance(v, (list, tuple)): return ', '.join(str(x) for x in v)
    return str(v)

# =============================================================================
# HEADER BAR  — KOMETA title + molecule name
# =============================================================================
