# kometa/gui/header_bar.py
import sys, os, json, traceback, subprocess, platform
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget, QDialog,
    QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QCheckBox,
    QFileDialog, QTextEdit, QProgressBar, QSplitter,
    QGroupBox, QScrollArea, QMessageBox, QFrame, QAction,
    QStatusBar, QDoubleSpinBox, QSpinBox, QSizePolicy
)
from PyQt5.QtCore  import Qt, QThread, pyqtSignal
from PyQt5.QtGui   import QFont, QColor, QPalette, QPixmap
import matplotlib
matplotlib.use('Qt5Agg')

try:
    from .style import *
except ImportError:
    from kometa.gui.style import *

class HeaderBar(QWidget):
    def __init__(self):
        super().__init__()
        self.setFixedHeight(68)
        self.setStyleSheet(f"""
            background: qlineargradient(x1:0,y1:0,x2:1,y2:0,
                stop:0 #071828, stop:0.5 #0D2240, stop:1 #071828);
            border-bottom: 2px solid {C_TEAL};
        """)
        lay = QHBoxLayout(self); lay.setContentsMargins(20, 6, 20, 6)

        # App name
        name_lbl = QLabel("KOMETA")
        name_lbl.setFont(QFont("Palatino Linotype", 26, QFont.Bold))
        name_lbl.setStyleSheet(f"color: {C_TEAL}; letter-spacing: 6px;")

        # Molecule input
        mol_lbl = QLabel("Molecule:")
        mol_lbl.setFont(QFont("Palatino Linotype", 13, QFont.Bold))
        mol_lbl.setStyleSheet(f"color: {C_SILVER};")

        self.mol_name = QLineEdit("PtOEP")
        self.mol_name.setFont(QFont("Palatino Linotype", 15, QFont.Bold))
        self.mol_name.setFixedWidth(200)
        self.mol_name.setFixedHeight(38)
        self.mol_name.setStyleSheet(f"""
            background: {C_CARD};
            color: {C_GREEN};
            border: 2px solid {C_GREEN2};
            border-radius: 5px;
            padding: 4px 14px;
        """)
        self.mol_name.setToolTip(
            "Molecule name — works for ANY organic molecule.\n"
            "PtOEP · Rubrene · DPA · TIPS-pentacene · etc.")

        # Status
        self.status_dot = QLabel("●")
        self.status_dot.setFont(QFont("Palatino Linotype", 20))
        self.status_dot.setStyleSheet(f"color: {C_SILVER};")
        self.status_txt = QLabel("Ready")
        self.status_txt.setFont(QFont("Palatino Linotype", 12))
        self.status_txt.setStyleSheet(f"color: {C_SILVER};")

        lay.addWidget(name_lbl)
        lay.addStretch()
        lay.addWidget(mol_lbl)
        lay.addSpacing(8)
        lay.addWidget(self.mol_name)
        lay.addSpacing(24)
        lay.addWidget(self.status_dot)
        lay.addWidget(self.status_txt)

    def set_running(self):
        self.status_dot.setStyleSheet(f"color: {C_GOLD};")
        self.status_txt.setText("Running…")
        self.status_txt.setStyleSheet(f"color: {C_GOLD};")

    def set_done(self):
        self.status_dot.setStyleSheet(f"color: {C_GREEN};")
        self.status_txt.setText("Done")
        self.status_txt.setStyleSheet(f"color: {C_GREEN};")

    def set_error(self):
        self.status_dot.setStyleSheet(f"color: {C_RED};")
        self.status_txt.setText("Error")
        self.status_txt.setStyleSheet(f"color: {C_RED};")

    def set_ready(self):
        self.status_dot.setStyleSheet(f"color: {C_SILVER};")
        self.status_txt.setText("Ready")
        self.status_txt.setStyleSheet(f"color: {C_SILVER};")


# =============================================================================
# SIM WORKER
# =============================================================================
