# kometa/gui/results_tab.py
import sys, os, json, traceback, subprocess, platform
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget, QDialog,
    QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QCheckBox,
    QFileDialog, QTextEdit, QProgressBar, QSplitter,
    QGroupBox, QScrollArea, QMessageBox, QFrame,
    QStatusBar, QDoubleSpinBox, QSpinBox, QSizePolicy
)
from PySide6.QtCore  import Qt, QThread, Signal
from PySide6.QtGui   import QFont, QColor, QPalette, QPixmap, QAction
import matplotlib
matplotlib.use('QtAgg')

try:
    from .style import *
except ImportError:
    from kometa.gui.style import *

class ResultsTab(QWidget):
    def __init__(self):
        super().__init__()
        lay = QVBoxLayout(self); lay.setContentsMargins(12, 12, 12, 12)
        hdr = QHBoxLayout()
        lbl = QLabel("Output Plots"); lbl.setObjectName("section")
        lbl.setFont(QFont("Palatino Linotype", 15, QFont.Bold))
        self.ref_btn = QPushButton("🔄  Refresh")
        self.ref_btn.setObjectName("browse_btn"); self.ref_btn.setFixedWidth(120)
        hdr.addWidget(lbl); hdr.addStretch(); hdr.addWidget(self.ref_btn)
        lay.addLayout(hdr)
        self.scroll = QScrollArea(); self.scroll.setWidgetResizable(True)
        self.inner  = QWidget(); self.grid = QGridLayout(self.inner)
        self.grid.setSpacing(14)
        self.scroll.setWidget(self.inner)
        lay.addWidget(self.scroll)

    def load_plots(self, folder):
        while self.grid.count():
            item = self.grid.takeAt(0)
            if item.widget(): item.widget().deleteLater()
        if not os.path.isdir(folder):
            l = QLabel(f"Folder '{folder}' not found.\nRun a simulation first.")
            l.setStyleSheet(f"color: {C_GOLD}; font-size: 13pt;")
            self.grid.addWidget(l, 0, 0); return
        pngs = sorted(f for f in os.listdir(folder) if f.endswith('.png'))
        if not pngs:
            l = QLabel("No plots found. Run a simulation first.")
            l.setStyleSheet(f"color: {C_GOLD}; font-size: 13pt;")
            self.grid.addWidget(l, 0, 0); return
        for idx, fname in enumerate(pngs):
            frame = QGroupBox(fname.replace('.png', ''))
            fl = QVBoxLayout(frame)
            img = QLabel()
            pix = QPixmap(os.path.join(folder, fname))
            pix = pix.scaledToWidth(400, Qt.SmoothTransformation)
            img.setPixmap(pix); img.setAlignment(Qt.AlignCenter)
            fl.addWidget(img)
            self.grid.addWidget(frame, idx // 3, idx % 3)


# =============================================================================
# PARAMETER TABS
# =============================================================================

