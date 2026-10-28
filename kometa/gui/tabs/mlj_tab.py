# kometa/gui/tabs/mlj_tab.py
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
    from ..style import *
except ImportError:
    from kometa.gui.style import *

class MLJTab(QWidget):
    def __init__(self):
        super().__init__()
        inner = QWidget(); vbox = QVBoxLayout(inner); vbox.setSpacing(18)
        vbox.addWidget(info_lbl(
            "Marcus-Levich-Jortner parameters — active when Rate model = <b>mlj</b>.<br>"
            "Required for T &lt; 200 K or when Huang-Rhys factor S &gt; 0.3.<br>"
            "MLJ adds Franck-Condon weighted quantum nuclear tunnelling to Marcus theory."))

        box, g = make_group("MLJ Parameters  (from ORCA normal-mode analysis)")
        self.S     = make_dspin(0.5, 0.0, 10.0, 3, 0.05,
                                "Huang-Rhys factor S = λᵢₙ/ħω\nFrom ORCA normal-mode analysis — S = λᵢₙ / ħω")
        self.hw    = make_dspin(0.17, 0.01, 1.0, 3, 0.01,
                                "Dominant quantum mode ħω (eV)\nC=C stretch ≈ 0.17 eV; C-H wag ≈ 0.37 eV")
        self.l_out = make_le("0.15, 0.18",
                             "Outer-sphere λ per layer (eV); λ_outer = λ_total − S·ℏω", 240)
        self.n_max = make_spin(20, 1, 100, "Number of vibrational levels in Franck-Condon sum")

        add_row(g, "Huang-Rhys factor S",
                self.S, "", "S = λ<sub>in</sub>/ℏω — from ORCA normal-mode analysis")
        add_row(g, "ℏω  (dominant mode)",
                self.hw, "eV", "C=C stretch ≈ 0.17 eV")
        add_row(g, "λ<sub>outer</sub>  per layer",
                self.l_out, "eV / layer",
                "λ<sub>outer</sub> = λ<sub>total</sub> − S·ℏω")
        add_row(g, "FC sum n<sub>max</sub>",
                self.n_max, "levels")
        vbox.addWidget(box); vbox.addStretch()
        lay = QVBoxLayout(self); lay.addWidget(scrolled(inner))

    def get_values(self):
        return {
            'huang_rhys_factor': self.S.value(),
            'hbar_omega'       : self.hw.value(),
            'lambda_outer'     : parse_list(self.l_out.text()),
            'mlj_n_max'        : self.n_max.value(),
        }

    def set_values(self, d):
        if 'huang_rhys_factor' in d: self.S.setValue(float(d['huang_rhys_factor']))
        if 'hbar_omega'        in d: self.hw.setValue(float(d['hbar_omega']))
        if 'lambda_outer'      in d: self.l_out.setText(l2s(d['lambda_outer']))
        if 'mlj_n_max'         in d: self.n_max.setValue(int(d['mlj_n_max']))


