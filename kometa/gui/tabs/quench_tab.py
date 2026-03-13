# kometa/gui/tabs/quench_tab.py
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
    from ..style import *
except ImportError:
    from kometa.gui.style import *

class QuenchTab(QWidget):
    def __init__(self):
        super().__init__()
        inner = QWidget(); vbox = QVBoxLayout(inner); vbox.setSpacing(18)
        vbox.addWidget(info_lbl(
            "Impurity quenching — Module 4B (DOI: 10.1103/PhysRevApplied.14.034050)<br><br>"
            "k<sub>quench</sub> is physically <b>separate</b> from k<sub>decay</sub>:<br>"
            "• k<sub>decay</sub> = natural phosphorescence on host molecule<br>"
            "• k<sub>quench</sub> = irreversible transfer to dark impurity state (typically 10<sup>5</sup>× faster)"))

        box, g = make_group("Impurity Quenching Parameters")
        self.enable  = QCheckBox("Enable impurity quenching  (Module 4B)")
        self.frac    = make_dspin(0.01, 0.0, 1.0, 4, 0.001,
                                  "Fraction of sites that are quenchers\n0.01 = 1% impurity concentration")
        self.k_q     = make_le("1e10",
                               "Quenching rate k_quench (s⁻¹)\nTypical range: 10⁹ – 10¹² s⁻¹", 180)
        self.sp_idx  = make_le("None",
                               "Species index of quenchers (integer) or None for random placement\n"
                               "None = random | 0, 1, 2 = specific species in blend", 180)
        self.on_hop  = QCheckBox(
            "Instantaneous quench on hop  (recommended: OFF — uses BKL event list)")

        re = g.rowCount(); g.addWidget(self.enable, re, 1, 1, 2)
        add_row(g, "Quencher fraction",        self.frac,   "",  "e.g. 0.01 = 1% impurity")
        add_row(g, "k<sub>quench</sub>",       self.k_q,    "s⁻¹")
        add_row(g, "Species index",            self.sp_idx, "",
                "None = random  |  0, 1, 2 = species in blend")
        roh = g.rowCount(); g.addWidget(self.on_hop, roh, 1, 1, 2)
        vbox.addWidget(box); vbox.addStretch()
        lay = QVBoxLayout(self); lay.addWidget(scrolled(inner))

    def get_values(self):
        s = self.sp_idx.text().strip()
        return {
            'enable_quenching'    : self.enable.isChecked(),
            'quencher_fraction'   : self.frac.value(),
            'k_quench'            : float(self.k_q.text()),
            'quencher_species_idx': None if s.lower() == 'none' else int(s),
            'quench_on_hop'       : self.on_hop.isChecked(),
        }

    def set_values(self, d):
        if 'enable_quenching'    in d: self.enable.setChecked(bool(d['enable_quenching']))
        if 'quencher_fraction'   in d: self.frac.setValue(float(d['quencher_fraction']))
        if 'k_quench'            in d: self.k_q.setText(str(d['k_quench']))
        if 'quencher_species_idx'in d:
            v = d['quencher_species_idx']
            self.sp_idx.setText('None' if v is None else str(v))
        if 'quench_on_hop'       in d: self.on_hop.setChecked(bool(d['quench_on_hop']))


# =============================================================================
# MAIN WINDOW
# =============================================================================
