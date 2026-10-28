# kometa/gui/tabs/blend_tab.py
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

class BlendTab(QWidget):
    def __init__(self):
        super().__init__()
        inner = QWidget(); vbox = QVBoxLayout(inner); vbox.setSpacing(18)
        box, g = make_group("Multi-Component Blend")
        self.blend   = QComboBox(); self.blend.addItems(['single', 'binary', 'ternary'])
        self.fracs   = make_le("1.0",
                               "Mole fractions — must sum to 1.0\n"
                               "Single: 1.0 | Binary: 0.9, 0.1 | Ternary: 0.8, 0.15, 0.05", 240)
        self.sp_T1   = make_le("1.8",
                               "T₁ per species (eV) — one value per component from ORCA TDDFT", 240)
        self.sp_kdec = make_le("1e5",
                               "k_decay per species (s⁻¹) — one per component from TRPL", 240)
        self.inter_h = QCheckBox("Allow inter-species triplet hopping")

        add_row(g, "Blend type",             self.blend)
        add_row(g, "Mole fractions",         self.fracs,   "Σ = 1")
        add_row(g, "T<sub>1</sub> / species",self.sp_T1,   "eV",
                "One value per species — from ORCA TDDFT")
        add_row(g, "k<sub>decay</sub> / species", self.sp_kdec, "s⁻¹",
                "One value per species — from TRPL")
        rh = g.rowCount(); g.addWidget(self.inter_h, rh, 1, 1, 2)
        vbox.addWidget(box); vbox.addStretch()
        lay = QVBoxLayout(self); lay.addWidget(scrolled(inner))

    def get_values(self):
        return {
            'blend_type'                : self.blend.currentText(),
            'blend_fractions'           : parse_list(self.fracs.text()),
            'species_triplet_energy'    : parse_list(self.sp_T1.text()),
            'species_k_decay'           : parse_list(self.sp_kdec.text()),
            'allow_inter_species_hopping': self.inter_h.isChecked(),
        }

    def set_values(self, d):
        if 'blend_type'                  in d: self.blend.setCurrentText(d['blend_type'])
        if 'blend_fractions'             in d: self.fracs.setText(l2s(d['blend_fractions']))
        if 'species_triplet_energy'      in d: self.sp_T1.setText(l2s(d['species_triplet_energy']))
        if 'species_k_decay'             in d: self.sp_kdec.setText(l2s(d['species_k_decay']))
        if 'allow_inter_species_hopping' in d: self.inter_h.setChecked(bool(d['allow_inter_species_hopping']))


