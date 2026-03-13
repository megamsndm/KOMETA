# kometa/gui/tabs/tta_tab.py
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

class TTATab(QWidget):
    def __init__(self):
        super().__init__()
        inner = QWidget(); vbox = QVBoxLayout(inner); vbox.setSpacing(18)
        vbox.addWidget(info_lbl(
            "Triplet energy T<sub>1</sub> and emission yield Φ<sub>em</sub> from ORCA TDDFT.<br>"
            "k<sub>decay</sub> = 1/τ<sub>triplet</sub> from time-resolved PL measurement."))

        box, g = make_group("Triplet & Dexter Parameters  (from ORCA TDDFT)")
        self.T1 = make_le("1.8, 1.75",
                          "Triplet energy T₁ per layer (eV) — TDDFT excitation energy", 240)
        self.k0 = make_le("1e12", "Dexter transfer prefactor k₀ (s⁻¹)\nTypical: 10¹¹ – 10¹³ s⁻¹", 180)
        self.L  = make_dspin(1.5, 0.1, 20.0, 2, 0.1, "Dexter decay length L (Å)\nTypical: 1–3 Å for organic molecules")

        add_row(g, "Triplet energy T<sub>1</sub>",
                self.T1, "eV / layer", "From ORCA TDDFT — T₁ excitation energy")
        add_row(g, "Dexter prefactor k<sub>0</sub>",
                self.k0, "s⁻¹")
        add_row(g, "Dexter decay length L",
                self.L, "Å", "k<sub>Dexter</sub> ∝ exp(−2r/L)")
        vbox.addWidget(box)

        box2, g2 = make_group("TTA Kinetics")
        self.k_TTA   = make_le("1e9", "Bimolecular TTA rate constant (s⁻¹)\nFrom intensity-dependent PL", 180)
        self.k_decay = make_le("1e5",
                               "Natural triplet decay rate k_decay = 1/τ (s⁻¹)\nFrom TRPL measurement", 180)
        self.f_s     = make_dspin(0.1111, 0.0, 1.0, 4, 0.001,
                                  "Singlet fraction fₛ = 1/9 (spin statistics for TTA)\n9 spin sub-levels: 1 singlet + 8 non-emissive")
        self.phi_em  = make_dspin(0.8, 0.0, 1.0, 3, 0.01,
                                  "Emission quantum yield Φₑₘ\nFrom TDDFT oscillator strength or PLQY measurement")
        self.n_pairs = make_spin(200, 10, 50000, "Number of triplet pairs to simulate")
        self.n_diff  = make_spin(100, 10, 10000,
                                 "Number of single excitons for D_T measurement")

        add_row(g2, "k<sub>TTA</sub>",
                self.k_TTA,   "s⁻¹")
        add_row(g2, "k<sub>decay</sub>  (= 1/τ)",
                self.k_decay, "s⁻¹", "From TRPL: k<sub>decay</sub> = 1/τ<sub>triplet</sub>")
        add_row(g2, "Singlet fraction f<sub>s</sub>",
                self.f_s,     "", "= 1/9 from spin statistics")
        add_row(g2, "Emission yield Φ<sub>em</sub>",
                self.phi_em)
        add_row(g2, "TTA triplet pairs",
                self.n_pairs, "pairs")
        add_row(g2, "Excitons for D<sub>T</sub>",
                self.n_diff,  "excitons")
        vbox.addWidget(box2); vbox.addStretch()
        lay = QVBoxLayout(self); lay.addWidget(scrolled(inner))

    def get_values(self):
        return {
            'triplet_energy'        : parse_list(self.T1.text()),
            'k_dexter_prefactor'    : float(self.k0.text()),
            'dexter_decay_length'   : self.L.value(),
            'k_TTA'                 : float(self.k_TTA.text()),
            'k_decay'               : float(self.k_decay.text()),
            'singlet_fraction'      : self.f_s.value(),
            'phi_emission'          : self.phi_em.value(),
            'num_triplet_pairs'     : self.n_pairs.value(),
            'diffusion_num_excitons': self.n_diff.value(),
        }

    def set_values(self, d):
        if 'triplet_energy'         in d: self.T1.setText(l2s(d['triplet_energy']))
        if 'k_dexter_prefactor'     in d: self.k0.setText(str(d['k_dexter_prefactor']))
        if 'dexter_decay_length'    in d: self.L.setValue(float(d['dexter_decay_length']))
        if 'k_TTA'                  in d: self.k_TTA.setText(str(d['k_TTA']))
        if 'k_decay'                in d: self.k_decay.setText(str(d['k_decay']))
        if 'singlet_fraction'       in d: self.f_s.setValue(float(d['singlet_fraction']))
        if 'phi_emission'           in d: self.phi_em.setValue(float(d['phi_emission']))
        if 'num_triplet_pairs'      in d: self.n_pairs.setValue(int(d['num_triplet_pairs']))
        if 'diffusion_num_excitons' in d: self.n_diff.setValue(int(d['diffusion_num_excitons']))


