# kometa/gui/tabs/morphology_tab.py
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

class MorphologyTab(QWidget):
    def __init__(self):
        super().__init__()
        inner = QWidget(); vbox = QVBoxLayout(inner); vbox.setSpacing(18)

        # mol_res kept as a hidden attribute to sync with header bar
        self.mol_res = make_le("PtOEP", "Synced with Molecule bar at top", 180)

        box2, g2 = make_group("Grid Geometry  (used when source = grid)")
        self.n_layers = make_spin(2, 1, 50, "Number of molecular layers")
        self.spl      = make_spin(50, 1, 5000, "Sites per layer")
        self.d_intra  = make_dspin(4.0, 0.5, 50.0, 2, 0.1, "Intra-layer distance (Å)")
        self.d_inter  = make_dspin(3.5, 0.5, 50.0, 2, 0.1, "Inter-layer stacking distance (Å)")
        self.pos_dis  = make_dspin(0.3, 0.0, 5.0, 2, 0.05, "Positional disorder σ (Å); 0 = off")
        add_row(g2, "Number of layers",              self.n_layers, "layers")
        add_row(g2, "Sites per layer",               self.spl,      "sites")
        add_row(g2, "Intra-layer distance",          self.d_intra,  "Å")
        add_row(g2, "Inter-layer distance",          self.d_inter,  "Å")
        add_row(g2, "Positional disorder σ",         self.pos_dis,  "Å")
        vbox.addWidget(box2)

        box3, g3 = make_group("Grain Boundaries")
        self.detect_gb = QCheckBox("Enable grain boundary detection  (S₂ order parameter)")
        self.detect_gb.setChecked(True)
        self.gb_cut = make_dspin(0.6, 0.0, 1.0, 2, 0.05, "S₂ cutoff — sites below are grain boundaries")
        self.gb_J   = make_dspin(0.3, 0.01, 2.0, 2, 0.05, "Transfer integral J multiplier at GB")
        self.gb_lam = make_dspin(1.5, 0.5, 5.0, 2, 0.1,   "Reorganization energy λ multiplier at GB")
        self.gb_sig = make_dspin(1.5, 0.5, 5.0, 2, 0.1,   "Disorder σ multiplier at GB")
        r0 = g3.rowCount(); g3.addWidget(self.detect_gb, r0, 1, 1, 2)
        add_row(g3, "S₂ order cutoff",            self.gb_cut)
        add_row(g3, "GB  J factor",               self.gb_J,   "× J")
        add_row(g3, "GB  λ factor",               self.gb_lam, "× λ")
        add_row(g3, "GB  σ factor",               self.gb_sig, "× σ")
        vbox.addWidget(box3)
        vbox.addStretch()
        lay = QVBoxLayout(self); lay.addWidget(scrolled(inner))



    def get_values(self):
        return {
            'morphology_source'        : 'grid',
            'molecule_name'            : self.mol_res.text(),
            'num_layers'               : self.n_layers.value(),
            'sites_per_layer'          : self.spl.value(),
            'intra_layer_distance'     : self.d_intra.value(),
            'inter_layer_distance'     : self.d_inter.value(),
            'positional_disorder_sigma': self.pos_dis.value(),
            'detect_grain_boundaries'  : self.detect_gb.isChecked(),
            'domain_order_cutoff'      : self.gb_cut.value(),
            'gb_J_factor'              : self.gb_J.value(),
            'gb_lambda_factor'         : self.gb_lam.value(),
            'gb_sigma_factor'          : self.gb_sig.value(),
        }

    def set_values(self, d):
        if 'molecule_name'             in d: self.mol_res.setText(str(d['molecule_name']))
        if 'num_layers'                in d: self.n_layers.setValue(int(d['num_layers']))
        if 'sites_per_layer'           in d: self.spl.setValue(int(d['sites_per_layer']))
        if 'intra_layer_distance'      in d: self.d_intra.setValue(float(d['intra_layer_distance']))
        if 'inter_layer_distance'      in d: self.d_inter.setValue(float(d['inter_layer_distance']))
        if 'positional_disorder_sigma' in d: self.pos_dis.setValue(float(d['positional_disorder_sigma']))
        if 'detect_grain_boundaries'   in d: self.detect_gb.setChecked(bool(d['detect_grain_boundaries']))
        if 'domain_order_cutoff'       in d: self.gb_cut.setValue(float(d['domain_order_cutoff']))
        if 'gb_J_factor'               in d: self.gb_J.setValue(float(d['gb_J_factor']))
        if 'gb_lambda_factor'          in d: self.gb_lam.setValue(float(d['gb_lambda_factor']))
        if 'gb_sigma_factor'           in d: self.gb_sig.setValue(float(d['gb_sigma_factor']))


