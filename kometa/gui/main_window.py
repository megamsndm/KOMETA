# kometa/gui/main_window.py
import sys, os, json, traceback, subprocess, platform, importlib.util
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
try:
    from .header_bar import HeaderBar
except ImportError:
    from kometa.gui.header_bar import HeaderBar
try:
    from .run_panel import RunPanel
except ImportError:
    from kometa.gui.run_panel import RunPanel
try:
    from .results_tab import ResultsTab
except ImportError:
    from kometa.gui.results_tab import ResultsTab
try:
    from .tabs.general_tab import GeneralTab
except ImportError:
    from kometa.gui.tabs.general_tab import GeneralTab
try:
    from .tabs.morphology_tab import MorphologyTab
except ImportError:
    from kometa.gui.tabs.morphology_tab import MorphologyTab
try:
    from .tabs.charge_tab import ChargeTab
except ImportError:
    from kometa.gui.tabs.charge_tab import ChargeTab
try:
    from .tabs.mlj_tab import MLJTab
except ImportError:
    from kometa.gui.tabs.mlj_tab import MLJTab
try:
    from .tabs.tta_tab import TTATab
except ImportError:
    from kometa.gui.tabs.tta_tab import TTATab
try:
    from .tabs.blend_tab import BlendTab
except ImportError:
    from kometa.gui.tabs.blend_tab import BlendTab
try:
    from .tabs.quench_tab import QuenchTab
except ImportError:
    from kometa.gui.tabs.quench_tab import QuenchTab
try:
    from ..io.config import save_config, load_config
except ImportError:
    from kometa.io.config import save_config, load_config

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("KOMETA")
        self.resize(1280, 880)
        self.setStyleSheet(QSS)
        # Menu bar — File and Help only
        self._build_menu()

        # Status bar
        self.sb = QStatusBar(); self.setStatusBar(self.sb)
        self.sb.showMessage(
"✔  KOMETA  ready")

        # Central layout
        central = QWidget(); self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # Header bar
        self.header = HeaderBar()
        root.addWidget(self.header)

        # Main splitter
        splitter = QSplitter(Qt.Horizontal)

        # ── LEFT: parameter tabs ─────────────────────────────────────────
        self.param_tabs = QTabWidget()
        self.tab_gen    = GeneralTab()
        self.tab_morph  = MorphologyTab()
        self.tab_charge = ChargeTab()
        self.tab_mlj    = MLJTab()
        self.tab_tta    = TTATab()
        self.tab_blend  = BlendTab()
        self.tab_quench = QuenchTab()

        self.param_tabs.addTab(self.tab_gen,    "⚙  General")
        self.param_tabs.addTab(self.tab_morph,  "🔬  Morphology")
        self.param_tabs.addTab(self.tab_charge, "⚡  Charge")
        self.param_tabs.addTab(self.tab_mlj,    "🌡  MLJ")
        self.param_tabs.addTab(self.tab_tta,    "✨  TTA")
        self.param_tabs.addTab(self.tab_blend,  "🧪  Blend")
        self.param_tabs.addTab(self.tab_quench, "☢  Quench")
        self.param_tabs.setMinimumWidth(560)

        # ── RIGHT: run panel + results ───────────────────────────────────
        right_tabs = QTabWidget()
        self.run_panel   = RunPanel(self)
        self.results_tab = ResultsTab()
        self.results_tab.ref_btn.clicked.connect(self.refresh_results)
        right_tabs.addTab(self.run_panel,   "▶  Run")
        right_tabs.addTab(self.results_tab, "📊  Results")
        right_tabs.setMinimumWidth(520)

        splitter.addWidget(self.param_tabs)
        splitter.addWidget(right_tabs)
        splitter.setSizes([600, 640])
        root.addWidget(splitter)

        # Sync molecule name: header bar ↔ morphology tab (block signals to avoid loops)
        self.header.mol_name.textChanged.connect(self._sync_mol_from_header)
        self.tab_morph.mol_res.textChanged.connect(self._sync_mol_from_morph)
        self._syncing = False

    def _sync_mol_from_header(self, t):
        if not self._syncing:
            self._syncing = True
            self.tab_morph.mol_res.setText(t)
            self._syncing = False

    def _sync_mol_from_morph(self, t):
        if not self._syncing:
            self._syncing = True
            self.header.mol_name.setText(t)
            self._syncing = False

    # ── Menu ──────────────────────────────────────────────────────────────
    def _build_menu(self):
        mb = self.menuBar()

        def act(name, slot, sc=None):
            a = QAction(name, self); a.triggered.connect(slot)
            if sc: a.setShortcut(sc)
            return a

        # File
        fm = mb.addMenu("&File")
        fm.addAction(act("&Save Configuration…", self.save_config, "Ctrl+S"))
        fm.addAction(act("&Load Configuration…", self.load_config, "Ctrl+O"))
        fm.addSeparator()
        fm.addAction(act("Open Output &Folder",  self._open_folder))
        fm.addSeparator()
        fm.addAction(act("E&xit", self.close, "Ctrl+Q"))

        # Help
        hm = mb.addMenu("&Help")
        hm.addAction(act("&Introduction",  self._show_intro))
        hm.addAction(act("&About KOMETA",  self._about))

    # ── Collect params ─────────────────────────────────────────────────────
    def collect_params(self):
        p = {}
        for tab in [self.tab_gen, self.tab_morph, self.tab_charge,
                    self.tab_mlj, self.tab_tta,   self.tab_blend,
                    self.tab_quench]:
            p.update(tab.get_values())
        mol = self.header.mol_name.text().strip() or 'molecule'
        p['molecule_name'] = mol
        # Auto-set output folder to molecule name if left blank
        if not p.get('output_folder', '').strip():
            p['output_folder'] = mol
        p.setdefault('site_energy_model',   'gaussian')
        p.setdefault('J_source',            'uniform')
        p.setdefault('tta_max_time',        1e-3)
        p.setdefault('tta_max_steps',       int(5e4))
        p.setdefault('diffusion_max_steps', int(1e5))
        p.setdefault('molecular_dipole',    2.0)
        return p

    # ── Results ───────────────────────────────────────────────────────────
    def refresh_results(self):
        p = self.collect_params()
        folder = p.get('output_folder', '') or p.get('molecule_name', 'molecule')
        self.results_tab.load_plots(folder)

    # ── Config ────────────────────────────────────────────────────────────
    def save_config(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save configuration", "kometa_config.json", "JSON (*.json)")
        if not path: return
        with open(path, 'w') as f:
            json.dump(self.collect_params(), f, indent=2, default=str)
        self.sb.showMessage(f"✔  Configuration saved: {path}", 6000)

    def load_config(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load configuration", "", "JSON (*.json)")
        if not path: return
        with open(path) as f: d = json.load(f)
        for tab in [self.tab_gen, self.tab_morph, self.tab_charge,
                    self.tab_mlj, self.tab_tta,   self.tab_blend,
                    self.tab_quench]:
            tab.set_values(d)
        if 'molecule_name' in d:
            self.header.mol_name.setText(str(d['molecule_name']))
        self.sb.showMessage(f"✔  Configuration loaded: {path}", 6000)

    def _open_folder(self):
        folder = self.tab_gen.output_folder.text()
        if os.path.isdir(folder):
            if platform.system() == 'Windows':  os.startfile(folder)
            elif platform.system() == 'Darwin': subprocess.Popen(['open', folder])
            else:                               subprocess.Popen(['xdg-open', folder])
        else:
            QMessageBox.information(self, "Folder not found",
                f"'{folder}' does not exist yet.\nRun a simulation first.")


    # ── Help dialogs ───────────────────────────────────────────────────────
    def _show_intro(self):
        dlg = QDialog(self); dlg.setWindowTitle("Introduction to KOMETA")
        dlg.resize(700, 600)
        lay = QVBoxLayout(dlg)
        txt = QTextEdit(); txt.setReadOnly(True)
        txt.setStyleSheet(f"""
            background: {C_PANEL}; color: {C_WHITE};
            font-family: 'Palatino Linotype', serif; font-size: 13pt;
            border: none; padding: 10px;
        """)
        txt.setHtml(f"""
        <h2 style='color:{C_TEAL}; letter-spacing:3px;'>KOMETA</h2>
        <p style='color:{C_SILVER}; font-size:11pt; font-style:italic;'>
        KMC Organic Molecular Exciton Transport Analysis — Version 1.0</p>

        <p style='margin-top:16px;'>
        <b style='color:{C_GOLD};'>Developer:</b>
        <span style='color:{C_WHITE};'> Dr. Mohsen Doust Mohammadi</span></p>

        <hr style='border:1px solid {C_BORDER}; margin:14px 0;'/>

        <h3 style='color:{C_TEAL};'>What is KOMETA?</h3>
        <p>KOMETA is a professional Kinetic Monte Carlo (KMC) simulation platform
        for organic molecular semiconductors and emitters. It integrates state-of-the-art
        theoretical methods to model charge transport, triplet exciton diffusion, and
        triplet-triplet annihilation upconversion (TTA-UC) in layered molecular systems.</p>

        <h3 style='color:{C_TEAL}; margin-top:16px;'>Physical Modules</h3>
        <ul style='line-height:1.9;'>
        <li><b style='color:{C_GOLD};'>Marcus / MLJ Rates</b> — Classical and quantum-corrected
            charge hopping rates. MLJ extends Marcus theory to low temperatures via
            Franck-Condon weighted nuclear tunnelling.</li>
        <li><b style='color:{C_GOLD};'>Dexter Energy Transfer</b> — Short-range triplet exciton
            hopping governed by orbital overlap.</li>
        <li><b style='color:{C_GOLD};'>BKL Algorithm</b> — Bortz-Kalos-Lebowitz rejection-free
            KMC for efficient event selection.</li>
        <li><b style='color:{C_GOLD};'>TTA Upconversion</b> — Two-body triplet-triplet annihilation
            kinetics, singlet yield (f<sub>s</sub> = 1/9), and Φ<sub>UC</sub> calculation.</li>
        <li><b style='color:{C_GOLD};'>Triplet Diffusion (D<sub>T</sub>)</b> — Single-exciton
            mean-square displacement and diffusion length L<sub>diff</sub>.</li>
        <li><b style='color:{C_GOLD};'>Grain Boundaries</b> — S₂ order parameter detection,
            domain-resolved mobility.</li>
        <li><b style='color:{C_GOLD};'>Impurity Quenching</b> — Module 4B, irreversible
            Dexter-type energy transfer to dark impurity states
            (DOI: 10.1103/PhysRevApplied.14.034050).</li>
        <li><b style='color:{C_GOLD};'>Multi-component Blends</b> — Binary and ternary systems
            with species-specific energetics.</li>
        </ul>

        <h3 style='color:{C_TEAL}; margin-top:16px;'>Workflow</h3>
        <ol style='line-height:1.9;'>
        <li>Enter your molecule name in the green box at the top.</li>
        <li>Run ORCA calculations externally to obtain λ, J, T₁, S, ℏω, Φ<sub>em</sub>.</li>
        <li>Enter ORCA values into the parameter tabs on the left.</li>
        <li>Save your configuration:  <i>File → Save Configuration</i>.</li>
        <li>Select a simulation mode and press <b>▶ RUN SIMULATION</b>.</li>
        <li>View output plots in the Results tab.</li>
        </ol>

        <p style='margin-top:16px; color:{C_SILVER}; font-size:11pt;'>
        📖  For detailed descriptions of every parameter, physical derivations,
        and interpretation of output plots, please refer to the
        <b style='color:{C_TEAL};'>KOMETA Manual</b>  (KOMETA_manual.pdf).</p>

        <hr style='border:1px solid {C_BORDER}; margin:14px 0;'/>
        <p style='color:{C_SILVER}; font-size:10pt; text-align:center;'>
        KOMETA  ·  Dr. Mohsen Doust Mohammadi  ·  All rights reserved</p>
        """)
        lay.addWidget(txt)
        close_btn = QPushButton("Close"); close_btn.clicked.connect(dlg.close)
        lay.addWidget(close_btn)
        dlg.exec_()

    def _about(self):
        QMessageBox.about(self, "About KOMETA",
            f"<h2 style='color:{C_TEAL}; letter-spacing:3px;'>KOMETA</h2>"
            f"<p style='color:{C_SILVER}; font-style:italic;'>"
            "KMC Organic Molecular Exciton Transport Analysis</p>"
            f"<p><b>Version:</b> 3.0</p>"
            f"<p><b>Developer:</b> Dr. Mohsen Doust Mohammadi</p>"
            f"<hr/>"
            f"<p>Modules: Marcus · MLJ · Dexter · BKL · TTA-UC ·"
            " Grain Boundaries · Impurity Quenching · Multi-component Blends</p>"
            f"<p>Works with <b>any organic molecule</b> — enter ORCA-computed"
            " parameters in the tabs.</p>"
            f"<p style='color:{C_SILVER}; font-size:10pt;'>"
            "📖  Refer to KOMETA_manual.pdf for full documentation.</p>")


# =============================================================================
# ENTRY POINT
# =============================================================================
def main():
    app = QApplication(sys.argv)
    app.setApplicationName("KOMETA")
    app.setStyle("Fusion")

    pal = QPalette()
    pal.setColor(QPalette.Window,          QColor(C_BG))
    pal.setColor(QPalette.WindowText,      QColor(C_WHITE))
    pal.setColor(QPalette.Base,            QColor(C_PANEL))
    pal.setColor(QPalette.AlternateBase,   QColor(C_BG))
