# kometa/gui/run_panel.py
import os, subprocess, platform
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QTextEdit, QProgressBar,
    QMessageBox, QSizePolicy
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui  import QFont

try:
    from .style import *
except ImportError:
    from kometa.gui.style import *
try:
    from .worker import SimWorker
except ImportError:
    from kometa.gui.worker import SimWorker


class RunPanel(QWidget):
    def __init__(self, parent_window):
        super().__init__()
        self.pw     = parent_window
        self.worker = None
        lay = QVBoxLayout(self)
        lay.setContentsMargins(20, 20, 20, 20)
        lay.setSpacing(16)

        # ── RUN / STOP BUTTONS ────────────────────────────────────────
        run_row = QHBoxLayout()
        self.run_btn = QPushButton("\u25b6   Run Simulation")
        self.run_btn.setObjectName("run_btn")
        self.run_btn.setCursor(Qt.PointingHandCursor)
        self.run_btn.setFixedHeight(48)
        self.stop_btn = QPushButton("\u25a0  Stop")
        self.stop_btn.setObjectName("stop_btn")
        self.stop_btn.setEnabled(False)
        self.stop_btn.setFixedWidth(110)
        run_row.addWidget(self.run_btn)
        run_row.addSpacing(12)
        run_row.addWidget(self.stop_btn)
        lay.addLayout(run_row)

        # ── ACTIVE CALCULATIONS (read-only info label) ────────────────
        self.active_lbl = QLabel("Calculations: —")
        self.active_lbl.setFont(QFont("Palatino Linotype", 11))
        self.active_lbl.setStyleSheet(f"color: {C_TEAL}; padding: 4px 0;")
        self.active_lbl.setWordWrap(True)
        lay.addWidget(self.active_lbl)

        # ── PROGRESS ──────────────────────────────────────────────────
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        self.prog_lbl = QLabel("")
        self.prog_lbl.setFont(QFont("Palatino Linotype", 12))
        self.prog_lbl.setStyleSheet(f"color: {C_GOLD};")
        pr = QHBoxLayout()
        pr.addWidget(self.progress); pr.addWidget(self.prog_lbl)
        lay.addLayout(pr)

        # ── CONSOLE ───────────────────────────────────────────────────
        con_lbl = QLabel("Console Output")
        con_lbl.setObjectName("section")
        con_lbl.setFont(QFont("Palatino Linotype", 14, QFont.Bold))
        lay.addWidget(con_lbl)
        self.console = QTextEdit()
        self.console.setReadOnly(True)
        self.console.setMinimumHeight(260)
        lay.addWidget(self.console, stretch=1)

        # ── BOTTOM BUTTONS ────────────────────────────────────────────
        bot = QHBoxLayout()
        self.save_btn = QPushButton("\U0001f4be  Save Config");  self.save_btn.setObjectName("browse_btn")
        self.load_btn = QPushButton("\U0001f4c2  Load Config");  self.load_btn.setObjectName("browse_btn")
        self.open_btn = QPushButton("\U0001f4c1  Open Results"); self.open_btn.setObjectName("browse_btn")
        bot.addWidget(self.save_btn); bot.addWidget(self.load_btn)
        bot.addWidget(self.open_btn); bot.addStretch()
        lay.addLayout(bot)

        # ── SIGNALS ───────────────────────────────────────────────────
        self.run_btn.clicked.connect(self._run)
        self.stop_btn.clicked.connect(self._stop)
        self.save_btn.clicked.connect(self.pw.save_config)
        self.load_btn.clicked.connect(self.pw.load_config)
        self.open_btn.clicked.connect(self._open_folder)

    # ------------------------------------------------------------------
    def _run(self):
        self.console.clear()
        params = self.pw.collect_params()

        mol  = params.get('molecule_name', 'molecule') or 'molecule'
        mode = params.get('simulation_mode', 'charge')

        # Build output folder: user override OR molecule name
        if not params.get('output_folder', '').strip():
            params['output_folder'] = mol
        params['molecule_name'] = mol

        # Update active calcs label
        calcs = []
        if params.get('calc_charge'):    calcs.append('Charge Transport')
        if params.get('calc_diffusion'): calcs.append('Triplet Diffusion')
        if params.get('calc_tta'):       calcs.append('TTA Upconversion')
        self.active_lbl.setText("Calculations: " + ("  |  ".join(calcs) if calcs else "none selected"))

        self.console.append(
            f"KOMETA  \u00b7  Molecule: {mol}  \u00b7  Mode: {mode}\n"
            f"Output folder: {params['output_folder']}\n"
            f"{'─'*55}\n"
        )
        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.progress.setVisible(True)
        self.prog_lbl.setText("Running\u2026")
        self.pw.header.set_running()

        self.worker = SimWorker(params)
        self.worker.log_signal.connect(self._log)
        self.worker.done_signal.connect(self._on_done)
        self.worker.start()

    def _stop(self):
        if self.worker and self.worker.isRunning():
            self.stop_btn.setEnabled(False)
            self.stop_btn.setText("Stopping\u2026")
            self.prog_lbl.setText("Stopping \u2014 please wait\u2026")
            self.worker.request_stop()
            if not self.worker.wait(8000):
                self.worker.quit(); self.worker.wait(2000)
            self.stop_btn.setText("\u25a0  Stop")

    def _log(self, text):
        self.console.append(text.rstrip())
        self.console.ensureCursorVisible()

    def _on_done(self, ok, msg):
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.progress.setVisible(False)
        if ok:
            self.prog_lbl.setText("\u2714  Done")
            self.prog_lbl.setStyleSheet(f"color: {C_GREEN}; font-weight: bold;")
            self.console.append(f"\n{'='*55}\n\u2714  {msg}\n")
            self.pw.header.set_done()
            self.pw.refresh_results()
        else:
            self.prog_lbl.setText("\u2718  " + msg[:70])
            self.prog_lbl.setStyleSheet(f"color: {C_RED}; font-weight: bold;")
            self.pw.header.set_error()

    def _open_folder(self):
        folder = self.pw.tab_gen.output_folder.text().strip()
        if not folder:
            mol = self.pw.header.mol_name.text().strip() or 'molecule'
            folder = mol
        if os.path.isdir(folder):
            if platform.system() == 'Windows':   os.startfile(folder)
            elif platform.system() == 'Darwin':  subprocess.Popen(['open', folder])
            else:                                 subprocess.Popen(['xdg-open', folder])
        else:
            QMessageBox.information(self, "Not found",
                f"Folder '{folder}' does not exist yet. Run a simulation first.")
