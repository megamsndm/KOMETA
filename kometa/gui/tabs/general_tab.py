# kometa/gui/tabs/general_tab.py
import os
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QCheckBox,
    QFileDialog, QGroupBox, QScrollArea, QFrame
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui  import QFont

try:
    from ..style import *
except ImportError:
    from kometa.gui.style import *


class GeneralTab(QWidget):
    def __init__(self):
        super().__init__()
        inner = QWidget(); vbox = QVBoxLayout(inner); vbox.setSpacing(18)

        # ── KMC Settings ──────────────────────────────────────────────
        box, g = make_group("KMC Settings")
        self.num_kmc_steps = make_spin(100000, 1000, 10**8, "KMC steps per run")
        self.num_runs      = make_spin(10, 1, 500, "Independent runs")
        self.kmc_algorithm = QComboBox()
        self.kmc_algorithm.addItems([
            'bkl', 'gillespie', 'frm', 'nrm',
            'random_selection', 'rejection', 'binary_tree', 'akmc', 'pegasus'
        ])
        self.kmc_algorithm.setToolTip(
            "bkl           O(N)     Standard BKL — default for all systems\n"
            "gillespie     O(N)     Identical to BKL, CME derivation basis\n"
            "frm           O(logN)  First Reaction Method — multi-particle\n"
            "nrm           O(logN)  Next Reaction Method — TTA / bipolar\n"
            "random_select O(1)     Global sampling — dense multi-particle\n"
            "rejection     O(1)     Fast only for uniform rate landscapes\n"
            "binary_tree   O(logN)  Fenwick tree — best for large systems\n"
            "akmc          Var.     Accelerated — deep traps / low T\n"
            "pegasus       O(1-N)   KOMETA original — ordered/biased systems"
        )
        self.temperature = make_dspin(300.0, 1.0, 1000.0, 1, 10.0, "Temperature (K)")
        self.random_seed = make_spin(42, 0, 99999, "Random seed (0 = random)")
        add_row(g, "KMC steps per run", self.num_kmc_steps, "steps")
        add_row(g, "Number of runs",    self.num_runs,      "runs")
        add_row(g, "KMC algorithm",     self.kmc_algorithm)
        add_row(g, "Temperature",       self.temperature,   "K")
        add_row(g, "Random seed",       self.random_seed)
        vbox.addWidget(box)

        # ── Calculation Types ──────────────────────────────────────────
        box2, g2 = make_group("Calculations to Run")
        self.calc_charge    = QCheckBox("Charge Transport  (Marcus / MLJ rates, mobility)")
        self.calc_diffusion = QCheckBox("Triplet Diffusion  (D_T, diffusion length L_diff)")
        self.calc_tta       = QCheckBox("TTA Upconversion  (Phi_UC, annihilation statistics)")
        self.calc_charge.setChecked(True)
        self.calc_diffusion.setChecked(False)
        self.calc_tta.setChecked(False)
        self.calc_charge.setFont(QFont("Palatino Linotype", 13))
        self.calc_diffusion.setFont(QFont("Palatino Linotype", 13))
        self.calc_tta.setFont(QFont("Palatino Linotype", 13))
        r0 = g2.rowCount(); g2.addWidget(self.calc_charge,    r0, 0, 1, 3)
        r1 = g2.rowCount(); g2.addWidget(self.calc_diffusion, r1, 0, 1, 3)
        r2 = g2.rowCount(); g2.addWidget(self.calc_tta,       r2, 0, 1, 3)
        vbox.addWidget(box2)

        # ── Output Plots to Generate ───────────────────────────────────
        box3, g3 = make_group("Output Figures")
        self.plt_dos       = QCheckBox("Plot 1  — Charge DOS  (crystalline vs grain boundary)")
        self.plt_traj      = QCheckBox("Plot 2  — Carrier trajectory  (displacement vs time)")
        self.plt_layer_occ = QCheckBox("Plot 3  — Layer occupancy  (time in each layer)")
        self.plt_energy    = QCheckBox("Plot 4  — Carrier energy vs time")
        self.plt_mob_dist  = QCheckBox("Plot 5  — Mobility distribution  (histogram)")
        self.plt_layer_res = QCheckBox("Plot 6  — Layer residence fractions  (bar chart)")
        self.plt_domain    = QCheckBox("Plot 6b — Domain-resolved mobility  (cryst vs GB)")
        self.plt_tta_pie   = QCheckBox("Plot 7  — TTA outcome pie chart")
        self.plt_tta_time  = QCheckBox("Plot 8  — Annihilation time distribution")
        self.plt_tta_layer = QCheckBox("Plot 9  — TTA interlayer / GB breakdown")
        self.plt_phi_uc    = QCheckBox("Plot 10 -- Phi_UC step-by-step breakdown")
        self.plt_diffusion = QCheckBox("Plot 11 — Triplet diffusion coefficient")
        self.plt_zerofield = QCheckBox("Plot 12 — Zero-field vs finite-field mobility")

        _charge_plots  = [self.plt_dos, self.plt_traj, self.plt_layer_occ,
                          self.plt_energy, self.plt_mob_dist, self.plt_layer_res,
                          self.plt_domain, self.plt_zerofield]
        _tta_plots     = [self.plt_tta_pie, self.plt_tta_time,
                          self.plt_tta_layer, self.plt_phi_uc]
        _diff_plots    = [self.plt_diffusion]

        # Default: all charge plots on, TTA/diffusion off
        for cb in _charge_plots:  cb.setChecked(True)
        for cb in _tta_plots:     cb.setChecked(False)
        for cb in _diff_plots:    cb.setChecked(False)

        all_plot_checks = [
            self.plt_dos, self.plt_traj, self.plt_layer_occ, self.plt_energy,
            self.plt_mob_dist, self.plt_layer_res, self.plt_domain,
            self.plt_tta_pie, self.plt_tta_time, self.plt_tta_layer,
            self.plt_phi_uc, self.plt_diffusion, self.plt_zerofield,
        ]
        for cb in all_plot_checks:
            cb.setFont(QFont("Palatino Linotype", 12))
            row = g3.rowCount(); g3.addWidget(cb, row, 0, 1, 3)

        # Select-all / none buttons
        sel_row = QHBoxLayout()
        btn_all  = QPushButton("Select all");  btn_all.setObjectName("browse_btn")
        btn_none = QPushButton("Select none"); btn_none.setObjectName("browse_btn")
        btn_all.clicked.connect(lambda: [cb.setChecked(True)  for cb in all_plot_checks])
        btn_none.clicked.connect(lambda: [cb.setChecked(False) for cb in all_plot_checks])
        sel_row.addWidget(btn_all); sel_row.addWidget(btn_none); sel_row.addStretch()
        rsel = g3.rowCount(); g3.addLayout(sel_row, rsel, 0, 1, 3)
        vbox.addWidget(box3)

        # ── Output Folder ──────────────────────────────────────────────
        box4, g4 = make_group("Output")
        self.output_folder   = make_le("", "Output folder — leave blank to auto-name from molecule", 280)
        self.save_trajectory = QCheckBox("Save trajectory data")
        self.save_trajectory.setChecked(True)
        browse = QPushButton("Browse…"); browse.setObjectName("browse_btn")
        browse.clicked.connect(self._browse)
        r = g4.rowCount()
        g4.addWidget(rich("Output folder"), r, 0)
        g4.addWidget(self.output_folder, r, 1)
        g4.addWidget(browse, r, 2)
        note = QLabel("  Leave blank → folder named after molecule automatically")
        note.setStyleSheet(f"color: {C_SILVER}; font-size: 11pt; font-style: italic;")
        rn = g4.rowCount(); g4.addWidget(note, rn, 1, 1, 2)
        rt = g4.rowCount(); g4.addWidget(self.save_trajectory, rt, 1)
        vbox.addWidget(box4)

        vbox.addStretch()
        lay = QVBoxLayout(self); lay.addWidget(scrolled(inner))

        # Keep reference lists for engine integration
        self._all_plot_checks = all_plot_checks

    def _browse(self):
        d = QFileDialog.getExistingDirectory(self, "Select output folder")
        if d: self.output_folder.setText(d)

    def _build_simulation_mode(self):
        """Derive simulation_mode string from checkboxes."""
        c = self.calc_charge.isChecked()
        d = self.calc_diffusion.isChecked()
        t = self.calc_tta.isChecked()
        if c and not d and not t:   return 'charge'
        if not c and not d and t:   return 'TTA'
        if not c and d and not t:   return 'diffusion'
        if c and t:                 return 'both'
        if d or t or c:             return 'full'
        return 'charge'  # fallback

    def get_values(self):
        return {
            'num_kmc_steps'     : self.num_kmc_steps.value(),
            'num_runs'          : self.num_runs.value(),
            'kmc_algorithm'     : self.kmc_algorithm.currentText(),
            'temperature'       : self.temperature.value(),
            'random_seed'       : self.random_seed.value(),
            'output_folder'     : self.output_folder.text().strip(),
            'save_trajectory'   : self.save_trajectory.isChecked(),
            # Calculation selection
            'calc_charge'       : self.calc_charge.isChecked(),
            'calc_diffusion'    : self.calc_diffusion.isChecked(),
            'calc_tta'          : self.calc_tta.isChecked(),
            'simulation_mode'   : self._build_simulation_mode(),
            # Plot selection — individual flags
            'plot_dos'          : self.plt_dos.isChecked(),
            'plot_trajectory'   : self.plt_traj.isChecked(),
            'plot_layer_occ'    : self.plt_layer_occ.isChecked(),
            'plot_energy'       : self.plt_energy.isChecked(),
            'plot_mob_dist'     : self.plt_mob_dist.isChecked(),
            'plot_layer_res'    : self.plt_layer_res.isChecked(),
            'plot_domain'       : self.plt_domain.isChecked(),
            'plot_tta_pie'      : self.plt_tta_pie.isChecked(),
            'plot_tta_time'     : self.plt_tta_time.isChecked(),
            'plot_tta_layer'    : self.plt_tta_layer.isChecked(),
            'plot_phi_uc'       : self.plt_phi_uc.isChecked(),
            'plot_diffusion'    : self.plt_diffusion.isChecked(),
            'plot_zerofield'    : self.plt_zerofield.isChecked(),
        }

    def set_values(self, d):
        if 'num_kmc_steps'   in d: self.num_kmc_steps.setValue(int(d['num_kmc_steps']))
        if 'num_runs'        in d: self.num_runs.setValue(int(d['num_runs']))
        if 'kmc_algorithm'   in d: self.kmc_algorithm.setCurrentText(d['kmc_algorithm'])
        if 'temperature'     in d: self.temperature.setValue(float(d['temperature']))
        if 'random_seed'     in d: self.random_seed.setValue(int(d['random_seed']))
        if 'output_folder'   in d: self.output_folder.setText(str(d['output_folder']))
        if 'save_trajectory' in d: self.save_trajectory.setChecked(bool(d['save_trajectory']))
        if 'calc_charge'     in d: self.calc_charge.setChecked(bool(d['calc_charge']))
        if 'calc_diffusion'  in d: self.calc_diffusion.setChecked(bool(d['calc_diffusion']))
        if 'calc_tta'        in d: self.calc_tta.setChecked(bool(d['calc_tta']))
        m = {
            'plot_dos': self.plt_dos, 'plot_trajectory': self.plt_traj,
            'plot_layer_occ': self.plt_layer_occ, 'plot_energy': self.plt_energy,
            'plot_mob_dist': self.plt_mob_dist, 'plot_layer_res': self.plt_layer_res,
            'plot_domain': self.plt_domain, 'plot_tta_pie': self.plt_tta_pie,
            'plot_tta_time': self.plt_tta_time, 'plot_tta_layer': self.plt_tta_layer,
            'plot_phi_uc': self.plt_phi_uc, 'plot_diffusion': self.plt_diffusion,
            'plot_zerofield': self.plt_zerofield,
        }
        for k, cb in m.items():
            if k in d: cb.setChecked(bool(d[k]))
