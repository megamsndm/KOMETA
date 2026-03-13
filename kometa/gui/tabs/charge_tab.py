# kometa/gui/tabs/charge_tab.py
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QCheckBox,
    QGroupBox, QScrollArea, QFrame, QStackedWidget
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui  import QFont

try:
    from ..style import *
except ImportError:
    from kometa.gui.style import *


# Model descriptions shown in the info box
_MODEL_INFO = {
    'marcus': (
        "Classical Marcus rate — gold standard for organic semiconductors.\n"
        "Requires: lambda (reorganization energy) and J (transfer integral) from ORCA.\n"
        "Formula:  k = (2pi/hbar)|J|^2 / sqrt(4pi*lambda*kT) * exp(-(dG+lambda)^2 / 4*lambda*kT)\n"
        "Best for: room temperature, moderate disorder, well-ordered films."
    ),
    'mlj': (
        "Marcus-Levich-Jortner — extends Marcus with quantum nuclear tunnelling.\n"
        "Adds Franck-Condon sum over vibrational levels via Huang-Rhys factor S.\n"
        "Requires: lambda_outer, S, hbar_omega from ORCA frequency calculation.\n"
        "Best for: low temperature, high-frequency C=C modes, triplet emitters."
    ),
    'miller_abrahams': (
        "Miller-Abrahams (1960) — simplest hopping model, no ORCA input needed.\n"
        "Rate determined by: spatial tunnelling exp(-2*alpha*r) and Boltzmann.\n"
        "No reorganization energy. Parameters: nu_0 (attempt freq) and alpha (localization).\n"
        "Best for: amorphous films, deeply localized states, benchmarking."
    ),
    'forster': (
        "Forster FRET — long-range dipole-dipole singlet exciton transfer.\n"
        "k_FRET = k_D * (R_0 / r)^6  where R_0 is the Forster radius.\n"
        "Requires: R_0 (Forster radius, nm) and k_D (singlet decay rate) from TDDFT.\n"
        "Best for: fluorescent OLEDs, host-guest singlet harvesting, TADF."
    ),
    'gdm': (
        "Gaussian Disorder Model (Bassler) — purely empirical, no ORCA needed.\n"
        "Rate = nu_0 * exp(-dE/kT) for uphill hops; nu_0 for downhill hops.\n"
        "Parameters: nu_0 (attempt frequency). Uses disorder_sigma for site energies.\n"
        "Best for: amorphous films, quick parameter scans, no quantum chemistry data."
    ),
    'cdm': (
        "Correlated Disorder Model — GDM with spatial exponential decay added.\n"
        "Rate = nu_0 * exp(-r/r_loc) * Boltzmann. Captures spatial correlations.\n"
        "Parameters: nu_0 and r_loc (localization radius, Ang).\n"
        "Best for: polar molecules with spatially correlated site energies."
    ),
    'small_polaron': (
        "Small Polaron / Holstein model — carrier self-trapping in soft lattice.\n"
        "Carrier distorts local geometry, lowering its own energy by E_b = S*hbar_omega.\n"
        "J is band-narrowed: J_eff = J * exp(-S/2). Rates use MLJ formula.\n"
        "Best for: high-lambda molecules (lambda > 0.5 eV), soft molecular crystals."
    ),
}

_MODELS = list(_MODEL_INFO.keys())
_MODEL_LABELS = [
    'marcus          — Marcus (classical)',
    'mlj             — Marcus-Levich-Jortner',
    'miller_abrahams — Miller-Abrahams',
    'forster         — Forster FRET (singlet)',
    'gdm             — Gaussian Disorder Model',
    'cdm             — Correlated Disorder Model',
    'small_polaron   — Small Polaron / Holstein',
]


class ChargeTab(QWidget):
    def __init__(self):
        super().__init__()
        inner = QWidget(); vbox = QVBoxLayout(inner); vbox.setSpacing(14)

        # ── Rate model selector ────────────────────────────────────────
        box0, g0 = make_group("Rate Model")
        self.rate_model = QComboBox()
        for lbl in _MODEL_LABELS:
            self.rate_model.addItem(lbl)
        self.rate_model.setFont(QFont("Palatino Linotype", 12))
        self.rate_model.setMinimumWidth(340)

        self.carrier = QComboBox()
        self.carrier.addItems(['hole', 'electron'])

        add_row(g0, "Rate model",   self.rate_model)
        add_row(g0, "Carrier type", self.carrier)

        # Dynamic info label
        self.model_info = QLabel()
        self.model_info.setWordWrap(True)
        self.model_info.setStyleSheet(
            f"color: {C_SILVER}; font-size: 11pt; font-style: italic;"
            f" background: {C_PANEL}; padding: 8px; border-radius: 4px;")
        self.model_info.setFont(QFont("Palatino Linotype", 11))
        ri = g0.rowCount(); g0.addWidget(self.model_info, ri, 0, 1, 3)
        vbox.addWidget(box0)
        self.rate_model.currentIndexChanged.connect(self._update_model_info)

        # ── Universal parameters (all models) ─────────────────────────
        box1, g1 = make_group("Site Energies & Disorder  (all models)")
        self.mean_E = make_le("0.0, -0.1",  "HOMO/LUMO per layer (eV)", 240)
        self.sigma  = make_le("0.07, 0.07", "Gaussian disorder sigma per layer (eV)", 240)
        self.E_field = make_dspin(1e-3, 0.0, 1.0, 5, 1e-4, "Applied electric field (V/Ang)")
        self.zero_f  = QCheckBox("Also compute zero-field mobility  (Einstein relation)")
        self.zero_f.setChecked(True)
        add_row(g1, "Mean site energy",  self.mean_E,  "eV / layer")
        add_row(g1, "Disorder \u03c3",      self.sigma,   "eV / layer")
        add_row(g1, "Electric field",    self.E_field, "V/\u00c5")
        rz = g1.rowCount(); g1.addWidget(self.zero_f, rz, 1, 1, 2)
        vbox.addWidget(box1)

        # ── Marcus / MLJ / Small-polaron parameters ────────────────────
        self.box_marcus, gm = make_group("Marcus / MLJ / Small-Polaron Parameters  (from ORCA)")
        self.lam_r   = make_le("0.2, 0.25",  "Reorganisation energy lambda per layer (eV)", 240)
        self.J_intra = make_le("0.05, 0.05", "Intra-layer transfer integral J (eV)", 240)
        self.J_inter = make_le("0.03",        "Inter-layer transfer integral J (eV)", 240)
        add_row(gm, "Reorganisation \u03bb", self.lam_r,   "eV / layer")
        add_row(gm, "J intra-layer",         self.J_intra, "eV / layer")
        add_row(gm, "J inter-layer",         self.J_inter, "eV")
        vbox.addWidget(self.box_marcus)

        # ── MLJ extra ─────────────────────────────────────────────────
        self.box_mlj, gmlj = make_group("MLJ / Small-Polaron Extra Parameters  (from ORCA freq calc)")
        self.lam_outer  = make_le("0.15, 0.18", "Outer-sphere lambda_out per layer (eV)", 240)
        self.huang_rhys = make_dspin(0.5, 0.0, 20.0, 2, 0.1, "Huang-Rhys factor S = lambda_in / hbar_omega")
        self.hbar_omega = make_dspin(0.18, 0.01, 1.0, 3, 0.01, "Dominant vibrational mode hbar_omega (eV)")
        self.mlj_nmax   = make_spin(20, 5, 100, "Max vibrational levels N in FC sum")
        add_row(gmlj, "lambda_outer",       self.lam_outer,  "eV / layer")
        add_row(gmlj, "Huang-Rhys S",       self.huang_rhys)
        add_row(gmlj, "hbar\u03c9",         self.hbar_omega, "eV")
        add_row(gmlj, "FC levels N",        self.mlj_nmax)
        vbox.addWidget(self.box_mlj)

        # ── Miller-Abrahams parameters ─────────────────────────────────
        self.box_ma, gma = make_group("Miller-Abrahams Parameters  (empirical, no ORCA needed)")
        self.ma_nu0   = make_dspin(1e12, 1e8, 1e15, 2, 1e11, "Attempt frequency nu_0 (s^-1)")
        self.ma_nu0.setDecimals(2)
        self.ma_alpha = make_dspin(2.0, 0.1, 10.0, 2, 0.1, "Inverse localization length alpha (Ang^-1)")
        add_row(gma, "nu_0 (attempt freq)", self.ma_nu0,   "s\u207b\u00b9")
        add_row(gma, "alpha (loc. length)", self.ma_alpha, "\u00c5\u207b\u00b9")
        vbox.addWidget(self.box_ma)

        # ── Forster parameters ─────────────────────────────────────────
        self.box_forster, gf = make_group("Forster FRET Parameters  (from TDDFT)")
        self.forster_R0    = make_dspin(3.0, 0.5, 20.0, 2, 0.1, "Forster radius R_0 (nm)")
        self.forster_kD    = make_dspin(1e8, 1e5, 1e12, 2, 1e7, "Donor singlet decay rate k_D (s^-1)")
        self.forster_kD.setDecimals(2)
        self.forster_kappa = make_dspin(0.6667, 0.0, 4.0, 4, 0.1, "Dipole orientation factor kappa^2 (2/3 isotropic)")
        add_row(gf, "Forster radius R_0",   self.forster_R0,    "nm")
        add_row(gf, "Donor decay k_D",      self.forster_kD,    "s\u207b\u00b9")
        add_row(gf, "kappa\u00b2 (orient.)", self.forster_kappa)
        vbox.addWidget(self.box_forster)

        # ── GDM / CDM parameters ───────────────────────────────────────
        self.box_gdm, gg = make_group("GDM / CDM Parameters  (empirical, no ORCA needed)")
        self.gdm_nu0  = make_dspin(1e12, 1e8, 1e15, 2, 1e11, "Attempt frequency nu_0 (s^-1)")
        self.gdm_nu0.setDecimals(2)
        self.gdm_rloc = make_dspin(1.0, 0.1, 20.0, 2, 0.1, "Localization radius r_loc (Ang) — CDM only")
        add_row(gg, "nu_0 (attempt freq)", self.gdm_nu0,  "s\u207b\u00b9")
        add_row(gg, "r_loc (CDM only)",    self.gdm_rloc, "\u00c5")
        vbox.addWidget(self.box_gdm)

        vbox.addStretch()
        lay = QVBoxLayout(self); lay.addWidget(scrolled(inner))

        self._update_model_info(0)

    # ------------------------------------------------------------------
    def _current_model_key(self):
        idx = self.rate_model.currentIndex()
        return _MODELS[idx] if idx < len(_MODELS) else 'marcus'

    def _update_model_info(self, _=None):
        key = self._current_model_key()
        self.model_info.setText(_MODEL_INFO.get(key, ''))

        needs_marcus  = key in ('marcus', 'mlj', 'small_polaron')
        needs_mlj     = key in ('mlj', 'small_polaron')
        needs_ma      = key == 'miller_abrahams'
        needs_forster = key == 'forster'
        needs_gdm     = key in ('gdm', 'cdm')

        self.box_marcus.setVisible(needs_marcus)
        self.box_mlj.setVisible(needs_mlj)
        self.box_ma.setVisible(needs_ma)
        self.box_forster.setVisible(needs_forster)
        self.box_gdm.setVisible(needs_gdm)

    # ------------------------------------------------------------------
    def get_values(self):
        key = self._current_model_key()
        d = {
            'charge_rate_model': key,
            'carrier_type'     : self.carrier.currentText(),
            'mean_site_energy' : parse_list(self.mean_E.text()),
            'disorder_sigma'   : parse_list(self.sigma.text()),
            'electric_field'   : self.E_field.value(),
            'zero_field_runs'  : self.zero_f.isChecked(),
            # Marcus/MLJ/polaron
            'lambda_reorg'     : parse_list(self.lam_r.text()),
            'J_intra'          : parse_list(self.J_intra.text()),
            'J_inter'          : parse_list(self.J_inter.text()),
            # MLJ/polaron
            'lambda_outer'     : parse_list(self.lam_outer.text()),
            'huang_rhys_factor': self.huang_rhys.value(),
            'hbar_omega'       : self.hbar_omega.value(),
            'mlj_n_max'        : self.mlj_nmax.value(),
            # Miller-Abrahams
            'ma_nu0'           : self.ma_nu0.value(),
            'ma_alpha'         : self.ma_alpha.value(),
            # Forster
            'forster_radius'   : self.forster_R0.value(),
            'k_singlet_decay'  : self.forster_kD.value(),
            'forster_kappa2'   : self.forster_kappa.value(),
            # GDM/CDM
            'gdm_nu0'          : self.gdm_nu0.value(),
            'gdm_r_loc'        : self.gdm_rloc.value(),
        }
        return d

    def set_values(self, d):
        if 'charge_rate_model' in d:
            key = str(d['charge_rate_model'])
            if key in _MODELS:
                self.rate_model.setCurrentIndex(_MODELS.index(key))
        if 'carrier_type'      in d: self.carrier.setCurrentText(d['carrier_type'])
        if 'mean_site_energy'  in d: self.mean_E.setText(l2s(d['mean_site_energy']))
        if 'disorder_sigma'    in d: self.sigma.setText(l2s(d['disorder_sigma']))
        if 'electric_field'    in d: self.E_field.setValue(float(d['electric_field']))
        if 'zero_field_runs'   in d: self.zero_f.setChecked(bool(d['zero_field_runs']))
        if 'lambda_reorg'      in d: self.lam_r.setText(l2s(d['lambda_reorg']))
        if 'J_intra'           in d: self.J_intra.setText(l2s(d['J_intra']))
        if 'J_inter'           in d: self.J_inter.setText(l2s(d['J_inter']))
        if 'lambda_outer'      in d: self.lam_outer.setText(l2s(d['lambda_outer']))
        if 'huang_rhys_factor' in d: self.huang_rhys.setValue(float(d['huang_rhys_factor']))
        if 'hbar_omega'        in d: self.hbar_omega.setValue(float(d['hbar_omega']))
        if 'mlj_n_max'         in d: self.mlj_nmax.setValue(int(d['mlj_n_max']))
        if 'ma_nu0'            in d: self.ma_nu0.setValue(float(d['ma_nu0']))
        if 'ma_alpha'          in d: self.ma_alpha.setValue(float(d['ma_alpha']))
        if 'forster_radius'    in d: self.forster_R0.setValue(float(d['forster_radius']))
        if 'k_singlet_decay'   in d: self.forster_kD.setValue(float(d['k_singlet_decay']))
        if 'forster_kappa2'    in d: self.forster_kappa.setValue(float(d['forster_kappa2']))
        if 'gdm_nu0'           in d: self.gdm_nu0.setValue(float(d['gdm_nu0']))
        if 'gdm_r_loc'         in d: self.gdm_rloc.setValue(float(d['gdm_r_loc']))
