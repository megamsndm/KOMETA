<div align="center">

# KOMETA

### KMC Organic Molecular Exciton Transport Analysis

*Kinetic Monte Carlo simulation platform for organic molecular semiconductors*

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue?logo=python&logoColor=white)](https://www.python.org)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey)]()
[![PyQt5](https://img.shields.io/badge/GUI-PyQt5-41CD52?logo=qt)](https://riverbankcomputing.com/software/pyqt/)

**Developed by Dr. Mohsen Doust Mohammadi**

[Overview](#overview) · [Features](#features) · [Installation](#installation) · [Usage](#usage) · [Rate Models](#rate-models) · [Parameters](#parameter-guide) · [Structure](#project-structure) · [Citation](#citation)

</div>

---

## Overview

KOMETA is a professional Kinetic Monte Carlo (KMC) simulation suite for modelling:

- **Charge transport** — drift mobility, diffusion, Einstein relation
- **Triplet exciton diffusion** — D_T, diffusion length, bimolecular annihilation rate
- **Triplet-Triplet Annihilation Upconversion (TTA-UC)** — Φ_UC, P(TTA), intensity dependence

It integrates quantum-chemistry parameters computed by **ORCA** (λ, J, T₁, S, ħω) directly into a full KMC engine, enabling first-principles-informed simulation of organic semiconductor thin films and OLED/upconversion devices.

**No molecular geometry input required.** KOMETA takes scalar parameters from ORCA and builds a model lattice internally.

```
Molecule → ORCA → λ, J, T₁, S, ħω → KOMETA → KMC simulation → Results
```

---

## Features

### 7 Rate Models
| Model | Key | Physics | Needs ORCA |
|:---|:---|:---|:---:|
| Marcus (classical) | `marcus` | Gold standard for organic semiconductors | ✔ |
| Marcus-Levich-Jortner | `mlj` | Quantum nuclear tunnelling, Franck-Condon sum | ✔ |
| Miller-Abrahams | `miller_abrahams` | Spatial tunnelling + Boltzmann, no λ needed | ✘ |
| Förster FRET | `forster` | Dipole-dipole singlet exciton transfer, r⁻⁶ | TDDFT |
| Gaussian Disorder Model | `gdm` | Bässler empirical model, purely energetic | ✘ |
| Correlated Disorder Model | `cdm` | GDM with spatial exponential decay | ✘ |
| Small Polaron (Holstein) | `small_polaron` | Band-narrowed J, carrier self-trapping | ✔ |

### 9 KMC Algorithms
| Algorithm | Complexity | Best for |
|:---|:---:|:---|
| BKL (Bortz-Kalos-Lebowitz) | O(N) | Default — all systems |
| Gillespie | O(N) | Chemical master equation basis |
| FRM (First Reaction Method) | O(log N) | Multi-particle simulations |
| NRM (Next Reaction Method) | O(log N) | TTA / bipolar transport |
| Random Selection | O(1) | Dense multi-particle |
| Rejection KMC | O(1) | Uniform rate landscapes |
| Binary Tree (Fenwick) | O(log N) | Large systems (N > 10⁴) |
| AKMC (Accelerated KMC) | Variable | Deep traps, low temperature |
| **PEGASUS** *(original)* | O(1–N) | Ordered/field-biased systems |

### Physical Modules
- Grain boundary morphology (S₂ order parameter, domain-resolved mobility)
- Multi-component molecular blends (binary and ternary)
- Impurity quenching — Dexter-type irreversible triplet quenching (Module 4B)
- Zero-field mobility — Einstein relation μ₀ = qD/kT
- Seebeck coefficient (thermoelectric characterization)
- Bipolar transport (hole + electron simultaneously)
- Time-of-flight (TOF) simulation
- Temperature scan μ(T)
- Excitation intensity scan Φ_UC(I)
- Anisotropic 3D diffusion tensor

### GUI — Deep Space Science Theme
- PyQt5 desktop application with dark navy/teal/gold palette
- 7 parameter tabs: General, Morphology, Charge, MLJ, TTA, Blend, Quenching
- Tick exactly which **calculations** and **output figures** you want — no wasted runs
- Live console output during simulation
- Output folder automatically named after your molecule
- Save/load JSON configurations

### Output
- **13 selectable figures** — choose only what you need
- **`KOMETA_Result.txt`** — dynamic summary of all computed quantities, adapts to what was run
- Output folder named after molecule automatically

---

## Installation

### Requirements
```
Python  >= 3.10
numpy   >= 1.20
scipy   >= 1.6
matplotlib >= 3.3
PyQt5   >= 5.12
h5py    >= 3.0
scikit-learn >= 0.24
```

### From source (recommended)
```bash
git clone https://github.com/megamsndm/KOMETA.git
cd KOMETA
pip install -r requirements.txt
python kometa_gui.py
```

### pip install
```bash
pip install git+https://github.com/megamsndm/KOMETA.git
kometa
```

### Windows standalone executable (no Python needed)
```
1. Install Python 3.10+ and add to PATH
2. Double-click  build_exe.bat
3. Find KOMETA.exe in  dist\
4. Copy dist\KOMETA.exe anywhere — runs on any Windows PC
```
Requires PyInstaller 6.x (auto-installed by the .bat file).

---

## Usage

### Launch GUI
```bash
python kometa_gui.py
```

### Step-by-step workflow

1. **Enter molecule name** in the green header bar — this names the output folder
2. **Run ORCA** separately to obtain: λ, J, T₁, S, ħω, Φ_em
3. **Fill in the parameter tabs:**
   - *General* — KMC steps, runs, algorithm, temperature; tick calculations and figures
   - *Morphology* — layers, sites per layer, distances, grain boundaries
   - *Charge* — rate model, λ, J, disorder σ, electric field
   - *TTA* — triplet energies, k_TTA, k_decay, Φ_emission
   - *Blend* — multi-component systems
   - *Quenching* — impurity quencher fraction
4. **Tick calculations** in the General tab: Charge Transport / Triplet Diffusion / TTA
5. **Tick output figures** — only the selected ones will be generated
6. **Press ▶ Run Simulation**
7. Results appear in a folder named after your molecule

### CLI / scripting
```python
from kometa.engine import run

results = run({
    'molecule_name'    : 'Alq3',
    'temperature'      : 300.0,
    'charge_rate_model': 'marcus',
    'calc_charge'      : True,
    'calc_tta'         : False,
    'num_runs'         : 20,
    'lambda_reorg'     : [0.21],
    'J_intra'          : [0.048],
    'disorder_sigma'   : [0.085],
    'electric_field'   : 1e-3,
})

print(f"Mobility = {results['transport']['mobility_mean']:.3e} cm²/(V·s)")
```

---

## Rate Models

### Marcus (classical)
```
k = (2π/ħ) |J|² / √(4πλkT) · exp(-(ΔG+λ)² / 4λkT)
```
Standard non-adiabatic Marcus theory. Requires λ and J from ORCA 4-point method and dimer FO calculation.

### Marcus-Levich-Jortner
Extends Marcus with a Franck-Condon sum over intramolecular quantum vibrational modes:
```
k_MLJ = prefactor · Σₙ [e⁻ˢ Sⁿ/n!] · exp(-(ΔG + λ_out + nħω)² / 4λ_out kT)
```
Essential at low temperature where C=C quantum modes dominate.

### Miller-Abrahams
```
k = ν₀ · exp(-2αr) · exp(-ΔE/kT)  if ΔE > 0
k = ν₀ · exp(-2αr)                  if ΔE ≤ 0
```
No λ needed. α = inverse localization length. Original 1960 hopping model.

### Förster FRET
```
k_FRET = k_D · κ² · (R₀/r)⁶
```
Long-range dipole-dipole singlet exciton transfer. R₀ from TDDFT spectral overlap integral.

### GDM / CDM
```
GDM:  k = ν₀ · Boltzmann(ΔE)
CDM:  k = ν₀ · exp(-r/r_loc) · Boltzmann(ΔE)
```
Fully empirical — no ORCA input needed. Bässler model for amorphous organic films.

### Small Polaron (Holstein)
MLJ with carrier self-trapping renormalization:
```
J_eff = J · exp(-S/2)      (band narrowing)
E_b   = S · ħω             (polaron binding energy, lowers site energy)
```
Applies MLJ rates with renormalized parameters. Important when λ > 0.5 eV.

---

## Parameter Guide

| Parameter | Symbol | Source | Typical range |
|:---|:---:|:---|:---:|
| Reorganization energy | λ | ORCA 4-point method | 0.10 – 0.50 eV |
| Transfer integral (intra-layer) | J_intra | ORCA dimer fragment orbital | 0.01 – 0.10 eV |
| Transfer integral (inter-layer) | J_inter | ORCA dimer fragment orbital | 0.005 – 0.05 eV |
| Triplet energy | T₁ | ORCA TDDFT / CASSCF | 1.50 – 2.50 eV |
| Huang-Rhys factor | S | ORCA frequency calculation | 0.1 – 2.0 |
| Dominant vibrational mode | ħω | ORCA frequency calculation | 0.15 – 0.20 eV |
| Energetic disorder | σ | Literature / ORCA ensemble | 0.05 – 0.12 eV |
| Outer-sphere reorganization | λ_out | ORCA (MLJ/polaron only) | 0.05 – 0.30 eV |
| Emission quantum yield | Φ_em | ORCA TDDFT | 0 – 1 |
| Förster radius | R₀ | TDDFT spectral overlap | 1 – 10 nm |
| Attempt frequency | ν₀ | Literature (MA/GDM/CDM) | 10¹¹ – 10¹³ s⁻¹ |

---

## PEGASUS Algorithm

KOMETA includes **PEGASUS** (Predictive Energy-Gradient Accelerated Stochastic Uniform Sampling) — an original KMC algorithm developed specifically for anisotropic, field-biased organic semiconductor physics.

PEGASUS identifies sequences of high-probability hops along the applied field gradient and executes them as compound cascade events — eliminating redundant random number generation without any statistical approximation.

**Expected speedup over BKL:**
- Ordered crystal (σ = 0): ~50×
- Weakly disordered film (σ = 0.05 eV): ~12×
- Typical OLED film (σ = 0.08 eV): ~7×

---

## Project Structure

```
KOMETA/
├── kometa_gui.py                  ← Entry point: python kometa_gui.py
│
├── kometa/
│   ├── engine.py                  ← Central simulation orchestrator
│   │
│   ├── core/
│   │   ├── constants.py           ← Physical constants (KB, HBAR, ...)
│   │   ├── rates.py               ← All 7 rate models + dispatcher
│   │   ├── kmc_algorithms.py      ← All 9 KMC algorithms incl. PEGASUS
│   │   ├── morphology.py          ← Grid lattice builder
│   │   ├── bkl.py                 ← BKL algorithm & run_kmc shim
│   │   ├── tta.py                 ← TTA-UC simulation
│   │   ├── diffusion.py           ← Triplet diffusion
│   │   ├── quenching.py           ← Impurity quenching (Module 4B)
│   │   └── advanced.py            ← Extended modules (bipolar, TOF, ...)
│   │
│   ├── analysis/
│   │   ├── mobility.py            ← Transport property extraction
│   │   └── plots.py               ← 13 selectable output figures
│   │
│   ├── gui/
│   │   ├── main_window.py         ← Main application window
│   │   ├── run_panel.py           ← Run / Stop / Console panel
│   │   ├── header_bar.py          ← Molecule name header
│   │   ├── results_tab.py         ← Results viewer
│   │   ├── worker.py              ← Background simulation thread
│   │   ├── style.py               ← Deep-space dark theme + Qt helpers
│   │   └── tabs/
│   │       ├── general_tab.py     ← KMC settings, calc & figure selection
│   │       ├── morphology_tab.py  ← Lattice, layers, grain boundaries
│   │       ├── charge_tab.py      ← All 7 rate models with dynamic panels
│   │       ├── mlj_tab.py         ← MLJ/polaron extra parameters
│   │       ├── tta_tab.py         ← TTA-UC parameters
│   │       ├── blend_tab.py       ← Multi-component blend
│   │       └── quench_tab.py      ← Impurity quenching
│   │
│   ├── io/
│   │   ├── params.py              ← All default parameters
│   │   ├── config.py              ← JSON config save/load
│   │   └── writer.py              ← KOMETA_Result.txt output
│   │
│   └── tests/
│       └── test_rates.py          ← Unit tests (pytest)
│
├── hooks/
│   ├── hook-kometa.py             ← PyInstaller hidden imports
│   └── rthook_kometa.py           ← PyInstaller runtime hook
│
├── KOMETA.spec                    ← PyInstaller build spec (Python 3.14 compatible)
├── build_exe.bat                  ← Windows .exe builder
├── setup.py
├── requirements.txt
├── README.md
├── LICENSE.txt
├── CITATION.txt
└── TERMS_OF_USE.txt
```

---

## Running Tests

```bash
pip install pytest
python -m pytest kometa/tests/ -v
```

---

## Citation

If you use KOMETA in your research, please cite:

```bibtex
@software{DoustMohammadi2026kometa,
  author  = {Doust Mohammadi, Mohsen},
  title   = {{KOMETA}: {KMC} Organic Molecular Exciton Transport Analysis},
  year    = {2026},
  version = {1.0},
  url     = {https://github.com/megamsndm/KOMETA}
}
```

See [CITATION.txt](CITATION.txt) for ACS and plain-text formats.

---

## License

All rights reserved.

For academic and non-commercial use: citation is required.
For commercial distribution with PyQt5: a Riverbank commercial license is required.

---

## Contact

**Dr. Mohsen Doust Mohammadi**
GitHub: [@megamsndm](https://github.com/megamsndm)
Repository: [github.com/megamsndm/KOMETA](https://github.com/megamsndm/KOMETA)
