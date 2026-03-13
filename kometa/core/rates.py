# =============================================================================
# KOMETA — kometa/core/rates.py
# =============================================================================
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

from .constants import KB, HBAR, ANG_TO_CM, E_CHARGE

from math import factorial

def calculate_marcus_rates(system, params):
    print("\n" + "="*60)
    print("  MODULE 3: Calculating Marcus charge-hopping rates...")
    print("="*60)

    k_B  = 8.617333e-5
    hbar = 6.582119e-16
    T    = params['temperature']
    kT   = k_B * T
    F    = params['electric_field']
    q_sign = +1 if params['carrier_type'] == 'hole' else -1

    total_sites   = system['total_sites']
    site_energies = system['site_energies']
    positions     = system['positions']
    neighbors     = system['neighbors']
    site_layer    = system['site_layer']
    gb_lambda     = system['gb_lambda']
    is_gb         = system['is_grain_boundary']

    # Load J map if available
    J_map = _load_J_map(params, total_sites)

    rate_matrix = [[] for _ in range(total_sites)]
    count = 0

    for i in range(total_sites):
        layer_i = site_layer[i]
        for (hop_type, j, layer_j) in neighbors[i]:

            # Transfer integral
            if J_map is not None and (i, j) in J_map:
                J = J_map[(i, j)]
            else:
                if hop_type == 'intra':
                    J = params['J_intra'][min(layer_i, len(params['J_intra'])-1)]
                else:
                    J = params['J_inter'][min(min(layer_i, layer_j),
                                            len(params['J_inter'])-1)]
                # Reduce J at grain boundaries
                if is_gb[i] or is_gb[j]:
                    J *= params.get('gb_J_factor', 0.3)

            # Reorganization energy (use per-site GB-modified lambda)
            lam = 0.5 * (gb_lambda[i] + gb_lambda[j])

            delta_G = ((site_energies[j] - site_energies[i])
                       - q_sign * F * (positions[j, 2] - positions[i, 2]))

            if lam <= 0 or J <= 0:
                continue

            prefactor = (2*np.pi/hbar) * J**2 / np.sqrt(4*np.pi*lam*kT)
            exponent  = -((delta_G + lam)**2) / (4*lam*kT)
            rate_matrix[i].append((prefactor * np.exp(exponent), j))
            count += 1

    print(f"  [OK] Marcus rates computed : {count}")
    print(f"  [OK] Thermal energy kT     : {kT:.4f} eV")
    return rate_matrix


def _load_J_map(params, total_sites):
    """Load transfer integral table from file if available."""
    if params.get('J_source', 'uniform') != 'orca_table':
        return None
    filepath = params.get('J_table_file', '')
    if not os.path.exists(filepath):
        print(f"  [WARNING] J table not found: {filepath}. Using uniform J.")
        return None
    J_map = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                i, j, J_val = int(parts[0]), int(parts[1]), float(parts[2])
                J_map[(i, j)] = J_val
                J_map[(j, i)] = J_val
    print(f"  [OK] J table loaded: {len(J_map)//2} pairs")
    return J_map


# =============================================================================
# MODULE 3B: MARCUS-LEVICH-JORTNER (MLJ) RATE CALCULATOR
# =============================================================================
# MLJ extends Marcus theory to include quantum nuclear tunnelling through
# high-frequency intramolecular vibrational modes (C=C stretches etc.).
# Essential at low temperatures where classical Marcus underestimates rates.
#
# MLJ RATE FORMULA:
#   k_MLJ = (2pi/hbar) * |J|^2 / sqrt(4*pi*lambda_out*kT)
#           * sum_{n=0}^{N} [ e^(-S) * S^n / n! ]
#             * exp( -(DeltaG + lambda_out + n*hbar_omega)^2 / (4*lambda_out*kT) )
#
# WHERE:
#   lambda_out : outer-sphere (classical) reorganization energy
#   S          : Huang-Rhys factor = lambda_in / hbar_omega
#   hbar_omega : dominant quantum vibrational mode frequency (eV)
#   n          : vibrational quantum number of the acceptor final state
#   N          : truncation of the sum (typically n_max = 20 is sufficient)
# =============================================================================



def calculate_mlj_rates(system, params):
    print("\n" + "="*60)
    print("  MODULE 3B: Calculating Marcus-Levich-Jortner rates...")
    print("="*60)

    k_B       = 8.617333e-5
    hbar      = 6.582119e-16
    T         = params['temperature']
    kT        = k_B * T
    F         = params['electric_field']
    q_sign    = +1 if params['carrier_type'] == 'hole' else -1
    S         = params['huang_rhys_factor']    # Huang-Rhys factor
    hbar_w    = params['hbar_omega']            # eV
    n_max     = params['mlj_n_max']             # truncation

    total_sites   = system['total_sites']
    site_energies = system['site_energies']
    positions     = system['positions']
    neighbors     = system['neighbors']
    site_layer    = system['site_layer']
    gb_lambda     = system['gb_lambda']
    is_gb         = system['is_grain_boundary']
    J_map         = _load_J_map(params, total_sites)

    # Precompute Franck-Condon factors e^(-S) * S^n / n!
    from math import factorial, exp
    fc_factors = np.array([np.exp(-S) * S**n / factorial(n) for n in range(n_max+1)])

    rate_matrix = [[] for _ in range(total_sites)]
    count = 0

    for i in range(total_sites):
        layer_i = site_layer[i]
        for (hop_type, j, layer_j) in neighbors[i]:

            if J_map and (i, j) in J_map:
                J = J_map[(i, j)]
            else:
                if hop_type == 'intra':
                    J = params['J_intra'][min(layer_i, len(params['J_intra'])-1)]
                else:
                    J = params['J_inter'][min(min(layer_i, layer_j),
                                            len(params['J_inter'])-1)]
                if is_gb[i] or is_gb[j]:
                    J *= params.get('gb_J_factor', 0.3)

            # Outer-sphere lambda (classical part only)
            lam_out_i = params['lambda_outer'][min(layer_i, len(params['lambda_outer'])-1)]
            lam_out_j = params['lambda_outer'][min(layer_j, len(params['lambda_outer'])-1)]
            lam_out = 0.5 * (lam_out_i + lam_out_j)

            if lam_out <= 0 or J <= 0:
                continue

            delta_G = ((site_energies[j] - site_energies[i])
                       - q_sign * F * (positions[j, 2] - positions[i, 2]))

            # Prefactor (uses only outer-sphere lambda)
            prefactor = (2*np.pi/hbar) * J**2 / np.sqrt(4*np.pi*lam_out*kT)

            # Sum over vibrational quantum numbers n
            rate = 0.0
            for n in range(n_max+1):
                fc  = fc_factors[n]
                arg = delta_G + lam_out + n * hbar_w
                exp_term = np.exp(-arg**2 / (4*lam_out*kT))
                rate += fc * exp_term

            rate_matrix[i].append((prefactor * rate, j))
            count += 1

    print(f"  [OK] MLJ rates computed    : {count}")
    print(f"  [OK] Huang-Rhys S          : {S:.3f}")
    print(f"  [OK] hbar_omega            : {hbar_w:.3f} eV")
    print(f"  [OK] Vibrational levels N  : {n_max}")
    return rate_matrix


# =============================================================================
# MODULE 4: DEXTER RATE CALCULATOR
# =============================================================================



def calculate_dexter_rates(system, params):
    print("\n" + "="*60)
    print("  MODULE 4: Calculating Dexter triplet hopping rates...")
    print("="*60)

    k_B = 8.617333e-5
    T   = params['temperature']
    kT  = k_B * T
    k0  = params['k_dexter_prefactor']
    L   = params['dexter_decay_length']

    total_sites      = system['total_sites']
    positions        = system['positions']
    triplet_energies = system['triplet_energies']
    neighbors        = system['neighbors']
    site_species     = system['site_species']
    is_gb            = system['is_grain_boundary']

    allow_inter = params.get('allow_inter_species_hopping', False)
    sp_k_decay  = params.get('species_k_decay', [params['k_decay']])

    dexter_matrix = [[] for _ in range(total_sites)]
    count = 0

    for i in range(total_sites):
        for (hop_type, j, layer_j) in neighbors[i]:

            # Block inter-species hopping if requested
            if not allow_inter and site_species[i] != site_species[j]:
                continue

            r_vec = positions[j] - positions[i]
            r     = np.linalg.norm(r_vec)
            if r < 0.1:
                continue

            dexter_spatial = np.exp(-2.0 * r / L)

            # Reduce rate at grain boundaries (slower diffusion through GB)
            if is_gb[i] or is_gb[j]:
                dexter_spatial *= params.get('gb_J_factor', 0.3)

            delta_E       = triplet_energies[j] - triplet_energies[i]
            thermal_factor = np.exp(-max(delta_E, 0.0) / kT)

            rate = k0 * dexter_spatial * thermal_factor
            dexter_matrix[i].append((rate, j))
            count += 1

    print(f"  [OK] Dexter rates computed : {count}")
    print(f"  [OK] k0 (prefactor)        : {k0:.2e} s^-1")
    print(f"  [OK] Decay length L        : {L:.2f} Ang")
    return dexter_matrix


# =============================================================================
# MODULE 4B: IMPURITY / QUENCHER SITE ASSIGNMENT
# =============================================================================
# Physical basis: DOI 10.1103/PhysRevApplied.14.034050
#
# Impurity molecules (e.g. non-emissive heavy-metal complexes, O2, radical
# species) act as irreversible triplet quenchers. When a diffusing triplet
# exciton encounters a quencher site it undergoes rapid Dexter-type energy
# transfer to a low-lying dark state on the impurity, with rate k_quench.
# This is PHYSICALLY DISTINCT from natural phosphorescence decay (k_decay):
#
#   k_decay  -- radiative + non-radiative decay on the HOST molecule
#   k_quench -- irreversible transfer onto IMPURITY dark state (no emission)
#
# Two BKL modes:
#   quench_on_hop = True  : triplet is quenched INSTANTLY upon hopping onto
#                           a quencher site (infinite-rate limit approximation)
#   quench_on_hop = False : k_quench enters the BKL event list as a competing
#                           rate when the triplet currently occupies a quencher
#                           site. This is the physically correct treatment.
#
# Returns: is_quencher -- bool array, shape (total_sites,)
# =============================================================================


# =============================================================================
# MODULE 3C: MILLER-ABRAHAMS RATE CALCULATOR
# =============================================================================
# The Miller-Abrahams (MA) rate is the simplest hopping model:
#   k_MA(i->j) = nu_0 * exp(-2*alpha*r_ij) * exp(-dE/kT)  if dE > 0
#              = nu_0 * exp(-2*alpha*r_ij)                  if dE <= 0
#
# No reorganization energy — transfer is purely determined by the spatial
# overlap (tunnelling factor) and energy difference (Boltzmann penalty).
# Originally derived for impurity hopping in doped semiconductors (1960).
#
# Parameters:
#   nu_0   : attempt frequency (s^-1), typically 1e12–1e13
#   alpha  : inverse localization length (Ang^-1), typically 1–3 Ang^-1
#
# Physical regime: amorphous / deeply localized organics where lambda >> J
# DOI: 10.1103/PhysRev.120.745
# =============================================================================

def calculate_miller_abrahams_rates(system, params):
    print("\n" + "="*60)
    print("  MODULE 3C: Calculating Miller-Abrahams hopping rates...")
    print("="*60)

    k_B    = 8.617333e-5
    T      = params['temperature']
    kT     = k_B * T
    F      = params['electric_field']
    q_sign = +1 if params['carrier_type'] == 'hole' else -1
    nu_0   = params.get('ma_nu0',   1e12)     # attempt frequency s^-1
    alpha  = params.get('ma_alpha', 2.0)      # inverse loc. length Ang^-1

    total_sites   = system['total_sites']
    site_energies = system['site_energies']
    positions     = system['positions']
    neighbors     = system['neighbors']
    is_gb         = system['is_grain_boundary']

    rate_matrix = [[] for _ in range(total_sites)]
    count = 0

    for i in range(total_sites):
        for (hop_type, j, layer_j) in neighbors[i]:
            r_vec = positions[j] - positions[i]
            r     = np.linalg.norm(r_vec)
            if r < 0.1:
                continue

            # GB penalty on tunnelling factor
            gb_fac = params.get('gb_J_factor', 0.3) if (is_gb[i] or is_gb[j]) else 1.0

            tunnel = np.exp(-2.0 * alpha * r) * gb_fac

            delta_G = ((site_energies[j] - site_energies[i])
                       - q_sign * F * (positions[j, 2] - positions[i, 2]))

            boltzmann = np.exp(-max(delta_G, 0.0) / kT)
            rate = nu_0 * tunnel * boltzmann
            rate_matrix[i].append((rate, j))
            count += 1

    print(f"  [OK] Miller-Abrahams rates : {count}")
    print(f"  [OK] nu_0                  : {nu_0:.2e} s^-1")
    print(f"  [OK] alpha                 : {alpha:.2f} Ang^-1")
    return rate_matrix


# =============================================================================
# MODULE 3D: FÖRSTER RATE CALCULATOR  (singlet exciton transfer)
# =============================================================================
# Förster Resonance Energy Transfer (FRET) — long-range dipole-dipole coupling
#   k_FRET(i->j) = k_D * (R_0 / r_ij)^6
#
# where:
#   k_D : donor decay rate (s^-1) = 1 / tau_D  (from TDDFT)
#   R_0 : Förster radius (Ang) — distance at which k_FRET = k_D
#       R_0^6 = (9 * kappa^2 * phi_D * J_overlap) / (128 * pi^5 * n^4)
#   kappa^2 : orientation factor (2/3 for isotropic, 4 for parallel dipoles)
#
# Physical regime: singlet exciton diffusion in fluorescent OLEDs,
# host-guest energy transfer, TADF systems.
# DOI: 10.1039/tf9600056441
# =============================================================================

def calculate_forster_rates(system, params):
    print("\n" + "="*60)
    print("  MODULE 3D: Calculating Forster singlet transfer rates...")
    print("="*60)

    k_B    = 8.617333e-5
    T      = params['temperature']
    kT     = k_B * T
    k_D    = params.get('k_singlet_decay', 1e8)   # s^-1  donor decay rate
    R0_ang = params.get('forster_radius', 3.0) * 10.0  # nm -> Ang
    r_min  = params.get('forster_r_min',  2.0)     # Ang  cutoff
    r_max  = params.get('forster_r_max',  R0_ang * 2.5)  # Ang  cutoff

    total_sites   = system['total_sites']
    positions     = system['positions']
    neighbors     = system['neighbors']
    site_energies = system['site_energies']    # here = singlet energies
    is_gb         = system['is_grain_boundary']

    rate_matrix = [[] for _ in range(total_sites)]
    count = 0

    for i in range(total_sites):
        for (hop_type, j, layer_j) in neighbors[i]:
            r_vec = positions[j] - positions[i]
            r     = np.linalg.norm(r_vec)
            if r < r_min or r > r_max:
                continue

            # Dipole orientation + GB penalty
            kappa2  = params.get('forster_kappa2', 0.6667)
            gb_fac  = params.get('gb_J_factor', 0.3) if (is_gb[i] or is_gb[j]) else 1.0

            rate = k_D * kappa2 * (R0_ang / r)**6 * gb_fac

            # Detailed balance: Boltzmann penalty for uphill transfer
            delta_E = site_energies[j] - site_energies[i]
            if delta_E > 0:
                rate *= np.exp(-delta_E / kT)

            rate_matrix[i].append((rate, j))
            count += 1

    print(f"  [OK] Forster rates computed: {count}")
    print(f"  [OK] R_0 (Forster radius)  : {R0_ang/10:.2f} nm")
    print(f"  [OK] k_D (donor decay)     : {k_D:.2e} s^-1")
    return rate_matrix


# =============================================================================
# MODULE 3E: GAUSSIAN DISORDER MODEL (GDM / Bassler model)
# =============================================================================
# The GDM describes hopping transport in amorphous organics with Gaussian
# DOS. The rate is a Boltzmann factor of site energy differences only —
# no explicit J or lambda required. The prefactor is an empirical attempt
# frequency fitted to experiment.
#
#   k_GDM(i->j) = nu_0 * exp(-(dE)^2 / (2*sigma^2)) if dE > 0
#               = nu_0                                 if dE <= 0
#
# In the CORRELATED DISORDER MODEL (CDM) variant, the exponential spatial
# decay adds a distance dependence:
#   k_CDM(i->j) = nu_0 * exp(-r/r_loc) * f_Boltzmann
#
# Parameters (no ORCA input needed — fully empirical):
#   nu_0  : hopping attempt frequency (s^-1)
#   sigma : energetic disorder width (eV) — same as disorder_sigma
#   r_loc : localization radius (Ang)     — CDM only
#
# DOI (GDM):  10.1016/0379-6779(93)90in
# DOI (CDM):  10.1103/PhysRevLett.70.2in
# =============================================================================

def calculate_gdm_rates(system, params):
    model = params.get('charge_rate_model', 'gdm')  # 'gdm' or 'cdm'
    label = 'CDM (Correlated Disorder)' if model == 'cdm' else 'GDM (Gaussian Disorder)'

    print("\n" + "="*60)
    print(f"  MODULE 3E: Calculating {label} rates...")
    print("="*60)

    k_B    = 8.617333e-5
    T      = params['temperature']
    kT     = k_B * T
    F      = params['electric_field']
    q_sign = +1 if params['carrier_type'] == 'hole' else -1
    nu_0   = params.get('gdm_nu0',  1e12)    # attempt frequency s^-1
    r_loc  = params.get('gdm_r_loc', 1.0)   # localization radius Ang (CDM)

    total_sites   = system['total_sites']
    site_energies = system['site_energies']
    positions     = system['positions']
    neighbors     = system['neighbors']
    site_layer    = system['site_layer']
    is_gb         = system['is_grain_boundary']

    rate_matrix = [[] for _ in range(total_sites)]
    count = 0

    for i in range(total_sites):
        layer_i = site_layer[i]
        sigma_i = params['disorder_sigma'][min(layer_i, len(params['disorder_sigma'])-1)]

        for (hop_type, j, layer_j) in neighbors[i]:
            r_vec = positions[j] - positions[i]
            r     = np.linalg.norm(r_vec)
            if r < 0.1:
                continue

            gb_fac = params.get('gb_J_factor', 0.3) if (is_gb[i] or is_gb[j]) else 1.0

            delta_G = ((site_energies[j] - site_energies[i])
                       - q_sign * F * (positions[j, 2] - positions[i, 2]))

            # GDM: Boltzmann factor only
            if delta_G > 0:
                boltzmann = np.exp(-delta_G / kT)
            else:
                boltzmann = 1.0

            # CDM: add spatial decay
            spatial = np.exp(-r / r_loc) if model == 'cdm' else 1.0

            rate = nu_0 * spatial * boltzmann * gb_fac
            rate_matrix[i].append((rate, j))
            count += 1

    print(f"  [OK] {label} rates : {count}")
    print(f"  [OK] nu_0          : {nu_0:.2e} s^-1")
    if model == 'cdm':
        print(f"  [OK] r_loc (CDM)   : {r_loc:.2f} Ang")
    return rate_matrix


# =============================================================================
# MODULE 3F: SMALL POLARON / HOLSTEIN RATE CALCULATOR
# =============================================================================
# The small polaron model accounts for carrier self-trapping: the carrier
# distorts the local molecular geometry, lowering its own site energy by
# the polaron binding energy E_b = lambda_in / 2.
#
# The renormalised Marcus rate is:
#   k_polaron(i->j) = (2pi/hbar) * J_eff^2 / sqrt(4pi*lambda_eff*kT)
#                     * exp(-(dG + lambda_eff)^2 / (4*lambda_eff*kT))
#
# where the effective parameters are renormalised by polaron formation:
#   J_eff        = J * exp(-S/2)          (polaron narrowing)
#   lambda_eff   = lambda_out             (only outer-sphere remains classical)
#   site_E_eff   = site_E - E_b           (self-trapping lowers site energy)
#   E_b          = S * hbar_omega         (polaron binding energy)
#
# This is the adiabatic small polaron limit. For non-adiabatic (J << kT),
# reduce to standard MLJ.
#
# Parameters: same as MLJ (S, hbar_omega, lambda_outer) — no new ORCA inputs
# DOI: 10.1103/PhysRevLett.70.2928
# =============================================================================

def calculate_small_polaron_rates(system, params):
    print("\n" + "="*60)
    print("  MODULE 3F: Calculating small polaron (Holstein) rates...")
    print("="*60)

    from math import factorial

    k_B    = 8.617333e-5
    hbar   = 6.582119e-16
    T      = params['temperature']
    kT     = k_B * T
    F      = params['electric_field']
    q_sign = +1 if params['carrier_type'] == 'hole' else -1
    S      = params['huang_rhys_factor']
    hbar_w = params['hbar_omega']
    n_max  = params.get('mlj_n_max', 20)

    # Polaron binding energy lowers site energies self-consistently
    E_b = S * hbar_w   # eV

    total_sites   = system['total_sites']
    site_energies = system['site_energies'] - E_b  # renormalised by polaron
    positions     = system['positions']
    neighbors     = system['neighbors']
    site_layer    = system['site_layer']
    gb_lambda     = system['gb_lambda']
    is_gb         = system['is_grain_boundary']
    J_map         = _load_J_map(params, total_sites)

    # Polaron band-narrowing: J_eff = J * exp(-S/2)
    polaron_factor = np.exp(-S / 2.0)

    # Franck-Condon factors (same as MLJ)
    fc_factors = np.array([np.exp(-S) * S**n / factorial(n) for n in range(n_max+1)])

    rate_matrix = [[] for _ in range(total_sites)]
    count = 0

    for i in range(total_sites):
        layer_i = site_layer[i]
        for (hop_type, j, layer_j) in neighbors[i]:

            if J_map and (i, j) in J_map:
                J = J_map[(i, j)]
            else:
                if hop_type == 'intra':
                    J = params['J_intra'][min(layer_i, len(params['J_intra'])-1)]
                else:
                    J = params['J_inter'][min(min(layer_i, layer_j),
                                            len(params['J_inter'])-1)]
                if is_gb[i] or is_gb[j]:
                    J *= params.get('gb_J_factor', 0.3)

            J_eff = J * polaron_factor   # band-narrowed coupling

            lam_out_i = params['lambda_outer'][min(layer_i, len(params['lambda_outer'])-1)]
            lam_out_j = params['lambda_outer'][min(layer_j, len(params['lambda_outer'])-1)]
            lam_out   = 0.5 * (lam_out_i + lam_out_j)

            if lam_out <= 0 or J_eff <= 0:
                continue

            delta_G = ((site_energies[j] - site_energies[i])
                       - q_sign * F * (positions[j, 2] - positions[i, 2]))

            prefactor = (2*np.pi/hbar) * J_eff**2 / np.sqrt(4*np.pi*lam_out*kT)

            rate = 0.0
            for n in range(n_max+1):
                fc  = fc_factors[n]
                arg = delta_G + lam_out + n * hbar_w
                rate += fc * np.exp(-arg**2 / (4*lam_out*kT))

            rate_matrix[i].append((prefactor * rate, j))
            count += 1

    print(f"  [OK] Small polaron rates   : {count}")
    print(f"  [OK] Polaron binding E_b   : {E_b:.4f} eV")
    print(f"  [OK] Band-narrowing factor : {polaron_factor:.4f}")
    print(f"  [OK] Huang-Rhys S          : {S:.3f}")
    return rate_matrix


# =============================================================================
# RATE DISPATCHER
# =============================================================================

def get_charge_rate_matrix(system, params):
    """
    Central dispatcher — call this instead of individual functions.
    Reads params['charge_rate_model'] and calls the appropriate calculator.

    Supported models:
        'marcus'        — Classical Marcus (default)
        'mlj'           — Marcus-Levich-Jortner (quantum nuclear tunnelling)
        'miller_abrahams'— Miller-Abrahams (amorphous, no lambda needed)
        'forster'       — Forster FRET (singlet exciton, dipole-dipole)
        'gdm'           — Gaussian Disorder Model (Bassler, empirical)
        'cdm'           — Correlated Disorder Model (GDM + spatial decay)
        'small_polaron' — Holstein small polaron (self-trapping, band-narrowing)
    """
    model = params.get('charge_rate_model', 'marcus').lower().strip()
    dispatch = {
        'marcus':          calculate_marcus_rates,
        'mlj':             calculate_mlj_rates,
        'miller_abrahams': calculate_miller_abrahams_rates,
        'forster':         calculate_forster_rates,
        'gdm':             calculate_gdm_rates,
        'cdm':             calculate_gdm_rates,   # same function, reads model key
        'small_polaron':   calculate_small_polaron_rates,
    }
    if model not in dispatch:
        print(f"  [WARNING] Unknown rate model '{model}' — falling back to Marcus.")
        model = 'marcus'
    return dispatch[model](system, params)
