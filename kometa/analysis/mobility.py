# =============================================================================
# KOMETA — kometa/analysis/mobility.py
# =============================================================================
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

from ..core.constants import KB, HBAR, ANG_TO_CM, E_CHARGE
from ..core.rates import calculate_marcus_rates, calculate_mlj_rates, get_charge_rate_matrix
from ..core.bkl import run_kmc

def calculate_zero_field_mobility(system, params):
    """
    Compute zero-field mobility by running KMC at F=0 and using Einstein.
    Returns zero-field mobility and diffusion coefficient.
    """
    print("\n" + "="*60)
    print("  MODULE 5B: Zero-field mobility (Einstein relation)...")
    print("="*60)

    # Temporarily set F = 0
    params_zf = dict(params)
    params_zf['electric_field'] = 0.0

    if params.get('charge_rate_model', 'marcus') == 'mlj':
        rate_matrix_zf = calculate_mlj_rates(system, params_zf)
    else:
        rate_matrix_zf = get_charge_rate_matrix(system, params_zf)

    k_B = 8.617333e-5
    T   = params['temperature']
    q   = 1.0  # elementary charge (hole)

    D_list = []
    n_runs = min(params['num_runs'], 5)  # fewer runs needed

    for run_idx in range(n_runs):
        result = run_kmc(system, rate_matrix_zf, params_zf,
                         run_index=run_idx+500, zero_field=True)
        t_tot = result['total_time']
        if t_tot > 0:
            msd = np.mean(result['displacement']**2) * 1e-16  # Ang^2 -> cm^2
            D = msd / (2 * t_tot)
            D_list.append(D)

    D_mean  = np.mean(D_list)
    mu0     = q * D_mean / (k_B * T)

    print(f"  [OK] D (zero-field)        : {D_mean:.4e} cm^2/s")
    print(f"  [OK] mu_0 (Einstein)       : {mu0:.4e} cm^2/(V.s)")

    return {'D_zero_field': D_mean, 'mu_zero_field': mu0}


# =============================================================================
# MODULE 5C: DOMAIN-RESOLVED MOBILITY ANALYSIS
# =============================================================================
# Separates carrier trajectory into crystalline-domain hops and
# grain-boundary hops. Reports mobility and residence time per domain type.
# This directly answers: is transport filamentary? how much do GBs slow it?
# =============================================================================

def analyze_domain_mobility(all_results, system, params):
    """
    Splits carrier trajectory into crystalline vs grain-boundary segments.
    Reports mobility and residence fraction per domain type.
    """
    print("\n" + "="*60)
    print("  MODULE 5C: Domain-resolved mobility analysis...")
    print("="*60)

    is_gb     = system['is_grain_boundary']
    domain_id = system['domain_id']
    positions = system['positions']

    k_B   = 8.617333e-5
    T     = params['temperature']
    F     = params['electric_field'] * 1e8   # V/Ang -> V/cm

    gb_time    = []
    cryst_time = []
    gb_disp    = []
    cryst_disp = []

    for result in all_results:
        sites      = result['site']
        times      = result['time']
        disps      = result['displacement']

        for step in range(len(sites)-1):
            dt   = times[step+1] - times[step]
            dz   = disps[step+1] - disps[step]
            site = sites[step]
            if is_gb[site]:
                gb_time.append(dt)
                gb_disp.append(dz)
            else:
                cryst_time.append(dt)
                cryst_disp.append(dz)

    total_t = sum(gb_time) + sum(cryst_time) + 1e-100

    f_gb    = sum(gb_time)    / total_t
    f_cryst = sum(cryst_time) / total_t

    # Mobility in each region (drift velocity / field)
    def region_mobility(disps_list, times_list):
        if not times_list or F == 0:
            return 0.0
        v_drift = (sum(disps_list) * 1e-8) / sum(times_list)  # cm/s
        return abs(v_drift / F)

    mu_gb    = region_mobility(gb_disp,    gb_time)
    mu_cryst = region_mobility(cryst_disp, cryst_time)

    n_gb_sites    = int(np.sum(is_gb))
    n_cryst_sites = int(np.sum(~is_gb))
    n_domains     = len(set(domain_id[~is_gb]))

    print(f"  [OK] Crystalline sites     : {n_cryst_sites} ({f_cryst*100:.1f}% of trajectory time)")
    print(f"  [OK] Grain boundary sites  : {n_gb_sites} ({f_gb*100:.1f}% of trajectory time)")
    print(f"  [OK] Number of domains     : {n_domains}")
    print(f"  [OK] mu (crystalline)      : {mu_cryst:.4e} cm^2/(V.s)")
    print(f"  [OK] mu (grain boundary)   : {mu_gb:.4e} cm^2/(V.s)")
    if mu_cryst > 0:
        ratio = mu_cryst / max(mu_gb, 1e-30)
        print(f"  [OK] Mobility ratio cryst/GB: {ratio:.1f}x")

    return {
        'f_crystalline'   : f_cryst,
        'f_grain_boundary': f_gb,
        'mu_crystalline'  : mu_cryst,
        'mu_grain_boundary': mu_gb,
        'n_domains'       : n_domains,
        'n_gb_sites'      : n_gb_sites,
        'n_cryst_sites'   : n_cryst_sites,
    }


# =============================================================================
# MODULE 6: KMC CORE -- TTA SIMULATION (two-body BKL)
# =============================================================================



def calculate_transport_properties(all_results, params, system):
    print("\n" + "="*60)
    print("  MODULE 8: Calculating charge transport properties...")
    print("="*60)

    k_B   = 8.617333e-5
    T     = params['temperature']
    F     = params['electric_field']
    F_Vcm = F * 1e8

    mobilities       = []
    diffusion_coeffs = []
    total_crossings  = np.zeros((system['num_layers'], system['num_layers']))
    layer_fractions  = np.zeros(system['num_layers'])

    for run_result in all_results:
        time_arr  = run_result['time']
        disp_arr  = run_result['displacement']
        layer_arr = run_result['layer']
        total_t   = run_result['total_time']

        if total_t <= 0:
            continue

        final_disp_cm = disp_arr[-1] * 1e-8
        mobility = abs(final_disp_cm / (F_Vcm * total_t))
        mobilities.append(mobility)

        msd = np.mean(disp_arr**2) * 1e-16
        D   = msd / (2 * total_t)
        diffusion_coeffs.append(D)

        for layer_idx in range(system['num_layers']):
            layer_fractions[layer_idx] += (
                np.sum(layer_arr == layer_idx) / len(layer_arr))

        total_crossings += run_result['interlayer_crossings']

    num_runs = max(len(all_results), 1)
    layer_fractions /= num_runs

    mobility_mean = np.mean(mobilities)
    mobility_std  = np.std(mobilities)
    D_mean        = np.mean(diffusion_coeffs)

    # Einstein relation check
    mu_einstein = 1.0 * D_mean / (k_B * T)
    ratio = mobility_mean / mu_einstein if mu_einstein > 0 else 0

    print(f"\n  Mobility mu  = {mobility_mean:.4e} +/- {mobility_std:.2e} cm^2/(V.s)")
    print(f"  Diffusion D  = {D_mean:.4e} cm^2/s")
    print(f"  mu(Einstein) = {mu_einstein:.4e} cm^2/(V.s)  [ratio: {ratio:.2f}]")
    for i, frac in enumerate(layer_fractions):
        print(f"  Layer {i+1} residence: {frac*100:.1f}%")

    return {
        'mobility_mean'  : mobility_mean,
        'mobility_std'   : mobility_std,
        'mobility_all'   : mobilities,
        'diffusion_mean' : D_mean,
        'mu_einstein'    : mu_einstein,
        'einstein_ratio' : ratio,
        'layer_fractions': layer_fractions,
        'total_crossings': total_crossings,
    }


# =============================================================================
# MODULE 9: PLOTTING (all modes)
# =============================================================================

