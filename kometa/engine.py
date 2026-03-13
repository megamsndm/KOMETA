# =============================================================================
# KOMETA — kometa/engine.py
# Central simulation orchestrator.
# Replaces the monolithic main() from engine.py.
# Called by the GUI (RunPanel) and by the CLI entry point.
# =============================================================================

import time
from typing import Optional
import numpy as np

from .core.morphology  import build_system
from .core.rates       import (calculate_marcus_rates, calculate_mlj_rates,
                               calculate_dexter_rates, get_charge_rate_matrix)
from .core.bkl         import run_kmc
from .core.quenching   import assign_quencher_sites
from .core.tta         import run_TTA_simulation, analyze_TTA
from .core.diffusion   import calculate_triplet_diffusion
from .core.advanced    import (
    apply_trap_dos, build_voronoi_neighbors,
    run_bipolar_simulation, run_multicarrier_kmc,
    run_tof_simulation, run_temperature_scan,
    run_singlet_kmc, run_tpq_simulation,
    compute_energy_correlation, compute_tresp_couplings,
    compute_diffusion_tensor, coarse_grain_rates,
    identify_traps, save_hdf5_state,
    compute_seebeck, run_intensity_scan,
)
from .analysis.mobility import (
    calculate_zero_field_mobility, analyze_domain_mobility,
    calculate_transport_properties,
)
from .analysis.plots   import plot_results, plot_extended_results
from .io.params        import get_input_parameters
from .io.writer        import save_results


# =============================================================================
# PUBLIC API
# =============================================================================

def run(params: Optional[dict] = None) -> dict:
    """
    Run a full KOMETA simulation.

    Parameters
    ----------
    params : dict, optional
        If None, loads defaults from io/params.py.
        Pass a dict (e.g. from the GUI) to override any parameter.

    Returns
    -------
    results : dict
        All computed quantities — transport, tta_stats, diffusion_stats, etc.
        Safe to call from GUI thread via SimWorker.
    """
    print("\n")
    print("+================================================================+")
    print("|  KOMETA — KMC Organic Molecular Exciton Transport         |")
    print("+================================================================+")

    t0 = time.time()
    p  = get_input_parameters()   # load defaults
    if params:
        p.update(params)          # GUI / caller overrides

    # ── Output folder = molecule name if not set ─────────────────────
    mol = p.get('molecule_name', 'molecule') or 'molecule'
    if not p.get('output_folder', '').strip():
        p['output_folder'] = mol

    # ── Derive which calculations to run from checkbox flags ─────────
    _mode = p.get('simulation_mode', 'charge')
    do_charge    = p.get('calc_charge',    _mode in ('charge', 'both', 'full'))
    do_diffusion = p.get('calc_diffusion', _mode in ('diffusion', 'both', 'full'))
    do_tta       = p.get('calc_tta',       _mode in ('TTA', 'both', 'full'))

    calcs = [c2 for c2, f in [('Charge', do_charge), ('Diffusion', do_diffusion), ('TTA', do_tta)] if f]
    print(f"\n  Calculations  : {' | '.join(calcs) if calcs else 'none'}")
    print(f"  Rate model    : {p['charge_rate_model'].upper()}")
    print(f"  Molecule      : {mol}")
    print(f"  Morphology    : {p['morphology_source']}")
    print(f"  Temperature   : {p['temperature']} K")
    print(f"  Layers        : {p['num_layers']}  |  Sites/layer: {p['sites_per_layer']}")
    print(f"  Output folder : {p['output_folder']}")

    # ── BUILD MORPHOLOGY ────────────────────────────────────────────────
    system = build_system(p)

    if p.get('use_voronoi_neighbors', False):
        system = build_voronoi_neighbors(system, p)
    system = apply_trap_dos(system, p)

    # ── RESULT HOLDERS ──────────────────────────────────────────────────
    all_charge_results = transport = tta_raw = tta_stats = None
    diffusion_stats = domain_stats = zero_field_stats = dexter_matrix = None
    bipolar_stats = multicarrier_stats = tof_stats = temp_scan = None
    singlet_stats = tpq_stats = energy_corr = diff_tensor = None
    seebeck_stats = intensity_scan = trap_id = tresp_stats = cg_stats = None

    is_quencher = assign_quencher_sites(system, p)

    # ── TRIPLET / DEXTER PATH ───────────────────────────────────────────
    if do_diffusion or do_tta:
        dexter_matrix   = calculate_dexter_rates(system, p)
        diffusion_stats = calculate_triplet_diffusion(
            system, dexter_matrix, p, is_quencher=is_quencher)
        diff_tensor = compute_diffusion_tensor(system, dexter_matrix, p)
        if p.get('enable_tpq', False):
            tpq_stats = run_tpq_simulation(system, dexter_matrix, p)

    # ── CHARGE TRANSPORT PATH ───────────────────────────────────────────
    if do_charge:
        print("\n  -- RUNNING CHARGE TRANSPORT --")
        rate_matrix = get_charge_rate_matrix(system, p)

        trap_id  = identify_traps(system, rate_matrix, p)
        cg_stats = coarse_grain_rates(system, rate_matrix, p)

        print("\n" + "="*60)
        print(f"  MODULE 5: Running {p['num_runs']} charge KMC runs...")
        print("="*60)
        all_charge_results = []
        for run_idx in range(p['num_runs']):
            print(f"  -> Run {run_idx+1}/{p['num_runs']} ...", end=' ', flush=True)
            t_run = time.time()
            result = run_kmc(system, rate_matrix, p, run_index=run_idx)
            all_charge_results.append(result)
            print(f"done in {time.time()-t_run:.2f}s")

        transport    = calculate_transport_properties(all_charge_results, p, system)
        domain_stats = analyze_domain_mobility(all_charge_results, system, p)
        if p.get('zero_field_runs', True):
            zero_field_stats = calculate_zero_field_mobility(system, p)
        seebeck_stats = compute_seebeck(all_charge_results, system, p)
        if p.get('enable_bipolar', False):
            bipolar_stats = run_bipolar_simulation(system, p)
        if p.get('enable_multicarrier', False):
            multicarrier_stats = run_multicarrier_kmc(system, rate_matrix, p)
        if p.get('enable_tof', False):
            tof_stats = run_tof_simulation(system, rate_matrix, p)
        if p.get('enable_temperature_scan', False):
            temp_scan = run_temperature_scan(system, p)

    # ── TTA PATH ────────────────────────────────────────────────────────
    if do_tta:
        print("\n  -- RUNNING TTA SIMULATION --")
        if dexter_matrix is None:
            dexter_matrix = calculate_dexter_rates(system, p)
        tta_raw   = run_TTA_simulation(system, dexter_matrix, p,
                                       is_quencher=is_quencher)
        tta_stats = analyze_TTA(tta_raw, p)
        intensity_scan = run_intensity_scan(system, dexter_matrix, p,
                                            is_quencher=is_quencher)

    # ── OPTIONAL MODULES ────────────────────────────────────────────────
    if p.get('enable_singlet_kmc', False):
        singlet_stats = run_singlet_kmc(system, p)

    energy_corr = compute_energy_correlation(system, p)
    tresp_stats = compute_tresp_couplings(system, p)

    # ── OUTPUT ──────────────────────────────────────────────────────────
    plot_results(all_charge_results, transport, system, p,
                 tta_raw, tta_stats, diffusion_stats, domain_stats,
                 zero_field_stats)
    plot_extended_results(p,
                          temp_scan=temp_scan, tof_stats=tof_stats,
                          energy_corr=energy_corr,
                          diffusion_tensor=diff_tensor,
                          seebeck=seebeck_stats,
                          intensity_scan=intensity_scan,
                          singlet_stats=singlet_stats,
                          tpq_stats=tpq_stats,
                          bipolar_stats=bipolar_stats)
    save_results(transport, p, system, tta_stats,
                 diffusion_stats, domain_stats, zero_field_stats)
    save_hdf5_state(system, p, transport, tta_stats, diffusion_stats)

    # ── SUMMARY ─────────────────────────────────────────────────────────
    elapsed = time.time() - t0
    _print_summary(elapsed, transport, zero_field_stats, domain_stats,
                   diffusion_stats, tta_stats, bipolar_stats,
                   tof_stats, temp_scan, seebeck_stats, diff_tensor,
                   singlet_stats, tpq_stats, energy_corr, p)

    return {
        'transport'       : transport,
        'zero_field_stats': zero_field_stats,
        'domain_stats'    : domain_stats,
        'diffusion_stats' : diffusion_stats,
        'tta_stats'       : tta_stats,
        'bipolar_stats'   : bipolar_stats,
        'tof_stats'       : tof_stats,
        'temp_scan'       : temp_scan,
        'seebeck_stats'   : seebeck_stats,
        'diffusion_tensor': diff_tensor,
        'singlet_stats'   : singlet_stats,
        'tpq_stats'       : tpq_stats,
        'energy_corr'     : energy_corr,
        'system'          : system,
        'params'          : p,
        'elapsed_s'       : elapsed,
    }


# =============================================================================
# INTERNAL
# =============================================================================

def _print_summary(elapsed, transport, zero_field_stats, domain_stats,
                   diffusion_stats, tta_stats, bipolar_stats,
                   tof_stats, temp_scan, seebeck_stats, diff_tensor,
                   singlet_stats, tpq_stats, energy_corr, p):
    print("\n+================================================================+")
    print(f"|  COMPLETE  |  Runtime: {elapsed:.1f}s")
    print("+================================================================+")
    if transport:
        print(f"|  mu (drift)       = {transport['mobility_mean']:.4e} cm²/(V·s)")
    if zero_field_stats:
        print(f"|  mu_0 (Einstein)  = {zero_field_stats['mu_zero_field']:.4e} cm²/(V·s)")
    if domain_stats:
        print(f"|  mu (crystalline) = {domain_stats['mu_crystalline']:.4e} cm²/(V·s)")
        print(f"|  mu (grain bndry) = {domain_stats['mu_grain_boundary']:.4e} cm²/(V·s)")
    if diffusion_stats:
        print(f"|  D_T (triplet)    = {diffusion_stats['D_T']:.4e} cm²/s")
        print(f"|  gamma_TTA        = {diffusion_stats['gamma_TTA']:.4e} cm³/s")
    if tta_stats:
        print(f"|  P(TTA)           = {tta_stats['P_TTA']:.4f}")
        print(f"|  Phi_UC           = {tta_stats['phi_UC']:.6f}")
    if bipolar_stats:
        print(f"|  mu_hole          = {bipolar_stats['mu_hole']:.4e} cm²/(V·s)")
        print(f"|  mu_electron      = {bipolar_stats['mu_electron']:.4e} cm²/(V·s)")
    if tof_stats and tof_stats.get('mu_tof', 0) > 0:
        print(f"|  mu_TOF           = {tof_stats['mu_tof']:.4e} cm²/(V·s)")
    if temp_scan:
        print(f"|  mu(T) regime     = {temp_scan['regime']}")
    if seebeck_stats:
        print(f"|  Seebeck S        = {seebeck_stats['S_uV_K']:.1f} µV/K")
    if diff_tensor:
        print(f"|  D_T (isotropic)  = {diff_tensor['D_iso']:.4e} cm²/s")
    if singlet_stats:
        print(f"|  L_S (singlet)    = {singlet_stats['L_S_Ang']:.1f} Å")
    if tpq_stats:
        print(f"|  TPQ fraction     = {tpq_stats['tpq_fraction']*100:.1f}%")
    if energy_corr:
        print(f"|  Corr. length ξ   = {energy_corr['xi']:.2f} Å")
    print(f"|  Results → '{p['output_folder']}/'")
    print("+================================================================+\n")
