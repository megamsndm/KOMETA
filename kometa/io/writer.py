# =============================================================================
# KOMETA — kometa/io/writer.py
# Output: summary.txt, HDF5, trajectory files
# =============================================================================
import numpy as np
import os

def save_results(transport, params, system, tta_stats=None,
                 diffusion_stats=None, domain_stats=None,
                 zero_field_stats=None):
    print("\n" + "="*60)
    print("  MODULE 10: Saving results to file...")
    print("="*60)

    os.makedirs(params['output_folder'], exist_ok=True)
    filepath = os.path.join(params['output_folder'], 'KOMETA_Result.txt')

    with open(filepath, 'w', encoding='utf-8') as f:
        # Dynamic header based on what was actually computed
        calcs_run = []
        if transport:       calcs_run.append("CHARGE TRANSPORT")
        if diffusion_stats: calcs_run.append("TRIPLET DIFFUSION")
        if tta_stats:       calcs_run.append("TTA UPCONVERSION")
        header_line = "  KOMETA RESULTS — " + " + ".join(calcs_run) if calcs_run else "  KOMETA RESULTS"
        f.write("="*65 + "\n")
        f.write(header_line + "\n")
        f.write(f"  Molecule: {params.get('molecule_name', '—')}\n")
        f.write("="*65 + "\n\n")

        f.write("INPUT PARAMETERS:\n")
        f.write("-"*40 + "\n")
        f.write(f"  Morphology source         : {params['morphology_source']}\n")
        f.write(f"  Blend type                : {params['blend_type']}\n")
        f.write(f"  Site energy model         : {params['site_energy_model']}\n")
        f.write(f"  Charge rate model         : {params['charge_rate_model']}\n")
        f.write(f"  Number of layers          : {params['num_layers']}\n")
        f.write(f"  Sites per layer           : {params['sites_per_layer']}\n")
        f.write(f"  Temperature               : {params['temperature']} K\n")
        f.write(f"  Mean site energies (eV)   : {params['mean_site_energy']}\n")
        f.write(f"  Disorder sigma (eV)       : {params['disorder_sigma']}\n")
        f.write(f"  Reorganization lambda (eV): {params['lambda_reorg']}\n")
        f.write(f"  J_intra (eV)              : {params['J_intra']}\n")
        f.write(f"  J_inter (eV)              : {params['J_inter']}\n")
        f.write(f"  Electric field (V/Ang)    : {params['electric_field']:.2e}\n")
        f.write(f"  Positional disorder (Ang) : {params['positional_disorder_sigma']}\n")
        if params['charge_rate_model'] == 'mlj':
            f.write(f"  Huang-Rhys S              : {params['huang_rhys_factor']}\n")
            f.write(f"  hbar_omega (eV)           : {params['hbar_omega']}\n")
        f.write(f"  Triplet energies T1 (eV)  : {params['triplet_energy']}\n")
        f.write(f"  Dexter k0 (s^-1)          : {params['k_dexter_prefactor']:.2e}\n")
        f.write(f"  Dexter decay length (Ang) : {params['dexter_decay_length']}\n")
        f.write(f"  k_TTA (s^-1)              : {params['k_TTA']:.2e}\n")
        f.write(f"  k_decay (s^-1)            : {params['k_decay']:.2e}\n")
        # Module 4B quenching summary
        if params.get('enable_quenching', False):
            f.write(f"\n  --- IMPURITY QUENCHING (Module 4B) ---\n")
            f.write(f"  enable_quenching          : True\n")
            f.write(f"  quencher_fraction         : {params.get('quencher_fraction',0)*100:.2f}%\n")
            f.write(f"  k_quench (s^-1)           : {params.get('k_quench',0):.2e}\n")
            f.write(f"  quench_on_hop             : {params.get('quench_on_hop',False)}\n")
            if tta_stats and 'n_quenched' in (tta_stats or {}):
                n_q   = tta_stats.get('n_quenched', 0)
                n_tot = params.get('num_triplet_pairs', 1)
                f.write(f"  Triplets quenched (TTA)   : {n_q}/{n_tot} ({100*n_q/n_tot:.1f}%)\n")
                f.write(f"  Quenching reduces Phi_UC by factor: "
                        f"{1 - n_q/max(n_tot,1):.3f}\n")
        f.write(f"  Singlet fraction          : {params['singlet_fraction']:.4f}\n")
        f.write(f"  Emission yield            : {params['phi_emission']:.4f}\n\n")

        f.write("SYSTEM PROPERTIES:\n")
        f.write("-"*40 + "\n")
        n_gb = int(np.sum(system['is_grain_boundary']))
        f.write(f"  Total sites               : {system['total_sites']}\n")
        f.write(f"  Grain boundary sites      : {n_gb} ({n_gb/system['total_sites']*100:.1f}%)\n")
        f.write(f"  Crystalline domains       : {len(set(system['domain_id'][~system['is_grain_boundary']]))}\n\n")

        if transport is not None:
            f.write("CHARGE TRANSPORT RESULTS:\n")
            f.write("-"*40 + "\n")
            f.write(f"  Mobility mu  = {transport['mobility_mean']:.6e} +/- {transport['mobility_std']:.2e} cm^2/(V.s)\n")
            f.write(f"  Diffusion D  = {transport['diffusion_mean']:.6e} cm^2/s\n")
            f.write(f"  mu(Einstein) = {transport['mu_einstein']:.6e} cm^2/(V.s)\n")
            f.write(f"  Einstein ratio = {transport['einstein_ratio']:.3f} (ideal = 1.0)\n")
            for i, frac in enumerate(transport['layer_fractions']):
                f.write(f"  Layer {i+1} residence = {frac*100:.2f}%\n")
            f.write("\n")

        if zero_field_stats is not None:
            f.write("ZERO-FIELD MOBILITY:\n")
            f.write("-"*40 + "\n")
            f.write(f"  mu_0 (Einstein) = {zero_field_stats['mu_zero_field']:.6e} cm^2/(V.s)\n")
            f.write(f"  D_0 (zero-field)= {zero_field_stats['D_zero_field']:.6e} cm^2/s\n\n")

        if domain_stats is not None:
            f.write("DOMAIN-RESOLVED MOBILITY:\n")
            f.write("-"*40 + "\n")
            f.write(f"  mu (crystalline)    = {domain_stats['mu_crystalline']:.6e} cm^2/(V.s)\n")
            f.write(f"  mu (grain boundary) = {domain_stats['mu_grain_boundary']:.6e} cm^2/(V.s)\n")
            f.write(f"  f (crystalline)     = {domain_stats['f_crystalline']*100:.1f}%\n")
            f.write(f"  f (grain boundary)  = {domain_stats['f_grain_boundary']*100:.1f}%\n")
            f.write(f"  Number of domains   = {domain_stats['n_domains']}\n\n")

        if diffusion_stats is not None:
            f.write("TRIPLET DIFFUSION:\n")
            f.write("-"*40 + "\n")
            f.write(f"  D_T              = {diffusion_stats['D_T']:.6e} cm^2/s\n")
            f.write(f"  Mean lifetime    = {diffusion_stats['mean_lifetime']:.3e} s\n")
            f.write(f"  L_diff           = {diffusion_stats['L_diff_ang']:.2f} Ang\n")
            f.write(f"  gamma_TTA (Smol) = {diffusion_stats['gamma_TTA']:.4e} cm^3/s\n\n")

        if tta_stats is not None:
            f.write("TTA UPCONVERSION RESULTS:\n")
            f.write("-"*40 + "\n")
            f.write(f"  P(TTA)           = {tta_stats['P_TTA']:.4f} +/- {tta_stats['P_TTA_err']:.4f}\n")
            f.write(f"  Phi_UC           = {tta_stats['phi_UC']:.6f}\n")
            f.write(f"  Max Phi_UC       = {tta_stats['phi_UC_max']:.6f}\n")
            f.write(f"  TTA events       = {tta_stats['n_TTA']}\n")
            f.write(f"  Decay events     = {tta_stats['n_decay']}\n")
            f.write(f"  Mean TTA time    = {tta_stats['mean_ann_time']:.3e} s\n")
            f.write(f"  Interlayer TTA   = {tta_stats['f_interlayer_TTA']*100:.1f}%\n")
            f.write(f"  TTA at GB        = {tta_stats['f_gb_TTA']*100:.1f}%\n")
            f.write(f"  Mean diff length = {tta_stats['mean_diff_length']:.2f} Ang\n")

    print(f"  [OK] Results saved to: {filepath}")


# =============================================================================
# MODULE 5D: BIPOLAR TRANSPORT — HOLES + ELECTRONS SIMULTANEOUSLY
# =============================================================================
# Simulates both carrier types on the same morphology under a common field.
# Relevant for OPV active layers and recombination studies.
# Uses BKL KMC with two separate rate matrices (hole Marcus, electron Marcus).
# Records: mu_hole, mu_electron, recombination events, Langevin rate.
# =============================================================================

