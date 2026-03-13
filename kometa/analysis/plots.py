# =============================================================================
# KOMETA — kometa/analysis/plots.py
# =============================================================================
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import warnings
warnings.filterwarnings('ignore')

def plot_results(all_results, transport, system, params,
                 tta_results=None, tta_stats=None,
                 diffusion_stats=None, domain_stats=None,
                 zero_field_stats=None):
    # Helper: should we generate this plot?
    def _want(key, default=True):
        return params.get(key, default)

    print("\n" + "="*60)
    print("  MODULE 9: Generating plots...")
    print("="*60)

    os.makedirs(params['output_folder'], exist_ok=True)
    plt.rcParams.update({'font.size': 12, 'axes.labelsize': 13,
                         'axes.titlesize': 14, 'figure.dpi': 150})
    colors = plt.cm.tab10(np.linspace(0, 1, max(system['num_layers'], 2)))
    gb_color   = '#e74c3c'
    crys_color = '#2ecc71'

    # ---------------------------------------------------------------
    # CHARGE TRANSPORT PLOTS
    # ---------------------------------------------------------------
    if all_results is not None and transport is not None:

        if _want('plot_dos'):
            # Plot 1: DOS + grain boundary overlay
            fig, ax = plt.subplots(figsize=(8, 5))
            for layer_idx in range(system['num_layers']):
                mask_layer = system['site_layer'] == layer_idx
                # Crystalline sites
                mask_cryst = mask_layer & ~system['is_grain_boundary']
                mask_gb    = mask_layer & system['is_grain_boundary']
                E_cryst = system['site_energies'][mask_cryst]
                E_gb    = system['site_energies'][mask_gb]
                if len(E_cryst) > 0:
                    ax.hist(E_cryst, bins=20, alpha=0.55, color=colors[layer_idx],
                            density=True, label=f'Layer {layer_idx+1} (cryst)')
                if len(E_gb) > 0:
                    ax.hist(E_gb, bins=15, alpha=0.4, color=gb_color,
                            density=True, label=f'Layer {layer_idx+1} (GB)',
                            histtype='step', linewidth=2, linestyle='--')
            ax.set_xlabel('Site Energy (eV)')
            ax.set_ylabel('Density of States')
            ax.set_title(f"Charge DOS  [{params.get('charge_rate_model','marcus').upper()}]")
            ax.legend(fontsize=10); ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(params['output_folder'], 'plot1_DOS.png'))
            plt.close()
            print("  [OK] Plot 1: DOS (with GB overlay)")

        if _want('plot_trajectory'):
            # Plot 2: Trajectory
            fig, ax = plt.subplots(figsize=(9, 5))
            for i, run in enumerate(all_results[:5]):
                ax.plot(run['time']*1e12, run['displacement'], alpha=0.7, lw=1,
                        label=f'Run {i+1}' if i < 3 else '')
            ax.set_xlabel('Time (ps)'); ax.set_ylabel('Displacement (Ang)')
            ax.set_title(f"Carrier Trajectory  [{params.get('charge_rate_model','marcus').upper()}]")
            ax.legend(); ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(params['output_folder'], 'plot2_trajectory.png'))
            plt.close()
            print("  [OK] Plot 2: Trajectory")

        if _want('plot_layer_occ'):
            # Plot 3: Layer occupancy
            fig, axes = plt.subplots(system['num_layers'], 1,
                                     figsize=(9, 3*system['num_layers']), sharex=True)
            if system['num_layers'] == 1:
                axes = [axes]
            run = all_results[0]
            for li in range(system['num_layers']):
                axes[li].fill_between(run['time']*1e12,
                                      (run['layer']==li).astype(float),
                                      alpha=0.7, color=colors[li])
                axes[li].set_ylabel(f'Layer {li+1}')
                axes[li].set_ylim(-0.05, 1.15)
                axes[li].grid(True, alpha=0.3)
            axes[-1].set_xlabel('Time (ps)')
            fig.suptitle('Layer Occupancy')
            plt.tight_layout()
            plt.savefig(os.path.join(params['output_folder'], 'plot3_layer_occupancy.png'))
            plt.close()
            print("  [OK] Plot 3: Layer occupancy")

        if _want('plot_energy'):
            # Plot 4: Energy vs time
            fig, ax = plt.subplots(figsize=(9, 5))
            run = all_results[0]
            ax.plot(run['time']*1e12, run['energy'], color='darkred', lw=0.8)
            ax.axhline(np.mean(run['energy']), color='k', ls='--', lw=1.5,
                       label=f"Mean = {np.mean(run['energy']):.3f} eV")
            ax.set_xlabel('Time (ps)'); ax.set_ylabel('Site Energy (eV)')
            ax.set_title('Carrier Energy vs Time')
            ax.legend(); ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(params['output_folder'], 'plot4_energy_vs_time.png'))
            plt.close()
            print("  [OK] Plot 4: Energy vs time")

        if _want('plot_mob_dist'):
            # Plot 5: Mobility distribution
            fig, ax = plt.subplots(figsize=(7, 5))
            ax.hist(transport['mobility_all'],
                    bins=max(5, len(transport['mobility_all'])//2),
                    color='steelblue', edgecolor='white', alpha=0.85)
            ax.axvline(transport['mobility_mean'], color='red', lw=2, ls='--',
                       label=f"Mean = {transport['mobility_mean']:.3e}")
            ax.set_xlabel('Mobility (cm^2/(V.s))'); ax.set_ylabel('Count')
            ax.set_title(f"Mobility Distribution  [{params.get('charge_rate_model','marcus').upper()}]")
            ax.legend(); ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(params['output_folder'], 'plot5_mobility_distribution.png'))
            plt.close()
            print("  [OK] Plot 5: Mobility distribution")

        if _want('plot_layer_res'):
            # Plot 6: Layer residence
            fig, ax = plt.subplots(figsize=(7, 5))
            fracs = transport['layer_fractions'] * 100
            bars = ax.bar([f'Layer {i+1}' for i in range(system['num_layers'])],
                          fracs, color=colors[:system['num_layers']], edgecolor='white')
            for bar, frac in zip(bars, fracs):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                        f'{frac:.1f}%', ha='center', fontsize=12)
            ax.set_ylabel('Residence Time (%)')
            ax.set_title('Layer Residence Fractions')
            ax.set_ylim(0, max(fracs)*1.2); ax.grid(True, axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(params['output_folder'], 'plot6_layer_residence.png'))
            plt.close()
            print("  [OK] Plot 6: Layer residence")

        if _want('plot_domain'):
            # Plot 6B: Domain-resolved mobility (if available)
            if domain_stats is not None:
                fig, axes = plt.subplots(1, 2, figsize=(12, 5))
                # Left: mobility comparison
                ax = axes[0]
                labels = ['Crystalline', 'Grain Boundary']
                vals   = [domain_stats['mu_crystalline'],
                          domain_stats['mu_grain_boundary']]
                bar_c  = [crys_color, gb_color]
                brs = ax.bar(labels, vals, color=bar_c, edgecolor='white', width=0.5)
                for bar, v in zip(brs, vals):
                    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()*1.05,
                            f'{v:.2e}', ha='center', fontsize=10)
                ax.set_ylabel('Mobility (cm^2/(V.s))')
                ax.set_title('Domain-Resolved Mobility')
                ax.set_yscale('log'); ax.grid(True, axis='y', alpha=0.3)
                # Right: trajectory fraction
                ax = axes[1]
                fracs_d = [domain_stats['f_crystalline']*100,
                           domain_stats['f_grain_boundary']*100]
                ax.bar(labels, fracs_d, color=bar_c, edgecolor='white', width=0.5)
                ax.set_ylabel('Fraction of Trajectory (%)')
                ax.set_title('Time Spent in Each Domain Type')
                ax.grid(True, axis='y', alpha=0.3)
                plt.tight_layout()
                plt.savefig(os.path.join(params['output_folder'],
                                         'plot6b_domain_mobility.png'))
                plt.close()
                print("  [OK] Plot 6b: Domain-resolved mobility")

    # ---------------------------------------------------------------
    # TTA PLOTS
    # ---------------------------------------------------------------
    if tta_results is not None and tta_stats is not None:

        outcomes    = tta_results['outcomes']
        ann_times   = tta_stats['ann_times']
        layer_pairs = tta_stats['layer_pairs']

        if _want('plot_tta_pie'):
            # Plot 7: Outcome pie
            fig, ax = plt.subplots(figsize=(7, 6))
            n_TTA     = tta_stats['n_TTA']
            n_decay   = tta_stats['n_decay']
            n_timeout = tta_stats['n_timeout']
            labels    = ['TTA (annihilation)', 'Natural decay', 'Timeout']
            sizes     = [n_TTA, n_decay, n_timeout]
            c_pie     = ['#2ecc71', '#e74c3c', '#95a5a6']
            non_zero  = [(l, s, c) for l, s, c in zip(labels, sizes, c_pie) if s > 0]
            if non_zero:
                labels_nz, sizes_nz, colors_nz = zip(*non_zero)
                wedges, texts, autotexts = ax.pie(
                    sizes_nz, labels=labels_nz, colors=colors_nz,
                    autopct='%1.1f%%', startangle=90,
                    wedgeprops={'edgecolor': 'white', 'linewidth': 1.5})
            ax.set_title(f'TTA Outcomes (N={tta_results["num_pairs"]} pairs)\n'
                         f'Phi_UC = {tta_stats["phi_UC"]:.5f}')
            plt.tight_layout()
            plt.savefig(os.path.join(params['output_folder'], 'plot7_TTA_outcomes.png'))
            plt.close()
            print("  [OK] Plot 7: TTA outcome pie")

        if _want('plot_tta_time'):
            # Plot 8: Annihilation time distribution
            if len(ann_times) > 5:
                fig, ax = plt.subplots(figsize=(8, 5))
                ax.hist(np.array(ann_times)*1e6, bins=30,
                        color='#e67e22', edgecolor='white', alpha=0.85)
                ax.axvline(tta_stats['mean_ann_time']*1e6, color='darkred',
                           lw=2, ls='--',
                           label=f'Mean = {tta_stats["mean_ann_time"]*1e6:.3f} mus')
                ax.set_xlabel('Time to Annihilation (mus)')
                ax.set_ylabel('Count')
                ax.set_title('Annihilation Time Distribution')
                ax.legend(); ax.grid(True, alpha=0.3)
                plt.tight_layout()
                plt.savefig(os.path.join(params['output_folder'],
                                         'plot8_TTA_annihilation_times.png'))
                plt.close()
                print("  [OK] Plot 8: Annihilation time distribution")

        if _want('plot_tta_layer'):
            # Plot 9: Interlayer + GB breakdown
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            ax = axes[0]
            bar_labels = ['Intralayer TTA', 'Interlayer TTA']
            bar_vals   = [tta_stats['f_intralayer_TTA']*100,
                          tta_stats['f_interlayer_TTA']*100]
            ax.bar(bar_labels, bar_vals, color=['#3498db', '#e74c3c'],
                   edgecolor='white', width=0.5)
            for i, v in enumerate(bar_vals):
                ax.text(i, v+0.5, f'{v:.1f}%', ha='center', fontweight='bold')
            ax.set_ylabel('Fraction of TTA Events (%)')
            ax.set_title('Intralayer vs Interlayer TTA')
            ax.set_ylim(0, 120); ax.grid(True, axis='y', alpha=0.3)
            ax = axes[1]
            gb_frac    = tta_stats['f_gb_TTA'] * 100
            cryst_frac = (1 - tta_stats['f_gb_TTA']) * 100
            ax.bar(['Crystalline sites', 'Grain boundary sites'],
                   [cryst_frac, gb_frac],
                   color=[crys_color, gb_color], edgecolor='white', width=0.5)
            for i, v in enumerate([cryst_frac, gb_frac]):
                ax.text(i, v+0.5, f'{v:.1f}%', ha='center', fontweight='bold')
            ax.set_ylabel('Fraction of TTA Events (%)')
            ax.set_title('TTA at Crystalline vs Grain Boundary Sites')
            ax.set_ylim(0, 120); ax.grid(True, axis='y', alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(params['output_folder'], 'plot9_TTA_layer_breakdown.png'))
            plt.close()
            print("  [OK] Plot 9: TTA layer + GB breakdown")

        if _want('plot_phi_uc'):
            # Plot 10: Phi_UC breakdown
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            ax = axes[0]
            steps     = ['P(TTA)', 'x f_sing', 'x phi_em', 'x 0.5 = Phi_UC']
            step_vals = [tta_stats['P_TTA'],
                         tta_stats['P_TTA']*params['singlet_fraction'],
                         tta_stats['P_TTA']*params['singlet_fraction']*params['phi_emission'],
                         tta_stats['phi_UC']]
            ax.bar(steps, step_vals,
                   color=['#3498db','#9b59b6','#e67e22','#2ecc71'],
                   edgecolor='white')
            for i, v in enumerate(step_vals):
                ax.text(i, v+0.002, f'{v:.4f}', ha='center', fontsize=9, fontweight='bold')
            ax.set_ylabel('Value')
            ax.set_title('Phi_UC Step-by-Step Breakdown')
            ax.set_ylim(0, max(step_vals)*1.3); ax.grid(True, axis='y', alpha=0.3)
            plt.setp(ax.get_xticklabels(), fontsize=9)
            ax = axes[1]
            for layer_idx in range(system['num_layers']):
                mask = system['site_layer'] == layer_idx
                t_E  = system['triplet_energies'][mask]
                ax.hist(t_E, bins=20, alpha=0.5, color=colors[layer_idx],
                        density=True, label=f'Layer {layer_idx+1}')
            ax.set_xlabel('Triplet Energy T1 (eV)')
            ax.set_ylabel('Density')
            ax.set_title('Triplet DOS')
            ax.legend(); ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(os.path.join(params['output_folder'],
                                     'plot10_TTA_quantum_yield.png'))
            plt.close()
            print("  [OK] Plot 10: Phi_UC summary")

    # ---------------------------------------------------------------
    # DIFFUSION COEFFICIENT PLOTS
    # ---------------------------------------------------------------
    if diffusion_stats is not None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        ax = axes[0]
        D_est = diffusion_stats['D_estimates']
        ax.hist(np.array(D_est)*1e9, bins=25,
                color='#8e44ad', edgecolor='white', alpha=0.85)
        ax.axvline(diffusion_stats['D_T']*1e9, color='k', lw=2, ls='--',
                   label=f'D_T = {diffusion_stats["D_T"]:.3e} cm^2/s')
        ax.set_xlabel('D_T (x10^-9 cm^2/s)')
        ax.set_ylabel('Count')
        ax.set_title('Triplet Diffusion Coefficient Distribution')
        ax.legend(); ax.grid(True, alpha=0.3)
        ax = axes[1]
        metrics = ['D_T\n(cm^2/s)', 'gamma_TTA\n(cm^3/s)', 'L_diff\n(Ang)']
        values  = [diffusion_stats['D_T'],
                   diffusion_stats['gamma_TTA'],
                   diffusion_stats['L_diff_ang']]
        bar_c   = ['#8e44ad', '#e74c3c', '#27ae60']
        for idx, (m, v, c) in enumerate(zip(metrics, values, bar_c)):
            ax2 = axes[1] if idx == 0 else axes[1].twinx() if idx == 1 else axes[1]
        # Simple text summary instead
        ax.axis('off')
        summary = (f"Triplet Diffusion Summary\n\n"
                   f"D_T = {diffusion_stats['D_T']:.4e} cm^2/s\n"
                   f"D_T std = {diffusion_stats['D_T_std']:.4e} cm^2/s\n\n"
                   f"Lifetime = {diffusion_stats['mean_lifetime']:.3e} s\n"
                   f"L_diff = {diffusion_stats['L_diff_ang']:.1f} Ang\n\n"
                   f"gamma_TTA = {diffusion_stats['gamma_TTA']:.4e} cm^3/s")
        ax.text(0.1, 0.5, summary, transform=ax.transAxes,
                fontsize=12, verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='#f0e6ff', alpha=0.8))
        ax.set_title('Diffusion Summary')
        plt.tight_layout()
        plt.savefig(os.path.join(params['output_folder'],
                                 'plot11_triplet_diffusion.png'))
        plt.close()
        print("  [OK] Plot 11: Triplet diffusion coefficient")

    # ---------------------------------------------------------------
    # ZERO-FIELD MOBILITY PLOT
    # ---------------------------------------------------------------
    if zero_field_stats is not None and transport is not None:
        fig, ax = plt.subplots(figsize=(7, 5))
        labels = ['Zero-field mu_0\n(Einstein)', 'Finite-field mu\n(drift)']
        vals   = [zero_field_stats['mu_zero_field'],
                  transport['mobility_mean']]
        colors_bar = ['#27ae60', '#2980b9']
        bars = ax.bar(labels, vals, color=colors_bar, edgecolor='white', width=0.5)
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()*1.05,
                    f'{v:.3e}', ha='center', fontsize=11, fontweight='bold')
        ax.set_ylabel('Mobility (cm^2/(V.s))')
        ax.set_title('Zero-Field vs Finite-Field Mobility Comparison')
        ax.set_yscale('log'); ax.grid(True, axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(params['output_folder'],
                                 'plot12_zero_field_mobility.png'))
        plt.close()
        print("  [OK] Plot 12: Zero-field vs finite-field mobility")

    print(f"\n  All plots saved to: '{params['output_folder']}/'")


# =============================================================================
# MODULE 10: SAVE RESULTS
# =============================================================================



def plot_extended_results(params, temp_scan=None, tof_stats=None,
                          energy_corr=None, diffusion_tensor=None,
                          seebeck=None, intensity_scan=None,
                          singlet_stats=None, tpq_stats=None,
                          bipolar_stats=None):
    """Generate additional plots for new modules."""
    import os, matplotlib.pyplot as plt
    os.makedirs(params['output_folder'], exist_ok=True)
    folder = params['output_folder']
    plotted = []

    # ── Temperature scan μ(T) ────────────────────────────────────────────────
    if temp_scan and len(temp_scan['temperatures']) > 1:
        fig, ax = plt.subplots(figsize=(6, 4))
        T_arr  = np.array(temp_scan['temperatures'])
        mu_arr = np.array(temp_scan['mobilities'])
        ax.semilogy(T_arr, mu_arr, 'o-', color='steelblue', lw=2, ms=7)
        ax.set_xlabel('Temperature (K)', fontsize=12)
        ax.set_ylabel('Mobility μ (cm²/(V·s))', fontsize=12)
        ax.set_title(f"μ(T) — {temp_scan['regime']}\nα = {temp_scan['alpha']:.2f}", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('#F8F9FA')
        plt.tight_layout()
        fname = f"{folder}/plot_temperature_scan.png"
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close()
        plotted.append(fname)

    # ── Site energy correlation C(r) ─────────────────────────────────────────
    if energy_corr is not None:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(energy_corr['r'], energy_corr['C_r'], color='darkorange', lw=2)
        ax.axhline(0, color='gray', lw=1, ls='--')
        ax.set_xlabel('Distance r (Å)', fontsize=12)
        ax.set_ylabel('C(r) = ⟨δε(0)δε(r)⟩/⟨δε²⟩', fontsize=12)
        ax.set_title(f"Site Energy Correlation Function\nξ = {energy_corr['xi']:.2f} Å", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('#F8F9FA')
        plt.tight_layout()
        fname = f"{folder}/plot_energy_correlation.png"
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close()
        plotted.append(fname)

    # ── Intensity scan Phi_UC(I) ──────────────────────────────────────────────
    if intensity_scan is not None:
        fig, ax = plt.subplots(figsize=(6, 4))
        n_arr   = np.array([r['n_pairs'] for r in intensity_scan['scan']])
        phi_arr = np.array([r['phi_UC']  for r in intensity_scan['scan']])
        ax.plot(n_arr, phi_arr, 'o-', color='mediumpurple', lw=2, ms=7)
        ax.set_xlabel('Triplet Pair Density (proxy for I)', fontsize=12)
        ax.set_ylabel('Φ_UC', fontsize=12)
        ax.set_title('Excitation Density Dependence Φ_UC(I)', fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_facecolor('#F8F9FA')
        plt.tight_layout()
        fname = f"{folder}/plot_intensity_scan.png"
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close()
        plotted.append(fname)

    # ── Diffusion tensor bar chart ────────────────────────────────────────────
    if diffusion_tensor is not None:
        fig, ax = plt.subplots(figsize=(5, 4))
        D = diffusion_tensor['D_tensor']
        labels = ['D_xx', 'D_yy', 'D_zz']
        values = [D[0,0], D[1,1], D[2,2]]
        colors = ['steelblue', 'darkorange', 'seagreen']
        ax.bar(labels, values, color=colors, edgecolor='black', linewidth=0.8)
        ax.set_ylabel('Diffusion Coefficient (cm²/s)', fontsize=12)
        ax.set_title('3D Anisotropic Diffusion Tensor\n(diagonal components)', fontsize=11)
        ax.set_facecolor('#F8F9FA')
        plt.tight_layout()
        fname = f"{folder}/plot_diffusion_tensor.png"
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close()
        plotted.append(fname)

    # ── Bipolar transport summary ─────────────────────────────────────────────
    if bipolar_stats is not None:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.bar(['Hole μ', 'Electron μ'],
               [bipolar_stats['mu_hole'], bipolar_stats['mu_electron']],
               color=['royalblue', 'tomato'], edgecolor='black', linewidth=0.8)
        ax.set_ylabel('Mobility (cm²/(V·s))', fontsize=12)
        ax.set_title('Bipolar Transport: Hole vs Electron\n'
                     f"Recombination: {bipolar_stats['recombination_fraction']*100:.1f}%", fontsize=11)
        ax.set_facecolor('#F8F9FA')
        plt.tight_layout()
        fname = f"{folder}/plot_bipolar.png"
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close()
        plotted.append(fname)

    if plotted:
        print(f"\n  [OK] Extended plots saved: {len(plotted)} figures")


# =============================================================================
# FUTURE DEVELOPMENT ROADMAP
# =============================================================================
# The following features are planned for KMC v4.0 and beyond.
# They are listed here as code comments so developers can track progress.
#
# ── PHYSICS ──────────────────────────────────────────────────────────────────
# TODO v4.0: Polaron formation & self-trapping (small polaron model)
#   Physics: Large lambda causes carrier to distort lattice around itself.
#   Method:  Self-consistent Marcus rate with site energy renormalisation.
#   Ref:     DOI 10.1103/PhysRevLett.70.2928
#
# TODO v4.0: Marcus-Hush-Chidsey (MHC) rates for electrode interfaces
#   Physics: Correct hopping rate at metal/organic contact (continuum DOS).
#   Method:  Numerical integral of Fermi-Dirac weighted Marcus rates.
#   Ref:     DOI 10.1021/acs.jpclett.0c03434
#
# TODO v4.0: Quantum coherence effects (Redfield / Lindblad master equation)
#   Physics: Partial coherence between hops when J ~ sigma (borderline).
#   Method:  Density matrix propagation on Hamiltonian with dephasing.
#   Ref:     DOI 10.1063/1.4794784
#
# TODO v4.0: Non-equilibrium hot carrier injection
#   Physics: Carriers injected above band edge relax into DOS.
#   Method:  Start carrier at elevated energy, track thermalisation.
#
# ── MORPHOLOGY ───────────────────────────────────────────────────────────────
# TODO v4.0: True 3D thin-film morphology (lateral disorder)
#   Current code is quasi-1D stack. Extend to full 3D grid.
#
# TODO v4.0: Log-normal crystalline domain size distribution
#   Current: single S2 cutoff -> all GB regions same size.
#   Needed: realistic polycrystalline film with varied domain sizes.
#
# TODO v4.0: ML force field integration (MACE, NequIP)
#   Direct link to Task 4.2 ML-MD morphologies without manual export.
#
# ── CHARGE TRANSPORT ─────────────────────────────────────────────────────────
# TODO v4.0: AC mobility / impedance spectroscopy simulation
#   Frequency-dependent mu(omega) for comparison to CELIV/impedance.
#
# TODO v4.0: Trap-filling kinetics (space-charge-limited conduction)
#   Multi-carrier with trap states -> SCLC J(V) with trap filling.
#
# ── TTA & PHOTOPHYSICS ───────────────────────────────────────────────────────
# TODO v4.0: Delayed fluorescence (TADF) emitter support
#   UC with TADF annihilators; rISC from T1 -> S1 competes with TTA.
#
# TODO v4.0: Photon reabsorption (inner filter effect)
#   Emitted upconverted photons reabsorbed by sensitizer in thick films.
#
# TODO v4.0: Oxygen quenching of triplets
#   k_O2 = k_q * [O2]; important for non-encapsulated devices.
#
# ── ML & DATA ────────────────────────────────────────────────────────────────
# TODO v4.0: ML-predicted transfer integrals (GNN on molecular geometry)
#   Replace ORCA dimer calculations with graph neural network J predictor.
#
# TODO v4.0: Bayesian optimisation of parameters
#   Auto-fit lambda, k0, L to experimental D_T and Phi_UC data.
#
# TODO v4.0: Surrogate model for Phi_UC(parameters)
#   Neural network emulator for high-throughput material screening.
#
# ── SOFTWARE ─────────────────────────────────────────────────────────────────
# TODO v4.0: GPU acceleration via CuPy (rate matrix build + KMC steps)
# TODO v4.0: Parallel trajectory execution (multiprocessing.Pool)
# TODO v4.0: Unit test suite (pytest) for all 20+ modules
# TODO v4.0: Web interface for parameter input and result visualisation
# TODO v4.0: Database of pre-computed lambda, J, T1 for common molecules
# =============================================================================


# =============================================================================
# MODULE 11: MAIN  (updated to call all new modules)
# =============================================================================

