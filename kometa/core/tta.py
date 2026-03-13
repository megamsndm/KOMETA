# =============================================================================
# KOMETA — kometa/core/tta.py
# =============================================================================
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

from .constants import KB, HBAR, ANG_TO_CM, E_CHARGE

def run_TTA_simulation(system, dexter_matrix, params, is_quencher=None):
    """
    Two-body BKL TTA simulation with optional impurity quenching (Module 4B).

    When enable_quenching=True, each step also includes a quenching event
    for whichever exciton currently sits on a quencher site.
    Outcome 'quench_A' or 'quench_B' means that exciton was irreversibly
    lost to an impurity dark state -- no photon emitted.
    """
    print("\n" + "="*60)
    print("  MODULE 6: Running TTA simulation...")
    print("="*60)

    k_decay   = params['k_decay']
    k_TTA     = params['k_TTA']
    num_pairs = params['num_triplet_pairs']
    max_steps = params['tta_max_steps']
    max_time  = params['tta_max_time']
    f_sing    = params['singlet_fraction']
    phi_em    = params['phi_emission']

    enable_q      = params.get('enable_quenching', False)
    k_quench      = params.get('k_quench', 1e10)
    quench_on_hop = params.get('quench_on_hop', False)

    total_sites  = system['total_sites']
    neighbors    = system['neighbors']
    positions    = system['positions']
    site_layer   = system['site_layer']
    site_species = system['site_species']
    is_gb        = system['is_grain_boundary']

    # Quencher array (all False if quenching disabled)
    if is_quencher is None:
        is_quencher = np.zeros(total_sites, dtype=bool)

    neighbor_set = [set() for _ in range(total_sites)]
    for i in range(total_sites):
        for (_, j, _) in neighbors[i]:
            neighbor_set[i].add(j)

    sp_k_decay = params.get('species_k_decay', [k_decay])

    outcomes           = []
    annihilation_times = []
    diffusion_lengths  = []
    layer_crossings_A  = []
    layer_crossings_B  = []
    tta_layer_pairs    = []
    tta_at_gb          = []   # NEW: was TTA at a grain boundary?

    np.random.seed(params['random_seed'] + 9999
                   if params['random_seed'] else None)

    print(f"  Simulating {num_pairs} triplet pairs...")

    for pair_idx in range(num_pairs):
        if (pair_idx + 1) % 50 == 0:
            print(f"  -> Pair {pair_idx+1}/{num_pairs} ...", end='\r')

        site_A = np.random.randint(0, total_sites)
        site_B = np.random.randint(0, total_sites)
        while site_B == site_A:
            site_B = np.random.randint(0, total_sites)

        start_pos_A = positions[site_A].copy()
        start_pos_B = positions[site_B].copy()

        sim_time  = 0.0
        outcome   = 'timeout'
        n_cross_A = 0
        n_cross_B = 0

        # Species-specific decay rates
        kd_A = sp_k_decay[min(site_species[site_A], len(sp_k_decay)-1)]
        kd_B = sp_k_decay[min(site_species[site_B], len(sp_k_decay)-1)]

        for step in range(max_steps):
            if sim_time > max_time:
                outcome = 'timeout'
                break

            events = []

            for (rate_hop, dest) in dexter_matrix[site_A]:
                if dest != site_B:
                    events.append((rate_hop, 'hop_A', dest))

            for (rate_hop, dest) in dexter_matrix[site_B]:
                if dest != site_A:
                    events.append((rate_hop, 'hop_B', dest))

            events.append((kd_A, 'decay_A', None))
            events.append((kd_B, 'decay_B', None))

            if site_B in neighbor_set[site_A]:
                # TTA rate enhanced at grain boundaries (higher local density)
                k_tta_local = k_TTA
                if is_gb[site_A] or is_gb[site_B]:
                    k_tta_local *= 2.0  # local concentration enhancement at GB
                events.append((k_tta_local, 'TTA', None))

            # ── MODULE 4B: impurity quenching events ──────────────────────
            # Added to BKL list only when exciton sits on a quencher site.
            # k_quench >> k_decay ensures rapid irreversible loss.
            if enable_q and not quench_on_hop:
                if is_quencher[site_A]:
                    events.append((k_quench, 'quench_A', None))
                if is_quencher[site_B]:
                    events.append((k_quench, 'quench_B', None))
            # ─────────────────────────────────────────────────────────────

            if not events:
                outcome = 'timeout'
                break

            rates_arr = np.array([ev[0] for ev in events])
            R_total   = np.sum(rates_arr)
            if R_total <= 0:
                outcome = 'timeout'
                break

            u1  = np.random.uniform(0, 1)
            cum = np.cumsum(rates_arr)
            idx = min(np.searchsorted(cum, u1*R_total), len(events)-1)
            chosen_rate, event_type, event_info = events[idx]

            u2 = max(np.random.uniform(0, 1), 1e-300)
            sim_time += -np.log(u2) / R_total

            if event_type == 'hop_A':
                old_layer = site_layer[site_A]
                site_A    = event_info
                kd_A = sp_k_decay[min(site_species[site_A], len(sp_k_decay)-1)]
                if site_layer[site_A] != old_layer:
                    n_cross_A += 1
                # instantaneous quenching mode: check on arrival
                if enable_q and quench_on_hop and is_quencher[site_A]:
                    outcome = 'quench_A'
                    break

            elif event_type == 'hop_B':
                old_layer = site_layer[site_B]
                site_B    = event_info
                kd_B = sp_k_decay[min(site_species[site_B], len(sp_k_decay)-1)]
                if site_layer[site_B] != old_layer:
                    n_cross_B += 1
                # instantaneous quenching mode: check on arrival
                if enable_q and quench_on_hop and is_quencher[site_B]:
                    outcome = 'quench_B'
                    break

            elif event_type == 'decay_A':
                outcome = 'decay_A'
                break

            elif event_type == 'decay_B':
                outcome = 'decay_B'
                break

            elif event_type == 'TTA':
                outcome = 'TTA'
                annihilation_times.append(sim_time)
                tta_layer_pairs.append((site_layer[site_A], site_layer[site_B]))
                tta_at_gb.append(is_gb[site_A] or is_gb[site_B])
                break

            elif event_type == 'quench_A':
                outcome = 'quench_A'
                break

            elif event_type == 'quench_B':
                outcome = 'quench_B'
                break

        dist_A = np.linalg.norm(positions[site_A] - start_pos_A)
        dist_B = np.linalg.norm(positions[site_B] - start_pos_B)
        diffusion_lengths.append((dist_A, dist_B))
        layer_crossings_A.append(n_cross_A)
        layer_crossings_B.append(n_cross_B)
        outcomes.append(outcome)

    n_tta     = outcomes.count('TTA')
    n_quench  = outcomes.count('quench_A') + outcomes.count('quench_B')
    n_decay   = outcomes.count('decay_A')  + outcomes.count('decay_B')
    n_timeout = outcomes.count('timeout')
    print(f"\n  [OK] TTA simulation complete. {num_pairs} pairs simulated.")
    if enable_q:
        print(f"  Quenching summary: TTA={n_tta}, Quenched={n_quench}, "
              f"Decay={n_decay}, Timeout={n_timeout}")
        print(f"  Quench fraction   : {100*n_quench/num_pairs:.1f}%  "
              f"(quencher conc = {params.get('quencher_fraction',0)*100:.1f}%)")

    return {
        'outcomes'           : outcomes,
        'annihilation_times' : annihilation_times,
        'diffusion_lengths'  : diffusion_lengths,
        'layer_crossings_A'  : layer_crossings_A,
        'layer_crossings_B'  : layer_crossings_B,
        'tta_layer_pairs'    : tta_layer_pairs,
        'tta_at_gb'          : tta_at_gb,
        'num_pairs'          : num_pairs,
        'f_singlet'          : f_sing,
        'phi_emission'       : phi_em,
        'n_quenched'         : n_quench,
        'quencher_fraction'  : params.get('quencher_fraction', 0.0),
    }


# =============================================================================
# MODULE 6B: SINGLE TRIPLET DIFFUSION + DIFFUSION COEFFICIENT
# =============================================================================
# Simulates individual triplet excitons (no partner, no TTA).
# Extracts triplet diffusion coefficient D_T from MSD = 6*D_T*t (3D).
# Also estimates bimolecular rate constant: gamma_TTA ~ 4*pi*R_cap*D_T
# This is the essential first step before full TTA simulation.
# =============================================================================

def calculate_triplet_diffusion(system, dexter_matrix, params, is_quencher=None):
    """
    Simulate single triplet diffusion and extract D_T.
    Also estimate bimolecular TTA rate constant gamma_TTA.
    Supports impurity quenching (Module 4B): quenched excitons recorded
    separately and shorten the effective diffusion length.
    """
    print("\n" + "="*60)
    print("  MODULE 6B: Single triplet diffusion coefficient...")
    print("="*60)

    k_decay   = params['k_decay']
    n_exc     = params.get('diffusion_num_excitons', 100)
    max_steps = params.get('diffusion_max_steps', int(1e5))

    enable_q      = params.get('enable_quenching', False)
    k_quench      = params.get('k_quench', 1e10)
    quench_on_hop = params.get('quench_on_hop', False)

    total_sites  = system['total_sites']
    positions    = system['positions']
    site_layer   = system['site_layer']
    neighbor_set = [set(j for (_, j, _) in system['neighbors'][i])
                    for i in range(total_sites)]

    if is_quencher is None:
        is_quencher = np.zeros(total_sites, dtype=bool)

    np.random.seed(params['random_seed'] + 12345
                   if params['random_seed'] else None)

    msd_list        = []
    lifetime_list   = []
    D_estimates     = []
    n_quenched_6b   = 0

    print(f"  Simulating {n_exc} single triplet excitons...")

    for exc_idx in range(n_exc):
        site = np.random.randint(0, total_sites)
        start_pos = positions[site].copy()

        sim_time   = 0.0
        trajectory = [positions[site].copy()]
        time_track = [0.0]
        quenched   = False

        for step in range(max_steps):
            events_full = [(r, 'hop', d) for (r, d) in dexter_matrix[site]]
            events_full.append((k_decay, 'decay', None))

            # Module 4B: add quenching event if on quencher site (BKL mode)
            if enable_q and not quench_on_hop and is_quencher[site]:
                events_full.append((k_quench, 'quench', None))

            rates_arr = np.array([e[0] for e in events_full])
            R_total   = np.sum(rates_arr)
            if R_total <= 0:
                break

            u1  = np.random.uniform(0, 1)
            cum = np.cumsum(rates_arr)
            idx = min(np.searchsorted(cum, u1*R_total), len(events_full)-1)
            ev  = events_full[idx]

            u2 = max(np.random.uniform(0, 1), 1e-300)
            sim_time += -np.log(u2) / R_total

            if ev[1] == 'decay':
                break
            elif ev[1] == 'quench':
                quenched = True
                n_quenched_6b += 1
                break
            else:
                site = ev[2]
                trajectory.append(positions[site].copy())
                time_track.append(sim_time)
                # Module 4B: instantaneous quenching on hop arrival
                if enable_q and quench_on_hop and is_quencher[site]:
                    quenched = True
                    n_quenched_6b += 1
                    break

        lifetime_list.append(sim_time)
        traj_arr = np.array(trajectory)
        if len(traj_arr) > 2 and sim_time > 0:
            disp_sq = np.sum((traj_arr[-1] - traj_arr[0])**2)  # Ang^2
            msd_list.append(disp_sq)
            D_estimates.append(disp_sq * 1e-16 / (6 * sim_time))

    D_T       = np.mean(D_estimates) if D_estimates else 0.0
    D_T_std   = np.std(D_estimates)  if D_estimates else 0.0
    mean_life = np.mean(lifetime_list)

    # Estimate diffusion length: L_diff = sqrt(6 * D_T * tau)
    L_diff    = np.sqrt(6 * D_T * mean_life) * 1e4  # cm -> um ... back to Ang
    L_diff_ang = np.sqrt(6 * D_T * mean_life) * 1e8  # cm -> Ang

    # Smoluchowski bimolecular rate: gamma = 4*pi*R_cap*D_T
    a_lattice = params['intra_layer_distance'] * 1e-8   # Ang -> cm
    gamma_TTA = 4 * np.pi * a_lattice * D_T

    print(f"  [OK] Excitons simulated    : {n_exc}")
    print(f"  [OK] Mean triplet lifetime : {mean_life:.3e} s")
    print(f"  [OK] D_T (triplet diff.)   : {D_T:.4e} cm^2/s")
    print(f"  [OK] D_T std               : {D_T_std:.4e} cm^2/s")
    print(f"  [OK] L_diff (diffusion)    : {L_diff_ang:.1f} Ang")
    print(f"  [OK] gamma_TTA (Smoluch.)  : {gamma_TTA:.4e} cm^3/s")
    if enable_q:
        print(f"  [OK] Quenched excitons     : {n_quenched_6b}/{n_exc} "
              f"({100*n_quenched_6b/n_exc:.1f}%)  -- impurity quenching active")

    return {
        'D_T'          : D_T,
        'D_T_std'      : D_T_std,
        'mean_lifetime' : mean_life,
        'L_diff_ang'   : L_diff_ang,
        'gamma_TTA'    : gamma_TTA,
        'msd_list'     : msd_list,
        'D_estimates'  : D_estimates,
    }


# =============================================================================
# MODULE 7: TTA ANALYSIS
# =============================================================================



def analyze_TTA(tta_results, params):
    print("\n" + "="*60)
    print("  MODULE 7: Analyzing TTA results...")
    print("="*60)

    outcomes    = tta_results['outcomes']
    num_pairs   = tta_results['num_pairs']
    f_singlet   = tta_results['f_singlet']
    phi_em      = tta_results['phi_emission']
    ann_times   = tta_results['annihilation_times']
    diff_lens   = tta_results['diffusion_lengths']
    layer_pairs = tta_results['tta_layer_pairs']
    tta_at_gb   = tta_results.get('tta_at_gb', [])

    n_TTA     = outcomes.count('TTA')
    n_decay_A = outcomes.count('decay_A')
    n_decay_B = outcomes.count('decay_B')
    n_timeout = outcomes.count('timeout')

    P_TTA   = n_TTA / num_pairs
    phi_UC  = 0.5 * P_TTA * f_singlet * phi_em
    phi_UC_max = 0.5 * 1.0 * f_singlet * phi_em

    mean_ann_time = np.mean(ann_times) if ann_times else 0.0
    std_ann_time  = np.std(ann_times)  if ann_times else 0.0

    n_interlayer_TTA = sum(1 for (la, lb) in layer_pairs if la != lb)
    f_interlayer_TTA = n_interlayer_TTA / n_TTA if n_TTA > 0 else 0.0

    n_gb_TTA  = sum(1 for x in tta_at_gb if x)
    f_gb_TTA  = n_gb_TTA / n_TTA if n_TTA > 0 else 0.0

    if diff_lens:
        all_dists = [d for pair in diff_lens for d in pair]
        mean_diff_length = np.mean(all_dists)
    else:
        mean_diff_length = 0.0

    # Binomial uncertainty on P(TTA)
    P_TTA_err = np.sqrt(P_TTA*(1-P_TTA)/num_pairs) if num_pairs > 0 else 0

    print(f"\n  +------------------------------------------------------+")
    print(f"  |                TTA RESULTS SUMMARY                   |")
    print(f"  +------------------------------------------------------+")
    print(f"  |  Pairs simulated           : {num_pairs:<6d}                |")
    print(f"  |  TTA annihilation events   : {n_TTA:<6d} ({P_TTA*100:.1f}%)        |")
    print(f"  |  Natural decay events      : {n_decay_A+n_decay_B:<6d} ({(n_decay_A+n_decay_B)/num_pairs*100:.1f}%)        |")
    print(f"  |  Timeout events            : {n_timeout:<6d}                |")
    print(f"  +------------------------------------------------------+")
    print(f"  |  P(TTA)         = {P_TTA:.4f} +/- {P_TTA_err:.4f}           |")
    print(f"  |  f_singlet      = {f_singlet:.4f}                          |")
    print(f"  |  phi_emission   = {phi_em:.4f}                          |")
    print(f"  |  Phi_UC         = {phi_UC:.6f}                        |")
    print(f"  |  Max Phi_UC     = {phi_UC_max:.6f}                        |")
    print(f"  +------------------------------------------------------+")
    print(f"  |  Mean TTA time  = {mean_ann_time:.3e} s                   |")
    print(f"  |  Interlayer TTA = {f_interlayer_TTA*100:.1f}%                          |")
    print(f"  |  TTA at GB      = {f_gb_TTA*100:.1f}% of all TTA events        |")
    print(f"  |  Mean diff len  = {mean_diff_length:.2f} Ang                       |")
    print(f"  +------------------------------------------------------+")

    return {
        'P_TTA'             : P_TTA,
        'P_TTA_err'         : P_TTA_err,
        'phi_UC'            : phi_UC,
        'phi_UC_max'        : phi_UC_max,
        'n_TTA'             : n_TTA,
        'n_decay'           : n_decay_A + n_decay_B,
        'n_timeout'         : n_timeout,
        'mean_ann_time'     : mean_ann_time,
        'std_ann_time'      : std_ann_time,
        'f_interlayer_TTA'  : f_interlayer_TTA,
        'f_intralayer_TTA'  : 1.0 - f_interlayer_TTA,
        'f_gb_TTA'          : f_gb_TTA,
        'mean_diff_length'  : mean_diff_length,
        'ann_times'         : ann_times,
        'layer_pairs'       : layer_pairs,
    }


# =============================================================================
# MODULE 8: CHARGE TRANSPORT PROPERTIES
# =============================================================================

