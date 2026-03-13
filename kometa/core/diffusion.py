# =============================================================================
# KOMETA — kometa/core/diffusion.py
# =============================================================================
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

from .constants import KB, HBAR, ANG_TO_CM, E_CHARGE

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

