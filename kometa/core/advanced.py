# =============================================================================
# KOMETA — kometa/core/advanced.py
# =============================================================================
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

from .constants import KB, HBAR, ANG_TO_CM, E_CHARGE

def run_bipolar_kmc(system, rate_matrix_hole, rate_matrix_electron, params):
    """Run a single bipolar KMC trajectory (one hole + one electron)."""
    total_sites  = system['total_sites']
    positions    = system['positions']
    site_layer   = system['site_layer']
    neighbors    = system['neighbors']
    k_rec        = params.get('k_recombination', 1e6)
    num_steps    = params['num_kmc_steps']
    field        = params['electric_field']
    T            = params['temperature']

    rng = np.random.default_rng(params['random_seed'])

    # Place hole and electron on different random sites
    site_h = rng.integers(0, total_sites)
    site_e = rng.integers(0, total_sites)
    while site_e == site_h:
        site_e = rng.integers(0, total_sites)

    z_h0, z_e0 = positions[site_h, 2], positions[site_e, 2]
    sim_time = 0.0
    outcome  = 'timeout'
    neighbor_set = [set(j for (_, j, _) in neighbors[i]) for i in range(total_sites)]

    for step in range(num_steps):
        events = []
        # Hole hops
        for (rate, dest, _) in neighbors[site_h]:
            if dest != site_e:
                events.append((rate_matrix_hole[site_h].get(dest, 0), 'hop_h', dest))
        # Electron hops
        for (rate, dest, _) in neighbors[site_e]:
            if dest != site_h:
                events.append((rate_matrix_electron[site_e].get(dest, 0), 'hop_e', dest))
        # Recombination if neighbours
        if site_e in neighbor_set[site_h]:
            events.append((k_rec, 'recombine', None))

        events = [(r, t, d) for r, t, d in events if r > 0]
        if not events:
            break

        rates = np.array([e[0] for e in events])
        R     = np.sum(rates)
        u1    = rng.uniform(0, 1)
        idx   = min(np.searchsorted(np.cumsum(rates), u1 * R), len(events) - 1)
        _, ev_type, ev_dest = events[idx]

        u2 = max(rng.uniform(0, 1), 1e-300)
        sim_time += -np.log(u2) / R

        if ev_type == 'hop_h':
            site_h = ev_dest
        elif ev_type == 'hop_e':
            site_e = ev_dest
        elif ev_type == 'recombine':
            outcome = 'recombine'
            break

    dz_h = positions[site_h, 2] - z_h0
    dz_e = positions[site_e, 2] - z_e0
    mu_h = (dz_h * 1e-8) / (field * 1e8 * sim_time) if sim_time > 0 else 0
    mu_e = (dz_e * 1e-8) / (field * 1e8 * sim_time) if sim_time > 0 else 0
    return {'mu_hole': abs(mu_h), 'mu_electron': abs(mu_e),
            'outcome': outcome, 'sim_time': sim_time}


def run_bipolar_simulation(system, params):
    """Run ensemble of bipolar KMC trajectories."""
    print("\n" + "="*60)
    print("  MODULE 5D: Bipolar Transport (hole + electron)...")
    print("="*60)

    # Build separate rate matrices for holes and electrons
    params_e = {**params,
                'mean_site_energy': params.get('electron_mean_site_energy', params['mean_site_energy']),
                'lambda_reorg'    : params.get('electron_lambda_reorg',     params['lambda_reorg']),
                'disorder_sigma'  : params.get('electron_disorder_sigma',   params['disorder_sigma']),
                'J_intra'         : params.get('electron_J_intra',          params['J_intra']),
                'J_inter'         : params.get('electron_J_inter',          params['J_inter']),
                'electric_field'  : -params['electric_field']}   # electrons move opposite

    rate_h = calculate_marcus_rates(system, params)
    rate_e = calculate_marcus_rates(system, params_e)

    # Convert rate matrices to dicts for fast lookup
    def to_dict(mat):
        return [{dest: rate for rate, dest in mat[i]} for i in range(len(mat))]

    rh_dict = to_dict(rate_h)
    re_dict = to_dict(rate_e)

    n_runs  = params.get('num_bipolar_runs', 5)
    results = []
    for ri in range(n_runs):
        r = run_bipolar_kmc(system, rh_dict, re_dict, params)
        results.append(r)
        print(f"  Run {ri+1}/{n_runs}: mu_h={r['mu_hole']:.3e}  mu_e={r['mu_electron']:.3e}  outcome={r['outcome']}")

    mu_h_mean = np.mean([r['mu_hole']      for r in results])
    mu_e_mean = np.mean([r['mu_electron']  for r in results])
    n_rec     = sum(1 for r in results if r['outcome'] == 'recombine')

    print(f"\n  [OK] mu_hole     = {mu_h_mean:.4e} cm^2/(V.s)")
    print(f"  [OK] mu_electron = {mu_e_mean:.4e} cm^2/(V.s)")
    print(f"  [OK] Recombination events: {n_rec}/{n_runs}")
    return {'mu_hole': mu_h_mean, 'mu_electron': mu_e_mean,
            'recombination_fraction': n_rec / n_runs, 'raw': results}


# =============================================================================
# MODULE 5E: MULTI-CARRIER KMC (space charge, N simultaneous carriers)
# =============================================================================
# Tracks N holes on the lattice simultaneously with site exclusion.
# No two carriers can occupy the same site (hard-core exclusion).
# Optional: Coulomb repulsion energy correction to hop rates.
# =============================================================================

def run_multicarrier_kmc(system, rate_matrix, params):
    """Run KMC with N simultaneous carriers and site exclusion."""
    print("\n" + "="*60)
    print("  MODULE 5E: Multi-Carrier KMC (space charge)...")
    print("="*60)

    N_c       = params.get('num_carriers', 3)
    total_s   = system['total_sites']
    positions = system['positions']
    num_steps = params['num_kmc_steps']
    field     = params['electric_field']
    use_coul  = not params.get('coulomb_exclusion_only', True)
    eps_r     = params.get('coulomb_epsilon_r', 3.5)
    rng       = np.random.default_rng(params['random_seed'] + 55555)

    # Place N carriers on distinct random sites
    sites = rng.choice(total_s, size=N_c, replace=False).tolist()
    occupied = set(sites)
    z0 = np.array([positions[s, 2] for s in sites])

    kB_eV = 8.617e-5
    T     = params['temperature']

    sim_time = 0.0
    for step in range(num_steps):
        events = []
        for ci, si in enumerate(sites):
            for (rate, dest, _) in system['neighbors'][si]:
                if dest in occupied:
                    continue
                r = rate_matrix[si]
                hop_rate = next((rt for rt, dt in r if dt == dest), 0.0) if isinstance(r, list) else r.get(dest, 0.0)
                if hop_rate <= 0:
                    continue
                # Optional Coulomb correction
                if use_coul:
                    dE_coul = 0.0
                    for cj, sj in enumerate(sites):
                        if cj == ci:
                            continue
                        r_ij = np.linalg.norm(positions[dest] - positions[sj]) * 1e-10
                        r_cur = np.linalg.norm(positions[si]   - positions[sj]) * 1e-10
                        if r_ij > 1e-12 and r_cur > 1e-12:
                            ke = 8.9875e9 * (1.602e-19)**2 / (eps_r * 1.602e-19)
                            dE_coul += ke * (1/r_ij - 1/r_cur)
                    if dE_coul > 0:
                        hop_rate *= np.exp(-dE_coul / (kB_eV * T))
                events.append((hop_rate, ci, dest))

        if not events:
            break

        rates = np.array([e[0] for e in events])
        R     = rates.sum()
        if R <= 0:
            break
        u1  = rng.uniform(0, 1)
        idx = min(np.searchsorted(np.cumsum(rates), u1 * R), len(events) - 1)
        _, ci_chosen, dest_chosen = events[idx]

        u2 = max(rng.uniform(0, 1), 1e-300)
        sim_time += -np.log(u2) / R

        occupied.discard(sites[ci_chosen])
        sites[ci_chosen] = dest_chosen
        occupied.add(dest_chosen)

    dz  = np.array([positions[s, 2] for s in sites]) - z0
    mu_list = [(dz[i] * 1e-8) / (field * 1e8 * sim_time) if sim_time > 0 else 0
               for i in range(N_c)]
    mu_mean = np.mean(np.abs(mu_list))
    print(f"  [OK] N_carriers = {N_c}  |  mu_mean = {mu_mean:.4e} cm^2/(V.s)")
    return {'mu_multicarrier': mu_mean, 'n_carriers': N_c}


# =============================================================================
# MODULE 5F: TIME-OF-FLIGHT (TOF) SIMULATION
# =============================================================================
# Simulates transient photocurrent: carrier injected at bottom layer,
# driven to top layer by applied voltage. Records transit time t_tr.
# mu = L^2 / (V * t_tr)  — directly comparable to TOF experiment.
# =============================================================================

def run_tof_simulation(system, rate_matrix, params):
    """Simulate time-of-flight mobility measurement."""
    print("\n" + "="*60)
    print("  MODULE 5F: Time-of-Flight (TOF) Simulation...")
    print("="*60)

    total_sites = system['total_sites']
    positions   = system['positions']
    site_layer  = system['site_layer']
    neighbors   = system['neighbors']
    n_carriers  = params.get('tof_num_carriers', 50)
    max_time    = params.get('tof_max_time', 1e-6)
    V           = params.get('tof_voltage', 10.0)
    rng         = np.random.default_rng(params['random_seed'] + 33333)

    # Film thickness from positions
    z_min = positions[:, 2].min()
    z_max = positions[:, 2].max()
    L_cm  = (z_max - z_min) * 1e-8   # Ang -> cm

    # Inject carriers at bottom layer (layer 0)
    start_sites = np.where(site_layer == 0)[0]
    transit_times = []

    for ci in range(n_carriers):
        site = rng.choice(start_sites)
        sim_time = 0.0
        transit  = None

        for step in range(params['num_kmc_steps']):
            if sim_time > max_time:
                break
            evts  = [(r, d) for r, d, _ in neighbors[site]]
            rates = np.array([r for r, _ in evts])
            if len(rates) == 0 or rates.sum() <= 0:
                break
            R   = rates.sum()
            u1  = rng.uniform(0, 1)
            idx = min(np.searchsorted(np.cumsum(rates), u1 * R), len(evts) - 1)
            site = evts[idx][1]
            u2 = max(rng.uniform(0, 1), 1e-300)
            sim_time += -np.log(u2) / R
            # Check if reached top layer
            if site_layer[site] == params['num_layers'] - 1:
                transit = sim_time
                break

        if transit:
            transit_times.append(transit)

    if transit_times:
        t_tr = np.mean(transit_times)
        mu_tof = L_cm**2 / (V * t_tr) if t_tr > 0 else 0.0
        print(f"  [OK] Transit carriers : {len(transit_times)}/{n_carriers}")
        print(f"  [OK] Mean t_tr        : {t_tr:.4e} s")
        print(f"  [OK] Film thickness L : {L_cm*1e7:.1f} nm")
        print(f"  [OK] mu_TOF           : {mu_tof:.4e} cm^2/(V.s)")
    else:
        mu_tof = 0.0
        t_tr   = 0.0
        print("  [WARN] No carriers reached top layer — increase tof_max_time or film size")

    return {'mu_tof': mu_tof, 'transit_time_mean': t_tr,
            'n_transit': len(transit_times), 'L_cm': L_cm}


# =============================================================================
# MODULE 5G: TEMPERATURE SCAN μ(T)
# =============================================================================
# Runs charge KMC at each temperature in temperature_scan_list.
# Fits power law mu ~ T^alpha to classify band-like vs hopping transport.
# =============================================================================

def run_temperature_scan(system, params):
    """Compute mobility at multiple temperatures."""
    print("\n" + "="*60)
    print("  MODULE 5G: Temperature scan μ(T)...")
    print("="*60)

    T_list  = params.get('temperature_scan_list', [200, 250, 300, 350, 400])
    results = []
    for T in T_list:
        p = {**params, 'temperature': T, 'num_runs': max(3, params['num_runs'] // 2)}
        if p['charge_rate_model'] == 'mlj':
            rm = calculate_mlj_rates(system, p)
        else:
            rm = calculate_marcus_rates(system, p)
        run_res = []
        for ri in range(p['num_runs']):
            run_res.append(run_kmc(system, rm, p, run_index=ri))
        tr = calculate_transport_properties(run_res, p, system)
        mu = tr['mobility_mean']
        results.append({'T': T, 'mu': mu})
        print(f"  T = {T:4d} K  |  mu = {mu:.4e} cm^2/(V.s)")

    T_arr  = np.array([r['T']  for r in results], dtype=float)
    mu_arr = np.array([r['mu'] for r in results], dtype=float)
    valid  = mu_arr > 0
    if valid.sum() > 1:
        alpha = np.polyfit(np.log(T_arr[valid]), np.log(mu_arr[valid]), 1)[0]
        regime = 'BAND-LIKE (alpha < 0)' if alpha < 0 else 'HOPPING (alpha > 0)'
        print(f"\n  [OK] Power law exponent alpha = {alpha:.2f}  →  {regime}")
    else:
        alpha = 0.0
        regime = 'undetermined'

    return {'temperatures': T_arr.tolist(), 'mobilities': mu_arr.tolist(),
            'alpha': alpha, 'regime': regime}


# =============================================================================
# MODULE 5H: TRAP-LIMITED TRANSPORT (exponential DOS tail / Urbach tail)
# =============================================================================
# Adds exponential trap states below the Gaussian DOS.
# Trap sites have deeper site energies drawn from exponential distribution.
# DOI: 10.1103/PhysRevB.63.085202 (Vissenberg & Matters trap model)
# =============================================================================

def apply_trap_dos(system, params):
    """
    Modify a fraction of site energies to follow an exponential trap DOS.
    Traps drawn from:  E_trap = E_bandedge - Exponential(E_U)
    Returns updated system dict with trap-modified site energies.
    """
    if not params.get('enable_trap_dos', False):
        return system

    print("\n  MODULE 5H: Applying exponential trap DOS (Urbach tail)...")

    rng          = np.random.default_rng(params['random_seed'] + 44444)
    total_sites  = system['total_sites']
    f_trap       = params.get('trap_fraction', 0.05)
    E_U          = params.get('urbach_energy', 0.03)
    trap_depth   = params.get('trap_depth_mean', 0.15)
    n_traps      = max(1, int(f_trap * total_sites))

    trap_indices = rng.choice(total_sites, size=n_traps, replace=False)
    site_energies = system['site_energies'].copy()

    for idx in trap_indices:
        # Trap energy: pulled below band edge by exponential distribution
        depth = rng.exponential(E_U) + trap_depth
        site_energies[idx] -= depth

    is_trap = np.zeros(total_sites, dtype=bool)
    is_trap[trap_indices] = True

    system = {**system, 'site_energies': site_energies, 'is_trap': is_trap}
    print(f"  [OK] {n_traps} trap sites added  |  E_U = {E_U:.3f} eV  |  "
          f"mean depth = {trap_depth:.3f} eV")
    return system


# =============================================================================
# MODULE 6C: SINGLET EXCITON KMC (Förster FRET)
# =============================================================================
# Simulates singlet exciton diffusion via Förster energy transfer.
# k_FRET = k0 * (R0/r)^6 * exp(-max(dE,0)/kBT)
# Much faster decay than triplets (k_S ~ 10^9 s^-1 vs k_T ~ 10^5 s^-1)
# Extracts singlet diffusion length L_S.
# =============================================================================

def run_singlet_kmc(system, params):
    """Simulate singlet exciton diffusion via Förster transfer."""
    print("\n" + "="*60)
    print("  MODULE 6C: Singlet Exciton KMC (Förster FRET)...")
    print("="*60)

    k0         = params.get('k_fret_prefactor', 1e12)
    R0_Ang     = params.get('forster_radius', 3.0) * 10   # nm -> Ang
    k_s_decay  = params.get('k_singlet_decay', 1e9)
    n_exc      = params.get('singlet_num_excitons', 100)
    max_steps  = params.get('singlet_max_steps', int(1e5))
    T          = params['temperature']
    kBT        = 8.617e-5 * T
    total_sites= system['total_sites']
    positions  = system['positions']
    neighbors  = system['neighbors']
    site_energies = system['site_energies']
    rng        = np.random.default_rng(params['random_seed'] + 66666)

    # Build Förster rate matrix
    print("  Building Förster rate matrix...")
    forster_matrix = [[] for _ in range(total_sites)]
    for i in range(total_sites):
        for (_, j, _) in neighbors[i]:
            r_ij = np.linalg.norm(positions[i] - positions[j])
            if r_ij < 1e-6:
                continue
            dE = site_energies[j] - site_energies[i]
            rate = k0 * (R0_Ang / r_ij)**6 * np.exp(-max(dE, 0) / kBT)
            forster_matrix[i].append((rate, j))

    D_estimates = []
    lifetime_list = []
    print(f"  Simulating {n_exc} singlet excitons...")

    for _ in range(n_exc):
        site      = rng.integers(0, total_sites)
        start_pos = positions[site].copy()
        sim_time  = 0.0

        for step in range(max_steps):
            evts  = [(r, d) for r, d in forster_matrix[site]]
            evts.append((k_s_decay, -1))
            rates = np.array([e[0] for e in evts])
            R     = rates.sum()
            if R <= 0:
                break
            u1  = rng.uniform(0, 1)
            idx = min(np.searchsorted(np.cumsum(rates), u1 * R), len(evts) - 1)
            u2  = max(rng.uniform(0, 1), 1e-300)
            sim_time += -np.log(u2) / R
            if evts[idx][1] == -1:
                break
            site = evts[idx][1]

        disp_sq = np.sum((positions[site] - start_pos)**2)
        lifetime_list.append(sim_time)
        if sim_time > 0:
            D_estimates.append(disp_sq * 1e-16 / (6 * sim_time))

    D_S      = np.mean(D_estimates) if D_estimates else 0.0
    tau_S    = np.mean(lifetime_list)
    L_S_Ang  = np.sqrt(6 * D_S * tau_S) * 1e8 if D_S > 0 else 0.0

    print(f"  [OK] D_S (singlet diff.)  : {D_S:.4e} cm^2/s")
    print(f"  [OK] tau_S (lifetime)     : {tau_S:.4e} s")
    print(f"  [OK] L_S (diff. length)   : {L_S_Ang:.1f} Ang")
    return {'D_S': D_S, 'tau_S': tau_S, 'L_S_Ang': L_S_Ang}


# =============================================================================
# MODULE 6D: TRIPLET-POLARON QUENCHING (TPQ)
# =============================================================================
# Triplet exciton + trapped charge carrier -> quenched (dark state).
# DOI: 10.1103/PhysRevLett.105.166601
# Polaron sites are pre-assigned; triplet quenched when visiting polaron site.
# =============================================================================

def run_tpq_simulation(system, dexter_matrix, params):
    """Simulate triplet diffusion with triplet-polaron quenching."""
    print("\n" + "="*60)
    print("  MODULE 6D: Triplet-Polaron Quenching (TPQ)...")
    print("="*60)

    k_tpq      = params.get('k_tpq', 1e8)
    f_pol      = params.get('polaron_fraction', 0.01)
    k_decay    = params['k_decay']
    n_exc      = params.get('diffusion_num_excitons', 100)
    max_steps  = params.get('diffusion_max_steps', int(1e5))
    total_sites= system['total_sites']
    positions  = system['positions']
    rng        = np.random.default_rng(params['random_seed'] + 77888)

    # Assign polaron sites
    n_pol     = max(1, int(f_pol * total_sites))
    pol_sites = set(rng.choice(total_sites, n_pol, replace=False).tolist())
    is_polaron = np.zeros(total_sites, dtype=bool)
    for s in pol_sites:
        is_polaron[s] = True

    n_quenched = 0
    D_estimates = []

    for _ in range(n_exc):
        site      = rng.integers(0, total_sites)
        start_pos = positions[site].copy()
        sim_time  = 0.0
        quenched  = False

        for step in range(max_steps):
            evts = [(r, d) for r, d in dexter_matrix[site]]
            evts.append((k_decay, -1))
            if is_polaron[site]:
                evts.append((k_tpq, -2))   # -2 = TPQ event
            rates = np.array([e[0] for e in evts])
            R     = rates.sum()
            if R <= 0:
                break
            u1  = rng.uniform(0, 1)
            idx = min(np.searchsorted(np.cumsum(rates), u1 * R), len(evts) - 1)
            u2  = max(rng.uniform(0, 1), 1e-300)
            sim_time += -np.log(u2) / R
            dest = evts[idx][1]
            if dest == -1:
                break
            elif dest == -2:
                quenched = True
                n_quenched += 1
                break
            site = dest

        disp_sq = np.sum((positions[site] - start_pos)**2)
        if sim_time > 0:
            D_estimates.append(disp_sq * 1e-16 / (6 * sim_time))

    D_T_tpq = np.mean(D_estimates) if D_estimates else 0.0
    print(f"  [OK] Polaron sites        : {n_pol} ({100*f_pol:.1f}%)")
    print(f"  [OK] TPQ quenched         : {n_quenched}/{n_exc} ({100*n_quenched/n_exc:.1f}%)")
    print(f"  [OK] D_T with TPQ         : {D_T_tpq:.4e} cm^2/s")
    return {'D_T_tpq': D_T_tpq, 'tpq_fraction': n_quenched / n_exc,
            'n_polaron_sites': n_pol}


# =============================================================================
# MODULE 7B: VORONOI NEIGHBOUR CONSTRUCTION
# =============================================================================
# Uses scipy.spatial.Voronoi to find natural neighbours of each site.
# More physical than a fixed cutoff for amorphous/melt morphologies.
# Falls back to cutoff method if scipy unavailable.
# =============================================================================

def build_voronoi_neighbors(system, params):
    """Replace cutoff neighbour list with Voronoi neighbours."""
    print("\n  MODULE 7B: Building Voronoi neighbour list...")
    try:
        from scipy.spatial import Voronoi
    except ImportError:
        print("  [WARN] scipy not available — falling back to cutoff neighbours")
        return system

    positions    = system['positions']
    total_sites  = system['total_sites']
    site_layer   = system['site_layer']
    intra_d      = params['intra_layer_distance']
    inter_d      = params['inter_layer_distance']

    vor = Voronoi(positions)
    neighbors_vor = [[] for _ in range(total_sites)]
    ridge_set     = set()

    for ridge in vor.ridge_points:
        i, j = int(ridge[0]), int(ridge[1])
        if (i, j) in ridge_set or (j, i) in ridge_set:
            continue
        ridge_set.add((i, j))
        r_ij     = np.linalg.norm(positions[i] - positions[j])
        hop_type = 'intra' if site_layer[i] == site_layer[j] else 'inter'
        neighbors_vor[i].append((0.0, j, hop_type))
        neighbors_vor[j].append((0.0, i, hop_type))

    system = {**system, 'neighbors': neighbors_vor}
    n_avg  = np.mean([len(neighbors_vor[i]) for i in range(total_sites)])
    print(f"  [OK] Voronoi neighbours built  |  avg neighbours/site = {n_avg:.1f}")
    return system


# =============================================================================
# MODULE 7C: SITE ENERGY CORRELATION FUNCTION
# =============================================================================
# Computes C(r) = <δε(0) δε(r)> / <δε²>
# Characterises the spatial correlation length of energetic disorder.
# A correlation length > lattice constant indicates correlated disorder.
# =============================================================================

def compute_energy_correlation(system, params):
    """Compute normalised spatial correlation function of site energies."""
    if not params.get('compute_energy_correlation', False):
        return None
    print("\n  MODULE 7C: Site energy correlation function C(r)...")

    positions    = system['positions']
    site_energies= system['site_energies']
    max_r        = params.get('correlation_max_r', 20.0)
    n_bins       = params.get('correlation_n_bins', 50)
    total_sites  = system['total_sites']

    de       = site_energies - site_energies.mean()
    var      = np.mean(de**2)
    if var < 1e-20:
        print("  [WARN] Zero variance — correlation undefined")
        return None

    bin_edges  = np.linspace(0, max_r, n_bins + 1)
    bin_center = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    C_sum      = np.zeros(n_bins)
    C_count    = np.zeros(n_bins)

    # Sample random pairs for speed (avoid O(N^2) for large systems)
    rng      = np.random.default_rng(params['random_seed'] + 11111)
    n_sample = min(total_sites * 20, 50000)
    idx_i    = rng.integers(0, total_sites, n_sample)
    idx_j    = rng.integers(0, total_sites, n_sample)

    for i, j in zip(idx_i, idx_j):
        r_ij = np.linalg.norm(positions[i] - positions[j])
        if r_ij > max_r:
            continue
        b = np.searchsorted(bin_edges, r_ij) - 1
        if 0 <= b < n_bins:
            C_sum[b]   += de[i] * de[j]
            C_count[b] += 1

    C_r = np.where(C_count > 0, C_sum / (C_count * var), 0.0)

    # Estimate correlation length from exponential fit
    valid = (bin_center > 1.0) & (C_r > 0)
    if valid.sum() > 3:
        try:
            from scipy.optimize import curve_fit
            def exp_decay(r, xi): return np.exp(-r / xi)
            popt, _ = curve_fit(exp_decay, bin_center[valid], C_r[valid],
                                p0=[3.0], maxfev=1000)
            xi = popt[0]
        except Exception:
            xi = 0.0
    else:
        xi = 0.0

    print(f"  [OK] Correlation computed  |  correlation length xi = {xi:.2f} Ang")
    return {'r': bin_center, 'C_r': C_r, 'xi': xi}


# =============================================================================
# MODULE 7D: TRAP IDENTIFICATION — GRAPH LAPLACIAN METHOD
# =============================================================================
# Identifies trap sites without running full KMC dynamics.
# Method: random-walk Laplacian L = D^{-1} A (row-normalised transition matrix)
# Small eigenvalues of L -> slow modes -> trap clusters
# Based on: Chen, van der Hoorn, Baumeier, arXiv 2024
# =============================================================================

def identify_traps(system, rate_matrix, params):
    """Identify trap sites using random-walk Laplacian eigenanalysis."""
    if not params.get('enable_trap_identification', False):
        return None
    print("\n  MODULE 7D: Trap identification (graph Laplacian)...")

    try:
        from scipy.sparse import csr_matrix
        from scipy.sparse.linalg import eigs
    except ImportError:
        print("  [WARN] scipy.sparse not available — skipping trap identification")
        return None

    total_sites = system['total_sites']
    n_eig       = params.get('trap_id_n_eigenvalues', 10)

    # Build transition matrix (row-normalised)
    rows, cols, vals = [], [], []
    for i in range(total_sites):
        row_rates = [r for r, _ in rate_matrix[i]] if rate_matrix[i] else []
        R_i = sum(row_rates)
        if R_i <= 0:
            continue
        for r, j in rate_matrix[i]:
            rows.append(i)
            cols.append(j)
            vals.append(r / R_i)

    T_mat = csr_matrix((vals, (rows, cols)), shape=(total_sites, total_sites))
    # Laplacian: L = I - T
    from scipy.sparse import eye
    L = eye(total_sites, format='csr') - T_mat

    try:
        eigenvalues, eigenvectors = eigs(L, k=min(n_eig, total_sites - 2),
                                         which='SM', tol=1e-6)
        eigenvalues = np.real(eigenvalues)
        eigenvectors = np.real(eigenvectors)
        # Sites with large components in smallest eigenvectors are traps
        trap_score = np.sum(eigenvectors[:, :3]**2, axis=1)
        trap_threshold = np.percentile(trap_score, 95)
        is_trap_graph  = trap_score > trap_threshold
        n_traps_found  = np.sum(is_trap_graph)
        print(f"  [OK] Eigenvalues computed  |  {n_traps_found} trap sites identified")
        return {'eigenvalues': eigenvalues, 'trap_score': trap_score,
                'is_trap': is_trap_graph, 'n_traps': int(n_traps_found)}
    except Exception as e:
        print(f"  [WARN] Eigenvalue computation failed: {e}")
        return None


# =============================================================================
# MODULE 7E: EXCITON COUPLING — TrESP APPROXIMATION
# =============================================================================
# Transition charge from electrostatic potentials (TrESP).
# Approximates BSE exciton coupling using atomic transition charges.
# J_TrESP = sum_{a in i} sum_{b in j} q_a * q_b / (4pi eps0 * r_ab)
# Faster than full GW-BSE; good for screening many molecular pairs.
# =============================================================================

def compute_tresp_couplings(system, params):
    """Compute exciton coupling matrix using TrESP approximation."""
    if not params.get('enable_tresp_coupling', False):
        return None
    print("\n  MODULE 7E: Exciton coupling via TrESP...")

    tresp_file = params.get('tresp_charges_file', 'tresp_charges.txt')
    try:
        data = np.loadtxt(tresp_file)
        # Expected format: atom_idx  mol_idx  q_trans  x  y  z
        mol_ids    = data[:, 1].astype(int)
        q_trans    = data[:, 2]
        atom_pos   = data[:, 3:6]   # Ang
    except Exception as e:
        print(f"  [WARN] Could not read TrESP charges from {tresp_file}: {e}")
        print("  [INFO] TrESP file format: atom_idx  mol_idx  q_trans  x  y  z")
        return None

    ke_eV_Ang = 14.3996   # eV*Ang/e^2 (Coulomb constant in convenient units)
    unique_mols = np.unique(mol_ids)
    n_mols      = len(unique_mols)
    J_matrix    = np.zeros((n_mols, n_mols))

    for mi, mol_i in enumerate(unique_mols):
        mask_i = mol_ids == mol_i
        qi     = q_trans[mask_i]
        ri     = atom_pos[mask_i]
        for mj, mol_j in enumerate(unique_mols):
            if mj <= mi:
                continue
            mask_j = mol_ids == mol_j
            qj     = q_trans[mask_j]
            rj     = atom_pos[mask_j]
            J_ij   = 0.0
            for a, (qa, ra) in enumerate(zip(qi, ri)):
                for b, (qb, rb) in enumerate(zip(qj, rj)):
                    r_ab = np.linalg.norm(ra - rb)
                    if r_ab > 1e-6:
                        J_ij += ke_eV_Ang * qa * qb / r_ab
            J_matrix[mi, mj] = J_ij
            J_matrix[mj, mi] = J_ij

    J_mean = np.abs(J_matrix[J_matrix != 0]).mean() if (J_matrix != 0).any() else 0
    print(f"  [OK] TrESP coupling computed for {n_mols} molecules")
    print(f"  [OK] Mean |J_TrESP| = {J_mean:.4f} eV")
    return {'J_matrix': J_matrix, 'J_mean': J_mean, 'n_molecules': n_mols}


# =============================================================================
# MODULE 8B: 3D ANISOTROPIC DIFFUSION TENSOR
# =============================================================================
# Computes full 3x3 diffusion tensor D_ab = <dr_a * dr_b> / (2 * t)
# where a,b ∈ {x, y, z}. Reveals anisotropy of exciton transport.
# Also computes principal axes and eigenvalues via diagonalisation.
# =============================================================================

def compute_diffusion_tensor(system, dexter_matrix, params):
    """Compute full 3D anisotropic triplet diffusion tensor."""
    if not params.get('compute_diffusion_tensor', False):
        return None
    print("\n  MODULE 8B: 3D Anisotropic Diffusion Tensor...")

    n_exc     = params.get('diffusion_tensor_excitons', 200)
    max_steps = params.get('diffusion_max_steps', int(1e5))
    k_decay   = params['k_decay']
    total_s   = system['total_sites']
    positions = system['positions']
    rng       = np.random.default_rng(params['random_seed'] + 99999)

    D_outer = np.zeros((3, 3))   # accumulates dr_a * dr_b
    n_valid = 0

    for _ in range(n_exc):
        site      = rng.integers(0, total_s)
        start_pos = positions[site].copy()
        sim_time  = 0.0

        for step in range(max_steps):
            evts  = [(r, d) for r, d in dexter_matrix[site]]
            evts.append((k_decay, -1))
            rates = np.array([e[0] for e in evts])
            R     = rates.sum()
            if R <= 0:
                break
            u1  = rng.uniform(0, 1)
            idx = min(np.searchsorted(np.cumsum(rates), u1 * R), len(evts) - 1)
            u2  = max(rng.uniform(0, 1), 1e-300)
            sim_time += -np.log(u2) / R
            if evts[idx][1] == -1:
                break
            site = evts[idx][1]

        if sim_time > 0:
            dr = (positions[site] - start_pos) * 1e-8   # Ang -> cm
            D_outer += np.outer(dr, dr) / (2 * sim_time)
            n_valid += 1

    if n_valid > 0:
        D_tensor = D_outer / n_valid
    else:
        D_tensor = np.zeros((3, 3))

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(D_tensor)
    idx_sort   = np.argsort(eigenvalues)[::-1]
    eigenvalues   = eigenvalues[idx_sort]
    eigenvectors  = eigenvectors[:, idx_sort]
    D_iso         = np.trace(D_tensor) / 3.0

    print(f"  [OK] D_xx = {D_tensor[0,0]:.4e}  D_yy = {D_tensor[1,1]:.4e}  D_zz = {D_tensor[2,2]:.4e}")
    print(f"  [OK] D_iso (trace/3) = {D_iso:.4e} cm^2/s")
    print(f"  [OK] Principal values: {eigenvalues[0]:.4e}, {eigenvalues[1]:.4e}, {eigenvalues[2]:.4e}")
    return {'D_tensor': D_tensor, 'eigenvalues': eigenvalues,
            'eigenvectors': eigenvectors, 'D_iso': D_iso}


# =============================================================================
# MODULE 8C: HDF5 STATE FILE OUTPUT
# =============================================================================
# Saves complete simulation state (parameters, system, results) to HDF5.
# Enables post-processing, restarts, and comparison across runs.
# Requires: h5py (pip install h5py)
# =============================================================================



def compute_seebeck(all_charge_results, system, params):
    """Compute Seebeck coefficient from carrier site energy statistics."""
    if not params.get('compute_seebeck', False):
        return None
    print("\n  MODULE 8D: Seebeck Coefficient...")

    site_energies = system['site_energies']
    T             = params['temperature']
    kB            = 8.617e-5   # eV/K

    mean_E_list = []
    for result in all_charge_results:
        # run_kmc stores site indices in 'site' array
        traj = result.get('site', result.get('trajectory', []))
        if traj is None or len(traj) == 0:
            continue
        energies = [site_energies[int(s)] for s in traj if 0 <= int(s) < len(site_energies)]
        if energies:
            mean_E_list.append(np.mean(energies))

    if not mean_E_list:
        print("  [WARN] No trajectory site data available for Seebeck calculation")
        return None

    E_mean = np.mean(mean_E_list)
    E_mean_eV = E_mean - np.mean(site_energies)   # relative to band centre
    S_uV_K    = -E_mean_eV / (T) * 1e6            # eV/K -> uV/K (q=1 electron)

    print(f"  [OK] <E_carrier> relative = {E_mean_eV*1000:.2f} meV")
    print(f"  [OK] Seebeck S            = {S_uV_K:.1f} uV/K")
    return {'S_uV_K': S_uV_K, 'E_carrier_mean_eV': E_mean_eV}


# =============================================================================
# MODULE 8E: EXCITATION DENSITY DEPENDENCE Phi_UC(I)
# =============================================================================
# Simulates TTA at multiple pair densities (proxy for excitation intensity).
# Produces threshold curve: Phi_UC vs I (relative intensity).
# At low density: Phi_UC ~ I (linear).
# At high density: Phi_UC ~ sqrt(I) (TTA-limited, quadratic -> sub-linear).
# =============================================================================



def run_intensity_scan(system, dexter_matrix, params, is_quencher=None):
    """Compute Phi_UC at multiple triplet pair numbers (excitation densities)."""
    if not params.get('enable_intensity_scan', False):
        return None
    print("\n  MODULE 8E: Excitation Density Dependence Phi_UC(I)...")

    pairs_list = params.get('intensity_scan_pairs', [20, 50, 100, 200, 500])
    results    = []

    for n_pairs in pairs_list:
        p = {**params, 'num_triplet_pairs': n_pairs}
        raw  = run_TTA_simulation(system, dexter_matrix, p, is_quencher=is_quencher)
        stat = analyze_TTA(raw, p)
        phi  = stat['phi_UC']
        pta  = stat['P_TTA']
        results.append({'n_pairs': n_pairs, 'P_TTA': pta, 'phi_UC': phi})
        print(f"  n_pairs = {n_pairs:5d}  |  P(TTA) = {pta:.4f}  |  Phi_UC = {phi:.5f}")

    # Estimate threshold: where Phi_UC starts to saturate
    phi_arr   = np.array([r['phi_UC']  for r in results])
    n_arr     = np.array([r['n_pairs'] for r in results], dtype=float)
    phi_max   = phi_arr.max()
    threshold = n_arr[np.argmax(phi_arr >= 0.5 * phi_max)] if phi_max > 0 else 0

    print(f"\n  [OK] Max Phi_UC = {phi_max:.5f}  |  50% threshold ~ {threshold:.0f} pairs")
    return {'scan': results, 'phi_UC_max': phi_max, 'threshold_pairs': threshold}


# =============================================================================
# MODULE 8F: COARSE-GRAINING (Iterative Boltzmann Inversion — IBI)
# =============================================================================
# Groups molecular sites into supercells and computes effective inter-cell
# hopping rates. Produces a coarse-grained rate map for large-scale modelling.
# Based on VOTCA-CSG methodology (Reith, Pütz, Müller-Plathe 2003).
# =============================================================================



def coarse_grain_rates(system, rate_matrix, params):
    """Compute coarse-grained effective rates between supercells."""
    if not params.get('enable_coarse_graining', False):
        return None
    print("\n  MODULE 8F: Coarse-Graining (IBI-inspired supercell averaging)...")

    cg_size     = params.get('cg_supercell_size', 3)
    total_sites = system['total_sites']
    positions   = system['positions']
    site_layer  = system['site_layer']

    # Assign each site to a supercell based on position grid
    pos_min = positions.min(axis=0)
    intra_d = params['intra_layer_distance']
    cg_idx  = np.floor((positions - pos_min) / (cg_size * intra_d)).astype(int)
    # Create unique supercell label
    max_cg  = cg_idx.max(axis=0) + 1
    cell_id = (cg_idx[:, 0] * max_cg[1] * max_cg[2] +
               cg_idx[:, 1] * max_cg[2] +
               cg_idx[:, 2])

    unique_cells = np.unique(cell_id)
    n_cells      = len(unique_cells)
    cell_map     = {c: i for i, c in enumerate(unique_cells)}

    # Compute effective inter-cell rate: sum all rates crossing cell boundary
    cg_rates = {}
    for i in range(total_sites):
        ci = cell_map[cell_id[i]]
        for rate, j in rate_matrix[i]:
            cj = cell_map[cell_id[j]]
            if ci != cj:
                key = (ci, cj)
                cg_rates[key] = cg_rates.get(key, 0.0) + rate

    # Normalise by number of sites in source cell
    cell_sizes = np.bincount([cell_map[c] for c in cell_id], minlength=n_cells)
    cg_rates_norm = {k: v / max(1, cell_sizes[k[0]]) for k, v in cg_rates.items()}

    mean_cg_rate = np.mean(list(cg_rates_norm.values())) if cg_rates_norm else 0
    print(f"  [OK] Supercell size   : {cg_size} sites/edge")
    print(f"  [OK] Coarse-grained cells  : {n_cells}")
    print(f"  [OK] CG inter-cell pairs   : {len(cg_rates_norm)}")
    print(f"  [OK] Mean CG rate          : {mean_cg_rate:.4e} s^-1")
    return {'n_cells': n_cells, 'cg_rates': cg_rates_norm,
            'cell_sizes': cell_sizes, 'mean_cg_rate': mean_cg_rate}


# =============================================================================
# MODULE 9B: EXTENDED PLOTTING (new modules)
# =============================================================================



def save_hdf5_state(system, params, transport, tta_stats, diffusion_stats):
    """Save all results to HDF5 state file."""
    if not params.get('save_hdf5', False):
        return
    try:
        import h5py
    except ImportError:
        print("  [WARN] h5py not installed — skipping HDF5 output (pip install h5py)")
        return

    import os
    os.makedirs(params['output_folder'], exist_ok=True)
    fname = os.path.join(params['output_folder'], params.get('hdf5_filename', 'kmc_state.h5'))

    with h5py.File(fname, 'w') as f:
        # System geometry
        grp = f.create_group('system')
        grp.create_dataset('positions',         data=system['positions'])
        grp.create_dataset('site_energies',     data=system['site_energies'])
        grp.create_dataset('site_layer',        data=system['site_layer'])
        grp.create_dataset('site_species',      data=system['site_species'])
        grp.create_dataset('is_grain_boundary', data=system['is_grain_boundary'])

        # Parameters
        pgrp = f.create_group('parameters')
        for k, v in params.items():
            try:
                pgrp.attrs[k] = v if not isinstance(v, list) else str(v)
            except Exception:
                pass

        # Charge transport results
        if transport:
            tgrp = f.create_group('charge_transport')
            for k, v in transport.items():
                try:
                    tgrp.create_dataset(k, data=v)
                except Exception:
                    tgrp.attrs[k] = str(v)

        # TTA results
        if tta_stats:
            agrp = f.create_group('tta')
            for k, v in tta_stats.items():
                try:
                    agrp.create_dataset(k, data=v)
                except Exception:
                    agrp.attrs[k] = str(v)

        # Diffusion
        if diffusion_stats:
            dgrp = f.create_group('diffusion')
            for k, v in diffusion_stats.items():
                try:
                    dgrp.create_dataset(k, data=v)
                except Exception:
                    dgrp.attrs[k] = str(v)

    print(f"  [OK] HDF5 state saved: {fname}")


# =============================================================================
# MODULE 8D: SEEBECK COEFFICIENT (THERMOPOWER)
# =============================================================================
# S = -<E_carrier> / (q * T)
# where <E_carrier> is the mean site energy of the carrier during transport.
# Physically: a carrier in a higher energy site carries more thermal energy.
# Seebeck gives the entropy per carrier, proportional to distance from Fermi.
# =============================================================================

