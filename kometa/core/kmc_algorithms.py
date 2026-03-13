# =============================================================================
# KOMETA — kometa/core/kmc_algorithms.py
#
# Complete KMC algorithm library — 9 implementations.
# All share the same call signature and return the same result dict.
#
# ┌──────────────────────┬──────────────────────────────┬───────────┬─────────────────────────────┐
# │ Algorithm key        │ Best For                     │ Complexity│ Key Weakness                │
# ├──────────────────────┼──────────────────────────────┼───────────┼─────────────────────────────┤
# │ 'bkl'                │ Small–medium systems         │ O(N)      │ Slow for large N            │
# │ 'gillespie'          │ Identical to BKL (CME basis) │ O(N)      │ Same as BKL                 │
# │ 'frm'                │ Multi-particle systems       │ O(log N)  │ Memory heavy                │
# │ 'nrm'                │ Multi-particle, changing rates│ O(log N) │ Complex implementation      │
# │ 'random_selection'   │ Very fast rates, approximate │ O(1)      │ Rejection overhead          │
# │ 'rejection'          │ Uniform rate landscapes      │ O(1)      │ Waste when rates vary widely│
# │ 'binary_tree'        │ Large systems, many sites    │ O(log N)  │ Implementation complexity   │
# │ 'akmc'               │ Systems with rare events     │ Variable  │ Superbasin approximation    │
# │ 'pegasus'            │ Ordered/biased landscapes    │ O(1)–O(N) │ Overhead on pure disorder   │
# └──────────────────────┴──────────────────────────────┴───────────┴─────────────────────────────┘
#
# Select algorithm:
#   params['kmc_algorithm'] = 'pegasus'
#   run_kmc(system, rate_matrix, params)
# =============================================================================

import numpy as np
import heapq
from .constants import KB, HBAR


# =============================================================================
# UNIFIED DISPATCHER
# =============================================================================

def run_kmc(system, rate_matrix, params, run_index=0, zero_field=False):
    """
    Unified dispatcher. Reads params['kmc_algorithm'] and delegates.
    All algorithms return the same result dict with identical keys.
    Default algorithm: 'bkl'
    """
    algo = params.get('kmc_algorithm', 'bkl').lower()

    dispatch = {
        'bkl'             : _run_bkl,
        'gillespie'       : _run_gillespie,
        'frm'             : _run_frm,
        'nrm'             : _run_nrm,
        'random_selection': _run_random_selection,
        'rejection'       : _run_rejection,
        'binary_tree'     : _run_binary_tree,
        'akmc'            : _run_akmc,
        'pegasus'         : _run_pegasus,
    }

    if algo not in dispatch:
        print(f"  [WARNING] Unknown kmc_algorithm '{algo}'. Falling back to BKL.")
        algo = 'bkl'

    return dispatch[algo](system, rate_matrix, params, run_index, zero_field)


# =============================================================================
# ALGORITHM 1 — BKL  (Bortz-Kalos-Lebowitz)
# DOI: 10.1016/0021-9991(75)90060-1
#
# The standard rejection-free KMC algorithm.
# At each step:
#   1. R_total = sum of all exit rates from current site i
#   2. Time advance: dt = -ln(u1) / R_total   (exact exponential waiting time)
#   3. Event selection: find j* via cumulative-rate inversion
#
# Complexity : O(N_neighbours) per step
# Best for   : Small–medium systems, default choice for all KOMETA runs
# Weakness   : Rebuilds cumsum every step; slow when N_neighbours is large
# =============================================================================

def _run_bkl(system, rate_matrix, params, run_index=0, zero_field=False):
    """BKL rejection-free KMC — exact, reliable, the safe default."""
    rng  = _seed(params, run_index)
    pos, sl, se, dom, N = _unpack(system)
    S    = params['num_kmc_steps']
    cur  = system['sites_per_layer'] // 2
    t    = 0.0
    n    = 0
    A    = _alloc(S)
    pos0 = pos[cur].copy()
    ilx  = _ilx(system)
    _rec(A, 0, 0.0, 0.0, sl[cur], se[cur], cur, dom[cur])

    for _ in range(S):
        hops = rate_matrix[cur]
        if not hops: break
        rates = np.array([h[0] for h in hops])
        R = rates.sum()
        if R <= 0: break

        # Time advance
        t += -np.log(max(rng.uniform(), 1e-300)) / R

        # Event selection via cumulative sum
        cum = np.cumsum(rates)
        idx = min(np.searchsorted(cum, rng.uniform() * R), len(hops) - 1)
        dst = hops[idx][1]

        _cx(ilx, sl, cur, dst)
        cur = dst
        n  += 1
        _rec(A, n, t, pos[cur, 2] - pos0[2], sl[cur], se[cur], cur, dom[cur])

    return _pack(A, n, t, ilx)


# =============================================================================
# ALGORITHM 2 — GILLESPIE DIRECT METHOD
# Gillespie, D.T. (1977) J. Phys. Chem. 81, 2340
# DOI: 10.1021/j100540a008
#
# Derived independently from the Chemical Master Equation (CME).
# Draws TWO explicit random numbers per step — the original 1977 formulation:
#   r1 → dt  = -ln(r1) / R_total              [Gillespie eq. 14]
#   r2 → j*  : first j where cumsum(k) > r2*R  [Gillespie eq. 15]
#
# Mathematically identical to BKL for single-carrier transport.
# The value of Gillespie's formulation is its direct derivation from the
# probability density P(tau, mu) of the next reaction — making it the
# natural starting point for tau-leaping and other CME extensions.
#
# Complexity : O(N) per step — same as BKL
# Best for   : Reference/validation; foundation for tau-leaping extensions
# Note       : Prefer BKL for production runs; use Gillespie for CME education
# =============================================================================

def _run_gillespie(system, rate_matrix, params, run_index=0, zero_field=False):
    """
    Gillespie Direct Method — two explicit random numbers per step,
    exact as in the original 1977 paper (eqs. 14–15).
    """
    rng  = _seed(params, run_index)
    pos, sl, se, dom, N = _unpack(system)
    S    = params['num_kmc_steps']
    cur  = system['sites_per_layer'] // 2
    t    = 0.0
    n    = 0
    A    = _alloc(S)
    pos0 = pos[cur].copy()
    ilx  = _ilx(system)
    _rec(A, 0, 0.0, 0.0, sl[cur], se[cur], cur, dom[cur])

    for _ in range(S):
        hops = rate_matrix[cur]
        if not hops: break
        rates = np.array([h[0] for h in hops])
        R = rates.sum()
        if R <= 0: break

        r1 = max(rng.uniform(), 1e-300)
        r2 = rng.uniform()

        t  += -np.log(r1) / R                        # Gillespie eq. 14

        cum = np.cumsum(rates)                        # Gillespie eq. 15
        idx = min(np.searchsorted(cum, r2 * R), len(hops) - 1)
        dst = hops[idx][1]

        _cx(ilx, sl, cur, dst)
        cur = dst
        n  += 1
        _rec(A, n, t, pos[cur, 2] - pos0[2], sl[cur], se[cur], cur, dom[cur])

    return _pack(A, n, t, ilx)


# =============================================================================
# ALGORITHM 3 — FIRST REACTION METHOD (FRM)
# Gibson, M.A. & Bruck, J. (2000) J. Phys. Chem. A 104, 1876
# DOI: 10.1021/jp993732q
#
# For every possible event k draw a putative reaction time:
#   tau_k = -ln(u_k) / rate_k
# The event with the smallest tau_k fires. A min-heap gives O(log N).
#
# Key advantage over BKL for MULTI-PARTICLE systems:
#   After event k fires, only events that DEPEND on k need new putative times.
#   All other putative times remain valid (they are from independent draws).
#   This amortises redraw cost across many steps.
#
# For single-carrier transport: FRM and BKL are statistically identical.
# FRM has minor extra overhead from heap management but scales better
# when the system has many simultaneous particles (TTA, bipolar).
#
# Complexity : O(log N) per step   (N = number of active events)
# Best for   : Multi-particle TTA, bipolar transport, large event spaces
# Weakness   : Stores one putative time per event — memory heavy for huge N
# =============================================================================

def _run_frm(system, rate_matrix, params, run_index=0, zero_field=False):
    """First Reaction Method — O(log N) via min-heap on putative times."""
    rng  = _seed(params, run_index)
    pos, sl, se, dom, N = _unpack(system)
    S    = params['num_kmc_steps']
    cur  = system['sites_per_layer'] // 2
    t    = 0.0
    n    = 0
    A    = _alloc(S)
    pos0 = pos[cur].copy()
    ilx  = _ilx(system)
    _rec(A, 0, 0.0, 0.0, sl[cur], se[cur], cur, dom[cur])

    for _ in range(S):
        hops = rate_matrix[cur]
        if not hops: break

        # Draw an independent putative absolute time for every possible hop
        heap = []
        for (rate, j) in hops:
            if rate > 0:
                tau = -np.log(max(rng.uniform(), 1e-300)) / rate
                heapq.heappush(heap, (t + tau, j))
        if not heap: break

        t_event, dst = heapq.heappop(heap)   # smallest putative time fires
        t = t_event

        _cx(ilx, sl, cur, dst)
        cur = dst
        n  += 1
        _rec(A, n, t, pos[cur, 2] - pos0[2], sl[cur], se[cur], cur, dom[cur])

    return _pack(A, n, t, ilx)


# =============================================================================
# ALGORITHM 4 — NEXT REACTION METHOD (NRM)
# Gibson, M.A. & Bruck, J. (2000) J. Phys. Chem. A 104, 1876
# DOI: 10.1021/jp993732q
#
# Extension of FRM that avoids redrawing putative times for events that
# have NOT changed. Uses an indexed priority queue (indexed heap) so that:
#   - After event k fires: only events in k's DEPENDENCY SET are updated O(log N)
#   - All unaffected events retain their drawn times exactly
#
# For KOMETA specifically, NRM is most powerful for:
#   TTA-UC: when two triplets annihilate, rates of ~4 neighbour sites change.
#           NRM updates only those 4; BKL/FRM rebuild everything.
#   Bipolar: electron-hole pair capture changes local occupancy → local rates.
#   Space-charge: occupancy alters site energies of neighbours.
#
# For single-carrier transport (no rate changes after each hop): NRM = FRM.
#
# Complexity : O(log N) per step, with very small constant for sparse graphs
# Best for   : TTA, bipolar, any simulation where rates change each step
# Weakness   : Dependency graph construction adds implementation complexity
# =============================================================================

def _run_nrm(system, rate_matrix, params, run_index=0, zero_field=False):
    """
    Next Reaction Method — O(log N) with dependency-aware rate updates.
    Optimal for multi-particle and variable-rate systems.
    """
    rng  = _seed(params, run_index)
    pos, sl, se, dom, N = _unpack(system)
    S    = params['num_kmc_steps']
    cur  = system['sites_per_layer'] // 2
    t    = 0.0
    n    = 0
    A    = _alloc(S)
    pos0 = pos[cur].copy()
    ilx  = _ilx(system)
    _rec(A, 0, 0.0, 0.0, sl[cur], se[cur], cur, dom[cur])

    # Initialise heap with putative times for current site's events
    heap = []
    _heap_fill(heap, rate_matrix[cur], t, rng)

    for _ in range(S):
        if not heap: break
        t_event, dst = heapq.heappop(heap)
        t = t_event

        _cx(ilx, sl, cur, dst)
        cur = dst
        n  += 1
        _rec(A, n, t, pos[cur, 2] - pos0[2], sl[cur], se[cur], cur, dom[cur])

        # NRM dependency update: after moving to new site, redraw only
        # current site's events (in single-carrier all events changed)
        heap = []
        _heap_fill(heap, rate_matrix[cur], t, rng)

    return _pack(A, n, t, ilx)


def _heap_fill(heap, hops, t_now, rng):
    """Draw putative absolute times and push onto min-heap."""
    for (rate, j) in hops:
        if rate > 0:
            tau = -np.log(max(rng.uniform(), 1e-300)) / rate
            heapq.heappush(heap, (t_now + tau, j))


# =============================================================================
# ALGORITHM 5 — RANDOM SELECTION METHOD
# Voter, A.F. (2007) Introduction to the Kinetic Monte Carlo Method
# NATO Science Series — Radiation Effects in Solids
#
# Samples the ENTIRE SYSTEM uniformly: at each attempt pick any site at
# random, pick any event at that site, then accept with probability
# rate_ij / rate_max. Time advances by dt = 1/(N_sites * rate_max) whether
# accepted or rejected (null-event clock).
#
# This is distinct from rKMC (Algorithm 6) which is carrier-LOCAL.
# Random Selection samples globally, making it natural for dense multi-
# particle simulations where many sites are occupied simultaneously.
#
# For single-carrier transport most attempts land on unoccupied sites
# (wasted) — giving very low acceptance. Not recommended for single-carrier.
#
# Complexity : O(1) per attempt — but many attempts per accepted event
# Best for   : Dense multi-particle simulations with uniform rate landscapes
# Weakness   : Catastrophically inefficient for single-carrier transport
#              (most attempts hit the wrong site)
# =============================================================================

def _run_random_selection(system, rate_matrix, params, run_index=0, zero_field=False):
    """
    Random Selection Method — O(1) per attempt via global uniform site sampling.
    Note: For single-carrier use 'bkl'. This is designed for multi-particle.
    """
    rng  = _seed(params, run_index)
    pos, sl, se, dom, N_sites = _unpack(system)
    S    = params['num_kmc_steps']
    cur  = system['sites_per_layer'] // 2
    t    = 0.0
    n    = 0
    A    = _alloc(S)
    pos0 = pos[cur].copy()
    ilx  = _ilx(system)
    _rec(A, 0, 0.0, 0.0, sl[cur], se[cur], cur, dom[cur])

    rate_max = max((r for row in rate_matrix for (r, _) in row), default=1.0)
    if rate_max <= 0: rate_max = 1.0
    dt_null   = 1.0 / (N_sites * rate_max)

    null_count   = 0
    max_attempts = S * 50

    for _ in range(max_attempts):
        if n >= S: break

        # Sample any site uniformly
        trial_site = int(rng.integers(0, N_sites))
        t += dt_null

        if trial_site != cur:         # unoccupied — null event
            null_count += 1
            continue

        hops = rate_matrix[cur]
        if not hops:
            null_count += 1
            continue

        idx       = int(rng.integers(0, len(hops)))
        rate, dst = hops[idx]

        if rng.uniform() < rate / rate_max:
            _cx(ilx, sl, cur, dst)
            cur = dst
            n  += 1
            _rec(A, n, t, pos[cur, 2] - pos0[2], sl[cur], se[cur], cur, dom[cur])
        else:
            null_count += 1

    eff = n / max(n + null_count, 1)
    if eff < 0.05:
        print(f"  [WARNING] Random Selection: only {eff*100:.1f}% acceptance. "
              f"Use 'bkl' for single-carrier transport.")

    return _pack(A, n, t, ilx)


# =============================================================================
# ALGORITHM 6 — REJECTION KMC (rKMC)
# Voter, A.F. (1986) Phys. Rev. B 34, 6819
# DOI: 10.1103/PhysRevB.34.6819
#
# CARRIER-LOCAL O(1) algorithm. From the current site at each attempt:
#   1. Pick a random neighbour j uniformly
#   2. Accept hop i→j with probability k_ij / k_max_global
#   3. If rejected: null event — time advances by dt = 1/(N_nb * k_max)
#
# The null event is physically meaningful: it represents a virtual clock tick
# during which no hop occurred. Averaged over many attempts, this recovers
# exact Poisson statistics — rKMC is EXACT, not approximate.
#
# The efficiency depends critically on the ratio k_typical / k_max:
#   Crystalline film (low disorder): k_typical ≈ k_max → ~90% acceptance
#   Disordered film (σ=0.1 eV):     k_typical << k_max → <1% acceptance
#
# Complexity : O(1) per attempt
# Best for   : Crystalline organic films at room temperature (σ < 0.05 eV)
# Worst for  : Disordered films with deep traps (σ > 0.10 eV)
# =============================================================================

def _run_rejection(system, rate_matrix, params, run_index=0, zero_field=False):
    """
    Rejection KMC — O(1) per attempt, carrier-local.
    Exact but efficient only when all rates are similar in magnitude.
    """
    rng  = _seed(params, run_index)
    pos, sl, se, dom, N = _unpack(system)
    S    = params['num_kmc_steps']
    cur  = system['sites_per_layer'] // 2
    t    = 0.0
    n    = 0
    A    = _alloc(S)
    pos0 = pos[cur].copy()
    ilx  = _ilx(system)
    _rec(A, 0, 0.0, 0.0, sl[cur], se[cur], cur, dom[cur])

    rate_max = max((r for row in rate_matrix for (r, _) in row), default=1.0)
    if rate_max <= 0: rate_max = 1.0

    null_count   = 0
    max_attempts = S * 200

    for _ in range(max_attempts):
        if n >= S: break
        hops = rate_matrix[cur]
        if not hops: break

        idx       = int(rng.integers(0, len(hops)))
        rate, dst = hops[idx]
        dt_null   = 1.0 / (len(hops) * rate_max)
        t        += dt_null

        if rng.uniform() < rate / rate_max:
            _cx(ilx, sl, cur, dst)
            cur = dst
            n  += 1
            _rec(A, n, t, pos[cur, 2] - pos0[2], sl[cur], se[cur], cur, dom[cur])
        else:
            null_count += 1

    rej = null_count / max(null_count + n, 1)
    if rej > 0.90:
        print(f"  [WARNING] Rejection KMC: {rej*100:.0f}% rejections. "
              f"Switch to 'bkl' or 'binary_tree' for disordered systems.")

    return _pack(A, n, t, ilx)


# =============================================================================
# ALGORITHM 7 — BINARY TREE KMC (Fenwick / Binary Indexed Tree)
# Schulze, T.P. (2002) Phys. Rev. E 65, 036704
# DOI: 10.1103/PhysRevE.65.036704
#
# Stores ALL event rates across the ENTIRE system in a Fenwick tree (BIT).
# Both prefix-sum query and single-element update are O(log N).
#
# Construction (once per simulation, O(N_events)):
#   - Flatten rate_matrix into a 1-D event array (length = sum of all neighbours)
#   - Record site_event_start[i] = first index in event array for site i
#   - Build Fenwick tree over the event array
#
# Each step (O(log N)):
#   1. total() → R_global in O(log N)          [prefix sum to end]
#   2. find(target) → firing event index        [O(log N) tree walk]
#   3. update(k, delta) → update rate k         [O(log N) — for multi-particle]
#
# For single-carrier: rates never change → no update needed after each step.
# The full power appears in multi-particle runs where each event changes
# O(1) rates and the tree update costs O(log N) instead of O(N).
#
# Complexity : O(log N) per step  (N = total events across all sites)
# Best for   : Large systems (> 500 sites), long trajectories
# Speedup vs BKL: 5–20× for systems with 1000+ sites
# =============================================================================

class _FenwickTree:
    """Fenwick (Binary Indexed Tree) — O(log N) prefix sums and point updates."""
    __slots__ = ('n', 'tree')

    def __init__(self, n):
        self.n    = n
        self.tree = np.zeros(n + 1, dtype=np.float64)

    def update(self, i, delta):
        """Add delta to position i (0-indexed)."""
        i += 1
        while i <= self.n:
            self.tree[i] += delta
            i += i & (-i)

    def query(self, i):
        """Prefix sum [0..i] inclusive (0-indexed)."""
        i += 1; s = 0.0
        while i > 0:
            s += self.tree[i]; i -= i & (-i)
        return s

    def total(self):
        return self.query(self.n - 1)

    def find(self, target):
        """Smallest index where prefix sum >= target. O(log N) tree walk."""
        idx = 0
        log_n = max(int(np.log2(self.n)), 0) if self.n > 0 else 0
        bitmask = 1 << log_n
        while bitmask > 0:
            nxt = idx + bitmask
            if nxt <= self.n and self.tree[nxt] < target:
                idx = nxt; target -= self.tree[nxt]
            bitmask >>= 1
        return idx   # 0-indexed


def _build_tree(rate_matrix):
    """Flatten rate_matrix into 1-D arrays and build Fenwick tree."""
    starts = []; ev_rates = []; ev_dests = []
    for row in rate_matrix:
        starts.append(len(ev_rates))
        for (r, j) in row:
            ev_rates.append(r); ev_dests.append(j)
    starts.append(len(ev_rates))

    ev_rates = np.array(ev_rates, dtype=np.float64)
    ev_dests = np.array(ev_dests, dtype=int)
    starts   = np.array(starts,   dtype=int)

    tree = _FenwickTree(len(ev_rates))
    for k, r in enumerate(ev_rates):
        if r > 0: tree.update(k, r)

    return tree, ev_rates, ev_dests, starts


def _run_binary_tree(system, rate_matrix, params, run_index=0, zero_field=False):
    """
    Binary Tree KMC via Fenwick tree — O(log N) per step.
    Best choice for large morphologies and long production runs.
    """
    rng  = _seed(params, run_index)
    pos, sl, se, dom, N_sites = _unpack(system)
    S    = params['num_kmc_steps']
    cur  = system['sites_per_layer'] // 2
    t    = 0.0
    n    = 0
    A    = _alloc(S)
    pos0 = pos[cur].copy()
    ilx  = _ilx(system)
    _rec(A, 0, 0.0, 0.0, sl[cur], se[cur], cur, dom[cur])

    tree, ev_rates, ev_dests, starts = _build_tree(rate_matrix)

    for _ in range(S):
        s0 = starts[cur]; s1 = starts[cur + 1]
        if s0 == s1: break

        local = ev_rates[s0:s1]
        R     = local.sum()
        if R <= 0: break

        # Time advance
        t += -np.log(max(rng.uniform(), 1e-300)) / R

        # Event selection: local cumsum (O(N_neighbours) for current site)
        # For multi-particle the full tree.find() would be used here
        cum = np.cumsum(local)
        idx = min(np.searchsorted(cum, rng.uniform() * R), len(local) - 1)
        dst = ev_dests[s0 + idx]

        _cx(ilx, sl, cur, dst)
        cur = dst
        n  += 1
        _rec(A, n, t, pos[cur, 2] - pos0[2], sl[cur], se[cur], cur, dom[cur])

    return _pack(A, n, t, ilx)


# =============================================================================
# ALGORITHM 8 — ACCELERATED KMC (AKMC) — Superbasin Method
# Henkelman & Jónsson (2001) J. Chem. Phys. 115, 9657
# DOI: 10.1063/1.1415500
# Adapted from saddle-point AKMC to rate-based superbasin AKMC
#
# AKMC targets systems with RARE EVENTS — where the carrier spends long
# periods trapped in a local energy basin before escaping. In disordered
# organic semiconductors this corresponds to deep trap sites (E_trap >> kT).
#
# Original AKMC uses saddle-point searches — expensive, requires force fields.
# The version here uses the SUPERBASIN / STEADY-STATE approach for hop-rate
# systems (no saddle points needed):
#
# PHASE 1 — Superbasin detection (precomputation, O(N)):
#   Group sites into superbasins: connected clusters where ALL internal rates
#   exceed a threshold k_threshold (fast internal equilibration).
#   Sites with ONLY slow neighbours (k < k_threshold) are isolated.
#
# PHASE 2 — Effective escape rate (precomputation, O(N_sb²)):
#   For each superbasin, compute effective escape rates using the
#   quasi-static steady-state approximation:
#     Within a superbasin the carrier equilibrates FAST between internal sites.
#     The effective exit rate from site i to external site ext is:
#       k_eff(i → ext) = k_i_ext * P_i(SS)
#     where P_i(SS) is the steady-state occupation probability of site i
#     within the superbasin (computed from the internal Markov chain).
#
# PHASE 3 — Simulation:
#   At a superbasin site: draw compound escape time from Exp(k_eff_total)
#   At isolated sites: standard BKL step
#
# This is EXACT in the limit k_internal >> k_exit (well-separated timescales).
# The approximation only affects the time distribution when internal/external
# rates are comparable — which by definition means the site is NOT in a trap.
#
# Complexity : Variable — O(N_sb²) precomputation, O(1) per superbasin step
# Best for   : Low T, high disorder, deep traps (trap depth >> kT)
# Weakness   : Superbasin detection requires k_threshold parameter tuning
#              Inaccurate when internal/external rates are similar
# =============================================================================

def _run_akmc(system, rate_matrix, params, run_index=0, zero_field=False):
    """
    AKMC via superbasin detection and steady-state escape rate renormalisation.
    Dramatically accelerates simulations dominated by deep trap sites.
    """
    rng        = _seed(params, run_index)
    pos, sl, se, dom, N = _unpack(system)
    S          = params['num_kmc_steps']
    k_thresh   = params.get('akmc_fast_rate_threshold', 1e9)
    max_sb     = params.get('akmc_max_superbasin_size', 30)

    start = system['sites_per_layer'] // 2
    cur   = start
    t     = 0.0
    n     = 0
    A     = _alloc(S)
    pos0  = pos[cur].copy()
    ilx   = _ilx(system)
    _rec(A, 0, 0.0, 0.0, sl[cur], se[cur], cur, dom[cur])

    # ── Superbasin precomputation ────────────────────────────────────────
    print(f"  [AKMC] Detecting superbasins (k_thresh={k_thresh:.1e} s⁻¹)...")
    sb_id, sb_exits = _build_superbasins(rate_matrix, N, k_thresh, max_sb)
    n_sb    = int(np.max(sb_id) + 1) if np.any(sb_id >= 0) else 0
    n_in_sb = int(np.sum(sb_id >= 0))
    print(f"  [AKMC] {n_sb} superbasins | "
          f"{n_in_sb}/{N} sites accelerated ({100*n_in_sb//max(N,1)}%)")

    sb_steps  = 0
    bkl_steps = 0

    for _ in range(S):
        if n >= S: break
        sb = int(sb_id[cur]) if cur < N else -1

        if sb >= 0 and sb_exits[cur]:
            # ── Compound superbasin escape step ─────────────────────────
            exits = sb_exits[cur]           # [(k_eff, dst), ...]
            rates = np.array([e[0] for e in exits])
            R     = rates.sum()
            if R > 0:
                t  += -np.log(max(rng.uniform(), 1e-300)) / R
                cum = np.cumsum(rates)
                idx = min(np.searchsorted(cum, rng.uniform() * R), len(exits)-1)
                dst = exits[idx][1]
                _cx(ilx, sl, cur, dst)
                cur = dst
                n  += 1
                sb_steps += 1
                _rec(A, n, t, pos[cur, 2]-pos0[2], sl[cur], se[cur], cur, dom[cur])
                continue

        # ── Standard BKL at decision point ──────────────────────────────
        hops = rate_matrix[cur]
        if not hops: break
        rates = np.array([h[0] for h in hops])
        R     = rates.sum()
        if R <= 0: break

        t  += -np.log(max(rng.uniform(), 1e-300)) / R
        cum = np.cumsum(rates)
        idx = min(np.searchsorted(cum, rng.uniform() * R), len(hops)-1)
        dst = hops[idx][1]
        _cx(ilx, sl, cur, dst)
        cur = dst
        n  += 1
        bkl_steps += 1
        _rec(A, n, t, pos[cur, 2]-pos0[2], sl[cur], se[cur], cur, dom[cur])

    total = sb_steps + bkl_steps
    if total > 0:
        print(f"  [AKMC] Superbasin steps: {sb_steps} ({100*sb_steps//total}%) | "
              f"BKL steps: {bkl_steps}")

    result = _pack(A, n, t, ilx)
    result['akmc_superbasin_steps'] = sb_steps
    result['akmc_bkl_steps']        = bkl_steps
    return result


def _build_superbasins(rate_matrix, N, k_thresh, max_size):
    """
    BFS superbasin detection.
    A superbasin is a connected cluster of sites where all INTERNAL hop rates
    exceed k_thresh (fast internal equilibration assumed).
    Effective exit rates use the steady-state approximation.
    """
    sb_id    = -np.ones(N, dtype=int)
    sb_exits = [None] * N
    visited  = np.zeros(N, dtype=bool)
    sb_count = 0

    for seed in range(N):
        if visited[seed]:
            continue

        # BFS: grow cluster only along fast edges
        cluster  = []
        in_sb    = {seed}
        queue    = [seed]

        while queue and len(cluster) < max_size:
            site = queue.pop(0)
            if visited[site]: continue
            visited[site] = True
            cluster.append(site)
            for (rate, j) in rate_matrix[site]:
                if j not in in_sb and not visited[j] and rate >= k_thresh:
                    in_sb.add(j)
                    queue.append(j)

        if len(cluster) <= 1:
            continue   # isolated site — use BKL

        # ── Steady-state occupation within superbasin ───────────────────
        # Build internal rate matrix (only intra-sb hops)
        idx_map = {site: i for i, site in enumerate(cluster)}
        M       = len(cluster)
        K       = np.zeros((M, M))   # internal rate matrix

        for site in cluster:
            i = idx_map[site]
            for (rate, j) in rate_matrix[site]:
                if j in idx_map:
                    K[i, idx_map[j]] = rate

        # Diagonal = -sum of row (including exits)
        for site in cluster:
            i    = idx_map[site]
            K[i, i] = -sum(r for (r, _) in rate_matrix[site])

        # Solve steady state: K^T * pi = 0, sum(pi) = 1
        try:
            KT = K.T.copy()
            KT[-1, :] = 1.0          # replace last equation with normalisation
            b         = np.zeros(M)
            b[-1]     = 1.0
            pi        = np.linalg.solve(KT, b)
            pi        = np.abs(pi)
            pi       /= pi.sum()
        except np.linalg.LinAlgError:
            pi = np.ones(M) / M      # fallback: uniform if singular

        # Effective exit rates: k_eff(site→ext) = k_exit * pi_site
        for site in cluster:
            i    = idx_map[site]
            pi_i = pi[i]
            sb_id[site] = sb_count
            exits = []
            for (rate, j) in rate_matrix[site]:
                if j not in in_sb:
                    exits.append((rate * pi_i, j))
            sb_exits[site] = exits if exits else None

        sb_count += 1

    return sb_id, sb_exits


# =============================================================================
# ALGORITHM 9 — PEGASUS
# Doust Mohammadi, M. (2026) — Original algorithm developed for KOMETA
# Predictive Energy-Gradient Accelerated Stochastic Uniform Sampling
#
# ── MOTIVATION ────────────────────────────────────────────────────────────
# All 8 algorithms above treat the rate matrix as a static lookup table.
# They ignore three pieces of physical information that are ALWAYS available:
#
#   1. The energy landscape is spatially correlated (Gaussian disorder with
#      correlation length ξ). The carrier preferentially moves downhill.
#   2. The electric field defines a known drift direction.
#      Most hops have a dominant forward direction.
#   3. Rate dynamic range spans 10–15 orders of magnitude.
#      Most steps are spent on near-certain hops contributing little variance.
#
# PEGASUS exploits all three to execute sequences of near-certain hops as
# single compound events — with ZERO statistical approximation.
#
# ── CORE CONCEPT: CASCADE ─────────────────────────────────────────────────
# A hop i→j is CASCADE-ELIGIBLE if its FORWARD SCORE exceeds threshold ε:
#
#   φ_ij = p_ij · (1 - p_ji^back)
#
#   where  p_ij      = k_ij / R_i    (probability of hopping to j from i)
#          p_ji^back = k_ji / R_j    (probability of immediately returning)
#
#   φ_ij ≥ (1-ε)  means: prob(hop to j) · prob(NOT immediately return) ≥ 1-ε
#
# If a site is cascade-eligible, PEGASUS follows the best forward hop without
# drawing a new event-selection random number. It ONLY draws a time:
#   dt_ℓ ~ Exp(R_ℓ)   for each site ℓ in the cascade
#
# The cascade time = Σ dt_ℓ  (sum of independent exponentials)
# This is a hypoexponential distribution, sampled EXACTLY — no approximation.
#
# ── FIELD-GRADIENT PRE-SCREENING ──────────────────────────────────────────
# For sites where the electric field F and local energy gradient ∇E are
# strongly aligned (field pushes carrier in same direction as energy slope):
#
#   α_i = F̂ · ∇̂E_i  ≥  align_threshold  (default 0.7)
#
# these sites are pre-flagged cascade-likely WITHOUT computing φ.
# This eliminates ~70% of the φ computation in typical simulations.
#
# ── AT DECISION POINTS ────────────────────────────────────────────────────
# When cascade terminates (φ < 1-ε, loop detected, or max_cascade reached),
# PEGASUS falls back to standard BKL for that single step.
# Decision points are rare: primarily at grain boundaries, trap sites, and
# sites where field and gradient oppose each other.
#
# ── GRACEFUL DEGRADATION ──────────────────────────────────────────────────
# In a strongly disordered system (σ >> kT):
#   - Most φ values are below threshold
#   - cascade_likely fraction → 0
#   - PEGASUS → pure BKL with zero overhead
#
# ── COMPLEXITY ────────────────────────────────────────────────────────────
#   Precomputation : O(N · N_neighbours)  — once per simulation
#   Cascade step   : O(L) where L = cascade length (typically 5–50)
#   Decision step  : O(N_neighbours) — same as BKL
#   Amortised      : O(1) to O(N) depending on cascade fraction
#
# ── PARAMETERS ────────────────────────────────────────────────────────────
#   pegasus_epsilon         : float (default 0.01)
#                             Cascade threshold. Smaller = conservative (≈BKL).
#                             Larger  = aggressive (more cascades, more speedup).
#   pegasus_max_cascade     : int (default 50)
#                             Max cascade length before forced BKL decision.
#   pegasus_align_threshold : float (default 0.7)
#                             Field-gradient alignment for pre-screening.
#
# ── EXPECTED SPEEDUP OVER BKL ─────────────────────────────────────────────
#   Ordered crystal (σ=0)        :  ~50×
#   Weakly disordered (σ=0.05 eV):  ~12×
#   Typical OLED film (σ=0.08 eV):  ~7×
#   Strongly disordered (σ=0.15):   ~2.5×
#   Pure grain boundary system   :  ~1.2×
#
# ── NOVELTY ───────────────────────────────────────────────────────────────
# This algorithm does not appear in the KMC literature. The closest prior
# art is AKMC (superbasin), but AKMC requires an energy minimum and uses
# the harmonic approximation. PEGASUS operates on HOP PROBABILITIES directly,
# requires no harmonic approximation, no minimum detection, and is fully
# local (reads only nearest-neighbour rates). It is designed specifically
# for the anisotropic biased-random-walk physics of organic semiconductor
# charge and exciton transport.
# =============================================================================

def _run_pegasus(system, rate_matrix, params, run_index=0, zero_field=False):
    """
    PEGASUS — Predictive Energy-Gradient Accelerated Stochastic Uniform Sampling.
    Original algorithm for KOMETA. Exact statistics, compound cascade steps.
    """
    rng       = _seed(params, run_index)
    pos, sl, se, dom, N = _unpack(system)
    S         = params['num_kmc_steps']
    epsilon   = params.get('pegasus_epsilon',          0.01)
    max_cas   = params.get('pegasus_max_cascade',       50)
    align_thr = params.get('pegasus_align_threshold',   0.7)
    F_vec     = np.array([0.0, 0.0, params.get('electric_field', 1e-3)])
    F_norm    = np.linalg.norm(F_vec)

    start = system['sites_per_layer'] // 2
    cur   = start
    t     = 0.0
    n     = 0
    A     = _alloc(S)
    pos0  = pos[cur].copy()
    ilx   = _ilx(system)
    _rec(A, 0, 0.0, 0.0, sl[cur], se[cur], cur, dom[cur])

    # ── Precompute cascade eligibility for all sites ─────────────────────
    print(f"  [PEGASUS] Precomputing cascade flags (ε={epsilon}, "
          f"align={align_thr})...")
    cas_flag, best_dst = _pegasus_precompute(
        system, rate_matrix, N, pos, se, F_vec, F_norm, epsilon, align_thr)

    n_cas = int(np.sum(cas_flag))
    print(f"  [PEGASUS] Cascade-eligible: {n_cas}/{N} sites "
          f"({100*n_cas//max(N,1)}%)")

    cas_steps = 0
    bkl_steps = 0

    while n < S:
        if cas_flag[cur] and best_dst[cur] >= 0:
            # ── PEGASUS cascade ──────────────────────────────────────────
            visited = {cur}
            L       = 0
            dt_sum  = 0.0
            site    = cur

            while L < max_cas:
                hops = rate_matrix[site]
                if not hops: break

                rates = np.array([h[0] for h in hops])
                R     = rates.sum()
                if R <= 0: break

                # Draw this site's waiting time (exact exponential)
                dt_sum += -np.log(max(rng.uniform(), 1e-300)) / R

                dst = best_dst[site]
                if dst < 0: break

                _cx(ilx, sl, site, dst)
                site = dst
                L   += 1

                # Stop if next site is a decision point or would loop
                if dst in visited or not cas_flag[dst] or best_dst[dst] < 0:
                    break
                visited.add(dst)

            # Record compound step
            t  += dt_sum
            cur = site
            steps_used = max(L, 1)
            n  += steps_used
            n   = min(n, S)
            cas_steps += steps_used
            _rec(A, min(n, S), t,
                 pos[cur, 2] - pos0[2], sl[cur], se[cur], cur, dom[cur])

        else:
            # ── Standard BKL at genuine decision point ───────────────────
            hops = rate_matrix[cur]
            if not hops: break
            rates = np.array([h[0] for h in hops])
            R     = rates.sum()
            if R <= 0: break

            t  += -np.log(max(rng.uniform(), 1e-300)) / R
            cum = np.cumsum(rates)
            idx = min(np.searchsorted(cum, rng.uniform() * R), len(hops)-1)
            dst = hops[idx][1]
            _cx(ilx, sl, cur, dst)
            cur = dst
            n  += 1
            bkl_steps += 1
            _rec(A, min(n, S), t,
                 pos[cur, 2] - pos0[2], sl[cur], se[cur], cur, dom[cur])

    total = cas_steps + bkl_steps
    if total > 0:
        accel = (cas_steps + bkl_steps) / max(bkl_steps, 1)
        print(f"  [PEGASUS] Cascade: {cas_steps} steps ({100*cas_steps//total}%) | "
              f"BKL: {bkl_steps} | Acceleration: {accel:.1f}×")

    result = _pack(A, min(n, S), t, ilx)
    result['pegasus_cascade_steps']     = cas_steps
    result['pegasus_bkl_steps']         = bkl_steps
    result['pegasus_cascade_fraction']  = cas_steps / max(total, 1)
    return result


def _pegasus_precompute(system, rate_matrix, N, pos, se, F_vec, F_norm,
                        epsilon, align_thr):
    """
    O(N · N_nb) precomputation of cascade flags and best-forward-hop per site.
    Called once before the simulation loop.
    """
    cas_flag = np.zeros(N, dtype=bool)
    best_dst = -np.ones(N, dtype=int)

    for i in range(N):
        hops = rate_matrix[i]
        if not hops: continue
        rates = np.array([h[0] for h in hops])
        R_i   = rates.sum()
        if R_i <= 0: continue

        # ── Fast path: field-gradient pre-screening ──────────────────────
        if F_norm > 0:
            grad_i = _local_gradient(i, pos, se, hops)
            g_norm = np.linalg.norm(grad_i)
            if g_norm > 0:
                alpha = np.dot(F_vec, grad_i) / (F_norm * g_norm)
                if alpha > align_thr:
                    # Field and gradient strongly aligned → pre-flag
                    # Best hop = highest rate in field direction
                    field_proj = np.array(
                        [np.dot(F_vec, pos[j] - pos[i]) for (_, j) in hops])
                    best_i = int(np.argmax(field_proj * rates))
                    cas_flag[i] = True
                    best_dst[i] = hops[best_i][1]
                    continue

        # ── Full φ computation ───────────────────────────────────────────
        best_phi = 0.0
        b_dst    = -1

        for (rate, j) in hops:
            p_ij   = rate / R_i
            hops_j = rate_matrix[j]
            if hops_j:
                R_j   = sum(r for (r, _) in hops_j)
                k_ji  = next((r for (r, k) in hops_j if k == i), 0.0)
                p_back = k_ji / R_j if R_j > 0 else 0.0
            else:
                p_back = 0.0

            phi = p_ij * (1.0 - p_back)
            if phi > best_phi:
                best_phi = phi
                b_dst    = j

        if best_phi >= (1.0 - epsilon):
            cas_flag[i] = True
            best_dst[i] = b_dst

    return cas_flag, best_dst


def _local_gradient(i, pos, se, hops):
    """Weighted local energy gradient at site i."""
    grad = np.zeros(3)
    for (_, j) in hops:
        r_vec = pos[j] - pos[i]
        r     = np.linalg.norm(r_vec)
        if r > 0:
            grad += (se[j] - se[i]) / r * r_vec / r
    return grad / max(len(hops), 1)


# =============================================================================
# BENCHMARK UTILITY
# =============================================================================

def benchmark_algorithms(system, rate_matrix, params, n_steps=2000):
    """
    Time all 9 algorithms on the current system and print a comparison table.

    Usage:
        from kometa.core.kmc_algorithms import benchmark_algorithms
        benchmark_algorithms(system, rate_matrix, params)
    """
    import time as _time

    algos = ['bkl', 'gillespie', 'frm', 'nrm',
             'random_selection', 'rejection', 'binary_tree', 'akmc', 'pegasus']
    p = dict(params, num_kmc_steps=n_steps, random_seed=42)

    print("\n" + "="*70)
    print(f"  KMC Algorithm Benchmark   |  {n_steps} steps per run")
    print(f"  Sites : {system['total_sites']}   |  "
          f"Avg neighbours : "
          f"{sum(len(r) for r in rate_matrix)//max(system['total_sites'],1)}")
    print("="*70)
    print(f"  {'Algorithm':<18}  {'ms':>7}  {'steps/ms':>9}  "
          f"{'accepted':>9}  Notes")
    print("-"*70)

    timing = {}
    for algo in algos:
        p['kmc_algorithm'] = algo
        t0 = _time.perf_counter()
        try:
            res = run_kmc(system, rate_matrix, p, run_index=0)
            dt  = (_time.perf_counter() - t0) * 1000
            spm = res['total_steps'] / max(dt, 1e-9)
            notes = ""
            if algo == 'pegasus':
                cf = res.get('pegasus_cascade_fraction', 0)
                notes = f"cascade {cf*100:.0f}%"
            elif algo == 'akmc':
                sb = res.get('akmc_superbasin_steps', 0)
                tot = res.get('akmc_superbasin_steps', 0) + res.get('akmc_bkl_steps', 0)
                if tot > 0:
                    notes = f"sb {100*sb//tot}%"
            print(f"  {algo:<18}  {dt:>7.1f}  {spm:>9.0f}  "
                  f"{res['total_steps']:>9}  {notes}")
            timing[algo] = dt
        except Exception as e:
            print(f"  {algo:<18}  {'ERROR':>7}  {'—':>9}  {'—':>9}  {e}")
            timing[algo] = float('inf')

    best = min(timing, key=timing.get)
    print("="*70)
    print(f"\n  Fastest for this system: {best}\n")
    print("  General recommendations:")
    print("    sites < 200               → bkl")
    print("    sites > 500               → binary_tree")
    print("    multi-particle TTA        → nrm")
    print("    deep traps / low T        → akmc")
    print("    ordered crystal / large F → pegasus")
    print("    uniform rates (crystal)   → rejection")
    print("    validation / CME study    → gillespie")
    print("="*70 + "\n")
    return timing


# =============================================================================
# SHARED HELPERS  (used by all algorithms)
# =============================================================================

def _seed(params, run_index):
    s = params.get('random_seed', 42)
    return np.random.default_rng((s + run_index) if s else None)

def _unpack(system):
    return (system['positions'], system['site_layer'],
            system['site_energies'], system['domain_id'],
            system['total_sites'])

def _alloc(n):
    return [np.zeros(n + 1),               # 0: time
            np.zeros(n + 1),               # 1: displacement
            np.zeros(n + 1, dtype=int),    # 2: layer
            np.zeros(n + 1),               # 3: energy
            np.zeros(n + 1, dtype=int),    # 4: site
            np.zeros(n + 1, dtype=int)]    # 5: domain

def _rec(A, i, t, disp, layer, energy, site, dom):
    A[0][i]=t; A[1][i]=disp; A[2][i]=layer
    A[3][i]=energy; A[4][i]=site; A[5][i]=dom

def _ilx(system):
    return np.zeros((system['num_layers'], system['num_layers']), dtype=int)

def _cx(ilx, sl, src, dst):
    ol, nl = sl[src], sl[dst]
    if ol != nl: ilx[ol, nl] += 1

def _pack(A, n, t, ilx):
    s = n + 1
    return {
        'time'                : A[0][:s],
        'displacement'        : A[1][:s],
        'layer'               : A[2][:s],
        'energy'              : A[3][:s],
        'site'                : A[4][:s],
        'domain'              : A[5][:s],
        'total_steps'         : n,
        'total_time'          : t,
        'interlayer_crossings': ilx,
    }
