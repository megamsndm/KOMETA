# =============================================================================
# KOMETA — kometa/core/quenching.py
# =============================================================================
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

from .constants import KB, HBAR, ANG_TO_CM, E_CHARGE

def assign_quencher_sites(system, params):
    """
    Randomly designate a fraction of molecular sites as quencher/impurity sites.

    Parameters
    ----------
    system : dict   -- system dict from Module 2 / 2B
    params : dict   -- input parameters from Module 1

    Returns
    -------
    is_quencher : np.ndarray, bool, shape (total_sites,)
    """
    if not params.get('enable_quenching', False):
        total_sites = system['total_sites']
        print("  MODULE 4B: Quenching disabled -- all sites are clean host.")
        return np.zeros(total_sites, dtype=bool)

    print("\n" + "="*60)
    print("  MODULE 4B: Assigning impurity / quencher sites...")
    print("="*60)

    total_sites     = system['total_sites']
    f_q             = params.get('quencher_fraction', 0.01)
    q_species_idx   = params.get('quencher_species_idx', None)
    site_species    = system['site_species']
    rng             = np.random.default_rng(
        params['random_seed'] + 77777 if params['random_seed'] else None
    )

    is_quencher = np.zeros(total_sites, dtype=bool)

    if q_species_idx is not None:
        # Quenchers are a specific species (e.g. impurity molecules in blend)
        candidate_sites = np.where(site_species == q_species_idx)[0]
        n_q = len(candidate_sites)
        print(f"  Quenchers = species {q_species_idx} ({n_q} sites, fraction {n_q/total_sites:.3f})")
        is_quencher[candidate_sites] = True
    else:
        # Random impurity distribution -- independent of species
        n_q = max(1, int(round(f_q * total_sites)))
        quencher_idx = rng.choice(total_sites, size=n_q, replace=False)
        is_quencher[quencher_idx] = True
        print(f"  Random quencher sites : {n_q} / {total_sites}  ({100*n_q/total_sites:.2f}%)")

    k_q = params.get('k_quench', 1e10)
    mode = "instantaneous (quench_on_hop)" if params.get('quench_on_hop', False) \
           else "BKL rate-based event"
    print(f"  k_quench              : {k_q:.2e} s^-1")
    print(f"  Quenching mode        : {mode}")
    print(f"  [OK] {np.sum(is_quencher)} quencher sites assigned.")
    return is_quencher


# =============================================================================
# MODULE 5: KMC CORE -- CHARGE TRANSPORT (BKL algorithm, single carrier)
# =============================================================================

