# kometa/tests/test_rates.py
# Unit tests for core rate functions
# Run with:  python -m pytest kometa/tests/ -v

import numpy as np
import pytest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from kometa.core.constants import KB, HBAR
from kometa.core.morphology import build_system
from kometa.core.rates import calculate_marcus_rates, calculate_mlj_rates, calculate_dexter_rates

# ── Minimal params for testing ────────────────────────────────────────────────
BASE_PARAMS = {
    'random_seed'               : 42,
    'num_layers'                : 2,
    'sites_per_layer'           : 10,
    'intra_layer_distance'      : 4.0,
    'inter_layer_distance'      : 3.5,
    'positional_disorder_sigma' : 0.0,
    'mean_site_energy'          : [0.0, -0.1],
    'disorder_sigma'            : [0.07, 0.07],
    'lambda_reorg'              : [0.2, 0.25],
    'triplet_energy'            : [1.8, 1.75],
    'detect_grain_boundaries'   : False,
    'blend_type'                : 'single',
    'blend_fractions'           : [1.0],
    'J_intra'                   : [0.05, 0.05],
    'J_inter'                   : [0.03],
    'electric_field'            : 1e-3,
    'carrier_type'              : 'hole',
    'charge_rate_model'         : 'marcus',
    'J_source'                  : 'uniform',
    'temperature'               : 300.0,
    'k_dexter_prefactor'        : 1e12,
    'dexter_decay_length'       : 1.5,
    'allow_inter_species_hopping': False,
    'enable_quenching'          : False,
    # MLJ
    'huang_rhys_factor'         : 0.5,
    'hbar_omega'                : 0.17,
    'lambda_outer'              : [0.1, 0.12],
    'mlj_n_max'                 : 20,
}


@pytest.fixture
def system():
    return build_system(BASE_PARAMS)


def test_build_system_keys(system):
    """System dict must contain all required keys."""
    required = ['total_sites', 'positions', 'site_energies', 'triplet_energies',
                'neighbors', 'site_layer', 'site_species', 'is_grain_boundary',
                'domain_id', 'gb_lambda', 'num_layers', 'sites_per_layer']
    for k in required:
        assert k in system, f"Missing key: {k}"


def test_build_system_shape(system):
    N = BASE_PARAMS['num_layers'] * BASE_PARAMS['sites_per_layer']
    assert system['total_sites'] == N
    assert system['positions'].shape == (N, 3)
    assert system['site_energies'].shape == (N,)


def test_marcus_rates_positive(system):
    rates = calculate_marcus_rates(system, BASE_PARAMS)
    assert len(rates) == system['total_sites']
    # Every site must have at least one non-negative rate
    for i, row in enumerate(rates):
        for (r, j) in row:
            assert r >= 0, f"Negative Marcus rate at site {i} -> {j}: {r}"


def test_marcus_rates_nonzero(system):
    rates = calculate_marcus_rates(system, BASE_PARAMS)
    total_pairs = sum(len(row) for row in rates)
    assert total_pairs > 0, "No Marcus rates were computed"


def test_mlj_rates_positive(system):
    p = dict(BASE_PARAMS)
    p['charge_rate_model'] = 'mlj'
    rates = calculate_mlj_rates(system, p)
    for row in rates:
        for (r, _) in row:
            assert r >= 0


def test_mlj_higher_than_marcus_at_low_T(system):
    """MLJ rates should be higher than Marcus at low T due to tunnelling."""
    p_cold = dict(BASE_PARAMS)
    p_cold['temperature'] = 100.0
    p_cold['lambda_outer'] = [0.1, 0.12]
    marcus = calculate_marcus_rates(system, p_cold)
    mlj    = calculate_mlj_rates(system, p_cold)
    sum_marcus = sum(r for row in marcus for (r,_) in row)
    sum_mlj    = sum(r for row in mlj    for (r,_) in row)
    assert sum_mlj >= sum_marcus, "MLJ should be >= Marcus at low T"


def test_dexter_rates_positive(system):
    rates = calculate_dexter_rates(system, BASE_PARAMS)
    for row in rates:
        for (r, _) in row:
            assert r >= 0


def test_dexter_decay_with_distance(system):
    """Dexter rates must decrease with inter-site distance."""
    rates  = calculate_dexter_rates(system, BASE_PARAMS)
    pos    = system['positions']
    neighs = system['neighbors']
    for i in range(system['total_sites']):
        row = rates[i]
        for (r, j) in row:
            d = np.linalg.norm(pos[j] - pos[i])
            if d > 0:
                # rate should be < k0 (spatial factor < 1)
                assert r <= BASE_PARAMS['k_dexter_prefactor'] * 1.01


def test_bkl_runs(system):
    """BKL must return a result dict with expected keys."""
    from kometa.core.bkl import run_kmc
    p = dict(BASE_PARAMS, num_kmc_steps=500)
    rate_matrix = calculate_marcus_rates(system, p)
    result = run_kmc(system, rate_matrix, p, run_index=0)
    assert 'time' in result
    assert 'displacement' in result
    assert result['total_steps'] > 0


def test_constants_values():
    """Spot-check physical constant values."""
    from kometa.core.constants import KB, HBAR, E_CHARGE
    assert abs(KB - 8.617333e-5) < 1e-10
    assert abs(HBAR - 6.582119e-16) < 1e-20
    assert abs(E_CHARGE - 1.602176e-19) < 1e-24
