# kometa/core/bkl.py
# Backward-compatible shim — imports from kmc_algorithms.py
# run_kmc() now dispatches based on params['kmc_algorithm']
# Default is 'bkl' so existing code is unaffected.

from .kmc_algorithms import run_kmc, benchmark_algorithms  # noqa: F401
