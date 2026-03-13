# =============================================================================
# KOMETA — kometa/io/config.py
# JSON configuration save / load
# =============================================================================
import json
import os
import numpy as np


def _convert(obj):
    """Make numpy types JSON-serialisable."""
    if isinstance(obj, (np.integer,)):  return int(obj)
    if isinstance(obj, (np.floating,)): return float(obj)
    if isinstance(obj, np.ndarray):     return obj.tolist()
    return str(obj)


def save_config(params: dict, path: str) -> None:
    """Save parameter dict to a JSON file."""
    with open(path, 'w') as f:
        json.dump(params, f, indent=2, default=_convert)
    print(f"  [OK] Config saved: {path}")


def load_config(path: str) -> dict:
    """Load parameter dict from a JSON file."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Config not found: {path}")
    with open(path) as f:
        params = json.load(f)
    # Convert lists back to numpy arrays where expected
    for key in ('mean_site_energy', 'disorder_sigma', 'lambda_reorg',
                'J_intra', 'J_inter', 'triplet_energy', 'lambda_outer',
                'blend_fractions', 'species_triplet_energy', 'species_k_decay'):
        if key in params and isinstance(params[key], list):
            params[key] = params[key]   # keep as list — downstream accepts both
    print(f"  [OK] Config loaded: {path}")
    return params
