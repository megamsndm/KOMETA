# =============================================================================
# KOMETA — kometa/core/morphology.py
# =============================================================================
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

from .constants import KB, HBAR, ANG_TO_CM, E_CHARGE

def build_system(params):
    print("\n" + "="*60)
    print("  MODULE 2: Building the molecular site lattice...")
    print("="*60)

    if params['random_seed'] is not None:
        np.random.seed(params['random_seed'])

    num_layers      = params['num_layers']
    sites_per_layer = params['sites_per_layer']
    total_sites     = num_layers * sites_per_layer

    positions = np.zeros((total_sites, 3))
    for layer_idx in range(num_layers):
        for site_idx in range(sites_per_layer):
            global_idx = layer_idx * sites_per_layer + site_idx
            positions[global_idx, 0] = site_idx * params['intra_layer_distance']
            positions[global_idx, 1] = 0.0
            positions[global_idx, 2] = layer_idx * params['inter_layer_distance']

    # Apply positional disorder
    pos_sigma = params.get('positional_disorder_sigma', 0.0)
    if pos_sigma > 0:
        positions += np.random.normal(0, pos_sigma, positions.shape)
        print(f"  [OK] Positional disorder applied: sigma = {pos_sigma:.2f} Ang")

    site_energies = np.zeros(total_sites)
    for layer_idx in range(num_layers):
        mean_E = params['mean_site_energy'][layer_idx]
        sigma  = params['disorder_sigma'][layer_idx]
        for site_idx in range(sites_per_layer):
            global_idx = layer_idx * sites_per_layer + site_idx
            site_energies[global_idx] = np.random.normal(mean_E, sigma)

    triplet_energies = np.zeros(total_sites)
    for layer_idx in range(num_layers):
        mean_T1 = params['triplet_energy'][layer_idx]
        sigma   = params['disorder_sigma'][layer_idx]
        for site_idx in range(sites_per_layer):
            global_idx = layer_idx * sites_per_layer + site_idx
            triplet_energies[global_idx] = np.random.normal(mean_T1, sigma)

    neighbors = [[] for _ in range(total_sites)]
    for layer_idx in range(num_layers):
        for site_idx in range(sites_per_layer):
            global_idx = layer_idx * sites_per_layer + site_idx
            if site_idx > 0:
                neighbors[global_idx].append(
                    ('intra', layer_idx * sites_per_layer + site_idx - 1, layer_idx))
            if site_idx < sites_per_layer - 1:
                neighbors[global_idx].append(
                    ('intra', layer_idx * sites_per_layer + site_idx + 1, layer_idx))
            if layer_idx > 0:
                neighbors[global_idx].append(
                    ('inter', (layer_idx-1) * sites_per_layer + site_idx, layer_idx-1))
            if layer_idx < num_layers - 1:
                neighbors[global_idx].append(
                    ('inter', (layer_idx+1) * sites_per_layer + site_idx, layer_idx+1))

    site_layer = np.zeros(total_sites, dtype=int)
    for layer_idx in range(num_layers):
        for site_idx in range(sites_per_layer):
            site_layer[layer_idx * sites_per_layer + site_idx] = layer_idx

    # Species assignment (for blend morphologies)
    site_species = _assign_species(total_sites, params)

    # Grain boundary detection
    is_grain_boundary = np.zeros(total_sites, dtype=bool)
    domain_id         = np.zeros(total_sites, dtype=int)
    if params.get('detect_grain_boundaries', False):
        is_grain_boundary, domain_id = _detect_grain_boundaries(
            positions, neighbors, params)
        n_gb = np.sum(is_grain_boundary)
        print(f"  [OK] Grain boundary sites  : {n_gb} ({n_gb/total_sites*100:.1f}%)")

    # Modify energies/lambda at grain boundaries
    gb_lambda = np.array(
        [params['lambda_reorg'][site_layer[i]] for i in range(total_sites)])
    if params.get('detect_grain_boundaries', False):
        for i in range(total_sites):
            if is_grain_boundary[i]:
                site_energies[i] += np.random.normal(
                    0, params['disorder_sigma'][site_layer[i]] *
                    (params['gb_sigma_factor'] - 1.0))
                gb_lambda[i] *= params['gb_lambda_factor']

    system = {
        'total_sites'       : total_sites,
        'positions'         : positions,
        'site_energies'     : site_energies,
        'triplet_energies'  : triplet_energies,
        'neighbors'         : neighbors,
        'site_layer'        : site_layer,
        'site_species'      : site_species,
        'is_grain_boundary' : is_grain_boundary,
        'domain_id'         : domain_id,
        'gb_lambda'         : gb_lambda,
        'num_layers'        : num_layers,
        'sites_per_layer'   : sites_per_layer,
    }

    print(f"  [OK] Total sites           : {total_sites}")
    print(f"  [OK] Layers                : {num_layers}")
    print(f"  [OK] Sites per layer       : {sites_per_layer}")
    print(f"  [OK] Charge energy range   : [{site_energies.min():.3f}, {site_energies.max():.3f}] eV")
    print(f"  [OK] Triplet energy range  : [{triplet_energies.min():.3f}, {triplet_energies.max():.3f}] eV")
    return system


def _assign_species(total_sites, params):
    """Assign molecular species to each site based on blend fractions."""
    blend_type = params.get('blend_type', 'single')
    fractions  = params.get('blend_fractions', [1.0])
    species    = np.zeros(total_sites, dtype=int)
    if blend_type == 'single':
        return species
    cum = np.cumsum(fractions)
    r   = np.random.uniform(0, 1, total_sites)
    for i in range(total_sites):
        for s, c in enumerate(cum):
            if r[i] <= c:
                species[i] = s
                break
    n_sp = [np.sum(species == s) for s in range(len(fractions))]
    print(f"  [OK] Species counts        : {n_sp}")
    return species


def _detect_grain_boundaries(positions, neighbors, params):
    """
    Detect grain boundaries using a local orientational order parameter.
    Sites with low S2 (orientational order) are flagged as grain boundaries.
    Uses the bond-orientational order approach: for each site, compute the
    mean cosine-squared of bond angles with neighbors relative to z-axis.
    """
    total_sites = len(positions)
    S2 = np.zeros(total_sites)
    cutoff = params.get('domain_order_cutoff', 0.6)

    for i in range(total_sites):
        angles = []
        for (_, j, _) in neighbors[i]:
            r_vec = positions[j] - positions[i]
            r_norm = np.linalg.norm(r_vec)
            if r_norm > 0:
                cos_theta = r_vec[2] / r_norm
                angles.append(cos_theta)
        if angles:
            cos_arr = np.array(angles)
            S2[i] = np.abs(np.mean(1.5 * cos_arr**2 - 0.5))
        else:
            S2[i] = 1.0

    is_gb = S2 < cutoff

    # Simple domain labelling: connected components of non-GB sites
    domain_id = -np.ones(total_sites, dtype=int)
    current_domain = 0
    for start in range(total_sites):
        if domain_id[start] == -1 and not is_gb[start]:
            stack = [start]
            while stack:
                node = stack.pop()
                if domain_id[node] != -1:
                    continue
                domain_id[node] = current_domain
                for (_, j, _) in neighbors[node]:
                    if not is_gb[j] and domain_id[j] == -1:
                        stack.append(j)
            current_domain += 1
    # GB sites get domain -1
    n_domains = current_domain
    print(f"  [OK] Crystalline domains   : {n_domains}")
    return is_gb, domain_id


# =============================================================================
# MODULE 2B: MD MORPHOLOGY IMPORTER
# =============================================================================



def import_md_morphology(params):
    """
    Imports molecular positions from an MD simulation output file.
    Supported formats: .xyz, .gro, .pdb
    Returns a system dict compatible with the rest of the code.

    USAGE:
      Set params['morphology_source'] = 'md_file'
      Set params['morphology_file'] to the path of your MD output.
      Set params['molecule_name'] to select the molecule of interest.

    The function:
      1. Reads all atomic positions
      2. Computes molecular centres of mass
      3. Builds a neighbour list based on distance cutoff
      4. Assigns layer indices from z-coordinate clustering
      5. Returns a system dict for use in rate calculation modules
    """
    print("\n" + "="*60)
    print("  MODULE 2B: Importing MD morphology...")
    print("="*60)

    filepath = params['morphology_file']
    if not os.path.exists(filepath):
        print(f"  [WARNING] Morphology file not found: {filepath}")
        print("  [WARNING] Falling back to grid morphology.")
        return build_system(params)

    ext = os.path.splitext(filepath)[1].lower()

    if ext == '.xyz':
        positions, mol_ids = _read_xyz_file(filepath, params['molecule_name'])
    elif ext == '.gro':
        positions, mol_ids = _read_gro_file(filepath, params['molecule_name'])
    elif ext == '.pdb':
        positions, mol_ids = _read_pdb_file(filepath, params['molecule_name'])
    else:
        print(f"  [WARNING] Unsupported format: {ext}. Falling back to grid.")
        return build_system(params)

    total_sites = len(positions)
    print(f"  [OK] Molecules loaded      : {total_sites}")

    # Assign layers by z-coordinate clustering
    site_layer = _assign_layers_by_z(positions, params['num_layers'])

    # Assign site energies from chosen model
    site_energies, triplet_energies = _compute_site_energies(
        positions, site_layer, params)

    # Build neighbour list with distance cutoff
    cutoff = params['inter_layer_distance'] * 2.5
    neighbors = _build_neighbor_list_cutoff(positions, cutoff)

    # Positional disorder already encoded in MD positions
    site_species = _assign_species(total_sites, params)
    is_grain_boundary, domain_id = (
        _detect_grain_boundaries(positions, neighbors, params)
        if params.get('detect_grain_boundaries', False)
        else (np.zeros(total_sites, dtype=bool), np.zeros(total_sites, dtype=int))
    )

    gb_lambda = np.array(
        [params['lambda_reorg'][min(site_layer[i], len(params['lambda_reorg'])-1)]
         for i in range(total_sites)])

    system = {
        'total_sites'      : total_sites,
        'positions'        : positions,
        'site_energies'    : site_energies,
        'triplet_energies' : triplet_energies,
        'neighbors'        : neighbors,
        'site_layer'       : site_layer,
        'site_species'     : site_species,
        'is_grain_boundary': is_grain_boundary,
        'domain_id'        : domain_id,
        'gb_lambda'        : gb_lambda,
        'num_layers'       : params['num_layers'],
        'sites_per_layer'  : total_sites // params['num_layers'],
    }
    return system


def _read_xyz_file(filepath, molecule_name):
    """Parse .xyz file, compute centre of mass per molecule."""
    positions = []
    mol_ids   = []
    with open(filepath, 'r') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        try:
            n_atoms = int(lines[i].strip())
        except ValueError:
            i += 1
            continue
        comment = lines[i+1].strip()
        mol_coords = []
        for j in range(n_atoms):
            parts = lines[i+2+j].split()
            if len(parts) >= 4:
                mol_coords.append([float(parts[1]),
                                   float(parts[2]),
                                   float(parts[3])])
        if mol_coords:
            com = np.mean(mol_coords, axis=0)
            positions.append(com)
            mol_ids.append(len(positions)-1)
        i += n_atoms + 2
    return np.array(positions), mol_ids


def _read_gro_file(filepath, molecule_name):
    """Parse GROMACS .gro file, compute centre of mass per residue."""
    from collections import defaultdict
    residue_coords = defaultdict(list)
    with open(filepath, 'r') as f:
        lines = f.readlines()
    for line in lines[2:-1]:
        if len(line) < 44:
            continue
        try:
            res_id   = int(line[0:5].strip())
            res_name = line[5:10].strip()
            if molecule_name.lower() not in res_name.lower():
                continue
            x = float(line[20:28]) * 10.0  # nm to Ang
            y = float(line[28:36]) * 10.0
            z = float(line[36:44]) * 10.0
            residue_coords[res_id].append([x, y, z])
        except (ValueError, IndexError):
            continue
    positions = [np.mean(coords, axis=0)
                 for coords in residue_coords.values()]
    return np.array(positions), list(residue_coords.keys())


def _read_pdb_file(filepath, molecule_name):
    """Parse .pdb file, compute centre of mass per residue."""
    from collections import defaultdict
    residue_coords = defaultdict(list)
    with open(filepath, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        res_name = line[17:20].strip()
        if molecule_name.lower() not in res_name.lower():
            continue
        try:
            res_id = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            residue_coords[res_id].append([x, y, z])
        except (ValueError, IndexError):
            continue
    positions = [np.mean(coords, axis=0)
                 for coords in residue_coords.values()]
    return np.array(positions), list(residue_coords.keys())


def _assign_layers_by_z(positions, num_layers):
    """Assign layer index to each site by clustering z-coordinates."""
    from scipy.cluster.vq import kmeans2
    z_vals = positions[:, 2].reshape(-1, 1)
    centroids, labels = kmeans2(z_vals, num_layers, seed=42)
    # Sort labels so layer 0 = smallest z
    order  = np.argsort(centroids[:, 0])
    remap  = {old: new for new, old in enumerate(order)}
    labels = np.array([remap[l] for l in labels])
    return labels


def _build_neighbor_list_cutoff(positions, cutoff):
    """Build neighbor list using spatial KD-tree and distance cutoff."""
    tree = cKDTree(positions)
    pairs = tree.query_pairs(cutoff)
    total_sites = len(positions)
    neighbors   = [[] for _ in range(total_sites)]
    for (i, j) in pairs:
        layer_i = 0  # will be updated after layer assignment
        layer_j = 0
        hop_type = 'intra'
        neighbors[i].append((hop_type, j, layer_j))
        neighbors[j].append((hop_type, i, layer_i))
    return neighbors


def _update_neighbor_hop_types(neighbors, site_layer):
    """Update hop types (intra/inter) in neighbor list after layer assignment."""
    total_sites = len(neighbors)
    updated = [[] for _ in range(total_sites)]
    for i in range(total_sites):
        for (_, j, _) in neighbors[i]:
            hop_type = 'intra' if site_layer[i] == site_layer[j] else 'inter'
            updated[i].append((hop_type, j, site_layer[j]))
    return updated


def _compute_site_energies(positions, site_layer, params):
    """Compute site energies based on the chosen model."""
    total_sites = len(positions)
    model = params.get('site_energy_model', 'gaussian')

    if model == 'qmm_file':
        qmm_file = params.get('qmm_energy_file', '')
        if os.path.exists(qmm_file):
            energies = np.loadtxt(qmm_file)
            if len(energies) == total_sites:
                print(f"  [OK] QM/MM energies loaded from {qmm_file}")
                site_energies = energies
            else:
                print("  [WARNING] QM/MM file size mismatch, using Gaussian")
                model = 'gaussian'
        else:
            print(f"  [WARNING] QM/MM file not found: {qmm_file}")
            model = 'gaussian'

    if model == 'dipole':
        site_energies = _compute_dipole_energies(positions, site_layer, params)
    elif model != 'qmm_file':
        # Default: Gaussian
        site_energies = np.zeros(total_sites)
        for i in range(total_sites):
            l = min(site_layer[i], len(params['mean_site_energy'])-1)
            site_energies[i] = np.random.normal(
                params['mean_site_energy'][l],
                params['disorder_sigma'][l])

    triplet_energies = np.zeros(total_sites)
    for i in range(total_sites):
        l = min(site_layer[i], len(params['triplet_energy'])-1)
        s = min(site_layer[i], len(params['disorder_sigma'])-1)
        triplet_energies[i] = np.random.normal(
            params['triplet_energy'][l], params['disorder_sigma'][s])

    return site_energies, triplet_energies


def _compute_dipole_energies(positions, site_layer, params):
    """
    Compute site energies from electrostatic dipole-dipole interactions.
    E_i = sum_{j != i} p_i . r_ij / |r_ij|^3
    Uses random dipole orientations drawn from uniform sphere distribution.
    This produces spatially correlated disorder, unlike independent Gaussian.
    """
    total_sites = len(positions)
    p_mag = params.get('molecular_dipole', 2.0) * 0.2082   # Debye -> eA
    # Random dipole directions on unit sphere
    theta = np.arccos(1 - 2*np.random.uniform(0, 1, total_sites))
    phi   = np.random.uniform(0, 2*np.pi, total_sites)
    dipoles = p_mag * np.column_stack([
        np.sin(theta)*np.cos(phi),
        np.sin(theta)*np.sin(phi),
        np.cos(theta)])

    site_energies = np.zeros(total_sites)
    cutoff_r = 15.0  # Ang -- truncation radius
    for i in range(total_sites):
        E_i = 0.0
        for j in range(total_sites):
            if i == j:
                continue
            r_vec = positions[j] - positions[i]
            r     = np.linalg.norm(r_vec)
            if r < cutoff_r and r > 0.5:
                E_i += np.dot(dipoles[j], r_vec) / r**3
        site_energies[i] = E_i

    # Add mean layer energy offset
    for i in range(total_sites):
        l = min(site_layer[i], len(params['mean_site_energy'])-1)
        site_energies[i] += params['mean_site_energy'][l]

    print(f"  [OK] Dipole disorder sigma : {site_energies.std():.4f} eV")
    return site_energies


# =============================================================================
# MODULE 3: MARCUS RATE CALCULATOR (charge transport)
# =============================================================================

