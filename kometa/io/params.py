# =============================================================================
# KOMETA — kometa/io/params.py
# All default parameters — the ONLY place you edit inputs
# =============================================================================
import numpy as np

def get_input_parameters():
    params = {

        # ------------------------------------------------------------------
        # GENERAL
        # ------------------------------------------------------------------
        'simulation_mode': 'both',
        # Options: 'charge', 'TTA', 'diffusion', 'both', 'full'

        # ------------------------------------------------------------------
        # MORPHOLOGY SOURCE
        # ------------------------------------------------------------------
        # 'grid'     : use simple regular grid (default, no MD needed)
        # 'md_file'  : import from MD simulation output file
        'morphology_source': 'grid',

        # Path to MD output file (used only if morphology_source = 'md_file')
        # Supported formats: .xyz, .gro, .pdb
        'morphology_file': 'morphology.xyz',

        # Which molecule/residue name to select from MD file
        'molecule_name': 'PtOEP',

        # ------------------------------------------------------------------
        # SYSTEM GEOMETRY (used when morphology_source = 'grid')
        # ------------------------------------------------------------------
        'num_layers'          : 2,
        'sites_per_layer'     : 50,
        'intra_layer_distance': 4.0,    # Ang
        'inter_layer_distance': 3.5,    # Ang

        # ------------------------------------------------------------------
        # GRAIN BOUNDARY SETTINGS
        # ------------------------------------------------------------------
        'detect_grain_boundaries'  : True,
        'domain_order_cutoff'      : 0.6,   # S2 order parameter threshold
        # Sites with S2 < cutoff are labelled grain boundary
        'gb_lambda_factor'         : 1.5,   # lambda multiplier at GB
        'gb_J_factor'              : 0.3,   # J multiplier at GB
        'gb_sigma_factor'          : 1.5,   # disorder multiplier at GB

        # ------------------------------------------------------------------
        # POSITIONAL DISORDER
        # ------------------------------------------------------------------
        'positional_disorder_sigma': 0.3,  # Ang -- std of position noise
        # Set to 0 to disable positional disorder

        # ------------------------------------------------------------------
        # MULTI-COMPONENT BLEND
        # ------------------------------------------------------------------
        # 'single'  : one molecular species (default)
        # 'binary'  : two species (e.g. PtOEP + host)
        # 'ternary' : three species
        'blend_type': 'single',

        # Fraction of each component (must sum to 1.0)
        # For binary:  [f_component1, f_component2]
        # For ternary: [f_comp1, f_comp2, f_comp3]
        'blend_fractions': [1.0],

        # Species-specific triplet energies (eV) -- one per species
        'species_triplet_energy': [1.8],

        # Species-specific decay rates (s^-1) -- one per species
        'species_k_decay': [1e5],

        # Whether triplets can hop between different species
        # False = triplets stay on same species (no inter-species transfer)
        'allow_inter_species_hopping': False,

        # ------------------------------------------------------------------
        # MODULE 4B: TRIPLET QUENCHING BY IMPURITIES
        # ------------------------------------------------------------------
        # Impurity/trap sites quench triplet excitons on contact.
        # Physical basis: DOI 10.1103/PhysRevApplied.14.034050
        # Quenching is a SEPARATE process from natural decay (k_decay).
        # When a triplet visits a quencher site it is immediately lost
        # with rate k_quench (Dexter-type irreversible transfer to a
        # low-lying dark state on the impurity).
        #
        # 'enable_quenching'    : turn on impurity quenching module
        # 'quencher_fraction'   : fraction of sites that are quenchers (0-1)
        # 'k_quench'            : quenching rate constant (s^-1)
        #                         Typical: 10^9 - 10^12 s^-1
        #                         (much faster than k_decay ~ 10^5)
        # 'quencher_species_idx': if blend, which species index are quenchers
        #                         None = random sites regardless of species
        # 'quench_on_hop'       : if True, quenching fires when triplet
        #                         HOPS ONTO a quencher site (instantaneous).
        #                         if False, quenching is a rate-based event
        #                         in the BKL event list (more physical).
        'enable_quenching'     : False,
        'quencher_fraction'    : 0.01,       # 1% impurity concentration
        'k_quench'             : 1e10,       # s^-1 (fast quenching)
        'quencher_species_idx' : None,       # None = random; int = species index
        'quench_on_hop'        : False,      # False = BKL event (recommended)

        # ------------------------------------------------------------------
        # SITE ENERGY MODEL
        # ------------------------------------------------------------------
        # 'gaussian' : standard Gaussian random disorder (default)
        # 'dipole'   : electrostatic dipole-dipole interaction model
        # 'qmm_file' : read site energies from external QM/MM output file
        'site_energy_model': 'gaussian',

        # Path to QM/MM site energy file (one energy per line, eV)
        'qmm_energy_file': 'site_energies_qmm.txt',

        # Molecular dipole moment magnitude (Debye) -- used for dipole model
        'molecular_dipole': 2.0,      # Debye

        # ------------------------------------------------------------------
        # TRANSFER INTEGRAL SOURCE
        # ------------------------------------------------------------------
        # 'uniform'   : single J value per layer pair (default)
        # 'orca_table': read pair-specific J values from file
        # 'ml_model'  : predict J from geometry using Gaussian Process
        'J_source': 'uniform',

        # Path to ORCA transfer integral table
        # Format: site_i  site_j  J_eV  (tab-separated)
        'J_table_file': 'transfer_integrals.txt',

        # ------------------------------------------------------------------
        # CHARGE TRANSPORT PARAMETERS
        # ------------------------------------------------------------------
        'mean_site_energy': [0.0, -0.1],      # eV per layer
        'disorder_sigma'  : [0.07, 0.07],     # eV per layer
        'lambda_reorg'    : [0.2, 0.25],       # eV per layer
        'J_intra'         : [0.05, 0.05],      # eV intra-layer coupling
        'J_inter'         : [0.03],            # eV inter-layer coupling

        'carrier_type'    : 'hole',            # 'hole' or 'electron'
        'electric_field'  : 1e-3,              # V/Ang
        'zero_field_runs' : True,              # also compute zero-field mobility

        # ------------------------------------------------------------------
        # RATE MODEL FOR CHARGE TRANSPORT
        # ------------------------------------------------------------------
        # 'marcus'   : classical Marcus theory (default, room temperature)
        # 'mlj'      : Marcus-Levich-Jortner (required for low T, quantum modes)
        'charge_rate_model': 'marcus',

        # MLJ parameters (only used when charge_rate_model = 'mlj')
        'huang_rhys_factor' : 0.5,    # S = lambda_in / hbar_omega
        'hbar_omega'        : 0.17,   # eV -- dominant quantum mode (C=C stretch)
        'lambda_outer'      : [0.15, 0.18],  # eV -- outer-sphere lambda per layer
        'mlj_n_max'         : 20,     # number of vibrational levels to sum

        # ------------------------------------------------------------------
        # TTA PARAMETERS
        # ------------------------------------------------------------------
        'triplet_energy'     : [1.8, 1.75],   # eV per layer (ORCA TDDFT)
        'k_dexter_prefactor' : 1e12,          # s^-1
        'dexter_decay_length': 1.5,           # Ang
        'k_TTA'              : 1e9,           # s^-1
        'k_decay'            : 1e5,           # s^-1 (= 1/lifetime)
        'singlet_fraction'   : 1/9,           # spin statistics
        'phi_emission'       : 0.8,           # fluorescence QY

        'num_triplet_pairs'  : 200,
        'tta_max_time'       : 1e-3,          # s
        'tta_max_steps'      : int(5e4),

        # ------------------------------------------------------------------
        # SINGLE TRIPLET DIFFUSION (Module 6B)
        # ------------------------------------------------------------------
        'diffusion_num_excitons' : 100,       # independent excitons to simulate
        'diffusion_max_steps'    : int(1e5),  # steps per exciton

        # ------------------------------------------------------------------
        # SHARED KMC SETTINGS
        # ------------------------------------------------------------------
        'temperature'     : 300.0,
        'num_kmc_steps'   : int(1e5),
        'num_runs'        : 10,

        # KMC algorithm selection
        # Options:
        #   'bkl'        -- Bortz-Kalos-Lebowitz  (default, O(N), reliable)
        #   'fenwick'    -- Fenwick Tree KMC       (O(log N), best for large systems)
        #   'frm'        -- First Reaction Method  (O(log N), good for multi-particle)
        #   'nrm'        -- Next Reaction Method   (O(log N), best for TTA/bipolar)
        #   'rejection'  -- Rejection KMC          (O(1), only for uniform rates)
        'kmc_algorithm'   : 'bkl',
        # ── Algorithm-specific parameters ─────────────────────────────────
        # AKMC (Accelerated KMC)
        'akmc_fast_rate_threshold' : 1e9,   # s⁻¹ — internal/external separator
        'akmc_max_superbasin_size' : 30,    # max sites per superbasin
        # PEGASUS
        'pegasus_epsilon'          : 0.01,  # cascade threshold (0.001–0.05)
        'pegasus_max_cascade'      : 50,    # max cascade length
        'pegasus_align_threshold'  : 0.7,   # field-gradient alignment cutoff
        'random_seed'     : 42,
        'output_folder'   : 'kmc_results',
        'save_trajectory' : True,

        # ------------------------------------------------------------------
        # MODULE 5D: BIPOLAR TRANSPORT (holes + electrons simultaneously)
        # ------------------------------------------------------------------
        # Simulates both carrier types moving under the same field.
        # Enables: recombination rate, Langevin recombination, OPV modelling.
        # Requires separate Marcus rate matrices for holes and electrons.
        'enable_bipolar'           : False,
        'electron_lambda_reorg'    : [0.2, 0.25],   # eV per layer
        'electron_mean_site_energy': [0.0, -0.1],   # eV per layer (LUMO)
        'electron_disorder_sigma'  : [0.07, 0.07],  # eV per layer
        'electron_J_intra'         : [0.04, 0.04],  # eV
        'electron_J_inter'         : [0.02],         # eV
        'k_recombination'          : 1e6,            # s^-1 Langevin recombination
        'num_bipolar_runs'         : 5,

        # ------------------------------------------------------------------
        # MODULE 5E: MULTI-CARRIER (space charge, N carriers simultaneously)
        # ------------------------------------------------------------------
        # Tracks N holes on the same lattice with Coulomb exclusion.
        # Simple exclusion: no two carriers on same site (Pauli-like).
        # Coulomb interaction: optional energy penalty for nearby carriers.
        'enable_multicarrier'      : False,
        'num_carriers'             : 3,              # number of simultaneous holes
        'coulomb_exclusion_only'   : True,           # False = add Coulomb energy
        'coulomb_epsilon_r'        : 3.5,            # relative permittivity

        # ------------------------------------------------------------------
        # MODULE 5F: TIME-OF-FLIGHT (TOF) SIMULATION
        # ------------------------------------------------------------------
        # Simulates transient photocurrent experiment:
        # Carrier injected at z=0, field drives it to z=L.
        # Records transit time t_tr -> mobility = L^2 / (V * t_tr)
        'enable_tof'               : False,
        'tof_voltage'              : 10.0,           # V (applied voltage)
        'tof_num_carriers'         : 50,             # carriers per TOF run
        'tof_max_time'             : 1e-6,           # s

        # ------------------------------------------------------------------
        # MODULE 5G: TEMPERATURE SCAN μ(T)
        # ------------------------------------------------------------------
        # Automatically runs charge KMC at multiple temperatures.
        # Fits power law: mu ~ T^alpha to classify transport regime.
        'enable_temperature_scan'  : False,
        'temperature_scan_list'    : [150, 200, 250, 300, 350, 400],  # K

        # ------------------------------------------------------------------
        # MODULE 5H: TRAP-LIMITED TRANSPORT (exponential DOS tail)
        # ------------------------------------------------------------------
        # Adds an exponential tail of trap states below the Gaussian DOS.
        # Common in disordered organic films (Urbach tail).
        # g_trap(E) ~ exp(E / E_U)  where E_U = Urbach energy
        'enable_trap_dos'          : False,
        'trap_fraction'            : 0.05,           # fraction of sites that are traps
        'urbach_energy'            : 0.03,           # eV (Urbach tail energy)
        'trap_depth_mean'          : 0.15,           # eV below band edge

        # ------------------------------------------------------------------
        # MODULE 6C: SINGLET EXCITON KMC
        # ------------------------------------------------------------------
        # Förster (FRET) + internal decay for singlet excitons.
        # Extracts singlet diffusion length L_S (short due to fast decay).
        'enable_singlet_kmc'       : False,
        'k_fret_prefactor'         : 1e12,           # s^-1
        'forster_radius'           : 3.0,            # nm (R_0)
        'k_singlet_decay'          : 1e9,            # s^-1 (fast radiative decay)
        'singlet_num_excitons'     : 100,
        'singlet_max_steps'        : int(1e5),

        # ------------------------------------------------------------------
        # MODULE 6D: TRIPLET-POLARON QUENCHING (TPQ)
        # ------------------------------------------------------------------
        # Triplet exciton encounters a trapped charge -> quenched.
        # Major loss channel in OLEDs and TTA films.
        # Rate: k_TPQ when triplet and trapped carrier are neighbours.
        'enable_tpq'               : False,
        'k_tpq'                    : 1e8,            # s^-1
        'polaron_fraction'         : 0.01,           # fraction of sites with trapped carriers

        # ------------------------------------------------------------------
        # MODULE 7B: VORONOI NEIGHBOUR CONSTRUCTION
        # ------------------------------------------------------------------
        # More physical than cutoff for amorphous morphologies.
        # Uses scipy.spatial.Voronoi to find natural neighbours.
        'use_voronoi_neighbors'    : False,          # True = Voronoi, False = cutoff

        # ------------------------------------------------------------------
        # MODULE 7C: SITE ENERGY CORRELATION FUNCTION
        # ------------------------------------------------------------------
        # Computes spatial correlation of site energies C(r) = <de(0)de(r)>.
        # Characterises energetic disorder correlation length.
        'compute_energy_correlation': False,
        'correlation_max_r'        : 20.0,           # Ang
        'correlation_n_bins'       : 50,

        # ------------------------------------------------------------------
        # MODULE 7D: TRAP IDENTIFICATION (graph Laplacian method)
        # ------------------------------------------------------------------
        # Identifies trap sites from the random-walk Laplacian eigenvalues.
        # Based on: DOI 10.1103/PhysRevApplied (Chen et al. 2024)
        'enable_trap_identification': False,
        'trap_id_n_eigenvalues'    : 10,             # eigenvalues to compute

        # ------------------------------------------------------------------
        # MODULE 7E: EXCITON COUPLING (TrESP approximation)
        # ------------------------------------------------------------------
        # Transition charge from electrostatic potentials (TrESP).
        # Fast approximation to BSE exciton coupling without full GW-BSE.
        'enable_tresp_coupling'    : False,
        'tresp_charges_file'       : 'tresp_charges.txt',  # partial charges per atom

        # ------------------------------------------------------------------
        # MODULE 8B: 3D ANISOTROPIC DIFFUSION TENSOR
        # ------------------------------------------------------------------
        # Computes full 3x3 diffusion tensor D_ij from MSD along all axes.
        # Reveals anisotropy: D_xx, D_yy, D_zz and off-diagonal terms.
        'compute_diffusion_tensor' : False,
        'diffusion_tensor_excitons': 200,

        # ------------------------------------------------------------------
        # MODULE 8C: HDF5 STATE FILE OUTPUT
        # ------------------------------------------------------------------
        # Saves all simulation state to HDF5 for post-processing / restart.
        'save_hdf5'                : False,
        'hdf5_filename'            : 'kmc_state.h5',

        # ------------------------------------------------------------------
        # MODULE 8D: SEEBECK COEFFICIENT (THERMOPOWER)
        # ------------------------------------------------------------------
        # Computes S = -<E_carrier> / (q*T) from average site energy
        # of the carrier during transport. Experimental comparison via
        # Seebeck measurement.
        'compute_seebeck'          : False,

        # ------------------------------------------------------------------
        # MODULE 8E: EXCITATION DENSITY DEPENDENCE Phi_UC(I)
        # ------------------------------------------------------------------
        # Varies num_triplet_pairs to simulate different excitation densities.
        # Produces threshold curve Phi_UC vs I (intensity).
        'enable_intensity_scan'    : False,
        'intensity_scan_pairs'     : [20, 50, 100, 200, 500, 1000],

        # ------------------------------------------------------------------
        # MODULE 8F: COARSE-GRAINING (Iterative Boltzmann Inversion)
        # ------------------------------------------------------------------
        # Derives effective coarse-grained rates from fine-grained KMC.
        # Groups sites into supercells and computes effective inter-cell rates.
        'enable_coarse_graining'   : False,
        'cg_supercell_size'        : 3,              # sites per supercell edge
    }
    return params


# =============================================================================
# MODULE 2: SYSTEM SETUP (original grid builder)
# =============================================================================

