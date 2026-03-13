# =============================================================================
# PyInstaller hook for the kometa package
# Lists all submodules EXPLICITLY instead of using collect_submodules().
# collect_submodules() triggers setuptools -> distutils on Python 3.14
# which conflicts with PyInstaller 6.19's internal hook-distutils.py.
# =============================================================================

hiddenimports = [
    'kometa',
    'kometa.engine',
    'kometa.core',
    'kometa.core.constants',
    'kometa.core.morphology',
    'kometa.core.rates',
    'kometa.core.bkl',
    'kometa.core.kmc_algorithms',
    'kometa.core.quenching',
    'kometa.core.tta',
    'kometa.core.diffusion',
    'kometa.core.advanced',
    'kometa.analysis',
    'kometa.analysis.mobility',
    'kometa.analysis.plots',
    'kometa.io',
    'kometa.io.params',
    'kometa.io.config',
    'kometa.io.writer',
    'kometa.gui',
    'kometa.gui.style',
    'kometa.gui.worker',
    'kometa.gui.header_bar',
    'kometa.gui.run_panel',
    'kometa.gui.results_tab',
    'kometa.gui.main_window',
    'kometa.gui.tabs',
    'kometa.gui.tabs.general_tab',
    'kometa.gui.tabs.morphology_tab',
    'kometa.gui.tabs.charge_tab',
    'kometa.gui.tabs.mlj_tab',
    'kometa.gui.tabs.tta_tab',
    'kometa.gui.tabs.blend_tab',
    'kometa.gui.tabs.quench_tab',
]

# No collect_data_files() call — source .py files are listed in datas in the spec
datas = []
