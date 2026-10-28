# =============================================================================
# KOMETA.spec — PyInstaller build specification
# Tested: PyInstaller 6.19, Python 3.14, Windows 11
# Usage:  pyinstaller KOMETA.spec
# =============================================================================
import os
import sys

block_cipher = None

def optional_file(src, dest):
    """Include a data file only if it exists — avoids build errors."""
    return [(src, dest)] if os.path.isfile(src) else []

# ---------------------------------------------------------------------------
# DATA FILES  — every .py source file included so frozen imports work
# ---------------------------------------------------------------------------
datas = [
    ('kometa/__init__.py',                'kometa'),
    ('kometa/engine.py',                  'kometa'),
    ('kometa/core/__init__.py',           'kometa/core'),
    ('kometa/core/constants.py',          'kometa/core'),
    ('kometa/core/rates.py',              'kometa/core'),
    ('kometa/core/morphology.py',         'kometa/core'),
    ('kometa/core/bkl.py',               'kometa/core'),
    ('kometa/core/kmc_algorithms.py',     'kometa/core'),
    ('kometa/core/quenching.py',          'kometa/core'),
    ('kometa/core/tta.py',               'kometa/core'),
    ('kometa/core/diffusion.py',          'kometa/core'),
    ('kometa/core/advanced.py',           'kometa/core'),
    ('kometa/analysis/__init__.py',       'kometa/analysis'),
    ('kometa/analysis/mobility.py',       'kometa/analysis'),
    ('kometa/analysis/plots.py',          'kometa/analysis'),
    ('kometa/io/__init__.py',             'kometa/io'),
    ('kometa/io/params.py',              'kometa/io'),
    ('kometa/io/config.py',              'kometa/io'),
    ('kometa/io/writer.py',              'kometa/io'),
    ('kometa/gui/__init__.py',            'kometa/gui'),
    ('kometa/gui/style.py',              'kometa/gui'),
    ('kometa/gui/worker.py',             'kometa/gui'),
    ('kometa/gui/header_bar.py',         'kometa/gui'),
    ('kometa/gui/run_panel.py',          'kometa/gui'),
    ('kometa/gui/results_tab.py',        'kometa/gui'),
    ('kometa/gui/main_window.py',        'kometa/gui'),
    ('kometa/gui/tabs/__init__.py',           'kometa/gui/tabs'),
    ('kometa/gui/tabs/general_tab.py',        'kometa/gui/tabs'),
    ('kometa/gui/tabs/morphology_tab.py',     'kometa/gui/tabs'),
    ('kometa/gui/tabs/charge_tab.py',         'kometa/gui/tabs'),
    ('kometa/gui/tabs/mlj_tab.py',            'kometa/gui/tabs'),
    ('kometa/gui/tabs/tta_tab.py',            'kometa/gui/tabs'),
    ('kometa/gui/tabs/blend_tab.py',          'kometa/gui/tabs'),
    ('kometa/gui/tabs/quench_tab.py',         'kometa/gui/tabs'),
]
datas += optional_file('KOMETA_manual.pdf', '.')
datas += optional_file('README.txt',        '.')
datas += optional_file('CITATION.txt',      '.')

a = Analysis(
    ['kometa_gui.py'],
    pathex=['.'],
    binaries=[],
    datas=datas,
    hiddenimports=[
        # Qt
        'PyQt5.sip',
        'PyQt5.QtCore',
        'PyQt5.QtGui',
        'PyQt5.QtWidgets',
        'PyQt5.QtPrintSupport',
        # matplotlib
        'matplotlib',
        'matplotlib.backends.backend_qt5agg',
        'matplotlib.backends.backend_agg',
        'matplotlib.figure',
        'matplotlib.pyplot',
        'matplotlib.style',
        # pyparsing (required by matplotlib >= 3.5 on Python 3.14)
        'pyparsing',
        'pyparsing.helpers',
        'pyparsing.core',
        'pyparsing.common',
        'pyparsing.exceptions',
        'pyparsing.results',
        'pyparsing.unicode',
        'pyparsing.actions',
        # numpy — NumPy 2.x ships with Python 3.14; numpy.core is REMOVED in 2.0
        # DO NOT list numpy.core.* — it will cause ImportError in NumPy 2.x
        'numpy',
        'numpy.lib',
        'numpy.lib.format',
        'numpy.linalg',
        'numpy.fft',
        'numpy.random',
        'numpy.polynomial',
        # scipy
        'scipy',
        'scipy.linalg',
        'scipy.linalg.blas',
        'scipy.linalg.lapack',
        'scipy.optimize',
        'scipy.special',
        'scipy.spatial',
        'scipy.stats',
        # optional
        'h5py',
        # KOMETA — all submodules explicit
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
    ],
    hookspath=['hooks'],
    runtime_hooks=['hooks/rthook_kometa.py'],
    # -----------------------------------------------------------------------
    # IMPORTANT — Python 3.14 + PyInstaller 6.19 KNOWN BUG:
    # Do NOT put 'distutils' in excludes.  PyInstaller's hook-distutils.py
    # calls add_alias_module() and crashes with:
    #   ValueError: Target module "distutils" already imported as
    #               "ExcludedModule('distutils',)."
    # Removing it from excludes lets PyInstaller manage distutils itself.
    # -----------------------------------------------------------------------
    excludes=[
        'tkinter',
        '_tkinter',
        'IPython',
        'jupyter',
        'notebook',
        'pytest',
        '_pytest',
        'pip',
    ],
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='KOMETA',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    onefile=True,
    icon=None,
)
