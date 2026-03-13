# =============================================================================
# KOMETA Runtime Hook — rthook_kometa.py
# Ensures the frozen package root is on sys.path so that
# "import kometa" and all subpackage imports work inside the .exe
# =============================================================================
import sys
import os

# When running as a PyInstaller .exe, sys._MEIPASS contains the
# extracted temp directory. We add it to sys.path so that
# "from kometa.core.rates import ..." etc. resolve correctly.
if hasattr(sys, '_MEIPASS'):
    meipass = sys._MEIPASS
    if meipass not in sys.path:
        sys.path.insert(0, meipass)

    # Also ensure matplotlib finds its data files
    os.environ.setdefault('MPLCONFIGDIR', os.path.join(meipass, 'matplotlib'))
