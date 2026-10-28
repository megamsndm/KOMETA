# kometa/gui/worker.py
import sys, traceback
from PySide6.QtCore import QThread, Signal

class SimWorker(QThread):
    log_signal  = Signal(str)
    done_signal = Signal(bool, str)

    def __init__(self, params: dict):
        super().__init__()
        self.params     = params
        self._stop_flag = False

    def request_stop(self):
        self._stop_flag = True

    def run(self):
        old_o, old_e = sys.stdout, sys.stderr
        worker_self  = self

        class Redirect:
            def __init__(self, sig): self._s = sig; self._b = ""
            def write(self, t):
                self._b += t
                if "\n" in t: self._s.emit(self._b); self._b = ""
            def flush(self): pass

        sys.stdout = sys.stderr = Redirect(self.log_signal)
        try:
            # Works both as package and inside PyInstaller frozen .exe
            try:
                from ..engine import run as engine_run
            except ImportError:
                from kometa.engine import run as engine_run
            engine_run(self.params)
            if self._stop_flag:
                self.done_signal.emit(False, "Stopped by user.")
            else:
                self.done_signal.emit(True, "Simulation completed successfully.")
        except Exception:
            tb = traceback.format_exc()
            if self._stop_flag:
                self.log_signal.emit("\n[INFO] Stopped by user.\n")
                self.done_signal.emit(False, "Stopped by user.")
            else:
                self.log_signal.emit(f"\n[ERROR]\n{tb}")
                self.done_signal.emit(False, tb.splitlines()[-1])
        finally:
            sys.stdout, sys.stderr = old_o, old_e
