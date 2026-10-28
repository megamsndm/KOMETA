#!/usr/bin/env python3
"""
KOMETA — KMC Organic Molecular Exciton Transport Analysis
Entry point. Run:  python kometa_gui.py
"""
import sys
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QColor, QPalette
from kometa.gui.main_window import MainWindow
from kometa.gui.style import C_BG, C_WHITE, C_TEAL, C_PANEL


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("KOMETA")
    app.setStyle("Fusion")

    pal = QPalette()
    pal.setColor(QPalette.Window,          QColor(C_BG))
    pal.setColor(QPalette.WindowText,      QColor(C_WHITE))
    pal.setColor(QPalette.Base,            QColor(C_PANEL))
    pal.setColor(QPalette.AlternateBase,   QColor(C_BG))
    pal.setColor(QPalette.Text,            QColor(C_WHITE))
    pal.setColor(QPalette.Button,          QColor(C_TEAL))
    pal.setColor(QPalette.ButtonText,      QColor(C_BG))
    pal.setColor(QPalette.Highlight,       QColor(C_TEAL))
    pal.setColor(QPalette.HighlightedText, QColor(C_BG))
    app.setPalette(pal)

    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
