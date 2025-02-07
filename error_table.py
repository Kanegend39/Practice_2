from PyQt6 import QtWidgets
from PyQt6 import uic
import sys

class Error_table(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setup_ui()
        self.button_ok.clicked.connect(self.exit)

    def setup_ui(self):
        try:
            uic.loadUi("error_table.ui", self)
        except Exception as e:
            print(f"Error loading UI: {e}")
            sys.exit(1)

    def exit(self):
        self.close()