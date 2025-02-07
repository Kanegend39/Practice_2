# by Kanegend :)
import sys
from PyQt6 import QtWidgets
from PyQt6 import uic


from layer import Layer
from n_layers import N_Layers


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setup_ui()
        self.version_1.clicked.connect(self.first_version)
        self.version_2.clicked.connect(self.second_version)

    def first_version(self):
        self.layer_window = Layer(self)
        self.layer_window.show()
        self.close()

    def second_version(self):
        self.n_layers_window = N_Layers(self)
        self.n_layers_window.show()
        self.close()

    def setup_ui(self):
        try:
            uic.loadUi("main_window.ui", self)
        except Exception as e:
            print(f"Error loading UI: {e}")
            sys.exit(1)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())
