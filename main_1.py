import sys
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6 import uic
import numpy as np

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from main import gauss


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self):
        super().__init__()
        uic.loadUi("untitled.ui", self)
        self.canvas_1 = MplCanvas_3D(self, width=10, height=8, dpi=100)
        self.toolbar_1 = NavigationToolbar(self.canvas_1, self)
        self.canvas_2 = MplCanvas_3D(self, width=10, height=8, dpi=100)
        self.toolbar_2 = NavigationToolbar(self.canvas_2, self)
        self.canvas_3 = MplCanvas_3D(self, width=10, height=8, dpi=100)
        self.toolbar_3 = NavigationToolbar(self.canvas_3, self)
        self.canvas_4 = MplCanvas_3D(self, width=10, height=8, dpi=100)
        self.toolbar_4 = NavigationToolbar(self.canvas_4, self)
        self.canvas_5 = MplCanvas_3D(self, width=10, height=8, dpi=100)
        self.toolbar_5 = NavigationToolbar(self.canvas_5, self)
        self.canvas_6 = MplCanvas_3D(self, width=10, height=8, dpi=100)
        self.toolbar_6 = NavigationToolbar(self.canvas_6, self)
        self.form1.addWidget(self.toolbar_1)
        self.form2.addWidget(self.canvas_1)
        self.form3.addWidget(self.toolbar_2)
        self.form4.addWidget(self.canvas_2)
        self.form5.addWidget(self.toolbar_3)
        self.form6.addWidget(self.canvas_3)
        self.form7.addWidget(self.toolbar_4)
        self.form8.addWidget(self.canvas_4)
        self.form9.addWidget(self.toolbar_5)
        self.form10.addWidget(self.canvas_5)
        self.form11.addWidget(self.toolbar_6)
        self.form12.addWidget(self.canvas_6)
        self.canvas_1.axes.set_zlabel("|Ex|")
        self.canvas_2.axes.set_zlabel("|Ey|")
        self.canvas_3.axes.set_zlabel("|Ez|")
        self.canvas_4.axes.set_zlabel("|Ex|")
        self.canvas_5.axes.set_zlabel("|Ey|")
        self.canvas_6.axes.set_zlabel("|Ez|")
        self.pushButton1.clicked.connect(self.create_graphs)

    def create_graphs(self):
        self.canvas_1.axes.cla()
        self.canvas_2.axes.cla()
        self.canvas_3.axes.cla()
        self.canvas_4.axes.cla()
        self.canvas_5.axes.cla()
        self.canvas_6.axes.cla()
        self.canvas_1.axes.set_xlabel("X")
        self.canvas_1.axes.set_ylabel("Y")
        self.canvas_2.axes.set_xlabel("X")
        self.canvas_2.axes.set_ylabel("Y")
        self.canvas_3.axes.set_xlabel("X")
        self.canvas_3.axes.set_ylabel("Y")
        self.canvas_4.axes.set_xlabel("X")
        self.canvas_4.axes.set_ylabel("Y")
        self.canvas_5.axes.set_xlabel("X")
        self.canvas_5.axes.set_ylabel("Y")
        self.canvas_6.axes.set_xlabel("X")
        self.canvas_6.axes.set_ylabel("Y")
        self.canvas_1.axes.set_zlabel("|Ex|")
        self.canvas_2.axes.set_zlabel("|Ey|")
        self.canvas_3.axes.set_zlabel("|Ez|")
        self.canvas_4.axes.set_zlabel("|Ex|")
        self.canvas_5.axes.set_zlabel("|Ey|")
        self.canvas_6.axes.set_zlabel("|Ez|")
        right = self.doubleSpinBox.value()
        left = self.doubleSpinBox2.value()
        step = self.doubleSpinBox3.value()
        z = self.doubleSpinBox4.value()
        w = self.doubleSpinBox5.value()
        A = gauss(left, right, step, z, w)
        X, Y, E_x, E_y, E_z, E_x_, E_y_, E_z_ = A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7]
        self.canvas_1.axes.plot_surface(X, Y, abs(E_x_), cmap='inferno')
        self.canvas_2.axes.plot_surface(X, Y, abs(E_y_), cmap='inferno')
        self.canvas_3.axes.plot_surface(X, Y, abs(E_z_), cmap='inferno')
        self.canvas_4.axes.plot_surface(X, Y, abs(E_x), cmap='inferno')
        self.canvas_5.axes.plot_surface(X, Y, abs(E_y), cmap='inferno')
        self.canvas_6.axes.plot_surface(X, Y, abs(E_z), cmap='inferno')
        self.canvas_1.draw()
        self.canvas_2.draw()
        self.canvas_3.draw()
        self.canvas_4.draw()
        self.canvas_5.draw()
        self.canvas_6.draw()


class MplCanvas_3D(FigureCanvas):
    def __init__(self, parent=None, width=10, height=8, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111, projection='3d')
        self.axes.set_facecolor('white')
        self.axes.grid(color='black', linewidth=0.5)
        self.axes.set_xlabel("X")
        self.axes.set_ylabel("Y")
        self.axes.grid(True)
        super(MplCanvas_3D, self).__init__(fig)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())
