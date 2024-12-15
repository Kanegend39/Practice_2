import sys
from PyQt6 import QtWidgets
from PyQt6 import uic

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from the_first_task import gaussian_beam_propagation_no_vector
from the_second_task import gaussian_beam_propagation_vector


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self):
        super().__init__()
        uic.loadUi("main_window.ui", self)
        self.canvas = MplCanvas_3D(self, width=10, height=8, dpi=100)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.form_for_tools.addWidget(self.toolbar)
        self.form_for_graph.addWidget(self.canvas)
        self.build_graph.clicked.connect(self.button)
        self.flag = False
        self.vector.currentTextChanged.connect(self.change)

    def change(self):
        if self.flag is False:
            self.angle.setEnabled(False)
            self.angle.setValue(0.0)
            self.graph_axis.clear()
            new_items = ["Transmitted field on x-axis", "Transmitted field on y-axis",
                         "Transmitted field on z-axis", "Reflected field on x-axis",
                         "Reflected field on y-axis", "Reflected field on z-axis",
                         "Initial field on x-axis", "Initial field on y-axis", "Initial field on z-axis"]
            self.graph_axis.addItems(new_items)
            self.flag = True
        else:
            self.angle.setEnabled(True)
            self.graph_axis.clear()
            new_items = ["Transmitted field (p - polarization)", "Transmitted field (s - polarization)",
                         "Reflected field (p - polarization)", "Reflected field (s - polarization)",
                         "Initial field"]
            self.graph_axis.addItems(new_items)
            self.flag = False

    def button(self):
        if self.vector.currentText() == "Do not consider the vector nature of the field":
            self.no_vector()
        else:
            self.yes_vector()

    def yes_vector(self):
        left = self.left_border.value()
        right = self.right_border.value()
        step = self.step.value()
        z = self.distance.value()
        w = self.wave_waist.value()
        d = self.plate_thickness.value()
        X, Y, E_x, E_y, E_z, E_x_T_, E_y_T_, E_z_T_, E_x_R_, E_y_R_, E_z_R_ = gaussian_beam_propagation_vector(left,
                                                                                                               right,
                                                                                                               step,
                                                                                                               z, w,
                                                                                                               d)
        self.create_graphs_2(X, Y, E_x, E_y, E_z, E_x_T_, E_y_T_, E_z_T_, E_x_R_, E_y_R_, E_z_R_, z)

    def no_vector(self):
        left = self.left_border.value()
        right = self.right_border.value()
        step = self.step.value()
        z = self.distance.value()
        w = self.wave_waist.value()
        d = self.plate_thickness.value()
        tetta = self.angle.value()
        X, Y, U0_x_y, U_x_y_T_p, U_x_y_T_s, U_x_y_R_p, U_x_y_R_s = gaussian_beam_propagation_no_vector(left, right,
                                                                                                       step, z, w,
                                                                                                       tetta, d)
        self.create_graphs_1(X, Y, U0_x_y, U_x_y_T_p, U_x_y_T_s, U_x_y_R_p, U_x_y_R_s, z)

    def create_graphs_2(self, X, Y, E_x, E_y, E_z, E_x_T_, E_y_T_, E_z_T_, E_x_R_, E_y_R_, E_z_R_, z):
        self.canvas.axes.cla()
        colormap = 'inferno'
        if z == 0 or self.graph_axis.currentText() == "Initial field on x-axis":
            self.canvas.axes.plot_surface(X, Y, abs(E_x), cmap=colormap)
            self.canvas.axes.set_zlabel("Ex")
        else:
            if self.graph_axis.currentText() == "Initial field on y-axis":
                self.canvas.axes.plot_surface(X, Y, abs(E_y), cmap=colormap)
                self.canvas.axes.set_zlabel("Ey")
            elif self.graph_axis.currentText() == "Initial field on z-axis":
                self.canvas.axes.plot_surface(X, Y, abs(E_z), cmap=colormap)
                self.canvas.axes.set_zlabel("Ez")
            elif self.graph_axis.currentText() == "Transmitted field on x-axis":
                self.canvas.axes.plot_surface(X, Y, abs(E_x_T_), cmap=colormap)
                self.canvas.axes.set_zlabel("Ex")
            elif self.graph_axis.currentText() == "Transmitted field on y-axis":
                self.canvas.axes.plot_surface(X, Y, abs(E_y_T_), cmap=colormap)
                self.canvas.axes.set_zlabel("Ey")
            elif self.graph_axis.currentText() == "Transmitted field on z-axis":
                self.canvas.axes.plot_surface(X, Y, abs(E_z_T_), cmap=colormap)
                self.canvas.axes.set_zlabel("Ez")
            elif self.graph_axis.currentText() == "Reflected field on x-axis":
                self.canvas.axes.plot_surface(X, Y, abs(E_x_R_), cmap=colormap)
                self.canvas.axes.set_zlabel("Ex")
            elif self.graph_axis.currentText() == "Reflected field on y-axis":
                self.canvas.axes.plot_surface(X, Y, abs(E_y_R_), cmap=colormap)
                self.canvas.axes.set_zlabel("Ey")
            elif self.graph_axis.currentText() == "Reflected field on z-axis":
                self.canvas.axes.plot_surface(X, Y, abs(E_z_R_), cmap=colormap)
                self.canvas.axes.set_zlabel("Ez")
        self.canvas.axes.set_xlabel("X")
        self.canvas.axes.set_ylabel("Y")
        self.canvas.draw()

    def create_graphs_1(self, X, Y, U0_x_y, U_x_y_T_p, U_x_y_T_s, U_x_y_R_p, U_x_y_R_s, z):
        self.canvas.axes.cla()
        colormap = 'inferno'
        if z == 0 or self.graph_axis.currentText() == "Initial field":
            self.canvas.axes.plot_surface(X, Y, abs(U0_x_y), cmap=colormap)
        else:
            if self.graph_axis.currentText() == "Transmitted field (p - polarization)":
                self.canvas.axes.plot_surface(X, Y, abs(U_x_y_T_p), cmap=colormap)
            elif self.graph_axis.currentText() == "Transmitted field (s - polarization)":
                self.canvas.axes.plot_surface(X, Y, abs(U_x_y_T_s), cmap=colormap)
            elif self.graph_axis.currentText() == "Reflected field (p - polarization)":
                self.canvas.axes.plot_surface(X, Y, abs(U_x_y_R_p), cmap=colormap)
            elif self.graph_axis.currentText() == "Reflected field (s - polarization)":
                self.canvas.axes.plot_surface(X, Y, abs(U_x_y_R_s), cmap=colormap)
        self.canvas.axes.set_zlabel("U")
        self.canvas.axes.set_xlabel("X")
        self.canvas.axes.set_ylabel("Y")
        self.canvas.draw()


class MplCanvas_3D(FigureCanvas):
    def __init__(self, parent=None, width=10, height=8, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111, projection='3d')
        self.axes.set_facecolor('white')
        self.axes.grid(color='black', linewidth=0.5)
        self.axes.grid(True)
        super(MplCanvas_3D, self).__init__(fig)


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())
