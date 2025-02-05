# by Kanegend :)
import sys
from PyQt6 import QtWidgets
from PyQt6 import uic

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from the_first_task import gaussian_beam_propagation_no_vector
from the_second_task import gaussian_beam_propagation_vector


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self):
        super().__init__()
        self.setup_ui()
        self.setup_graphs()
        self.setup_controls()
        self.build_graph.clicked.connect(self.button)
        self.vector.currentTextChanged.connect(self.change)
        self.distance.valueChanged.connect(self.upd)

    def setup_ui(self):
        try:
            uic.loadUi("main_window.ui", self)
        except Exception as e:
            print(f"Error loading UI: {e}")
            sys.exit(1)

    def setup_graphs(self):
        self.canvases = []
        self.toolbars = []
        self.graph_layouts = [
            (self.form_for_tools_1, self.form_for_graph_1),
            (self.form_for_tools_2, self.form_for_graph_2),
            (self.form_for_tools_3, self.form_for_graph_3)
        ]

        for tool_layout, graph_layout in self.graph_layouts:
            canvas = MplCanvas_3D(self, width=10, height=8, dpi=100)
            toolbar = NavigationToolbar(canvas, self)

            tool_layout.addWidget(toolbar)
            graph_layout.addWidget(canvas)

            self.canvases.append(canvas)
            self.toolbars.append(toolbar)

    def setup_controls(self):
        self.plate_thickness.setEnabled(False)
        self.plate_thickness.setValue(0.00000)
        self.epsilons = [self.epsilon_1, self.epsilon_2, self.epsilon_3]
        for epsilon in self.epsilons:
            epsilon.setEnabled(False)
            epsilon.setValue(0.00000)
        self.flag_1 = False
        self.flag_2 = False

    def upd(self):
        if self.distance.value() == 0 and self.vector.currentText() == "Consider the vector nature of the field":
            self.graph_axis.addItems(["Reflected field"])
        elif self.distance.value() != 0 and self.vector.currentText() == "Consider the vector nature of the field":
            index = self.graph_axis.findText("Reflected field")
            if index != -1:
                self.graph_axis.removeItem(index)
        if self.distance.value() >= self.plate_thickness.value() and \
                self.vector.currentText() == "Consider the vector nature of the field" and self.flag_2 is False:
            self.graph_axis.addItems(["Transmitted field"])
            self.flag_2 = True
        elif self.distance.value() < self.plate_thickness.value() and \
                self.vector.currentText() == "Consider the vector nature of the field" and self.flag_2 is True:
            index = self.graph_axis.findText("Transmitted field")
            if index != -1:
                self.graph_axis.removeItem(index)
            self.flag_2 = False

    def change(self):
        if self.flag_1 is False:
            self.configure_parameters(True, d_min=0.00001, d_max=1000, d_standart=0.00675, epsilon_min=0.0001,
                                      epsilon_max=100, epsilon_1=1, epsilon_2=12, epsilon_3=1)
            self.graph_axis.clear()
            new_items = ["Propagation in free space"]
            self.graph_axis.addItems(new_items)
            if self.distance.value() == 0:
                self.graph_axis.addItems(["Reflected field"])
            elif self.distance.value() >= self.plate_thickness.value() and self.flag_2 is False:
                self.graph_axis.addItems(["Transmitted field"])
                self.flag_2 = True
            self.flag_1 = True
        else:
            self.flag_2 = False
            self.configure_parameters(False, d_min=0, d_max=1, d_standart=0, epsilon_min=0,
                                      epsilon_max=100, epsilon_1=0, epsilon_2=0, epsilon_3=0)
            self.graph_axis.clear()
            new_items = ["Propagation in free space"]
            self.graph_axis.addItems(new_items)
            self.flag_1 = False

    def configure_parameters(self, flag, d_min, d_max, d_standart, epsilon_min, epsilon_max, epsilon_1, epsilon_2,
                             epsilon_3):
        self.plate_thickness.setRange(d_min, d_max)
        self.plate_thickness.setEnabled(flag)
        self.plate_thickness.setValue(d_standart)
        self.epsilon_1.setRange(epsilon_min, epsilon_max)
        self.epsilon_1.setEnabled(flag)
        self.epsilon_1.setValue(epsilon_1)
        self.epsilon_2.setRange(epsilon_min, epsilon_max)
        self.epsilon_2.setEnabled(flag)
        self.epsilon_2.setValue(epsilon_2)
        self.epsilon_3.setRange(epsilon_min, epsilon_max)
        self.epsilon_3.setEnabled(flag)
        self.epsilon_3.setValue(epsilon_3)

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
        epsilon_1 = self.epsilon_1.value()
        epsilon_2 = self.epsilon_2.value()
        epsilon_3 = self.epsilon_3.value()
        X, Y, E_x, E_y, E_z, E_x_T_, E_y_T_, E_z_T_, E_x_R_, E_y_R_, E_z_R_ = \
            gaussian_beam_propagation_vector(left, right, step, z, w, d, epsilon_1, epsilon_2, epsilon_3)
        self.create_graphs_2(X, Y, E_x, E_y, E_z, E_x_T_, E_y_T_, E_z_T_, E_x_R_, E_y_R_, E_z_R_)

    def no_vector(self):
        left = self.left_border.value()
        right = self.right_border.value()
        step = self.step.value()
        z = self.distance.value()
        w = self.wave_waist.value()
        X, Y, E_x, E_y, E_z = gaussian_beam_propagation_no_vector(left, right, step, z, w)
        self.create_graphs_1(X, Y, E_x, E_y, E_z)

    def create_graphs_2(self, X, Y, E_x, E_y, E_z, E_x_T_, E_y_T_, E_z_T_, E_x_R_, E_y_R_, E_z_R_):
        colormap = 'inferno'
        z_axis = ["|Ex|", "|Ey|", "|Ez|"]
        graphs_free = [abs(E_x), abs(E_y), abs(E_z)]
        graphs_transmitted = [abs(E_x_T_), abs(E_y_T_), abs(E_z_T_)]
        graphs_reflected = [abs(E_x_R_), abs(E_y_R_), abs(E_z_R_)]
        for i in range(len(self.canvases)):
            canvas = self.canvases[i]
            canvas.axes.cla()
            if self.graph_axis.currentText() == "Propagation in free space":
                canvas.axes.plot_surface(X, Y, graphs_free[i], cmap=colormap)
                canvas.axes.set_zlabel(z_axis[i])
            elif self.graph_axis.currentText() == "Transmitted field":
                canvas.axes.plot_surface(X, Y, graphs_transmitted[i], cmap=colormap)
                canvas.axes.set_zlabel(z_axis[i])
            elif self.graph_axis.currentText() == "Reflected field":
                canvas.axes.plot_surface(X, Y, graphs_reflected[i], cmap=colormap)
                canvas.axes.set_zlabel(z_axis[i])
            canvas.axes.set_xlabel("X")
            canvas.axes.set_ylabel("Y")
            canvas.draw()

    def create_graphs_1(self, X, Y, E_x, E_y, E_z):
        colormap = 'inferno'
        z_axis = ["|Ex|", "|Ey|", "|Ez|"]
        graphs_free = [abs(E_x), abs(E_y), abs(E_z)]
        for i in range(len(self.canvases)):
            canvas = self.canvases[i]
            canvas.axes.cla()
            canvas.axes.plot_surface(X, Y, graphs_free[i], cmap=colormap)
            canvas.axes.set_xlabel("X")
            canvas.axes.set_ylabel("Y")
            canvas.axes.set_zlabel(z_axis[i])
            canvas.draw()


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
