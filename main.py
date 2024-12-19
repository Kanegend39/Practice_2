# by Kanegend :)
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
        try:
            uic.loadUi("main_window.ui", self)
        except Exception as e:
            print(f"Error loading UI: {e}")
            sys.exit(1)
        self.canvas = MplCanvas_3D(self, width=10, height=8, dpi=100)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.form_for_tools.addWidget(self.toolbar)
        self.form_for_graph.addWidget(self.canvas)
        self.build_graph.clicked.connect(self.button)
        self.plate_thickness.setEnabled(False)
        self.plate_thickness.setValue(0.00000)
        self.flag_1 = False
        self.flag_2 = False
        self.vector.currentTextChanged.connect(self.change)
        self.distance.valueChanged.connect(self.upd)

    def upd(self):
        if self.distance.value() == 0 and self.vector.currentText() == "Consider the vector nature of the field":
            self.graph_axis.addItems(["Reflected field on x-axis", "Reflected field on y-axis",
                                      "Reflected field on z-axis"])
        elif self.distance.value() != 0 and self.vector.currentText() == "Consider the vector nature of the field":
            index = self.graph_axis.findText("Reflected field on x-axis")
            if index != -1:
                self.graph_axis.removeItem(index)
            index = self.graph_axis.findText("Reflected field on y-axis")
            if index != -1:
                self.graph_axis.removeItem(index)
            index = self.graph_axis.findText("Reflected field on z-axis")
            if index != -1:
                self.graph_axis.removeItem(index)
        if self.distance.value() >= self.plate_thickness.value() and self.vector.currentText() == "Consider the vector nature of the field" and self.flag_2 is False:
            self.graph_axis.addItems(["Transmitted field on x-axis", "Transmitted field on y-axis",
                                      "Transmitted field on z-axis"])
            self.flag_2 = True
        elif self.distance.value() < self.plate_thickness.value() and self.vector.currentText() == "Consider the vector nature of the field" and self.flag_2 is True:
            index = self.graph_axis.findText("Transmitted field on x-axis")
            if index != -1:
                self.graph_axis.removeItem(index)
            index = self.graph_axis.findText("Transmitted field on y-axis")
            if index != -1:
                self.graph_axis.removeItem(index)
            index = self.graph_axis.findText("Transmitted field on z-axis")
            if index != -1:
                self.graph_axis.removeItem(index)
            self.flag_2 = False

    def change(self):
        if self.flag_1 is False:
            self.plate_thickness.setRange(0.00001, 1)
            self.plate_thickness.setEnabled(True)
            self.plate_thickness.setValue(0.00675)
            self.graph_axis.clear()
            new_items = ["Free propagation x-axis", "Free propagation y-axis",
                         "Free propagation z-axis"]
            self.graph_axis.addItems(new_items)
            if self.distance.value() == 0:
                self.graph_axis.addItems(["Reflected field on x-axis", "Reflected field on y-axis",
                                          "Reflected field on z-axis"])
            elif self.distance.value() >= self.plate_thickness.value() and self.flag_2 is False:
                self.graph_axis.addItems(["Transmitted field on x-axis", "Transmitted field on y-axis",
                                          "Transmitted field on z-axis"])
                self.flag_2 = True
            self.flag_1 = True
        else:
            self.flag_2 = False
            self.plate_thickness.setRange(0.00000, 1)
            self.plate_thickness.setEnabled(False)
            self.plate_thickness.setValue(0.00000)
            self.graph_axis.clear()
            new_items = ["x-axis field",
                         "y-axis field",
                         "z-axis field"]
            self.graph_axis.addItems(new_items)
            self.flag_1 = False

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
        self.canvas.axes.cla()
        colormap = 'inferno'
        if self.graph_axis.currentText() == "Free propagation x-axis":
            self.canvas.axes.plot_surface(X, Y, abs(E_x), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ex|")
        elif self.graph_axis.currentText() == "Free propagation y-axis":
            self.canvas.axes.plot_surface(X, Y, abs(E_y), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ey|")
        elif self.graph_axis.currentText() == "Free propagation z-axis":
            self.canvas.axes.plot_surface(X, Y, abs(E_z), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ez|")
        elif self.graph_axis.currentText() == "Transmitted field on x-axis":
            self.canvas.axes.plot_surface(X, Y, abs(E_x_T_), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ex|")
        elif self.graph_axis.currentText() == "Transmitted field on y-axis":
            self.canvas.axes.plot_surface(X, Y, abs(E_y_T_), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ey|")
        elif self.graph_axis.currentText() == "Transmitted field on z-axis":
            self.canvas.axes.plot_surface(X, Y, abs(E_z_T_), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ez|")
        elif self.graph_axis.currentText() == "Reflected field on x-axis":
            self.canvas.axes.plot_surface(X, Y, abs(E_x_R_), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ex|")
        elif self.graph_axis.currentText() == "Reflected field on y-axis":
            self.canvas.axes.plot_surface(X, Y, abs(E_y_R_), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ey|")
        elif self.graph_axis.currentText() == "Reflected field on z-axis":
            self.canvas.axes.plot_surface(X, Y, abs(E_z_R_), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ez|")
        self.canvas.axes.set_xlabel("X")
        self.canvas.axes.set_ylabel("Y")
        self.canvas.draw()

    def create_graphs_1(self, X, Y, E_x, E_y, E_z):
        self.canvas.axes.cla()
        colormap = 'inferno'
        if self.graph_axis.currentText() == "x-axis field":
            self.canvas.axes.plot_surface(X, Y, abs(E_x), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ex|")
        elif self.graph_axis.currentText() == "y-axis field":
            self.canvas.axes.plot_surface(X, Y, abs(E_y), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ey|")
        elif self.graph_axis.currentText() == "z-axis field":
            self.canvas.axes.plot_surface(X, Y, abs(E_z), cmap=colormap)
            self.canvas.axes.set_zlabel("|Ez|")
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
