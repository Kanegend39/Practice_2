from PyQt6 import QtWidgets
from PyQt6 import uic
import sys
import time
import psutil
import os
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from the_second_task import gaussian_beam_propagation_vector
from Canvas_3D import MplCanvas_3D
from error_table import Error_table


class N_Layers(QtWidgets.QWidget):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.setup_ui()
        self.setup_graphs()
        self.setup_controls()
        self.exit_button.clicked.connect(self.exit)
        self.plot_graphs.clicked.connect(self.graph)
        self.number_of_layers.valueChanged.connect(self.change_number_of_layers)
        self.left_border.valueChanged.connect(self.change_borders)
        self.right_border.valueChanged.connect(self.change_borders)
        self.plot_section.clicked.connect(self.section)
        self.x_section.toggled.connect(self.on_x_section)
        self.y_section.toggled.connect(self.on_y_section)
        self.export_section.clicked.connect(self.export_2d)
        self.button_export_3d.clicked.connect(self.export_3d)
        self.button_add_data.clicked.connect(self.add_data)
        self.number_of_file = 1

    def setup_ui(self):
        try:
            uic.loadUi("n_layers.ui", self)
        except Exception as e:
            print(f"Error loading UI: {e}")
            sys.exit(1)

    def setup_controls(self):
        self.x_const_section.setRange(self.left_border.value(), self.right_border.value())
        self.y_const_section.setRange(self.left_border.value(), self.right_border.value())
        self.y_const_section.setEnabled(False)
        self.export_section.setEnabled(False)
        self.button_export_3d.setEnabled(False)

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

    def change_number_of_layers(self):
        self.table.setRowCount(self.number_of_layers.value())

    def change_borders(self):
        self.x_const_section.setRange(self.left_border.value(), self.right_border.value())
        self.y_const_section.setRange(self.left_border.value(), self.right_border.value())

    def read_data(self):
        left = self.left_border.value()
        right = self.right_border.value()
        step = self.step.value()
        z = self.distance.value()
        w = self.wave_waist.value()
        n = self.number_of_layers.value()
        d_i = np.empty(self.table.rowCount(), dtype="float64")
        epsilons = np.empty(self.table.rowCount() * 2 + 1, dtype="complex64")
        epsilons[0] = 1 + 0j
        for i in range(self.table.rowCount()):
            d_i[i] = float(self.table.item(i, 0).text().replace(",", "."))
            epsilons[i * 2 + 1] = float(self.table.item(i, 1).text().replace(",", ".")) + 1j * float(
                self.table.item(i, 2).text().replace(",", "."))
            epsilons[(i + 1) * 2] = 1 + 0j
        return gaussian_beam_propagation_vector(left, right, step, z, w, d_i, epsilons, n)

    def graph(self):
        try:
            start_time = time.time()
            process = psutil.Process(os.getpid())
            memory_usage_before = process.memory_info().rss
            X, Y, E_x, E_y, E_z, E_x_T_, E_y_T_, E_z_T_, E_x_R_, E_y_R_, E_z_R = self.read_data()
            colormap = 'inferno'
            z_axis = ["|Ex|", "|Ey|", "|Ez|"]
            graphs_transmitted = [abs(E_x_T_), abs(E_y_T_), abs(E_z_T_)]
            for i in range(len(self.canvases)):
                canvas = self.canvases[i]
                canvas.axes.cla()
                canvas.axes.plot_surface(X, Y, graphs_transmitted[i], cmap=colormap)
                canvas.axes.set_zlabel(z_axis[i])
                canvas.axes.set_xlabel("X")
                canvas.axes.set_ylabel("Y")
                canvas.draw()
            self.button_export_3d.setEnabled(True)
            memory_usage_after = process.memory_info().rss
            end_time = time.time()
            execution_time = end_time - start_time
            print(f"Время выполнения метода run: {execution_time} секунд")
            print(f"Использование памяти: {(memory_usage_after - memory_usage_before) / (1024 * 1024)} МБ")
        except Exception:
            self.error_table_window = Error_table()
            self.error_table_window.show()

    def find_section(self):
        X, Y, E_x, E_y, E_z, E_x_T_, E_y_T_, E_z_T_, E_x_R_, E_y_R_, E_z_R = self.read_data()
        if self.x_section.isChecked():
            y = np.empty(len(X))
            E_T_1 = np.empty(len(X), dtype='complex64')
            E_T_2 = np.empty(len(X), dtype='complex64')
            E_T_3 = np.empty(len(X), dtype='complex64')
            for i in range(len(X)):
                for j in range(len(X)):
                    if abs(X[i][j] - self.x_const_section.value()) < 1e-9:
                        y[i] = Y[i][j]
                        E_T_1[i] = E_x_T_[i][j]
                        E_T_2[i] = E_y_T_[i][j]
                        E_T_3[i] = E_z_T_[i][j]
            return y, E_T_1, E_T_2, E_T_3
        elif self.y_section.isChecked():
            for i in range(len(Y)):
                if abs(Y[i][i] - self.y_const_section.value()) < 1e-9:
                    return X[i], E_x_T_[i], E_y_T_[i], E_z_T_[i]

    def section(self):
        start_time = time.time()
        process = psutil.Process(os.getpid())
        memory_usage_before = process.memory_info().rss
        plt.close('all')
        if self.x_section.isChecked():
            try:
                y, E_T_1, E_T_2, E_T_3 = self.find_section()
                self.export_section.setEnabled(True)
                plt.figure()
                plt.plot(y, abs(E_T_1))
                plt.grid(True)
                plt.xlabel('Y')
                plt.ylabel('|Ex|')
                plt.figure()
                plt.plot(y, abs(E_T_2))
                plt.grid(True)
                plt.xlabel('Y')
                plt.ylabel('|Ey|')
                plt.figure()
                plt.plot(y, abs(E_T_3))
                plt.grid(True)
                plt.xlabel('Y')
                plt.ylabel('|Ez|')
                plt.show()
            except Exception:
                self.error_table_window = Error_table()
                self.error_table_window.show()
        elif self.y_section.isChecked():
            try:
                x, E_T_1, E_T_2, E_T_3 = self.find_section()
                self.export_section.setEnabled(True)
                plt.figure()
                plt.plot(x, abs(E_T_1))
                plt.grid(True)
                plt.xlabel('X')
                plt.ylabel('|Ex|')
                plt.figure()
                plt.plot(x, abs(E_T_2))
                plt.grid(True)
                plt.xlabel('X')
                plt.ylabel('|Ey|')
                plt.figure()
                plt.plot(x, abs(E_T_3))
                plt.grid(True)
                plt.xlabel('X')
                plt.ylabel('|Ez|')
                plt.show()
            except Exception:
                self.error_table_window = Error_table()
                self.error_table_window.show()
        memory_usage_after = process.memory_info().rss
        end_time = time.time()
        execution_time = end_time - start_time
        print(f"Время выполнения метода run: {execution_time} секунд")
        print(f"Использование памяти: {(memory_usage_after - memory_usage_before) / (1024 * 1024)} МБ")

    def on_x_section(self):
        self.y_const_section.setEnabled(False)
        self.x_const_section.setEnabled(True)

    def on_y_section(self):
        self.y_const_section.setEnabled(True)
        self.x_const_section.setEnabled(False)

    def export_2d(self):
        try:
            start_time = time.time()
            process = psutil.Process(os.getpid())
            memory_usage_before = process.memory_info().rss
            A, B, C, D = self.find_section()
            if self.x_section.isChecked():
                name_1 = f"x={self.x_const_section.value()}_section_{self.number_of_file}.txt"
                name_2 = f"README_x={self.x_const_section.value()}_section_{self.number_of_file}.txt"
                my_file_1 = open(name_1, "w+")
                my_file_2 = open(name_2, "w+")
                my_file_2.write(f"left_border={self.left_border.value()}, right_border={self.right_border.value()}, "
                                f"step={self.step.value()}, distance={self.distance.value()}, wave_waist={self.wave_waist.value()}, "
                                f"number_of_layers={self.number_of_layers.value()}\n")
                for i in range(self.table.rowCount()):
                    my_file_2.write(f"thickness={self.table.item(i, 0).text().replace(',', '.')}, epsilon={float(self.table.item(i, 1).text().replace(',', '.')) + 1j * float(self.table.item(i, 2).text().replace(',', '.'))}\n")
                my_file_2.write("y  Re(Ex)  Im(Ex)  Re(Ey)  Im(Ey)  Re(Ez)  Im(Ez)\n")
                for i in range(len(A)):
                    my_file_1.write(
                        f"{A[i]}  {B[i].real}  {B[i].imag}  {C[i].real}  {C[i].imag}  {D[i].real}  {D[i].imag}\n")
                my_file_1.close()
                my_file_2.close()
            elif self.y_section.isChecked():
                name_1 = f"y={self.y_const_section.value()}_section_{self.number_of_file}.txt"
                name_2 = f"README_y={self.y_const_section.value()}_section_{self.number_of_file}.txt"
                my_file_1 = open(name_1, "w+")
                my_file_2 = open(name_2, "w+")
                my_file_2.write(f"left_border={self.left_border.value()}, right_border={self.right_border.value()}, "
                                f"step={self.step.value()}, distance={self.distance.value()}, wave_waist={self.wave_waist.value()}, "
                                f"number_of_layers={self.number_of_layers.value()}\n")
                for i in range(self.table.rowCount()):
                    my_file_2.write(
                        f"thickness={self.table.item(i, 0).text().replace(',', '.')}, epsilon={float(self.table.item(i, 1).text().replace(',', '.')) + 1j * float(self.table.item(i, 2).text().replace(',', '.'))}\n")
                my_file_2.write("x  Re(Ex)  Im(Ex)  Re(Ey)  Im(Ey)  Re(Ez)  Im(Ez)\n")
                for i in range(len(A)):
                    my_file_1.write(
                        f"{A[i]}  {B[i].real}  {B[i].imag}  {C[i].real}  {C[i].imag}  {D[i].real}  {D[i].imag}\n")
                my_file_1.close()
                my_file_2.close()
            self.number_of_file += 1
            memory_usage_after = process.memory_info().rss
            end_time = time.time()
            execution_time = end_time - start_time
            print(f"Время выполнения метода run: {execution_time} секунд")
            print(f"Использование памяти: {(memory_usage_after - memory_usage_before) / (1024 * 1024)} МБ")
        except Exception:
            self.error_table_window = Error_table()
            self.error_table_window.show()

    def export_3d(self):
        try:
            start_time = time.time()
            process = psutil.Process(os.getpid())
            memory_usage_before = process.memory_info().rss
            X, Y, E_x, E_y, E_z, E_x_T_, E_y_T_, E_z_T_, E_x_R_, E_y_R_, E_z_R = self.read_data()
            name_1 = f"graph_3d_{self.number_of_file}.txt"
            name_2 = f"README_graph_3d_{self.number_of_file}.txt"
            my_file_1 = open(name_1, "w+")
            my_file_2 = open(name_2, "w+")
            my_file_2.write(f"left_border={self.left_border.value()}, right_border={self.right_border.value()}, "
                            f"step={self.step.value()}, distance={self.distance.value()}, wave_waist={self.wave_waist.value()}, "
                            f"number_of_layers={self.number_of_layers.value()}\n")
            for i in range(self.table.rowCount()):
                my_file_2.write(
                    f"thickness={self.table.item(i, 0).text().replace(',', '.')}, epsilon={float(self.table.item(i, 1).text().replace(',', '.')) + 1j * float(self.table.item(i, 2).text().replace(',', '.'))}\n")
            my_file_2.write("x  y  Re(Ex)  Im(Ex)  Re(Ey)  Im(Ey)  Re(Ez)  Im(Ez)\n")
            for i in range(len(X)):
                for j in range(len(X)):
                    my_file_1.write(
                        f"{X[i][j]}  {Y[i][j]}  {E_x_T_[i][j].real}  {E_x_T_[i][j].imag}  {E_y_T_[i][j].real}  {E_y_T_[i][j].imag} "
                        f"{E_z_T_[i][j].real}  {E_z_T_[i][j].imag}\n")
            my_file_1.close()
            my_file_2.close()
            self.number_of_file += 1
            memory_usage_after = process.memory_info().rss
            end_time = time.time()
            execution_time = end_time - start_time
            print(f"Время выполнения метода run: {execution_time} секунд")
            print(f"Использование памяти: {(memory_usage_after - memory_usage_before) / (1024 * 1024)} МБ")
        except Exception:
            self.error_table_window = Error_table()
            self.error_table_window.show()

    def add_data(self):
        H = self.plate_thickness_h.value()
        L = self.layer_thickness_l.value()
        n = self.number_of_layers_n.value()
        epsilon_0_shtrih = self.epsilon_0_shtrih.value()
        epsilon_0_shtrih_shtrih = self.epsilon_0_shtrih_shtrih.value()
        epsilon_1 = self.epsilon_1.value()
        epsilon_2 = self.epsilon_2.value()
        self.number_of_layers.setValue(n + 1)
        for i in range(self.table.rowCount()):
            self.table.setItem(i, 0, QtWidgets.QTableWidgetItem(str(L / n)))
            self.table.setItem(i, 1, QtWidgets.QTableWidgetItem(str(epsilon_0_shtrih + epsilon_1 * np.exp(- (i * (L / n)) / (L / n)))))
            self.table.setItem(i, 2, QtWidgets.QTableWidgetItem(str(epsilon_0_shtrih_shtrih + epsilon_2 * np.exp(- (i * (L / n)) / (L / n)))))
        self.table.setItem(n, 0, QtWidgets.QTableWidgetItem(str(H - L)))
        self.table.setItem(n, 1, QtWidgets.QTableWidgetItem(str(12.)))
        self.table.setItem(n, 2, QtWidgets.QTableWidgetItem(str(0.)))

    def exit(self):
        self.close()
        self.parent.show()
