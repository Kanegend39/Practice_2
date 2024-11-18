# by Kanegend :)
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft2, fftfreq, ifft2

import gaussian_decomposition
from matrix_method import matrix_method


def gauss(left, right, step, z, w):
    # CONSTANTS #
    c = 3 * 1e10  # cm/s
    f = 170 * 1e9  # Gh
    k0 = (2 * math.pi * f) / c  # cm^-1
    lambda_ = c / f  # cm
    epsilon_1 = 1  # dielectric permittivity
    epsilon_2 = 12
    epsilon_3 = 1
    tetta = (math.pi / 2) - 1.31  # angle between the horizontal axis and the plate
    d_2 = 0.0675  # layer diameter in cm
    # GAUSSIAN BEAM
    x = np.arange(left, right, step)
    y = np.arange(left, right, step)
    N = (right - left) / step
    X, Y = np.meshgrid(x, y)
    v0 = 1j * (math.pi * (w ** 2)) / lambda_
    E_x = (1 / v0) * np.exp(-1j * k0 * ((X ** 2 + Y ** 2) / 2) * (1 / v0))
    E_y = 0 * X
    E_z = (X / (v0 ** 2)) * np.exp(-1j * k0 * ((X ** 2 + Y ** 2) / 2) * (1 / v0))

    # DECOMPOSITION OF A GAUSSIAN BEAM INTO A FOURIER SPECTRUM #
    gd = gaussian_decomposition.decomposition(E_x, x, y, k0)
    k_x, k_y, k_z = gd[1], gd[2], gd[3]
    unit_vector_E_s_x = np.empty((len(k_x), len(k_x)), dtype="complex64")
    unit_vector_E_s_y = np.empty((len(k_x), len(k_x)), dtype="complex64")
    unit_vector_E_s_z = np.empty((len(k_x), len(k_x)), dtype="complex64")
    for i in range(len(k_x)):
        for j in range(len(k_x)):
            if k_x[i][j] == k_y[i][j] == 0:
                unit_vector_E_s_x[i][j] = 1
                unit_vector_E_s_y[i][j] = 0
                unit_vector_E_s_z[i][j] = 0
            else:
                unit_vector_E_s_x[i][j] = (k_y[i][j] / np.sqrt(abs(k_x[i][j]) ** 2 + abs(k_y[i][j]) ** 2))
                unit_vector_E_s_y[i][j] = - (k_x[i][j] / np.sqrt(abs(k_x[i][j]) ** 2 + abs(k_y[i][j]) ** 2))
                unit_vector_E_s_z[i][j] = 0
    unit_vector_E_p_x = np.empty((len(k_x), len(k_x)), dtype="complex64")
    unit_vector_E_p_y = np.empty((len(k_x), len(k_x)), dtype="complex64")
    unit_vector_E_p_z = np.empty((len(k_x), len(k_x)), dtype="complex64")
    for i in range(len(k_x)):
        for j in range(len(k_x)):
            if k_x[i][j] == k_y[i][j] == 0:
                unit_vector_E_p_x[i][j] = 0
                unit_vector_E_p_y[i][j] = 1
                unit_vector_E_p_z[i][j] = 0
            else:
                unit_vector_E_p_x[i][j] = (unit_vector_E_s_y[i][j] * k_z[i][j]) / \
                                          np.sqrt(abs((unit_vector_E_s_y[i][j] * k_z[i][j])) ** 2 +
                                                  abs((unit_vector_E_s_x[i][j] * k_z[i][j])) ** 2 +
                                                  abs((unit_vector_E_s_x[i][j] * k_y[i][j] - unit_vector_E_s_y[i][j] *
                                                       k_x[i][j])) ** 2)
                unit_vector_E_p_y[i][j] = - (unit_vector_E_s_x[i][j] * k_z[i][j]) / \
                                          np.sqrt(abs((unit_vector_E_s_y[i][j] * k_z[i][j])) ** 2 +
                                                  abs((unit_vector_E_s_x[i][j] * k_z[i][j])) ** 2 +
                                                  abs((unit_vector_E_s_x[i][j] * k_y[i][j] - unit_vector_E_s_y[i][j] *
                                                       k_x[i][j])) ** 2)
                unit_vector_E_p_z[i][j] = (unit_vector_E_s_x[i][j] * k_y[i][j] - unit_vector_E_s_y[i][j] * k_x[i][j]) / \
                                          np.sqrt( abs((unit_vector_E_s_y[i][j] * k_z[i][j])) ** 2 +
                                                  abs((unit_vector_E_s_x[i][j] * k_z[i][j])) ** 2 +
                                                  abs((unit_vector_E_s_x[i][j] * k_y[i][j] - unit_vector_E_s_y[i][j] *
                                                       k_x[i][j])) ** 2)
    D = abs(abs(unit_vector_E_s_x) ** 2 + abs(unit_vector_E_s_y) ** 2 + abs(unit_vector_E_s_z) ** 2)
    E = abs(abs(unit_vector_E_p_x) ** 2 + abs(unit_vector_E_p_y) ** 2 + abs(unit_vector_E_p_z) ** 2)
    H = abs(
        unit_vector_E_s_x * unit_vector_E_p_x + unit_vector_E_s_y * unit_vector_E_p_y + unit_vector_E_s_z * unit_vector_E_p_z)
    C_s = fft2(
        E_x * unit_vector_E_s_x.conjugate() + E_y * unit_vector_E_s_y.conjugate() + E_z * unit_vector_E_s_z.conjugate()) / (
                      4 * math.pi ** 2)
    C_p = fft2(
        E_x * unit_vector_E_p_x.conjugate() + E_y * unit_vector_E_p_y.conjugate() + E_z * unit_vector_E_p_z.conjugate()) / (
                      4 * math.pi ** 2)
    # SEARCH FOR TRANSMISSION AND REFRACTION COEFFICIENT #
    T_p = np.empty((len(k_z), len(k_z)), dtype='complex64')
    R_p = np.empty((len(k_z), len(k_z)), dtype='complex64')
    T_s = np.empty((len(k_z), len(k_z)), dtype='complex64')
    R_s = np.empty((len(k_z), len(k_z)), dtype='complex64')
    n = np.array([-math.cos(tetta), 0, -math.sin(tetta)])
    for i in range(len(k_x)):
        for j in range(len(k_y)):
            if tetta != 0:
                cos_alpha = (n[0] * k_x[i][j] + n[1] * k_y[i][j] + n[2] * k_z[i][j]) / (
                        np.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2) * np.sqrt(
                    k_x[i][j] ** 2 + k_y[i][j] ** 2 + k_z[i][j] ** 2))
                alpha = np.sqrt(1 - cos_alpha ** 2)
            else:
                alpha = k_x[i][j] / k0
            T_p[i][j] = matrix_method(f, epsilon_1, epsilon_2, epsilon_3, d_2, alpha)[0]
            R_p[i][j] = matrix_method(f, epsilon_1, epsilon_2, epsilon_3, d_2, alpha)[1]
            T_s[i][j] = matrix_method(f, epsilon_1, epsilon_2, epsilon_3, d_2, alpha)[2]
            R_s[i][j] = matrix_method(f, epsilon_1, epsilon_2, epsilon_3, d_2, alpha)[3]
    E_x_T = ifft2((T_p * C_p * unit_vector_E_p_x + T_s * C_s * unit_vector_E_s_x) * np.exp(1j * k_z * z))
    E_x_R = ifft2((R_p * C_p * unit_vector_E_p_x + R_s * C_s * unit_vector_E_s_x) * np.exp(1j * k_z * z))
    E_x_ = ifft2(C_p * unit_vector_E_p_x + C_s * unit_vector_E_s_x)
    E_y_ = ifft2(C_p * unit_vector_E_p_y + C_s * unit_vector_E_s_y)
    E_z_ = ifft2(C_p * unit_vector_E_p_z + C_s * unit_vector_E_s_z)
    E_z_T = ifft2(C_p * unit_vector_E_p_z + C_s * unit_vector_E_s_z)
    # FOLDING THE FOURIER SPECTRUM #
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.set_xlabel("X")
    # ax.set_ylabel("Y")
    # ax.set_zlabel("|Ex|")
    # ax.plot_surface(X, Y, abs(E_x_), color='orange')
    # ax.plot_surface(X, Y, abs(E_x_), color='blue')
    # ax.plot_surface(X, Y, abs(E_y_), color='green')
    # ax.plot_surface(X, Y, abs(E_x_), color='yellow')
    # ax.plot_surface(X, Y, abs(E_x_), color='blue')
    # ax.plot_surface(X, Y, abs(E_z), color='green')
    # ax.plot_surface(X, Y, abs(E_x), color='green')
    # plt.show()
    return X, Y, E_x, E_y, E_z, E_x_, E_y_, E_z_, k_x, k_y
