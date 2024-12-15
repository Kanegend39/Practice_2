# Propagation of a Gaussian beam at an arbitrary angle to the plate normal #
# with considering the vector nature of the field #

import numpy as np
from scipy.fft import fft2, fftfreq, ifft2
from matrix_method import matrix_method


def gaussian_beam_propagation_vector(left, right, step, z, w, d_2, epsilon_1=1, epsilon_2=12, epsilon_3=1):
    # CONSTANTS #
    c = 3e10  # cm/s
    f = 170e9  # Gh
    k0 = (2 * np.pi * f) / c  # cm^-1
    x = np.arange(left, right, step)
    y = np.arange(left, right, step)
    X, Y = np.meshgrid(x, y)
    v0 = -1j * k0 * (w ** 2)
    E_x = (1 / v0) * (-1j * k0 * (w ** 2)) * np.exp(1j * k0 * ((X ** 2 + Y ** 2) / 2) * (1 / v0))
    E_y = 0 * X
    E_z = (1 / (v0 ** 2)) * (1j * X * k0 * (w ** 2)) * np.exp(1j * k0 * ((X ** 2 + Y ** 2) / 2) * (1 / v0))

    # DECOMPOSITION OF A GAUSSIAN BEAM INTO A FOURIER SPECTRUM #
    fx = fftfreq(len(x), d=x[1] - x[0])
    fy = fftfreq(len(y), d=y[1] - y[0])
    FX, FY = np.meshgrid(fx, fy)
    k_x = 2 * np.pi * FX
    k_y = 2 * np.pi * FY
    k_z = np.where(k0 ** 2 >= k_x ** 2 + k_y ** 2, np.sqrt(k0 ** 2 - k_x ** 2 - k_y ** 2, dtype='complex64'),
                   1j * np.sqrt(k_x ** 2 + k_y ** 2 - k0 ** 2, dtype='complex64'))
    unit_vector_E_s_x = np.empty((len(k_x), len(k_x)), dtype="complex64")
    unit_vector_E_s_y = np.empty((len(k_x), len(k_x)), dtype="complex64")
    unit_vector_E_s_z = np.empty((len(k_x), len(k_x)), dtype="complex64")
    for i in range(len(k_x)):
        for j in range(len(k_x)):
            if k_x[i][j] ** 2 + k_y[i][j] ** 2 == 0:
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
            if k0 ** 2 >= k_x[i][j] ** 2 + k_y[i][j] ** 2:
                unit_vector_E_p_x[i][j] = (unit_vector_E_s_y[i][j] * k_z[i][j]) / k0
                unit_vector_E_p_y[i][j] = - (unit_vector_E_s_x[i][j] * k_z[i][j]) / k0
                unit_vector_E_p_z[i][j] = (unit_vector_E_s_x[i][j] * k_y[i][j] - unit_vector_E_s_y[i][j] * k_x[i][
                    j]) / k0
            else:
                unit_vector_E_p_x[i][j] = 0
                unit_vector_E_p_y[i][j] = 0
                unit_vector_E_p_z[i][j] = 1
    C_s = fft2(E_x) * unit_vector_E_s_x + fft2(E_y) * unit_vector_E_s_y + fft2(E_z) * unit_vector_E_s_z
    C_p = fft2(E_x) * unit_vector_E_p_x + fft2(E_y) * unit_vector_E_p_y + fft2(E_z) * unit_vector_E_p_z

    # SEARCH FOR TRANSMISSION AND REFRACTION COEFFICIENT #
    T_p = np.empty(k_x.shape, dtype='complex64')
    R_p = np.empty(k_x.shape, dtype='complex64')
    T_s = np.empty(k_x.shape, dtype='complex64')
    R_s = np.empty(k_x.shape, dtype='complex64')
    for i in range(len(k_x)):
        for j in range(len(k_y)):
            alpha = k_x[i][j] / k0
            result = matrix_method(f, epsilon_1, epsilon_2, epsilon_3, d_2, alpha)
            T_p[i][j], R_p[i][j], T_s[i][j], R_s[i][j] = result

    # GAUSSIAN BEAM AFTER THE INVERSE FOURIER TRANSFORM CONSIDERING ITS PASSAGE THROUGH THE PLATE #
    E_x = ifft2(C_p * unit_vector_E_p_x + C_s * unit_vector_E_s_x)
    E_y = ifft2(C_p * unit_vector_E_p_y + C_s * unit_vector_E_s_y)
    E_z = ifft2(C_p * unit_vector_E_p_z + C_s * unit_vector_E_s_z)
    C_s_T = T_s * C_s * np.exp(1j * k_z * z)
    C_p_T = T_p * C_p * np.exp(1j * k_z * z)
    E_x_T_ = ifft2(C_p_T * unit_vector_E_p_x + C_s_T * unit_vector_E_s_x)
    E_y_T_ = ifft2(C_p_T * unit_vector_E_p_y + C_s_T * unit_vector_E_s_y)
    E_z_T_ = ifft2(C_p_T * unit_vector_E_p_z + C_s_T * unit_vector_E_s_z)
    C_s_R = R_s * C_s * np.exp(1j * k_z * z)
    C_p_R = R_p * C_p * np.exp(1j * k_z * z)
    E_x_R_ = ifft2(C_p_R * unit_vector_E_p_x + C_s_R * unit_vector_E_s_x)
    E_y_R_ = ifft2(C_p_R * unit_vector_E_p_y + C_s_R * unit_vector_E_s_y)
    E_z_R_ = ifft2(C_p_R * unit_vector_E_p_z + C_s_R * unit_vector_E_s_z)
    return X, Y, E_x, E_y, E_z, E_x_T_, E_y_T_, E_z_T_, E_x_R_, E_y_R_, E_z_R_
