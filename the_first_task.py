# Propagation of a Gaussian beam at an arbitrary angle to the plate normal #
# without considering the vector nature of the field #

import numpy as np
from scipy.fft import fft2, fftfreq, ifft2
from matrix_method import matrix_method


def gaussian_beam_propagation_no_vector(left, right, step, z, w, tetta, d_2, epsilon_1=1, epsilon_2=12, epsilon_3=1):
    # CONSTANTS #
    c = 3e10  # cm/s
    f = 170e9  # Gh
    k0 = (2 * np.pi * f) / c  # cm^-1
    x = np.arange(left, right, step)
    y = np.arange(left, right, step)
    X, Y = np.meshgrid(x, y)
    v0 = -1j * k0 * (w ** 2)
    U0_x_y = (1 / v0) * (-1j * k0 * (w ** 2)) * np.exp(1j * k0 * ((X ** 2 + Y ** 2) / 2) * (1 / v0))

    # DECOMPOSITION OF A GAUSSIAN BEAM INTO A FOURIER SPECTRUM #
    E_kx_ky = fft2(U0_x_y)
    fx = fftfreq(len(x), d=x[1] - x[0])
    fy = fftfreq(len(y), d=y[1] - y[0])
    FX, FY = np.meshgrid(fx, fy)
    k_x = 2 * np.pi * FX
    k_y = 2 * np.pi * FY
    k_z = np.where(k0 ** 2 >= k_x ** 2 + k_y ** 2, np.sqrt(k0 ** 2 - k_x ** 2 - k_y ** 2, dtype='complex64'),
                   1j * np.sqrt(k_x ** 2 + k_y ** 2 - k0 ** 2, dtype='complex64'))

    # SEARCH FOR TRANSMISSION AND REFRACTION COEFFICIENT #
    T_p = np.empty(k_x.shape, dtype='complex64')
    R_p = np.empty(k_x.shape, dtype='complex64')
    T_s = np.empty(k_x.shape, dtype='complex64')
    R_s = np.empty(k_x.shape, dtype='complex64')
    n = np.array([-np.cos(tetta), 0, -np.sin(tetta)])
    for i in range(len(k_x)):
        for j in range(len(k_y)):
            if not np.isclose(tetta, 0):
                cos_alpha = (n[0] * k_x[i][j] + n[1] * k_y[i][j] + n[2] * k_z[i][j]) / (
                        np.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2) * np.sqrt(
                    k_x[i][j] ** 2 + k_y[i][j] ** 2 + k_z[i][j] ** 2))
                alpha = np.sqrt(1 - cos_alpha ** 2)
            else:
                alpha = k_x[i][j] / k0
            result = matrix_method(f, epsilon_1, epsilon_2, epsilon_3, d_2, alpha)
            T_p[i][j], R_p[i][j], T_s[i][j], R_s[i][j] = result

    # GAUSSIAN BEAM AFTER THE INVERSE FOURIER TRANSFORM CONSIDERING ITS PASSAGE THROUGH THE PLATE #
    U_x_y_T_p = ifft2(E_kx_ky * T_p * np.exp(1j * k_z * z))
    U_x_y_T_s = ifft2(E_kx_ky * T_s * np.exp(1j * k_z * z))
    U_x_y_R_p = ifft2(E_kx_ky * R_p * np.exp(1j * k_z * z))
    U_x_y_R_s = ifft2(E_kx_ky * R_s * np.exp(1j * k_z * z))
    return X, Y, U0_x_y, U_x_y_T_p, U_x_y_T_s, U_x_y_R_p, U_x_y_R_s
