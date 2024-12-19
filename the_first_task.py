# Propagation of a Gaussian beam in free space (scalar problem) #

import numpy as np
from scipy.fft import fft2, fftfreq, ifft2


def gaussian_beam_propagation_no_vector(left, right, step, z, w):
    # CONSTANTS #
    c = 3e10  # cm/s
    f = 170e9  # Gh
    k0 = (2 * np.pi * f) / c  # cm^-1
    x = np.arange(left, right, step)
    y = np.arange(left, right, step)
    X, Y = np.meshgrid(x, y)
    v0 = -1j * k0 * (w ** 2)
    U0_x = (1 / v0) * (-1j * k0 * (w ** 2)) * np.exp(1j * k0 * ((X ** 2 + Y ** 2) / 2) * (1 / v0))
    U0_y = X * 0
    U0_z = (1 / (v0 ** 2)) * (1j * X * k0 * (w ** 2)) * np.exp(1j * k0 * ((X ** 2 + Y ** 2) / 2) * (1 / v0))

    # DECOMPOSITION OF A GAUSSIAN BEAM INTO A FOURIER SPECTRUM #
    F_U_x = fft2(U0_x)
    F_U_y = fft2(U0_y)
    F_U_z = fft2(U0_z)
    fx = fftfreq(len(x), d=x[1] - x[0])
    fy = fftfreq(len(y), d=y[1] - y[0])
    FX, FY = np.meshgrid(fx, fy)
    k_x = 2 * np.pi * FX
    k_y = 2 * np.pi * FY
    k_z = np.where(k0 ** 2 >= k_x ** 2 + k_y ** 2, np.sqrt(k0 ** 2 - k_x ** 2 - k_y ** 2, dtype='complex64'),
                   1j * np.sqrt(k_x ** 2 + k_y ** 2 - k0 ** 2, dtype='complex64'))

    # GAUSSIAN BEAM AFTER THE INVERSE FOURIER TRANSFORM#
    E_x = ifft2(F_U_x * np.exp(1j * k_z * z))
    E_y = ifft2(F_U_y * np.exp(1j * k_z * z))
    E_z = ifft2(F_U_z * np.exp(1j * k_z * z))
    return X, Y, E_x, E_y, E_z