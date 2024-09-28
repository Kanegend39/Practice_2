import math
import numpy as np
from scipy.fft import fft2, fftfreq, ifft2


def decomposition(U0, N, T, k):
    k_x_y = fftfreq(N, T)
    k_z = np.empty((len(k_x_y), len(k_x_y)), dtype='complex64')
    E = fft2(U0)
    for i in range(len(k_x_y)):
        for j in range(len(k_x_y)):
            if k ** 2 >= k_x_y[i] ** 2 + k_x_y[j] ** 2:
                k_z[i][j] = math.sqrt(k ** 2 - (k_x_y[j] ** 2) - (k_x_y[i] ** 2))
            else:
                k_z[i][j] = 1j * math.sqrt(- (k ** 2) + (k_x_y[j] ** 2) + (k_x_y[i] ** 2))
    return E, k_x_y, k_z


def folding(J, N, T, k_z, z):
    R = J * np.exp(1j * N * T * (math.pi / 5) * math.pi * k_z * z)
    return ifft2(R)
