import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, fftfreq, ifft2


def decomposition(U0, x, y, k):
    fx = fftfreq(len(x), d=x[1] - x[0])
    fy = fftfreq(len(y), d=y[1] - y[0])
    FX, FY = np.meshgrid(fx, fy)
    E = fft2(U0)
    k_x = 2 * math.pi * FX
    k_y = 2 * math.pi * FY
    k_z = np.where(k ** 2 >= k_x ** 2 + k_y ** 2, np.sqrt(k ** 2 - k_x ** 2 - k_y ** 2, dtype='complex64'),
                   1j * np.sqrt(k_x ** 2 + k_y ** 2 - k ** 2, dtype='complex64'))
    return E, k_x, k_y, k_z


def folding(J, k_z, z):
    R = J * np.exp(1j * k_z * z)
    return ifft2(R)
