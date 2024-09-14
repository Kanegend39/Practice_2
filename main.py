# by Kanegend
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, fftfreq, ifftn, ifft2


def grafik_abs_U_X_Y(x, y, k, q, z):
    X, Y = np.meshgrid(x, y)
    q += z
    U = (1 / q) * np.exp(-1j * k * ((X ** 2 + Y ** 2) / 2) * (1 / q))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, abs(U))
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("U")


def grafik_kx_ky_F(x, y, k, q, z):
    X, Y = np.meshgrid(x, y)
    N = 2048
    T = 0.01
    U = (1 / q) * np.exp(-1j * k * ((X ** 2 + Y ** 2) / 2) * (1 / q))
    k_x_y = fftfreq(N, T)
    F = fft2(U)
    A, B = np.meshgrid(k_x_y, k_x_y)
    # ifft2 doesnt work, need algorithm for fourier transform with sigma
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.plot_surface(abs(k_x_y), abs(k_x_y), 1.0 / N * np.abs(F))
    # ax.set_xlabel("|kx|")
    # ax.set_ylabel("|ky|")
    # ax.set_zlabel("Fourie")
    ax.plot_surface(X, Y, np.abs(ifft2(F)))
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("U")


def main():
    f = 170 * 1e9  # в герцах
    c = 3 * 1e10  # см/с
    k = (2 * math.pi * f) / c
    w = 0.5  # при z = 0 в см
    lambda_ = (2 * math.pi) / k  # в см
    x = np.arange(-10, 10.48, 0.01)
    y = np.arange(-10, 10.48, 0.01)
    z = 10  # в см
    q = 1j * ((math.pi * (w ** 2)) / lambda_)  # R(z=0) → ∞
    grafik_kx_ky_F(x, y, k, q, z)
    grafik_abs_U_X_Y(x, y, k, q, z)
    plt.show()


main()
