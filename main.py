# by Kanegend :)
import math
import cmath
import numpy as np
import gaussian_decomposition
import matplotlib.pyplot as plt
from matrix_method import matrix_method

# CONSTANTS #
c = 3 * 1e10  # cm/s
f = 170 * 1e9  # Gh
w = 0.5  # cm
k0 = (2 * math.pi * f) / c  # cm^-1
lambda_ = c / f  # cm
epsilon_1 = 1  # dielectric permittivity
epsilon_2 = 12
epsilon_3 = 1
d_2 = 0.0675   # layer diameter in cm

# GAUSSIAN BEAM
x = np.arange(-2, 2, 0.01)
y = np.arange(-2, 2, 0.01)
N = 400  # total number of points
T = 0.01  # step
z = 10  # cm
X, Y = np.meshgrid(x, y)
v0 = 1j * (math.pi * (w ** 2)) / lambda_
v = v0 + z
U0 = (1 / v0) * np.exp(-1j * k0 * ((X ** 2 + Y ** 2) / 2) * (1 / v0))
U1 = (1 / v) * np.exp(-1j * k0 * ((X ** 2 + Y ** 2) / 2) * (1 / v))

# DECOMPOSITION OF A GAUSSIAN BEAM INTO A FOURIER SPECTRUM #
spectrum, k_x_y, k_z = gaussian_decomposition.decomposition(U0, N, T, k0)[0], \
                       gaussian_decomposition.decomposition(U0, N, T, k0)[1], \
                       gaussian_decomposition.decomposition(U0, N, T, k0)[2]

# SEARCH FOR TRANSMISSION AND REFRACTION COEFFICIENT #
T_p = np.empty((len(k_x_y), len(k_x_y)), dtype='complex64')
R_p = np.empty((len(k_x_y), len(k_x_y)), dtype='complex64')
T_s = np.empty((len(k_x_y), len(k_x_y)), dtype='complex64')
R_s = np.empty((len(k_x_y), len(k_x_y)), dtype='complex64')
for i in range(len(k_x_y)):
    for j in range(len(k_x_y)):
        if k_z[i][j] == 0:
            alpha = math.pi / 2
        else:
            alpha = cmath.atan(k_x_y[i] / k_z[i][j])
        T_p[i][j] = matrix_method(f, epsilon_1, epsilon_2, epsilon_3, d_2, alpha)[0]
        R_p[i][j] = matrix_method(f, epsilon_1, epsilon_2, epsilon_3, d_2, alpha)[1]
        T_s[i][j] = matrix_method(f, epsilon_1, epsilon_2, epsilon_3, d_2, alpha)[2]
        R_s[i][j] = matrix_method(f, epsilon_1, epsilon_2, epsilon_3, d_2, alpha)[3]
B = spectrum * T_p
C = spectrum * R_p
D = spectrum * T_s
G = spectrum * R_s

# FOLDING THE FOURIER SPECTRUM #
U = gaussian_decomposition.folding(B, N, T, k_z, z)
F = gaussian_decomposition.folding(C, N, T, k_z, z)
H = gaussian_decomposition.folding(D, N, T, k_z, z)
L = gaussian_decomposition.folding(G, N, T, k_z, z)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("U")
ax.plot_surface(X, Y, abs(U1))
ax.plot_surface(X, Y, abs(U))
plt.show()
