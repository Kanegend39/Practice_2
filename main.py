# by Kanegend :)
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, fftfreq, ifft2

f = 170 * 1e9
c = 3 * 1e10
k = (2 * math.pi * f) / c
lambda_ = c / f
a = 0.5
v0 = 1j * (math.pi * (a ** 2)) / lambda_
z = 10
v = v0 + z
x = np.arange(-10, 10, 0.01)
y = np.arange(-10, 10, 0.01)
X, Y = np.meshgrid(x, y)
U0 = (1 / v0) * np.exp(-1j * k * ((X ** 2 + Y ** 2) / 2) * (1 / v0))
U = (1 / v) * np.exp(-1j * k * ((X ** 2 + Y ** 2) / 2) * (1 / v))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, abs(U0))
ax.plot_surface(X, Y, abs(U))
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("U")
N = 2000
T = 0.01
k_x_y = fftfreq(N, T)
k_z = np.empty((len(k_x_y), len(k_x_y)), dtype='complex64')
for i in range(len(k_x_y)):
    for j in range(len(k_x_y)):
        if k ** 2 >= k_x_y[i] ** 2 + k_x_y[j] ** 2:
            k_z[i][j] = math.sqrt(k ** 2 - (k_x_y[j] ** 2) - (k_x_y[i] ** 2))
        else:
            k_z[i][j] = 1j * math.sqrt(- (k ** 2) + (k_x_y[j] ** 2) + (k_x_y[i] ** 2))
E = np.exp(1j * z * k_z) * fft2(U)
ax.plot_surface(X, Y, abs(ifft2(E)))
plt.show()
