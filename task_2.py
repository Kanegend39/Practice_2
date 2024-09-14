import math
import matplotlib.pyplot as plt
import numpy as np

f = 170 * 10e9  # GH
epsilon_1 = 1
epsilon_2 = 12
epsilon_3 = 1
n_1 = np.sqrt(epsilon_1)
n_2 = np.sqrt(epsilon_2)
n_3 = np.sqrt(epsilon_3)
c = 3 * 10e8
d_2 = 0.675 * 10e-3  # meters
k_0 = (2 * math.pi * f) / c
alpha = np.arange(0., 3.14 / 2, 0.01)
k_x = k_0 * np.sin(alpha)
N_1 = np.sqrt(epsilon_1 - ((k_x ** 2) / (k_0 ** 2)))
N_2 = np.sqrt(epsilon_2 - ((k_x ** 2) / (k_0 ** 2)))
N_3 = np.sqrt(epsilon_3 - ((k_x ** 2) / (k_0 ** 2)))
phi_2 = np.exp(1j * N_2 * k_0 * d_2)
M_21_s = np.array([(N_2 + N_1) / (2 * N_2), (N_1 - N_2) / (2 * N_2),  # transfer matrix, s - polarization
                   (N_1 - N_2) / (2 * N_2), (N_2 + N_1) / (2 * N_2)])
M_21_p = np.array([((epsilon_1 * N_2 + epsilon_2 * N_1) / (2 * N_2 * epsilon_1 * epsilon_2)) * n_1 * n_2,
                   # transfer matrix, p - polarization
                   ((epsilon_1 * N_2 - epsilon_2 * N_1) / (2 * N_2 * epsilon_1 * epsilon_2)) * n_1 * n_2,
                   ((epsilon_1 * N_2 - epsilon_2 * N_1) / (2 * N_2 * epsilon_1 * epsilon_2)) * n_1 * n_2,
                   ((epsilon_1 * N_2 + epsilon_2 * N_1) / (2 * N_2 * epsilon_1 * epsilon_2)) * n_1 * n_2])
M_32_s = np.array([(N_3 + N_2) / (2 * N_3), (N_2 - N_3) / (2 * N_3),  # transfer matrix, s - polarization
                   (N_2 - N_3) / (2 * N_3), (N_3 + N_2) / (2 * N_3)])
M_32_p = np.array([((epsilon_2 * N_3 + epsilon_3 * N_2) / (2 * N_3 * epsilon_2 * epsilon_3)) * n_2 * n_3,
                   # transfer matrix, p - polarization
                   ((epsilon_2 * N_3 - epsilon_3 * N_2) / (2 * N_3 * epsilon_2 * epsilon_3)) * n_2 * n_3,
                   ((epsilon_2 * N_3 - epsilon_3 * N_2) / (2 * N_3 * epsilon_2 * epsilon_3)) * n_2 * n_3,
                   ((epsilon_2 * N_3 + epsilon_3 * N_2) / (2 * N_3 * epsilon_2 * epsilon_3)) * n_2 * n_3])
PHI_2 = [phi_2, 0, 0, phi_2 ** (-1)]  # propagation matrix
A_s = [M_32_s[0] * PHI_2[0] + M_32_s[1] * PHI_2[2], M_32_s[0] * PHI_2[2] + M_32_s[1] * PHI_2[3],
       M_32_s[2] * PHI_2[0] + M_32_s[3] * PHI_2[2], M_32_s[2] * PHI_2[2] + M_32_s[3] * PHI_2[3]]  # M32 * PHI2
Tw_s = [A_s[0] * M_21_s[0] + A_s[1] * M_21_s[2], A_s[0] * M_21_s[2] + A_s[1] * M_21_s[3],
        A_s[2] * M_21_s[0] + A_s[3] * M_21_s[2], A_s[2] * M_21_s[2] + A_s[3] * M_21_s[3]]
A_p = [M_32_p[0] * PHI_2[0] + M_32_p[1] * PHI_2[2], M_32_p[0] * PHI_2[2] + M_32_p[1] * PHI_2[3],
       M_32_p[2] * PHI_2[0] + M_32_p[3] * PHI_2[2], M_32_p[2] * PHI_2[2] + M_32_p[3] * PHI_2[3]]  # M32 * PHI2
Tw_p = [A_p[0] * M_21_p[0] + A_p[1] * M_21_p[2], A_p[0] * M_21_p[2] + A_p[1] * M_21_p[3],
        A_p[2] * M_21_p[0] + A_p[3] * M_21_p[2], A_p[2] * M_21_p[2] + A_p[3] * M_21_p[3]]
R_s = - (Tw_s[2] / Tw_s[3])
R_p = - (Tw_p[2] / Tw_p[3])
T_p = Tw_p[0] + R_p * Tw_p[1]
T_s = Tw_s[0] + R_s * Tw_s[1]
plt.plot(alpha, abs(R_p) ** 2, color='blue', label='p - polarization')
plt.plot(alpha, abs(R_s) ** 2, color='green', label='s - polarization')
# plt.plot(alpha, abs(T_p) ** 2, color='blue', label='p - polarization')
# plt.plot(alpha, abs(T_s) ** 2, color='green', label='s - polarization')
plt.legend(loc='best')
plt.xlabel('alpha')
plt.ylabel('|R|²')
plt.show()
