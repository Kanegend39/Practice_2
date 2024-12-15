import numpy as np


def matrix_method(f, epsilon_1, epsilon_2, epsilon_3, d_2, alpha):
    n_1 = np.sqrt(epsilon_1)
    n_2 = np.sqrt(epsilon_2)
    n_3 = np.sqrt(epsilon_3)
    c = 3 * 1e10
    k_0 = (2 * np.pi * f) / c
    k_x = k_0 * alpha
    N_1 = np.sqrt(epsilon_1 - ((k_x ** 2) / (k_0 ** 2)), dtype='complex64')
    N_2 = np.sqrt(epsilon_2 - ((k_x ** 2) / (k_0 ** 2)), dtype='complex64')
    N_3 = np.sqrt(epsilon_3 - ((k_x ** 2) / (k_0 ** 2)), dtype='complex64')
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
    return T_p, R_p, T_s, R_s
