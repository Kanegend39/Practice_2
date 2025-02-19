import numpy as np


def matrix_method(f, epsilons, d_i, n, k_x, k_y):
    n += 2
    n_i = np.arange(0, n, 1., dtype='complex64')
    N_i = np.arange(0, n, 1., dtype='complex64')
    phi_i = np.arange(0, n - 2, 1., dtype='complex128')
    M_i_s = []
    M_i_p = []
    PHI_i = []
    c = 3 * 1e10
    k_0 = (2 * np.pi * f) / c
    for i in range(n):
        n_i[i] = np.sqrt(epsilons[i], dtype='complex64')
        N_i[i] = np.sqrt(epsilons[i] - ((k_x ** 2 + k_y ** 2) / (k_0 ** 2)), dtype='complex64')
    for i in range(n - 2):
        phi_i[i] = np.exp(1j * N_i[i + 1] * k_0 * d_i[i], dtype='complex128')
        PHI_i.append(np.array([[phi_i[i], 0], [0, phi_i[i] ** (-1)]], dtype='complex128'))
    for i in range(n - 1):
        M_i_s.append(np.array([[(N_i[i + 1] + N_i[i]) / (2 * N_i[i + 1]), (N_i[i] - N_i[i + 1]) / (2 * N_i[i + 1])],
                            [(N_i[i] - N_i[i + 1]) / (2 * N_i[i + 1]), (N_i[i + 1] + N_i[i]) / (2 * N_i[i + 1])]], dtype='complex128'))
        M_i_p.append(np.array([[((epsilons[i] * N_i[i + 1] + epsilons[i + 1] * N_i[i]) / (2 * N_i[i + 1] * epsilons[i] * epsilons[i + 1])) * n_i[i] * n_i[i + 1],
                       ((epsilons[i] * N_i[i + 1] - epsilons[i + 1] * N_i[i]) / (2 * N_i[i + 1] * epsilons[i] * epsilons[i + 1])) * n_i[i] * n_i[i + 1]],
                       [((epsilons[i] * N_i[i + 1] - epsilons[i + 1] * N_i[i]) / (2 * N_i[i + 1] * epsilons[i] * epsilons[i + 1])) * n_i[i] * n_i[i + 1],
                       ((epsilons[i] * N_i[i + 1] + epsilons[i + 1] * N_i[i]) / (2 * N_i[i + 1] * epsilons[i] * epsilons[i + 1])) * n_i[i] * n_i[i + 1]]], dtype='complex128'))
    W_s = M_i_s[n - 2]
    for i in range(n - 3, -1, -1):
        W_s = W_s @ PHI_i[i] @ M_i_s[i]
    W_p = M_i_p[n - 2]
    for i in range(n - 3, -1, -1):
        W_p = W_p @ PHI_i[i] @ M_i_p[i]
    R_s = - (W_s[1][0] / W_s[1][1])
    R_p = - (W_p[1][0] / W_p[1][1])
    T_p = W_p[0][0] + R_p * W_p[0][1]
    T_s = W_s[0][0] + R_s * W_s[0][1]
    return T_p, R_p, T_s, R_s
