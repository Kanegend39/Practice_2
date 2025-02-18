import numpy as np
from scipy.fft import fft2, fftfreq, ifft2, ifft
from matrix_method import matrix_method


def gaussian_beam_propagation_vector_angle(left, right, step, w, phi, z_shtrih, d_i, epsilons, n):
    # CONSTANTS #
    right += step
    c = 3e10  # cm/s
    f = 170e9  # Gh
    k0 = (2 * np.pi * f) / c  # cm^-1
    x = np.arange(left, right, step)
    y = np.arange(left, right, step)
    X_shtrih, Y_shtrih = np.meshgrid(x, y)
    v0 = -1j * k0 * (w ** 2)
    E_x_shtrih = (1 / v0) * (-1j * k0 * (w ** 2)) * np.exp(1j * k0 * ((X_shtrih ** 2 + Y_shtrih ** 2) / 2) * (1 / v0))
    E_y_shtrih = 0 * X_shtrih
    E_z_shtrih = (1 / (v0 ** 2)) * (1j * X_shtrih * k0 * (w ** 2)) * np.exp(1j * k0 * ((X_shtrih ** 2 + Y_shtrih ** 2) / 2) * (1 / v0))
    E_x = np.empty(E_x_shtrih.shape, dtype="complex64")
    E_y = E_y_shtrih
    E_z = np.empty(E_z_shtrih.shape, dtype="complex64")
    X = X_shtrih
    Y = Y_shtrih
    rotation_matrix = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
    for i in range(len(E_x_shtrih)):
        for j in range(len(E_x_shtrih)):
            A = rotation_matrix @ np.array([[E_x_shtrih[i][j]], [E_z_shtrih[i][j]]])
            E_x[i][j] = A[0][0]
            E_z[i][j] = A[1][0]

    # DECOMPOSITION OF A GAUSSIAN BEAM INTO A FOURIER SPECTRUM #
    fx = fftfreq(len(x), d=x[1] - x[0])
    fy = fftfreq(len(y), d=y[1] - y[0])
    FX, FY = np.meshgrid(fx, fy)
    k_x = 2 * np.pi * FX
    k_y = 2 * np.pi * FY
    k_z = np.where(k0 ** 2 >= k_x ** 2 + k_y ** 2, np.sqrt(k0 ** 2 - k_x ** 2 - k_y ** 2, dtype='complex64'),
                   1j * np.sqrt(k_x ** 2 + k_y ** 2 - k0 ** 2, dtype='complex64'))
    unit_vector_E_s_x = np.empty(k_x.shape, dtype="complex64")
    unit_vector_E_s_y = np.empty(k_x.shape, dtype="complex64")
    unit_vector_E_s_z = np.empty(k_x.shape, dtype="complex64")
    for i in range(len(k_x)):
            if k_x[i][j] ** 2 + k_y[i][j] ** 2 == 0:
                unit_vector_E_s_x[i][j] = np.cos(phi)
                unit_vector_E_s_y[i][j] = 0
                unit_vector_E_s_z[i][j] = -np.sin(phi)
            else:
                unit_vector_E_s_x[i][j] = (k_y[i][j] / np.sqrt(abs(k_x[i][j]) ** 2 + abs(k_y[i][j]) ** 2))
                unit_vector_E_s_y[i][j] = - (k_x[i][j] / np.sqrt(abs(k_x[i][j]) ** 2 + abs(k_y[i][j]) ** 2))
                unit_vector_E_s_z[i][j] = 0
    unit_vector_E_p_x = np.empty(k_x.shape, dtype="complex64")
    unit_vector_E_p_y = np.empty(k_x.shape, dtype="complex64")
    unit_vector_E_p_z = np.empty(k_x.shape, dtype="complex64")
    for i in range(len(k_x)):
        for j in range(len(k_x)):
            if k0 ** 2 >= k_x[i][j] ** 2 + k_y[i][j] ** 2:
                unit_vector_E_p_x[i][j] = (unit_vector_E_s_y[i][j] * k_z[i][j]) / k0
                unit_vector_E_p_y[i][j] = - (unit_vector_E_s_x[i][j] * k_z[i][j]) / k0
                unit_vector_E_p_z[i][j] = (unit_vector_E_s_x[i][j] * k_y[i][j] - unit_vector_E_s_y[i][j] * k_x[i][
                    j]) / k0
            else:
                unit_vector_E_p_x[i][j] = np.sin(phi)
                unit_vector_E_p_y[i][j] = 0
                unit_vector_E_p_z[i][j] = np.cos(phi)
    C_s = fft2(E_x) * unit_vector_E_s_x + fft2(E_y) * unit_vector_E_s_y + fft2(E_z) * unit_vector_E_s_z
    C_p = fft2(E_x) * unit_vector_E_p_x + fft2(E_y) * unit_vector_E_p_y + fft2(E_z) * unit_vector_E_p_z

    T_p = np.empty(k_x.shape, dtype='complex128')
    R_p = np.empty(k_x.shape, dtype='complex128')
    T_s = np.empty(k_x.shape, dtype='complex128')
    R_s = np.empty(k_x.shape, dtype='complex128')
    for i in range(len(k_x)):
        for j in range(len(k_y)):
            alpha = k_x[i][j] / k0
            result = matrix_method(f, epsilons, d_i, alpha, n)
            T_p[i][j], R_p[i][j], T_s[i][j], R_s[i][j] = result

    # GAUSSIAN BEAM AFTER THE INVERSE FOURIER TRANSFORM CONSIDERING ITS PASSAGE THROUGH THE PLATE #
    z = z_shtrih * np.cos(phi)
    E_x = ifft2(C_p * unit_vector_E_p_x * np.exp(1j * k_z * z) + C_s * unit_vector_E_s_x * np.exp(1j * k_z * z))
    E_y = ifft2(C_p * unit_vector_E_p_y * np.exp(1j * k_z * z) + C_s * unit_vector_E_s_y * np.exp(1j * k_z * z))
    E_z = ifft2(C_p * unit_vector_E_p_z * np.exp(1j * k_z * z) + C_s * unit_vector_E_s_z * np.exp(1j * k_z * z))
    C_s_T = T_s * C_s * np.exp(1j * k_z * z)
    C_p_T = T_p * C_p * np.exp(1j * k_z * z)
    E_x_T_ = ifft2(C_p_T * unit_vector_E_p_x + C_s_T * unit_vector_E_s_x)
    E_y_T_ = ifft2(C_p_T * unit_vector_E_p_y + C_s_T * unit_vector_E_s_y)
    E_z_T_ = ifft2(C_p_T * unit_vector_E_p_z + C_s_T * unit_vector_E_s_z)
    for i in range(len(E_x_shtrih)):
        for j in range(len(E_x_shtrih)):
            A = np.linalg.inv(rotation_matrix) @ np.array([[E_x[i][j]], [E_z[i][j]]])
            E_x_shtrih[i][j] = A[0][0]
            E_z_shtrih[i][j] = A[1][0]
    E_y_shtrih = E_y
    E_x_T_shtrih = np.empty(E_x_shtrih.shape, dtype="complex64")
    E_y_T_shtrih = E_y_T_
    E_z_T_shtrih = np.empty(E_z_shtrih.shape, dtype="complex64")
    for i in range(len(E_x_shtrih)):
        for j in range(len(E_x_shtrih)):
            A = np.linalg.inv(rotation_matrix) @ np.array([[E_x_T_[i][j]], [E_z_T_[i][j]]])
            E_x_T_shtrih[i][j] = A[0][0]
            E_z_T_shtrih[i][j] = A[1][0]
            X[i][j] = X_shtrih[i][j] * np.cos(phi) + z_shtrih * np.sin(phi)
    return X, Y, E_x_shtrih, E_y_shtrih, E_z_shtrih, E_x_T_, E_y_T_, E_z_T_
