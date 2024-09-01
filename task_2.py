import math
import cmath

f = 170 * (10 ** 9)  # GH
epsilon_1 = 1
epsilon_2 = 12
epsilon_3 = 1
n_1 = math.sqrt(epsilon_1)
n_2 = math.sqrt(epsilon_2)
n_3 = math.sqrt(epsilon_3)
c = 3 * (10 ** 8)
d_2 = 0.3  # meters
k_0 = (2 * math.pi * f) / c
alpha = math.pi / 50  # radian (3.6 degrees)
k_x = k_0 * math.sin(alpha)
N_1 = math.sqrt(epsilon_1 - ((k_x ** 2) / (k_0 ** 2)))
N_2 = math.sqrt(epsilon_2 - ((k_x ** 2) / (k_0 ** 2)))
N_3 = math.sqrt(epsilon_3 - ((k_x ** 2) / (k_0 ** 2)))
phi_2 = cmath.exp(1j * N_2 * k_0 * d_2)
M_21_s = [(N_2 + N_1) / (2 * N_2), (N_2 - N_1) / (2 * N_2),  # transfer matrix, s - polarization
          (N_2 - N_1) / (2 * N_2), (N_2 + N_1) / (2 * N_2)]
M_21_p = [(epsilon_1 * N_2 + epsilon_2 * N_1) / (2 * N_2 * n_2 * n_1),  # transfer matrix, p - polarization
          (epsilon_1 * N_2 - epsilon_2 * N_1) / (2 * N_2 * n_2 * n_1),
          (epsilon_1 * N_2 - epsilon_2 * N_1) / (2 * N_2 * n_2 * n_1),
          (epsilon_1 * N_2 + epsilon_2 * N_1) / (2 * N_2 * n_2 * n_1)]
M_32_s = [(N_3 + N_2) / (2 * N_3), (N_3 - N_2) / (2 * N_3),  # transfer matrix, s - polarization
          (N_3 - N_2) / (2 * N_3), (N_3 + N_2) / (2 * N_3)]
M_32_p = [(epsilon_2 * N_3 + epsilon_3 * N_2) / (2 * N_3 * n_3 * n_2),  # transfer matrix, p - polarization
          (epsilon_2 * N_3 - epsilon_3 * N_2) / (2 * N_3 * n_3 * n_2),
          (epsilon_2 * N_3 - epsilon_3 * N_2) / (2 * N_3 * n_3 * n_2),
          (epsilon_2 * N_3 + epsilon_3 * N_2) / (2 * N_3 * n_3 * n_2)]
PHI_2 = [phi_2, 0, 0, phi_2 ** (-1)]  # propagation matrix
A_s = [M_32_s[0] * PHI_2[0] + M_32_s[1] * PHI_2[2], M_32_s[0] * PHI_2[2] + M_32_s[1] * PHI_2[3],
       M_32_s[2] * PHI_2[0] + M_32_s[3] * PHI_2[2], M_32_s[2] * PHI_2[2] + M_32_s[3] * PHI_2[3]]  # M32 * PHI2
T_s = [A_s[0] * M_21_s[0] + A_s[1] * M_21_s[2], A_s[0] * M_21_s[2] + A_s[1] * M_21_s[3],
       A_s[2] * M_21_s[0] + A_s[3] * M_21_s[2], A_s[2] * M_21_s[2] + A_s[3] * M_21_s[3]]
A_p = [M_32_p[0] * PHI_2[0] + M_32_p[1] * PHI_2[2], M_32_p[0] * PHI_2[2] + M_32_p[1] * PHI_2[3],
       M_32_p[2] * PHI_2[0] + M_32_p[3] * PHI_2[2], M_32_p[2] * PHI_2[2] + M_32_p[3] * PHI_2[3]]  # M32 * PHI2
T_p = [A_p[0] * M_21_p[0] + A_p[1] * M_21_p[2], A_p[0] * M_21_p[2] + A_p[1] * M_21_p[3],
       A_p[2] * M_21_p[0] + A_p[3] * M_21_p[2], A_p[2] * M_21_p[2] + A_p[3] * M_21_p[3]]
R_s = - (T_s[2] / T_s[3])
R_p = - (T_p[2] / T_p[3])
T_p = cmath.sqrt(1 - (R_p ** 2))
T_s = cmath.sqrt(1 - (R_s ** 2))
print(R_p, R_s)
print(T_p, T_s)
print(R_p.real, T_p.real)
print(R_s.real, T_s.real)
#(0.820406304605871+0.1438778225369187j) (-0.8215540897979996-0.14363824110227494j)
#(0.6196170777085491-0.19050196798767824j) (0.6181612905403462-0.1908993433183326j)
#0.820406304605871 0.6196170777085491
#-0.8215540897979996 0.6181612905403462
