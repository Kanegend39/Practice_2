import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import rfft, rfftfreq
from scipy.signal import hilbert


def spectrum(A):    #функция рисует график спектра
    S = rfft(A)
    f = rfftfreq(5000, 1 / 1000) # 1 число - количество точек; 2 число - 1 / (количество точек в 1 секунде)
    plt.plot(f, np.abs(S), color='red')
    plt.xlabel("f, ГГц")
    plt.ylabel('A.U.')


t = np.arange(5, 10, 0.001)  # нс
f_0 = 1  # ГГц
w_0 = 2 * math.pi * f_0
A_0 = 1
A = A_0 * np.cos(w_0 * t)
plt.xlabel("t, нс")
plt.ylabel('A.U.')
plt.grid(color='black', linewidth=0.5)
plt.grid(True)
plt.plot(t, A, color='orange')
F = hilbert(A)
F[0], F[-1] = 0, 0
plt.plot(t, np.sqrt(F.real ** 2 + F.imag ** 2), color='green')
#spectrum(A)
plt.show()
