from math import *
import math
import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
# from scipy.interpolate import griddata
import sympy as sp


# from mpl_toolkits.mplot3d import Axes3D


def plot_elips(x0, y0, a, b, alpha_rad, color, name_curve):
    u = x0  # x-position of the center
    v = y0  # y-position of the center
    t_rot = alpha_rad  # rotation angle

    t = np.linspace(0, 2 * pi, 100)
    Ell = np.array([a * np.cos(t), b * np.sin(t)])
    # u,v removed to keep the same center location
    R_rot = np.array([[cos(t_rot), -sin(t_rot)], [sin(t_rot), cos(t_rot)]])
    # 2-D rotation matrix

    Ell_rot = np.zeros((2, Ell.shape[1]))
    for i in range(Ell.shape[1]):
        Ell_rot[:, i] = np.dot(R_rot, Ell[:, i])
    plt.plot(u + Ell_rot[0, :], v + Ell_rot[1, :], color, label=name_curve)  # rotated ellipse

def rotate_ellise(A,B,C,D,E,F):
    b_array = np.array([[A, B], [B, C]])
    d_array = np.array([[A, B, D], [B, C, E], [D, E, F]])
    det_b = np.linalg.det(b_array)
    det_d = np.linalg.det(d_array)
    S = A + C
    print("S",S, "det_b", det_b,"det_d", det_d)
    A1, C1, F1 = sp.symbols('A1 C1 F1')
    eq1 = sp.Eq(A1 + C1, S)
    eq2 = sp.Eq(A1 * C1, det_b)
    eq3 = sp.Eq(A1 * C1 * F1, det_d)
    solution = sp.solve((eq1, eq2, eq3), (A1, C1, F1))
    A1 = solution[0][0]
    C1 = solution[0][1]
    F1 = solution[0][2]
    print("A1=", A1, " C1=", C1, " F1=", F1)
    alpha = math.atan((A1 - A) / B)
    print("alpha=", math.degrees(alpha))
    res=[]
    res.append(solution[0])
    res.append(alpha)
    return res

fig, ax = plt.subplots(figsize=(7, 6))

plt.grid(color='lightgray', linestyle='--')
ax.set_title("Критерии прочности")
ax.set_xlabel("SIGMA1")
ax.set_ylabel("SIGMA2")

# Сопротивления растяжению в осях ортотропии
X = 1.2
Y = 0.6
Z = 1
# Сопротивления срезу в осях ортотропии
R = 0.6
S = 0.6
T = 0.7

F = (-1 / X ** 2 + 1 / Y ** 2 + 1 / Z ** 2) / 2
G = (1 / X ** 2 - 1 / Y ** 2 + 1 / Z ** 2) / 2
H = (1 / X ** 2 + 1 / Y ** 2 - 1 / Z ** 2) / 2
L = 1 / (2 * R ** 2)
M = 1 / (2 * S ** 2)
N = 1 / (2 * T ** 2)
print(F)
print(G)
print(H)
# print(L)
# print(M)
# print(N)

# sigma1=Symbol("SIGMA1")
# sigma2=Symbol("SIGMA2")
# sigma3=Symbol("SIGMA3")
sigma1 = 0
sigma2 = 1
sigma3 = 0
print("sigma1= ", sigma1, " sigma2= ", sigma2)
expr_Mises_Hill = F * (sigma2 - sigma3) ** 2 + G * (sigma3 - sigma1) ** 2 + H * (sigma1 - sigma2) ** 2
expr_Mises = (sigma2 - sigma3) ** 2 + (sigma3 - sigma1) ** 2 + (sigma1 - sigma2) ** 2
print("Mises-Hill intens=", expr_Mises_Hill)
print("Mises intens=", expr_Mises)
print("Mises-Hill:")
A = H + G
B = -H
C = F + H
F = -1
D, E, = 0, 0
print("A=", A, "B= ", B, "C=", C)
b_array = np.array([[A, B], [B, C]])
d_array = np.array([[A, B, D], [B, C, E], [D, E, F]])

# print(b_array)
# print(d_array)

det_b = np.linalg.det(b_array)
det_d = np.linalg.det(d_array)
S = A + C
print("S",S, "det_b", det_b,"det_d", det_d)

A1, C1, F1 = sp.symbols('A1 C1 F1')
eq1 = sp.Eq(A1 + C1, S)
eq2 = sp.Eq(A1 * C1, det_b)
eq3 = sp.Eq(A1 * C1 * F1, det_d)
solution = sp.solve((eq1, eq2, eq3), (A1, C1, F1))
A1 = solution[0][0]
C1 = solution[0][1]
F1 = solution[0][2]
print(solution)
print("A1=", A1, " C1=", C1, " F1=", F1)
alpha = math.atan((A1 - A) / B)
print("alpha=", math.degrees(alpha))

# полуоси эллипса:
a = math.sqrt(-F1 / A1)
b = math.sqrt(-F1 / C1)
print("a=", a, " b=", b)
# Центр эллипса

x, y = sp.symbols('x y')
eq1 = sp.Eq(A * x + B * y + D, 0)
eq2 = sp.Eq(B * x + C * y + E, 0)
solution = sp.solve((eq1, eq2), (x, y))
print("x0= ", solution[x], "y0= ", solution[y])

plot_elips(solution[x], solution[y], a, b, alpha, "darkorange", "Mises-Hill")
#plot_elips(0, 0, a, b, alpha, "darkorange", "Mises-Hill")
# Чистый Мизис
print("Mises:")
A = 2
B = -1
C = 2
F = -1
D, E, = 0, 0

b_array = np.array([[A, B], [B, C]])
d_array = np.array([[A, B, D], [B, C, E], [D, E, F]])
det_b = np.linalg.det(b_array)
det_d = np.linalg.det(d_array)
# print(det_b, det_d)
S = A + C

A1, C1, F1 = sp.symbols('A1 C1  F1')
eq1 = sp.Eq(A1 + C1, S)
eq2 = sp.Eq(A1 * C1, det_b)
eq3 = sp.Eq(A1 * C1 * F1, det_d)
solution = sp.solve((eq1, eq2, eq3), (A1, C1, F1))
A1 = solution[0][0]
C1 = solution[0][1]
F1 = solution[0][2]
print("A1=", A1, " C1=", C1, " F1=", F1)
alpha = math.atan((A1 - A) / B)
print("alpha=", math.degrees(alpha))

# полуоси эллипса:
a = math.sqrt(-F1 / A1)/0.707
b = math.sqrt(-F1 / C1)/0.707
print("a=", a, " b=", b)
# Центры эллипса
eq1 = sp.Eq(A * x + B * y + D, 0)
eq2 = sp.Eq(B * x + C * y + E, 0)
solution = sp.solve((eq1, eq2), (x, y))
print("x0= ", solution[x], "y0= ", solution[y])

plot_elips(solution[x], solution[y], a, b, alpha, "darkblue", "Mises")
#plot_elips(0, 0, a, b, alpha, "darkblue", "Mises")

# Тсаи-Хилл
print("Tsai-Hill:")
A = 1/X**2
B = -0.5/max(abs(X), abs(Y))**2
C = 1/Y**2
F = -1
D, E, = 0, 0
print("A=", A, "B=", B,"C=",C)
res = rotate_ellise(A,B,C,D,E,F)
A1 = res[0][0]
C1 = res[0][1]
F1 = res[0][2]
alpha=res[1]
# полуоси эллипса:
a = math.sqrt(-F1 / A1)
b = math.sqrt(-F1 / C1)
print("a=", a, " b=", b)
# Центры эллипса
eq1 = sp.Eq(A * x + B * y + D, 0)
eq2 = sp.Eq(B * x + C * y + E, 0)
solution = sp.solve((eq1, eq2), (x, y))
print("x0= ", solution[x], "y0= ", solution[y])

plot_elips(solution[x], solution[y], a, b, alpha, "red", "Tsai-Hill")
ax.set_aspect('equal', adjustable='box')
fig.legend(loc="lower center")
plt.show()
