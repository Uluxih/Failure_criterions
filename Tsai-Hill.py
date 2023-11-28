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
    plt.plot(u + Ell_rot[0, :], v + Ell_rot[1, :], color, label=name_curve, linestyle='dashed')  # rotated ellipse

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

def plot_tsai_hill(X,Y,color,text):
    x, y = sp.symbols('x y')
    # Тсаи-Хилл
    print("Tsai-Hill c-c:")
    A = 1/X**2
    B = -0.5/max(abs(X), abs(Y))**2
    C = 1/Y**2
    F = -1
    D, E, = 0, 0
    print("A=", A, "B=", B,"C=",C)
    res = rotate_ellise(A, B, C, D, E, F)
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

    plot_elips(solution[x], solution[y], a, b, alpha, color,text )

fig, ax = plt.subplots(figsize=(7, 6))

plt.grid(color='lightgray', linestyle='--')
ax.set_title("Tsai-Hill failure criterion")
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

plot_tsai_hill(X,Y,"blue", "Tsai-Hill compression-compression")
plot_tsai_hill(X/2,Y/2,"red", "Tsai-Hill tension-tension")
plot_tsai_hill(X/2,Y,"violet", "Tsai-Hill tension-compression")
plot_tsai_hill(X,Y/2,"orange", "Tsai-Hill compression-tension")
ax.set_aspect('equal', adjustable='box')
fig.legend(loc="lower center")
plt.show()
