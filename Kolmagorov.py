from math import *
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
import numpy as np
# import pandas as pd
# from scipy.interpolate import griddata
import sympy as sp
from functools import reduce
import operator
from scipy.interpolate import make_interp_spline, BSpline

sigma_f=750
sigma_f02=500
tau_b=1000
#(sigma1-sigma2)/2

sigma1 = 1.001
sigma2 = 1
sigma3 = 0

if(sigma_f/sigma_f02<=1.32):
    m = (log(sigma_f/sigma_f02,e)+0.056)/3.44
else:
    m = (log(sigma_f/sigma_f02,e)+0.216)/4.78
a = (tau_b*sqrt(3)/sigma_f)**(1/m)
SIGMAs=[]
SIGMAsX=[]
SIGMAsY=[]
# a=1.74
# m=0.136
sigma_i = 1 / sqrt(2) * sqrt((sigma1 - sigma2) ** 2 + (sigma2 - sigma3) ** 2 + (sigma3 - sigma1) ** 2)
sigma0 = (sigma1 + sigma2 + sigma3) / 3
sigma_eq = sigma_i / ((a * exp(-3 * log(a) * sigma0 / sigma_i)) ** m)
print(sigma_eq)
while (sigma1 >= -1):
    sigma2=1.001
    while (sigma2 >= -1):
        sigma_i = 1 / sqrt(2) * sqrt((sigma1 - sigma2) ** 2 + (sigma2 - sigma3) ** 2 + (sigma3 - sigma1) ** 2)
        sigma0 = (sigma1 + sigma2 + sigma3) / 3
        sigma_eq = sigma_i / ((a * exp(-3 * log(a) * sigma0 / sigma_i)) ** m)
        print("sigma1 =", round(sigma1,2), "\t sigma2 =", round(sigma2,2), "\t s0 =", round(sigma0,2),
              "\t si =", round(sigma_i,2), "\t s eq = ", round(sigma_eq,2), "\t sx =", round(sigma1/sigma_eq, 2),
              round(sigma2/sigma_eq, 2))
        SIGMAsX.append(sigma1/sigma_eq)
        SIGMAsY.append(sigma2/sigma_eq)
        sigs=[sigma1/sigma_eq, sigma2/sigma_eq]
        SIGMAs.append(sigs)
        sigma2 = sigma2 - 0.1
    sigma1 = sigma1 - 0.1

#сортирую точки по часавой стрелке
points=[]
coords = SIGMAs
center = tuple(map(operator.truediv, reduce(lambda x, y: map(operator.add, x, y), coords), [len(coords)] * 2))
points = (sorted(coords, key=lambda coord: (-135 - math.degrees(math.atan2(*tuple(map(operator.sub, coord, center))[::-1]))) % 360))
print(points)
nSigX=[]
nSigY=[]
for i in points:
    nSigX.append(i[0])
    nSigY.append(i[1])
nSigX.append(points[0][0])
nSigY.append(points[0][1])

#plt.scatter(SIGMAsX,SIGMAsY, color="black", linewidths=0.5, marker="v")

fig, ax = plt.subplots(figsize=(7, 6))
plt.grid(color='lightgray', linestyle='--')
ax.set_title(f"Критерий Колмогорова τ={tau_b}")
ax.set_xlabel("σ_x") #/σ_Fx
ax.set_ylabel("σ_y") #/σ_Fy
plt.plot(nSigX,nSigY, color="black")
ax.set_aspect('equal', adjustable='box')

# ax.set_xlim(-2, 2)
# ax.set_ylim(-2, 2)
ax.xaxis.set_major_locator(MultipleLocator(0.5))
# ax.yaxis.set_major_locator(MultipleLocator(10))
plt.show()
# print(a*sigma0/sigma_i)
# print(-3*log(a*sigma0/sigma_i, e))
