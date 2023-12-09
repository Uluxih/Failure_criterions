from math import *
import math
import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
# from scipy.interpolate import griddata
import sympy as sp
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)

sigma1=2
sigma2=3
sigma3=0
tau_nu_ravn= 1 / 3 * sqrt((sigma1 - sigma2) ** 2 + (sigma2 - sigma3) ** 2 + (sigma3 - sigma1) ** 2)
#направляющий вектор:
i=1
j=1
k=1
#длина вектора:
v_l=sqrt(i**2+j**2+k**2)
print("длина вектора:", v_l)
#направляющие косинусы:
l=i/v_l
m=j/v_l
n=k/v_l
print("l,m,n", l,m,n)
s_nu2= (sigma1 * l) ** 2 + (sigma2 * m) ** 2 + (sigma3 * n) ** 2
sigma_nu=(sigma1*l**2+sigma2*m**2+sigma3*n**2)
tau_nu=sqrt(s_nu2-sigma_nu**2)
print("tau_nu_ravn", tau_nu_ravn)
print("tau_nu", tau_nu)
