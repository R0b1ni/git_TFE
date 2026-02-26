import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import sys
from scipy.optimize import minimize, fsolve
import matplotlib.cm as cm
from scipy.interpolate import make_interp_spline
from scipy.interpolate import PchipInterpolator
from matplotlib.patches import FancyArrowPatch

mu_0 = 4 * np.pi * 1e-7
mu_r = 1.05
p = 12 # poles
L = 24 * 1e-3
lg = 10 * 1e-3
R = (50 - (lg * 1e3)/2) * 1e-3
N  = 120
f = (p * 2000) / 120 # 6000 rpm
E = 10.26
k_w = 0.85
Br = 1.1
C_teta = 0.685

phi_g = E/(4.44 * f * N * k_w)
Rg = (lg * p) / (2 * np.pi * mu_0 * R * L)
hm = C_teta * (2 * np.pi * R) / p
Rm = lambda lm : lm  / (mu_0 * mu_r * hm * L)

lm = fsolve(lambda x: x - (phi_g * (Rm(x) + Rg) * mu_0)/(Br), 0.1)[0]
k_f = 1 + lg/hm

print("PM length: ", lm * 1e3, "mm")
print("PM height: ", hm * 1e3, "mm")
print("k_f: ", k_f)
print("PM reluctance: ", Rm(lm), "A/Wb")
print("Air gap reluctance: ", Rg, "A/Wb")
print("Air gap flux: ", phi_g, "Wb")