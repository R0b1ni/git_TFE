import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import sys
from scipy.optimize import minimize
import matplotlib.cm as cm


# Data
P = 1000 # W
c_x = 1.5 # m/s (axial speed)
rho_w = 1000 # kg/m3
R_i = 0.02 # m (inner radius)
R_o = 0.05 # m (outer radius)
g = 9.81 # m/s2
H = 1 # m (height)
N_q = 250

Q = c_x * np.pi * (R_o**2 - R_i**2)
n = (H)**(5/4) * (N_q) / (P**(1/2))
Re = (rho_w * c_x * (R_o)) / (1e-3) # dynamic viscosity of water ~ 1e-3 Pa.s

#velocity triangles
U = 2 * np.pi * n * (R_i + R_o) / 2

c_1x = c_x
c_1 = c_1x # assuming axial inflow
c_2x = c_x # mass conservation

w1 = np.sqrt((c_1x - U)**2)
c_2teta = P/(rho_w * Q * U)
c_2 = np.sqrt(c_2x**2 + c_2teta**2)
w2 = np.sqrt((c_2teta + U)**2 + c_2x**2)
alpha_1 = 0.0
alpha_2 = np.degrees(np.arctan2(c_2teta, c_2x))
beta_1 = np.degrees(np.arctan2(U, c_1x))
beta_2 = np.degrees(np.arctan2(c_2teta + U, c_2x))
rotor_deviation = beta_1 - beta_2

print(f"Specific speed n: {n:.2f} t/s")
print(f"mean blade speed U: {U:.2f} m/s")
print(f"Inlet absolute velocity c1: {c_1:.2f} m/s")
print(f"Outlet absolute velocity c2: {c_2:.2f} m/s")
print(f"Inlet relative velocity w1: {w1:.2f} m/s")
print(f"Outlet relative velocity w2: {w2:.2f} m/s")
print(f"Inlet flow angle alpha1: {alpha_1:.2f} degrees")
print(f"Outlet flow angle alpha2: {alpha_2:.2f} degrees")
print(f"Inlet relative flow angle beta1: {beta_1:.2f} degrees")
print(f"Outlet relative flow angle beta2: {beta_2:.2f} degrees")
print(f"Rotor deviation: {rotor_deviation:.2f} degrees")

def c_teta(r):
    A = ((R_i + R_o) / 2) * c_2teta
    return A / r

def rotor_deviation_at_r(r):
    c_teta_r = c_teta(r)
    beta_2_r = np.degrees(np.arctan2(c_teta_r + U, c_2x))
    return beta_1 - beta_2_r

print("\nRotor deviation at different radii:")
for r in [R_i, (R_i + R_o) / 2, R_o]:
    dev = rotor_deviation_at_r(r)
    print(f"Radius {r:.3f} m: Deviation {dev:.2f} degrees")

twist_angle = rotor_deviation_at_r(R_i) - rotor_deviation_at_r(R_o)
print(f"\nTwist angle from hub to tip: {twist_angle:.2f} degrees")

def U(r):
    return 2 * np.pi * n * r