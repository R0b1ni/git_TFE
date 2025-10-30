import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import sys
from scipy.optimize import minimize
import matplotlib.cm as cm
from scipy.interpolate import make_interp_spline
from scipy.interpolate import PchipInterpolator

# Data
P = 1000 # W
c_x = 2 # m/s (axial speed)
rho_w = 1000 # kg/m3
R_o = 0.155 # m (outer radius)
g = 9.81 # m/s2
H = 1 # m (height)
R_i = 0.4 * R_o # m (inner radius) (from fig 15 of 2022 Abeykoon)
sigma = 2.0 # blockage ratio (from fig 15 of 2022 Abeykoon)

Q = c_x * np.pi * (R_o**2 - R_i**2)
N = sigma * (2 * g * H)**(3/4) / (2*np.sqrt(np.pi * Q))
n = N * 60 # rpm


Re = (rho_w * c_x * (R_o)) / (1e-3) # dynamic viscosity of water ~ 1e-3 Pa.s

print(f"Rotationnal speed N: {n :.2f} rpm <=> {N :.2f} t/s \n")

#velocity triangles
U_m = 2 * np.pi * N * (R_i + R_o) / 2 # mean blade speed

Radii = np.linspace(R_i, R_o, 5)

c_1x = c_x
c_1 = c_1x # assuming axial inflow
c_2x = c_x # mass conservation

# 1D design
c_2teta = P/(rho_w * Q * U_m)

def c_teta(r):
    A = ((R_i + R_o) / 2) * c_2teta
    return A / r # free vortex distribution = constant work along the blade ?

w1 = []
w2 = []
w_inf = []
c2 = []
alpha_2 = []
beta_1 = []
beta_2 = []
beta_inf = []
chord_to_pitch = [1.30, 1.163, 1.025, 0.888, 0.750] # à vérifier, valeurs du livre
pitch = []
chord = []

for i in range (len(Radii)):
    r = Radii[i]
    U = 2 * np.pi * N * r
    pitch.append( 2 * np.pi * r / 4) # assuming 4 blades
    chord.append( chord_to_pitch[i] * pitch[i] )
    w1.append( np.sqrt(c_1x**2 + U**2) )
    c2.append( np.sqrt(c_2x**2 + c_teta(r)**2) )
    w2.append( np.sqrt((c_teta(r) + U)**2 + c_2x**2) )
    w_inf.append( (w1[i] + w2[i]) / 2 )
    alpha_2.append( np.degrees(np.arctan2(c_2teta, c_2x)) )
    beta_1.append( np.degrees(np.arctan2(U, c_1x)) +90)
    beta_2.append( np.degrees(np.arctan2(c_teta(r) + U, c_2x)) +90 )
    beta_inf.append(( beta_1[i] + beta_2[i] ) /2 )
    
eta_r = P / (rho_w * g * Q * H)
print(f"Hydraulic efficiency eta: {eta_r :.2f} \n")
print("Hydraulic Power P_hydraulic: {:.2f} W \n".format(rho_w * g * Q * H))
df = pd.DataFrame({
    "radius_m": Radii,
    "U_m_s": [2 * np.pi * N * r for r in Radii],
    "w1_m_s": w1,
    "w2_m_s": w2,
    "beta1_deg": beta_1,
    "beta2_deg": beta_2,
    "c2_m_s": c2,
    "alpha2_deg": alpha_2,
    "pitch_m": pitch,
    "chord_m": chord,
    "chord_to_pitch": chord_to_pitch
})
pd.options.display.float_format = "{:,.3f}".format
print(df.to_string(index=False))
print("\n")

#calculation of the x, y, z coordinates of the blades
x_1 = []
y_1 = []
z_1 = []
x_2 = []
y_2 = []
z_2 = []
x_inf = []
y_inf = []
z_inf = []

for i in range (len(Radii)):
    D = 2 * Radii[i]
    
    y_1.append( np.sin( np.radians(beta_1[i] - 90)) * ((D * np.sin(np.radians((chord[i]/D )* (180/np.pi))))/(2 * np.cos(np.radians(beta_inf[i] - beta_1[i])))) )
    z_1.append(y_1[i] / np.tan( np.radians(beta_1[i] - 90) ) )
    x_1.append((D/2) - D*np.sin( np.arcsin(2*y_1[i]/D)/2 )**2 )
    
    y_2.append(-np.sin( np.radians(beta_2[i] - 90)) * ((D * np.sin(np.radians((chord[i]/D )* (180/np.pi))))/(2 * np.cos(np.radians(beta_2[i] - beta_inf[i])))) )
    z_2.append(y_2[i] / np.tan( np.radians(beta_2[i] - 90) ) )
    x_2.append((D/2) - D*np.sin( np.arcsin(2*y_2[i]/D)/2 )**2 )
    
    x_inf.append( D/2)
    y_inf.append( 0)
    z_inf.append( 0 )
    
df = pd.DataFrame({
    "Diameter": Radii * 2,
    "x1_m": x_1,
    "y1_m": y_1,
    "z1_m": z_1,
    "x2_m": x_2,
    "y2_m": y_2,
    "z2_m": z_2,
    "xinf_m": x_inf,
    "yinf_m": y_inf,
    "zinf_m": z_inf
})
pd.options.display.float_format = "{:,.3f}".format
print(df.to_string(index=False))
print("\n")

# Plotting the blade shape - Top view

plt.figure(figsize=(10,6))
for i in range(len(Radii)):
    color = cm.viridis(i / len(Radii))
    plt.plot(y_1[i],z_1[i], 'o-',color = color)
    plt.plot(y_2[i],z_2[i], 's--',color = color)
    plt.plot(y_inf[i],z_inf[i],  'x',color = color)

plt.title("Blade shape Turbine 2 - Top view")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.axis('equal')
plt.grid(True)
plt.show()

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, projection='3d')

for i in range(len(Radii)):
    color = cm.viridis(i / len(Radii))
    
    x_points = [x_1[i], x_inf[i], x_2[i]]
    y_points = [y_1[i], y_inf[i], y_2[i]]
    z_points = [z_1[i], z_inf[i], z_2[i]]
    
    # interpolation spline
    t = np.linspace(0, 1, len(x_points))
    t_fine = np.linspace(0, 1, 50)  # plus de points pour une courbe douce

    # Interpolation PCHIP (préserve forme, évite oscillations)
    pchip_x = PchipInterpolator(t, x_points)(t_fine)
    pchip_y = PchipInterpolator(t, y_points)(t_fine)
    pchip_z = PchipInterpolator(t, z_points)(t_fine)
    
    # Tracé
    ax.plot(pchip_x, pchip_y, pchip_z, '-', color=color)
    ax.scatter(x_points, y_points, z_points, color=color, marker='o')

# Interpolation des lignes reliant tous les points 1, 2 et inf
x1 = np.array(x_1)
y1 = np.array(y_1)
z1 = np.array(z_1)
x2 = np.array(x_2)
y2 = np.array(y_2)
z2 = np.array(z_2)
xinf = np.array(x_inf)
yinf = np.array(y_inf)
zinf = np.array(z_inf)

t = np.linspace(0, 1, len(x1))
t_fine = np.linspace(0, 1, 200)

spline_x1 = make_interp_spline(t, x1, k=2)(t_fine)
spline_y1 = make_interp_spline(t, y1, k=2)(t_fine)
spline_z1 = make_interp_spline(t, z1, k=2)(t_fine)
spline_x2 = make_interp_spline(t, x2, k=2)(t_fine)
spline_y2 = make_interp_spline(t, y2, k=2)(t_fine)
spline_z2 = make_interp_spline(t, z2, k=2)(t_fine)
spline_xinf = make_interp_spline(t, xinf, k=2)(t_fine)
spline_yinf = make_interp_spline(t, yinf, k=2)(t_fine)
spline_zinf = make_interp_spline(t, zinf, k=2)(t_fine)

ax.plot(spline_x1, spline_y1, spline_z1, 'black', linewidth=2.0)
ax.plot(spline_x2, spline_y2, spline_z2, 'black', linewidth=2.0)
ax.plot(spline_xinf, spline_yinf, spline_zinf, 'black', linewidth=2.0)

# Paramètres du rotor (cylindre intérieur)
height = (max(z_1) - min(z_2) ) * 1.2  # hauteur approximative
theta = np.linspace(0, 2*np.pi, 60)
z_cyl = np.linspace(-height, height, 40)
theta, z_cyl = np.meshgrid(theta, z_cyl)

# Coordonnées du cylindre
x_cyl = R_i * np.cos(theta)
y_cyl = R_i * np.sin(theta)
z_cyl = z_cyl

# Tracé du cylindre
ax.plot_surface(x_cyl, y_cyl, z_cyl, color = 'black', alpha=0.5, linewidth=0, shade=True)



ax.set_title("Blade shape Turbine 3D view")
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.set_zlabel("Z (m)")
ax.set_xlim(-0.15, 0.15)
ax.set_ylim(-0.15, 0.15)
ax.set_zlim(-0.15, 0.15)
plt.show()