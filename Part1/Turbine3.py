import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import sys
from scipy.optimize import minimize
import matplotlib.cm as cm
from scipy.interpolate import make_interp_spline
from scipy.interpolate import PchipInterpolator
from matplotlib.patches import FancyArrowPatch

def draw_arrow(ax, x0, y0, x1, y1, color, label=None,
               lw=2.0, text_offset=(0, 0), mscale=12):
    arrow = FancyArrowPatch(
        (x0, y0), (x1, y1),
        arrowstyle='->',
        linewidth=lw,
        color=color,
        mutation_scale=mscale   # ⬅️ clé de la visibilité
    )
    ax.add_patch(arrow)

    if label is not None:
        ax.text(
            (x0 + x1) / 2 + text_offset[0],
            (y0 + y1) / 2 + text_offset[1],
            label,
            color=color,
            fontsize=12,
            ha='center',
            va='center'
        )

# Data
c_x_pipe = 2.5 # m/s (axial speed before the turbine)

rho_w = 1000 # kg/m3
g = 9.81 # m/s2
R_o = 0.1 # maximum diameter =  D20 = input diameter old value = 0.155
R_i = 0.4 * R_o # m (inner radius) (from fig 15 of 2022 Abeykoon)

#the fluid comes frome a pipe of R_o at speed c_x_pipe => by mass conservation
c_x = c_x_pipe * (R_o**2)/(R_o**2 - R_i**2) # en fonction de l'élargissement
print(f"Axial speed c_x: {c_x :.2f} m/s\n")

sigma = 1.6 # blockage ratio (from fig 15 of 2022 Abeykoon)
eta_r = 0.92 #theorical maximal value for a Kaplan turbine

Q = c_x * np.pi * (R_o**2 - R_i**2)
H = 1.0
P = rho_w * g * Q * H * eta_r
print(f"height of the water column: {H :.4f} m \n")
print (f"Power P: {P :.2f} W \n")

print(f"Volumetric Flow rate Q: {Q * 3600 :.4f} m3/h <=> {Q :.4f} m3/s \n")
#N = sigma * (2 * g * H)**(3/4) / (2*np.sqrt(np.pi * Q))
#N = 0.35 * N #(correction factor)
#n = N * 60 # rpm
n = 800
N = n / 60 # t/s

N_s = N * np.sqrt(Q) / (H**(3/4))

print(f"Hydraulic efficiency eta: {eta_r :.2f} \n")
print("Hydraulic Power P_hydraulic: {:.2f} W \n".format(rho_w * g * Q * H))
print(f"Rotationnal speed N: {n :.2f} rpm <=> {N :.2f} t/s \n")
print(f"Specific Rotationnal speed N_s: {N_s :.2f} \n")


Re = (rho_w * c_x * (R_o)) / (1e-3) # dynamic viscosity of water ~ 1e-3 Pa.s


#velocity triangles
U_m = 2 * np.pi * N * (R_i + R_o) / 2 # mean blade speed

Radii = np.linspace(R_i, R_o, 5)

c_1x = c_x
c1 = c_1x # assuming axial inflow
c_2x = c_x # mass conservation
c_3x = c_x
c3 = c_3x # assuming axial outflow

# 1D design
c_2teta = P/(rho_w * Q * U_m)

def c_teta(r):
    A = ((R_i + R_o) / 2) * c_2teta
    return A / r # free vortex distribution = constant work along the blade ?


w1 = []
w2 = []
w3 = []
w_inf = []
c2 = []
alpha_1 = []
alpha_2 = []
beta_1 = []
beta_2 = []
beta_3 = []
beta_inf = []
chord_to_pitch = [1.30, 1.163, 1.025, 0.888, 0.750] # à vérifier, valeurs du livre
pitch = []
chord = []

for i in range (len(Radii)):
    r = Radii[i]
    U = 2 * np.pi * N * r
    pitch.append( 2 * np.pi * r / 4) # assuming 4 blades
    chord.append( chord_to_pitch[i] * pitch[i] )
    c2.append( np.sqrt(c_2x**2 + c_teta(r)**2) )
    w1.append( np.sqrt(c_1x**2 + U**2) )
    w2.append( np.sqrt((U - c_teta(r))**2 + c_2x**2) )
    w3.append( np.sqrt(c_3x**2 + U**2) )
    w_inf.append( (w3[i] + w2[i]) / 2 )
    
    alpha_1.append(89)
    alpha_2.append( 90 + np.degrees(np.arctan2(c_2teta, c_2x)) )
    beta_1.append( np.degrees(np.arctan2(U, c_1x)) +90)
    beta_2.append( np.degrees(np.arctan2( np.abs(U - c_teta(r)), c_2x)) +90 )
    beta_3.append( np.degrees(np.arctan2(U, c_3x)) +90)
    beta_inf.append(( beta_3[i] + beta_2[i] ) /2 )
    if r == (R_i + R_o) / 2:
        c_2x = float(c_2x)
        c_2y = float(c_teta(r))
        c_3x = float(c_3x)
        c_3y = 0.0

        color1 = 'tab:blue'
        color2 = 'tab:gray'
        color3 = 'lightskyblue'

        fig, ax = plt.subplots(figsize=(10,8))

        # V2
        draw_arrow(ax, 0, 0, -c_2y, -c_2x, color1, 'C2', lw=2.5, text_offset=(0, 0.15))

        # U (entrée)
        draw_arrow(ax, -c_2y + U_m, -c_2x, -c_2y, -c_2x,
                color2, 'U', lw=2.0, text_offset=(0, 0.15))

        # W2
        draw_arrow(ax, 0, 0, -c_2y + U_m, -c_2x,
                color3, 'W2', lw=2.5, text_offset=(0, -0.2))

        # V3
        draw_arrow(ax, 0, 0, c_3y, -c_3x, color1, 'C3', lw=2.5, text_offset=(0.15, 0))

        # U (sortie)
        draw_arrow(ax, c_3y + U_m, -c_3x, c_3y, -c_3x,
                color2, 'U', lw=2.0, text_offset=(0, 0.15))

        # W3
        draw_arrow(ax, 0, 0, c_3y + U_m, -c_3x,
                color3, 'W3', lw=2.5, text_offset=(0.2, 0))

        ax.set_title('Velocity Triangles', fontsize=14)
        ax.scatter([0], [0], color='black', zorder=5)

        #ax.set_aspect('equal')

        #ax.set_xlim(-0.25, 9)
        # ax.set_ylim(-2.5, 0)

        ax.grid(True, linestyle='--', alpha=0.4)
        plt.tight_layout()
        #plt.savefig('Part1/results/velocity_triangles.pdf')
        plt.show()



df = pd.DataFrame({
    "radius_m": Radii,
    "U_m_s": [2 * np.pi * N * r for r in Radii],
    "w1_m_s": w1,
    "w2_m_s": w2,
    "beta1_deg": beta_1,
    "beta2_deg": beta_2,
    "beta3_deg": beta_3,
    "c2_m_s": c2,
    "alpha1_deg": alpha_1,
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

xs_1 = []
ys_1 = []
zs_1 = []
xs_2 = []
ys_2 = []
zs_2 = []
xs_inf = []
ys_inf = []
zs_inf = []
alpha_inf = []

xs_1_new = []
ys_1_new = []
zs_1_new = []
xs_2_new = []
ys_2_new = []
zs_2_new = []
xs_inf_new = []
ys_inf_new = []
zs_inf_new = []

def delta_r (z):
    alpha = 0.1
    return alpha * z

for i in range (len(Radii)):
    D = 2 * Radii[i]
    
    y_1.append( np.sin( np.radians(beta_2[i] - 90)) * ((D * np.sin(np.radians((chord[i]/D )* (180/np.pi))))/(2 * np.cos(np.radians(beta_inf[i] - beta_2[i])))) )
    z_1.append(y_1[i] / np.tan( np.radians(beta_2[i] - 90) ) )
    x_1.append((D/2) - D*np.sin( np.arcsin(2*y_1[i]/D)/2 )**2 )
    
    y_2.append(-np.sin( np.radians(beta_3[i] - 90)) * ((D * np.sin(np.radians((chord[i]/D )* (180/np.pi))))/(2 * np.cos(np.radians(beta_3[i] - beta_inf[i])))) )
    z_2.append(y_2[i] / np.tan( np.radians(beta_3[i] - 90) ) )
    x_2.append((D/2) - D*np.sin( np.arcsin(2*y_2[i]/D)/2 )**2 )
    
    x_inf.append( D/2)
    y_inf.append( 0)
    z_inf.append( 0 )
    
    
    offset = 0.1 # m, distance between the stator and the rotor in height
    alpha_inf.append((alpha_1[i] + alpha_2[i]) / 2)
    
    ys_1.append( -np.sin( np.radians(alpha_1[i] - 90)) * ((D * np.sin(np.radians((chord[i]/D )* (180/np.pi))))/(2 * np.cos(np.radians(alpha_inf[i] - alpha_1[i])))) )
    zs_1.append(-ys_1[i] / np.tan( np.radians(alpha_1[i] - 90) ) +offset)
    xs_1.append((D/2) - D*np.sin( np.arcsin(2*ys_1[i]/D)/2 )**2 )
    
    ys_2.append(np.sin( np.radians(alpha_2[i] - 90)) * ((D * np.sin(np.radians((chord[i]/D )* (180/np.pi))))/(2 * np.cos(np.radians(alpha_2[i] - alpha_inf[i])))) )
    zs_2.append(-ys_2[i] / np.tan( np.radians(alpha_2[i] - 90) ) + offset)
    xs_2.append((D/2) - D*np.sin( np.arcsin(2*ys_2[i]/D)/2 )**2 )
    
    xs_inf.append( D/2)
    ys_inf.append( 0)
    zs_inf.append(offset)
    
    r_old = Radii[i]
    r_new = np.sqrt(r_old**2 + np.abs(delta_r(zs_2[i]) * (2 * Radii[0] - delta_r(zs_2[i])) ))
    
    xs_inf_new.append(r_new)
    ys_inf_new.append( 0)
    zs_inf_new.append(offset)
    
    print(zs_1[i],z_inf[i], zs_2[i])
    r_new = np.sqrt(r_old**2 + np.abs(delta_r(zs_1[i] + zs_2[i]) * (2 * Radii[0] + delta_r(zs_1[i] + zs_2[i]))) )
    teta = np.arctan2(xs_1[i], ys_1[i])
    xs_1_new.append( (r_new) * np.sin(teta) )
    ys_1_new.append( (r_new) * np.cos(teta) )
    zs_1_new.append(zs_1[i])
    
    teta = np.arctan2(xs_2[i], ys_2[i])
    xs_2_new.append( (Radii[i]) * np.sin(teta) )
    ys_2_new.append( (Radii[i]) * np.cos(teta) )
    zs_2_new.append(zs_2[i])


xs_1 = xs_1_new
ys_1 = ys_1_new
zs_1 = zs_1_new
xs_2 = xs_2_new
ys_2 = ys_2_new
zs_2 = zs_2_new
xs_inf = xs_inf_new
ys_inf = ys_inf_new
zs_inf = zs_inf_new    
# Plotting the blade shape - Top view

'''plt.figure(figsize=(10,6))
for i in range(len(Radii)):
    color = cm.viridis(i / len(Radii))
    plt.plot([y_1[i], y_inf[i]], [z_1[i], z_inf[i]], color='black', linewidth=0.2)
    plt.plot([y_inf[i], y_2[i]], [z_inf[i], z_2[i]], color='black', linewidth=0.2)
    plt.plot(y_1[i],z_1[i], 'o-',color = 'red')
    plt.plot(y_2[i],z_2[i], 's--',color = 'red')
    plt.plot(y_inf[i],z_inf[i],  'x',color = 'red')

plt.title("Blade shape Turbine 2 - Top view")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.axis('equal')
plt.grid(False)
#plt.savefig('Part1/results/Twist_view.pdf')
plt.show()'''

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, projection='3d')

for i in range(len(Radii)):
    color = cm.viridis(i / len(Radii))
    
    x_points = [x_1[i], x_inf[i], x_2[i]]
    y_points = [y_1[i], y_inf[i], y_2[i]]
    z_points = [z_1[i], z_inf[i], z_2[i]]
    
    xs_points = [xs_1[i], xs_inf[i], xs_2[i]]
    ys_points = [ys_1[i], ys_inf[i], ys_2[i]]
    zs_points = [zs_1[i], zs_inf[i], zs_2[i]]
    
    # interpolation spline
    t = np.linspace(0, 1, len(x_points))
    t_fine = np.linspace(0, 1, 50)  # plus de points pour une courbe douce

    # Interpolation PCHIP (préserve forme, évite oscillations)
    pchip_x = PchipInterpolator(t, x_points)(t_fine)
    pchip_y = PchipInterpolator(t, y_points)(t_fine)
    pchip_z = PchipInterpolator(t, z_points)(t_fine)
    
    pchip_xs = PchipInterpolator(t, xs_points)(t_fine)
    pchip_ys = PchipInterpolator(t, ys_points)(t_fine)
    pchip_zs = PchipInterpolator(t, zs_points)(t_fine)
    
    # Tracé
    ax.plot(pchip_x, pchip_y, pchip_z, '-', color='black')
    ax.plot(pchip_xs, pchip_ys, pchip_zs, '-', color='black')
    ax.scatter(x_points, y_points, z_points, color='red', marker='o', s=30, alpha=0.7)
    ax.scatter(xs_points, ys_points, zs_points, color='red', marker='o', s=30, alpha=0.7)

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

xs1 = np.array(xs_1)
ys1 = np.array(ys_1)
zs1 = np.array(zs_1)
xs2 = np.array(xs_2)
ys2 = np.array(ys_2)
zs2 = np.array(zs_2)
xsinf = np.array(xs_inf)
ysinf = np.array(ys_inf)
zsinf = np.array(zs_inf)

spline_xs1 = make_interp_spline(t, xs1, k=2)(t_fine)
spline_ys1 = make_interp_spline(t, ys1, k=2)(t_fine)
spline_zs1 = make_interp_spline(t, zs1, k=2)(t_fine)
spline_xs2 = make_interp_spline(t, xs2, k=2)(t_fine)
spline_ys2 = make_interp_spline(t, ys2, k=2)(t_fine)
spline_zs2 = make_interp_spline(t, zs2, k=2)(t_fine)
spline_xsinf = make_interp_spline(t, xsinf, k=2)(t_fine)
spline_ysinf = make_interp_spline(t, ysinf, k=2)(t_fine)
spline_zsinf = make_interp_spline(t, zsinf, k=2)(t_fine)

ax.plot(spline_xs1, spline_ys1, spline_zs1, 'black', linewidth=2.0)
ax.plot(spline_xs2, spline_ys2, spline_zs2, 'black', linewidth=2.0)
ax.plot(spline_xsinf, spline_ysinf, spline_zsinf, 'black', linewidth=2.0)

# Paramètres du rotor (cylindre intérieur)
height = (max(zs_1) - min(z_2) ) * 1.2  # hauteur approximative
theta = np.linspace(0, 2*np.pi, 60)
z_cyl = np.linspace(-height/4, height, 40)
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
ax.set_xlim(-0.1, 0.1)
ax.set_ylim(-0.1, 0.1)
ax.set_zlim(-0.1, 0.1)
#plt.savefig('Part1/results/3D_Model.pdf')
plt.show()
