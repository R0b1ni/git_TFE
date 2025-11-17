import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# ======================================================
def NACA4(m, p, t, c, x):
    yt = 5 * t * c * (
        0.2969 * np.sqrt(x/c)
        - 0.1260 * (x/c)
        - 0.3516 * (x/c)**2
        + 0.2843 * (x/c)**3
        - 0.1015 * (x/c)**4
    )

    # Ligne de cambrure
    yc = np.where(
        x < p*c,
        m/p**2 * (2*p*(x/c) - (x/c)**2) * c,
        m/(1 - p)**2 * ((1 - 2*p) + 2*p*(x/c) - (x/c)**2) * c
    )

    # Pente de la ligne de cambrure
    dyc_dx = np.where(
        x < p*c,
        2*m/p**2 * (p - x/c),
        2*m/(1 - p)**2 * (p - x/c)
    )
    theta = np.arctan(dyc_dx)

    # Coordonnées des surfaces supérieure et inférieure
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)

    return xu, yu, xl, yl

c = (0.155 - 0.062) * 1.1    # corde
m = 0.04   # cambrure max (4%)
p = 0.4    # position de la cambrure max (40%)
t = 0.12   # épaisseur max (12%)
x = np.linspace(0, c, 200)

xu, yu, xl, yl = NACA4(m, p, t, c, x)

plt.figure(figsize=(8, 4))
plt.plot(xu, yu, label='Extrados')
plt.plot(xl, yl, label='Intrados')
plt.title(f'NACA {int(m*100):01d}{int(p*10):01d}{int(t*100):02d}')
plt.xlabel('x (Corde)')
plt.ylabel('y')
plt.axis('equal')
plt.grid(False)
plt.legend()
plt.tight_layout()
plt.show()
# ======================================================
span = (0.155 - 0.062) * 1.1      # envergure
twist = -3.2    # torsion totale (degrés)
twist_axis = 0.25 * c  # quart de corde comme axe de torsion

z = np.linspace(0, span, 20)
x_3D, y_3D,z_3D  = [], [], []

for zi in z:
    # Angle de torsion local
    twist_angle = np.radians(twist * (zi / span)) - np.radians(73.73)  # on soustrait l'angle de sortie du fluide

    # On translate pour mettre l'axe de torsion à x=0 avant rotation
    xu_shifted = xu - twist_axis
    xl_shifted = xl - twist_axis

    # Rotation autour de l'axe de torsion (x=0.25*c)
    x_upper = xu_shifted * np.cos(twist_angle) - yu * np.sin(twist_angle)
    y_upper = xu_shifted * np.sin(twist_angle) + yu * np.cos(twist_angle)
    x_lower = xl_shifted * np.cos(twist_angle) - yl * np.sin(twist_angle)
    y_lower = xl_shifted * np.sin(twist_angle) + yl * np.cos(twist_angle)

    # On retransfère le profil à sa position originale
    x_upper += twist_axis
    x_lower += twist_axis

    # Empilage des points extrados + intrados
    x_section = np.concatenate((x_upper, x_lower[::-1]))
    y_section = np.concatenate((y_upper, y_lower[::-1]))
    z_section = np.full_like(x_section, zi)

    x_3D.append(x_section)
    y_3D.append(y_section)
    z_3D.append(z_section)

X = np.array(x_3D)
Y = np.array(y_3D)
Z = np.array(z_3D)
# Visualisation 3D
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(Z, Y, -X, rstride=1, cstride=1, color='royalblue', alpha=0.85, edgecolor='none')



# Paramètres du rotor (cylindre intérieur)
height = (0.155 - 0.062 )  # hauteur approximative
theta = np.linspace(0, 2*np.pi, 60)
z_cyl = np.linspace(-height, height, 40)
theta, z_cyl = np.meshgrid(theta, z_cyl)

# Coordonnées du cylindre
x_cyl = 0.062 * np.cos(theta) - 0.032
y_cyl = 0.062 * np.sin(theta) - 0.025
z_cyl = z_cyl + 0.031
# Tracé du cylindre
ax.plot_surface(x_cyl, y_cyl, -z_cyl, color = 'black', alpha=0.7, linewidth=0, shade=True)

ax.set_title('Torsion géométrique : corde réelle constante, projection variable')
ax.set_xlabel('X (Corde)')
ax.set_ylabel('Y')
ax.set_zlabel('Z (Envergure)')
#ax.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(Z)))
ax.set_xlim(-0.15, 0.15)
ax.set_ylim(-0.15, 0.15)
ax.set_zlim(-0.15, 0.15)
plt.tight_layout()
plt.show()