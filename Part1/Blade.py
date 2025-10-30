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

c = 1.0    # corde
m = 0.02   # cambrure max (4%)
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
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
# ======================================================
span = 2      # envergure
twist = 2.46    # torsion totale (degrés)
twist_axis = 0.25 * c  # quart de corde comme axe de torsion

z = np.linspace(0, span, 20)
x_3D, y_3D,z_3D  = [], [], []

for zi in z:
    # Angle de torsion local
    twist_angle = np.radians(twist * (zi / span))

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
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, color='royalblue', alpha=0.85, edgecolor='none')
ax.set_title('Torsion géométrique : corde réelle constante, projection variable')
ax.set_xlabel('X (Corde)')
ax.set_ylabel('Y')
ax.set_zlabel('Z (Envergure)')
ax.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(Z)))
plt.tight_layout()
plt.show()