import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import sys
from scipy.optimize import minimize
import matplotlib.cm as cm

#NACA 4 digit airfoil generator
def NACA4(m,p,t,c,x):
    """
    NACA 4 digit airfoil generator
    Inputs:
    m : maximum camber (as fraction of chord)
    p : location of maximum camber (as fraction of chord)
    t : maximum thickness (as fraction of chord)
    c : chord length
    x : x-coordinates along the chord (array)

    Outputs:
    yu : y-coordinates of upper surface
    yl : y-coordinates of lower surface
    """

    # Thickness distribution
    yt = 5 * t * c * (0.2969 * np.sqrt(x/c) - 0.1260 * (x/c) - 0.3516 * (x/c)**2 + 0.2843 * (x/c)**3 - 0.1015 * (x/c)**4)

    # Camber line and its slope
    yc = np.where(x < p*c,
                  m/p**2 * (2*p*(x/c) - (x/c)**2) * c,
                  m/(1 - p)**2 * ((1 - 2*p) + 2*p*(x/c) - (x/c)**2) * c)

    dyc_dx = np.where(x < p*c,
                      2*m/p**2 * (p - x/c),
                      2*m/(1 - p)**2 * (p - x/c))

    theta = np.arctan(dyc_dx)

    # Upper and lower surface coordinates
    yu = yc + yt * np.cos(theta)
    yl = yc - yt * np.cos(theta)

    return yu, yl


    
c = 1.0  # Chord length
m = 0.02  # Maximum camber (2% of chord)
p = 0.4   # Location of maximum camber (40% of chord)
t = 0.12  # Maximum thickness (12% of chord)

x = np.linspace(0, c, 100)
yu, yl = NACA4(m, p, t, c, x)

plt.figure()
plt.plot(x, yu, label='Upper Surface')
plt.plot(x, yl, label='Lower Surface')
plt.title('NACA 4-Digit Airfoil: NACA {}{}{}'.format(int(m*100), int(p*10), int(t*100)))
plt.xlabel('x (Chord Length)')
plt.ylabel('y')
plt.axis('equal')
plt.grid(True)
plt.legend()
plt.show()

#3D plot of the airfoil extruded in the spanwise direction
from mpl_toolkits.mplot3d import Axes3D

span = 2  # Spanwise length
twist = 15  # Twist angle in degrees
z = np.linspace(0, span, 10)
x_3D,y_3D = [],[]

for zi in z:
    twist_angle = np.radians(twist * (zi/span))
    x_upper = x * np.cos(twist_angle) - yu * np.sin(twist_angle)
    y_upper = x * np.sin(twist_angle) + yu * np.cos(twist_angle)
    x_lower = x * np.cos(twist_angle) - yl * np.sin(twist_angle)
    y_lower = x * np.sin(twist_angle) + yl * np.cos(twist_angle)
    x_3D.append(np.concatenate((x_upper, x_lower[::-1])))
    y_3D.append(np.concatenate((y_upper, y_lower[::-1])))
x_3D = np.array(x_3D)
y_3D = np.array(y_3D)
X, Z = np.meshgrid(np.concatenate((x, x[::-1])), z)
Y = y_3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, color='b', alpha=0.7)
ax.set_title('3D View of NACA 4-Digit Airfoil with Twist')
ax.set_xlabel('X (Chord Length)')
ax.set_ylabel('Y')
ax.set_zlabel('Z (Spanwise Length)')
# Make 3D axes have equal scale
try:
    ax.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(Z)))
except Exception:
    x_range = X.max() - X.min()
    y_range = Y.max() - Y.min()
    z_range = Z.max() - Z.min()
    max_range = max(x_range, y_range, z_range)
    x_mid = 0.5 * (X.max() + X.min())
    y_mid = 0.5 * (Y.max() + Y.min())
    z_mid = 0.5 * (Z.max() + Z.min())
    ax.set_xlim(x_mid - max_range/2, x_mid + max_range/2)
    ax.set_ylim(y_mid - max_range/2, y_mid + max_range/2)
    ax.set_zlim(z_mid - max_range/2, z_mid + max_range/2)
plt.show()