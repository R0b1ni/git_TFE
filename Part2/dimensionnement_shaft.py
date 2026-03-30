"""
MÉMOIRE: DIMENSIONNEMENT SHAFT

Created on Mon Feb  9 17:08:39 2026

@author: H. Nishio et R. Thonon

(IA utilisée pour la génération des plots)
"""

import numpy as np
import math as m
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


# Paramètres globaux
w   = 900       #[RPM]  vitesse de rotation de l'hélice
S_y = 210e6     #[Pa]   sigma_yield acier inox

#SF  = 1         #[-]    facteur de sécurité


# Shaft Satellite-Carrier  (STATIQUE)------------------------------------------

a_sc        = 23.5*1e-3   #[m]   distance hélice-roulement de gauche
b_sc        = 38*1e-3     #[m]   distance entre les 2 roulements
c_sc        = 18*1e-3     #[m]   distance entre le roulement de gauche et le satellite-carrier
T_h         = 13          #[Nm]  couple hélice
m_pale      = 25          #[g]   poids d'une pale
m_coeur     = 343         #[g]   poids du coeur de l'hélice
m_sc        = 115         #[g]   poids satellite-carrier

# Forces [N]
F_a_sc = (4*m_pale+m_coeur)*1e-3*9.81 
F_d_sc = m_sc*1e-3*9.81

F_c_sc = ((b_sc+c_sc)*F_d_sc-a_sc*F_a_sc)/b_sc
F_b_sc = F_a_sc + F_d_sc - F_c_sc

M_max_sc = a_sc*F_a_sc

k_sc = ((32*M_max_sc/np.pi)**2) + 3*((16*T_h/np.pi)**2)

d_min_sc = (m.sqrt(k_sc)/(S_y))**(1/3)

#print("Diamètre minimum shaft SC (statique): ", d_min_sc*1e3, " mm.")

#Shaft Satellite-Carrier - DIAGRAMME EFFORTS TRANCHANTS------------------------

#position
pos_final = a_sc+b_sc+c_sc
x_sc = [0, 0, pos_final-b_sc-c_sc, pos_final-c_sc, pos_final]

#effort tranchant 
F1 = 0
F2 = F_a_sc
F3 = F_a_sc - F_b_sc
F4 = F3 - F_c_sc
F5 = 0
y_sc = [F1, F2, F3, F4, F5] 

fig, ax = plt.subplots(figsize=(8, 4))

# plt.step trace les lignes horizontales et verticales.
# where='post' indique que la valeur y reste constante jusqu'au prochain point x, puis "saute" verticalement.
ax.step(x_sc, y_sc, where='post', color='blue', linewidth=2)

# Remplissage sous la courbe 
ax.fill_between(x_sc, y_sc, step='post', alpha=0.3, color='blue')

# Trace la ligne neutre à y=0
ax.axhline(0, color='black', linewidth=1.5)

ax.set_xticks(x_sc)
ax.set_xticklabels(['A', 'A', 'B', 'C', 'D'], fontsize=12, fontweight='bold')
ax.set_yticks([F2, 0, F3, F1])
#ax.set_yticklabels(['F2', '0', 'F3', 'F1'], fontsize=12, fontweight='bold')
ax.invert_yaxis()

ax.set_title("Shearing diagram", pad=15, fontsize=14) 
ax.set_ylabel("Shearing [N]")
ax.set_xlabel("Position along the runner-driven shaft")

plt.tight_layout()

#plt.savefig("tranchantF_sc.pdf")
plt.show()

#Shaft SC - DIAGRAMME MOMENTS--------------------------------------------------

M_sc = [0, 0, M_max_sc, 0.01, 0]

plt.plot(x_sc, M_sc, color='purple')
plt.axhline(0, color='black', linewidth=1.5)
plt.fill_between(x_sc, M_sc, alpha=0.3, color='purple')
plt.title("Bending diagram")
plt.xlabel("Position along the runner-driven shaft")
plt.xticks(x_sc, ['A', 'A', 'B', 'C', 'D'], fontsize=12, fontweight='bold')
plt.ylabel("Moments [Nm]")
plt.gca().invert_yaxis()

#plt.savefig("moment_sc.pdf")
plt.show()

# Shaft Satellite-Carrier  (FATIGUE)-------------------------------------------
 
S_u      = 634e6    #[Pa]
C_L      = 0.58
C_G      = 1
C_S      = 0.76
C_T      = 1
C_R      = 0.814
Sprime_N = 0.5*S_u
S_n      = Sprime_N*C_L*C_G*C_S*C_T*C_R

q    = 0.72 

K_tb_sc = 1.7
K_fb_sc = 1 + (K_tb_sc - 1)*q
sigma_ea_sc = lambda d : 32*M_max_sc*K_fb_sc/(np.pi*d**3)

K_tt_sc = 1.4
K_ft_sc = 1 + (K_tt_sc - 1)*q
sigma_em_sc = lambda d : 16*T_h*K_ft_sc/(np.pi*d**3)

f_sigma_sc = lambda x : 32*M_max_sc*K_fb_sc*x/(16*T_h*K_ft_sc)
f_SnSu_sc = lambda x : S_n*(1-x/S_u)

f_diff_sc = lambda x : f_sigma_sc(x) - f_SnSu_sc(x)
guess_xP_sc = 100e6
x_P_sc = fsolve(f_diff_sc, guess_xP_sc)[0]
y_P_sc = f_sigma_sc(x_P_sc)

d_sc_fatigue = (K_ft_sc*16*T_h/(m.pi*x_P_sc))**(1/3)

#print("Diamètre minimum shaft SC (fatigue): ", d_sc_fatigue*1e3, " mm")

plt.figure(figsize=(8, 6))

# Droite 1 : relie (0, S_n) à (S_u, 0)
plt.plot([0, S_u], [S_n, 0], color='red')

# Droite 2 : relie (0, C) à (C, 0)
plt.plot([0, S_y], [S_y, 0], color='blue')

# Droite 3 : passe par l'origine avec une pente m
# Pour la tracer joliment, on la fait aller de x=0 jusqu'au point x le plus éloigné (B ou C)
pente_sc = 32*M_max_sc*K_fb_sc/(16*T_h*K_ft_sc)
x_max = max(S_u, S_y) 
plt.plot([0, x_max], [0, pente_sc * x_max], color='green')

# Tracer le point (x_p, y_p)
plt.plot(x_P_sc, y_P_sc, 'ko', markersize=8, zorder=10)
# Tracer la ligne pontillée verticale (vers l'axe X)
plt.plot([x_P_sc, x_P_sc], [0, y_P_sc], 'k--', linewidth=1, alpha=0.7)
# Tracer la ligne pointillée horizontale (vers l'axe Y)
plt.plot([0, x_P_sc], [y_P_sc, y_P_sc], 'k--', linewidth=1, alpha=0.7)

# Tracer le point (sigma_em, sigma_ea)
plt.plot(sigma_em_sc(d_min_sc), sigma_ea_sc(d_min_sc), 'ko', markersize=8, zorder=10)
# Tracer la ligne pointillée verticale (vers l'axe X)
plt.plot([sigma_em_sc(d_min_sc), sigma_em_sc(d_min_sc)], [0, sigma_ea_sc(d_min_sc)], 'k--', linewidth=1, alpha=0.7)
# Tracer la ligne pointillée horizontale (vers l'axe Y)
plt.plot([0, sigma_em_sc(d_min_sc)], [sigma_ea_sc(d_min_sc), sigma_ea_sc(d_min_sc)], 'k--', linewidth=1, alpha=0.7)

# --- 3. Mise en forme ---
coord_x_sc = [0, sigma_em_sc(d_min_sc),S_y, x_P_sc, S_u]
coord_y_sc = [0, sigma_ea_sc(d_min_sc), S_y, y_P_sc, S_n]
graduation_x_sc = [0, r"$\sigma_{em}$",r"$S_y$", r"$x_p$", r"$S_u$"]
graduation_y_sc = [" ", r"$\sigma_{ea}$", r"$S_y$", r"$y_p$", r"$S_n$"]
plt.title("Goodman diagram of the runner-driven shaft")
plt.xticks(coord_x_sc, graduation_x_sc, fontsize=12, fontweight='bold')
plt.yticks(coord_y_sc, graduation_y_sc, fontsize=12, fontweight='bold')
plt.xlabel("Mean stress")
plt.ylabel("Alternating stress")

# On force l'affichage du point (0,0) en bas à gauche pour bien voir l'origine
plt.xlim(left=0)
plt.ylim(bottom=0)

plt.grid(True, linestyle='--', alpha=0.6)

#plt.savefig("goodman_sc.pdf")
plt.show()

# Shaft Rotor  (STATIQUE)------------------------------------------------------

a_r        = 21.5*1e-3            #[m]      distance soleil-roulement gauche
b_r        = 41.75*1e-3           #[m]      distance centre générateur-roulement gauche
c_r        = b_r                  #[m]      distance centre générateur-roulement droit
T_s        = 0.95*T_h/4           #[Nm]     couple soleil
L_m        = 6.15*1e-3            #[m]      longueur aimant
l_m        = 1.5*1e-3             #[m]      largeur aimant
h_m        = 0.5*1e-3             #[m]      épaisseur aimant
V_m        = L_m*l_m*h_m          #[m³]     volume aimant
rho        = 7.4                  #[kg/m³]  densité neodyme
m_aimant   = V_m*rho*9.81*1e3     #[g]      poids d'un aimant
m_rotor    = 343                  #[g]      poids carcasse rotor
m_soleil    = 115                  #[g]      poids soleil

# Forces [N]
F_a_r = m_soleil*1e-3*9.81 
F_c_r = (12*m_aimant+m_rotor)*1e-3*9.81

F_d_r = (b_r*F_c_r - a_r*F_a_r)/(b_r+c_r)
F_b_r = F_a_r + F_c_r - F_d_r

M_max_r = c_r*F_d_r

k_r = ((32*M_max_r/np.pi)**2) + 3*((16*T_s/np.pi)**2)

d_min_r = (m.sqrt(k_r)/(S_y))**(1/3)

#print("Diamètre minimum shaft rotor (statique): ", d_min_r*1e3, " mm.")

#Shaft rotor - DIAGRAMME EFFORTS TRANCHANTS------------------------------------

#position
pos_final = a_r+b_r+c_r
x_r = [0, 0, pos_final-b_r-c_r, pos_final-c_r, pos_final]

#effort tranchant 
F1 = 0
F2 = F_a_r
F3 = F_a_r - F_b_r
F4 = F3 + F_c_r
F5 = 0
y_r = [F1, F2, F3, F4, F5] 

fig, ax = plt.subplots(figsize=(8, 4))

# plt.step trace les lignes horizontales et verticales.
# where='post' indique que la valeur y reste constante jusqu'au prochain point x, puis "saute" verticalement.
ax.step(x_r, y_r, where='post', color='blue', linewidth=2)

# Remplissage sous la courbe 
ax.fill_between(x_r, y_r, step='post', alpha=0.3, color='blue')

# Trace la ligne neutre à y=0
ax.axhline(0, color='black', linewidth=1.5)

ax.set_xticks(x_r)
ax.set_xticklabels(['A', 'A', 'B', 'C', 'D'], fontsize=12, fontweight='bold')
ax.set_yticks([F2, 0, F3, F1])
#ax.set_yticklabels(['F2', '0', 'F3', 'F1'], fontsize=12, fontweight='bold')
ax.invert_yaxis()

ax.set_title("Shearing diagram", pad=15, fontsize=14)
ax.set_ylabel("Shearing [N]")
ax.set_xlabel("Position along the rotoric shaft")

plt.tight_layout()

#plt.savefig("tranchant_r.pdf")
plt.show()

#Shaft rotor - DIAGRAMME MOMENTS-----------------------------------------------
M_b_r = a_r*F_a_r
M_r = [0, 0, M_b_r, M_max_r, 0]

plt.plot(x_r, M_r, color='purple')
plt.axhline(0, color='black', linewidth=1.5)
plt.title("Bending diagram")
plt.xlabel("Position along the rotoric shaft")
plt.xticks(x_r, ['A', 'A', 'B', 'C', 'D'], fontsize=12, fontweight='bold')
plt.ylabel("Moments [Nm]")
plt.fill_between(x_r, M_r, alpha=0.3, color='purple')
plt.gca().invert_yaxis()
#plt.savefig("moment_r.pdf")
plt.show()

# Shaft Rotor  (FATIGUE)-------------------------------------------------------

K_tb_r = 2.05
K_fb_r = 1 + (K_tb_r - 1)*q
sigma_ea_r = lambda d : 32*M_max_r*K_fb_r/(np.pi*d**3)

K_tt_r = 1.63
K_ft_r = 1 + (K_tt_r - 1)*q
sigma_em_r = lambda d : 16*T_s*K_ft_r/(np.pi*d**3)

f_sigma_r = lambda x : 32*M_max_r*K_fb_r*x/(16*T_s*K_ft_r)
f_SnSu_r = lambda x : S_n*(1-x/S_u)

f_diff_r = lambda x : f_sigma_r(x) - f_SnSu_r(x)
guess_xP_r = 100e6
x_P_r = fsolve(f_diff_r, guess_xP_r)[0]
y_P_r = f_sigma_r(x_P_r)
d_r_fatigue = (K_ft_r*16*T_h/(m.pi*x_P_r))**(1/3)

#print("Diamètre minimum shaft rotor (fatigue): ", d_r_fatigue*1e3, " mm")

plt.figure(figsize=(8, 6))

# Droite 1 : relie (0, S_n) à (S_u, 0)
plt.plot([0, S_u], [S_n, 0], color='red')

# Droite 2 : relie (0, C) à (C, 0)
plt.plot([0, S_y], [S_y, 0], color='blue')

# Droite 3 : passe par l'origine avec une pente m
# Pour la tracer joliment, on la fait aller de x=0 jusqu'au point x le plus éloigné (B ou C)
pente_r = 32*M_max_r*K_fb_r/(16*T_s*K_ft_r)
x_max = max(S_u, S_y) 
plt.plot([0, x_max], [0, pente_r * x_max], color='green')

# Tracer le point (x_p, y_p)
plt.plot(x_P_r, y_P_r, 'ko', markersize=8, zorder=10)
# Tracer la ligne pointillée verticale (vers l'axe X)
plt.plot([x_P_r, x_P_r], [0, y_P_r], 'k--', linewidth=1, alpha=0.7)
# Tracer la ligne pointillée horizontale (vers l'axe Y)
plt.plot([0, x_P_r], [y_P_r, y_P_r], 'k--', linewidth=1, alpha=0.7)

# Tracer le point (sigma_em, sigma_ea)
plt.plot(sigma_em_r(d_min_r), sigma_ea_r(d_min_r), 'ko', markersize=8, zorder=10)
# Tracer la ligne pointillée verticale (vers l'axe X)
plt.plot([sigma_em_r(d_min_r), sigma_em_r(d_min_r)], [0, sigma_ea_r(d_min_r)], 'k--', linewidth=1, alpha=0.7)
# Tracer la ligne pointillée horizontale (vers l'axe Y)
plt.plot([0, sigma_em_r(d_min_r)], [sigma_ea_r(d_min_r), sigma_ea_r(d_min_r)], 'k--', linewidth=1, alpha=0.7)

# --- 3. Mise en forme ---
coord_x_r = [0, sigma_em_r(d_min_r),S_y, x_P_r, S_u]
coord_y_r = [0, sigma_ea_r(d_min_r), S_y, y_P_r, S_n]
graduation_x_r = [0, r"$\sigma_{em}$",r"$S_y$", r"$x_p$", r"$S_u$"]
graduation_y_r = [0, r"$\sigma_{ea}$", r"$S_y$", r"$y_p$", r"$S_n$"]
plt.title("Goodman diagram of the rotoric shaft")
plt.xticks(coord_x_r, graduation_x_r, fontsize=12, fontweight='bold')
plt.yticks(coord_y_r, graduation_y_r, fontsize=12, fontweight='bold')
plt.xlabel("Mean stress")
plt.ylabel("Alternating stress")

# On force l'affichage du point (0,0) en bas à gauche pour bien voir l'origine
plt.xlim(left=0)
plt.ylim(bottom=0)

plt.grid(True, linestyle='--', alpha=0.6)

#plt.savefig("goodman_r.pdf")
plt.show()

#CHOIX FINAL ------------------------------------------------------------------

d_sc = max(d_min_sc, d_sc_fatigue)

if d_min_sc > d_sc_fatigue : print("\nDiamètre minimum shaft SC    : ", round(d_sc*1e3, 3), " mm (point déterminant = PLASTICITÉ (statique))")
else : print("\nDiamètre minimum shaft SC    : ", round(d_sc*1e3, 3), " mm (mode déterminant = FATIGUE)")

d_r = max(d_min_r, d_r_fatigue)

if d_min_r > d_r_fatigue : print("Diamètre minimum shaft rotor : ", round(d_r*1e3, 3), " mm (point déterminant = PLASTICITÉ (statique))")
else : print("Diamètre minimum shaft rotor : ", round(d_r*1e3, 3), " mm (mode déterminant = FATIGUE)")
 
