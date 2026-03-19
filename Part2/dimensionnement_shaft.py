"""
MÉMOIRE: DIMENSIONNEMENT SHAFT

Created on Mon Feb  9 17:08:39 2026

@author: H. Nishio et R. Thonon
"""

import numpy as np
import math as m

# Paramètres globaux
w   = 900       #[RPM]  vitesse de rotation de l'hélice
S_y = 210e6     #[Pa]   sigma_yield acier inox

#SF  = 1         #[-]    facteur de sécurité


# Shaft Satellite-Carrier -----------------------------------------------------

a_sc        = 23.5*1e-3   #[m]   distance hélice-roulement de gauche
b_sc        = 38*1e-3     #[m]   distance entre les 2 roulements
c_sc        = 18*1e-3     #[m]   distance entre le roulement de gauche et le satellite-carrier
T_h         = 13          #[Nm]  couple hélice
m_pale      = 25          #[g]   poids d'une pale
m_coeur     = 343         #[g]   poids du coeur de l'hélice
m_SC        = 115         #[g]   poids satellite-carrier

# Forces [N]
F_a = (4*m_pale+m_coeur)*1e-3*9.81 
F_d = m_SC*1e-3*9.81

F_c = ((b_sc+c_sc)*F_d-a_sc*F_a)/b_sc
F_b = F_a + F_d - F_c

M_max = a_sc*F_a

k = ((32*M_max/np.pi)**2) + 3*((16*T_h/np.pi)**2)

d_min = (m.sqrt(k)/(S_y))**(1/3)

print("Diamètre minimum shaft SC: ", d_min*1e3, " mm.")
 

# Shaft Rotor -----------------------------------------------------------------

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
F_a = m_soleil*1e-3*9.81 
F_c = (12*m_aimant+m_rotor)*1e-3*9.81

F_d = (b_r*F_c - a_r*F_a)/(b_r+c_r)
F_b = F_a + F_c - F_d

M_max = c_r*F_d

k = ((32*M_max/np.pi)**2) + 3*((16*T_s/np.pi)**2)

d_min = (m.sqrt(k)/(S_y))**(1/3)

print("Diamètre minimum saft rotor: ", d_min*1e3, " mm.")
 
