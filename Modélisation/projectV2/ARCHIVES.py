import numpy as np
from numba import njit
from initialisation.positions import random_pos, pos_cristal2D
from initialisation.vitesses import vit_temp, random_vit
from dyna.dynam import verlet, LJ_walls
import time


#########################################################################################################
#####################ANCIEN CALCUL DE DLJPot pour une particule et pour n particules#####################
#########################################################################################################

@njit
def dLJpot_i(r, i, sig, eps, cut_off, nb_part, D):
    
    """
    Gradient du potentiel de Lennard-Jones pour une particule

    Args:
        r (array): positions des particules
        i (float): index de la particule
        sig (int): paramètre du puit de LJ
        eps (float): profondeur du puit de LJ
        cut_off (int): distance a partir de laquelle la force est nulle
        nb_part (int): nombre de particules
        D (int): nombre de dimensions

    Returns:
        dLJP (array): Gradient du potentiel de la particule
    """

    incr = 0
    dpart_cut = np.empty((nb_part-1,D), dtype=float)
    dpart_calc_cut = np.empty((nb_part-1,1), dtype=float)

    for j in range(nb_part):
        if j != i:
            dist_xy = r[j] - r[i]
            dist = np.sqrt(dist_xy[0]**2 + dist_xy[1]**2)
            if dist <= cut_off:
                dpart_cut[incr] = dist_xy   #on ajoute cet élément dans l'array
                dpart_calc_cut[incr] = dist
                incr +=1

    dpart_cut = dpart_cut[:incr]
    dpart_calc_cut = dpart_calc_cut[:incr]

    r8 =  (sig**6)*(1.0/dpart_calc_cut)**8
    r14 = 2.0*(sig**12)*(1.0/dpart_calc_cut)**14
    r814 = r14-r8
    r814v = dpart_cut*r814
    dLJP_i = -24.0*eps*np.sum(r814v,axis=0)
    
    return dLJP_i



@njit
def dLJpot(r0, sig, eps, cutoff, nb_part, D):
    """
    Gradient du potentiel de Lennard-Jones pour n particules

    Args:
        r (array): positions des particules
        sig (int): paramètre du puit de LJ
        eps (float): profondeur du puit de LJ
        cut_off (int): distance a partir de laquelle la force est nulle
        nb_part (int): nombre de particules
        D (int): nombre de dimensions

    Returns:
        dLJP (array): Gradient du potentiel des particules
    """
    dLJP = np.empty((nb_part,D), dtype=float)
    for i in range(nb_part):
        dLJP[i] = dLJpot_i(r0, i, sig, eps, cutoff, nb_part, D)
    return dLJP






#################################################################################################
#####################TEST DU TEMPS POUR L'ALGO VECTORIEL DE DLJP ET L'ANCIEN#####################
#################################################################################################


# Constantes
L_box = 30  #bord boite en nm
D = 2 #dimension
nb_part = 2  #nombre de particules
dt = 0.00002  #pas de temps en ps
T = 1000 #température initiale en Kelvin
m_part = 39.95  #masse particules en ua
nb_pas = 1000000

# Paramètres du potentiel Lennard-Jones
sig = 0.34 #paramètres de distance du potentiel en nm
Kb = 1.38e-23
Kb_scaled = 0.007831  #Kb en ua.nm^2.ps^-2.T
eps = 119.8*Kb_scaled   #paramètre décrivant la profondeur du puit de Lennard-Jones, énergie
cutoff = 3.2*sig

# Initialisation des positions et des vitesses
r, nb_part = pos_cristal2D(100, L_box)
v = vit_temp(nb_part, T, Kb_scaled, m_part)
r0 = r
r1 = r
v0 = v
v1 = v


@njit
def dLJP_test(r0, sig, eps, cutoff, nb_part, D):

    dLJP = np.empty((nb_part,D), dtype=float)
    cutoff2 = cutoff**2
    
    for j in range(nb_part):
        dist_xy = r0 - r0[j]
        dist2 = dist_xy.T[0]**2 + dist_xy.T[1]**2
        mask = ((dist2 <= cutoff2) & (dist2 > 0))  # Mask for particles within cutoff, excluding self
        valid_indices = np.where(mask)[0]
        dpart_cut = dist_xy[valid_indices]  # Assign valid distances
        dpart_calc_cut = dist2[valid_indices]  # Assign valid calculated distances
        dLJP[j] = -24.0*eps*np.dot(dpart_cut.T, 2.0*(sig**12)*(1.0/dpart_calc_cut)**7-(sig**6)*(1.0/dpart_calc_cut)**4)
    
    return dLJP


@njit
def cutoff_test(r0, cutoff, nb_part, D):
    dpart_cut = np.empty((nb_part,nb_part-1,D), dtype=float)
    dpart_calc_cut = np.empty((nb_part,nb_part-1,D), dtype=float)
    
    for i in range(nb_part):
        incr = 0
        dpart_cut_i = np.empty((nb_part-1,D), dtype=float)
        dpart_calc_cut_i = np.empty((nb_part-1,1), dtype=float)
        
        for j in range(nb_part):
            if j != i:
                dist_xy = r0[j] - r0[i]
                dist = np.sqrt(dist_xy[0]**2 + dist_xy[1]**2)
                if dist <= cutoff:
                    dpart_cut_i[incr] = dist_xy   #on ajoute cet élément dans l'array
                    dpart_calc_cut_i[incr] = dist
                    incr +=1
        dpart_cut[i] = dpart_cut_i
        dpart_calc_cut[i] = dpart_calc_cut_i
    return dpart_cut, dpart_calc_cut


@njit
def calc_dLJP_test(dpart_cut, dpart_calc_cut, sig, eps):
    dLJP = np.empty((nb_part,D), dtype=float)
    for i in range(nb_part):
        dLJP[i] = -24.0*eps*np.sum(dpart_cut[i]*(2.0*(sig**12)*(1.0/np.power(dpart_calc_cut[i], 14))-(sig**6)*(1.0/np.power(dpart_calc_cut[i], 8))), axis=0)
    return dLJP




for i in range(5):
    
    t = time.time()
    force_LJ0 = dLJP_test(r0, sig, eps, cutoff, nb_part, D)
    timeDLJP = time.time()-t
    print(f'temps new {i}: {timeDLJP}')
    force_wall0, force_wall_tot0 = LJ_walls(r, nb_part, sig, eps, L_box, D)
    force0 = force_LJ0+force_wall0
    r0, v0 = verlet(r0, v0, force0, dt, m_part)
    
    
    t = time.time()
    force_LJ1 = dLJpot(r1, sig, eps, cutoff, nb_part, D)
    timeDLJP = time.time()-t
    print(f'temps dLJP old {i}: {timeDLJP}')
    force_wall1, force_wall_tot0 = LJ_walls(r, nb_part, sig, eps, L_box, D)
    force1 = force_LJ1+force_wall1
    r1, v1 = verlet(r0, v1, force1, dt, m_part)
    
print(r0[:5])
print(r1[:5])