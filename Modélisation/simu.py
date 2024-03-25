import numpy as np
import matplotlib.pyplot as plt
from numba import njit

#Constants
L_box = 20  #bord boite
D = 2 #dimension
nb_part = 10  #nombre de particules
dt = 0.1  #pas de temps en ps
box = np.zeros([D])+L_box   #taille système (3D)
m_part = 1  #masse particules
nb_pas = 1000

# Paramètres du potentiel Lennard-Jones
sig = 0.34 #paramètres de distance du potentiel en Angstrom
Kb = 1.38e-23
eps = 120*Kb   #paramètre décrivant la profondeur du puit de Lennard-Jones, énergie
cutoff = 3.2*sig



def random_pos(n):
    """Initialise les positions de n particules aléatoires dans la boîte

    Args:
        n (int): nb de particules

    Returns:
        r_0 (D-d array): positions des n particules à D-dim
    """
    r_0 = L_box*np.random.rand(n,D)
    return r_0



def random_vit(n):
    """Initialise les vitesses aléatoires de particules

    Args:
        n (int): nb de particules

    Returns:
        v_0 (D-d array): vitesses nes n particules à D-dim
    """
    v_0 = L_box/20*(np.random.rand(n,D)-0.5)
    return v_0


@njit
def dLJpot(r, i, sig, eps, cut_off):
    
    """
    Gradient du potentiel de Lennard-Jones pour une particule

    Args:
        r (array): Positions des particules
        i (float): Index de la particule
        sig (float): Taille des particules
        eps (float): Profondeur du puit

    Returns:
        dLJP (array): Gradient du potentiel de la particule
    """
    
    dpart = r-r[i] #distance entre la i-eme particule et les autres
    dpart = np.delete(dpart,2*i)
    dpart = np.delete(dpart,2*i)  #on retire le i-eme element car pas d'interaction d'une particule avec elle-même
    dpart = dpart.reshape(nb_part-1,2)
    dpart_cut = np.empty((nb_part-1,D), dtype=float)
    dpart_calc_cut = np.empty((nb_part-1,1), dtype=float)

    incr = 0
    for elem in dpart :
        distance_calc = np.sqrt(elem[0]**2 + elem[1]**2)
    
        if distance_calc <= cut_off:
            dpart_cut[incr] = elem   #on ajoute cet élément dans l'array
            dpart_calc_cut[incr] = distance_calc
            incr +=1
    
    dpart_cut = dpart_cut[:incr]
    dpart_calc_cut = dpart_calc_cut[:incr]

    r8 =  (sig**6)*(1.0/dpart_calc_cut)**8
    r14 = 2.0*(sig**12)*(1.0/dpart_calc_cut)**14
    r814 = r14-r8
    r814v = (dpart_cut)*r814
    dLJP = 24.0*eps*np.sum(r814v,axis=0)
    
    return dLJP


@njit
def verlet(r0, v0, force, pas, m) :
    
    v1_2 = v0 + force*pas/(2*m)
    r1 = r0 + v1_2*pas/2
    v1 = v1_2 + force*pas/(2*m)

    return r1, v1



def update(r_0, v_0, pas, m) :
    
    force = -np.array([dLJpot(r0, j, sig, eps, cutoff) for j in range(nb_part)])
    r_1, v_1 = verlet(r_0, v_0, force, pas, m)
    
    
    return r_1, v_1


r0 = random_pos(nb_part)
v0 = random_vit(nb_part)
print("r0 =")
print(r0)

r1, v1 = update(r0, v0, dt, m_part)
print(r1)



