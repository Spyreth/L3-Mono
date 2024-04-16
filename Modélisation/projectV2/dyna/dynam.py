import numpy as np
from numba import njit
from dyna.walls import *



@njit
def verlet(r0, v0, force, pas, m) :
    """Calcule les positions et vitesses à t+1 via verlet

    Args:
        r0 (array): positions des particules à t0
        v0 (array): vitesses des particules à t0
        force (array): forces s'appliquant sur les particules à t0
        pas (int): pas de temps de calcul
        m (int): masse des particules

    Returns:
        r1 (array): positions des particules à t1
        v1 (array): vitesses des particules à t1
    """
    
    v1_2 = v0 + force*pas/(2*m)
    r1 = r0 + v1_2*pas
    v1 = v1_2 + force*pas/(2*m)

    return r1, v1




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




@njit
def update(r_0, v_0, pas, m, nb_part, sig, eps, cutoff, D, L_box) :
    """
    Calcule les positions et vitesses à t+1 à partir de t en utilisant verlet

    Args:
        r_0 (array): positions des particules à t0
        v_0 (array): vitesses des particules à t0
        pas (int): pas de temps d'intégration
        m (int): masse des particules
        sig (int): paramètre du puit de LJ
        eps (float): profondeur du puit de LJ
        cutoff (int): distance a partir de laquelle la force est nulle
        D (int): nombre de dimensions
        L_box (int): taille de la boîte

    Returns:
        r_1 (array): positions des particules à t1
        v_1 (array): vitesses des particules à t1
    """
    force_LJ = dLJpot(r_0, sig, eps, cutoff, nb_part, D)
    force_wall = LJ_walls(r_0, nb_part, sig, eps, L_box, D)
    #force_wall = testPot(r_0, nb_part, D, sig, L_box, 10000)
    #force_wall = 0

    force = force_LJ+force_wall

    r_1, v_1 = verlet(r_0, v_0, force, pas, m)
    return r_1, v_1



if __name__ == "__main__":
    print("Ceci n'est pas un script mais un package.")