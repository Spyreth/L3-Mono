import numpy as np
from numba import njit
from dyna.walls import reflectBC, LJ_walls



@njit
def verlet_LJ(r0, v0, force, pas, m, nb_part, sig, eps, cutoff, D, L_box) :
    """Calcule les positions et vitesses à t+1 via verlet

    Args:
        r0 (array): positions des particules à t0
        v0 (array): vitesses des particules à t0
        force (array): forces s'appliquant sur les particules à t0
        pas (float): pas de temps de calcul
        m (float): masse des particules

    Returns:
        r1 (array): positions des particules à t1
        v1 (array): vitesses des particules à t1
    """
    
    v1_2 = v0 + force*pas/(2*m)
    r1 = r0 + v1_2*pas
    force_LJ = dLJP(r0, sig, eps, cutoff, nb_part, D)
    force_wall, force_wall_tot = LJ_walls(r0, nb_part, sig, eps, L_box, D)
    force1 = force_LJ+force_wall
    v1 = v1_2 + force1*pas/(2*m)

    return r1, v1, force1, force_wall_tot


@njit
def verlet_billard(r0, v0, pas, nb_part, m, L_box, D, rayon) :
    r1 = r0 + v0*pas
    v1 = v0
    r1, v1, delta_p = reflectBC(r1, v1, nb_part, m, L_box, D, rayon)
    return r1, v1, delta_p



@njit
def dLJP(r, sig, eps, cutoff, nb_part, D):
    """Calcule les forces dûes aux interactions de Lennard-Jones
    Args:
        r (array): positions des particules
        sig (float): paramètre du puit de LJ
        eps (float): profondeur du puit de LJ
        cutoff (float): distance a partir de laquelle la force est nulle
        nb_part (int): nb de particules
        D (int): nombre de dimensions     
    Returns:
        (array): forces dûes à LJ s'appliquant sur les particules
    """

    dLJP = np.empty((nb_part,D), dtype=float)
    cutoff2 = cutoff**2
    
    for j in range(nb_part):
        
        #calcule des distances relatives à la particule j
        dist_xy = r - r[j]
        dist2 = dist_xy.T[0]**2 + dist_xy.T[1]**2 #norme au carrée
        
        #on garde en memoire les particules pour lesquelles dist<cutoff
        mask = ((dist2 <= cutoff2) & (dist2 > 0))
        valid_indices = np.where(mask)[0]
        dpart_cut = dist_xy[valid_indices]
        dpart_calc_cut = dist2[valid_indices]
        
        #calcule des forces s'appliquant sur la particule j dûes à LJ
        dLJP[j] = -24.0*eps*np.dot(dpart_cut.T, 2.0*(sig**12)*(1.0/dpart_calc_cut)**7-(sig**6)*(1.0/dpart_calc_cut)**4)
    
    return dLJP



if __name__ == "__main__":
    print("Ceci n'est pas un script mais un package.")