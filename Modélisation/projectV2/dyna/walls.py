import numpy as np
from numba import njit


@njit
def LJ_walls(r, nb_part, sig, eps, L, D):
    """Calcul les forces associées pour des murs avec un potentiel LJ sur toutes les particules
    Args:
        r (array): positions des particules
        sig (float): paramètre Lennard-Jones
        eps (float): profondeur du puit de LJ
        L (float): Taille de la boîte
        D (int): Nb dimensions
    Returns:
        force, force_tot (array, int): forces sur les particules, force totale sur le mur pour le calcul de pression
    """

    cut = 2**(1/6)*sig
    cut2 = L-cut
    force = np.empty((nb_part, D), dtype=float)
    incr = 0

    #dans chaque dimension et pour chaque particule, on vérifie si on ne dépasse pas le cut et on applique le potentiel LJ corrigé
    for elem in r:
        for d in range(D):
            if elem[d] < cut:
                force[incr][d] = -24*eps*((sig**6)*(1/(elem[d]**7)) - 2*(sig**12)*(1/elem[d]**13)) + eps
            elif elem[d] > cut2:
                force[incr][d] = 24*eps*((sig**6)*(1/((L-elem[d])**7)) - 2*(sig**12)*((1/(L-elem[d])**13))) + eps
            else:
                force[incr][d] = 0
        incr += 1

    force_tot = np.sum(np.abs(force))

    return force, force_tot


@njit
def reflectBC(r_0, v_0, nb_part, m_part, L_box, D, rayon):
    """Applique les conditions limites d'une boîte solide (rebond type billard)
    Args:
        r_0 (array): positions des particules à t0
        r_1 (array): positions des particules à t1
    Returns:
        r_0 (array): positions des particules à t0 corrigée
        r_1 (array): positions des particules à t1 corrigée
        delta_p (float): différence de qté de mouvement des particules
    """
    
    delta_p = 0
    r0 = r_0
    v0 = v_0
    
    for i in range(nb_part):  #pour chaque particule
        for j in range(D):   #dans chaque dimension
            #si r1 sort de la boite en 0 ou en L, on change les positions et les vitesses dans la dimension associée
            if r0[i][j]<rayon:
                r0[i][j] = -r0[i][j]+2*rayon
                v0[i][j] = -v0[i][j]
                delta_p += 2*m_part*v0[i][j]  #pour le calcul de pression
            if r0[i][j]>(L_box-rayon):
                r0[i][j] = 2.0*L_box-r0[i][j]-2*rayon
                v0[i][j] = -v0[i][j]
                delta_p += -2*m_part*v0[i][j]  #pour le calcul de pression
    return r0, v0, delta_p



if __name__ == "__main__":
    print("Ceci n'est pas un script mais un package.")