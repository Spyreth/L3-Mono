import numpy as np
from numba import njit


@njit
def LJ_walls(r, nb_part, sig, eps, L, D):
    """Calcul les forces associées pour des murs avec un potentiel LJ sur toutes les particules

    Args:
        r (array): Positions des particules
        sig (int): para LJ
        eps (int): para LJ
        L (int): Taille boîte
        D (int): Nb dimensions

    Returns:
        f (array): forces sur les particules
    """

    cut = 2**(1/6)*sig
    cut2 = L-cut
    force = np.empty((nb_part, D), dtype=float)
    incr = 0

    for elem in r:
        for d in range(D):
            if elem[d] < cut:
                force[incr][d] = -24*eps*((sig**6)*(1/(elem[d]**7)) - 2*(sig**12)*(1/elem[d]**13)) + eps
            elif elem[d] > cut2:
                force[incr][d] = 24*eps*((sig**6)*(1/((L-elem[d])**7)) - 2*(sig**12)*((1/(L-elem[d])**13))) + eps
            else:
                force[incr][d] = 0
        incr += 1

    return force


@njit
def reflectBC(r_0, v_0, nb_part, m_part, L_box, D, rayon):
    """
    Applique les conditions limites d'une boîte solide (rebond type billard)

    Args:
        r_0 (array): positions des particules à t0
        r_1 (array): positions des particules à t1

    Returns:
        r_0 (array): positions des particules à t0 corrigée
        r_1 (array): positions des particules à t1 corrigée
    """
    
    delta_p = 0
    r0 = r_0
    v0 = v_0
    
    for i in range(nb_part):  #pour chaque particule
        for j in range(D):   #dans chaque dimension
            #si r1 sort de la boite en 0 ou en L, on inverse les positions
            if r0[i][j]<rayon:
                r0[i][j] = -r0[i][j]+2*rayon
                v0[i][j] = -v0[i][j]
                delta_p += 2*m_part*v0[i][j]
            if r0[i][j]>(L_box-rayon):
                r0[i][j] = 2.0*L_box-r0[i][j]-2*rayon
                v0[i][j] = -v0[i][j]
                delta_p += -2*m_part*v0[i][j]
    return r0, v0, delta_p


@njit
def testPot(r, nb_part, D, sig, L, a):
    """
    Force linéaire test pour les murs, dérivant d'un potentiel quadratique

    Args:
        r (array): position des particules
        nb_part (int): nombre de particules
        D (int): nombre de dimensions
        sig (int): paramètre du puit de LJ
        L (int): taille de la boîte
        a (int): paramètre de la force linéaire en a*x

    Returns:
        array: force appliquée sur les particules dûe aux murs
    """

    force = np.empty((nb_part, D), dtype=float)
    incr = 0
    cut = sig
    cut2 = L - sig
    for elem in r:
        for d in range(D):
            if elem[d] < cut:
                force[incr][d] = -a*elem[d] + a*sig
            elif elem[d] > cut2:
                force[incr][d] = a*(L-elem[d]) - a*sig
            else:
                force[incr][d] = 0
        incr += 1
    return force



if __name__ == "__main__":
    print("Ceci n'est pas un script mais un package.")