import numpy as np
from numba import njit

@njit
def random_pos(n, L_box, D):
    """Initialise les positions de n particules aléatoires dans la boîte

    Args:
        n (int): nb de particules
        L_box (int): taille de la boite
        D (int): dimension

    Returns:
        r_0 (D-d array): positions des n particules à D-dim
    """
    r_0 = L_box*np.random.rand(n,D)
    return r_0


@njit
def pos_cristal2D(nb_particules, len_box):
    """
    Initialise les positions des particules dans la boîte en forme de cristal cubique

    Args:
        nb_particules (float): nb de particules par arête d'un axe du cristal cubique
        len_box (float): longueur de la boîte (un axe)

    Returns:
        r0 (array): Array de positions initiales
        n_particules: nb de particules dans le cristal
    """
    
    r0 = np.empty((nb_particules**2, 2), dtype=float)
    a = len_box/(nb_particules+1) #distance interparticulaire
    
    for i in range(nb_particules):
        for j in range(nb_particules):
            x = (i+1)*a
            y = (j+1)*a
            r0[i*nb_particules+j*1] = [x,y]
    n_particules = len(r0)

    return r0, n_particules

if __name__ == "__main__":
    positions_initiales, nb_part = pos_cristal2D(4,20)
    print(positions_initiales)
    print(nb_part)