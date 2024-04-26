import numpy as np
from numba import njit

@njit
def random_vit(n, L_box, D):
    """Initialise les vitesses aléatoires de particules

    Args:
        n (int): nb de particules
        L_box (int): taille de la boite
        D (int): dimension

    Returns:
        v_0 (D-d array): vitesses nes n particules à D-dim
    """
    v_0 = L_box/10*(np.random.rand(n,D)-0.5)
    return v_0


@njit
def vit_temp(n, T, kb, m):
    """Initialise les vitesses des particules en fonction de la température

    Args:
        n (int): nb de particules
        T (int): température initiale du système
        kb (int): constante de Boltzmann dans les bonnes unités
        m (int): masse des particules
        D (int): dimension

    Returns:
        v_0 (D-d array): vitesses nes n particules à D-dim
    """
    angles = np.random.uniform(0, 2*np.pi, n)
    x = np.cos(angles)
    y = np.sin(angles)
    v_0 = np.column_stack((x,y))*np.sqrt(2*kb*T/m)
    return v_0



if __name__ == "__main__":
    L_box = 100
    D = 2
    n = 10
    vitesses_initiales = random_vit(n, L_box, D)
    print(vitesses_initiales)