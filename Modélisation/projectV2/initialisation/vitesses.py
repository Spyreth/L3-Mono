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
    v_0 = L_box/2*(np.random.rand(n,D)-0.5)
    return v_0


if __name__ == "__main__":
    L_box = 100
    D = 2
    n = 10
    vitesses_initiales = random_vit(n, L_box, D)
    print(vitesses_initiales)