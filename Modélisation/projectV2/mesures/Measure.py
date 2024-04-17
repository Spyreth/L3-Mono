import numpy as np
from numba import njit


@njit
def LJpot(r, i, sig, eps, cut_off, nb_part, trace):
    """
    Calcul du potentiel de Lennard-Jones pour une particule

    Args:
        r (array): Positions des particules
        i (float): Index de la particule
        sig (float): Taille des particules
        eps (float): Paramètre de Lennard-Jones
        trace (array): array de bool gardant une trace des interactions pour ne pas les compter deux fois

    Returns:
        LJP (float): Energie potentielle de Lennard-Jones de la i-eme particule
    """
    
    incr = 0
    dpart_calc_cut = np.empty((nb_part-1,1), dtype=float)

    for j in range(nb_part):
        if j != i and (trace[i, j]==0):
            dist_xy = r[j] - r[i]
            dist = np.sqrt(dist_xy[0]**2 + dist_xy[1]**2)
            if dist <= cut_off:
                dpart_calc_cut[incr] = dist
                incr +=1
            trace[i,j] = 1

    dpart_calc_cut = dpart_calc_cut[:incr]
            
    #calcul du potentiel
    r6 = (sig/dpart_calc_cut)**6
    r12 = (sig/dpart_calc_cut)**12
    LJP = 4.0*eps*(r12-r6)
    LJP = np.sum(LJP)  #somme des potentiels
    
    return LJP



@njit
def LJpotWalls(r, sig, eps, nb_part, L):
    cut = 2**(1/6)*sig
    cut2 = L-cut
    E_walls = np.empty((nb_part), dtype=float)
    incr = 0
    for elem in r:
        E_LJi = np.empty(2, dtype=float)
        for d in range(2):
            if elem[d] < cut:
                E_LJi[d] = 4*eps*((sig/elem[d])**12-(sig/elem[d])**6) + eps
            elif elem[d] > cut2:
                E_LJi[d] = 4*eps*((sig**12)*((1/(L-elem[d])**12)) - (sig**6)*(1/((L-elem[d])**6))) + eps
            else:
                E_LJi[d] = 0
        E_walls[incr] = E_LJi[0] + E_LJi[1]
        incr += 1
    return E_walls


@njit
def sumLJpot(r, sig, eps, cut_off, nb_part):
    trace = np.zeros((nb_part,nb_part), dtype=np.int32)
    E_LJ = np.empty((nb_part, 1, 1), np.float64)
    for i in range(nb_part):
        E_LJ[i] = LJpot(r, i, sig, eps, cut_off, nb_part, trace)
    sum = np.sum(E_LJ)
    return sum


@njit
def sumEC(v, m, nb_part):
    vit_abs2 = np.empty((nb_part, 1, 1), np.float64)
    for i in range(nb_part):
        vit_abs2[i] = v[i][0]**2 + v[i][1]**2
    sum = m*np.sum(vit_abs2)/2
    return sum


@njit
def sumLJwalls(r, sig, eps, nb_part, L):
    E_walls = LJpotWalls(r, sig, eps, nb_part, L)
    E = np.sum(E_walls)
    return E




if __name__ == "__main__":
    print("Ceci n'est pas un script mais un package.")