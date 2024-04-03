import numpy as np
from numba import njit


@njit
def update(r_0, v_0, pas, m, nb_part, sig, eps, cutoff, D, L_box) :

    force = np.empty((nb_part,D), dtype=float)
    force_wall = LJ_walls(r_0, nb_part, sig, eps, L_box, D)
    #force_wall = 0
    for i in range(nb_part):
        force[i] = dLJpotnumb(r_0, i, sig, eps, cutoff, nb_part, D)

    force = force+force_wall

    r_1, v_1 = verlet(r_0, v_0, force, pas, m)
    
    return r_1, v_1


@njit
def LJ_walls(r, nb_part, sig, eps, L, D):
    """Calcul les forces associées aux murs

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
                force[incr][d] = -24*eps*((sig**6)*(1/(elem[d]**7)) - 2*(sig**12)*(1/elem[d]**13))
            elif elem[d] > cut2:
                force[incr][d] = 24*eps*((sig**6)*(1/((L-elem[d])**7)) - 2*(sig**12)*((1/(L-elem[d])**13)))
            else:
                force[incr][d] = 0
        incr += 1
    return force





@njit
def verlet(r0, v0, force, pas, m) :
    
    v1_2 = v0 + force*pas/(2*m)
    r1 = r0 + v1_2*pas
    v1 = v1_2 + force*pas/(2*m)

    return r1, v1


@njit
def dLJpot(r, i, sig, eps, cut_off, nb_part, D):
    
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
def dLJpotnumb(r, i, sig, eps, cut_off, nb_part, D):
    
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
    r814v = (dpart_cut)*r814
    dLJP = 24.0*eps*np.sum(r814v,axis=0)
    
    return dLJP



if __name__ == "__main__":

    #####TESTS#####

    L = 20
    D = 2
    nb_arete = 2
    sig = 3.4
    eps = 120
    cut_off = 3*sig

    r0 = np.empty((nb_arete**2, 2), dtype=float)
    a = L/(nb_arete+1) #distance interparticulaire
    
    for i in range(nb_arete):
        for j in range(nb_arete):
            x = (i+1)*a
            y = (j+1)*a
            r0[i*nb_arete+j*1] = [x,y]
    n_particules = len(r0)
    v_0 = 40/20*(np.random.rand(n_particules,D)-0.5)
    print(r0)
    print(v_0)