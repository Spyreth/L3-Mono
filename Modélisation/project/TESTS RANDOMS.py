import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import collections
from numba import njit
from initialisation.positions import *
from initialisation.vitesses import *
from dyna.dynam import *
from mesures.Measure import *



sig = 3.4
L = 20
eps = 120
cut_off = 3.2*sig

dLJP1 = np.empty((20,2), dtype=np.float64)
dLJP2 = np.empty((20,2), dtype=np.float64)
r1 = np.empty((20,1), dtype=np.float64)
r2 = np.empty((20,1), dtype=np.float64)

#r = random_pos(nb_part, L_box, D)
#r, nb_part = pos_cristal2D(2, L_box)
for i in range(20):

    ri = np.array([[10-((i+1) + 20)/10,10],[10+((i+1) + 20)/10,10]])
    nb_part = 2
    #v = random_vit(nb_part, L_box, D)
    v = np.array([[0,0], [0,0]])

    dLJP1[i] = dLJpotnumb(ri, 0, sig, eps, cut_off, nb_part, 2)
    dLJP2[i] = dLJpotnumb(ri, 1, sig, eps, cut_off, nb_part, 2)
    r1[i] = ri[0,0]
    r2[i] = ri[1,0]


print(dLJP1)
print(dLJP2)
plt.figure(figsize=(10,8))
plt.plot(r1, dLJP1.T[0], 'b-', label='particule 1 (gauche)')
plt.plot(r1, dLJP2.T[0], 'r-', label='particule 2 (droite)')
plt.xlabel('r')
plt.ylabel('force')
plt.legend()
plt.show()
