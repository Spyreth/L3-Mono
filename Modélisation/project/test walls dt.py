import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import collections
from numba import njit
from initialisation.positions import *
from initialisation.vitesses import *
from dyna.dynam import *
from mesures.Measure import *

# Constantes
L_box = 10  #bord boite en Angstrom
D = 2 #dimension
nb_part = 1  #nombre de particules
m_part = 20  #masse particules en ua

# Paramètres du potentiel Lennard-Jones
sig = 3.4 #paramètres de distance du potentiel en Angstrom
Kb = 1.38e-23
eps = 120   #paramètre décrivant la profondeur du puit de Lennard-Jones, énergie
cutoff = 3.2*sig

# Paramètres de l'animation
rayon = 2
save_interval = 10



r0 = np.array([[5,5]])
nb_part = 1
v0 = np.array([[40,0]])
nb_test = 1000
EC_i = sumEC(v0, m_part, nb_part)
print(EC_i)

dt_pas = 0.0000001  #pas de temps en ps
tot_time = 1.5
dt_arr = np.array([dt_pas*(i+1) for i in range(nb_test)])
nb_pas_arr = np.array([int(tot_time/dt_arr[i]) for i in range(nb_test)])
diff_EC = np.empty((nb_test), np.float64)

for j in range(nb_test):
    
    dt = dt_arr[j]
    nb_pas = nb_pas_arr[j]
    r = r0
    v = v0

    E_LJ = np.empty((nb_pas), np.float64)
    E_LJwalls = np.empty((nb_pas), np.float64)
    E_C = np.empty((nb_pas), np.float64)
    t = np.empty((nb_pas), np.float64)
    t2 = np.empty((nb_pas), np.float64)

    for i in range(nb_pas):
        r, v = update(r, v, dt, m_part, nb_part, sig, eps, cutoff, D, L_box)
        E_LJ[i] = sumLJpot(r, sig, eps, cutoff, nb_part)
        E_LJwalls[i] = sumLJwalls(r, sig, eps, nb_part, L_box)
        E_C[i] = sumEC(v, m_part, nb_part)
        t2[i] = i*dt

    Esum = E_C
    Esum = Esum - Esum[0]
    diff = Esum[-1]-Esum[0]
    relative_diff = diff/EC_i
    

    if relative_diff < 1000 and relative_diff>-1000:
        diff_EC[j] = relative_diff
    else:
        diff_EC[j] = 0

    if j % 10 == 0:
        progress = round(j / nb_pas * 100, 2)
        print(f'\rAvancement calculs: {progress}%')
print(f'\rAvancement calculs: Fin')


plt.plot(dt_arr, diff_EC)
plt.show()