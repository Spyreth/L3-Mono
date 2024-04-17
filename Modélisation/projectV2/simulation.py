import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import animation
from numba import njit
from initialisation.positions import *
from initialisation.vitesses import *
from dyna.dynam import *
from dyna.walls import *
from filemanager.write import *

# Constantes
L_box = 20  #bord boite en Angstrom
D = 2 #dimension
nb_part = 1  #nombre de particules
dt = 0.0001  #pas de temps en ps
m_part = 20  #masse particules en ua
nb_pas = 200000

# Paramètres du potentiel Lennard-Jones
sig = 3.4 #paramètres de distance du potentiel en Angstrom
Kb = 1.38e-23
eps = 120   #paramètre décrivant la profondeur du puit de Lennard-Jones, énergie
cutoff = 3.2*sig

# Paramètres de l'animation
rayon = 1
save_interval = 1000
script_directory = os.path.dirname(os.path.abspath(__file__))
save_folder = os.path.dirname(os.path.abspath(__file__)) + r'\Resultats'
results_name = r'\test0'

# Initialisation des positions et des vitesses
r = np.array([[10,10]])
v = np.array([[0,10]])

# Initialisation des fichiers de sauvegarde
csv_init(save_folder, results_name, 1, D)
save_parameters(save_folder, results_name, L_box, D, nb_part, dt, m_part, nb_pas, sig, eps, cutoff, rayon, save_interval)


progress_affichage = nb_pas/100 #pour afficher le progrès tous les %

for i in range(nb_pas):
    if i % save_interval == 0:
        datasave(save_folder, results_name, r, v, i*dt, D)


    if i % progress_affichage == 0:
        progress = round(i / nb_pas * 100, 2)
        print(f'\rAvancement calculs: {progress}%')

    r, v = update(r, v, dt, m_part, nb_part, sig, eps, cutoff, D, L_box)
print(f'\rAvancement calculs: Fin')
