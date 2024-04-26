import os
from initialisation.positions import random_pos, pos_cristal2D
from initialisation.vitesses import vit_temp, random_vit
from dyna.dynam import update
from filemanager.write import csv_init, save_parameters, datasave, pressureSaveBillard



# Constantes
L_box = 10  #bord boite en nm
D = 2 #dimension
nb_part = 2  #nombre de particules
dt = 0.0001  #pas de temps en ps
T = 300 #température en Kelvin
m_part = 39.95  #masse particules en ua
nb_pas = 200000

# Paramètres du potentiel Lennard-Jones
sig = 0.34 #paramètres de distance du potentiel en nm
Kb = 1.38e-23
Kb_scaled = 0.007831  #Kb en ua.nm^2.ps^-2.T
eps = 119.8*Kb_scaled   #paramètre décrivant la profondeur du puit de Lennard-Jones, énergie
cutoff = 3.2*sig

# Paramètres de l'animation
rayon = 0.1
save_interval = 1000
script_directory = os.path.dirname(os.path.abspath(__file__))
save_folder = os.path.dirname(os.path.abspath(__file__)) + r'\Resultats'
results_name = r'\testscale3'

# Initialisation des positions et des vitesses
r, nb_part = pos_cristal2D(5, L_box)
v = vit_temp(nb_part, T, Kb_scaled, m_part)

# Initialisation des fichiers de sauvegarde
csv_init(save_folder, results_name, 1, D)
save_parameters(save_folder, results_name, L_box, D, nb_part, dt, m_part, nb_pas, sig, eps, cutoff, rayon, save_interval, 0, Kb_scaled)


progress_affichage = nb_pas/100 #pour afficher le progrès tous les %

for i in range(nb_pas):
    if i % save_interval == 0:
        datasave(save_folder, results_name, r, v, i*dt, D)


    if i % progress_affichage == 0:
        progress = round(i / nb_pas * 100)
        print(f'\rAvancement calculs: {progress}%')

    r, v = update(r, v, dt, m_part, nb_part, sig, eps, cutoff, D, L_box)
print(f'\rAvancement calculs: Fin')
