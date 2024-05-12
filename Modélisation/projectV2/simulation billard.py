import os
from initialisation.positions import random_pos, pos_cristal2D
from initialisation.vitesses import random_vit, vit_temp
from dyna.dynam import verlet_billard
from filemanager.write import csv_init, save_parameters, datasave, pressureSave



# Constantes
L_box = 15  #bord boite en nm
D = 2 #dimension
nb_part = 2  #nombre de particules, remplacé dans la suite si on initialise les positions avec pos_cristal2D
dt = 0.00001  #pas de temps en ps
T = 300 #température initiale en Kelvin
m_part = 39.95  #masse particules en ua
nb_pas = 10_000_000

# Paramètres du potentiel Lennard-Jones
sig = 0.34 #paramètres de distance du potentiel en nm
Kb = 1.38e-23
Kb_scaled = 0.007831  #Kb en ua.nm^2.ps^-2.T
eps = 119.8*Kb_scaled   #paramètre décrivant la profondeur du puit de Lennard-Jones, énergie
cutoff = 3.2*sig

# Paramètres de l'animation et des mesures
rayon = 0.1
save_interval = 100_000
pressure_calc_interval = 1_000_000
script_directory = os.path.dirname(os.path.abspath(__file__))
save_folder = os.path.dirname(os.path.abspath(__file__)) + r'/Resultatsbillard'
results_name = r'/testGP'

# Initialisation des positions et des vitesses
r, nb_part = pos_cristal2D(12, L_box)
v = vit_temp(nb_part, T, Kb_scaled, m_part)

# Initialisation des fichiers de sauvegarde
csv_init(save_folder, results_name)
save_parameters(save_folder, results_name, L_box, D, nb_part, dt, m_part, nb_pas, sig, eps, cutoff, rayon, save_interval, pressure_calc_interval, Kb_scaled)


progress_affichage = nb_pas/100 #pour afficher le progrès tous les %
delta_p_tot = 0 #pour calculer la pression


for i in range(nb_pas):
    
    # Sauvegarde des données tous les save_interval pas de temps
    if i % save_interval == 0:
        datasave(save_folder, results_name, r, v, i*dt, D)

    # Affichage de l'avancement
    if i % progress_affichage == 0:
        progress = round(i / nb_pas * 100)
        print(f'Avancement calculs: {progress}%')

    # Application de verlet
    r, v, delta_p = verlet_billard(r, v, dt, nb_part, m_part, L_box, D, rayon)
    delta_p_tot += delta_p  #ajout de la qté de mvt échangée avec les murs sur ce pas de temps pour le calcul de pression

    if i % pressure_calc_interval == 0:
        pressure = delta_p_tot/(pressure_calc_interval*dt*4*L_box)
        pressureSave(save_folder, results_name, pressure)
        delta_p_tot = 0

print(f'Avancement calculs: Fin')
