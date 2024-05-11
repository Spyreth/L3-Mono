import os
from initialisation.positions import random_pos, pos_cristal2D
from initialisation.vitesses import vit_temp, random_vit
from dyna.dynam import update
from filemanager.write import csv_init, save_parameters, datasave, pressureSave



# Constantes
L_box = 15  #bord boite en nm
D = 2 #dimension
dt = 0.00001  #pas de temps en ps
T = [i*30 for i in range(30)] #température initiale en Kelvin
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
save_folder = os.path.dirname(os.path.abspath(__file__)) + r'/Resultats/GP2_12x12_L15'
results_name = r'/GP2_12x12_L15'

for temp in T:
    # Initialisation des positions et des vitesses
    r, nb_part = pos_cristal2D(12, L_box)
    v = vit_temp(nb_part, temp, Kb_scaled, m_part)

    # Initialisation des fichiers de sauvegarde
    csv_init(save_folder, results_name+f'_Ti={temp}', 1, D)
    save_parameters(save_folder, results_name+f'_Ti={temp}', L_box, D, nb_part, dt, m_part, nb_pas, sig, eps, cutoff, rayon, save_interval, pressure_calc_interval, Kb_scaled)


    progress_affichage = nb_pas/100 #pour afficher le progrès tous les %
    f_sum = 0

    for i in range(nb_pas):
        
        if i % save_interval == 0:
            datasave(save_folder, results_name+f'_Ti={temp}', r, v, i*dt, D)

        if i % progress_affichage == 0:
            progress = round(i / nb_pas * 100)
            print(f'Avancement calculs pour T={temp}: {progress}%')

        r, v, force_wall_tot = update(r, v, dt, m_part, nb_part, sig, eps, cutoff, D, L_box)
        f_sum += force_wall_tot
    
        if i % pressure_calc_interval == 0:
            pressure = f_sum/(pressure_calc_interval*4*L_box)
            pressureSave(save_folder, results_name+f'_Ti={temp}', pressure)
            f_sum = 0

print(f'Avancement calculs: Fin')
