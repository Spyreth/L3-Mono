import numpy as np
import matplotlib.pyplot as plt
import os
from filemanager.read import get_param, get_posvittime, get_pressure
from mesures.Measure import calcTemp
from scipy.optimize import curve_fit


T = [i*30 for i in range(30)] #températures initiales en Kelvin
save_folder = os.path.dirname(os.path.abspath(__file__)) + r'/Resultatsbillard/GP2_12x12_L15'
results_name = r'/GP2_12x12_L15'

param = get_param(save_folder+results_name+f'_Ti={T[0]}'+r'/param.txt')
#Utilisation des paramètres récupérés
L_box = float(param["L"])
D = int(param["D"])
nb_part = int(param["nb_part"])
dt = float(param["dt"])
m_part = float(param["m_part"])
nb_pas = int(param["nb_pas"])
sig = float(param["sig"])
eps = float(param["eps"])
cutoff = float(param["cutoff"])
rayon = float(param["rayon"])
save_interval = int(param["save_interval"])
pressure_calc_interval = int(param["pressure_calc_interval"])
kb = float(param['Kb'])

#Initialisation des variables de calcul de moyennes de T et p
T_calc = np.empty(len(T))
p_calc = np.empty(len(T))
i = 0


#Boucles sur chaque simulation faite avec simu GP1 billard
for temp in T:

    #Récupération des positions, des vitesses, des temps et des pressions
    r, v, t = get_posvittime((save_folder+results_name+f'_Ti={temp}'), D, nb_part, nb_pas, save_interval)
    pressure = get_pressure((save_folder+results_name+f'_Ti={temp}'))

    #Calcul des températures sur la simulation
    len = np.size(r, axis=0)
    T_i = np.empty((len), np.float64)
    for j in range(len):
        T_i[j] = calcTemp(v[j], m_part, kb)
    
    #Calcul des moyennes de T et p sur les 70 derniers pourcents de la simulation pour être sûrs d'être à l'équilibre thermo
    indice_T = int(0.3*np.size(T_i))
    indice_p = int(0.3*np.size(pressure))
    T_moy = np.mean(T_i[indice_T:])
    p_moy = np.mean(pressure[indice_p:])
    
    #Ajout des moyennes dans les array
    T_calc[i] = T_moy
    p_calc[i] = p_moy
    i += 1
    
#Modèle ax : loi des GP
def linear(x,a):
    return a*x
popt, pcov = curve_fit(linear, T_calc, p_calc)
slope = popt[0]
slope_err = np.sqrt(np.diag(pcov))[0]


#Pour l'affichage
T_plot = np.linspace(T[0], T[-1], 200)   

#Pente théorique : nR/V
slope_theorique = np.float64(nb_part*kb/(L_box**2)).round(6)

#Graphe de la régression
plt.figure(figsize=(12,8))
plt.scatter(T_calc, p_calc, label='Données simulées')
plt.plot(T_plot, T_plot*slope, 'r-', label='Ajustement linéaire')
plt.text(610, 1.6, f'pente = {slope.round(6)} +- {slope_err.round(6)}', fontsize=12)
plt.text(610, 1.3, f'Nkb/V = {slope_theorique}', fontsize=12)
plt.xlabel('T (K)')
plt.ylabel('p (unités du projet)')
plt.title('Vérification de la loi des gaz parfaits')
plt.grid()
plt.legend()
plt.savefig(save_folder+r'/Loi des gaz parfaits.png')
    