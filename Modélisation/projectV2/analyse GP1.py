import numpy as np
import matplotlib.pyplot as plt
import os
from filemanager.read import get_param, get_posvittime, get_pressure
from mesures.Measure import calcTemp
from scipy.optimize import curve_fit


T = [i*30 for i in range(20)] #températures initiales en Kelvin
save_folder = os.path.dirname(os.path.abspath(__file__)) + r'/Resultats/GP2_12x12_L15'
results_name = r'/GP2_12x12_L15'

param = get_param(save_folder+results_name+f'_Ti={T[0]}'+r'/param.txt')
# Utilisation des paramètres récupérés
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

    #Récupération des positions, des vitesses et des pressions
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

#Modèle ax + b : équation d'état de van der waals
def vdwpressure(x, a, b):
    return a*x+b
popt2, pcov2 = curve_fit(vdwpressure, T_calc[1:], p_calc[1:])
slope2 = popt2[0]
intercept2 = popt2[1]
slope_err2 = np.sqrt(np.diag(pcov2))[0]
intercept_err2 = np.sqrt(np.diag(pcov2))[1]

#Pour l'affichage
T_plot = np.linspace(T[0], T[-1], 200)

#Pente théorique modèle ax: nR/V
slope_theorique1 = np.float64(nb_part*kb/(L_box**2)).round(5)

#Pente théorique modèle ax+b : nR/(V-nb)
n_moles = nb_part/(6.022*(10**23))
b_theorique = 0.03219*(10**24)  #valeur trouvée pour l'argon dans les unités du projet
slope_theorique2 = np.float64(nb_part*kb/(L_box**2 - (n_moles*b_theorique))).round(5)
#Intercepts théorique modèle ax+b : -an^2/(V^2)
a_theorique = 0.1363*6.022*(10**47)  #valeur trouvée pour l'argon dans les unités du projet
intercept_theorique2 = np.float64(-a_theorique*(n_moles**2)/(L_box**4)).round(5)


#Graphe de la régression modèle ax
plt.figure(figsize=(12,8))
plt.scatter(T_calc, p_calc, label='p en fonction de T')
plt.plot(T_plot, T_plot*slope, 'r-', label='ajustement linéaire')
plt.text(400, 0.7, f'pente = {slope.round(5)} +- {slope_err.round(5)}', fontsize=12)
plt.text(400, 0.5, f'Nkb/V = {slope_theorique1}', fontsize=12)
plt.xlabel('T (K)')
plt.ylabel('p (unités du projet)')
plt.title('Vérification de la loi des gaz parfaits')
plt.grid()
plt.legend()
plt.savefig(save_folder+r'/Loi des gaz parfaits.png')


#Graphe de la régression modèle ax+b
plt.figure(figsize=(12,8))
plt.scatter(T_calc[1:], p_calc[1:], label='p en fonction de T')
plt.plot(T_plot, T_plot*slope2+intercept2, 'r-', label='ajustement affine')
plt.text(400, 0.9, f'pente = {slope2.round(5)} +- {slope_err2.round(5)}', fontsize=12)
plt.text(400, 0.7, f'intercept = {intercept2.round(5)} +- {intercept_err2.round(5)}', fontsize=12)
plt.text(400, 0.5, f'Nkb/(V-nb) = {slope_theorique2}', fontsize=12)
plt.text(400, 0.3, f'-an^2/V^2 = {intercept_theorique2}', fontsize=12)
plt.xlabel('T (K)')
plt.ylabel('p (unités du projet)')
plt.title('Equation d\'état de Van der Waals')
plt.grid()
plt.legend()
plt.savefig(save_folder+r'/Loi de VdW.png')
    