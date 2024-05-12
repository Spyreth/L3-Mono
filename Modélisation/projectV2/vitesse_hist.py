import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.animation import FuncAnimation
from filemanager.read import get_param, get_posvittime
from mesures.Measure import sumEC, sumLJpotsyst, sumLJwalls, calcTemp


save_folder = os.path.dirname(os.path.abspath(__file__)) + r'/Resultats'
results_name = r'/Vitesses_Maxwell'


param = get_param(save_folder+results_name+r'\param.txt')
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
kb = float(param['Kb'])


r, v, t = get_posvittime((save_folder+results_name), D, nb_part, nb_pas, save_interval)

#V est un tableau avec comme colone vx et vy, et comme lignes t, nous allons pour t0 et tfin calculer les normes de v pour chaque particule

length = np.size(r, axis=0)

def normvanim(v, num):
    result=np.empty((length), np.float64)
    for i in range(length):
         result[i] = np.sqrt(v[num][i][0]**2+v[num][i][1]**2)
    return result



fig, ax = plt.subplots(figsize=(15,15))
ax.set_title('Evolution de la distribution de vitesse')
ax.set_xlabel('Vitesse (nm/ps)')
ax.set_ylabel('Nombre de particules')

def animation(i):
    ax.clear()
    ax.hist(normvanim(v, i),range=(0,max(normvanim(v, i))), bins =100)

ani = FuncAnimation(fig,animation,length,interval=1)
ani.save(save_folder+results_name+'\gifhisto.gif',writer='pillow')
