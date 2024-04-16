import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from numba import njit
from initialisation.positions import *
from initialisation.vitesses import *
from dyna.dynam import *
from dyna.walls import *
from filemanager.write import *
from filemanager.read import *


save_folder = r'C:\GIT REPOS\L3-Mono\Modélisation\projectV2\Resultats'
results_name = r'\test0'


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

r, v, t = get_posvittime(r"C:\GIT REPOS\L3-Mono\Modélisation\projectV2\Resultats\test0", D, nb_part, nb_pas, save_interval)



#################################################
######## ANIMATION ###################################
def init():
    global r, rayon, nb_part, ax, particles
    
    for i in range(nb_part):
        circle = plt.Circle((-10, -10), rayon,facecolor='red', edgecolor='black', lw=1)
        ax.add_patch(circle)
        particles.append(circle)
    return particles 

def animate(i):
    global r, particles

    for j in range(nb_part):
        particles[j].center = (r[i,j, 0], r[i,j, 1])

    if i % 50 == 0:
        progress = round(i / nb_pas * save_interval * 100, 1)
        print(f'\rAvancement animation: {progress}%')

    return particles 



fig, ax = plt.subplots()
ax.set_xlim(0, L_box)
ax.set_ylim(0, L_box)
ax.set_aspect('equal')
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)  

particles = []

anim = animation.FuncAnimation(fig, animate, frames=int(nb_pas/save_interval), interval=10, blit=True, init_func=init, repeat=False)


writergif = animation.PillowWriter(fps=20)
anim.save('C:\GIT REPOS\L3-Mono\Modélisation\projectV2\Resultats' + results_name + '\gaz.gif', writer=writergif)
#writergif = animation.FFMpegWriter(fps=20)
#anim.save('gaz.mp4', fps=20, dpi=200)
