import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import animation
from filemanager.write import *
from filemanager.read import *
from mesures.Measure import *


save_folder = os.path.dirname(os.path.abspath(__file__)) + r'\Resultatsbillard'
results_name = r'\testpressure2'


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
pressure_calc_interval = int(param["pressure_calc_interval"])
kb = float(param['Kb'])


r, v, t = get_posvittime((save_folder+results_name), D, nb_part, nb_pas, save_interval)
pressure = get_pressure_billard(save_folder+results_name)

len = np.size(r, axis=0)
E_C = np.empty((len), np.float64)
T = np.empty((len), np.float64)


for i in range(len):
    E_C[i] = sumEC(v[i], m_part, nb_part)
    T[i] = calcTemp(v[i], m_part, kb)

E_tot = E_C
pV = pressure*(L_box**2)
NkbT = nb_part*kb*T

plt.figure(figsize=(12,8))
plt.plot(t, E_tot, 'r-', label='E_tot')
plt.xlabel('t (ps)')
plt.ylabel('E')
plt.legend()
plt.savefig(save_folder+results_name+r'\total energy.png')

plt.figure(figsize=(12,8))
plt.plot(t, T, 'r-', label='Température')
plt.xlabel('t (ps)')
plt.ylabel('T (K)')
plt.legend()
plt.savefig(save_folder+results_name+r'\Temperature.png')

plt.figure(figsize=(12,8))
plt.plot((np.linspace(0, nb_pas*dt, int(nb_pas/pressure_calc_interval))), pressure, 'r-', label='Pression')
plt.xlabel('t (ps)')
plt.ylabel('p')
plt.legend()
plt.savefig(save_folder+results_name+r'\Pression.png')

plt.figure(figsize=(12,8))
plt.plot(t, NkbT, 'r-', label='NKbT')
plt.plot((np.linspace(0, nb_pas*dt, int(nb_pas/pressure_calc_interval))), pV, 'b-', label='pV')
plt.xlabel('t (ps)')
plt.ylabel('pV, NKbT')
plt.legend()
plt.savefig(save_folder+results_name+r'\Loi des gaz parfaits.png')










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
anim.save(save_folder + results_name + '\gaz.gif', writer=writergif)
#writergif = animation.FFMpegWriter(fps=20)
#anim.save('gaz.mp4', fps=20, dpi=200)
