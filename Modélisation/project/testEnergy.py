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
L_box = 20  #bord boite en Angstrom
D = 2 #dimension
nb_part = 1  #nombre de particules
dt = 0.000001  #pas de temps en ps
m_part = 20  #masse particules en ua
nb_pas = 1000000

# Paramètres du potentiel Lennard-Jones
sig = 3.4 #paramètres de distance du potentiel en Angstrom
Kb = 1.38e-23
eps = 120   #paramètre décrivant la profondeur du puit de Lennard-Jones, énergie
cutoff = 3.2*sig

# Paramètres de l'animation
rayon = 1
save_interval = 10000




#r = random_pos(nb_part, L_box, D)
#r, nb_part = pos_cristal2D(2, L_box)
r = np.array([[10-2,10],[10+2,10]])
nb_part = 2
#v = random_vit(nb_part, L_box, D)
v = np.array([[0,0], [0,0]])

position = np.empty((nb_pas, nb_part, D), np.float64)
vitesse = np.empty((nb_pas, nb_part, D), np.float64)
E_LJ = np.empty((nb_pas), np.float64)
E_LJwalls = np.empty((nb_pas), np.float64)
E_C = np.empty((nb_pas), np.float64)
t = np.empty((nb_pas), np.float64)
t2 = np.empty((nb_pas), np.float64)



increment = 0
for i in range(nb_pas):
    r, v = update(r, v, dt, m_part, nb_part, sig, eps, cutoff, D, L_box)
    E_LJ[i] = sumLJpot(r, sig, eps, cutoff, nb_part)
    E_LJwalls[i] = sumLJwalls(r, sig, eps, nb_part, L_box)
    E_C[i] = sumEC(v, m_part, nb_part)
    t2[i] = i*dt
    #r, v = reflectBC(r,v,nb_part,L_box,D,rayon)
    if i % save_interval == 0:
        position[increment] = r
        vitesse[increment] = v
        t[increment] = i*dt
        increment += 1
    if i % 1000 == 0:
        progress = round(i / nb_pas * 100, 2)
        print(f'\rAvancement calculs: {progress}%')
print(f'\rAvancement calculs: Fin')


position = position[:increment]
vitesse = vitesse[:increment]
moy_vit = np.mean(vitesse, axis=1)
moy_vit2 = np.array([(np.sqrt(moy_vit[i][0]**2 + moy_vit[i][1]**2)) for i in range(increment)])
Esum = E_C + E_LJwalls + E_LJ
#Esum = Esum - Esum[0]
diff = Esum[-1]-Esum[0]
t = t[:increment]

plt.figure(figsize=(12,8))
plt.plot(t2, Esum, 'b-', label='E tot')
plt.plot(t2, E_C, 'r-', label='E cin')
plt.plot(t2, E_LJwalls, 'g-', label='E pot LJ walls')
plt.plot(t2, E_LJ, 'p-', label='E pot LJ')
plt.xlabel('time')
plt.ylabel('energy')
plt.legend()
print(f'energy difference : {diff}')
plt.show()




#################################################
######## ANIMATION ###################################

def init():
    global position, rayon, nb_part, ax, particles
    
    for i in range(nb_part):
        circle = plt.Circle((-10, -10), rayon,facecolor='red', edgecolor='black', lw=1)
        ax.add_patch(circle)
        particles.append(circle)
    return particles 

def animate(i):
    global position, particles

    for j in range(nb_part):
        particles[j].center = (position[i,j, 0], position[i,j, 1])

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
anim.save('gaz.gif', writer=writergif)
#writergif = animation.FFMpegWriter(fps=20)
#anim.save('gaz.mp4', fps=20, dpi=200)
