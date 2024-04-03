import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

#################################################
######## ANIMATION ###################################

def init():
    global position, rayon, N, ax, particles
    
    for i in range(N):
        circle = plt.Circle((-10, -10), rayon[i],facecolor='red', edgecolor='black', lw=1)
        ax.add_patch(circle)
        particles.append(circle)
    return particles 

def animate(i):
    global position, particles

    for j in range(N):
        particles[j].center = (position[i,j, 0], position[i,j, 1])

    return particles 

# N particules, N frames, coordonnées 2D:
position = np.empty((Nframes, N, 2), np.float64)

# SIMULATION qui remplit position...

# création de l'animation à partir du tableau position:

fig, ax = plt.subplots()
ax.set_xlim(0, L)
ax.set_ylim(0, L)
ax.set_aspect('equal')
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)  

particles = []

anim = animation.FuncAnimation(fig, animate, frames=Nframes, interval=10, blit=True, init_func=init, repeat=False)

#writergif = animation.PillowWriter(fps=20)
#anim.save('gaz.gif', writer=writergif)
writergif = animation.FFMpegWriter(fps=20)
anim.save('gaz.mp4', fps=20, dpi=200)
#plt.show()
