import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animate


def calc(n):
    x = np.random.rand(n)
    y = np.random.rand(n)
    dist = np.sqrt(x**2 + y**2)
    i = 0

    for d in dist:
        if d < 1:
            i += 1

    pi = 4*i/n
    return pi


def evol(start, end, pas):

    arr_pi = np.array([])
    arr_n = np.array(range(start,end,pas))
    i = 0

    for elem in arr_n :
        i += 1
        arr_pi = np.append(arr_pi, calc(elem))
        if i*pas/end*100 %5 == 0:
            print(f"{i*pas/end*100}% done")
    
    print(f"{100}% done")
    return arr_n, arr_pi


time = 200
fps = 30
nb_frames = time*fps
interval_vid = time/nb_frames
n_start = 100
n_end = 100000
step = (n_end-n_start)/nb_frames

arr_n, arr_pi = evol(int(n_start), int(n_end), int(step))
plt.xscale('log')
plt.xlabel('nb de points aléatoires')
plt.ylabel('valeur approchée de pi')
plt.scatter(arr_n, arr_pi)
plt.show()

def def_fig(n_start, n_end):
    fig, ax = plt.subplots(1,1, figsize = (20,4), dpi = 300)
    fig.set_facecolor("black")
    ax.set_facecolor("black")
    ax.set_xlim(n_start, n_end)
    ax.set_ylim(np.min(arr_pi), np.max(arr_pi))

    return fig, ax

fig, ax = def_fig(int(n_start), int(n_end))


def animate_graph(n_current):
    plt.scatter(arr_n[:np.size(arr_n)-n_current], arr_pi[:np.size(arr_pi)-n_current])
    plt.show()

"""
anim = animate.FuncAnimation(
    fig,
    animate_graph,
    frames = nb_frames,
    interval = interval_vid
)
anim.save("test.mp4")
"""