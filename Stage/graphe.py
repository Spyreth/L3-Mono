import numpy as np
import matplotlib.pyplot as plt



### CONSTANTES
gamma = np.array([1.16, 1.32, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66])
E_S = np.array([0.16, 0.28, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60])
E_N = np.array([0.50, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40])
E_R = 0.05
E_ES12 = 0.10
E_ES23 = 0
E_ES34 = 0
E_ES45 = 0
E_ESRESTE = 3
E_ES = np.array([[0.0, E_ES12, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE],
                [E_ES12, 0.0, E_ES23, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE],
                [E_ESRESTE, E_ES23, 0.0, E_ES34, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE],
                [E_ESRESTE, E_ESRESTE, E_ES34, 0.0, E_ES45, E_ESRESTE, E_ESRESTE, E_ESRESTE],
                [E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ES45, 0.0, E_ESRESTE, E_ESRESTE, E_ESRESTE],
                [E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, 0.0, E_ESRESTE, E_ESRESTE],
                [E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, 0.0, E_ESRESTE],
                [E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, E_ESRESTE, 0.0],
                ])
### CONSTANTES


def calc_maxmin(h):

    shape = h.shape[0]
    max = np.zeros(shape+1)
    min = np.zeros(shape)
    maxmin = np.zeros(2*shape+1)

    for i in range(shape):
        hprev = h[(i-1)%shape]
        hcurr = h[i]
        hnext = h[(i+1)%shape]

        if (hcurr<hnext and hcurr>hprev) or (hcurr>hnext and hcurr<hprev):
            min[i] = -gamma[hcurr]-E_N[hcurr]-E_R
        elif hcurr<hnext and hcurr>hprev:
            min[i] = -gamma[hcurr]-2*E_N[hcurr]
        elif hcurr<hnext or hcurr<hprev:
            min[i] = -gamma[hcurr]-E_N[hcurr]
        elif hcurr>hnext or hcurr>hprev:
            min[i] = -gamma[hcurr]-E_R
        else :
            min[i] = -gamma[hcurr]
        
        if i == 0:
            if hcurr<hprev:
                max[i] = -gamma[hprev]+E_S[hprev]+E_ES[hprev,hcurr]
            elif hcurr>hprev:
                max[i] = -gamma[hcurr]+E_S[hcurr]+E_ES[hprev,hcurr]
            else:
                max[i] = -gamma[i] + E_S[i]

        if hcurr<hnext:
            max[i+1] = -gamma[hnext]+E_S[hnext]+E_ES[hnext,hcurr]
        elif hcurr>hnext:
            max[i+1] = -gamma[hcurr]+E_S[hcurr]+E_ES[hnext,hcurr]
        else :
            max[i+1] = -gamma[hcurr]+E_S[hcurr]
    
    maxmin[0::2] = max
    maxmin[1::2] = min
    return maxmin


def calc_coeffsin(max,min):
    a = (max-min)/2
    b = (max+min)/2
    return a,b


def plot_func(h):

    x = np.linspace(-0.5, h.size-0.5, h.size*2+1)
    maxmin = calc_maxmin(h)
    hmax = np.max(h)
    potmin = np.min(maxmin)
    potmax = np.max(maxmin)

    fig, (ax1,ax2) = plt.subplots(2,1,figsize = (12,10), dpi=100)
    ax1.grid(color='black', linewidth=0.4, alpha = 0.8)
    #ax1.axis("off")
    ax1.set_ylabel("eV")
    ax1.set_xlabel("position")
    ax1.set_xlim(x[0],x[-1])
    ax1.set_ylim(potmin-0.02,potmax+0.02)

    aspect_ratio = (hmax+0.6)/h.size
    ax2.set_box_aspect(aspect_ratio)
    ax2.set_xlim(x[0],x[-1])
    ax2.set_ylim(-0.3,hmax+0.3)
    ax2.set_xlabel("position")
    ax2.set_yticks([])
    #ax2.axis("off")
    ax2.plot(x,np.zeros(x.size), color='black')

    z = np.linspace(-np.pi/2,np.pi/2,200)
    sin = np.sin(z)

    incr = 0
    for elem in h:
        for i in range(elem):
            circle = plt.Circle((incr,i+0.5), 0.50, edgecolor='black', facecolor='black')
            ax2.add_patch(circle)
        incr += 1


    for i in range(2*h.size):

        if i%2 == 0:
            max = maxmin[i]
            min = maxmin[i+1]
            a, b = calc_coeffsin(max,min)
            arr = np.linspace(x[i], x[i+1], 200)
            ax1.plot(arr,-a*sin + b, "black")

        else:
            min = maxmin[i]
            max = maxmin[i+1]
            a, b = calc_coeffsin(max,min)
            arr = np.linspace(x[i], x[i+1], 200)
            ax1.plot(arr,a*sin + b, "black")

    plt.show()


h = np.array([0,0,0,2,1,1,1,2,1,0,0,0])
plot_func(h)
