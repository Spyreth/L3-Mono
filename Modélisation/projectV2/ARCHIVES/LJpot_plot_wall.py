import numpy as np
import matplotlib.pyplot as plt


eps = 1
sig = 1
cutoff = 2**(1/6)*sig

def force(x):
    if x<cutoff:
        return -24*eps*((sig**6)*(1/(x**7)) - 2*(sig**12)*(1/x**13)) 
    else:
        return 0

def pot(x):
    if x<cutoff:
        return 4*eps*((sig/x)**12-(sig/x)**6) + eps
    else:
        return 0

x = np.linspace(0.8, 1.5, 1500)
f = np.empty((np.size(x)), dtype=np.float64)
p = np.empty((np.size(x)), dtype=np.float64)
incr=0
for elem in x:
    f[incr] = force(x[incr])
    p[incr] = pot(x[incr])
    incr += 1

plt.figure(figsize=(12,8))
plt.xlim(0.8, 1.5)
plt.ylim(-1, 4)
plt.plot(x, f, 'r-', label='force')
plt.plot(x, p, 'b-', label='potentiel')
plt.xlabel('r (unités arbitraires)')
plt.ylabel('E (unités arbitraires)')
plt.legend()
plt.grid()
plt.title("Potentiel de Lennard-Jones tronqué et force associée pour eps = 1 et sig = 1")
plt.show()