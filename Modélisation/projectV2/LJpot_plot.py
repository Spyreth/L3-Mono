import numpy as np
import matplotlib.pyplot as plt


eps = 1
sig = 1

def force(x):
    return -24*eps*((sig**6)*(1/(x**7)) - 2*(sig**12)*(1/x**13))

def pot(x):
    return 4*eps*((sig/x)**12-(sig/x)**6)

x = np.linspace(0, 4, 1500)
f = np.empty((np.size(x)), dtype=np.float64)
p = np.empty((np.size(x)), dtype=np.float64)
incr=0
for elem in x:
    f[incr] = force(x[incr])
    p[incr] = pot(x[incr])
    incr += 1

plt.figure(figsize=(12,8))
plt.xlim(0, 4)
plt.ylim(-4, 5)
plt.plot(x, f, 'r-', label='force')
plt.plot(x, p, 'b-', label='potentiel')
plt.xlabel('r (unités arbitraires)')
plt.ylabel('E (unités arbitraires)')
plt.legend()
plt.title("Potentiel de Lennard-Jones et force associée pour eps = 1 et sig = 1")
plt.grid()
plt.show()