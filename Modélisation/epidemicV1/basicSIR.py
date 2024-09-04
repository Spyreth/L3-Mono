import numpy as np
import matplotlib.pyplot as plt

N = 8700000
I0 = 10
R0 = 0
S0 = N - I0 - R0

state = np.array([S0, I0, R0])
t_i = 0
t_f = 200
n_step = 200

beta = 0.3
gamma = 0.2

def derivatives(state):
  dS = -beta*state[0]*state[1]/N #dY/dt = -aS(t)I(t)
  dI = beta*state[0]*state[1]/N - gamma*state[1]  #dI/dt = aS(t)I(t)-bI(t)
  dR = gamma*state[1]  #dR/dt = bI(t)
  return np.array([dS, dI, dR])


def rk4(f, state0, t0, tf, n):
  t = np.linspace(t0, tf, n+1) #x grid
  s = np.array((n+1)*[state0]) #array of the state of the sytem for each x
  h = t[1] - t[0] #stepsize
  for i in range(n):
    k0 = h * f(s[i]) 
    k1= h * f(s[i] + 0.5 * k0)
    k2 = h * f(s[i] + 0.5 * k1)
    k3 = h * f(s[i] + k2)
    s[i+1] = s[i] + (k0 + 2*(k1 + k2) + k3) / 6.
  return t, s

t, s = rk4(derivatives, state, t_i, t_f, n_step)
s = s.T

"""  #BASIC PLOT
plt.figure(figsize=(12,10))
plt.title("SIR model")
plt.plot(t, s[0], 'b-', label='Susceptible')
plt.plot(t, s[1], 'r-', label='Infected')
plt.plot(t, s[2], 'g-', label='Recovered with immunity')
plt.grid()
plt.xlabel("Time, $t$ [s]")
plt.ylabel("Numbers of individuals")
plt.legend(loc = "best")
plt.show()
"""

#STACKPLOT
plt.figure(figsize=(12,10))
plt.title("SIR model")
names = ['S', 'I', 'R']
plt.stackplot(t, s[0], s[1], s[2], labels=names)
plt.xlabel("Time, $t$ [s]")
plt.ylabel("Numbers of individuals")
plt.legend(loc = "best")
plt.show()