import numpy as np
import matplotlib.pyplot as plt

N = 8700000
I0 = 10
R0 = 0
D0 = 0
S0 = N - I0 - R0 - D0

state_i = np.array([S0, I0, R0, D0])
t_i = 0
t_f = 200
n_step = 200
t = np.linspace(t_i, t_f, n_step+1)
dt = t[1] - t[0]

beta = 0.3
gamma_r = 0.08
gamma_d = 0.02


def derivatives(state):
  dS = -beta*state[0]*state[1]/N
  dI = beta*state[0]*state[1]/N - gamma_r*state[1] - gamma_d*state[1]
  dR = gamma_r*state[1]
  dD = gamma_d*state[1]
  return np.array([dS, dI, dR, dD])


def rk4(f, state0, dt):
  k0 = dt * f(state0) 
  k1= dt * f(state0 + 0.5 * k0)
  k2 = dt * f(state0 + 0.5 * k1)
  k3 = dt * f(state0 + k2)
  state1 = state0 + (k0 + 2*(k1 + k2) + k3) / 6.
  return state1


state = np.empty((n_step+1, 4))
state[0] = state_i

for i in range(n_step):
  state[i+1] = rk4(derivatives, state[i], dt)
state = state.T

#BASIC PLOT
plt.figure(figsize=(12,10))
plt.title("SIR model")
plt.plot(t, state[0], 'b-', label='Susceptible')
plt.plot(t, state[1], 'r-', label='Infected')
plt.plot(t, state[2], 'g-', label='Recovered with immunity')
plt.plot(t, state[3], 'b.', label='Deceased')
plt.grid()
plt.xlabel("Time, $t$ [s]")
plt.ylabel("Numbers of individuals")
plt.legend(loc = "best")
plt.show()



#STACKPLOT
plt.figure(figsize=(12,10))
plt.title("SIR model")
names = ['Susceptible', 'Infected', 'Recovered with immunity', 'Deceased']
plt.stackplot(t, state[0], state[1], state[2], state[3], labels=names)
plt.xlabel("Time, $t$ [s]")
plt.ylabel("Numbers of individuals")
plt.legend(loc = "best")
plt.show()
