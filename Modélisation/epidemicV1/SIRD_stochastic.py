import numpy as np
import matplotlib.pyplot as plt

N = 200
I0 = 10
R0 = 0
D0 = 0
S0 = N - I0 - R0 - D0

state_0 = np.array([S0, I0, R0, D0])
t_0 = 0
n_step = 100000

state = np.empty((n_step+1, 4))
t = np.empty(n_step+1)
state[0] = state_0
t[0] = t_0

beta = 0.3
gamma_r = 0.08
gamma_d = 0.10


def simulation(state, t, beta, gamma_r, gamma_d):
    p_IS = beta*state[0]*state[1]/N
    p_SR = gamma_r*state[1]
    p_SD = gamma_d*state[1]

    p_sum = p_IS + p_SR + p_SD
    tau = 1/p_sum*np.log(1/np.random.rand())
    t = t + tau

    rand = np.random.rand()

    if rand*p_sum <= p_IS:
        state += np.array([-1, 1, 0, 0])

    elif rand*p_sum > p_IS and rand*p_sum <= p_IS+p_SR:
        state += np.array([0, -1, 1, 0])

    else:
        state += np.array([0, -1, 0, 1])

    return state, t



for i in range(n_step):
  state[i+1], t[i+1] = simulation(state[i], t[i], beta, gamma_r, gamma_d)
state = state.T

plt.figure(figsize=(12,10))
plt.title("SIR model")
plt.plot(t, state[0], 'b-', label='Susceptible')
plt.plot(t, state[1], 'r-', label='Infected')
plt.plot(t, state[2], 'g-', label='Recovered with immunity')
plt.plot(t, state[3], 'b.', label='Deceased')
plt.grid()
plt.xlabel("Time [d]")
plt.ylabel("Numbers of individuals")
plt.legend(loc = "best")
plt.show()


#STACKPLOT
plt.figure(figsize=(12,10))
plt.title("SIR model")
names = ['Susceptible', 'Infected', 'Recovered with immunity', 'Deceased']
plt.stackplot(t, state[0], state[1], state[2], state[3], labels=names)
plt.xlabel("Time [d]")
plt.ylabel("Numbers of individuals")
plt.legend(loc = "best")
plt.show()
