import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def scientific_to_decimal(scientific_str):
    return float(scientific_str)

<<<<<<< HEAD
df = pd.read_csv(r'D:\GIT REPOS\L3 Mono\L3-Mono\Fil Chaud\deuxiememesure.csv', converters={'CH1': scientific_to_decimal})
=======
df = pd.read_csv(r'C:\GIT REPOS\L3-Mono\Fil Chaud\A0-28U32-3.csv', converters={'CH1': scientific_to_decimal})
>>>>>>> 705bc1e090a08b0a8edf47c2fe8472ee9d4f184b

t = df['X']
U = df['CH1']
dt = 10

U1 = U[:-dt]
U2 = U[dt:]

<<<<<<< HEAD
colors = np.linspace(0, 1, len(U1))
plt.figure(figsize=(15,10))
plt.scatter(U1, U2, c=colors, cmap='viridis')
=======

# Generate colors for each dot
colors = np.linspace(0, 1, len(U1))

# Plot the scatter plot with changing colors
plt.figure(figsize=(15,10))

plt.scatter(U1, U2, c=colors, cmap='viridis')
plt.colorbar(label='Time')
#plt.plot(U1,U2)

>>>>>>> 705bc1e090a08b0a8edf47c2fe8472ee9d4f184b
plt.xlabel('U1')
plt.ylabel('U2')
plt.title('Attractor')
#plt.yticks([-0.3,-0.2,-0.1,0,0.1,0.2,0.3])  # Remplacez ces valeurs par celles de votre choix
plt.grid()
plt.show()