import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def scientific_to_decimal(scientific_str):
    return float(scientific_str)

df = pd.read_csv(r'C:\GIT REPOS\L3-Mono\Fil Chaud\A0-28U32-3.csv', converters={'CH1': scientific_to_decimal})

t = df['X']
U = df['CH1']
dt = 10

U1 = U[:-dt]
U2 = U[dt:]


# Generate colors for each dot
colors = np.linspace(0, 1, len(U1))

# Plot the scatter plot with changing colors
plt.figure(figsize=(15,10))

plt.scatter(U1, U2, c=colors, cmap='viridis')
plt.colorbar(label='Time')
#plt.plot(U1,U2)

plt.xlabel('U1')
plt.ylabel('U2')
plt.title('Attractor')
#plt.yticks([-0.3,-0.2,-0.1,0,0.1,0.2,0.3])  # Remplacez ces valeurs par celles de votre choix
plt.grid()
plt.show()