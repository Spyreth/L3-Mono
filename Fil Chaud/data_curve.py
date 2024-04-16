import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def scientific_to_decimal(scientific_str):
    return float(scientific_str)

df = pd.read_csv(r'C:\GIT REPOS\L3-Mono\Fil Chaud\A0-28U32-3.csv', converters={'CH1': scientific_to_decimal})

t = df['X']
U = df['CH1']
print(U)

plt.figure(figsize=(15,10))
plt.plot(t, U)
plt.xlabel('X')
plt.ylabel('CH1')
plt.title('Plot des données X vs CH1')
plt.yticks([-0.3,-0.2,-0.1,0,0.1,0.2,0.3])  # Remplacez ces valeurs par celles de votre choix
plt.grid()
plt.show()