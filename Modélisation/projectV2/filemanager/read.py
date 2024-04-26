import pandas as pd
import numpy as np


def get_posvittime(folder, D, nb_part, nb_pas, save_interval):
    nb_pas_saved = int(nb_pas/save_interval)

    t = pd.read_csv((folder + r"\time.csv"), header=None, skiprows=1).values

    r = np.empty((nb_pas_saved, nb_part, D), np.float64)
    v = np.empty((nb_pas_saved, nb_part, D), np.float64)

    for d in range(D):
        r_data = pd.read_csv((folder + fr"\pos_{d}.csv"), header=None, skiprows=1).values
        r[:, :, d] = r_data.reshape(nb_pas_saved, nb_part)
        v_data = pd.read_csv((folder + fr"\vit_{d}.csv"), header=None, skiprows=1).values
        v[:, :, d] = v_data.reshape(nb_pas_saved, nb_part)

    return r,v,t


def get_param(folder):
    param = {}
    with open(folder, 'r') as f:
        for ligne in f:
            name, value = ligne.strip().split(": ")
            param[name] = value

    return param


def get_pressure_billard(folder):
    pressure = pd.read_csv((folder + r"\pressure.csv"), header=None).values
    return pressure



if __name__ == "__main__":
    print("Ceci n'est pas un script mais un package.")
    