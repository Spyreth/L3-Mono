import csv
import os


def csv_init(folder, name):
    """Crée les dossier pour les résultats"""
    if not os.path.exists(folder + name):
        os.makedirs(folder + name)



def write_to_csv(file, data):
    """Ecrit la data spécifiée sur une nouvelle ligne du fichier csv spécifié"""
    with open(file, 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(data)



def datasave(folder, name, r, v, t, D):
    """Utilise write_to_csv pour sauvegarder les données d'un pas de temps dans les fichiers csv
    Args:
        r (array): positions des particules
        v (array): vitesses des particules
        t (float): temps
        D (int): nb de dimensions
    """
    write_to_csv((folder + f"/{name}" + r"/time.csv"), [t])
    for d in range(D):
        write_to_csv((folder + f"/{name}" + fr"/pos_{d}.csv"), r.T[d])
        write_to_csv((folder + f"/{name}" + fr"/vit_{d}.csv"), v.T[d])


def pressureSave(folder, name, pressure):
    """Sauvegarde la pression dans un fichier csv avec write_to_csv"""
    write_to_csv((folder + f"/{name}" + r"/pressure.csv"), [pressure])


def save_parameters(folder, sim_name, L, D, nb_part, dt, m_part, nb_pas, sig, eps, cutoff, rayon, save_interval, pressure_calc_interval, kb):
    """Sauvegarde les paramètres de la simlulation dans un fichier txt à l'endroit spécifié"""
    with open((folder+sim_name+r"/param.txt"), 'w') as f:
        # Écriture des paramètres dans le fichier
        f.write(f"L: {L}\n")
        f.write(f"D: {D}\n")
        f.write(f"nb_part: {nb_part}\n")
        f.write(f"dt: {dt}\n")
        f.write(f"m_part: {m_part}\n")
        f.write(f"nb_pas: {nb_pas}\n")
        f.write(f"sig: {sig}\n")
        f.write(f"eps: {eps}\n")
        f.write(f"cutoff: {cutoff}\n")
        f.write(f"rayon: {rayon}\n")
        f.write(f"save_interval: {save_interval}\n")
        f.write(f"pressure_calc_interval: {pressure_calc_interval}\n")
        f.write(f"Kb: {kb}\n")




if __name__ == "__main__":
    print("Ceci n'est pas un script mais un package.")