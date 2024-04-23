import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Données
all_fishing_gear_diameters = [0.148, 0.165, 0.205, 0.235, 0.260, 0.285, 0.310, 0.330, 0.370, 0.405, 0.407, 0.520, 0.570, 0.620]
log_afg_d = np.log(np.array(all_fishing_gear_diameters))
all_fishing_gear_masses = [0.907, 1.361, 1.814, 2.268, 2.722, 3.629, 4.082, 6.350, 6.804, 7.711, 11.340, 14.515, 14.969, 19.958]
log_afg_m = np.log(np.array(all_fishing_gear_masses))

teklon_gold_advanced_masses = [3.1, 3.82, 4.78, 6.04, 7.06, 7.90, 8.94, 11, 15.2, 17.16]
log_tga_m = np.log(np.array(teklon_gold_advanced_masses))
teklon_gold_advanced_diameters = [0.172, 0.194, 0.222, 0.248, 0.272, 0.294, 0.314, 0.352, 0.416, 0.444]
log_tga_d = np.log(np.array(teklon_gold_advanced_diameters))

savage_gear_silencer_masses = [13.8, 4.19, 11.92, 8.97, 2.69, 15.56, 7.17, 3.33, 6.15, 1.8, 5.23]
log_sgs_m = np.log(np.array(savage_gear_silencer_masses))
savage_gear_silencer_diameters = [0.435, 0.235, 0.405, 0.35, 0.18, 0.465, 0.31, 0.20, 0.285, 0.15, 0.26]
log_sgs_d = np.log(np.array(savage_gear_silencer_diameters))

savage_xl_strong_g2_masses = [1.9, 2.2, 2.8, 3.3, 4.4, 5.4, 6.6, 7.7]
log_x1_m = np.log(np.array(savage_xl_strong_g2_masses))
savage_xl_strong_g2_diameters = [0.14, 0.16, 0.18, 0.2, 0.23, 0.25, 0.28, 0.3]
log_x1_d = np.log(np.array(savage_xl_strong_g2_diameters))

sufix_advance_g2_clear_diameters = [0.3, 0.35, 0.4]
log_g2_d = np.log(np.array(sufix_advance_g2_clear_diameters))
sufix_advance_g2_clear_masses = [7.1, 9.1, 11.3]  # Ajout de valeurs None pour faire correspondre les tailles des listes
log_g2_m = np.log(np.array(sufix_advance_g2_clear_masses))


# Tracé pour chaque modèle
plt.figure(figsize=(12,8))
plt.scatter(log_afg_d, log_afg_m, label='All Fishing Gear', color='blue')
plt.scatter(log_tga_d, log_tga_m, label='Teklon Gold Advanced', color='red')
plt.scatter(log_sgs_d, log_sgs_m, label='Savage Gear Silencer', color='green')
plt.scatter(log_x1_d, log_x1_m, label='Savage XL Strong G2', color='purple')

# Étiquettes d'axe et titre
plt.xlabel('log(Diamètre)')
plt.ylabel('log(Masse de rupture)')
plt.title('Graphe log-log de la masse de rupture en fonction du diamètre des fils de pêche')

# Légende
plt.legend()
plt.grid(True, linestyle=':', linewidth='0.5')

#linreg
slope_afg, intercept_afg, r_value, p_value, std_err = linregress(log_afg_d, log_afg_m)
slope_tga, intercept_tga, r_value, p_value, std_err = linregress(log_tga_d, log_tga_m)
slope_sgs, intercept_sgs, r_value, p_value, std_err = linregress(log_sgs_d, log_sgs_m)
slope_x1, intercept_x1, r_value, p_value, std_err = linregress(log_x1_d, log_x1_m)
slope_g2, intercept_g2, r_value, p_value, std_err = linregress(log_g2_d, log_g2_m)

# print linreg
models = {
    'All Fishing Gear': [slope_afg,intercept_afg],
    'Teklon Gold Advanced': [slope_tga,intercept_tga],
    'Savage Gear Silencer': [slope_sgs,intercept_sgs],
    'Savage XL Strong G2': [slope_x1,intercept_x1],
    'Sufix Advance G2 Clear': [slope_g2,intercept_g2]
}

# Afficher les coefficients de régression
for model_name, model in models.items():
    print(f"Marque: {model_name}")
    print(f"Coefficient de pente: {model[0]}")
    print(f"Constante (intercept): {model[1]}")
    print("\n")

plt.plot(np.linspace(log_afg_d.min(), log_afg_d.max(), 100), np.linspace(log_afg_d.min(), log_afg_d.max(), 100)*slope_afg + intercept_afg, color='blue', linestyle='--')
plt.plot(np.linspace(log_tga_d.min(), log_tga_d.max(), 100), np.linspace(log_tga_d.min(), log_tga_d.max(), 100)*slope_tga + intercept_tga, color='red', linestyle='--')
plt.plot(np.linspace(log_sgs_d.min(), log_sgs_d.max(), 100), np.linspace(log_sgs_d.min(), log_sgs_d.max(), 100)*slope_sgs + intercept_sgs, color='green', linestyle='--')
plt.plot(np.linspace(log_x1_d.min(), log_x1_d.max(), 100), np.linspace(log_x1_d.min(), log_x1_d.max(), 100)*slope_x1 + intercept_x1, color='purple', linestyle='--')
plt.show()


log_all_d = np.concatenate((log_afg_d, log_tga_d, log_sgs_d, log_x1_d))
log_all_m = np.concatenate((log_afg_m, log_tga_m, log_sgs_m, log_x1_m))
slope_all, intercept_all, r_value, p_value, std_err = linregress(log_all_d, log_all_m)


print(f"Tous les fils :")
print(f"Coefficient de pente: {slope_all}")
print(f"Constante (intercept): {intercept_all}")
print("\n")


plt.figure(figsize=(12,8))
plt.grid(True, linestyle=':', linewidth='0.5')
plt.xlabel('log(Diamètre)')
plt.ylabel('log(Masse de rupture)')
plt.scatter(log_all_d, log_all_m, label='All fishing lines', color='blue', s=10)
plt.plot(np.linspace(log_all_d.min(), log_all_d.max(), 100), np.linspace(log_all_d.min(), log_all_d.max(), 100)*slope_all + intercept_all, color='red', linestyle='--', label ='linear fit')
plt.title('Graphe log-log de la masse de rupture en fonction du diamètre des fils de pêche')
plt.legend()
plt.show()


