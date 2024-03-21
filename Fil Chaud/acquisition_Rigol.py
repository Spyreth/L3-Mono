# -*- coding: utf-8 -*-

"""

Created on Sat Nov 28 11:39:45 2020



@author: Romain Long

"""





import pyvisa as visa

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt

#import scipy.optimize as opt # Pour l’ajustement

#from scipy.fft import fft, fftfreq



rm = visa.ResourceManager()

ls = rm.list_resources()

print(ls)

oscillo = rm.open_resource(ls[0])           # Ouvre l'oscilloscope généralement le premier de la liste ( ls[0] )



oscillo.write(':WAV:SOURCE CHAN1')          # Selection de la voie 1

oscillo.write(':WAV:FORMAT ASCII')          # Données en format texte

answer1 = oscillo.query(':WAV:DATA?')        # Récupère les données

data1 = answer1[11:]                          # élimine l'entête (les 11 premiers caractères)

y1 = np.float32( data1.split(',') )           # découpe les données et les convertit en float





# oscillo.write(':WAV:SOURCE CHAN2')          # Selection de la voie 1

# oscillo.write(':WAV:FORMAT ASCII')          # Données en format texte

# answer2 = oscillo.query(':WAV:DATA?')        # Récupère les données

# data2 = answer2[11:]                          # élimine l'entête (les 11 premiers caractères)

# y2 = np.float32( data2.split(',') )           # découpe les données et les convertit en float



answer = oscillo.query(':TIM:SCALE?')       # Récupère la base de temps (durée par division)

TimeScale = float(answer)

t = np.linspace(0, TimeScale * 12, len(y1) ) # Création du tableau des abscisses de même taille que y 

                                            # l'écran de l'oscillo fait 12 divisions

                                            

plt.figure()

plt.grid(axis='both', alpha=0.75)

plt.title('Oscillateurs couplés')

plt.xlabel("Temps (ms)")

plt.ylabel("Tension (V)")

plt.plot(1000*t, y1, label='Tension condensateur')

#plt.scatter(1000*t_util, y_util, label='experimental data')

#plt.plot( 1000*t, y2, label= 'Tension_generateur')

#plt.plot( 1000*t_util, y_fit, 'red', label= 'Ajustement')

plt.legend()

plt.show()

nval=np.shape(t)[0]

print('nval=',nval)







#data = pd.DataFrame({'Time (s)': t , 'CH1 (V)': y1 , 'CH2 (V)': y2 })   # Mets en forme les données avec entête



                                     

'''                                        

data = pd.DataFrame({'Time (s)': t , 'CH1 (V)': y1  })   # Mets en forme les données avec entête

data.to_csv('RC.csv', index = False)     # Sauvegarde les données dans un fichier CSV





#data = pd.read_csv('RC7.csv')     # Lecture d'un fichier pour réanalyse

                                            







t = data['Time (s)']

y1 = data['CH1 (V)']

'''

#y2 = data['CH2 (V)']

#t_util = t[365:665]-t[365] # Sélection des données : l’exponentielle croissante

#y_util = y1[365:665]





# def f(t,V,V0,tau) : # Fonction d’ajustement (charge d’un condensateur)

#     return V + (V0 - V ) * np.exp( -t/tau )



# init_param = [5 , 0 ,4.0e3 * 100e-9] # Valeurs initiales V=5V, V0=0V tau=RC

# final_param , var = opt.curve_fit(f, t_util, y_util, init_param)



# print('V: ', final_param[0],'\n V0: ', final_param[1], '\n tau:  ', final_param[2])



#y_fit = f(t_util,*final_param)





# plt.figure()

# plt.grid(axis='both', alpha=0.75)

# plt.title('Charge condensateur_100micros')

# plt.xlabel("Temps (ms)")

# plt.ylabel("Tension condensateur (V)")

# plt.plot(1000*t_util, y_util, label='experimental data')

# #plt.scatter(1000*t_util, y_util, label='experimental data')

# plt.plot( 1000*t_util, y_fit, 'red', label= 'Ajustement V= {:.2f}+ ( {:.2f} - {:.2f})*exp(-t/{:.5f})'.format(final_param[0], final_param[1],final_param[0],final_param[2]))

# #plt.plot( 1000*t_util, y_fit, 'red', label= 'Ajustement')

# plt.legend()

# plt.show()

# plt.savefig('Charge condensateur_100micros.png')

# plt.savefig('Charge condensateur_100micros.pdf')
