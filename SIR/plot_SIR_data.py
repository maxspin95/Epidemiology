import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) > 1:
    command_str = sys.argv[1]
else:
    print("Mentionner le nom du fichier csv après py plot_SIR_data.py.")

#importation des données
data_array = np.loadtxt(command_str, delimiter=",")
T = data_array[:,0]
S = data_array[:,1]
I = data_array[:,2]
R = data_array[:,3]
population_size = data_array[:,4]

#graphique des compartiments
plt.figure()
plt.scatter(T,S,color='b',marker='+',s=25,label=r"$S(t)$")
plt.scatter(T,I,color='r',marker='+',s=25,label=r"$I(t)$")
plt.scatter(T,R,color='k',marker='+',s=25,label=r"$R(t)$")
plt.xlabel(r"$t$")
plt.legend()
plt.show()
#graphique de l'erreur relative sur la taille de la population
plt.figure()
plt.scatter(T,100*np.abs(population_size-population_size[0])/population_size[0],color='k',marker='+',s=25)
plt.xlabel(r"$t$")
plt.ylabel(r"Erreur relative (%)")
plt.tight_layout()
plt.show()