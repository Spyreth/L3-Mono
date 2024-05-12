PROJET PHYNUM :

Les simulations peuvent prendre beaucoup de temps (jusqua une vingtaine d'heures sur nos ordinateurs pour le script simu GP1.py), 
nous avons sélectionné les données des simulations pertinentes qui sont incluses dans les dossiers "Résultats" et "Résultatsbillard".

Résultats:
- testEnergie : conservation de l'energie (section 5.1)
- GP2_12x12_L15 : tous les dossiers pour la loi des gp et l'équation d'état de vdw avec interactions
- Vitesses_Maxwell : distribution des vitesses gif (section 6)
Résultatsbillard:
- GP2_12x12_L15 : tous les dossiers pour la loi des gp sans interactions

Les scripts d'analyse sont configurés pour analyser ces dossiers directement sans changer de paramètre :
- analyse.py -> configuré pour testEnergie
- analyse GP1.py -> configuré pour l'analyse GP et vdw de GP2_12x12_L15 avec interactions
- analyse GP1 billard.py -> configuré pour l'analyse GP de GP2_12x12_L15 sans interactions

Pour lancer une nouvelle simulation, il est recommandé de se référer au rapport afin de comprendre le fonctionnement de chaque script et fonction.
Il suffit ensuite de modifier les paramètres de simulation du script que l'on veut lancer, et de modifier la variable "results_name". Les données de position,
de vitesse, de temps et de pression seront sauvegardées dans un nouveau dossier portant le nom de la variable, se trouvant dans le dossier "save_folder" à
l'intérieur du projet.

Pour lancer les analyses souhaitées, il faut ensuite modifier la variable "results_name" dans le script d'analyse voulu et lancer le script. Les diférents graphes
produit par le script sont alors sauvegardés dans le dossier "results_name" correspondant.
