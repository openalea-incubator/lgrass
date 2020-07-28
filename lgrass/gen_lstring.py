# coding: utf8
# '''
# Created on 02/06/2020
#
# @author: modelisation - TR
# '''

import os
import pickle
import pandas as pd
import lgrass

lstring_dir = "outputs/Sauvegardes_axiom/lstring/"
var_dir = "outputs/Sauvegardes_axiom/variables/"


# Sauvegarder une lstring dans un répertoire lstring_dir et les variables lgrass dans un répertoire var_dir
def save_lstring(lstring, lsystem):
    print lstring
    # la lstring
    if not os.path.isdir(lstring_dir):
        os.makedirs(lstring_dir)
    for mod in range(len(lstring)):
        # on créé un répertoire par module
        repertory = lstring_dir + str(mod)
        if not os.path.isdir(repertory):
            os.makedirs(repertory)
        # fichier avec le nom du module
        open(os.path.join(repertory, 'name' + str(mod) + '.txt'), 'w').write(lstring[mod].name)
        for i in range(len(lstring[mod])):
            # un fichier par paramètre présent dans le module (quel que soit le type)
            pickle.dump(lstring[mod][i], open(os.path.join(repertory, str(i) + '.txt'), 'w'))
    # les variables lgrass nécessaires à la récupération complète de la simulation
    if not os.path.isdir(var_dir):
        os.makedirs(var_dir)
    pickle.dump(lsystem.Biomasse_aerienne, open(os.path.join(var_dir, 'Biomasse_aerienne.txt'), 'w'))
    pickle.dump(lsystem.Biomasse_racinaire, open(os.path.join(var_dir, 'Biomasse_racinaire.txt'), 'w'))
    pickle.dump(lsystem.PourcentageRootGrowthRealized,
                open(os.path.join(var_dir, 'PourcentageRootGrowthRealized.txt'), 'w'))
    pickle.dump(lsystem.Reserve, open(os.path.join(var_dir, 'Reserve.txt'), 'w'))
    pickle.dump(lsystem.Demande_feuille, open(os.path.join(var_dir, 'Demande_feuille.txt'), 'w'))
    pickle.dump(lsystem.PourcentageFeuilGrowthRealized,
                open(os.path.join(var_dir, 'PourcentageFeuilGrowthRealized.txt'), 'w'))
    pickle.dump(lsystem.les_diametres, open(os.path.join(var_dir, 'les_diametres.txt'), 'w'))
    pickle.dump(lsystem.Seuil, open(os.path.join(var_dir, 'Seuil.txt'), 'w'))
    pickle.dump(lsystem.Taille_finale_gaine, open(os.path.join(var_dir, 'Taille_finale_gaine.txt'), 'w'))
    pickle.dump(lsystem.nb_talle, open(os.path.join(var_dir, 'nb_talle.txt'), 'w'))
    pickle.dump(lsystem.surface_foliaire_emergee, open(os.path.join(var_dir, 'surface_foliaire_emergee.txt'), 'w'))
    pickle.dump(lsystem.Surfoliaireencours, open(os.path.join(var_dir, 'Surfoliaireencours.txt'), 'w'))
    pickle.dump(lsystem.les_feuilles, open(os.path.join(var_dir, 'les_feuilles.txt'), 'w'))


# Charger et récupérer une lstring d'un répertoire dir
def load_lstring():
    modulenamelist = []  # contient tous les noms des modules
    moduleparam = []     # contient la liste des paramètres de chaque module
    # pour chaque dossier créé dans la simu sauvegardée
    for name in range(len([i for i in os.listdir(lstring_dir)])):
        df = pd.read_csv(os.path.join(lstring_dir + str(name), 'name' + str(name) + '.txt'))
        modulenamelist.append(df.columns[0])
        cpt = 0
        for mod in range(len([i for i in os.listdir(lstring_dir + str(name))])):
            cpt += 1
        if cpt == 1:
            moduleparam.append(None)
        else:
            param = []
            for j in range(cpt - 1):
                param.append(pickle.load(open(os.path.join(lstring_dir + str(name), str(j) + '.txt'), 'r')))
            moduleparam.append(param)
    return modulenamelist, moduleparam


# Récupérer les variables lgrass nécessaires à la récupération complète de la simulation
def load_variables():
    Biomasse_aerienne = pickle.load(open(os.path.join(var_dir, 'Biomasse_aerienne.txt'), 'r'))
    Biomasse_racinaire = pickle.load(open(os.path.join(var_dir, 'Biomasse_racinaire.txt'), 'r'))
    PourcentageRootGrowthRealized = pickle.load(open(os.path.join(var_dir, 'PourcentageRootGrowthRealized.txt'), 'r'))
    Reserve = pickle.load(open(os.path.join(var_dir, 'Reserve.txt'), 'r'))
    Demande_feuille = pickle.load(open(os.path.join(var_dir, 'Demande_feuille.txt'), 'r'))
    PourcentageFeuilGrowthRealized = pickle.load(open(os.path.join(var_dir, 'PourcentageFeuilGrowthRealized.txt'), 'r'))
    les_diametres = pickle.load(open(os.path.join(var_dir, 'les_diametres.txt'), 'r'))
    Seuil = pickle.load(open(os.path.join(var_dir, 'Seuil.txt'), 'r'))
    Taille_finale_gaine = pickle.load(open(os.path.join(var_dir, 'Taille_finale_gaine.txt'), 'r'))
    nb_talle = pickle.load(open(os.path.join(var_dir, 'nb_talle.txt'), 'r'))
    surface_foliaire_emergee = pickle.load(open(os.path.join(var_dir, 'surface_foliaire_emergee.txt'), 'r'))
    Surfoliaireencours = pickle.load(open(os.path.join(var_dir, 'Surfoliaireencours.txt'), 'r'))
    les_feuilles = pickle.load(open(os.path.join(var_dir, 'les_feuilles.txt'), 'r'))
    return Biomasse_aerienne, Biomasse_racinaire, PourcentageRootGrowthRealized, Reserve, Demande_feuille, PourcentageFeuilGrowthRealized, les_diametres, Seuil, Taille_finale_gaine, nb_talle, surface_foliaire_emergee, Surfoliaireencours, les_feuilles
