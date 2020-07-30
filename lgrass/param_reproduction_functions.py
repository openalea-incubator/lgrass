# coding: utf8
# '''
# Created on 05/04/2020
#
# @author: modelisation - TR
# '''

import random
import numpy as np
import pandas as pd
import math
import shutil
import os
import flowering_functions


# Créer la matrice de croisement des plantes pour établir une nouvelle génération de nb_plantes
def create_seeds(lstring, nb_plantes, opt_repro, cutting_freq, ParamP):
    matrix = np.zeros((nb_plantes, nb_plantes))
    seeds = []
    elected_seeds = []

    # construction des graines et de leurs parents
    mothers = []
    for mod in lstring:
        if mod.name in ('apex',):
            if mod[0].final_spikelet_number is not None:
                mothers.append((mod[0].id_plante, mod[0].final_spikelet_number))

    # Méthode de calcul du nombre de graines via les regressions graines/tiges du rapport de Pierre Guinard (2012)
    if opt_repro == "SPPR_2012":
        if cutting_freq <= 21:  # Coupes fréquentes
            nb_seeds = int(9.73*len(mothers) - 38.19)
        else:   # Coupes peu fréquentes
            nb_seeds = int(6.35*len(mothers) - 13.79)
        if nb_seeds >= 1:
            for i in range(nb_seeds):
                id_mother = random.Random().choice(mothers)[0]
                fathers = [j for j in range(nb_plantes)]
                fathers.remove(id_mother)
                if fathers is not None:
                    id_father = random.Random().choice(fathers)  # séléction aléatoire du père
                    seeds.append((id_mother, id_father))
                else:
                    raise NameError("La génération ne comporte qu'une seule plante, la reproduction est impossible.")
        else:
            raise NameError("Il n'y a pas eu suffisamment de graines produites pour établir une nouvelle génération.")

    # Méthode de calcul via un nombre de graines par épillet
    elif opt_repro == "spikelets":
        for k in range(len(mothers)):
            for i in range(int(mothers[k][1] * int(ParamP[mothers[k][0]]['seeds_by_spikelet']))):   # nb_epillets de la plante * nb_graines/épillet
                id_mother = mothers[k][0]  # séléction de la mère
                fathers = [j for j in range(nb_plantes)]
                fathers.remove(id_mother)
                if fathers is not None:
                    id_father = random.Random().choice(fathers)  # séléction aléatoire du père
                    seeds.append((id_mother, id_father))
                else:
                    raise NameError("La génération ne comporte qu'une seule plante, la reproduction est impossible.")

    # sélection aléatoire des graines pour la génération suivante et création de la matrice de croisement
    if len(seeds) < nb_plantes:
        raise NameError("Il n'y a pas eu suffisamment de graines produites pour établir une nouvelle génération.")

    for rand in range(nb_plantes):
        seed = random.Random().choice(seeds)
        elected_seeds.append(seed)
        seeds.remove(seed)
        # id_mere/ligne et id_pere/colonne
        matrix[elected_seeds[rand][0], elected_seeds[rand][1]] += 1
    return matrix


# Lire le fichier de données génétique (.r) et le formater en dataframe
def get_genet_file(in_genet_file=None):
    if in_genet_file is None:
        in_genet_file = 'modelgenet/ped.r'

    infile = pd.read_csv(in_genet_file, header=None)
    df = pd.DataFrame(columns=['geno', 'B', 'G', 'D', 'E', 'F', 'C', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'])
    data = []
    for i in range(len(infile)):
        line = []
        k = 0
        while k < len(infile[0][i]):
            if infile[0][i][k] != ' ':
                chain = ''
                while infile[0][i][k] != ' ':
                    chain += infile[0][i][k]
                    k += 1
                    if k == len(infile[0][i]):
                        break
                line.append(chain)
            else:
                k += 1
        data.append(line)
    for i in range(len(data)):
        df.loc[i] = data[i]
    return df


# La formule pour convertir la donnée génétique en paramètre C
def calculate_C(x):
    # La version de cette fonction est incomplète (Thibault Raquet), voir Leopoldo Sanchez Rodriguez pour compléter la formule
    # se base sur une observation des sorties génétiques dont le paramètre variait entre -4 et 4, à convertir en [0.8, 1.6]
    return 0.1 * x + 1.2


# Création du fichier de paramètres d'entrée pour chaque plante et configuration du planteur lgrass
def define_param(in_param_file=None, in_genet_file=None, out_param_file=None, id_gener=1, opt_repro=None):
    if in_param_file is None:
        in_param_file = 'inputs/liste_plantes.csv'
    if out_param_file is None:
        out_param_file = 'outputs/Simulation' + str(id_gener) + '.csv'
    if opt_repro != "False":
        infile = get_genet_file(in_genet_file=in_genet_file)
        genet_data = infile.loc[infile['D'] == str(id_gener), :]
    data = pd.read_csv(in_param_file)
    # Création du fichier de paramètres d'entrée de chaque plante
    param_init = open(out_param_file, 'w')
    param_name = list(data.columns)
    for par in range(len(param_name)):
        L = [param_name[par]]
        for i in range(len(data)):
            L.append(str(data[param_name[par]].iloc[i]))
        param_init.write(";".join(L) + "\n")
    param_init.close()

    # lecture du fichier, creation de ParamP et des parametres de floraison
    param_plante = pd.read_csv(out_param_file, sep=";", header=None, index_col=0)
    # Conversion du paramètre génétique en valeur de C
    if opt_repro != 'False':
        if len(param_plante.columns) < len(genet_data):
            print("Attention, des plantes sont manquantes dans le fichier de paramètres pour l'application du modèle génétique")
            print("Des copies de la dernière plante enregistrée ont été ajoutées pour completer la simulation")
            for line in range(len(param_plante.columns), len(genet_data)):
                param_plante[line+1] = param_plante[line]
        for i in range(len(genet_data)):
            param_plante.loc['C'].iloc[i] = calculate_C(float(genet_data['C'].iloc[i]))
            param_plante.loc['geno'].iloc[i] = int(genet_data['geno'].iloc[i])
        param_plante.to_csv(out_param_file, sep=';')

    flowering_param = flowering_functions.FloweringFunctions()
    flowering_param.param.__dict__.update(
        (k, list(param_plante.loc[k, :])) for k in param_plante.index.intersection(flowering_param.param.__dict__))

    ParamP = list(dict(zip(param_plante.index, param_plante.iloc[:, col])) for col in range(len(param_plante.columns)))
    # Creation des matrices d'identifiant des plantes et de leur genotype
    nb_plantes = len(ParamP)
    NBlignes = int(math.ceil(np.sqrt(nb_plantes)))
    NBcolonnes = int(math.floor(np.sqrt(nb_plantes)))
    posPlante = [[i, j] for i, j in zip(sorted(list(range(NBlignes)) * NBcolonnes), list(range(NBcolonnes)) * NBlignes)]
    Plantes = np.arange(nb_plantes).reshape(NBlignes, NBcolonnes)
    Genotypes = np.array([i for i in param_plante.loc['geno']]).reshape(NBlignes, NBcolonnes)
    return ParamP, nb_plantes, NBlignes, NBcolonnes, posPlante, Plantes, Genotypes, flowering_param


# Remplir le fichier d'entree du modèle génétique avec la matrice de croisement et executer le modèle
def rungenet(src, dst, exe, mat):
    if mat is None:
        shutil.copy(src, dst)
    else:
        source = open(src, "r").readlines()
        destination = open(os.path.join(dst, 'insim.txt'), "w")
        b = False
        for line in range(len(source)):
            if b:
                b = False
                continue
            else:
                destination.write(source[line])
                if source[line] == "*status_gener \n":
                    destination.write("1\n")
                    b = True
                if source[line] == "*mating_design \n":
                    break
        for i in mat:
            for j in i:
                # écriture de la matrice de croisement ligne par ligne
                destination.write(str(int(j)) + '\n')
        destination.close()
    os.chdir(dst)
    os.startfile(exe)
    os.chdir('..')
