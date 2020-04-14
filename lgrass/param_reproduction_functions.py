# coding: utf8
import random
import numpy as np
import pandas as pd
import math
import shutil
import os
from lgrass import flowering_functions


# Créer la matrice de croisement des plantes pour établir une nouvelle génération de nb_plantes
def create_seeds(lstring, param, nb_plantes):
    seed_number = int(param.loc[param['name'] == 'seeds_by_spikelet', :]['value'])
    matrix = np.zeros((nb_plantes, nb_plantes))
    seeds = []
    elected_seeds = []

    # construction des graines et de leurs parents
    for mod in lstring:
        if mod.name in ('apex',):
            if mod[0].final_spikelet_number is not None:
                for i in range(int(mod[0].final_spikelet_number * seed_number)):
                    id_mere = mod[0].id_plante  # séléction de la mère
                    peres = [j for j in range(nb_plantes)]
                    peres.remove(id_mere)
                    id_pere = random.Random().choice(peres)  # séléction aléatoire du père
                    seeds.append((id_mere, id_pere))
    print ("seeds : ", seeds)
    # sélection aléatoire des graines pour la génération suivante et création de la matrice de croisement
    if len(seeds) < nb_plantes:
        raise NameError("Il n'y a pas eu suffisamment de graines produites pour établir une nouvelle génération.")

    for rand in range(nb_plantes):
        seed = random.Random().choice(seeds)
        print("tirage :", seed)
        elected_seeds.append(seed)
        seeds.remove(seed)
        matrix[elected_seeds[rand][0], elected_seeds[rand][1]] += 1
    print(matrix)
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


# Récupérer les paramètres de plante constants
def get_param(in_param_file=None, sheet='ParamP'):
    if in_param_file is None:
        in_param_file = 'inputs/Parametre_plante_Lgrass.xls'
    return pd.read_excel(in_param_file, sheet_name=sheet)


# La formule pour convertir la donnée génétique en paramètre C
def calculate_C(n):

    ### A traiter ###

    return str(1.5)


# Création du fichier de paramètres d'entrée pour chaque plante et configuration du planteur lgrass
def define_param(in_genet_file=None, out_param_file=None, param=None, id_gener=1, opt_repro=False):
    if in_genet_file is None:
        in_genet_file = 'inputs/donnees_C.csv'
    if out_param_file is None:
        out_param_file = 'inputs/Simulation_G' + str(id_gener) + '.csv'
    if param is None:
        param = get_param()
    if opt_repro:
        infile = get_genet_file(in_genet_file=in_genet_file)
        data = infile.loc[infile['D'] == str(id_gener), :]
    else:
        data = pd.read_csv(in_genet_file)

    param_init = open(out_param_file, 'w')
    param_name = list(param['name'])
    param_value = list(param['value'])

    Geno = [param_name[0]]
    C = [param_name[1]]
    for i in range(len(data)):
        Geno.append(str(data['geno'].iloc[i]))
        if opt_repro:
            C.append(calculate_C(float(data['C'].iloc[i])))
        else:
            C.append(str(data['C'].iloc[i]))
    param_init.write(";".join(Geno) + "\n")
    param_init.write(";".join(C) + "\n")

    for par in range(2, len(param_name)):
        L = [param_name[par]]
        for i in range(len(data)):
            L.append(str(param_value[par]))
        param_init.write(";".join(L) + "\n")
    param_init.close()

    # lecture du fichier, creation de ParamP et des parametres de floraison
    param_plante = pd.read_csv(out_param_file, sep=";", header=None, index_col=0)

    flowering_param = flowering_functions.FloweringFunctions()
    flowering_param.param.__dict__.update(
        (k, list(param_plante.loc[k, :])) for k in param_plante.index.intersection(flowering_param.param.__dict__))

    ParamP = list(dict(zip(param_plante.index, param_plante.iloc[:, col])) for col in range(len(param_plante.columns)))

    # Creation des matrices d'identifiant des plantes et de leur genotype
    nb_plantes = len(ParamP)
    NBlignes = int(math.ceil(np.sqrt(nb_plantes)))
    NBcolonnes = int(math.floor(np.sqrt(nb_plantes)))
    posPlante = [[i, j] for i, j in zip(sorted(range(NBlignes) * NBcolonnes), range(NBcolonnes) * NBlignes)]
    Plantes = np.arange(nb_plantes).reshape(NBlignes, NBcolonnes)
    Genotypes = np.array([i for i in data['geno']]).reshape(NBlignes, NBcolonnes)
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
                destination.write(str(int(j)) + '\n')
        destination.close()
    os.chdir(dst)
    os.startfile(exe)
    os.chdir('..')
