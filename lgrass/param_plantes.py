# coding: utf8
import pandas as pd
import numpy as np
import math
from lgrass import flowering_functions


def get_param(in_param_file='inputs/Parametre_plante_Lgrass.xls', sheet='ParamP'):
    return pd.read_excel(in_param_file, sheet_name=sheet)


def define_param(in_genet_file='inputs/donnees_C.csv', out_param_file='inputs/Simulation_1.csv', param=None):
    value_C = pd.read_csv(in_genet_file)

    param_init = open(out_param_file, 'w')
    param_name = list(param['name'])
    param_value = list(param['value'])

    # Ecriture des parametres de plante variables
    Geno = [param_name[0]]
    C = [param_name[1]]
    for geno in range(len(value_C['geno'])):
        Geno.append(str(value_C['geno'].iloc[geno]))
        C.append(str(value_C['C'].iloc[geno]))
    param_init.write(";".join(Geno) + "\n")
    param_init.write(";".join(C) + "\n")

    # Ecriture des parametres de plante constants
    for par in range(2, len(param_name)):
        L = [param_name[par]]
        for geno in range(len(value_C['geno'])):
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
    Genotypes = np.array([i for i in value_C['geno']]).reshape(NBlignes, NBcolonnes)

    return ParamP, nb_plantes, NBlignes, NBcolonnes, posPlante, Plantes, Genotypes, flowering_param
