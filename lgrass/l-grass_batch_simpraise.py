# '''
# Created on 19/03/2020
#
# @author: modelisation - TR
# '''


# import the modules necessary to initiate the L-systems
import time
import re
from openalea.lpy import *
import math
from openalea.plantgl.all import *
#from alinea.caribu.CaribuScene import CaribuScene
import multiprocessing
import sys
import os
# from generateScene import *
import numpy as np
from lgrass import flowering_functions
from lgrass import meteo_ephem
import caribu
from openalea.lpy import Lsystem
import itertools
import pandas as pd
import random


TableParamP = pd.read_excel('inputs/Parametre_plante_Lgrass.xls', sheet_name='ParamP')
simul_conditions = pd.read_csv(os.path.join("inputs/plan_simulation.csv"))
value_C = pd.read_csv(os.path.join("inputs/donnees_C.csv"))
param_name = list(TableParamP['name'])
param_value = list(TableParamP['value'])
nb_graines_par_epillet = 3

# Cree la liste de L-systems et liste des noms
testsim = {}
names = []


def runlsystem(id_scenario = 0, opt_caribu=False, opt_reproduction=False, option_tiller_regression=False, option_mophogenetic_regulation_by_carbone=False):

    # Gestion de la recuperation des variables par le lsystem
    row = simul_conditions.iloc[id_scenario]
    name = str(row["name"])
    print(name)
    names.append(name)
    lpy_filename = os.path.join('lgrass.lpy')
    testsim[name] = Lsystem(lpy_filename)

    # Gestion des options
    testsim[name].option_tiller_regression = option_tiller_regression
    testsim[name].option_mophogenetic_regulation_by_carbone = option_mophogenetic_regulation_by_carbone

    # Ecriture du fichier des parametres de plante d'entree
    path_param = 'inputs/' + row[0] + '.csv'
    param_init = open(path_param, 'w')

    # Ecriture des parametres de plante variables
    Geno = []
    C = []
    Geno.append(param_name[0])
    C.append(param_name[1])
    for geno in range(len(value_C['geno'])):
        Geno.append(str(value_C['geno'].iloc[geno]))
        C.append(str(value_C['C'].iloc[geno]))
    param_init.write(";".join(Geno) + "\n")
    param_init.write(";".join(C) + "\n")

    # Ecriture des parametres de plante constants
    for par in range(2, len(param_name)):
        L = []
        L.append(param_name[par])
        for geno in range(len(value_C['geno'])):
            L.append(str(param_value[par]))
        param_init.write(";".join(L) + "\n")
    param_init.close()

    # lecture du fichier, creation de ParamP et des parametres de floraison
    param_plante = pd.read_csv(os.path.join(path_param), sep=";", header=None, index_col=0)

    flowering_param = flowering_functions.FloweringFunctions()
    flowering_param.param.__dict__.update((k, list(param_plante.loc[k,:])) for k in param_plante.index.intersection(flowering_param.param.__dict__))

    testsim[name].ParamP = list(dict(zip(param_plante.index, param_plante.iloc[:, col])) for col in range(len(param_plante.columns)))

    # Creation des matrices d'identifiant des plantes et de leur genotype
    testsim[name].nb_plantes = len(testsim[name].ParamP)
    testsim[name].NBlignes = int(math.ceil(np.sqrt(testsim[name].nb_plantes)))
    testsim[name].NBcolonnes = int(math.floor(np.sqrt(testsim[name].nb_plantes)))
    testsim[name].posPlante = [[i, j] for i, j in
                               zip(sorted(range(testsim[name].NBlignes) * testsim[name].NBcolonnes),
                                   range(testsim[name].NBcolonnes) * testsim[name].NBlignes)]
    testsim[name].Plantes = np.arange(testsim[name].nb_plantes).reshape(testsim[name].NBlignes,
                                                                        testsim[name].NBcolonnes)
    testsim[name].Genotypes = np.array([i for i in value_C['geno']]).reshape(testsim[name].NBlignes,
                                                                             testsim[name].NBcolonnes)
    # Donnees meteo
    sowing_date = "2018_11_13"
    site = 'treatment_1'
    testsim[name].meteo = meteo_ephem.import_meteo_data('inputs/meteo_controlled_conditions.csv', sowing_date, site)

    # Gestion caribu
    if opt_caribu:
        dico_caribu = {}
        dico_caribu['azimuths'] = 4
        dico_caribu['zeniths'] = 5
        dico_caribu['diffuse_model'] = 'soc'
        dico_caribu['scene_unit'] = 'mm'
        dico_caribu['RUE'] = 2
        dico_caribu['period_considered_tiller_regression'] = 10
        dico_caribu['radiation_threshold'] = 100000
        dico_caribu['Ray'] = [0.] * (testsim[name].nb_plantes)
        dico_caribu['radiation_interception'] = pd.DataFrame()
        dico_caribu['meteo'] = testsim[name].meteo
        dico_caribu['option_tiller_regression'] = testsim[name].option_tiller_regression
        dico_caribu['option_mophogenetic_regulation_by_carbone'] = testsim[
            name].option_mophogenetic_regulation_by_carbone

    # Parametres de simulation
    testsim[name].derivationLength = int(row["derivationLength"])
    testsim[name].option_tallage = row["option_tallage"]
    testsim[name].meteo_path = os.path.join(row["meteo_path"])
    testsim[name].sowing_date = row["sowing_date"]
    testsim[name].site = row["site"]
    testsim[name].flowering_model = flowering_param
    testsim[name].output_induction_file_name = name + '_' + 'induction'
    testsim[name].output_organ_lengths_file_name = name + '_' + 'organ_lengths'
    testsim[name].cutting_dates = [] if pd.isna(row["cutting_dates"]) \
        else [row["cutting_dates"]] if isinstance(row["cutting_dates"], int) \
        else [int(i) for i in row["cutting_dates"].split("_")]

    # Gestion de la creation d une nouvelle generation
    if opt_reproduction:
        Epillets = []

    # Lancement du lsystem
    lstring = testsim[name].axiom
    for dd in range(testsim[name].derivationLength):
        try:
            day = testsim[name].current_day
        except:
            day = 1
        lstring = testsim[name].derive(lstring, dd, 1)
        lscene = testsim[name].sceneInterpretation(lstring)
        if opt_caribu:
            testsim[name].BiomProd, dico_caribu['radiation_interception'], dico_caribu[
                'Ray'] = caribu.apply_caribu_lgrass(lstring, lscene, testsim[name].TPS,
                                                    testsim[name].current_day,
                                                    testsim[name].tiller_appearance,
                                                    testsim[name].nb_plantes, dico_caribu, day)
        Viewer.display(lscene)

    # Gestion des graines
    for mod in lstring:
        if mod.name in ('apex',):
            if mod[0].phenological_state == 'reproductive' and mod[0].final_spikelet_number is not None:
                Epillets.append((mod[0].id_plante, mod[0].id_talle, mod[0].final_spikelet_number))
    nb_tot_graines = []
    for tupple in Epillets:
        for i in range(int(tupple[2]) * nb_graines_par_epillet):
            nb_tot_graines.append(tupple)
    Graines = []
    for rand in range(testsim[name].nb_plantes):
        if len(nb_tot_graines) != 0:
            Graines.append(nb_tot_graines[random.randint(0, len(nb_tot_graines) - 1)])


    testsim[name].clear()
    print(''.join((name, " - done")))

#for x in range(simul_conditions.shape[0]):
runlsystem(id_scenario=0,opt_caribu=True, opt_reproduction=True, option_tiller_regression=True, option_mophogenetic_regulation_by_carbone=False)
