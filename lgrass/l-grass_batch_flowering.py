# '''
# Created on 11/10/2018
#
# @author: modelisation - SR
# '''
# batch pour L-grass with flowering

# import the modules necessary to initiate the L-systems
import time

import re
from openalea.lpy import *
import math
from openalea.plantgl.all import *
from alinea.caribu.CaribuScene import CaribuScene
import multiprocessing
import sys
import os
# from generateScene import *
import numpy as np
from lgrass import flowering_functions
from lgrass import meteo_ephem
from lgrass import caribu
from openalea.lpy import Lsystem
import itertools
import pandas as pd


# path_ = r'D:\Simon\Python\lgrass\lgrass'

INPUTS_DIRPATH = 'inputs'
OUTPUTS_DIRPATH = r'outputs'
# date_string = datetime.datetime.now() .strftime('%Y_%m_%d_%Hh%M')
# output_batch_path = open(os.path.join(OUTPUTS_DIRPATH, 'output_batch_path', date_string), "wb")
# output_batch = csv.writer(output_batch_path)
# output_batch.writerow(['vai', 'vbee', 'sldl', 'absolute_max_leaf_number', 'absolute_min_leaf_number'])


# path_param = os.path.join(INPUTS_DIRPATH,'Parametre_plante_Lgrass.xls') # Fichier contenant les parametres
# path_param = 'D:\Simon\Python\lgrass\lgrass\inputs\Parametre_plante_Lgrass.xls'

TableParamP = pd.read_excel('inputs/Parametre_plante_Lgrass.xls', sheet_name='ParamP')
simul_conditions = pd.read_csv(os.path.join("inputs/plan_simulation.csv"))
value_C = pd.read_csv(os.path.join("inputs/donnees_C.csv"))
param_name = list(TableParamP['name'])
param_value = list(TableParamP['value'])

# dict(zip(TableParamP['name'], TableParamP['value']))

# cree la liste de L-systems et liste des noms
testsim = {}
names = []


combination_site_sowing_date = np.array([("treatment_1", "2018_11_13")])
# combination_site_sowing_date = np.array([("LUSIGNAN", "2017_10_15")])
# combination_site_sowing_date = np.array([("LUSIGNAN", "2012_10_15"),
#                                          ("LUSIGNAN", "2017_10_15"),
#                                          ("THEIX", "2007_10_15"),
#                                          ("THEIX", "2012_10_15"),
#                                          ("THEIX", "2017_10_15"),
#                                          ("PLOUDANIEL", "2007_10_15"),
#                                          ("PLOUDANIEL", "2012_10_15"),
#                                          ("PLOUDANIEL", "2017_10_15")])
# value_C = np.array([0.7, 1.2])
Premiecroiss = np.array([90, 110, 150, 170])


# column_names = ["name", "temp_vern_min", "temp_vern_inter", "temp_vern_max", "daily_vern_rate",
#                 "basic_vern_rate", "photoperiod_min", "photoperiod_max", "max_photo_ind_rate",
#                 "coeff_primordia_emission_vegetative", "coeff_primordia_emission_reproductive",
#                 "derivationLength", "option_tallage", "meteo_path", "sowing_date", "site",
#                 "output_induction_file_name", "cutting_dates", "value_C", "Premiecroiss]

# simul_conditions = pd.DataFrame(columns=column_names)

# simul_conditions = pd.read_csv(r"D:\Simon\Python\lgrass\lgrass\inputs\plan_simulation.csv")

def runlsystem(opt_caribu):
    for x in range(simul_conditions.shape[0]):
        row = simul_conditions.iloc[x]
        name = str(row["name"])
        print(name)
        names.append(name)
        lpy_filename = os.path.join('lgrass.lpy')
        testsim[name] = Lsystem(lpy_filename)
        # ecriture du fichier des parametres d'entree
        path_param = 'inputs/' + row[0] + '.csv'
        param_init = open(path_param, 'w')
        # ecriture des parametres variables
        Geno = []
        C = []
        Geno.append(param_name[0])
        C.append(param_name[1])
        for geno in range(len(value_C['geno'])):
            Geno.append(str(value_C['geno'].iloc[geno]))
            C.append(str(value_C['C'].iloc[geno]))
        param_init.write(";".join(Geno) + "\n")
        param_init.write(";".join(C) + "\n")
        # ecriture des parametres constants
        for par in range(2, len(param_name)):
            L = []
            L.append(param_name[par])
            for geno in range(len(value_C['geno'])):
                L.append(str(param_value[par]))
            param_init.write(";".join(L) + "\n")
        param_init.close()
        # lecture du fichier, creation de ParamP et des parametres de floraison
        param_plante = pd.read_csv(os.path.join(path_param), sep=";", header=None)
        testsim[name].ParamP = []
        flowering_param = flowering_functions.FloweringFunctions()
        for par in range(1, len(param_plante.columns)):
            L = dict(zip(param_plante.iloc[:, 0], param_plante.iloc[:, par]))
            testsim[name].ParamP.append(L)
            flowering_param.param.temp_vern_min.append(L["temp_vern_min"])
            flowering_param.param.temp_vern_inter.append(L["temp_vern_inter"])
            flowering_param.param.temp_vern_max.append(L["temp_vern_max"])
            flowering_param.param.daily_vern_rate.append(L["daily_vern_rate"])
            flowering_param.param.basic_vern_rate.append(L["basic_vern_rate"])
            flowering_param.param.photoperiod_min.append(L["photoperiod_min"])
            flowering_param.param.photoperiod_max.append(L["photoperiod_max"])
            flowering_param.param.max_photo_ind_rate.append(L["max_photo_ind_rate"])
            flowering_param.param.coeff_primordia_emission_vegetative.append(L["coeff_primordia_emission_vegetative"])
            flowering_param.param.coeff_primordia_emission_reproductive.append(L["coeff_primordia_emission_reproductive"])
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
        testsim[name].meteo = meteo_ephem.import_meteo_data('inputs/meteo_controlled_conditions.csv', sowing_date,site)
        # Gestion caribu
        if (opt_caribu == 'On'):
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
            dico_caribu['option_mophogenetic_regulation_by_carbone'] = testsim[name].option_mophogenetic_regulation_by_carbone
        # Parametres de simulation
        # testsim[name].axiom = "CouvertVegetal(2)"
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
        # Lancement du lsystem
        lstring = testsim[name].axiom
        for dd in range(testsim[name].derivationLength):
            try :
                day = testsim[name].current_day
            except :
                day = 1
            lstring = testsim[name].derive(lstring, dd, 1)
            lscene = testsim[name].sceneInterpretation(lstring)
            if (opt_caribu == 'On'):
                testsim[name].BiomProd, dico_caribu['radiation_interception'], dico_caribu['Ray'] = caribu.apply_caribu_lgrass(lstring, lscene, testsim[name].TPS,
                                                                    testsim[name].current_day,
                                                                    testsim[name].tiller_appearance,
                                                                    testsim[name].nb_plantes, dico_caribu, day)
            Viewer.display(lscene)

        # testsim[names[n]].derive() # permet le declenchement de la fonction "End" du script lpy
        # testsim[names[x]].animate()
        # print(testsim[names[n]].output_dict)
        # with open(os.path.join(OUTPUTS_DIRPATH, 'sortie_test_path', str(n) + '.csv'), 'wb') as sortie_test_path:
        #     sortie_test = csv.writer(sortie_test_path)
        #     sortie_test.writerows(testsim[names[n]].output_dict.items())
        testsim[names[x]].clear()
        print(''.join((names[x], " - done")))

        # print(testsim[names[n]].output_dict)

# try:
runlsystem('On')
# except:
#    print(simul_conditions.iloc[i]["name"] + " " + "failed")


# for x in itertools.product(temp_vern_min_list, temp_vern_inter_list, temp_vern_max_list, daily_vern_rate_list,
#                            basic_vern_rate_list, photoperiod_min_list, photoperiod_max_list, max_photo_ind_rate_list,
#                            coeff_primordia_emission_vegetative_list, coeff_primordia_emission_reproductive_list,
#                            combination_site_sowing_date, value_C, Premiecroiss):
#     print(x)
#
#     nb_simul += 1
#     flowering_param = flowering_functions.FloweringFunctions()
#     flowering_param.param.temp_vern_min = x[0]
#     flowering_param.param.temp_vern_inter = x[1]
#     flowering_param.param.temp_vern_max = x[2]
#     flowering_param.param.daily_vern_rate = x[3]
#     flowering_param.param.basic_vern_rate = x[4]
#     flowering_param.param.photoperiod_min = x[5]
#     flowering_param.param.photoperiod_max = x[6]
#     flowering_param.param.max_photo_ind_rate = x[7]
#     flowering_param.param.coeff_primordia_emission_vegetative = x[8]
#     flowering_param.param.coeff_primordia_emission_reproductive = x[9]
#
#     name = str(x[10][0]) + "_" + str(x[10][1]) + "_C_" + str(x[11]) + "_" + str(x[12])
#     print(name)
#     names.append(name)
#     lpy_filename = os.path.join('lgrass.lpy')
#     testsim[name] = Lsystem(lpy_filename)
#     testsim[name].derivationLength = 3000
#     testsim[name].option_tallage = False
#     testsim[name].meteo_path = os.path.join('inputs/meteo_controlled_conditions_bloc_1.csv')
#     testsim[name].sowing_date = x[10][1]
#     testsim[name].site = x[10][0]
#     testsim[name].flowering_model = flowering_param
#     testsim[name].output_induction_file_name = 'induction_' + name
#     testsim[name].output_organ_lengths_file_name = 'organ_lengths_' + name
#     testsim[name].cutting_dates = []
#     testsim[name].ParamP[0]["C"] = x[11]
#     testsim[name].ParamP[0]["Premiecroiss"] = x[12]

#    simul_conditions = simul_conditions.append(pd.DataFrame([a, ], columns=column_names))


# for i in range(nb_simul):
#     runlsystem(i)

# #run the L-systems
# if __name__ == '__main__':
#     multiprocessing.freeze_support()
#     CPUnb = multiprocessing.cpu_count()-1 #nombre de processeurs, moins un par prudence. (et pour pouvoir faire d'autres choses en meme temps)
#     print 'nb CPU: ' + str(CPUnb)
#     pool = multiprocessing.Pool(processes=CPUnb)
#     for i in range(int(nb_simul)):
#         pool.apply_async(runlsystem, args=(i,)) #Lance CPUnb simulations en meme temps, lorsqu'une simulation se termine elle est immediatement remplacee par la suivante
#     pool.close()
#     pool.join()
