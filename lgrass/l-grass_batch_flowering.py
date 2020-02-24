'''
Created on 11/10/2018

@author: modelisation - SR
'''
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
from openalea.lpy import Lsystem
import itertools
import pandas as pd

#path_ = r'D:\Simon\Python\lgrass\lgrass'

INPUTS_DIRPATH = 'inputs'
OUTPUTS_DIRPATH = r'outputs'
# date_string = datetime.datetime.now() .strftime('%Y_%m_%d_%Hh%M')
# output_batch_path = open(os.path.join(OUTPUTS_DIRPATH, 'output_batch_path', date_string), "wb")
# output_batch = csv.writer(output_batch_path)
# output_batch.writerow(['vai', 'vbee', 'sldl', 'absolute_max_leaf_number', 'absolute_min_leaf_number'])


# path_param = os.path.join(INPUTS_DIRPATH,'Parametre_plante_Lgrass.xls') # Fichier contenant les parametres
#path_param = 'D:\Simon\Python\lgrass\lgrass\inputs\Parametre_plante_Lgrass.xls'
path_param = 'inputs/Parametre_plante_Lgrass.xls'

onglet1 = 'FL'
onglet2 = 'FC'
TableParamP1 = pd.read_excel(path_param, sheetname=onglet1)
TableParamP2 = pd.read_excel(path_param, sheetname=onglet2)
paramP1 = dict(zip(TableParamP1['name'], TableParamP1['value']))
paramP2 = dict(zip(TableParamP2['name'], TableParamP2['value']))

# cree la liste de L-systems et liste des noms
testsim = {}
names = []

nb_simul = 0
temp_vern_min_list = np.arange(0, 1, 10)
temp_vern_inter_list = np.arange(8, 9, 10)
temp_vern_max_list = np.arange(13, 18, 10)
daily_vern_rate_list = np.arange(0.001, 0.002, 10)
basic_vern_rate_list = np.arange(0.01, 0.11, 10)
photoperiod_min_list = np.arange(10, 11, 10)
photoperiod_max_list = np.arange(16, 17, 10)
max_photo_ind_rate_list = np.arange(1, 2, 10)
coeff_primordia_emission_vegetative_list = np.arange(1, 2, 10)
coeff_primordia_emission_reproductive_list = np.arange(1, 6, 10)
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
#value_C = np.array([0.7, 1.2])
Premiecroiss = np.array([90, 110, 150, 170])

# column_names = ["name", "temp_vern_min", "temp_vern_inter", "temp_vern_max", "daily_vern_rate",
#                 "basic_vern_rate", "photoperiod_min", "photoperiod_max", "max_photo_ind_rate",
#                 "coeff_primordia_emission_vegetative", "coeff_primordia_emission_reproductive",
#                 "derivationLength", "option_tallage", "meteo_path", "sowing_date", "site",
#                 "output_induction_file_name", "cutting_dates", "value_C", "Premiecroiss]

# simul_conditions = pd.DataFrame(columns=column_names)

#simul_conditions = pd.read_csv(r"D:\Simon\Python\lgrass\lgrass\inputs\plan_simulation.csv")
simul_conditions = pd.read_csv(os.path.join("inputs/plan_simulation.csv"))
value_C = pd.read_csv(os.path.join("inputs/donnees_C.csv"))

for x in range(simul_conditions.shape[0]):
    row = simul_conditions.iloc[x]
    print(row)

    flowering_param = flowering_functions.FloweringFunctions()
    flowering_param.param.temp_vern_min = row["temp_vern_min"]
    flowering_param.param.temp_vern_inter = row["temp_vern_inter"]
    flowering_param.param.temp_vern_max = row["temp_vern_max"]
    flowering_param.param.daily_vern_rate = row["daily_vern_rate"]
    flowering_param.param.basic_vern_rate = row["basic_vern_rate"]
    flowering_param.param.photoperiod_min = row["photoperiod_min"]
    flowering_param.param.photoperiod_max = row["photoperiod_max"]
    flowering_param.param.max_photo_ind_rate = row["max_photo_ind_rate"]
    flowering_param.param.coeff_primordia_emission_vegetative = row["coeff_primordia_emission_vegetative"]
    flowering_param.param.coeff_primordia_emission_reproductive = row["coeff_primordia_emission_reproductive"]

    name = str(row["name"])
    print(name)
    names.append(name)
    lpy_filename = os.path.join('lgrass.lpy')
    testsim[name] = Lsystem(lpy_filename)
    #testsim[name].axiom = "CouvertVegetal(2)"
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
    testsim[name].ParamP = []
    testsim[name].ParamP.append(paramP1)
    testsim[name].ParamP.append(paramP2)
    testsim[name].nb_plantes = len(testsim[name].ParamP)
    testsim[name].NBlignes = int(math.ceil(np.sqrt(testsim[name].nb_plantes)))
    testsim[name].NBcolonnes = int(math.floor(np.sqrt(testsim[name].nb_plantes)))
    testsim[name].posPlante = [[i, j] for i, j in zip(sorted(range(testsim[name].NBlignes) * testsim[name].NBcolonnes), range(testsim[name].NBcolonnes) * testsim[name].NBlignes)]
    # Creation des matrices d'identifiant des plantes et de leur genotype
    testsim[name].Plantes = np.arange(testsim[name].nb_plantes).reshape(testsim[name].NBlignes, testsim[name].NBcolonnes)
    testsim[name].Genotypes = np.array([i for i in value_C['geno']]).reshape(testsim[name].NBlignes, testsim[name].NBcolonnes)
    for i in range(len(value_C["geno"])):
        testsim[name].ParamP[i]["genotype"] = value_C['geno'].iloc[i]
        testsim[name].ParamP[i]["C"] = value_C['C'].iloc[i]
        testsim[name].ParamP[i]["Premiecroiss"] = row["Premiecroiss"]


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

# function to run an L-system from the 'testsim' dictionnary
def runlsystem(n):
    testsim[names[n]].derive()  # permet le declenchement de la fonction "End" du script lpy
    # print(testsim[names[n]].output_dict)
    # with open(os.path.join(OUTPUTS_DIRPATH, 'sortie_test_path', str(n) + '.csv'), 'wb') as sortie_test_path:
    #     sortie_test = csv.writer(sortie_test_path)
    #     sortie_test.writerows(testsim[names[n]].output_dict.items())
    testsim[names[n]].clear()
    print(''.join((names[n], " - done")))

    # print(testsim[names[n]].output_dict)


for i in range(simul_conditions.shape[0]):
    #try:
    runlsystem(i)
    #except:
    #    print(simul_conditions.iloc[i]["name"] + " " + "failed")


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
