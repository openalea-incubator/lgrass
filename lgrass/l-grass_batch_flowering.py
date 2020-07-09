import os
from lgrass import flowering_functions
from openalea.lpy import Lsystem
import pandas as pd


INPUTS_DIRPATH = 'inputs'
OUTPUTS_DIRPATH = r'outputs'

# cree la liste de L-systems et liste des noms
testsim = {}
names = []

simul_conditions = pd.read_csv(os.path.join("inputs/plan_simulation.csv"))

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
    testsim[name].derivationLength = int(row["derivationLength"])
    testsim[name].option_tallage = row["option_tallage"]
    testsim[name].meteo_path = os.path.join(row["meteo_filename"])
    testsim[name].sowing_date = row["sowing_date"]
    testsim[name].site = row["site"]
    testsim[name].flowering_model = flowering_param
    testsim[name].output_induction_file_name = name + '_' + 'induction'
    testsim[name].output_organ_lengths_file_name = name + '_' + 'organ_lengths'
    testsim[name].cutting_dates = [] if pd.isna(row["cutting_dates"]) \
        else [row["cutting_dates"]] if isinstance(row["cutting_dates"], int) \
        else [int(i) for i in row["cutting_dates"].split("_")]
    testsim[name].ParamP[0]["C"] = row["value_C"]
    testsim[name].ParamP[0]["Premiecroiss"] = row["Premiecroiss"]


# function to run an L-system from the 'testsim' dictionnary
def runlsystem(n):
    testsim[names[n]].derive()  # permet le declenchement de la fonction "End" du script lpy
    testsim[names[n]].clear()
    print(''.join((names[n], " - done")))


for i in range(simul_conditions.shape[0]):
    try:
        runlsystem(i)
    except:
        print(simul_conditions.iloc[i]["name"] + " " + "failed")