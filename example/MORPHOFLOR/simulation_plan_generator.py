import itertools
import numpy as np
import pandas as pd
import os


temp_vern_min_list = [0]
temp_vern_inter_list = [6]
temp_vern_max_list = [13]
daily_vern_rate_list = [0.00135]
basic_vern_rate_list = [0]
photoperiod_min_list = [10]
photoperiod_max_list = [16]
max_photo_ind_rate_list = [1]
coeff_primordia_emission_vegetative_list = [1]
coeff_primordia_emission_reproductive_list = [2]
derivationLength_list = [4400]
option_tallage_list = ["False"]
option_senescence_list = ["False"]
option_floraison_list = ["True"]
meteo_filename_list = ["environnement_morphoflor.csv"]
sowing_date_list = ["2019_01_16"]
site_list = ["treatment_" + "4" + "_all_blocks"]
value_C_list = np.round(np.arange(7.25, 8, 0.25), decimals=2)
Premiecroiss_list = np.round(np.arange(45, 75, 5), decimals=2)
PS_compensation_point_list = [14.1]
leaf_primary_induction_coeff_list = [1]
leaf_secondary_induction_coeff_list = [1]
increase_growth_Premiecroiss_list = [1]
increase_growth_C_list = [1]


simulation_plan = pd.DataFrame(columns=["temp_vern_min",
                           "temp_vern_inter",
                           "temp_vern_max",
                           "daily_vern_rate",
                           "basic_vern_rate",
                           "photoperiod_min",
                           "photoperiod_max",
                           "max_photo_ind_rate",
                           "coeff_primordia_emission_vegetative",
                           "coeff_primordia_emission_reproductive",
                           "derivationLength",
                           "option_tallage",
                           "option_senescence",
                           "option_floraison",
                           "meteo_filename",
                           "sowing_date",
                           "site",
                           "value_C",
                           "Premiecroiss",
                           "PS_compensation_point",
                           "leaf_primary_induction_coeff",
                           "leaf_secondary_induction_coeff",
                           "increase_growth_Premiecroiss",
                           "increase_growth_C"])

for i in itertools.product(temp_vern_min_list,
                           temp_vern_inter_list,
                           temp_vern_max_list,
                           daily_vern_rate_list,
                           basic_vern_rate_list,
                           photoperiod_min_list,
                           photoperiod_max_list,
                           max_photo_ind_rate_list,
                           coeff_primordia_emission_vegetative_list,
                           coeff_primordia_emission_reproductive_list,
                           derivationLength_list,
                           option_tallage_list,
                           option_senescence_list,
                           option_floraison_list,
                           meteo_filename_list,
                           sowing_date_list,
                           site_list,
                           value_C_list,
                           Premiecroiss_list,
                           PS_compensation_point_list,
                           leaf_primary_induction_coeff_list,
                           leaf_secondary_induction_coeff_list,
                           increase_growth_Premiecroiss_list,
                           increase_growth_C_list):
    simulation_plan = simulation_plan.append(pd.DataFrame([i], columns=["temp_vern_min",
                           "temp_vern_inter",
                           "temp_vern_max",
                           "daily_vern_rate",
                           "basic_vern_rate",
                           "photoperiod_min",
                           "photoperiod_max",
                           "max_photo_ind_rate",
                           "coeff_primordia_emission_vegetative",
                           "coeff_primordia_emission_reproductive",
                           "derivationLength",
                           "option_tallage",
                           "option_senescence",
                           "option_floraison",
                           "meteo_filename",
                           "sowing_date",
                           "site",
                           "value_C",
                           "Premiecroiss",
                           "PS_compensation_point",
                           "leaf_primary_induction_coeff",
                           "leaf_secondary_induction_coeff",
                           "increase_growth_Premiecroiss",
                           "increase_growth_C"]))


id_first_simulation = 314
simulation_plan['Scenario'] = [i for i in range(id_first_simulation, simulation_plan.shape[0] + id_first_simulation)]
simulation_plan['name'] = 'simulation_' + simulation_plan['Scenario'].astype(str)

simulation_plan.to_csv(os.path.join('inputs', 'plan_simulation_veg.csv'), index=False)