import itertools
import numpy as np
import pandas as pd


temp_vern_min_list = [0]
temp_vern_inter_list = [6]
temp_vern_max_list = [13]
daily_vern_rate_list = [0.00135]
basic_vern_rate_list = [0.01]
photoperiod_min_list = [10]
photoperiod_max_list = [16]
max_photo_ind_rate_list = [1]
coeff_primordia_emission_vegetative_list = [1]
coeff_primordia_emission_reproductive_list = [2]
derivationLength_list = [3000]
option_tallage_list = ["False"]
meteo_filename_list = ["environnement_morphoflor.csv"]
sowing_date_list = ["2019_01_16"]
treatment_bloc = list(itertools.product([4], [1, 2, 3, 4, 5, 6, 7, 8, 9]))
site_list = ["treatment_" + str(n[0]) + "_bloc_" + str(n[1]) for n in treatment_bloc]
cutting_dates_list = [np.nan]
value_C_list = np.arange(0.8, 1.75, 0.05)
Premiecroiss_list = [90]
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
                           "meteo_filename",
                           "sowing_date",
                           "site",
                           "cutting_dates",
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
                           meteo_filename_list,
                           sowing_date_list,
                           site_list,
                           cutting_dates_list,
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
                           "meteo_filename",
                           "sowing_date",
                           "site",
                           "cutting_dates",
                           "value_C",
                           "Premiecroiss",
                           "PS_compensation_point",
                           "leaf_primary_induction_coeff",
                           "leaf_secondary_induction_coeff",
                           "increase_growth_Premiecroiss",
                           "increase_growth_C"]))

simulation_plan['name'] = 'simulation_' + 'C_' + simulation_plan['value_C'].astype(str) + '_' + simulation_plan['site']
simulation_plan['Scenario'] = [i for i in range(simulation_plan.shape[0])]

simulation_plan.to_csv('plan_simulation.csv', index=False)