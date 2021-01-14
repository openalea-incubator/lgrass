import itertools
import numpy as np
import pandas as pd
import os

plantes_Bronsyn = pd.read_csv(
     'D:/Simon/Objectifs_modelisation/Reproductif_temps_eff/meilleures_simulations_V1.csv')
plantes_Bronsyn = plantes_Bronsyn[["C", "VR1", "LNM", "PPRM", "rho1", "L_ind", "kC", "kPremierCroiss", "PremierCroiss"]]
plantes_Bronsyn.rename(columns={'C': 'value_C'}, inplace=True)
plantes_Bronsyn.rename(columns={'VR1': 'daily_vern_rate'}, inplace=True)
plantes_Bronsyn.rename(columns={'LNM': 'leaf_number_max'}, inplace=True)
plantes_Bronsyn.rename(columns={'PPR': 'PPRM'}, inplace=True)
plantes_Bronsyn.rename(columns={'rho1': 'coeff_primordia_emission_reproductive'}, inplace=True)
plantes_Bronsyn.rename(columns={'L_ind': 'leaf_secondary_induction_coeff'}, inplace=True)
plantes_Bronsyn.rename(columns={'kC': 'increase_growth_C'}, inplace=True)
plantes_Bronsyn.rename(columns={'kPremierCroiss': 'increase_growth_Premiecroiss'}, inplace=True)
plantes_Bronsyn.rename(columns={'PremierCroiss': 'Premiecroiss'}, inplace=True)
plantes_Bronsyn['variete'] = 'Bronsyn'
plantes_Bronsyn = plantes_Bronsyn.drop_duplicates()

annees_a_simuler = pd.read_csv('D:/Simon/Objectifs_modelisation/GEVES/dates_epiaison_mesurees_GEVES.csv')
annees_a_simuler.annee = annees_a_simuler.annee.apply(str)
annees_a_simuler['sowing_date'] = annees_a_simuler.annee + '_09_15'
annees_a_simuler = annees_a_simuler[["site", "variete", "sowing_date"]]
plantes_Bronsyn = plantes_Bronsyn.merge(annees_a_simuler, on=["variete"])


temp_vern_min_list = [0]
temp_vern_inter_list = [6]
temp_vern_max_list = [13]
daily_vern_rate_list = [0.1]
basic_vern_rate_list = [0]
photoperiod_min_list = [10]
photoperiod_max_list = [16]
coeff_primordia_emission_vegetative_list = [1]
derivationLength_list = [5136]
option_tallage_list = ["False"]
option_senescence_list = ["False"]
option_floraison_list = ["True"]
cutting_height_list = [None]
cutting_end_list = [None]
PS_compensation_point_list = [14.1]
meteo_filename_list = ["Meteo_GEVES.csv"]
cutting_dates_list = ['']
leaf_primary_induction_coeff_list = [1]


simulation_plan = pd.DataFrame(columns=["temp_vern_min",
                                        "temp_vern_inter",
                                        "temp_vern_max",
                                        "daily_vern_rate",
                                        "basic_vern_rate",
                                        "photoperiod_min",
                                        "photoperiod_max",
                                        "coeff_primordia_emission_vegetative",
                                        "derivationLength",
                                        "option_tallage",
                                        "option_senescence",
                                        "option_floraison",
                                        "cutting_height",
                                        "cutting_end",
                                        "meteo_filename",
                                        "cutting_dates",
                                        "PS_compensation_point",
                                        "leaf_primary_induction_coeff",
                                        "value_C",
                                        "daily_vern_rate",
                                        "leaf_number_max",
                                        "PPRM",
                                        "coeff_primordia_emission_reproductive",
                                        "leaf_secondary_induction_coeff",
                                        "increase_growth_C",
                                        "increase_growth_Premiecroiss",
                                        "Premiecroiss",
                                        "variete",
                                        "site",
                                        "sowing_date"])

for i in itertools.product(temp_vern_min_list,
                           temp_vern_inter_list,
                           temp_vern_max_list,
                           daily_vern_rate_list,
                           basic_vern_rate_list,
                           photoperiod_min_list,
                           photoperiod_max_list,
                           coeff_primordia_emission_vegetative_list,
                           derivationLength_list,
                           option_tallage_list,
                           option_senescence_list,
                           option_floraison_list,
                           cutting_height_list,
                           cutting_end_list,
                           meteo_filename_list,
                           cutting_dates_list,
                           PS_compensation_point_list,
                           leaf_primary_induction_coeff_list):
    for r in [tuple(x) for x in plantes_Bronsyn.to_records(index=False)]:
        simulation_plan = simulation_plan.append(pd.DataFrame([i + r], columns=["temp_vern_min",
                                                                                "temp_vern_inter",
                                                                                "temp_vern_max",
                                                                                "daily_vern_rate",
                                                                                "basic_vern_rate",
                                                                                "photoperiod_min",
                                                                                "photoperiod_max",
                                                                                "coeff_primordia_emission_vegetative",
                                                                                "derivationLength",
                                                                                "option_tallage",
                                                                                "option_senescence",
                                                                                "option_floraison",
                                                                                "cutting_height",
                                                                                "cutting_end",
                                                                                "meteo_filename",
                                                                                "cutting_dates",
                                                                                "PS_compensation_point",
                                                                                "leaf_primary_induction_coeff",
                                                                                "value_C",
                                                                                "daily_vern_rate",
                                                                                "leaf_number_max",
                                                                                "PPRM",
                                                                                "coeff_primordia_emission_reproductive",
                                                                                "leaf_secondary_induction_coeff",
                                                                                "increase_growth_C",
                                                                                "increase_growth_Premiecroiss",
                                                                                "Premiecroiss",
                                                                                "variete",
                                                                                "site",
                                                                                "sowing_date"]))


simulation_plan['Scenario'] = [i for i in range(1, simulation_plan.shape[0] + 1)]
simulation_plan['name'] = 'simulation_' + simulation_plan['Scenario'].astype(str)

simulation_plan.to_csv(os.path.join('inputs', 'plan_simulation_Bronsyn_GEVES.csv'), index=False)
