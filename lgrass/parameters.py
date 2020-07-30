# coding: utf8
class ParametersValues:
    def __init__(self):
        self.temp_vern_min = []  # minimal vernalisation temperature (°C)
        self.temp_vern_inter = []  # intermediate vernalisation temperature (°C)
        self.temp_vern_max = []  # maximal vernalisation temperature (°C)
        self.daily_vern_rate = []  # daily vernalisation rate (°C^-1)
        self.basic_vern_rate = []  # basic vernalisation rate
        self.photoperiod_min = []  # minimal efficient photoperiod (h)
        self.photoperiod_max = []  # maximal efficient photoperiod (h)
        self.max_photo_ind_rate = []  # photoperiodic induction ratedaily_vern_rate
        # -----------------------------------------------------
        self.coeff_primordia_emission_vegetative = 1  # number of primordia producted when one new leaf appear on a vegetative tiller
        self.coeff_primordia_emission_reproductive = 2  # number of primordia producted when one new leaf appear on a vegetative tiller
        # -----------------------------------------------------
        self.leaf_primary_induction_coeff = 1
        self.leaf_secondary_induction_coeff = 1
        self.increase_growth_Premiecroiss = 1
        self.increase_growth_C = 1
