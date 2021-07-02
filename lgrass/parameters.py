# coding: utf8
class ParametersValues:
    def __init__(self):
        self.temp_vern_min = 0  # minimal vernalisation temperature (°C)
        self.temp_vern_inter = 8  # intermediate vernalisation temperature (°C)
        self.temp_vern_max = 17  # maximal vernalisation temperature (°C)
        self.daily_vern_increment = 0.001  # daily vernalisation increment (dimensionless)
        self.photoperiod_min = 10.  # minimal efficient photoperiod (h)
        self.photoperiod_max = 16.  # maximal efficient photoperiod (h)
        self.leaf_number_max = 14  # number of leaf produced before the secondary induction which maximizes secondary induction speed
        self.PPRM = 10  # maximal daily increment of the secondary induction (dimensionless)
        # -----------------------------------------------------
        self.coeff_primordia_emission_vegetative = 1  # number of primordia producted when one new leaf appear on a vegetative tiller
        self.coeff_primordia_emission_reproductive = 2  # number of primordia producted when one new leaf appear on a vegetative tiller
        # -----------------------------------------------------
        self.leaf_primary_induction_coeff = 1
        self.leaf_secondary_induction_coeff = 1
        self.increase_growth_Premiecroiss = 1
        self.increase_growth_C = 1
