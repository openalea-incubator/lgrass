# coding: utf8
class ParametersValues:
    def __init__(self):
        self.temp_vern_min = 0  # minimal vernalisation temperature (°C)
        self.temp_vern_inter = 8  # intermediate vernalisation temperature (°C)
        self.temp_vern_max = 17  # maximal vernalisation temperature (°C)
        self.daily_vern_rate = 0.001  # daily vernalisation rate (°C^-1)
        self.basic_vern_rate = 0.01  # basic vernalisation rate
        self.photoperiod_min = 10  # minimal efficient photoperiod (h)
        self.photoperiod_max = 16  # maximal efficient photoperiod (h)
        self.max_photo_ind_rate = 1  # photoperiodic induction ratedaily_vern_rate
        # -----------------------------------------------------
        self.coeff_primordia_emission_vegetative = 1  # number of primordia producted when one new leaf appear on a vegetative tiller
        self.coeff_primordia_emission_reproductive = 2  # number of primordia producted when one new leaf appear on a vegetative tiller
