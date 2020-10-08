# This file contains all fonctions associated with flowering
from lgrass import parameters


class FloweringFunctions:
    def __init__(self):
        self.param = parameters.ParametersValues()

    def vernalisation_function(self, temperature,i):
        """

        :param temperature:
        :return:
        """
        if isinstance(self.param.temp_vern_min, list):
            if self.param.temp_vern_min[i] <= temperature <= self.param.temp_vern_inter[i]:
                primary_induction_increment = self.param.daily_vern_rate[i] * temperature + self.param.basic_vern_rate[i]
            elif self.param.temp_vern_inter[i] < temperature <= self.param.temp_vern_max[i]:
                #primary_induction_increment = (self.param.daily_vern_rate * temperature + self.param.basic_vern_rate) * (self.param.temp_vern_max - temperature) / (self.param.temp_vern_max - self.param.temp_vern_inter)
                primary_induction_increment = (-temperature + self.param.temp_vern_max[i]) * (self.param.daily_vern_rate[i] * self.param.temp_vern_inter[i] + self.param.basic_vern_rate[i]) / (self.param.temp_vern_max[i] - self.param.temp_vern_inter[i])
            else:
                primary_induction_increment = 0
            return primary_induction_increment
        else:
            if self.param.temp_vern_min <= temperature <= self.param.temp_vern_inter:
                primary_induction_increment = self.param.daily_vern_rate * temperature + self.param.basic_vern_rate
            elif self.param.temp_vern_inter < temperature <= self.param.temp_vern_max:
                #primary_induction_increment = (self.param.daily_vern_rate * temperature + self.param.basic_vern_rate) * (self.param.temp_vern_max - temperature) / (self.param.temp_vern_max - self.param.temp_vern_inter)
                primary_induction_increment = (-temperature + self.param.temp_vern_max) * (self.param.daily_vern_rate * self.param.temp_vern_inter + self.param.basic_vern_rate) / (self.param.temp_vern_max - self.param.temp_vern_inter)
            else:
                primary_induction_increment = 0
            return primary_induction_increment

    def photoperiod_induction_function(self, daylength,i):
        """

        :param daylength:
        :return:
        """
        if isinstance(self.param.temp_vern_min, list):
            if daylength < self.param.photoperiod_min[i]:
                secondary_induction_increment = 0
            else:
                secondary_induction_increment = (daylength - self.param.photoperiod_min[i]) / (self.param.photoperiod_max[i] - self.param.photoperiod_min[i])
            return secondary_induction_increment
        else:
            if daylength < self.param.photoperiod_min:
                secondary_induction_increment = 0
            else:
                secondary_induction_increment = (daylength - self.param.photoperiod_min) / (self.param.photoperiod_max - self.param.photoperiod_min)
            return secondary_induction_increment

#    def final_phytomer_number



    #
    # def approx_final_leaf_number_calculation(self, daylength, potential_leaf_number):
    #     """
    #
    #     :param daylength:
    #     :param potential_leaf_number:
    #     :return:
    #     """
    #     approx_final_leaf_number = min(potential_leaf_number, potential_leaf_number + self.param.sldl * (self.param.saturation_daylength - daylength))
    #     return approx_final_leaf_number
    #
    #
    # def potential_leaf_number_calculation(self, primary_induction_index):
    #     """
    #
    #     :self.param primary_induction_index:
    #     :return:
    #     """
    #     potential_leaf_number = self.param.absolute_max_leaf_number - (self.param.absolute_max_leaf_number - self.param.absolute_min_leaf_number) * primary_induction_index
    #     return potential_leaf_number
    #
    # def test_vernalisation(self, primary_induction_index, potential_leaf_number, primordia_number):
    #     """
    #
    #     :self.param primary_induction_index:
    #     :self.param potential_leaf_number:
    #     :self.param primordia_number:
    #     :return:
    #     """
    #     if(primary_induction_index == 1 or potential_leaf_number == self.param.absolute_min_leaf_number or potential_leaf_number <= primordia_number):
    #         return True