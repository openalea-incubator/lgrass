# This file contains all fonctions associated with flowering
import parameters


class FloweringFunctions:
    def __init__(self):
        self.param = parameters.ParametersValues()

    def vernalisation_function(self, temperature):
        """

        :param temperature:
        :return:
        """
        if self.param.temp_vern_min <= temperature <= self.param.temp_vern_inter:
            primary_induction_increment = self.param.vai * temperature + self.param.vbee
        elif self.param.temp_vern_inter < temperature <= self.param.temp_vern_max:
            primary_induction_increment = (self.param.vai * temperature + self.param.vbee) * (self.param.temp_vern_max - temperature) / (self.param.temp_vern_max - self.param.temp_vern_inter)
        else:
            primary_induction_increment = 0
        return primary_induction_increment


    def approx_final_leaf_number_calculation(self, daylength, potential_leaf_number):
        """

        :param daylength:
        :param potential_leaf_number:
        :return:
        """
        approx_final_leaf_number = min(potential_leaf_number, potential_leaf_number + self.param.sldl * (self.param.saturation_daylength - daylength))
        return approx_final_leaf_number


    def potential_leaf_number_calculation(self, primary_induction_index):
        """

        :self.param primary_induction_index:
        :return:
        """
        potential_leaf_number = self.param.absolute_max_leaf_number - (self.param.absolute_max_leaf_number - self.param.absolute_min_leaf_number) * primary_induction_index
        return potential_leaf_number

    def test_vernalisation(self, primary_induction_index, potential_leaf_number, primordia_number):
        """

        :self.param primary_induction_index:
        :self.param potential_leaf_number:
        :self.param primordia_number:
        :return:
        """
        if(primary_induction_index == 1 or potential_leaf_number == self.param.absolute_min_leaf_number or potential_leaf_number <= primordia_number):
            return True