# This file contains all fonctions associated with flowering
from lgrass import parameters


class FloweringFunctions:
    def __init__(self):
        self.param = parameters.ParametersValues()

    def vernalisation_function(self, temperature):
        """ Daily increment of primary induction. Triangluar function of daily temperature.

        :param float temperature: Daily mean air temperature (°C)

        :return: primary_induction_increment (dimensionless)
        :rtype: float
        """

        if self.param.temp_vern_min <= temperature <= self.param.temp_vern_inter:
            primary_induction_increment = (self.param.daily_vern_increment /
                                           (self.param.temp_vern_inter - self.param.temp_vern_min)) * \
                                          (temperature - self.param.temp_vern_min)
        elif self.param.temp_vern_inter < temperature <= self.param.temp_vern_max:
            primary_induction_increment = (self.param.daily_vern_increment /
                                           (self.param.temp_vern_inter - self.param.temp_vern_max)) * \
                                          (temperature - self.param.temp_vern_max)
        else:
            primary_induction_increment = 0
        return primary_induction_increment

    def PPR_function(self, leaf_number):
        """ Calculate daily secondary induction increment as a function of leaf number.

        :param int leaf_number: number of visible leaves on the tiller

        :return: max_photo_ind_rate (dimensionless)
        :rtype: float
        """

        if leaf_number <= self.param.leaf_number_max:
            max_photo_ind_increment = leaf_number * self.param.PPRM / self.param.leaf_number_max
        else:
            max_photo_ind_increment = self.param.PPRM
        return max_photo_ind_increment

    def photoperiod_induction_function(self, daylength, leaf_number):
        """ Daily increment of secondary induction if photoperiod longer than a threshold. Daily increment regulated by leaf number.

        :param float daylength: photoperiod of the day
        :param int leaf_number: number of visible leaves on the tiller

        :return: secondary_induction_increment (dimensionless)
        :rtype: float
        """

        if daylength < self.param.photoperiod_min:
            secondary_induction_increment = 0
        else:
            secondary_induction_increment = min(self.PPR_function(leaf_number), self.PPR_function(leaf_number) * (
                        daylength - self.param.photoperiod_min) / (
                                                            self.param.photoperiod_max - self.param.photoperiod_min))
        return secondary_induction_increment
