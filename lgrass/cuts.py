# coding: utf8
# '''
# Created on 20/04/2020
#
# @author: modelisation - TR
# '''


# Fonction déterminant les jours de coupes programmés en fonction de la fréquence des tontes en jour
def define_cutting_dates(weather, max_degrees, cutting_freq):
    cutting_dates = []
    derivation_length = max_degrees
    degrees = 0
    do_cut = False
    cut = 0
    i = 0
    while degrees <= max_degrees:
        degrees += weather['mean_temperature'].iloc[i]
        i += 1
        if i == cut + cutting_freq:
            do_cut = True
        if do_cut & (degrees <= max_degrees):
            cut = i
            cutting_dates.append(cut)
            derivation_length += 3  # une tonte nécessite 3 itérations supplémentaires (x+0.25,x+0.5,x+0.75)
            do_cut = False
    return cutting_dates, derivation_length
