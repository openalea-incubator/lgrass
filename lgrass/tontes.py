# coding: utf8


# Fonction déterminant les jours de coupes programmés en fonction de la fréquence des tontes en jour
def define_cutting_dates(weather, max_degrees, cutting_freq):
    cutting_dates = []
    degrees = 0
    do_cut = True
    cut = 0
    i = 0
    while degrees <= max_degrees:
        while (weather['daylength'].iloc[i] < 13) & (degrees <= max_degrees):
            degrees += weather['mean_temperature'].iloc[i]
            i += 1
        degrees += weather['mean_temperature'].iloc[i]
        i += 1
        if i == cut + cutting_freq:
            do_cut = True
        if do_cut & (degrees <= max_degrees):
            cut = i
            cutting_dates.append(cut)
            do_cut = False
    return cutting_dates


