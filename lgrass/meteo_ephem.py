# coding: utf8

import ephem
import pandas as pd
import numpy as np
from geopy.geocoders import Nominatim
from datetime import datetime, timedelta


def import_meteo_data(meteo_path, sowing_date, site):
    meteo_data = pd.read_csv(meteo_path, sep=',')
    meteo_data = meteo_data[meteo_data.site == site]
    meteo_data = thermal_time_calculation(meteo_data, sowing_date)
    meteo_data['experimental_day'] = list(range(1, len(meteo_data) + 1))
    if 'daylength' in meteo_data.columns:
        return meteo_data
    else:
        meteo_data = daylength_series(meteo_data)
        return meteo_data


def thermal_time_calculation(meteo_data, sowing_date):
    """
    :param meteo_data: dataframe with 2 columns:   date: format 'YYYY_mm_dd'
                                                   temperature: daily mean temperature (float in °C)
    :param sowing_date: date: format 'YYYY_mm_dd'
    :return:
    """
    meteo_data.date = pd.to_datetime(meteo_data.date, format='%Y_%m_%d')
    meteo_data = meteo_data.sort_values(by=['date'])
    meteo_data.loc[meteo_data.mean_temperature < 0, 'mean_temperature'] = 0
    if len(meteo_data[meteo_data.mean_temperature.isnull()].index) > 0:
        for id_missing_value in meteo_data[meteo_data.mean_temperature.isnull()].index:
            if id_missing_value == min(meteo_data.index) or id_missing_value == max(meteo_data.index):
                meteo_data.mean_temperature[meteo_data.index == id_missing_value] = 0
            else:
                prev_value = meteo_data.mean_temperature[meteo_data.index == id_missing_value - 1].item()
                next_value = meteo_data.mean_temperature[meteo_data.index == id_missing_value + 1].item()
                meteo_data.loc[id_missing_value, 'mean_temperature'] = np.mean([prev_value, next_value])
    meteo_data = meteo_data.loc[meteo_data['date'] >= pd.to_datetime(sowing_date, format='%Y_%m_%d')]
    meteo_data['thermal_time_cumul'] = meteo_data.mean_temperature.cumsum()
    return meteo_data


def set_observer(address):
    """
    This method creates a object of class 'ephem.Observer'
    Arguments:
    - address:
        address of the geographical location of the site to be simulated
        type: str
    """
    geolocator = Nominatim(user_agent="rouet")
    location = geolocator.geocode(address)
    latitude = location.latitude
    longitude = location.longitude
    elev = location.altitude
    obs = ephem.Observer()
    obs.lon = str(longitude)
    obs.lat = str(latitude)
    obs.elev = elev
    return obs


# obs.horizon = '-0:34'
# We relocate the horizon to get twilight times
# obs.horizon = '-6' #-6=civil twilight, -12=nautical, -18=astronomical
# beg_twilight=obs.previous_rising(ephem.Sun(), use_center=True) #Begin civil twilight
# end_twilight=fred.next_setting   (ephem.Sun(), use_center=True) #End civil twilight
# https://stackoverflow.com/questions/2637293/calculating-dawn-and-sunset-times-using-pyephem


def daylength_for_a_date(date, observer):
    """
    This method calculates daylength for one date and one object of class 'ephem.Observer'
    Arguments:
    - date:
        date of the day as a string : 'YYYY_mm_dd'
        type: str
    - observer:

    return a daylength in hours
    """
    date_at_noon = date + timedelta(hours=12)
    observer.date = date_at_noon
    sunrise = observer.previous_rising(ephem.Sun())                              # GTM hour
    sunrise = datetime.strptime(str(sunrise), '%Y/%m/%d %H:%M:%S')               # GTM hour
    sunset = observer.next_setting(ephem.Sun())                                  # GTM hour
    sunset = datetime.strptime(str(sunset), '%Y/%m/%d %H:%M:%S')                 # GTM hour
    dl = sunset - sunrise
    daylength = dl.seconds/3600.                    # Change from unit 'second' to unit 'hour'
    return daylength


def daylength_series(data):
    """
    This method creates a pandas.DataFrame with dates and associated daylengths for a location
    Arguments:
    - data:
        dataframe with a column named date 'YYYY_mm_dd'
        type: str
    - address:
        address of the geographical location of the site to be simulated
        type: str
    """
    observer = set_observer(data.iloc[0].site)
    data['daylength'] = data.apply(lambda x: daylength_for_a_date(x['date'], observer), axis=1).tolist()
    return data
