import pandas as pd
from lgrass import meteo_ephem

meteo_data = pd.read_csv('D:/Simon/Python/lgrass/example/GEVES/Meteo_GEVES_sans_photoperiode.csv', sep=',')
meteo_data.date = pd.to_datetime(meteo_data.date, format='%Y_%m_%d')

new_meteo_data = pd.DataFrame()

for site in meteo_data.site.unique():
    meteo_site = meteo_data.loc[meteo_data['site'] == site]
    a = meteo_ephem.daylength_series(meteo_site)
    new_meteo_data = new_meteo_data.append(a)

new_meteo_data.date = new_meteo_data.date.astype('string')
new_meteo_data['date'] = new_meteo_data['date'].str.replace('-', '_')

new_meteo_data.to_csv('D:/Simon/Python/lgrass/example/GEVES/Meteo_GEVES.csv', index=False)