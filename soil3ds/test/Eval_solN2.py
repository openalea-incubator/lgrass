from scipy import *
from rpy_options import set_options
set_options(RHOME='c:/progra~1/R/R-2.12.1')
from rpy import r

#import sys
#path_ = r'H:\devel\grassland\grassland\luzerne'
#path2_ = r'H:\devel\grassland\grassland\L-gume'
#sys.path.insert(0, path_)
#sys.path.insert(0, path2_)

from soil_moduleN import * #! renommer car dans nouvelle version Lpy, mot module est reserve et fait planter!
from soil_modulevisu import *
from IOxls import *





####################
## Test sol Nu: usm Lusig99
## 3 compartiments de 30cm, profil sol lusignan99,  meteo 1999-2000, ZESX=60cm, CFES=5
## Nitrification et denitrification desactivees
####################

## lecture fichier meteo
#faire une fonction??
meteo_path = r'H:\Travail\devel\grassland\grassland\L-gume\meteo_exemple.xls'
ongletM = 'SolNu_Lusig99'#'testSolNu'#'exemple'#'competiluz'#
met = xlrd.open_workbook(meteo_path)
meteo = IOtable.conv_dataframe(get_xls_col(met.sheet_by_name(ongletM)))
for k in ['year', 'month', 'day', 'DOY','Coupe']: meteo[k] = list(map(int, meteo[k]))


## sol
pattern8 = [[-50,-50], [50,50]]#[[-2.5,-2.5], [5,5]]#[[-22.5,-22.5], [22.5,22.5]]#[[-35,-35], [35,35]]# #pattern.8 equivalent (cm)
Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
surfsolref = Lsol*largsol #m2
dz_sol = 30. #cm
ncouches_sol = 3

ZESX = 0.6 #m
CFES = 5.#coefficient de forme sans dimension [0.5 - 5.]
ACLIMc=26.#site de lusignan
ARGIs=18.3#lusignan99
CALCs = 0.2 #lusignan99
Norg = 1.1 #gN.kg-1 sol (*10 de Norg STICS exprime en %)
q0 = 9.464 #(mm de sol)
pH= 7.1#lusignan99
Tsol=15. #degresC
concrr = 0.02/10000. #kg N.m-2.mm-1 (concentration en N de la pluie)


### pour calculs d'humidite
#SN.HR()

CN0_30 = 9.5238 # (1/Wh de STICS = 0.105)
CN30_60 = 9.5238 # (1/Wh de STICS = 0.105)
CN60_90 = 6.51 #sans dim; moyenne competiluz

MO0_30 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)
MO30_60 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)
MO60_90 = mean([4.6, 4.8])#(g.kg-1) TransfertN

ARGIs0_30 = 18.3#lusignan99
ARGIs30_60 = 18.3#lusignan99
ARGIs60_90 = mean([30.4, 34.0])#(%) TransfertN


#vecteur d'initialisation
vDA = [1.46, 1.46, 1.41] #densite apparente de sol
vCN = [CN0_30]*1 + [CN30_60]*1 + [CN60_90]*1 #maxi 90cm en strates de 5cm
vMO = [MO0_30]*1 + [MO30_60]*1 + [MO60_90]*1 #maxi 90cm en strates de 5cm
vARGIs = [ARGIs0_30]*1 + [ARGIs30_60]*1 + [ARGIs60_90]*1 #maxi 90cm
vCALCs = [CALCs]*3
vNH4 = [7.7, 1.4, 0.9] #lusignan99 #!! kg d'N.ha-1 (entree de STICS)
vNO3 = [19.4, 3., 2.3] #lusignan99 # kg d'N.ha-1 (entree de STICS)






#recalcule a partir du fichier test sol nu qui reprend un sol de lusignan
par_sol = {'1':{'soil number': '1', 'soil type': 'H1 - Lusign99', 'teta_sat': 0.35, 'teta_fc': 0.3358, 'teta_wp': 0.1168, 'teta_ad': 0.011,   'WCST': '0.35',   'gamma_theo': '0.07'}}
par_sol['2'] = {'soil number': '2', 'soil type': 'H2 - Lusign99', 'teta_sat': 0.41, 'teta_fc': 0.3404, 'teta_wp': 0.1184, 'teta_ad': 0.026,   'WCST': '0.41',   'gamma_theo': '0.056'}
par_sol['3'] = {'soil number': '3', 'soil type': 'H3 - Lusign99', 'teta_sat': 0.35, 'teta_fc': 0.3243, 'teta_wp': 0.1734, 'teta_ad': 0.054,   'WCST': '0.35',   'gamma_theo': '0.038'}


#parametre sol complementaires pour l'N
par_SN = {}
par_SN['FMIN1G'] = 0.0006 #(day-1) (p145)
par_SN['FMIN2G'] = 0.0272 #(%ARGIS-1) para pot rate min a ARGIs (p145)
par_SN['FMIN3G'] = 0.0167 #(% CALC-1)para pot rate min a CALCss (p145)

par_SN['FINERTG'] = 0.65 #0.65 = default fraction of stable pool (p145) should be smaller in grassland & forest
par_SN['PROFHUMs'] = 35. # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220

par_SN['HMinMg'] = 0.3 #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
par_SN['HoptMg'] = 1. #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)

par_SN['TRefg'] = 15. #reference temperature
par_SN['FTEMHAg'] = 25. #asymptotic value of FTH (seen as a logistic response)
par_SN['FTEMHg'] = 0.120 #(K-1)
par_SN['FTEMHB'] = 145.

par_SN['FNXg'] = 0.5 #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
par_SN['PHMinNITg'] = 3. #pH min de nitrification (prop of Field Capacity) #value p149
par_SN['PHMaxNITg'] = 5.5 #pH max de nitrification (prop of Field Capacity) #value p149
par_SN['HMinNg'] = 0.67 #Humidite min de nitrification #value p149
par_SN['HoptNg'] = 1. #Humidite opt de nitrification  #value p149
par_SN['TNITMINg'] = 5. #Temperature min de nitrification #  (degreC)value p151
par_SN['TNITOPTg'] = 20. #Temperature opt de nitrification #  (degreC)value p151
par_SN['TNITMAXg'] = 45. #Temperature max de nitrification #  (degreC)value p151
par_SN['RATIONITs'] = 0. #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
par_SN['DIFNg'] = 0.018  #N diffusion coefficient at field capacity (cm2.day-1, p 161)



## soil initialisation
#S = SoilN(par_sol, soil_number = [1,2,3], dxyz = [[Lsol], [largsol], [dz_sol/100.]*ncouches_sol], vDA=DA,ZESX=ZESX , CFES=CFES, ARGIs = ARGIs, pH=pH, CALCs=CALCs)
S = SoilN(par_sol, par_SN, soil_number = [1,2,3], dxyz = [[Lsol], [largsol], [dz_sol/100.]*ncouches_sol], vDA=vDA, vCN=vCN,vMO=vMO, vARGIs = vARGIs,vNO3=vNO3,vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol,pH=pH, ZESX=ZESX , CFES=CFES)


#print S.m_soil_vox


#A inclure dans un bilan
S.Corg[0,0,0]*10000*0.35/0.3 #Corg initial en kg de C.ha-1 sur PROFHUM
S.Norg[0,0,0]*10000*0.35/0.3 #Corg initial en kg de N.ha-1 sur PROFHUM
S.InertCorg[0,0,0]*10000*0.35/0.3 #C inert initial en kg de C.ha-1 sur PROFHUM
S.InertNorg[0,0,0]*10000*0.35/0.3 #N inert initial en kg de C.ha-1 sur PROFHUM
S.Corg[0,0,0]/S.Norg[0,0,0]

S.m_NH4*10000
S.m_NO3*10000


##
#initialisation des humidite de sol -> mettre a jour la fonction! -> OK a partir des HRp (pas des HRv)
HRpinit = [11.8, 15.6, 18.3]
S.init_asw(HRp_init=HRpinit)

intialWC = sum3(S.tsw_t)

Uval = q0*0.1*sum(S.m_QH20fc[0])*surfsolref / (S.dxyz[2][0]*100.)#(epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
stateEV = [0.,0.,0.] #pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
b= bEV(ACLIMc, ARGIs, HXs=0.261)#1.#valeur empirique tres proche#0.1#0.63#0.63



#debut, fin de simulation
DOY_deb, DOY_fin = 239,623


##vegetation sans racine ni LAI
R1 = vert_roots(S.dxyz, [0.000000001,0.,0.,0.]) #pas zero sinon buf FTSW
ls_roots = [R1]
ls_epsi = [0.]


##bouble journaliere
cumEV,cumET0, cumPP, cumD, profH20, cumTransp = [],[],[], [], [], []
vlix, azomes = [], []
for DOY in range(DOY_deb, DOY_fin):
    meteo_j = extract_dataframe(meteo, ['TmoyDay','I0','Et0','Precip','Irrig','Coupe','FertNO3','FertNH4','Tsol'], 'DOY', val=DOY)
    for k in list(meteo_j.keys()): meteo_j[k]=meteo_j[k][0]

    #entrees eau
    #Precip = meteo_j['Precip']+meteo_j['Irrig']
    Rain = meteo_j['Precip']
    Irrig = meteo_j['Irrig']

    #entrees N
    #map_N = 0.*S.m_1[0,:,:]
    mapN_Rain = 1.*S.m_1[0,:,:] * Rain * concrr #Nmin de la pluie
    mapN_Irrig = 1.*S.m_1[0,:,:] * Irrig * concrr #Nmin de l'eau d'irrigation
    mapN_fertNO3 = 1.*S.m_1[0,:,:] * meteo_j['FertNO3'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    mapN_fertNH4 = 1.*S.m_1[0,:,:] * meteo_j['FertNH4'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel

    S.updateTsol(meteo_j['Tsol'])#(meteo_j['TmoyDay']) #Tsol forcee comme dans STICS


    ls_transp, evapo_tot, Drainage, stateEV,  m_frac_transpi, m_frac_evap, ls_ftsw = S.stepWBmc(meteo_j['Et0']*surfsolref, ls_roots, ls_epsi, Rain*surfsolref, Irrig*surfsolref, stateEV, ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    HAx = S.HRp()
    S.stepNB(par_SN)
    S.stepNitrif(par_SN)
    S.stepNINFILT(mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, Drainage, opt=1)


    #sorties
    print(DOY, S.tsw_t[0,0,0], evapo_tot) #sum3(S.tsw_t)
    cumEV.append(evapo_tot)
    cumTransp.append(sum(ls_transp))
    cumET0.append(meteo_j['Et0']*surfsolref)
    cumPP.append(meteo_j['Precip']*surfsolref)
    cumD.append(Drainage[-1][0][0])
    profH20.append([DOY]+HAx[:,0,0].tolist()+[evapo_tot, Drainage[-1][0][0]])

    vlix.append(S.lixiNO3*10000)
    azomes.append(sum3(S.m_NH4+S.m_NO3)*10000)


##termes du bilan hydrique global
intialWC
sum(cumPP)
intialWC+sum(cumPP)

sum3(S.tsw_t)
sum(cumEV)
sum(cumTransp)
sum(cumD)
sum3(S.tsw_t)+sum(cumEV)+sum(cumD)

sum(cumET0)

S.CloseWbalance() #-> equilibre
S.CloseCbalance() #-> equilibre
S.CloseNbalance() #-> equilibre




##termes du bilan N
S.Corg[0,0,0]*10000*0.35/0.3 #Corg initial en kg de C.ha-1 sur PROFHUM
S.Norg[0,0,0]*10000*0.35/0.3 #Corg initial en kg de N.ha-1 sur PROFHUM
S.InertCorg[0,0,0]*10000*0.35/0.3 #C inert initial en kg de C.ha-1 sur PROFHUM
S.InertNorg[0,0,0]*10000*0.35/0.3 #N inert initial en kg de C.ha-1 sur PROFHUM
S.Corg[0,0,0]/S.Norg[0,0,0]

S.m_NH4*10000
S.m_NO3*10000
S.N2ONitrif*10000
S.lixiNO3*10000

r.plot(vlix, ylab='', xlab='DOY')
r.plot(azomes, ylab='', xlab='DOY')

##ecriture sorties
#f = file(r'H:\simul\Valid_sol\Lusig99_profH2O_.csv', 'w')
#IOtable.ecriture_csv(profH20, f)
#f.close()



#A faire:

# lecture directe des parametres dans fichier soil.XML de STICS?
#minearalisation CN sous estimee: liee a temperature de sol -> hypothese de Tsol=Tair tes brute!
#introduire le calcul de Tcult, puis flux de chaleur dans le sol
# en forcant Tsol avec valeur de STICS, encore un ecart de mineralisation -> lie a difference dans dynamique d'humidite du sol (plus sec au moment ou mineralise le plus?)
#ordre de grandeur OK cependant





####################
## Test sol Nu: usm Lusig99
## 3 compartiments de 30cm, profil sol lusignan99,  meteo 1999-2000, ZESX=60cm, CFES=5
## discretise en couches de 10cm
## Nitrification et denitrification desactivees
####################

## lecture fichier meteo
#faire une fonction??
meteo_path = r'H:\Travail\devel\grassland\grassland\L-gume\meteo_exemple.xls'
ongletM = 'SolNu_Lusig99'#'testSolNu'#'exemple'#'competiluz'#
met = xlrd.open_workbook(meteo_path)
meteo = IOtable.conv_dataframe(get_xls_col(met.sheet_by_name(ongletM)))
for k in ['year', 'month', 'day', 'DOY','Coupe']: meteo[k] = list(map(int, meteo[k]))


## sol
pattern8 = [[-50,-50], [50,50]]#[[-2.5,-2.5], [5,5]]#[[-22.5,-22.5], [22.5,22.5]]#[[-35,-35], [35,35]]# #pattern.8 equivalent (cm)
Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
surfsolref = Lsol*largsol #m2
dz_sol = 10. #cm
ncouches_sol = 9

ZESX = 0.6 #m
CFES = 5.#coefficient de forme sans dimension [0.5 - 5.]
ACLIMc=26.#site de lusignan
ARGIs=18.3#lusignan99
CALCs = 0.2 #lusignan99
Norg = 1.1 #gN.kg-1 sol (*10 de Norg STICS exprime en %)
q0 = 9.464 #(mm de sol)
pH= 7.1#lusignan99
Tsol=15. #degresC
concrr = 0.02/10000. #kg N.m-2.mm-1 (concentration en N de la pluie)


### pour calculs d'humidite
#SN.HR()

CN0_30 = 9.5238 # (1/Wh de STICS = 0.105)
CN30_60 = 9.5238 # (1/Wh de STICS = 0.105)
CN60_90 = 6.51 #sans dim; moyenne competiluz

MO0_30 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)
MO30_60 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)
MO60_90 = mean([4.6, 4.8])#(g.kg-1) TransfertN

ARGIs0_30 = 18.3#lusignan99
ARGIs30_60 = mean([20.6, 21.1, 19.9])#(%) TransfertN
ARGIs60_90 = mean([30.4, 34.0])#(%) TransfertN


#vecteur d'initialisation
vDA = [1.46, 1.46, 1.46, 1.48, 1.48, 1.48, 1.41, 1.41, 1.41] #densite apparente de sol
vCN = [CN0_30]*3 + [CN30_60]*3 + [CN60_90]*3 #maxi 90cm en strates de 5cm
vMO = [MO0_30]*3 + [MO30_60]*3 + [MO60_90]*3 #maxi 90cm en strates de 5cm
vARGIs = [ARGIs0_30]*3 + [ARGIs30_60]*3 + [ARGIs60_90]*3 #maxi 90cm
vCALCs = [CALCs]*9
vNH4 = [7.7/3, 7.7/3, 7.7/3, 1.4/3, 1.4/3, 1.4/3, 0.9/3, 0.9/3, 0.9/3]#lusignan99 #!! kg d'N.ha-1 (entree de STICS)
vNO3 = [19.4/3, 19.4/3, 19.4/3, 3./3, 3./3, 3./3, 2.3/3, 2.3/3, 2.3/3] #lusignan99 # kg d'N.ha-1 (entree de STICS)



#recalcule a partir du fichier test sol nu qui reprend un sol de lusignan
par_sol = {'1':{'soil number': '1', 'soil type': 'H1 - Lusign99', 'teta_sat': 0.35, 'teta_fc': 0.3358, 'teta_wp': 0.1168, 'teta_ad': 0.011,   'WCST': '0.35',   'gamma_theo': '0.07'}}
par_sol['2'] = {'soil number': '2', 'soil type': 'H2 - Lusign99', 'teta_sat': 0.41, 'teta_fc': 0.3404, 'teta_wp': 0.1184, 'teta_ad': 0.026,   'WCST': '0.41',   'gamma_theo': '0.056'}
par_sol['3'] = {'soil number': '3', 'soil type': 'H3 - Lusign99', 'teta_sat': 0.35, 'teta_fc': 0.3243, 'teta_wp': 0.1734, 'teta_ad': 0.054,   'WCST': '0.35',   'gamma_theo': '0.038'}


#parametre sol complementaires pour l'N
par_SN = {}
par_SN['FMIN1G'] = 0.0006 #(day-1) (p145)
par_SN['FMIN2G'] = 0.0272 #(%ARGIS-1) para pot rate min a ARGIs (p145)
par_SN['FMIN3G'] = 0.0167 #(% CALC-1)para pot rate min a CALCss (p145)

par_SN['FINERTG'] = 0.65 #0.65 = default fraction of stable pool (p145) should be smaller in grassland & forest
par_SN['PROFHUMs'] = 35. # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220

par_SN['HMinMg'] = 0.3 #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
par_SN['HoptMg'] = 1. #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)

par_SN['TRefg'] = 15. #reference temperature
par_SN['FTEMHAg'] = 25. #asymptotic value of FTH (seen as a logistic response)
par_SN['FTEMHg'] = 0.120 #(K-1)
par_SN['FTEMHB'] = 145.

par_SN['FNXg'] = 0.5 #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
par_SN['PHMinNITg'] = 3. #pH min de nitrification (prop of Field Capacity) #value p149
par_SN['PHMaxNITg'] = 5.5 #pH max de nitrification (prop of Field Capacity) #value p149
par_SN['HMinNg'] = 0.67 #Humidite min de nitrification #value p149
par_SN['HoptNg'] = 1. #Humidite opt de nitrification  #value p149
par_SN['TNITMINg'] = 5. #Temperature min de nitrification #  (degreC)value p151
par_SN['TNITOPTg'] = 20. #Temperature opt de nitrification #  (degreC)value p151
par_SN['TNITMAXg'] = 45. #Temperature max de nitrification #  (degreC)value p151
par_SN['RATIONITs'] = 0. #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
par_SN['DIFNg'] = 0.018  #N diffusion coefficient at field capacity (cm2.day-1, p 161)





## soil initialisation
#S = SoilN(par_sol, soil_number = [1,2,3], dxyz = [[Lsol], [largsol], [dz_sol/100.]*ncouches_sol], vDA=DA,ZESX=ZESX , CFES=CFES, ARGIs = ARGIs, pH=pH, CALCs=CALCs)
S = SoilN(par_sol, par_SN, soil_number = [1,1,1,2,2,2,3,3,3], dxyz = [[Lsol/2,Lsol/2], [largsol/2, largsol/2], [dz_sol/100.]*ncouches_sol], vDA=vDA, vCN=vCN,vMO=vMO, vARGIs = vARGIs,vNO3=vNO3,vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol,pH=pH, ZESX=ZESX , CFES=CFES)

##
#initialisation des humidite de sol -> mettre a jour la fonction! -> OK a partir des HRp (pas des HRv)
HRpinit = [11.8, 11.8, 11.8, 15.6, 15.6, 15.6, 18.3, 18.3, 18.3]
S.init_asw(HRp_init=HRpinit)


Uval = q0*0.1*sum(S.m_QH20fc[0])*surfsolref / (S.dxyz[2][0]*100.)#(epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
stateEV = [0.,0.,0.] #pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
b= bEV(ACLIMc, ARGIs, HXs=0.261)#1.#valeur empirique tres proche#0.1#0.63#0.63



#debut, fin de simulation
DOY_deb, DOY_fin = 239,623


##vegetation sans racine ni LAI
R1 = vert_roots(S.dxyz, [0.000000001,0.,0.,0.,0.,0.,0.,0.,0.]) #pas zero sinon buf FTSW
ls_roots = [R1]
ls_epsi = [0.]


##bouble journaliere
cumEV,cumET0, cumPP, cumD, profH20, cumTransp = [],[],[], [], [], []
vlix, azomes = [], []
for DOY in range(DOY_deb, DOY_fin):
    meteo_j = extract_dataframe(meteo, ['TmoyDay','I0','Et0','Precip','Irrig','Coupe','FertNO3','FertNH4','Tsol'], 'DOY', val=DOY)
    for k in list(meteo_j.keys()): meteo_j[k]=meteo_j[k][0]

    #entrees eau
    #Precip = meteo_j['Precip']+meteo_j['Irrig']
    Rain = meteo_j['Precip']
    Irrig = meteo_j['Irrig']

    #entree N
    #map_N = 0.*S.m_1[0,:,:]
    mapN_Rain = 1.*S.m_1[0,:,:] * Rain * concrr #Nmin de la pluie
    mapN_Irrig = 1.*S.m_1[0,:,:] * Irrig * concrr #Nmin de l'eau d'irrigation
    mapN_fertNO3 = 1.*S.m_1[0,:,:] * meteo_j['FertNO3'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    mapN_fertNH4 = 1.*S.m_1[0,:,:] * meteo_j['FertNH4'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel

    S.updateTsol(meteo_j['Tsol'])#(meteo_j['TmoyDay']) #Tsol forcee comme dans STICS


    ls_transp, evapo_tot, Drainage, stateEV,  m_frac_transpi, m_frac_evap, ls_ftsw = S.stepWBmc(meteo_j['Et0']*surfsolref, ls_roots, ls_epsi, Rain*surfsolref, Irrig*surfsolref, stateEV, ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    HAx = S.HRp()
    S.stepNB(par_SN)
    S.stepNitrif(par_SN)
    S.stepNINFILT(mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, Drainage, opt=1)


    #sorties
    print(DOY, S.tsw_t[0,0,0], evapo_tot) #sum3(S.tsw_t)
    cumEV.append(evapo_tot)
    cumTransp.append(sum(ls_transp))
    cumET0.append(meteo_j['Et0']*surfsolref)
    cumPP.append(meteo_j['Precip']*surfsolref)
    cumD.append(Drainage[-1][0][0])
    profH20.append([DOY]+HAx[:,0,0].tolist()+[evapo_tot, Drainage[-1][0][0]])

    vlix.append(S.lixiNO3*10000)
    azomes.append(sum3(S.m_NH4+S.m_NO3)*10000)


S.CloseWbalance() #-> equilibre
S.CloseCbalance() #-> equilibre
S.CloseNbalance() #-> equilibre






####################
## Test sol Nu (3): usm Lusig99Res
## 3 compartiments de 30cm, profil sol lusignan99,  meteo 1999-2000, ZESX=60cm, CFES=5
## Nitrification et denitrification desactivees
## ajout de 20T d'un residu a DOY 300
####################

## lecture fichier meteo
#faire une fonction??
meteo_path = r'H:\Travail\devel\grassland\grassland\L-gume\meteo_exemple.xls'
ongletM = 'SolNu_Lusig99'#'testSolNu'#'exemple'#'competiluz'#
met = xlrd.open_workbook(meteo_path)
meteo = IOtable.conv_dataframe(get_xls_col(met.sheet_by_name(ongletM)))
for k in ['year', 'month', 'day', 'DOY','Coupe']: meteo[k] = list(map(int, meteo[k]))


## sol
pattern8 = [[-50,-50], [50,50]]#[[-2.5,-2.5], [5,5]]#[[-22.5,-22.5], [22.5,22.5]]#[[-35,-35], [35,35]]# #pattern.8 equivalent (cm)
Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
surfsolref = Lsol*largsol #m2
dz_sol = 30. #cm
ncouches_sol = 3

ZESX = 0.6 #m
CFES = 5.#coefficient de forme sans dimension [0.5 - 5.]
ACLIMc=26.#site de lusignan
ARGIs=18.3#lusignan99
CALCs = 0.2 #lusignan99
Norg = 1.1 #gN.kg-1 sol (*10 de Norg STICS exprime en %)
q0 = 9.464 #(mm de sol)
pH= 7.1#lusignan99
Tsol=15. #degresC
concrr = 0.02/10000. #kg N.m-2.mm-1 (concentration en N de la pluie)


### pour calculs d'humidite
#SN.HR()

CN0_30 = 9.5238 # (1/Wh de STICS = 0.105)
CN30_60 = 9.5238 # (1/Wh de STICS = 0.105)
CN60_90 = 6.51 #sans dim; moyenne competiluz

MO0_30 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)
MO30_60 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)
MO60_90 = mean([4.6, 4.8])#(g.kg-1) TransfertN

ARGIs0_30 = 18.3#lusignan99
ARGIs30_60 = 18.3#lusignan99
ARGIs60_90 = mean([30.4, 34.0])#(%) TransfertN


#vecteur d'initialisation
vDA = [1.46, 1.46, 1.41] #densite apparente de sol
vCN = [CN0_30]*1 + [CN30_60]*1 + [CN60_90]*1 #maxi 90cm en strates de 5cm
vMO = [MO0_30]*1 + [MO30_60]*1 + [MO60_90]*1 #maxi 90cm en strates de 5cm
vARGIs = [ARGIs0_30]*1 + [ARGIs30_60]*1 + [ARGIs60_90]*1 #maxi 90cm
vCALCs = [CALCs]*3
vNH4 = [7.7, 1.4, 0.9] #lusignan99 #!! kg d'N.ha-1 (entree de STICS)
vNO3 = [19.4, 3., 2.3] #lusignan99 # kg d'N.ha-1 (entree de STICS)

#1 residu = listes de 1 element
vAmount= [20.]# T Fresh Weight.ha-1 (equivalent QRES)
vCNRESt = [16.] #g.g-1 (equivalent CsurNres)
Vprop1 = [1.]+2*[0.] #distribution dans les horizons
vProps= [Vprop1]#[Vprop1]#[Vprop1, Vprop1, Vprop1]
vWC=[0.7]# fraction d'eau des residu frais (equivalent de Crespc /(%)/100)
vCC=[0.42]# fraction de C des residus sec (equivalent de Crespc /(%)/100)
vNmires = [0.00197]# fraction de poids frais residu en azote mineral (equivalent de Nminres(%)/100)



#recalcule a partir du fichier test sol nu qui reprend un sol de lusignan
par_sol = {'1':{'soil number': '1', 'soil type': 'H1 - Lusign99', 'teta_sat': 0.35, 'teta_fc': 0.3358, 'teta_wp': 0.1168, 'teta_ad': 0.011,   'WCST': '0.35',   'gamma_theo': '0.07'}}
par_sol['2'] = {'soil number': '2', 'soil type': 'H2 - Lusign99', 'teta_sat': 0.41, 'teta_fc': 0.3404, 'teta_wp': 0.1184, 'teta_ad': 0.026,   'WCST': '0.41',   'gamma_theo': '0.056'}
par_sol['3'] = {'soil number': '3', 'soil type': 'H3 - Lusign99', 'teta_sat': 0.35, 'teta_fc': 0.3243, 'teta_wp': 0.1734, 'teta_ad': 0.054,   'WCST': '0.35',   'gamma_theo': '0.038'}


#parametre sol complementaires pour l'N
par_SN = {}
par_SN['FMIN1G'] = 0.0006 #(day-1) (p145)
par_SN['FMIN2G'] = 0.0272 #(%ARGIS-1) para pot rate min a ARGIs (p145)
par_SN['FMIN3G'] = 0.0167 #(% CALC-1)para pot rate min a CALCss (p145)

par_SN['FINERTG'] = 0.65 #0.65 = default fraction of stable pool (p145) should be smaller in grassland & forest
par_SN['PROFHUMs'] = 35. # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220

par_SN['HMinMg'] = 0.3 #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
par_SN['HoptMg'] = 1. #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)

par_SN['TRefg'] = 15. #reference temperature
par_SN['FTEMHAg'] = 25. #asymptotic value of FTH (seen as a logistic response)
par_SN['FTEMHg'] = 0.120 #(K-1)
par_SN['FTEMHB'] = 145.

par_SN['FNXg'] = 0.5 #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
par_SN['PHMinNITg'] = 3. #pH min de nitrification (prop of Field Capacity) #value p149
par_SN['PHMaxNITg'] = 5.5 #pH max de nitrification (prop of Field Capacity) #value p149
par_SN['HMinNg'] = 0.67 #Humidite min de nitrification #value p149
par_SN['HoptNg'] = 1. #Humidite opt de nitrification  #value p149
par_SN['TNITMINg'] = 5. #Temperature min de nitrification #  (degreC)value p151
par_SN['TNITOPTg'] = 20. #Temperature opt de nitrification #  (degreC)value p151
par_SN['TNITMAXg'] = 45. #Temperature max de nitrification #  (dereC)value p151
par_SN['RATIONITs'] = 0. #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
par_SN['DIFNg'] = 0.018  #N diffusion coefficient at field capacity (cm2.day-1, p 161)





## soil initialisation
#S = SoilN(par_sol, soil_number = [1,2,3], dxyz = [[Lsol], [largsol], [dz_sol/100.]*ncouches_sol], vDA=DA,ZESX=ZESX , CFES=CFES, ARGIs = ARGIs, pH=pH, CALCs=CALCs)
S = SoilN(par_sol, par_SN, soil_number = [1,2,3], dxyz = [[Lsol], [largsol], [dz_sol/100.]*ncouches_sol], vDA=vDA, vCN=vCN,vMO=vMO, vARGIs = vARGIs,vNO3=vNO3,vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol,pH=pH, ZESX=ZESX , CFES=CFES)


##
#initialisation des humidite de sol -> mettre a jour la fonction! -> OK a partir des HRp (pas des HRv)
HRpinit = [11.8, 15.6, 18.3]
S.init_asw(HRp_init=HRpinit)

Uval = q0*0.1*sum(S.m_QH20fc[0])*surfsolref / (S.dxyz[2][0]*100.)#(epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
stateEV = [0.,0.,0.] #pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
b= bEV(ACLIMc, ARGIs, HXs=0.261)#1.#valeur empirique tres proche#0.1#0.63#0.63


#debut, fin de simulation
DOY_deb, DOY_fin = 239,623
DOYres = 300 #jour d'ajout des residus

##vegetation sans racine ni LAI
R1 = vert_roots(S.dxyz, [0.000000001,0.,0.,0.]) #pas zero sinon buf FTSW
ls_roots = [R1]
ls_epsi = [0.]


##bouble journaliere
cumEV,cumET0, cumPP, cumD, profH20, cumTransp = [],[],[], [], [], []
vlix, azomes = [], []
cresi, cbio = [],[]
minNH4, minNO3 = [], []#pour verif que pas de valeurs negatives
for DOY in range(DOY_deb, DOY_fin):
    meteo_j = extract_dataframe(meteo, ['TmoyDay','I0','Et0','Precip','Irrig','Coupe','FertNO3','FertNH4','Tsol'], 'DOY', val=DOY)
    for k in list(meteo_j.keys()): meteo_j[k]=meteo_j[k][0]

    if DOY == DOYres:
        S.init_residues(vCNRESt, vAmount, vProps, vWC, vCC)

    #entree eau
    #Precip = meteo_j['Precip']+meteo_j['Irrig']
    Rain = meteo_j['Precip']
    Irrig = meteo_j['Irrig']

    #eentrees N
    #map_N = 0.*S.m_1[0,:,:]
    mapN_Rain = 1.*S.m_1[0,:,:] * Rain * concrr #Nmin de la pluie
    mapN_Irrig = 1.*S.m_1[0,:,:] * Irrig * concrr #Nmin de l'eau d'irrigation
    mapN_fertNO3 = 1.*S.m_1[0,:,:] * meteo_j['FertNO3'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    mapN_fertNH4 = 1.*S.m_1[0,:,:] * meteo_j['FertNH4'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel

    S.updateTsol(meteo_j['Tsol'])#(meteo_j['TmoyDay']) #Tsol forcee comme dans STICS


    ls_transp, evapo_tot, Drainage, stateEV,  m_frac_transpi, m_frac_evap, ls_ftsw = S.stepWBmc(meteo_j['Et0']*surfsolref, ls_roots, ls_epsi, Rain*surfsolref, Irrig*surfsolref, stateEV, ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    HAx = S.HRp()
    S.stepNB(par_SN)
    if DOY>=DOYres:
        S.stepResidueMin(par_SN)
        minNH4.append(S.m_NH4.min())
        minNO3.append(S.m_NO3.min())
        S.stepMicrobioMin(par_SN)
        cresi.append(sum(S.ls_CRES)*10000)
        cbio.append(sum(S.ls_CBio)*10000)
    else:
        cresi.append(0.)
        cbio.append(0.)

    S.stepNitrif(par_SN)
    S.stepNINFILT(mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, Drainage, opt=1)

    #sorties
    print(DOY, S.tsw_t[0,0,0], evapo_tot) #sum3(S.tsw_t)
    cumEV.append(evapo_tot)
    cumTransp.append(sum(ls_transp))
    cumET0.append(meteo_j['Et0']*surfsolref)
    cumPP.append(meteo_j['Precip']*surfsolref)
    cumD.append(Drainage[-1][0][0])
    profH20.append([DOY]+HAx[:,0,0].tolist()+[evapo_tot, Drainage[-1][0][0]])

    vlix.append(S.lixiNO3*10000)
    azomes.append(sum3(S.m_NH4+S.m_NO3)*10000)


##termes du bilan hydrique global

S.CloseWbalance() #-> equilibre
S.CloseCbalance() #-> equilibre
S.CloseNbalance() #-> equilibre






####################
## Test sol Nu: usm Lusig99
## 3 compartiments de 30cm, profil sol lusignan99,  meteo 1999-2000, ZESX=60cm, CFES=5
## discretise en couches de 10cm
## Nitrification et denitrification desactivees
## ajout de 20T d'un residu a DOY 300
####################

## lecture fichier meteo
#faire une fonction??
meteo_path = r'H:\Travail\devel\grassland\grassland\L-gume\meteo_exemple.xls'
ongletM = 'SolNu_Lusig99'#'testSolNu'#'exemple'#'competiluz'#
met = xlrd.open_workbook(meteo_path)
meteo = IOtable.conv_dataframe(get_xls_col(met.sheet_by_name(ongletM)))
for k in ['year', 'month', 'day', 'DOY','Coupe']: meteo[k] = list(map(int, meteo[k]))


## sol
pattern8 = [[-50,-50], [50,50]]#[[-2.5,-2.5], [5,5]]#[[-22.5,-22.5], [22.5,22.5]]#[[-35,-35], [35,35]]# #pattern.8 equivalent (cm)
Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
surfsolref = Lsol*largsol #m2
dz_sol = 10. #cm
ncouches_sol = 9

ZESX = 0.6 #m
CFES = 5.#coefficient de forme sans dimension [0.5 - 5.]
ACLIMc=26.#site de lusignan
ARGIs=18.3#lusignan99
CALCs = 0.2 #lusignan99
Norg = 1.1 #gN.kg-1 sol (*10 de Norg STICS exprime en %)
q0 = 9.464 #(mm de sol)
pH= 7.1#lusignan99
Tsol=15. #degresC
concrr = 0.02/10000. #kg N.m-2.mm-1 (concentration en N de la pluie)


### pour calculs d'humidite
#SN.HR()

CN0_30 = 9.5238 # (1/Wh de STICS = 0.105)
CN30_60 = 9.5238 # (1/Wh de STICS = 0.105)
CN60_90 = 6.51 #sans dim; moyenne competiluz

MO0_30 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)
MO30_60 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)
MO60_90 = mean([4.6, 4.8])#(g.kg-1) TransfertN

ARGIs0_30 = 18.3#lusignan99
ARGIs30_60 = mean([20.6, 21.1, 19.9])#(%) TransfertN
ARGIs60_90 = mean([30.4, 34.0])#(%) TransfertN


#vecteur d'initialisation
vDA = [1.46, 1.46, 1.46, 1.48, 1.48, 1.48, 1.41, 1.41, 1.41] #densite apparente de sol
vCN = [CN0_30]*3 + [CN30_60]*3 + [CN60_90]*3 #maxi 90cm en strates de 5cm
vMO = [MO0_30]*3 + [MO30_60]*3 + [MO60_90]*3 #maxi 90cm en strates de 5cm
vARGIs = [ARGIs0_30]*3 + [ARGIs30_60]*3 + [ARGIs60_90]*3 #maxi 90cm
vCALCs = [CALCs]*9
vNH4 = [7.7/3, 7.7/3, 7.7/3, 1.4/3, 1.4/3, 1.4/3, 0.9/3, 0.9/3, 0.9/3]#lusignan99 #!! kg d'N.ha-1 (entree de STICS)
vNO3 = [19.4/3, 19.4/3, 19.4/3, 3./3, 3./3, 3./3, 2.3/3, 2.3/3, 2.3/3] #lusignan99 # kg d'N.ha-1 (entree de STICS)


#1 residu = listes de 1 element
vAmount= [20.]# T Fresh Weight.ha-1 (equivalent QRES)
vCNRESt = [16.] #g.g-1 (equivalent CsurNres)
Vprop1 = [1./3., 1./3., 1./3.]+6*[0.] #distribution dans les horizons
vProps= [Vprop1]#[Vprop1]#[Vprop1, Vprop1, Vprop1]
vWC=[0.7]# fraction d'eau des residu frais (equivalent de Crespc /(%)/100)
vCC=[0.42]# fraction de C des residus sec (equivalent de Crespc /(%)/100)
vNmires = [0.00197]# fraction de poids frais residu en azote mineral (equivalent de Nminres(%)/100)





#recalcule a partir du fichier test sol nu qui reprend un sol de lusignan
par_sol = {'1':{'soil number': '1', 'soil type': 'H1 - Lusign99', 'teta_sat': 0.35, 'teta_fc': 0.3358, 'teta_wp': 0.1168, 'teta_ad': 0.011,   'WCST': '0.35',   'gamma_theo': '0.07'}}
par_sol['2'] = {'soil number': '2', 'soil type': 'H2 - Lusign99', 'teta_sat': 0.41, 'teta_fc': 0.3404, 'teta_wp': 0.1184, 'teta_ad': 0.026,   'WCST': '0.41',   'gamma_theo': '0.056'}
par_sol['3'] = {'soil number': '3', 'soil type': 'H3 - Lusign99', 'teta_sat': 0.35, 'teta_fc': 0.3243, 'teta_wp': 0.1734, 'teta_ad': 0.054,   'WCST': '0.35',   'gamma_theo': '0.038'}


#parametre sol complementaires pour l'N
par_SN = {}
par_SN['FMIN1G'] = 0.0006 #(day-1) (p145)
par_SN['FMIN2G'] = 0.0272 #(%ARGIS-1) para pot rate min a ARGIs (p145)
par_SN['FMIN3G'] = 0.0167 #(% CALC-1)para pot rate min a CALCss (p145)

par_SN['FINERTG'] = 0.65 #0.65 = default fraction of stable pool (p145) should be smaller in grassland & forest
par_SN['PROFHUMs'] = 35. # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220

par_SN['HMinMg'] = 0.3 #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
par_SN['HoptMg'] = 1. #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)

par_SN['TRefg'] = 15. #reference temperature
par_SN['FTEMHAg'] = 25. #asymptotic value of FTH (seen as a logistic response)
par_SN['FTEMHg'] = 0.120 #(K-1)
par_SN['FTEMHB'] = 145.

par_SN['FNXg'] = 0.5 #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
par_SN['PHMinNITg'] = 3. #pH min de nitrification (prop of Field Capacity) #value p149
par_SN['PHMaxNITg'] = 5.5 #pH max de nitrification (prop of Field Capacity) #value p149
par_SN['HMinNg'] = 0.67 #Humidite min de nitrification #value p149
par_SN['HoptNg'] = 1. #Humidite opt de nitrification  #value p149
par_SN['TNITMINg'] = 5. #Temperature min de nitrification #  (degreC)value p151
par_SN['TNITOPTg'] = 20. #Temperature opt de nitrification #  (degreC)value p151
par_SN['TNITMAXg'] = 45. #Temperature max de nitrification #  (degreC)value p151
par_SN['RATIONITs'] = 0. #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
par_SN['DIFNg'] = 0.018  #N diffusion coefficient at field capacity (cm2.day-1, p 161)





## soil initialisation
#S = SoilN(par_sol, soil_number = [1,2,3], dxyz = [[Lsol], [largsol], [dz_sol/100.]*ncouches_sol], vDA=DA,ZESX=ZESX , CFES=CFES, ARGIs = ARGIs, pH=pH, CALCs=CALCs)
S = SoilN(par_sol, par_SN, soil_number = [1,1,1,2,2,2,3,3,3], dxyz = [[Lsol/2,Lsol/2], [largsol/2, largsol/2], [dz_sol/100.]*ncouches_sol], vDA=vDA, vCN=vCN,vMO=vMO, vARGIs = vARGIs,vNO3=vNO3,vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol,pH=pH, ZESX=ZESX , CFES=CFES)

##
#initialisation des humidite de sol -> mettre a jour la fonction! -> OK a partir des HRp (pas des HRv)
HRpinit = [11.8, 11.8, 11.8, 15.6, 15.6, 15.6, 18.3, 18.3, 18.3]
S.init_asw(HRp_init=HRpinit)


Uval = q0*0.1*sum(S.m_QH20fc[0])*surfsolref / (S.dxyz[2][0]*100.)#(epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
stateEV = [0.,0.,0.] #pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
b= bEV(ACLIMc, ARGIs, HXs=0.261)#1.#valeur empirique tres proche#0.1#0.63#0.63



#debut, fin de simulation
DOY_deb, DOY_fin = 239,623
DOYres = 300 #jour d'ajout des residus


##vegetation sans racine ni LAI
R1 = vert_roots(S.dxyz, [0.000000001,0.,0.,0.,0.,0.,0.,0.,0.]) #pas zero sinon buf FTSW
ls_roots = [R1]
ls_epsi = [0.]


##bouble journaliere
cumEV,cumET0, cumPP, cumD, profH20, cumTransp = [],[],[], [], [], []
vlix, azomes = [], []
for DOY in range(DOY_deb, DOY_fin):
    meteo_j = extract_dataframe(meteo, ['TmoyDay','I0','Et0','Precip','Irrig','Coupe','FertNO3','FertNH4','Tsol'], 'DOY', val=DOY)
    for k in list(meteo_j.keys()): meteo_j[k]=meteo_j[k][0]

    if DOY == DOYres:
        S.init_residues(vCNRESt, vAmount, vProps, vWC, vCC)

    #entrees eau
    Rain = meteo_j['Precip']
    Irrig = meteo_j['Irrig']

    #entrees N
    #map_N = 0.*S.m_1[0,:,:]
    mapN_Rain = 1.*S.m_1[0,:,:] * Rain * concrr #Nmin de la pluie
    mapN_Irrig = 1.*S.m_1[0,:,:] * Irrig * concrr #Nmin de l'eau d'irrigation
    mapN_fertNO3 = 1.*S.m_1[0,:,:] * meteo_j['FertNO3'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    mapN_fertNH4 = 1.*S.m_1[0,:,:] * meteo_j['FertNH4'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel

    S.updateTsol(meteo_j['Tsol'])#(meteo_j['TmoyDay']) #Tsol forcee comme dans STICS
    ls_transp, evapo_tot, Drainage, stateEV,  m_frac_transpi, m_frac_evap, ls_ftsw = S.stepWBmc(meteo_j['Et0']*surfsolref, ls_roots, ls_epsi, Rain*surfsolref, Irrig*surfsolref, stateEV, ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    HAx = S.HRp()
    S.stepNB(par_SN)
    if DOY>=DOYres:
        S.stepResidueMin(par_SN)
        S.stepMicrobioMin(par_SN)

    S.stepNitrif(par_SN)
    S.stepNINFILT(mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, Drainage, opt=1)

    #sorties
    print(DOY)


S.CloseWbalance() #-> equilibre
S.CloseCbalance() #-> equilibre
S.CloseNbalance() #-> equilibre

#tourne, et OK: tous les bilans equilibres independamment de la discretisation!
# petit effet de la discretisation sur la mineralisation des residus





#####################
## Test sol Nu: usm solnu_Llzir82 et solnu_Llzsc82(These JLD)
## 4 compartiments, profil sol lusignan,  meteo 1982, ZESX=39cm , CFES=3.5
## Tsol forcees comme stics
##
####################

## lecture fichier meteo
#faire une fonction??
meteo_path = r'H:\Travail\devel\grassland\grassland\L-gume\meteo_exemple.xls'
ongletM = 'Llzsc82'#'Llzir82'#
met = xlrd.open_workbook(meteo_path)
meteo = IOtable.conv_dataframe(get_xls_col(met.sheet_by_name(ongletM)))
for k in ['year', 'month', 'day', 'DOY','Coupe']: meteo[k] = list(map(int, meteo[k]))




## sol
pattern8 = [[-50,-50], [50,50]]#pattern.8 equivalent (cm)
Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
surfsolref = Lsol*largsol #m2
dz_sol = [0.25, 0.30, 0.35, 0.30] #m
ncouches_sol = 4

ZESX = 0.39 #m
CFES = 3.5#coefficient de forme sans dimension [0.5 - 5.]
ACLIMc=26.#site de lusignan
ARGIs=22.#lusignan
CALCs = 1. #lusignan
Norg = 1.2 #gN.kg-1 sol (*10 de Norg STICS exprime en %)
q0 = 9. #(mm de sol)
pH= 7.#lusignan
Tsol=15. #degresC
concrr = 0.02/10000. #kg N.m-2.mm-1 (concentration en N de la pluie)

### pour calculs d'humidite

CN0_30 = 9.5238 # (1/Wh de STICS = 0.105)
MO0_30 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)


#vecteur d'initialisation
vDA = [1.16, 1.29, 1.15, 0.91] #densite apparente de sol
vCN = [CN0_30]*ncouches_sol #maxi 90cm en strates de 5cm
vMO = [MO0_30]*ncouches_sol #maxi 90cm en strates de 5cm
vARGIs = [ARGIs]*ncouches_sol #maxi 90cm
vCALCs = [CALCs]*ncouches_sol
vNH4 = [0., 0., 0., 0.] # #!! kg d'N.ha-1 (entree de STICS)
vNO3 = [32., 12., 9., 0.] # kg d'N.ha-1 (entree de STICS)
HRpinit = [22.5, 24.4, 25.8, 28.7]



#recalcule a partir du fichier sol de lusignan
par_sol = {'1':{'soil number': '1', 'soil type': 'H1 - Lusignan', 'teta_sat': 0.35, 'teta_fc': 0.261, 'teta_wp': 0.0986, 'teta_ad': 0.011,   'WCST': '0.35',   'gamma_theo': '0.07'}}
par_sol['2'] = {'soil number': '2', 'soil type': 'H2 - Lusignan', 'teta_sat': 0.41, 'teta_fc': 0.3277, 'teta_wp': 0.1496, 'teta_ad': 0.026,   'WCST': '0.41',   'gamma_theo': '0.056'}
par_sol['3'] = {'soil number': '3', 'soil type': 'H3 - Lusignan', 'teta_sat': 0.35, 'teta_fc': 0.2967, 'teta_wp': 0.1759, 'teta_ad': 0.054,   'WCST': '0.35',   'gamma_theo': '0.038'}
par_sol['4'] = {'soil number': '3', 'soil type': 'H4 - Lusignan', 'teta_sat': 0.35, 'teta_fc': 0.2611, 'teta_wp': 0.1729, 'teta_ad': 0.054,   'WCST': '0.35',   'gamma_theo': '0.038'}
#Huv = Hup*DA


#parametre sol complementaires pour l'N
par_SN = {}
par_SN['FMIN1G'] = 0.0006 #(day-1) (p145)
par_SN['FMIN2G'] = 0.0272 #(%ARGIS-1) para pot rate min a ARGIs (p145)
par_SN['FMIN3G'] = 0.0167 #(% CALC-1)para pot rate min a CALCss (p145)

par_SN['FINERTG'] = 0.65 #0.65 = default fraction of stable pool (p145) should be smaller in grassland & forest
par_SN['PROFHUMs'] = 35. # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220

par_SN['HMinMg'] = 0.3 #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
par_SN['HoptMg'] = 1. #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)

par_SN['TRefg'] = 15. #reference temperature
par_SN['FTEMHAg'] = 25. #asymptotic value of FTH (seen as a logistic response)
par_SN['FTEMHg'] = 0.120 #(K-1)
par_SN['FTEMHB'] = 145.

par_SN['FNXg'] = 0.5 #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
par_SN['PHMinNITg'] = 3. #pH min de nitrification (prop of Field Capacity) #value p149
par_SN['PHMaxNITg'] = 5.5 #pH max de nitrification (prop of Field Capacity) #value p149
par_SN['HMinNg'] = 0.67 #Humidite min de nitrification #value p149
par_SN['HoptNg'] = 1. #Humidite opt de nitrification  #value p149
par_SN['TNITMINg'] = 5. #Temperature min de nitrification #  (degreC)value p151
par_SN['TNITOPTg'] = 20. #Temperature opt de nitrification #  (degreC)value p151
par_SN['TNITMAXg'] = 45. #Temperature max de nitrification #  (degreC)value p151
par_SN['RATIONITs'] = 0. #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
par_SN['DIFNg'] = 0.018  #N diffusion coefficient at field capacity (cm2.day-1, p 161)


#correctif de profhum pour tenir compte des changements de densite entre horizons
par_SN['PROFHUMs'] = 25. + 10.* vDA[0]/vDA[1]


## soil initialisation
S = SoilN(par_sol, par_SN, soil_number = [1,2,3,4], dxyz = [[Lsol], [largsol], dz_sol], vDA=vDA, vCN=vCN,vMO=vMO, vARGIs = vARGIs,vNO3=vNO3,vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol,pH=pH, ZESX=ZESX , CFES=CFES)
S.init_asw(HRp_init=HRpinit)




Uval = q0*0.1*sum(S.m_QH20fc[0])*surfsolref / (S.dxyz[2][0]*100.)#(epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
stateEV = [0.,0.,0.] #pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
b= bEV(ACLIMc, ARGIs, HXs=0.261)#1.#valeur empirique tres proche#0.1#0.63#0.63


#debut, fin de simulation
DOY_deb, DOY_fin = 1,360


##vegetation sans racine ni LAI
R1 = vert_roots(S.dxyz, [0.000000001,0.,0.,0.]) #pas zero sinon buf FTSW
ls_roots = [R1]
ls_epsi = [0.]



##boucle journaliere
vlix, azomes, resmes = [], [], []
for DOY in range(DOY_deb, DOY_fin):

    meteo_j = extract_dataframe(meteo, ['TmoyDay','I0','Et0','Precip','Irrig','Coupe','FertNO3','FertNH4','Tsolnu'], 'DOY', val=DOY)
    for k in list(meteo_j.keys()): meteo_j[k]=meteo_j[k][0]

    #Entrees eau
    Rain = meteo_j['Precip']
    Irrig = meteo_j['Irrig']

    #entrees N
    #map_N = 1.*S.m_1[0,:,:] * Precip * concrr
    mapN_Rain = 1.*S.m_1[0,:,:] * Rain * concrr #Nmin de la pluie
    mapN_Irrig = 1.*S.m_1[0,:,:] * Irrig * concrr #Nmin de l'eau d'irrigation
    mapN_fertNO3 = 1.*S.m_1[0,:,:] * meteo_j['FertNO3'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    mapN_fertNH4 = 1.*S.m_1[0,:,:] * meteo_j['FertNH4'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    S.updateTsol(meteo_j['Tsolnu'])#(meteo_j['TmoyDay']) #Tsol forcee comme dans STICS


    ls_transp, evapo_tot, Drainage, stateEV,  m_transpi, m_evap, ls_ftsw = S.stepWBmc(meteo_j['Et0']*surfsolref, ls_roots, ls_epsi, Rain*surfsolref, Irrig*surfsolref, stateEV, ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    S.stepNB(par_SN)
    S.stepNitrif(par_SN)
    ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N, idmin = S.stepNuptakePlt(par_SN)
    #S.stepNINFILT(map_N, Drainage, opt=1)
    S.stepNINFILT(mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, Drainage, opt=1)

    #sorties
    print(DOY)
    vlix.append(S.lixiNO3*10000)
    azomes.append(sum3(S.m_NH4+S.m_NO3)*10000)
    resmes.append(sum(S.tsw_t))


#r.plot(range(DOY_deb, DOY_fin), resmes, ylab='resmes', xlab='DOY', type='l')
#r.plot(range(DOY_deb, DOY_fin), azomes, ylab='azomes', xlab='DOY', type='l')
#r.plot(range(DOY_deb, DOY_fin), vlix, ylab='lix', xlab='DOY', type='l')








#####################
## Test Uptake plante: usm test_Llzir82 et test_Llzsc82(These JLD) et leur equivalent sans fixation test_Llzir82_0fix et test_Llzsc82_0fix
## 4 compartiments, profil sol lusignan,  meteo 1982, ZESX=39cm , CFES=3.5
## Tsol forcees comme stics
##
####################

def critN (MS, a=4.8, b=-0.33):
    """ courbe critique de dilution de l'N """
    return min(6.5, a*MS**b) #en %


## lecture fichier meteo
#faire une fonction??
meteo_path = r'H:\Travail\devel\grassland\grassland\L-gume\meteo_exemple.xls'
ongletM = 'Llzir82_0fix'#'Llzsc82_0fix'#'Llzsc82'#'Llzir82'#
met = xlrd.open_workbook(meteo_path)
meteo = IOtable.conv_dataframe(get_xls_col(met.sheet_by_name(ongletM)))
for k in ['year', 'month', 'day', 'DOY','Coupe']: meteo[k] = list(map(int, meteo[k]))




## sol
pattern8 = [[-50,-50], [50,50]]#pattern.8 equivalent (cm)
Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
surfsolref = Lsol*largsol #m2
dz_sol = [0.25, 0.30, 0.35, 0.30] #m
ncouches_sol = 4

ZESX = 0.39 #m
CFES = 3.5#coefficient de forme sans dimension [0.5 - 5.]
ACLIMc=26.#site de lusignan
ARGIs=22.#lusignan
CALCs = 1. #lusignan
Norg = 1.2 #gN.kg-1 sol (*10 de Norg STICS exprime en %)
q0 = 9. #(mm de sol)
pH= 7.#lusignan
Tsol=15. #degresC
concrr = 0.02/10000. #kg N.m-2.mm-1 (concentration en N de la pluie)

### pour calculs d'humidite

CN0_30 = 9.5238 # (1/Wh de STICS = 0.105)
MO0_30 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)


#vecteur d'initialisation
vDA = [1.16, 1.29, 1.15, 0.91] #densite apparente de sol
vCN = [CN0_30]*ncouches_sol #maxi 90cm en strates de 5cm
vMO = [MO0_30]*ncouches_sol #maxi 90cm en strates de 5cm
vARGIs = [ARGIs]*ncouches_sol #maxi 90cm
vCALCs = [CALCs]*ncouches_sol
vNH4 = [0., 0., 0., 0.] # #!! kg d'N.ha-1 (entree de STICS)
vNO3 = [32., 12., 9., 0.] # kg d'N.ha-1 (entree de STICS)
HRpinit = [22.5, 24.4, 25.8, 28.7]



#recalcule a partir du fichier sol de lusignan
par_sol = {'1':{'soil number': '1', 'soil type': 'H1 - Lusignan', 'teta_sat': 0.35, 'teta_fc': 0.261, 'teta_wp': 0.0986, 'teta_ad': 0.011,   'WCST': '0.35',   'gamma_theo': '0.07'}}
par_sol['2'] = {'soil number': '2', 'soil type': 'H2 - Lusignan', 'teta_sat': 0.41, 'teta_fc': 0.3277, 'teta_wp': 0.1496, 'teta_ad': 0.026,   'WCST': '0.41',   'gamma_theo': '0.056'}
par_sol['3'] = {'soil number': '3', 'soil type': 'H3 - Lusignan', 'teta_sat': 0.35, 'teta_fc': 0.2967, 'teta_wp': 0.1759, 'teta_ad': 0.054,   'WCST': '0.35',   'gamma_theo': '0.038'}
par_sol['4'] = {'soil number': '3', 'soil type': 'H4 - Lusignan', 'teta_sat': 0.35, 'teta_fc': 0.2611, 'teta_wp': 0.1729, 'teta_ad': 0.054,   'WCST': '0.35',   'gamma_theo': '0.038'}
#Huv = Hup*DA


#parametre sol complementaires pour l'N
par_SN = {}
par_SN['FMIN1G'] = 0.0006 #(day-1) (p145)
par_SN['FMIN2G'] = 0.0272 #(%ARGIS-1) para pot rate min a ARGIs (p145)
par_SN['FMIN3G'] = 0.0167 #(% CALC-1)para pot rate min a CALCss (p145)

par_SN['FINERTG'] = 0.65 #0.65 = default fraction of stable pool (p145) should be smaller in grassland & forest
par_SN['PROFHUMs'] = 35. # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220

par_SN['HMinMg'] = 0.3 #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
par_SN['HoptMg'] = 1. #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)

par_SN['TRefg'] = 15. #reference temperature
par_SN['FTEMHAg'] = 25. #asymptotic value of FTH (seen as a logistic response)
par_SN['FTEMHg'] = 0.120 #(K-1)
par_SN['FTEMHB'] = 145.

par_SN['FNXg'] = 0.5 #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
par_SN['PHMinNITg'] = 3. #pH min de nitrification (prop of Field Capacity) #value p149
par_SN['PHMaxNITg'] = 5.5 #pH max de nitrification (prop of Field Capacity) #value p149
par_SN['HMinNg'] = 0.67 #Humidite min de nitrification #value p149
par_SN['HoptNg'] = 1. #Humidite opt de nitrification  #value p149
par_SN['TNITMINg'] = 5. #Temperature min de nitrification #  (degreC)value p151
par_SN['TNITOPTg'] = 20. #Temperature opt de nitrification #  (degreC)value p151
par_SN['TNITMAXg'] = 45. #Temperature max de nitrification #  (degreC)value p151
par_SN['RATIONITs'] = 0. #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
par_SN['DIFNg'] = 0.018  #N diffusion coefficient at field capacity (cm2.day-1, p 161)


#correctif de profhum pour tenir compte des changements de densite entre horizons
par_SN['PROFHUMs'] = 25. + 10.* vDA[0]/vDA[1]


# parametre plante pour l'uptake d'N
ParamP = [{}]
ParamP[0]['Vmax1'] = 0.0018 #(alfalfa.plt) micromole.cm-1.h-1
ParamP[0]['Kmax1'] = 50. #(alfalfa.plt) micromole.L-1
ParamP[0]['Vmax2'] = 0.05 #(alfalfa.plt) micromole.cm-1.h-1
ParamP[0]['Kmax2'] = 25000. #(alfalfa.plt) micromole.L-1

#initialise teneur en N des plantes
Npc = 6.5


## soil initialisation
S = SoilN(par_sol, par_SN, soil_number = [1,2,3,4], dxyz = [[Lsol], [largsol], dz_sol], vDA=vDA, vCN=vCN,vMO=vMO, vARGIs = vARGIs,vNO3=vNO3,vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol,pH=pH, ZESX=ZESX , CFES=CFES)
S.init_asw(HRp_init=HRpinit)




Uval = q0*0.1*sum(S.m_QH20fc[0])*surfsolref / (S.dxyz[2][0]*100.)#(epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
stateEV = [0.,0.,0.] #pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
b= bEV(ACLIMc, ARGIs, HXs=0.261)#1.#valeur empirique tres proche#0.1#0.63#0.63


#debut, fin de simulation
DOY_deb, DOY_fin = 1,360


##vegetation sans racine ni LAI
droot = 0.5 #cm.cm-3
R1 = S.m_1 * droot* S.m_soil_vol * 100**3 / 100. #(m root.voxel-1) vert_roots(S.dxyz, [0.5,0.5,0.5,0.5]) #pas zero sinon buf FTSW -> longueur! pas densite de longueur!
ls_roots = [R1]
ls_epsi = [0.]



##boucle journaliere
vlix, azomes, resmes, demplt = [], [], [], []
for DOY in range(DOY_deb, DOY_fin):
    #meteo
    meteo_j = extract_dataframe(meteo, ['TmoyDay','I0','Et0','Precip','Irrig','Coupe','FertNO3','FertNH4','Tsol','LAI','MS','dMS', 'Npc'], 'DOY', val=DOY)
    for k in list(meteo_j.keys()): meteo_j[k]=meteo_j[k][0]

    #preparation des entrees eau
    Rain = meteo_j['Precip']
    Irrig = meteo_j['Irrig']
    ls_epsi = [1.-exp(-0.8*meteo_j['LAI'])]

    #preparation des entrees azote
    #map_N = 1.*S.m_1[0,:,:] * (Rain+Irrig) * concrr #Nmin de la pluie + irrigation
    mapN_Rain = 1.*S.m_1[0,:,:] * Rain * concrr #Nmin de la pluie
    mapN_Irrig = 1.*S.m_1[0,:,:] * Irrig * concrr #Nmin de l'eau d'irrigation
    mapN_fertNO3 = 1.*S.m_1[0,:,:] * meteo_j['FertNO3'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    mapN_fertNH4 = 1.*S.m_1[0,:,:] * meteo_j['FertNH4'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel

    S.updateTsol(meteo_j['Tsol'])#(meteo_j['TmoyDay']) #Tsol forcee comme dans STICS

    QN = meteo_j['MS'] * Npc/100.*1000#kg N.ha-1 #%N libre
    #QN = meteo_j['MS'] * meteo_j['Npc']/100.*1000#kg N.ha-1 #%N force
    PotN = (meteo_j['MS'] + meteo_j['dMS']) * critN (meteo_j['MS'] + meteo_j['dMS'])/100.*1000 #kg N.ha-1
    demande_N_plt =  max(PotN - QN, 0.) #kg N.ha-1
    ls_demandeN = [demande_N_plt/10000.] #kg N.surface de sol


    ls_transp, evapo_tot, Drainage, stateEV,  ls_m_transpi, m_evap, ls_ftsw = S.stepWBmc(meteo_j['Et0']*surfsolref, ls_roots, ls_epsi, Rain*surfsolref, Irrig*surfsolref, stateEV, ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    S.stepNB(par_SN)
    S.stepNitrif(par_SN)

    ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N, idmin = S.stepNuptakePlt(par_SN, ParamP, ls_roots, ls_m_transpi, ls_demandeN)
    #PotUpNtot, ls_Pot_Nuptake_plt, idmin = Distrib_Potential_Nuptake_Plt(S, par_SN, ParamP, ls_roots, ls_m_transpi)
    #ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N = Actual_Nuptake_plt(S, ls_Pot_Nuptake_plt, ls_demandeN)

    #S.stepNINFILT(map_N, Drainage, opt=1)
    S.stepNINFILT(mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, Drainage, opt=1)

    #mise a jour de Npc
    QN = QN + sum(ActUpNtot)/S.soilSurface()*10000 #kg N.ha-1
    Npc = (QN/1000.)/(meteo_j['MS']+meteo_j['dMS'])*100 #%

    #sorties
    print(DOY)
    vlix.append(S.lixiNO3*10000)
    azomes.append(sum3(S.m_NH4+S.m_NO3)/S.soilSurface()*10000)
    resmes.append(sum(S.tsw_t))
    demplt.append(sum(ActUpNtot)/S.soilSurface()*10000)
    #demplt.append(Npc)

S.CloseWbalance()
S.CloseNbalance()
#r.plot(range(DOY_deb, DOY_fin), demplt, ylab='uptake plt', xlab='DOY', type='l')

##ecriture sorties
#f = file(r'H:\simul\Valid_sol\Lusig0fix_.csv', 'w')
#IOtable.ecriture_csv(IOtable.t_list([azomes,resmes,demplt]), f)
#f.close()











#####################
## Test Aschyd- sol nu
## 30 * 9 compartiments ZESX=60cm , CFES=3.5
## Tsol forcees comme mesuree
##
####################

## lecture fichier meteo
#faire une fonction??
meteo_path = r'H:\Travail\devel\grassland\grassland\L-gume\meteo_exemple.xls'
ongletM = 'ASCHYD11'
met = xlrd.open_workbook(meteo_path)
meteo = IOtable.conv_dataframe(get_xls_col(met.sheet_by_name(ongletM)))
for k in ['year', 'month', 'day', 'DOY','Coupe']: meteo[k] = list(map(int, meteo[k]))




## sol
pattern8 = [[-49.5/2.,-22.3/2.], [49.5/2.,22.3/2.]]#pattern.8 equivalent (cm)
Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
surfsolref = Lsol*largsol #m2
dz_sol = 0.05 #m
ncouches_sol = 30

ZESX = 0.60 #m
CFES = 3.5#coefficient de forme sans dimension [0.5 - 5.]
ACLIMc= 3.#serre
ARGIs=22.#lusignan
CALCs = 1. #lusignan
Norg = 0.92 #gN.kg-1 sol (moyenne Analyse sol)
q0 = 0. #(mm de sol)
pH= 7.#lusignan
Tsol=15. #degresC
concrr = 0.02/10000. #kg N.m-2.mm-1 (concentration en N de la pluie)

### pour calculs d'humidite

CN0_30 = 11.97 # moyenne Analyse sol
MO0_30 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)


#vecteur d'initialisation
vDA = [1.81]+[1.31]*3+[1.37]*13+[1.42]*13 #densite apparente de sol (mesure pesees initial aschyd11)
vCN = [CN0_30]*ncouches_sol #maxi 90cm en strates de 5cm
vMO = [MO0_30]*ncouches_sol #maxi 90cm en strates de 5cm
vARGIs = [ARGIs]*ncouches_sol #maxi 90cm
vCALCs = [CALCs]*ncouches_sol
vNH4 = [2.]*ncouches_sol # #!! kg d'N.ha-1 (entree de STICS)

coeff = 0.15#coeff perte ressuyage -> a ajuster pour avoir environ 600 kg N.ha-1
vNO3 = [91.*coeff]*ncouches_sol # kg d'N.ha-1 (entree de STICS)
HRpinit = [25.5,26.,25.,25.5,26.,26.,26.,26.5,26.5,27.,27.,27.,27.5,27.5,27.5,27.5,27.5,29,29,29,29,29,29,29,29,30,30,30,30,30]#-> mesures ahscyd au jour 195 (140711) -> init sol nu

vsoilnumbers = [1]+[2]*3+[3]*13+[4]*13 #numeros de sol du profil -> mesures acsyd11



#recalcule a partir du fichier sol de lusignan
#Huv = Hup*DA
par_sol = {'1':{'soil number': '1', 'soil type': 'H5cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.479, 'teta_wp': 0.380, 'teta_ad': 0.340,   'WCST': '0.35',   'gamma_theo': '0.07'}}
par_sol['2'] = {'soil number': '2', 'soil type': 'H10cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.347, 'teta_wp': 0.275, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}
par_sol['3'] = {'soil number': '3', 'soil type': 'H25cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.377, 'teta_wp': 0.301, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}
par_sol['4'] = {'soil number': '4', 'soil type': 'H90cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.426, 'teta_wp': 0.312, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}


#parametre sol complementaires pour l'N
par_SN = {}
par_SN['FMIN1G'] = 0.0006 #(day-1) (p145)
par_SN['FMIN2G'] = 0.0272 #(%ARGIS-1) para pot rate min a ARGIs (p145)
par_SN['FMIN3G'] = 0.0167 #(% CALC-1)para pot rate min a CALCss (p145)

par_SN['FINERTG'] = 0.65 #0.65 = default fraction of stable pool (p145) should be smaller in grassland & forest
par_SN['PROFHUMs'] = 150. # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220

par_SN['HMinMg'] = 0.3 #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
par_SN['HoptMg'] = 1. #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)

par_SN['TRefg'] = 15. #reference temperature
par_SN['FTEMHAg'] = 25. #asymptotic value of FTH (seen as a logistic response)
par_SN['FTEMHg'] = 0.120 #(K-1)
par_SN['FTEMHB'] = 145.

par_SN['FNXg'] = 0.5 #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
par_SN['PHMinNITg'] = 3. #pH min de nitrification (prop of Field Capacity) #value p149
par_SN['PHMaxNITg'] = 5.5 #pH max de nitrification (prop of Field Capacity) #value p149
par_SN['HMinNg'] = 0.67 #Humidite min de nitrification #value p149
par_SN['HoptNg'] = 1. #Humidite opt de nitrification  #value p149
par_SN['TNITMINg'] = 5. #Temperature min de nitrification #  (degreC)value p151
par_SN['TNITOPTg'] = 20. #Temperature opt de nitrification #  (degreC)value p151
par_SN['TNITMAXg'] = 45. #Temperature max de nitrification #  (degreC)value p151
par_SN['RATIONITs'] = 0. #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
par_SN['DIFNg'] = 0.018  #N diffusion coefficient at field capacity (cm2.day-1, p 161)


#correctif de profhum pour tenir compte des changements de densite entre horizons
##par_SN['PROFHUMs'] = 25. + 10.* vDA[0]/vDA[1]


# parametre plante pour l'uptake d'N
ParamP = [{}]
ParamP[0]['Vmax1'] = 0.0018 #(alfalfa.plt) micromole.cm-1.h-1
ParamP[0]['Kmax1'] = 50. #(alfalfa.plt) micromole.L-1
ParamP[0]['Vmax2'] = 0.05 #(alfalfa.plt) micromole.cm-1.h-1
ParamP[0]['Kmax2'] = 25000. #(alfalfa.plt) micromole.L-1

#initialise teneur en N des plantes
Npc = 6.5


## soil initialisation
S = SoilN(par_sol, par_SN, soil_number = vsoilnumbers, dxyz = [[Lsol/9.]*9, [largsol], [dz_sol]*ncouches_sol], vDA=vDA, vCN=vCN,vMO=vMO, vARGIs = vARGIs,vNO3=vNO3,vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol,pH=pH, ZESX=ZESX , CFES=CFES)
S.init_asw(HRp_init=HRpinit)




Uval = q0*0.1*sum(S.m_QH20fc[0])*surfsolref / (S.dxyz[2][0]*100.)#(epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
stateEV = [0.,0.,0.] #pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
b= bEV(ACLIMc, ARGIs, HXs=0.261)#1.#valeur empirique tres proche#0.1#0.63#0.63



#lecture et construction du systeme racinaire
rac_path1 = r'H:\Travail\Lusignan\Analyses\ASCHYD11\matrices densites racinaires apparente RH1.xls'
rac_path2 = r'H:\Travail\Lusignan\Analyses\ASCHYD11\matrices densites racinaires apparente RH2.xls'

ls_datesR = ['080711','190711','290711','030811','050811','120811','190811','250811','290811','020911','050912']
onglet = '080711'#ls_datesR[9]#'050811'#'190811'#'050911'#


def read_rootmap(rac_path1, onglet, S):
    """ lit fichier excell contenant une map de densite de racines et le met dans sol 2D """
    #rac_path1 = r'H:\Lusignan\Analyses\ASCHYD11\matrices densites racinaires apparente RH1.xls'

    rac1 = xlrd.open_workbook(rac_path1)
    rac_i = get_xls_row(rac1.sheet_by_name(onglet))
    #remet dans sol au connes dimensions
    R1 = S.m_1*0.
    for i in range(len(rac_i)):
        R1[i,:,0] = array(rac_i[i])

    return R1
    #generaliser pour 3D?

R1 = read_rootmap(rac_path1, onglet, S)
R2 = read_rootmap(rac_path2, onglet, S)
Rcum = R1+R2 #somme des deux systemes
ls_roots = [Rcum]
ls_epsi = [0.]




#plot

#Mascene = S.plot_soil_properties(1.-R1/R1.max(), Scene(), 4)# #densite relative de racines
Mascene = plot_soil_properties(S, 1.-Rcum/Rcum.max(), Scene(), 4)# #densite relative
Monviewer = Viewer
Monviewer.display(Mascene)

Mascene = plot_soil_properties(S, S.ftsw_t, Scene(), 5)# ftsw
Monviewer = Viewer
Monviewer.display(Mascene)




#verif les numeros de date de RH2 -> pas forcement meme decalage que Rh1! (date 19?)

#calcul des humidite journalieres par horizon et des DA
#gerer les irrigations

#A faire -> calage des paramtres de sol (borne humidite



#debut, fin de simulation
DOY_deb, DOY_fin = 195,251 #187-251 -> 1ere pousse; 195=initialisation sol

#id couches sorties
id_out = [0,1,4,11,17,25]# 5,10,25,60,90,130 cm, id de voxel dans le sol
out_HR = [['DOY','HP5', 'HP10', 'HP25', 'HP60', 'HP90', 'HP130']]

##boucle journaliere
for DOY in range(DOY_deb, DOY_fin):
    #meteo
    meteo_j = extract_dataframe(meteo, ['TmoyDay','I0','Et0','Precip','Irrig','Coupe','FertNO3','FertNH4','irrig_Rh1N','systRac','ETPserre_cor'], 'DOY', val=DOY)
    for k in list(meteo_j.keys()): meteo_j[k]=meteo_j[k][0]

    #preparation du ls_roots
    onglet = str(meteo_j['systRac'])
    R1 = read_rootmap(rac_path1, onglet, S)
    R2 = read_rootmap(rac_path2, onglet, S)
    ls_roots = [R1+R2]#1 seul qui est la somme des deux systemes


    #preparation des entrees eau
    Rain = meteo_j['Precip']
    Irrig = meteo_j['irrig_Rh1N']#R1N = sol_nu
    ##ls_epsi = [1.-exp(-0.8*meteo_j['LAI'])]

    #preparation des entrees azote
    mapN_Rain = 1.*S.m_1[0,:,:] * Rain * concrr #Nmin de la pluie
    mapN_Irrig = 1.*S.m_1[0,:,:] * Irrig * concrr #Nmin de l'eau d'irrigation
    mapN_fertNO3 = 1.*S.m_1[0,:,:] * meteo_j['FertNO3'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    mapN_fertNH4 = 1.*S.m_1[0,:,:] * meteo_j['FertNH4'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    S.updateTsol(meteo_j['TmoyDay'])#(meteo_j['Tsol'])# #Tsol forcee comme dans STICS

    ##QN = meteo_j['MS'] * Npc/100.*1000#kg N.ha-1 #%N libre
    ##PotN = (meteo_j['MS'] + meteo_j['dMS']) * critN (meteo_j['MS'] + meteo_j['dMS'])/100.*1000 #kg N.ha-1
    ##demande_N_plt =  max(PotN - QN, 0.) #kg N.ha-1
    ##ls_demandeN = [demande_N_plt/10000.] #kg N.surface de sol


    ls_transp, evapo_tot, Drainage, stateEV,  ls_m_transpi, m_evap, ls_ftsw = S.stepWBmc(meteo_j['ETPserre_cor']*surfsolref, ls_roots, ls_epsi, Rain*surfsolref, Irrig*surfsolref, stateEV, ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    S.stepNB(par_SN)
    S.stepNitrif(par_SN)

    ##PotUpNtot, ls_Pot_Nuptake_plt, idmin = Distrib_Potential_Nuptake_Plt(S, par_SN, ParamP, ls_roots, ls_m_transpi)
    ##ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N = Actual_Nuptake_plt(S, ls_Pot_Nuptake_plt, ls_demandeN)

    S.stepNINFILT(mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, Drainage, opt=1)

    #sorties
    print(DOY)
    out_HR.append([DOY]+mean(S.HRp(), axis=1)[id_out,0].tolist())

S.CloseWbalance()
S.CloseNbalance()
#fichier contenant les profils d'humidite
#f = file(r'H:\simul\Valid_sol\out_HR.csv', 'w')
#IOtable.ecriture_csv(out_HR, f)
#f.close()







#####################
## Test Aschyd- sol avec racine
## 30 * 9 compartiments ZESX=60cm , CFES=3.5
## Tsol forcees comme mesuree
##
####################


## lecture fichier meteo
#faire une fonction??
meteo_path = r'H:\Travail\devel\grassland\grassland\L-gume\meteo_exemple.xls'
ongletM = 'ASCHYD11'
met = xlrd.open_workbook(meteo_path)
meteo = IOtable.conv_dataframe(get_xls_col(met.sheet_by_name(ongletM)))
for k in ['year', 'month', 'day', 'DOY','Coupe']: meteo[k] = list(map(int, meteo[k]))




## sol
pattern8 = [[-49.5/2.,-22.3/2.], [49.5/2.,22.3/2.]]#pattern.8 equivalent (cm)
Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
surfsolref = Lsol*largsol #m2
dz_sol = 0.05 #m
ncouches_sol = 30

ZESX = 0.60 #m
CFES = 3.5#coefficient de forme sans dimension [0.5 - 5.]
ACLIMc= 3.#serre
ARGIs=22.#lusignan
CALCs = 1. #lusignan
Norg = 0.92 #gN.kg-1 sol (moyenne Analyse sol)
q0 = 0. #(mm de sol)
pH= 7.#lusignan
Tsol=15. #degresC
concrr = 0.02/10000. #kg N.m-2.mm-1 (concentration en N de la pluie)

### pour calculs d'humidite

CN0_30 = 11.97 # moyenne Analyse sol
MO0_30 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)


#vecteur d'initialisation
vDA = [1.81]+[1.31]*3+[1.37]*13+[1.42]*13 #densite apparente de sol (mesure pesees initial aschyd11)
vCN = [CN0_30]*ncouches_sol #maxi 90cm en strates de 5cm
vMO = [MO0_30]*ncouches_sol #maxi 90cm en strates de 5cm
vARGIs = [ARGIs]*ncouches_sol #maxi 90cm
vCALCs = [CALCs]*ncouches_sol
vNH4 = [2.]*ncouches_sol # #!! kg d'N.ha-1 (entree de STICS)

coeff = 0.15#coeff perte ressuyage -> a ajuster pour avoir environ 600 kg N.ha-1
vNO3 = [91.*coeff]*ncouches_sol # kg d'N.ha-1 (entree de STICS)
HRpinit = [25.5,26.,25.,25.5,26.,26.,26.,26.5,26.5,27.,27.,27.,27.5,27.5,27.5,27.5,27.5,29,29,29,29,29,29,29,29,30,30,30,30,30]#-> mesures ahscyd au jour 195 (140711) -> init sol nu

vsoilnumbers = [1]+[2]*3+[3]*13+[4]*13 #numeros de sol du profil -> mesures acsyd11



#recalcule a partir du fichier sol de lusignan
#Huv = Hup*DA
par_sol = {'1':{'soil number': '1', 'soil type': 'H5cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.479, 'teta_wp': 0.380, 'teta_ad': 0.340,   'WCST': '0.35',   'gamma_theo': '0.07'}}
par_sol['2'] = {'soil number': '2', 'soil type': 'H10cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.347, 'teta_wp': 0.275, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}
par_sol['3'] = {'soil number': '3', 'soil type': 'H25cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.377, 'teta_wp': 0.301, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}
par_sol['4'] = {'soil number': '4', 'soil type': 'H90cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.426, 'teta_wp': 0.312, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}


#parametre sol complementaires pour l'N
par_SN = {}
par_SN['FMIN1G'] = 0.0006 #(day-1) (p145)
par_SN['FMIN2G'] = 0.0272 #(%ARGIS-1) para pot rate min a ARGIs (p145)
par_SN['FMIN3G'] = 0.0167 #(% CALC-1)para pot rate min a CALCss (p145)

par_SN['FINERTG'] = 0.65 #0.65 = default fraction of stable pool (p145) should be smaller in grassland & forest
par_SN['PROFHUMs'] = 150. # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220

par_SN['HMinMg'] = 0.3 #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
par_SN['HoptMg'] = 1. #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)

par_SN['TRefg'] = 15. #reference temperature
par_SN['FTEMHAg'] = 25. #asymptotic value of FTH (seen as a logistic response)
par_SN['FTEMHg'] = 0.120 #(K-1)
par_SN['FTEMHB'] = 145.

par_SN['FNXg'] = 0.5 #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
par_SN['PHMinNITg'] = 3. #pH min de nitrification (prop of Field Capacity) #value p149
par_SN['PHMaxNITg'] = 5.5 #pH max de nitrification (prop of Field Capacity) #value p149
par_SN['HMinNg'] = 0.67 #Humidite min de nitrification #value p149
par_SN['HoptNg'] = 1. #Humidite opt de nitrification  #value p149
par_SN['TNITMINg'] = 5. #Temperature min de nitrification #  (degreC)value p151
par_SN['TNITOPTg'] = 20. #Temperature opt de nitrification #  (degreC)value p151
par_SN['TNITMAXg'] = 45. #Temperature max de nitrification #  (degreC)value p151
par_SN['RATIONITs'] = 0. #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
par_SN['DIFNg'] = 0.018  #N diffusion coefficient at field capacity (cm2.day-1, p 161)


#correctif de profhum pour tenir compte des changements de densite entre horizons
##par_SN['PROFHUMs'] = 25. + 10.* vDA[0]/vDA[1]


# parametre plante pour l'uptake d'N
ParamP = [{}]
ParamP[0]['Vmax1'] = 0.0018 #(alfalfa.plt) micromole.cm-1.h-1
ParamP[0]['Kmax1'] = 50. #(alfalfa.plt) micromole.L-1
ParamP[0]['Vmax2'] = 0.05 #(alfalfa.plt) micromole.cm-1.h-1
ParamP[0]['Kmax2'] = 25000. #(alfalfa.plt) micromole.L-1

#initialise teneur en N des plantes
Npc = 6.5


## soil initialisation
S = SoilN(par_sol, par_SN, soil_number = vsoilnumbers, dxyz = [[Lsol/9.]*9, [largsol], [dz_sol]*ncouches_sol], vDA=vDA, vCN=vCN,vMO=vMO, vARGIs = vARGIs,vNO3=vNO3,vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol,pH=pH, ZESX=ZESX , CFES=CFES)
S.init_asw(HRp_init=HRpinit)




Uval = q0*0.1*sum(S.m_QH20fc[0])*surfsolref / (S.dxyz[2][0]*100.)#(epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
stateEV = [0.,0.,0.] #pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
b= bEV(ACLIMc, ARGIs, HXs=0.261)#1.#valeur empirique tres proche#0.1#0.63#0.63



#lecture et construction du systeme racinaire
rac_path1 = r'H:\Travail\Lusignan\Analyses\ASCHYD11\matrices densites racinaires apparente RH1.xls'
rac_path2 = r'H:\Travail\Lusignan\Analyses\ASCHYD11\matrices densites racinaires apparente RH2.xls'

ls_datesR = ['080711','190711','290711','030811','050811','120811','190811','250811','290811','020911','050912']
onglet = '080711'#ls_datesR[9]#'050811'#'190811'#'050911'#


def read_rootmap(rac_path1, onglet, S):
    """ lit fichier excell contenant une map de densite de racines et le met dans sol 2D """
    #rac_path1 = r'H:\Lusignan\Analyses\ASCHYD11\matrices densites racinaires apparente RH1.xls'

    rac1 = xlrd.open_workbook(rac_path1)
    rac_i = get_xls_row(rac1.sheet_by_name(onglet))
    #remet dans sol au connes dimensions
    R1 = S.m_1*0.
    for i in range(len(rac_i)):
        R1[i,:,0] = array(rac_i[i])

    return R1
    #generaliser pour 3D?

R1 = read_rootmap(rac_path1, onglet, S)
R2 = read_rootmap(rac_path2, onglet, S)
Rcum = R1+R2 #somme des deux systemes
ls_roots = [Rcum]
ls_epsi = [0.]




#plot

#Mascene = S.plot_soil_properties(1.-R1/R1.max(), Scene(), 4)# #densite relative de racines
Mascene = plot_soil_properties(S, 1.-Rcum/Rcum.max(), Scene(), 4)# #densite relative
Monviewer = Viewer
Monviewer.display(Mascene)

Mascene = plot_soil_properties(S, S.ftsw_t, Scene(), 5)# ftsw
Monviewer = Viewer
Monviewer.display(Mascene)




#verif les numeros de date de RH2 -> pas forcement meme decalage que Rh1! (date 19?)

#calcul des humidite journalieres par horizon et des DA
#gerer les irrigations

#A faire -> calage des paramtres de sol (borne humidite



#debut, fin de simulation
DOY_deb, DOY_fin = 195,251 #187-251 -> 1ere pousse; 195=initialisation sol

#id couches sorties
id_out = [0,1,4,11,17,25]# 5,10,25,60,90,130 cm, id de voxel dans le sol
out_HR = [['DOY','HP5', 'HP10', 'HP25', 'HP60', 'HP90', 'HP130']]

##boucle journaliere
for DOY in range(DOY_deb, DOY_fin):
    #meteo
    meteo_j = extract_dataframe(meteo, ['TmoyDay','I0','Et0','Precip','Irrig','Coupe','FertNO3','FertNH4','irrig_Rh1S','systRac','ETPserre_cor','LAI'], 'DOY', val=DOY)
    for k in list(meteo_j.keys()): meteo_j[k]=meteo_j[k][0]

    #preparation du ls_roots
    onglet = str(meteo_j['systRac'])
    R1 = read_rootmap(rac_path1, onglet, S)
    R2 = read_rootmap(rac_path2, onglet, S)
    corpiv = S.m_1*0.
    corpiv[0,:,0] = sum(R1+R2)*0.8/9.
    #corpiv[1,:,0] = sum(R1+R2)*0.1/9.
    ls_roots = [R1+R2+corpiv]#1 seul qui est la somme des deux systemes


    #preparation des entrees eau
    Rain = meteo_j['Precip']
    Irrig = meteo_j['irrig_Rh1S']#R1S = asso
    ls_epsi = [1.-exp(-0.8*meteo_j['LAI'])]

    #preparation des entrees azote
    mapN_Rain = 1.*S.m_1[0,:,:] * Rain * concrr #Nmin de la pluie
    mapN_Irrig = 1.*S.m_1[0,:,:] * Irrig * concrr #Nmin de l'eau d'irrigation
    mapN_fertNO3 = 1.*S.m_1[0,:,:] * meteo_j['FertNO3'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    mapN_fertNH4 = 1.*S.m_1[0,:,:] * meteo_j['FertNH4'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    S.updateTsol(meteo_j['TmoyDay'])#(meteo_j['Tsol'])# #Tsol forcee comme dans STICS

    ##QN = meteo_j['MS'] * Npc/100.*1000#kg N.ha-1 #%N libre
    ##PotN = (meteo_j['MS'] + meteo_j['dMS']) * critN (meteo_j['MS'] + meteo_j['dMS'])/100.*1000 #kg N.ha-1
    ##demande_N_plt =  max(PotN - QN, 0.) #kg N.ha-1
    ##ls_demandeN = [demande_N_plt/10000.] #kg N.surface de sol


    ls_transp, evapo_tot, Drainage, stateEV,  ls_m_transpi, m_evap, ls_ftsw = S.stepWBmc(meteo_j['ETPserre_cor']*surfsolref, ls_roots, ls_epsi, Rain*surfsolref, Irrig*surfsolref, stateEV, ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    S.stepNB(par_SN)
    S.stepNitrif(par_SN)

    ##PotUpNtot, ls_Pot_Nuptake_plt, idmin = Distrib_Potential_Nuptake_Plt(S, par_SN, ParamP, ls_roots, ls_m_transpi)
    ##ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N = Actual_Nuptake_plt(S, ls_Pot_Nuptake_plt, ls_demandeN)

    S.stepNINFILT(mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, Drainage, opt=1)

    #sorties
    print(DOY)
    out_HR.append([DOY]+mean(S.HRp(), axis=1)[id_out,0].tolist())

    Mascene = plot_soil_properties(S, S.ftsw_t, Scene(), 5)# ftsw
    Monviewer = Viewer
    Monviewer.display(Mascene)


S.CloseWbalance()
S.CloseNbalance()
#fichier contenant les profils d'humidite
#f = file(r'H:\simul\Valid_sol\out_HR2.csv', 'w')
#IOtable.ecriture_csv(out_HR, f)
#f.close()







#####################
## Test Aschyd- sol nu
## 30 * 9 compartiments ZESX=60cm , CFES=3.5
## Avec ajout de residus en continu
##
####################

## lecture fichier meteo
#faire une fonction??
meteo_path = r'H:\Travail\devel\grassland\grassland\L-gume\meteo_exemple.xls'
ongletM = 'ASCHYD11'
met = xlrd.open_workbook(meteo_path)
meteo = IOtable.conv_dataframe(get_xls_col(met.sheet_by_name(ongletM)))
for k in ['year', 'month', 'day', 'DOY','Coupe']: meteo[k] = list(map(int, meteo[k]))




## sol
pattern8 = [[-49.5/2.,-22.3/2.], [49.5/2.,22.3/2.]]#pattern.8 equivalent (cm)
Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
surfsolref = Lsol*largsol #m2
dz_sol = 0.05 #m
ncouches_sol = 30

ZESX = 0.60 #m
CFES = 3.5#coefficient de forme sans dimension [0.5 - 5.]
ACLIMc= 3.#serre
ARGIs=22.#lusignan
CALCs = 1. #lusignan
Norg = 0.92 #gN.kg-1 sol (moyenne Analyse sol)
q0 = 0. #(mm de sol)
pH= 7.#lusignan
Tsol=15. #degresC
concrr = 0.02/10000. #kg N.m-2.mm-1 (concentration en N de la pluie)

### pour calculs d'humidite

CN0_30 = 11.97 # moyenne Analyse sol
MO0_30 = Norg * CN0_30 * 1.72 #gC.kg-1 sol (1.72: facteur de passage de Corg a MO)


#vecteur d'initialisation
vDA = [1.81]+[1.31]*3+[1.37]*13+[1.42]*13 #densite apparente de sol (mesure pesees initial aschyd11)
vCN = [CN0_30]*ncouches_sol #maxi 90cm en strates de 5cm
vMO = [MO0_30]*ncouches_sol #maxi 90cm en strates de 5cm
vARGIs = [ARGIs]*ncouches_sol #maxi 90cm
vCALCs = [CALCs]*ncouches_sol
vNH4 = [2.]*ncouches_sol # #!! kg d'N.ha-1 (entree de STICS)

coeff = 0.15#coeff perte ressuyage -> a ajuster pour avoir environ 600 kg N.ha-1
vNO3 = [91.*coeff]*ncouches_sol # kg d'N.ha-1 (entree de STICS)
HRpinit = [25.5,26.,25.,25.5,26.,26.,26.,26.5,26.5,27.,27.,27.,27.5,27.5,27.5,27.5,27.5,29,29,29,29,29,29,29,29,30,30,30,30,30]#-> mesures ahscyd au jour 195 (140711) -> init sol nu

vsoilnumbers = [1]+[2]*3+[3]*13+[4]*13 #numeros de sol du profil -> mesures acsyd11



#1 residu = listes de 1 element
vAmount= [20.]# T Fresh Weight.ha-1 (equivalent QRES)
vCNRESt = [16.] #g.g-1 (equivalent CsurNres)
Vprop1 = [1./3., 1./3., 1./3.]+27*[0.] #distribution dans les horizons
vProps= [Vprop1]#[Vprop1]#[Vprop1, Vprop1, Vprop1]
vWC=[0.7]# fraction d'eau des residu frais (equivalent de Crespc /(%)/100)
vCC=[0.42]# fraction de C des residus sec (equivalent de Crespc /(%)/100)
vNmires = [0.00197]# fraction de poids frais residu en azote mineral (equivalent de Nminres(%)/100)




#recalcule a partir du fichier sol de lusignan
#Huv = Hup*DA
par_sol = {'1':{'soil number': '1', 'soil type': 'H5cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.479, 'teta_wp': 0.380, 'teta_ad': 0.340,   'WCST': '0.35',   'gamma_theo': '0.07'}}
par_sol['2'] = {'soil number': '2', 'soil type': 'H10cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.347, 'teta_wp': 0.275, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}
par_sol['3'] = {'soil number': '3', 'soil type': 'H25cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.377, 'teta_wp': 0.301, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}
par_sol['4'] = {'soil number': '4', 'soil type': 'H90cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.426, 'teta_wp': 0.312, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}


#parametre sol complementaires pour l'N
par_SN = {}
par_SN['FMIN1G'] = 0.0006 #(day-1) (p145)
par_SN['FMIN2G'] = 0.0272 #(%ARGIS-1) para pot rate min a ARGIs (p145)
par_SN['FMIN3G'] = 0.0167 #(% CALC-1)para pot rate min a CALCss (p145)

par_SN['FINERTG'] = 0.65 #0.65 = default fraction of stable pool (p145) should be smaller in grassland & forest
par_SN['PROFHUMs'] = 150. # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220

par_SN['HMinMg'] = 0.3 #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
par_SN['HoptMg'] = 1. #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)

par_SN['TRefg'] = 15. #reference temperature
par_SN['FTEMHAg'] = 25. #asymptotic value of FTH (seen as a logistic response)
par_SN['FTEMHg'] = 0.120 #(K-1)
par_SN['FTEMHB'] = 145.

par_SN['FNXg'] = 0.5 #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
par_SN['PHMinNITg'] = 3. #pH min de nitrification (prop of Field Capacity) #value p149
par_SN['PHMaxNITg'] = 5.5 #pH max de nitrification (prop of Field Capacity) #value p149
par_SN['HMinNg'] = 0.67 #Humidite min de nitrification #value p149
par_SN['HoptNg'] = 1. #Humidite opt de nitrification  #value p149
par_SN['TNITMINg'] = 5. #Temperature min de nitrification #  (degreC)value p151
par_SN['TNITOPTg'] = 20. #Temperature opt de nitrification #  (degreC)value p151
par_SN['TNITMAXg'] = 45. #Temperature max de nitrification #  (degreC)value p151
par_SN['RATIONITs'] = 0. #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
par_SN['DIFNg'] = 0.018  #N diffusion coefficient at field capacity (cm2.day-1, p 161)


#correctif de profhum pour tenir compte des changements de densite entre horizons
##par_SN['PROFHUMs'] = 25. + 10.* vDA[0]/vDA[1]


# parametre plante pour l'uptake d'N
ParamP = [{}]
ParamP[0]['Vmax1'] = 0.0018 #(alfalfa.plt) micromole.cm-1.h-1
ParamP[0]['Kmax1'] = 50. #(alfalfa.plt) micromole.L-1
ParamP[0]['Vmax2'] = 0.05 #(alfalfa.plt) micromole.cm-1.h-1
ParamP[0]['Kmax2'] = 25000. #(alfalfa.plt) micromole.L-1

#initialise teneur en N des plantes
Npc = 6.5


## soil initialisation
S = SoilN(par_sol, par_SN, soil_number = vsoilnumbers, dxyz = [[Lsol/9.]*9, [largsol], [dz_sol]*ncouches_sol], vDA=vDA, vCN=vCN,vMO=vMO, vARGIs = vARGIs,vNO3=vNO3,vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol,pH=pH, ZESX=ZESX , CFES=CFES)
S.init_asw(HRp_init=HRpinit)
#S.init_residues(vCNRESt, vAmount, vProps, vWC, vCC)



Uval = q0*0.1*sum(S.m_QH20fc[0])*surfsolref / (S.dxyz[2][0]*100.)#(epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
stateEV = [0.,0.,0.] #pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
b= bEV(ACLIMc, ARGIs, HXs=0.261)#1.#valeur empirique tres proche#0.1#0.63#0.63


#debut, fin de simulation
DOY_deb, DOY_fin = 195,400 #187-251 -> 1ere pousse; 195=initialisation sol
DOYres = 300 #jour d'ajout des residus
lsDOYaj = [330,360,390] #jour  ou de aouts supplementaires du meme reidu


##vegetation sans racine ni LAI
R1 = vert_roots(S.dxyz, [0.000000001,0.,0.,0.,0.,0.,0.,0.,0.]+[0.]*21) #pas zero sinon buf FTSW
ls_roots = [R1]
ls_epsi = [0.]


##bouble journaliere
cumEV,cumET0, cumPP, cumD, profH20, cumTransp = [],[],[], [], [], []
vlix, azomes = [], []
for DOY in range(DOY_deb, DOY_fin):
    meteo_j = extract_dataframe(meteo, ['TmoyDay','I0','Et0','Precip','Irrig','Coupe','FertNO3','FertNH4','Tsol'], 'DOY', val=DOY)
    for k in list(meteo_j.keys()): meteo_j[k]=meteo_j[k][0]

    if DOY == DOYres:
        S.init_residues(vCNRESt, vAmount, vProps, vWC, vCC)

    if DOY in lsDOYaj:#ajout des residus dans pools existats
        mat_res = 0. * S.m_1
        mat_res[0,:,:] = 10.
        S.mixResMat(mat_res, 0, vCC[0])
        print ('add res')

    #entrees eau
    Rain = meteo_j['Precip']
    Irrig = meteo_j['Irrig']

    #entrees N
    #map_N = 0.*S.m_1[0,:,:]
    mapN_Rain = 1.*S.m_1[0,:,:] * Rain * concrr #Nmin de la pluie
    mapN_Irrig = 1.*S.m_1[0,:,:] * Irrig * concrr #Nmin de l'eau d'irrigation
    mapN_fertNO3 = 1.*S.m_1[0,:,:] * meteo_j['FertNO3'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel
    mapN_fertNH4 = 1.*S.m_1[0,:,:] * meteo_j['FertNH4'] * S.m_vox_surf[0,:,:]/10000. #kg N par voxel

    S.updateTsol(meteo_j['Tsol'])#(meteo_j['TmoyDay']) #Tsol forcee comme dans STICS
    ls_transp, evapo_tot, Drainage, stateEV,  m_frac_transpi, m_frac_evap, ls_ftsw = S.stepWBmc(meteo_j['Et0']*surfsolref, ls_roots, ls_epsi, Rain*surfsolref, Irrig*surfsolref, stateEV, ZESX, leafAlbedo=0.15, U=Uval, b=b, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1)
    HAx = S.HRp()
    S.stepNB(par_SN)
    if DOY>=DOYres:
        S.stepResidueMin(par_SN)
        S.stepMicrobioMin(par_SN)

    S.stepNitrif(par_SN)
    S.stepNINFILT(mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, Drainage, opt=1)

    #sorties
    print(DOY)


S.CloseWbalance() #-> equilibre
S.CloseCbalance() #-> equilibre
S.CloseNbalance() #-> equilibre


#bilan equilibre mais valeur negative pour Resid. Mineralisation dans mineral N balance
#generation de sorie de mineralisation par type de residu 'NminfromNres'
#sum(S.bilanN['cumNRes1']), S.bilanN['NminfromNresCum'],S.bilanN['NminfromNres'] #pour un seul residu
