
from scipy import *
from soil3ds import soil_moduleN as solN




def init_sol_test(pattern8 = [[-50.,-50.], [50.,50.]], dz=5., size=[10,10,30], ):
    """ ceation d'un sol test (manip Ashyd)"""

    ## sol
    #pattern8 = [[-50.,-50.], [50.,50.]]#pattern.8 equivalent (cm)
    Lsol, largsol = (pattern8[1][0]-pattern8[0][0])/100., (pattern8[1][1]-pattern8[0][1])/100. #m
    surfsolref = Lsol*largsol #m2
    dz_sol = dz/100.#0.05 #m
    ncouches_sol = size[2]#30

    ZESX = 0.30 #m
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
    vNO3 = [10.]*ncouches_sol # kg d'N.ha-1 (entree de STICS)
    HRpinit = [25.5,26.,25.,25.5,26.,26.,26.,26.5,26.5,27.,27.,27.,27.5,27.5,27.5,27.5,27.5,29,29,29,29,29,29,29,29,30,30,30,30,30]#-> mesures ahscyd au jour 195 (140711) -> init sol nu

    vsoilnumbers = [1]+[2]*3+[3]*13+[4]*13 #numeros de sol du profil -> mesures acsyd11



    #recalcule a partir du fichier sol de lusignan
    #Huv = Hup*DA
    par_sol = {'1':{'soil number': '1', 'soil type': 'H5cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.479, 'teta_wp': 0.380, 'teta_ad': 0.340,   'WCST': '0.35',   'gamma_theo': '0.07'}}
    par_sol['2'] = {'soil number': '2', 'soil type': 'H10cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.347, 'teta_wp': 0.275, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}
    par_sol['3'] = {'soil number': '3', 'soil type': 'H25cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.377, 'teta_wp': 0.301, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}
    par_sol['4'] = {'soil number': '4', 'soil type': 'H90cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.426, 'teta_wp': 0.312, 'teta_ad': 0.248,   'WCST': '0.35',   'gamma_theo': '0.07'}


    #parametre sol complementaires pour l'N
    par_SN = solN.default_parSN()

    # par_SN = {}
    # par_SN['FMIN1G'] = 0.0006 #(day-1) (p145)
    # par_SN['FMIN2G'] = 0.0272 #(%ARGIS-1) para pot rate min a ARGIs (p145)
    # par_SN['FMIN3G'] = 0.0167 #(% CALC-1)para pot rate min a CALCss (p145)
    #
    # par_SN['FINERTG'] = 0.65 #0.65 = default fraction of stable pool (p145) should be smaller in grassland & forest
    # par_SN['PROFHUMs'] = 150. # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220
    #
    # par_SN['HMinMg'] = 0.3 #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
    # par_SN['HoptMg'] = 1. #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)
    #
    # par_SN['TRefg'] = 15. #reference temperature
    # par_SN['FTEMHAg'] = 25. #asymptotic value of FTH (seen as a logistic response)
    # par_SN['FTEMHg'] = 0.120 #(K-1)
    # par_SN['FTEMHB'] = 145.
    #
    # par_SN['FNXg'] = 0.5 #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
    # par_SN['PHMinNITg'] = 3. #pH min de nitrification (prop of Field Capacity) #value p149
    # par_SN['PHMaxNITg'] = 5.5 #pH max de nitrification (prop of Field Capacity) #value p149
    # par_SN['HMinNg'] = 0.67 #Humidite min de nitrification #value p149
    # par_SN['HoptNg'] = 1. #Humidite opt de nitrification  #value p149
    # par_SN['TNITMINg'] = 5. #Temperature min de nitrification #  (degreC)value p151
    # par_SN['TNITOPTg'] = 20. #Temperature opt de nitrification #  (degreC)value p151
    # par_SN['TNITMAXg'] = 45. #Temperature max de nitrification #  (degreC)value p151
    # par_SN['RATIONITs'] = 0. #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
    # par_SN['DIFNg'] = 0.018  #N diffusion coefficient at field capacity (cm2.day-1, p 161)




    ## soil initialisation
    S = solN.SoilN(par_sol, par_SN, soil_number = vsoilnumbers, dxyz = [[Lsol/size[0]]*size[0], [largsol/size[1]]*size[1], [dz_sol]*ncouches_sol], vDA=vDA, vCN=vCN,vMO=vMO, vARGIs = vARGIs,vNO3=vNO3,vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol, obstarac=None, pattern8=pattern8)
    S.init_asw(HRp_init=HRpinit)

    return S








