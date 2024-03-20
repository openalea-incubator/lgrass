#### initialisation functions for l-egume

import os
import sys

try:
    import legume

    path_ = os.path.dirname(os.path.abspath(legume.__file__))  # local absolute path of L-egume
except:
    path_ = r'C:\devel\l-egume\legume'  # r'C:\devel\grassland'


sys.path.insert(0, path_)

import numpy as np
#from scipy import *
from copy import deepcopy
import string
import time
import pandas as pd

from soil3ds import soil_moduleN as solN
from riri5 import RIRI5 as riri

import RootDistrib as rtd
import RootMorpho as rt
import ShootMorpho as sh
import IOtable
import IOxls



def init_glob_variables_simVGL(meteo, mng, DOYdeb, path_station, ongletSta):
    """ Initialise global variables used within the L-egume L-system """

    ## station
    station = IOxls.read_plant_param(path_station, ongletSta)  # ou lire ds fichier inis

    ## meteo du jour
    DOY = DOYdeb
    meteo_j = IOxls.extract_dataframe(meteo, ['TmoyDay', 'RG', 'Et0', 'Precip', 'Tmin', 'Tmax', 'Tsol'], 'DOY', val=DOY)
    meteo_j['I0'] = [0.48 * meteo_j['RG'][0] * 10000 / (3600 * 24)]  # flux PAR journalier moyen en W.m-2 / RG en j.cm-2
    mng_j = IOxls.extract_dataframe(mng, ['Coupe', 'Irrig', 'FertNO3', 'FertNH4', 'Hcut'], 'DOY', val=DOY)
    for k in list(meteo_j.keys()): meteo_j[k] = meteo_j[k][0]
    for k in list(mng_j.keys()): mng_j[k] = mng_j[k][0]
    meteo_j['durjour'] = sh.DayLength(station['latitude'], sh.DecliSun(DOY % 365))

    TT = 0
    TTsol = 0
    STEPS_ = meteo_j['TmoyDay'] - 5.  # dTT(meteo_j['TmoyDay'], [ParamP[0]['Tdev']])# #variable remise  ajour chaque jour
    STEPSsol_ = meteo_j['Tsol'] - 5.
    ls_epsi = [0.]  # [0.4, 0.4]##!! correspondance avec les nb de root systems!

    ##coupe
    TT_repousse = 0  # TT de la derniere coupev #resoud pb: TT utilise pour LAI pour NI revenait jamais a zero
    isTTcut = False
    wasTTcut = False
    # TTcutFreq = 18.*32#15.*32 #phyllochron
    isRegrowth = False  # indicateur: est on a la pousse initiale ou bien plus tard?
    Hcut = 1.  # 3.#simple initialisation : est passe en lecture fichier management
    cutNB = 0

    ## divers
    start_time, past_time = time.time(), 0.  # pour recuperer temps de calcul

    # distribution des retard a levee
    # deltalevmoy = 30 #degre.jours
    # deltalevsd = 15
    # pour gerer un decalage a la levee/reprise
    # tir1, tir2 = [], []
    # for i in range(1000):
    #  tir1.append(max(0,random.gauss(deltalevmoy,deltalevsd)))#test_retard.append(random.uniform(0,60))
    #  tir2.append(max(0,random.gauss(deltalevmoy,deltalevsd)))#test_retard.append(random.uniform(0,60))

    # test_retard = [tir1, tir2] #pour gerer deux especes

    return DOY, TT, TTsol, meteo_j, mng_j, STEPS_, STEPSsol_, ls_epsi, TT_repousse, isTTcut, wasTTcut, isRegrowth, cutNB, Hcut, start_time, past_time, station




def init_sol_fromLpy(inis, meteo_j, par_sol, par_SN, discret_solXY, dz_sol, pattern8, opt_residu, obstarac=None):
    """ 3DS soil initialisation from L-py - 3 couches de sol avec propirete differentes maxi"""

    #Tsol lu dans meteo
    Tsol = meteo_j['Tsol']  # 15. #degresC

    # pattern8 en cm!
    Lsol = max((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
    largsol = min((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m

    # vecteurs d'initialisation du sol (pour 3 couches maxi)
    num_nb = list(map(int, inis['num_nb']))  # [6,6,18] #nbr de couche de chaque num de sol
    ncouches_sol = int(inis['ncouches_sol'])#num_nb[0] + num_nb[1] + num_nb[2]
    vsoilnumbers = [1] * num_nb[0] + [2] * num_nb[1] + [3] * num_nb[2]  # convention autorise 3 types d'horizon max
    # vDA = [par_SN['DA'][0]]*num_nb[0] + [par_SN['DA'][1]]*num_nb[1] + [par_SN['DA'][2]]*num_nb[2] #densite apparente de sol
    vCN = [par_SN['CN0_30']] * num_nb[0] + [par_SN['CN30_60']] * num_nb[1] + [par_SN['CN60_90']] * num_nb[2]  # maxi 3 horizons
    vMO = [par_SN['MO0_30']] * num_nb[0] + [par_SN['MO30_60']] * num_nb[1] + [par_SN['MO60_90']] * num_nb[2]  # maxi 3 horizons
    vARGIs = [par_SN['ARGIs0_30']] * num_nb[0] + [par_SN['ARGIs30_60']] * num_nb[1] + [par_SN['ARGIs60_90']] * num_nb[2]
    vCALCs = [par_SN['CALCs']] * ncouches_sol
    vNH4 = inis['NH4']  # [2.]*ncouches_sol # #!! kg d'N.ha-1 (entree de STICS)
    vNO3 = inis['NO3']  # [0.]*ncouches_sol
    HRpinit = inis['HRp']  # []
    if min(HRpinit) < 0:  # code -1 pour pas d'initialisation
        HRpinit = []

    vDA = []
    for i in vsoilnumbers:
        vDA.append(par_sol[str(i)]['DA'])

    # vsoilnumbers = [1]+[2]*3+[3]*13+[4]*13 #numeros de sol du profil -> mesures acsyd11
    # vDA = [1.81]+[1.31]*3+[1.37]*13+[1.42]*13 #densite apparente de sol (mesure pesees initial aschyd11)
    # vCN = [par_SN['CN0_30']]*ncouches_sol #maxi 90cm en strates de 5cm
    # vMO = [par_SN['MO0_30']]*ncouches_sol #maxi 90cm en strates de 5cm
    # vARGIs = [par_SN['ARGIs']]*ncouches_sol #maxi 90cm
    # vCALCs = [par_SN['CALCs']]*ncouches_sol
    # vNH4 = [2.]*ncouches_sol # #!! kg d'N.ha-1 (entree de STICS)
    # coeff = 0.#0.09#coeff perte ressuyage -> a ajuster pour avoir environ 600 kg N.ha-1
    # vNO3 = [91.*coeff]*ncouches_sol # kg d'N.ha-1 (entree de STICS)
    # vNO3 = array([16.96, 16.07, 15.17, 33.92, 33.92, 33.92, 33.92, 62.49, 82.13, 89.27, 76.77, 107.13, 124.98, 142.84, 124.98, 142.84, 160.69, 151.76, 151.76, 142.84, 178.55, 133.91, 98.20, 89.27, 83.92, 89.27, 73.20, 89.27, 87.45, 62.49])*coeff #issu du profil en sol nu
    # HRpinit = [25.5,26.,25.,25.5,26.,26.,26.,26.5,26.5,27.,27.,27.,27.5,27.5,27.5,27.5,27.5,29,29,29,29,29,29,29,29,30,30,30,30,30]#-> mesures ahscyd au jour 195 (140711) -> init sol nu

    ## soil initialisation
    S = solN.SoilN(par_sol, par_SN, soil_number=vsoilnumbers,
                   dxyz=[[Lsol / discret_solXY[0]] * discret_solXY[0], [largsol / discret_solXY[1]] * discret_solXY[1],
                         [dz_sol / 100.] * ncouches_sol], vDA=vDA, vCN=vCN, vMO=vMO, vARGIs=vARGIs, vNO3=vNO3,
                   vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol, obstarac=obstarac, pattern8=pattern8)

    ## initialise humidite si fournit
    if HRpinit != []:  # initialise humidite si un vecteur est fourni
        S.init_asw(HRp_init=HRpinit)

    # Uval = 0.9*2.61#(epaisseur de sol* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)
    # Uval = par_SN['q0']*0.1*sum(S.m_QH20fc[0])*surfsolref / (S.dxyz[2][0]*100.)#(epaisseur de sol (cm)* mm d'eau dans 1cm) #U quantite d'eau dans une couche superieure en mm (5 par default)

    #Uval = par_SN['q0']
    #stateEV = [0., 0., 0.]  # pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
    #HXs = par_sol[str(vsoilnumbers[0])]['teta_fc']  # humidite a la capacite au champ de l'horizon de surface
    #b_ = solN.bEV(par_SN['ACLIMc'], par_SN['ARGIs'], HXs)  # HXs=0.261)#1.#valeur empirique tres proche#0.1#0.63#0.63
    ## print (par_SN['q0'],'Uval', ' b ', Uval, b_, sum(S.m_QH20fc[0]), surfsolref , (S.dxyz[2][0]*100.))

    return S, Tsol#, Uval, stateEV, b_




def init_scene_fromLpy(ParamP, inis, cote, nbcote, station, lsidP, type='damier8'):
    # initialoise la scene L-egume: arrangement des plantes (carto), discretisation souterraine, discretisation aerienne
    # 1) CARTO
    distplantes = cote / nbcote  # 1. #cm
    carto = sh.planter_coordinates(type, cote, nbcote)

    # reduit a une espece si veut simul separee
    if type == 'damier8_sp1' or type == 'damier8_sp2' or type == 'damier16_sp1' or type == 'damier16_sp2' or 'row4_sp1' or 'row4_sp2':
        carto = sh.reduce_carto(carto, lsidP)

    # 2) definition du pattern et discretisation sol
    pattern8 = [[0, 0], [cote, cote]]
    # pattern8 =[[min(xxx)-dp,min(yyy)-dp], [max(xxx)+dp,max(yyy)+dp]]
    # [[-2.5,-2.5], [2.5,2.5]]#[[-5.,-5.],[5.,5.]]#[[-2.5,-2.5], [5,5]]#[[-12.5,-12.5],[12.5,12.5]]
    # pattern8 = [[-22.3/2.,-49.5/2.], [22.3/2.,49.5/2.]]#pattern.8 rhizotron equivalent (cm)

    Lsol = max((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
    largsol = min((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
    surfsolref = Lsol * largsol  # m2
    dz_sol = inis['dz_sol']  # 4.#5. #cm
    ncouches_sol = int(inis['ncouches_sol'])  # 4#10#30
    prof_sol_max = ncouches_sol * dz_sol  # 80.

    discret_solXY = list(map(int, inis['discret_solXY']))  # [10,10]# nb de discretisation du sol en X et en Y
    lims_sol = rtd.lims_soil(pattern8, dxyz=[[Lsol / discret_solXY[0]] * discret_solXY[0], [largsol / discret_solXY[1]] * discret_solXY[1], [dz_sol / 100.] * ncouches_sol])

    # 3) discretisation au niveau aerien
    # creation des grid3D pour calcul de rayonnement
    ls_gammagroup = list(map(int, IOxls.get_lsparami(ParamP, 'gammagroup')))
    setp = list(set(ls_gammagroup))  # set equivalent fonction r.unique!
    n_gamagroup = len(setp)

    for nump in range(len(ParamP)):  # ajout a chaque plante de son id grille dans ParamP
        ParamP[nump]['id_grid'] = setp.index(int(ParamP[nump]['gammagroup']))  # pour savoir id de grille de la plante

    na, dxyz, lims_aer, origin_grid, surf_refVOX = riri.def_na_lims(pattern8, station['dz_aerien'], station['Hmaxcouv'], opt=station['opt1D3D'])  # '1D'
    # na, dxyz, lims_aer, origin_grid, surf_refVOX = riri.def_na_lims(pattern8, dz_aerien, Hmaxcouv,opt='1D') #version 1D, comme avant

    m_lais = np.zeros([n_gamagroup, na[2], na[1], na[0]])  # ngamma, Z,Y,,X
    m_lais_construct = deepcopy(m_lais)  # pour contsruction du m_lais a t+1
    m_laiPlt = np.zeros([len(ParamP), na[2], na[1], na[0]])  # Distrib lai des plantes indivs: nump, Z,Y,,X
    triplets = riri.get_ls_triplets(m_lais[0], opt=station['sky'])  # 'VXpXmYpYm')#opt='V')#

    # liste des k_teta (coeff extinction directionnels) par entite
    ls_dif = []
    # prepa des k_teta par entite
    for i in setp:
        nump = ls_gammagroup.index(i)  # retrouve numero de plante de premiere occurence de gammagroup
        ls_dif.append(ParamP[nump]['k_teta_distf'])

    # print('triplets', len(triplets))
    # print('ls_dif', len(ls_dif), ls_dif)

    res_trans, res_abs_i = [], []
    res_rfr = []

    return carto, distplantes, pattern8, Lsol, largsol, surfsolref, dz_sol, ncouches_sol, prof_sol_max, discret_solXY, lims_sol, ls_gammagroup, setp, n_gamagroup, na, dxyz, lims_aer, origin_grid, surf_refVOX, m_lais, m_lais_construct, m_laiPlt, triplets, ls_dif, res_trans, res_abs_i, res_rfr
    # par logique, passer cote et nbcote dans ini?



def init_plant_residues_fromParamP(S, opt_residu, ParamP, par_SN):
    """ separate initialisation of plant residue from plant files / lystem"""
    # initialise soil residues from a ParamP of plants + update ParamP[nump]['CNRES']

    if opt_residu == 1:  # initialisatio de residus

        # number of groupes?
        ls_groupres = list(map(int, IOxls.get_lsparami(ParamP, 'groupe_resid')))
        setg = list(set(ls_groupres))  # set equivalent fonction r.unique!
        n_groupres = len(setg)
        nbplantes = len(ls_groupres)
        # recup 1 jeu de param pour chaque groupe
        CNRES, CC, WC, Nmires = [], [], [], []
        for i in range(len(setg)):
            for nump in range(nbplantes):
                if ParamP[nump]['groupe_resid'] == setg[i]:
                    ParamP[nump]['CNRES'] = [ParamP[nump]['CNRESlf'], ParamP[nump]['CNRESst'], ParamP[nump]['CNRESr'], ParamP[nump]['CNRESpiv']]
                    CNRES = CNRES + ParamP[nump]['CNRES']  # ajout des 4 classes pardefaut
                    CC = CC + ParamP[nump]['CC']
                    WC = WC + ParamP[nump]['WC']
                    Nmires = Nmires + ParamP[nump]['Nmires']
                    # print ('par',CNRES, CC, WC, Nmires)
                    break  # s'arrete a premiere plante de ce groupe

        if len(setg) == 1:  # si 1 seul grope, met qd meme un deuxieme residu de meme type pour pas planter
            for nump in range(nbplantes):
                if ParamP[nump]['groupe_resid'] == setg[0]:
                    ParamP[nump]['CNRES'] = [ParamP[nump]['CNRESlf'], ParamP[nump]['CNRESst'], ParamP[nump]['CNRESr'],
                                             ParamP[nump]['CNRESpiv']]
                    CNRES = CNRES + ParamP[nump]['CNRES']
                    CC = CC + ParamP[nump]['CC']
                    WC = WC + ParamP[nump]['WC']
                    Nmires = Nmires + ParamP[nump]['Nmires']
                    # print ('par',CNRES, CC, WC, Nmires)
                    break
        # Nmires not used???

        # distrib dans le sol en dur!
        nb_res = len(CNRES)  # 4 types de residus par espece (4 compatiment du papier) * 2 especes #pas utilise jusque la? (force donc cycles boucle pas? ou  ajuster a 1 moment?)
        vAmount = [0.1] * nb_res  # [20.]# T Fresh Weight.ha-1 (equivalent QRES)
        Vprop1 = [1. / 3., 1. / 3., 1. / 3.] + 50 * [0.]  # distribution dans les horizons #-> change 27 en 50 pour etre sur d'avoir le nb d'horizon-> a adapter selon le vrai nbr d'horizons!!!
        vProps = [Vprop1] * nb_res  # [Vprop1]#[Vprop1, Vprop1, Vprop1]

        # S.init_residues(vCNRESt, vAmount, vProps, vWC, vCC)

        print('soil init', CNRES, vAmount, vProps, WC, CC)
        S.init_residues(par_SN, CNRES, vAmount, vProps, WC, CC)

        print('ls_CRES', np.shape(S.ls_CRES))
        print('ls_CBio', np.shape(S.ls_CBio))
        print('parResi', S.parResi)

        return CC



#def init_ParamP(path_plante, ongletP, ongletPvois, nbcote, deltalevmoy, deltalevsd, seed_=0, type='homogeneous', opt=4, ongletScenar1='default', ongletScenar2='default', idscenar1=1, idscenar2=1, mn_sc=None):
def init_ParamP_VGL_old(path_plante, ongletP, ongletPvois, nbcote, deltalevmoy, deltalevsd, Plt_seed, seed_=0, type='homogeneous', opt=4, ongletScenar1='default', ongletScenar2='default', idscenar1=1, idscenar2=1, mn_sc=None, opt_sd=0, opt_covar=0, path_variance_geno=None, path_variance_matrix=None, idscenar1_sd=None, idscenar2_sd=None):
    """ """

    # 1) cree liste des paramtres plante (1 dico par plante)
    # nbcote = nombre de plante sur un cote en supposant repartition homogene
    g4 = IOxls.read_plant_param(path_plante, ongletP)
    g4 = IOxls.modif_param(g4, ongletP, ongletScenar1, idscenar1, mn_sc=mn_sc)
    g5 = IOxls.read_plant_param(path_plante, ongletPvois)
    g5 = IOxls.modif_param(g5, ongletPvois, ongletScenar2, idscenar2, mn_sc=mn_sc)
    if type == 'homogeneous':  # cas d'un couvert monospe homogene
        ParamP = [g4] * nbcote * nbcote
    elif type == 'damier8' or type == 'random8' or type == 'damier8_sp1' or type == 'damier8_sp2':  # damier binaire 64 plantes
        if nbcote == 8:
            ParamP = sh.damier8(g4, g5, opt=opt)
        else:
            # if opt_verbose==1:
            print('Error! :' + type + ' option is for a 64 plant design')
    elif type == 'damier16' or type == 'random16' or type == 'damier16_sp1' or type == 'damier16_sp2':  # damier binaire 256 plantes
        if nbcote == 16:
            ParamP = sh.damier16(g4, g5, opt=opt)
        else:
            # if opt_verbose==1:
            print('Error! :' + type + ' option is for a 256 plant design')
    elif type == 'row4' or 'row4_sp1' or 'row4_sp2':  # 4 rangs - 500pl.m-2
        ParamP, cart_ = sh.row4(g4, g5, nbprow=nbcote, opt=opt)
    else:  # defautl= force nb plante comme nbcote
        ParamP = [g4] * nbcote
        # ParamP = [g4]*10#*2#*30#[g6, g4, g6, g4, g6, g4, g6]#[g6,g4,g4,g4,g4,g4,g4,g4,g4,g4,g4,g4,g4]#[g4]
        # [g1, g2, g3]#[g4, g5, g5, g5, g5, g5, g5]

    # 2) ajout variabilite sd si opt_sd==1 (variabilite intra) ; possible seulement si pas analyse de sensibilite (onglet scenar=default)
    # test pour esp 1, Len avec sd=0.5
    # ls_sdpar = [0.5] #ecart type parametre - a passer via un fichier d'entree comme scenar? autrement (multivarie ou directement dans fichier parametre plante?)
    # ls_parname = ['Len'] #liste a recuperer via un fichier d'entree

    if opt_sd == 2:  # new version: lecture des CV et calcul des SD a partir des parametres moyens

        sd_fichier_g4 = pd.read_excel(path_variance_geno, sheet_name=ongletP)
        ls_parname_g4 = list(sd_fichier_g4.columns)[1:]  # liste les noms de colonne a  partir de la deuxieme
        CVs4 = sd_fichier_g4.loc[sd_fichier_g4["id_scenario"] == idscenar1_sd][ls_parname_g4][0:]
        CVs4 = CVs4.iloc[0]
        # calcul des sigma a partir des cv et moy
        moyP4 = []
        for p in ls_parname_g4:
            moyP4.append(g4[p])

        ls_sdpar_g4 = abs(CVs4 * np.array(moyP4))

        sd_fichier_g5 = pd.read_excel(path_variance_geno, sheet_name=ongletPvois)
        ls_parname_g5 = list(sd_fichier_g5.columns)[1:]  # liste les noms de colonne a  partir de la deuxieme
        CVs5 = sd_fichier_g5.loc[sd_fichier_g5["id_scenario"] == idscenar2_sd][ls_parname_g5][0:]
        CVs5 = CVs5.iloc[0]
        # calcul des sigma a partir des cv et moy
        moyP5 = []
        for p in ls_parname_g5:
            moyP5.append(g5[p])

        ls_sdpar_g5 = abs(CVs5 * np.array(moyP5))

        # print(CVs5, CVs4)
        # print(ls_sdpar_g5, ls_sdpar_g4)
        # print(ls_parname_g5, ls_parname_g4)

        # ParamP = IOxls.modif_ParamP_sd(ParamP, g4, ls_parname= ['Len'], ls_sdpar= [0.5])
        # ParamP = IOxls.modif_ParamP_sd(ParamP, g5, ls_parname= ['Len'], ls_sdpar= [0.5])
        # ParamP = IOxls.modif_ParamP_sd(ParamP, g4, ls_parname=  ls_parname_g4, ls_sdpar=ls_sdpar_g4)
        # ParamP = IOxls.modif_ParamP_sd(ParamP, g5, ls_parname=  ls_parname_g5, ls_sdpar=ls_sdpar_g5)
        # print(df1, df2)

        if opt_covar == 1:
            # lecture de deux matrices de correlation -> bon onlet / scenario (id colonne dedie) dans bon fichier "corr_Matrix" -> a faire!
            # idscenar1_sd et idscenar2_sd: meme id que pour opt_sd!
            covar_fichier_g4 = pd.read_excel(path_variance_matrix, sheet_name=ongletP)
            cor_g4 = np.array(covar_fichier_g4[covar_fichier_g4['id_scenario'] == idscenar1_sd]["correlation"])
            size_mat4 = int(np.sqrt(len(cor_g4)))
            cor_mat4 = cor_g4.reshape((size_mat4, size_mat4))

            covar_fichier_g5 = pd.read_excel(path_variance_matrix, sheet_name=ongletPvois)
            cor_g5 = np.array(covar_fichier_g5[covar_fichier_g5['id_scenario'] == idscenar2_sd]["correlation"])
            size_mat5 = int(np.sqrt(len(cor_g5)))
            cor_mat5 = cor_g5.reshape((size_mat5, size_mat5))

            # print('covar g4', cor_mat4, cor_mat5)
            ParamP, df1 = IOxls.modif_ParamP_sdMulti(ParamP, g4, ls_parname=ls_parname_g4, ls_sdpar=ls_sdpar_g4, corrmatrix=cor_mat4)
            ParamP, df2 = IOxls.modif_ParamP_sdMulti(ParamP, g5, ls_parname=ls_parname_g5, ls_sdpar=ls_sdpar_g5, corrmatrix=cor_mat5)

        else:
            # defaut = no matrix of covariance = independant
            ParamP, df1 = IOxls.modif_ParamP_sdMulti(ParamP, g4, ls_parname=ls_parname_g4, ls_sdpar=ls_sdpar_g4, corrmatrix=None)
            ParamP, df2 = IOxls.modif_ParamP_sdMulti(ParamP, g5, ls_parname=ls_parname_g5, ls_sdpar=ls_sdpar_g5, corrmatrix=None)

    elif opt_sd == 1:  # remet old file = lecture directe des valeurs de SD independament des valeurs moyennes

        sd_fichier_g4 = pd.read_excel(path_variance_geno, sheet_name=ongletP)
        ls_parname_g4 = list(sd_fichier_g4.columns)[1:]  # liste les noms de colonne a  partir de la deuxieme
        ls_sdpar_g4 = sd_fichier_g4.loc[sd_fichier_g4["id_scenario"] == idscenar1_sd][ls_parname_g4][0:]
        ls_sdpar_g4 = ls_sdpar_g4.iloc[0]

        sd_fichier_g5 = pd.read_excel(path_variance_geno, sheet_name=ongletPvois)
        ls_parname_g5 = list(sd_fichier_g5.columns)[1:]  # liste les noms de colonne a  partir de la deuxieme
        ls_sdpar_g5 = sd_fichier_g5.loc[sd_fichier_g5["id_scenario"] == idscenar2_sd][ls_parname_g5][0:]
        ls_sdpar_g5 = ls_sdpar_g5.iloc[0]

        if opt_covar == 1:
            covar_fichier_g4 = pd.read_excel(path_variance_matrix, sheet_name=ongletP)
            cor_g4 = np.array(covar_fichier_g4[covar_fichier_g4['id_scenario'] == idscenar1_sd]["correlation"])
            size_mat4 = int(np.sqrt(len(cor_g4)))
            cor_mat4 = cor_g4.reshape((size_mat4, size_mat4))

            covar_fichier_g5 = pd.read_excel(path_variance_matrix, sheet_name=ongletPvois)
            cor_g5 = np.array(covar_fichier_g5[covar_fichier_g5['id_scenario'] == idscenar2_sd]["correlation"])
            size_mat5 = int(np.sqrt(len(cor_g5)))
            cor_mat5 = cor_g5.reshape((size_mat5, size_mat5))

            ParamP, df1 = IOxls.modif_ParamP_sdMulti(ParamP, g4, ls_parname=ls_parname_g4, ls_sdpar=ls_sdpar_g4, corrmatrix=cor_mat4)
            ParamP, df2 = IOxls.modif_ParamP_sdMulti(ParamP, g5, ls_parname=ls_parname_g5, ls_sdpar=ls_sdpar_g5, corrmatrix=cor_mat5)

        else:
            # defaut = no matrix of covariance = independant
            ParamP, df1 = IOxls.modif_ParamP_sdMulti(ParamP, g4, ls_parname=ls_parname_g4, ls_sdpar=ls_sdpar_g4, corrmatrix=None)
            ParamP, df2 = IOxls.modif_ParamP_sdMulti(ParamP, g5, ls_parname=ls_parname_g5, ls_sdpar=ls_sdpar_g5, corrmatrix=None)

    # 3) ajout de parametre 'recalcule'
    # roots
    for nump in range(len(ParamP)):
        # update des parametre racinaire
        rt.update_root_params(ParamP[nump])  # 'lsDrac', 'nb_ordre_rac', 'lsVrac', 'lsDemanDRac', 'LDs'
        # print 'LDs', ParamP[nump]['LDs2'], ParamP[nump]['LDs3'], ParamP[nump]['lsDrac'], ParamP[nump]['GDs2'], ParamP[nump]['GDs3']

        ParamP[nump]['profilRoot'] = rt.rootTropism(ParamP[nump]['IncRoot0'], ParamP[nump]['g_root'], segment=0.3, Long=300.)
        # ParamP[nump]['profilRoot'] = rt.rootTropism_df(ParamP[nump]['IncRoot0'], ParamP[nump]['g_root'], segment=0.3, Long=300.) #plus lent!

    # shoot profiles
    ParamP = sh.update_shoot_params(ParamP)

    # 4) random number generators per plant and test_retard vector
    # seed_=1
    # creation d'une liste d'objet random pour chaque plante (tirage random par plante et plus pour toutes les plantes)
    ls_seeds = []
    for nump in range(len(ParamP)):
        ls_seeds.append(Plt_seed(nump + seed_))
        # ls_seeds.append(Plt_seed(seed_))

    test_retard = []
    for nump in range(len(ParamP)):
        test_retard.append(max(0, ls_seeds[nump].randgen.normal(deltalevmoy, deltalevsd)))
        # max(0,random.gauss(deltalevmoy,deltalevsd))

    # !! ls_seeds pas passe pour tirage multivarie opt_sd? -> non, tirage multivarie avec Rssed pour ttes les plantes ensemble

    # 5) reduit a une espece si veut simul separee
    lsidP = range(len(ParamP))  # default value = all
    if type == 'damier8_sp1' or type == 'damier16_sp1' or type == 'row4_sp1':
        ParamP, lsidP = sh.reduce_ParamP(ParamP, ongletP)
        test_retard = sh.reduce_carto(test_retard, lsidP)
        ls_seeds = sh.reduce_carto(ls_seeds, lsidP)
    elif type == 'damier8_sp2' or type == 'damier16_sp2' or type == 'row4_sp2':
        ParamP, lsidP = sh.reduce_ParamP(ParamP, ongletPvois)
        test_retard = sh.reduce_carto(test_retard, lsidP)
        ls_seeds = sh.reduce_carto(ls_seeds, lsidP)

    # print('test_retard', test_retard)

    nbplantes = len(ParamP)
    return ParamP, nbplantes, ls_seeds, lsidP, test_retard

#pour liste espece
def init_ParamP_VGL(path_plante, ls_Spe, nbcote, deltalevmoy, deltalevsd, Plt_seed, seed_=0, type='homogeneous', opt=4, opt_scenar=0, ls_idscenar=[1, 1], mn_sc=None, opt_sd=0, opt_covar=0, path_variance_geno=None, path_variance_matrix=None, ls_idscenar_sd=[None, None], opt_shuffle=0):
    """ """
    # nbcote = nombre de plante sur un cote en supposant repartition homogene

    # 1) cree liste des paramtres plante (1 dico par plante)
    ls_g =[] # liste des parametrage d'espece
    for i in range(len(ls_Spe)):
        ongletP = ls_Spe[i]
        g = IOxls.read_plant_param(path_plante, ongletP)
        #g = default_paramp() # to test default_paramp() parameter set -> OK
        if opt_scenar!=0: #0:'default' -> pas de changement
            ongletScenar = ongletP # same name by convention
            idscenar = ls_idscenar[i] # ordre des ls_Spe by convention
            g = IOxls.modif_param(g, ongletP, ongletScenar, idscenar, mn_sc=mn_sc)

        ls_g.append(g)


    ParamP = sh.planter_order_ParamP(ls_g, type, nbcote, opt, opt_shuffle)

    # 2) modif ParamP et ajout variabilite sd si opt_sd==1 (variabilite intra) ; possible seulement si pas analyse de sensibilite (onglet scenar=default)
    # test pour esp 1, Len avec sd=0.5
    # ls_sdpar = [0.5] #ecart type parametre - a passer via un fichier d'entree comme scenar? autrement (multivarie ou directement dans fichier parametre plante?)
    # ls_parname = ['Len'] #liste a recuperer via un fichier d'entree

    if opt_sd == 2:  # new version: lecture des CV et calcul des SD a partir des parametres moyens

        #recup nom des parmetres et sd par sp
        ls_sdpar_g, ls_parname_g = [], []
        for i in range(len(ls_Spe)):
            ongletP = ls_Spe[i] # same name by convention
            sd_fichier_g4 = pd.read_excel(path_variance_geno, sheet_name=ongletP)
            ls_parname_g_i = list(sd_fichier_g4.columns)[1:]  # liste les noms de colonne a  partir de la deuxieme
            idscenar_sd = ls_idscenar_sd[i] # ordre des ls_Spe by convention
            CVs4 = sd_fichier_g4.loc[sd_fichier_g4["id_scenario"] == idscenar_sd][ls_parname_g_i][0:]
            CVs4 = CVs4.iloc[0]
            # calcul des sigma a partir des cv et moy
            moyP4 = []
            for p in ls_parname_g_i:
                g = ls_g[i]
                moyP4.append(g[p])

            ls_sdpar_g_i = abs(CVs4 * np.array(moyP4))

            ls_sdpar_g.append(ls_sdpar_g_i)
            ls_parname_g.append(ls_parname_g_i)


        # print(CVs5, CVs4)
        # print(ls_sdpar_g5, ls_sdpar_g4)
        # print(ls_parname_g5, ls_parname_g4)

        # ParamP = IOxls.modif_ParamP_sd(ParamP, g4, ls_parname= ['Len'], ls_sdpar= [0.5])
        # ParamP = IOxls.modif_ParamP_sd(ParamP, g5, ls_parname= ['Len'], ls_sdpar= [0.5])
        # ParamP = IOxls.modif_ParamP_sd(ParamP, g4, ls_parname=  ls_parname_g4, ls_sdpar=ls_sdpar_g4)
        # ParamP = IOxls.modif_ParamP_sd(ParamP, g5, ls_parname=  ls_parname_g5, ls_sdpar=ls_sdpar_g5)
        # print(df1, df2)

        if opt_covar == 1:
            #recup matrice de covariance
            ls_cor_mat = []
            for i in range(len(ls_Spe)):
                ongletP = ls_Spe[i] # same name by convention
                covar_fichier_g4 = pd.read_excel(path_variance_matrix, sheet_name=ongletP)
                idscenar_sd = ls_idscenar_sd[i]
                cor_g4 = np.array(covar_fichier_g4[covar_fichier_g4['id_scenario'] == idscenar_sd]["correlation"])
                size_mat4 = int(np.sqrt(len(cor_g4)))
                cor_mat_i = cor_g4.reshape((size_mat4, size_mat4))

                ls_cor_mat.append(cor_mat_i)


            ls_df = []
            for i in range(len(ls_Spe)):
                ParamP, dfi = IOxls.modif_ParamP_sdMulti(ParamP, ls_g[i], ls_parname=ls_parname_g[i], ls_sdpar=ls_sdpar_g[i], corrmatrix=ls_cor_mat[i])
                ls_df.append(dfi)

        else:
            # defaut = no matrix of covariance = independant
            ls_df = []
            for i in range(len(ls_Spe)):
                ParamP, dfi = IOxls.modif_ParamP_sdMulti(ParamP, ls_g[i], ls_parname=ls_parname_g[i], ls_sdpar=ls_sdpar_g[i], corrmatrix= None)
                ls_df.append(dfi)

    elif opt_sd == 1:  # remet old file = lecture directe des valeurs de SD independament des valeurs moyennes

        # recup nom des parmetres et sd par sp
        ls_sdpar_g, ls_parname_g = [], []
        for i in range(len(ls_Spe)):
            ongletP = ls_Spe[i]  # same name by convention
            sd_fichier_g4 = pd.read_excel(path_variance_geno, sheet_name=ongletP)
            ls_parname_g_i = list(sd_fichier_g4.columns)[1:]  # liste les noms de colonne a  partir de la deuxieme
            idscenar_sd = ls_idscenar_sd[i] # ordre des ls_Spe by convention
            ls_sdpar_g_i = sd_fichier_g4.loc[sd_fichier_g4["id_scenario"] == idscenar_sd][ls_parname_g_i][0:]
            ls_sdpar_g_i = ls_sdpar_g_i.iloc[0]

            ls_sdpar_g.append(ls_sdpar_g_i)
            ls_parname_g.append(ls_parname_g_i)


        if opt_covar == 1:
            # recup matrice de covariance
            ls_cor_mat = []
            for i in range(len(ls_Spe)):
                ongletP = ls_Spe[i]  # same name by convention

                covar_fichier_g4 = pd.read_excel(path_variance_matrix, sheet_name=ongletP)
                idscenar_sd = ls_idscenar_sd[i]
                cor_g4 = np.array(covar_fichier_g4[covar_fichier_g4['id_scenario'] == idscenar_sd]["correlation"])
                size_mat4 = int(np.sqrt(len(cor_g4)))
                cor_mat_i = cor_g4.reshape((size_mat4, size_mat4))

                ls_cor_mat.append(cor_mat_i)


            ls_df = []
            for i in range(len(ls_Spe)):
                ParamP, dfi = IOxls.modif_ParamP_sdMulti(ParamP, ls_g[i], ls_parname=ls_parname_g[i], ls_sdpar=ls_sdpar_g[i], corrmatrix=ls_cor_mat[i])
                ls_df.append(dfi)

        else:
            # defaut = no matrix of covariance = independant
            ls_df = []
            for i in range(len(ls_Spe)):
                ParamP, dfi = IOxls.modif_ParamP_sdMulti(ParamP, ls_g[i], ls_parname=ls_parname_g[i], ls_sdpar=ls_sdpar_g[i], corrmatrix=None)
                ls_df.append(dfi)


    # 3) ajout de parametre 'recalcule'
    # roots
    for nump in range(len(ParamP)):
        # update des parametre racinaire
        rt.update_root_params(ParamP[nump])  # 'lsDrac', 'nb_ordre_rac', 'lsVrac', 'lsDemanDRac', 'LDs'
        # print 'LDs', ParamP[nump]['LDs2'], ParamP[nump]['LDs3'], ParamP[nump]['lsDrac'], ParamP[nump]['GDs2'], ParamP[nump]['GDs3']

        ParamP[nump]['profilRoot'] = rt.rootTropism(ParamP[nump]['IncRoot0'], ParamP[nump]['g_root'], segment=0.3, Long=300.)
        # ParamP[nump]['profilRoot'] = rt.rootTropism_df(ParamP[nump]['IncRoot0'], ParamP[nump]['g_root'], segment=0.3, Long=300.) #plus lent!

    # shoot profiles
    ParamP = sh.update_shoot_params(ParamP)

    # 4) random number generators per plant and test_retard vector
    # seed_=1
    # creation d'une liste d'objet random pour chaque plante (tirage random par plante et plus pour toutes les plantes)
    ls_seeds = []
    for nump in range(len(ParamP)):
        ls_seeds.append(Plt_seed(nump + seed_))
        # ls_seeds.append(Plt_seed(seed_))

    test_retard = []
    for nump in range(len(ParamP)):
        test_retard.append(max(0, ls_seeds[nump].randgen.normal(deltalevmoy, deltalevsd)))
        # max(0,random.gauss(deltalevmoy,deltalevsd))

    # !! ls_seeds pas passe pour tirage multivarie opt_sd? -> non, tirage multivarie avec Rssed pour ttes les plantes ensemble

    # 5) reduit a une espece si veut simul separee
    lsidP = range(len(ParamP))  # default value = all
    if type == 'damier8_sp1' or type == 'damier16_sp1' or type == 'row4_sp1':
        ParamP, lsidP = sh.reduce_ParamP(ParamP, ls_Spe[0])
        test_retard = sh.reduce_carto(test_retard, lsidP)
        ls_seeds = sh.reduce_carto(ls_seeds, lsidP)
    elif type == 'damier8_sp2' or type == 'damier16_sp2' or type == 'row4_sp2':
        ParamP, lsidP = sh.reduce_ParamP(ParamP, ls_Spe[1])
        test_retard = sh.reduce_carto(test_retard, lsidP)
        ls_seeds = sh.reduce_carto(ls_seeds, lsidP)

    # print('test_retard', test_retard)

    nbplantes = len(ParamP)
    return ParamP, nbplantes, ls_seeds, lsidP, test_retard


def init_variables_plantes(ParamP, nbplantes, na):
    global epsilon
    epsilon = 10e-10

    # 1) INVAR
    # dico des variables interne instantanes par plante  utilises dans differents calculs ou preparant des sorties

    invar = init_invar(ParamP, nbplantes, epsilon)

    # 2) INVAR_SC: variables a autres echelle que plante (un dico par echelle des surface, surface verte, PARa cumules)
    ## 3 echelles: plante = 'plt', tige ramifiee ='sh', Axe = 'ax'
    ## Surf:
    ## SurfVerte:
    ## PARaF:
    ## MaxPiv: dictionnaire par cle d'axe de biomasse cumulee par pivot
    ## DiampivMax: dictionnaire par cle d'axe de diametre max de pivot
    ## AgePiv: dictionnaire par cle d'axe d'age des pivot en TT

    invar_sc = {'plt': {}, 'sh': {}, 'ax': {}}
    # reorganiser invar sur la base de ces 3 echelles??
    invar_sc['ax']['MaxPiv'] = {}
    invar_sc['ax']['DiampivMax'] = {}
    invar_sc['ax']['AgePiv'] = {}
    invar_sc['ax']['DemCRac'] = {}
    invar_sc['ax']['OfrCRac'] = {}
    invar_sc['ax']['QDCRac'] = {}
    invar_sc['ax']['QDCmoyRac'] = {}  # QD moyen integre dans le temps par pivot
    invar_sc['ax']['StressHRac'] = {}
    invar_sc['ax']['PonderStressHRac'] = {}
    invar_sc['ax']['StressHmoyRac'] = {}  # stress hydrique moyen des racines par pivot integre dans le temps
    invar_sc['ax']['NRac'] = {}  # liste de nb d'apex par ordre pour chaque pivot
    invar_sc['ax']['dlRac'] = {}  # delta de longueur des racines par ordre
    invar_sc['ax']['cumlRac'] = {}  # cumul de longueur de racine par ordre

    # pour stocker NI max
    invar_sc['sh']['MaxNI'] = {}  # dev max des axes primaires -> pour gerer vitesse de redemarrage des pivots

    # 3) Autres variables globales utilisees ds calcul
    ## lsAxes: liste d'apex I(axes) actifs (utilise dans calcLeafStemRatio)
    ## lsApex: liste d'apex I et II actifs (utilise dans calcNB_NI et cumul_lenIN)
    ## lsApexStop : liste d'apex I et II a l'arret
    ## lsApexAll : liste de tous les apex I et II
    ## lsOrgans: liste des organes (Lf/Stp/In/Pet) sur tous les axes (utilise dans cumul_lenIN, calcOffreC, calcDemandeC) (ancien lsActiveAxes)
    lsAxes = []
    lsApex = []
    lsApexStop = []
    lsApexAll = []
    lsOrgans = [['TT', 'organ', 'nump', 'nsh', 'rank', 'rankp', 'strate', 'surf', 'PARaF', 'statut', 'age', 'ordre', 'l','Long', 'DOY', 'cutNB', 'Larg']]
    savelsOrgans = []
    lsFeuilBilanR = [['nump', 'nsh', 'rank', 'rankp', 'status', 'surf', 'id_grid', 'X', 'Y', 'Z', 'Vox2', 'Vox1', 'Vox0', 'sVox', 'paraF']]

    # for i in range(nbplantes): lsOrgans.append([]) #faire une liste d'organe par plante??

    # 4)#initialisation des ls_systrac (#pour recuperer les enveloppes de racine par plante) et ls_roots_prev
    ls_systrac = {}
    for i in range(nbplantes): ls_systrac[i] = []
    ls_roots_prev = []  # to keep ls_roots in memory

    # 5) initialisation des indices de stress par plante a 1. (devrait passer dans invar)
    ls_ftswStress = {'WaterTreshExpSurf': [], 'WaterTreshDevII': [], 'WaterTreshDevI': [], 'WaterTreshFix': [], 'WaterTreshRUE': []}
    ls_NNIStress = {'NTreshRUE': [], 'NTreshExpSurf': [], 'NTreshDev': [], 'NTreshDevII': []}
    ls_TStress = {'stressTRUE': []}
    for i in range(nbplantes):
        ls_ftswStress['WaterTreshExpSurf'].append(1.)
        ls_ftswStress['WaterTreshDevII'].append(1.)
        ls_ftswStress['WaterTreshDevI'].append(1.)
        ls_ftswStress['WaterTreshFix'].append(1.)
        ls_ftswStress['WaterTreshRUE'].append(1.)
        ls_NNIStress['NTreshRUE'].append(1.)
        ls_NNIStress['NTreshExpSurf'].append(1.)
        ls_NNIStress['NTreshDev'].append(1.)
        ls_NNIStress['NTreshDevII'].append(1.)
        ls_TStress['stressTRUE'].append(1.)

    # 6) Variables de profil plante (surface/eaclairement/N)

    LAIprofil, SurfprofilPlant = {}, []
    for i in range(0, na[2]):  LAIprofil[i] = 0.  # initialise variables globales de profils
    for i in range(nbplantes): SurfprofilPlant.append(deepcopy(LAIprofil))  # liste de LAIprofil par plante

    deltaI_I0 = 0.05  # 0.1 #delta entre classe d'aclairement relatif
    nbI_I0 = int(1. / deltaI_I0)  # nb classes d'eclairement relatif
    I_I0Classes = np.arange(deltaI_I0 / 2., 1. + deltaI_I0 / 2., deltaI_I0)  # eclairememnt relatif moyen par classe
    I_I0profilLfPlant = []  # liste de surface de feuille par classe d'eclairement relatif
    for i in range(nbplantes): I_I0profilLfPlant.append(np.zeros(nbI_I0))
    I_I0profilPetPlant = deepcopy(I_I0profilLfPlant)  # liste de longueur cumulee de petiole par classe d'eclairement relatif
    I_I0profilInPlant = deepcopy(I_I0profilLfPlant)  # liste de longueur cumulee d'entre-noeuds par classe d'eclairement relatif

    NaClasses = ParamP[0]['Na0'] * sh.Na_N0(I_I0Classes)  # N0(INN=1.)*Na_N0(I_I0Classes)
    NlClasses = ParamP[0]['NL0Pet'] * sh.Na_N0(I_I0Classes)
    NlinClasses = ParamP[0]['NL0Sh'] * sh.Na_N0(I_I0Classes)
    # actuellement pas utilise -> serait a retirer proprement

    # variable de profil racine (a passer en dico?)
    res_root = []

    # RLProfil, RprospectProfil, rp0, rpp0 = [],[], {}, [] #liste de root length profil par horizon de sol; liste de profil des rayons de la racine primaire;rp0 et rpp0 sont les profils initiaux d'une plante, utilise pour faciliter l'instanciation de la liste de plantes
    # for i in range(0, ncouches_sol): rp0[i]=0.; rpp0.append(0.)
    # for i in range(nbplantes): RLProfil.append(deepcopy(rp0)); RprospectProfil.append(deepcopy(rpp0))

    return invar, invar_sc, lsAxes, lsApex, lsApexStop, lsApexAll, lsOrgans, savelsOrgans, lsFeuilBilanR, ls_systrac, ls_ftswStress, ls_NNIStress, ls_TStress, LAIprofil, SurfprofilPlant, deltaI_I0, nbI_I0, I_I0Classes, I_I0profilLfPlant, I_I0profilPetPlant, I_I0profilInPlant, NaClasses, NlClasses, NlinClasses, res_root, ls_roots_prev, epsilon


def mef_res_sd(ParamP, ls_Spe, path_variance_geno, test_retard, carto, opt_sd):
    """ mise en forme d'un tableau avec position, retard et parametre des plantes indiv"""

    if opt_sd != 0:
        # ls_parname = ['name']+['Len'] # a recuperer=la bonne liste
        sd_fichier_g4 = pd.read_excel(path_variance_geno, sheet_name=ls_Spe[0])
        ls_parname = ['name'] + list(sd_fichier_g4.columns)[
                                1:]  # liste les noms de colonne a  partir de la deuxieme; suppose la meme pour les 2 sp
        nbp = len(IOxls.get_lsparami(ParamP, 'name'))
        res_sd = {'nump': range(0, nbp)}
        res_sd['retard'] = test_retard[0:nbp]
        for p in ls_parname:
            res_sd[p] = IOxls.get_lsparami(ParamP, p)

        # ajout des coord x,y des plantes
        xcarto, ycarto = [], []
        for i in range(len(carto)):
            xcarto.append(carto[i][0]);
            ycarto.append(carto[i][1])

        res_sd['x'] = xcarto
        res_sd['y'] = ycarto

    else:  # si pas opt_sd
        nbp = len(IOxls.get_lsparami(ParamP, 'name'))
        res_sd = {'nump': range(0, nbp)}
        res_sd['retard'] = test_retard[0:nbp]
        # ajout des coord x,y des plantes
        xcarto, ycarto = [], []
        for i in range(len(carto)):
            xcarto.append(carto[i][0]);
            ycarto.append(carto[i][1])

        res_sd['x'] = xcarto
        res_sd['y'] = ycarto

    return res_sd





def init_outputs(ParamP, nbplantes, ncouches_sol, surfsolref):
    """

    :param ParamP:
    :param nbplantes:
    :param ncouches_sol:
    :return:
    """

    # 1) OUTVAR
    outvar = init_outvar(ParamP, nbplantes, surfsolref)

    # 2) Variables dynamique localisee du sol
    # id couches sorties sol (grand rhizotron)
    id_out = list(range(0,
                        ncouches_sol))  # tous les horizons verticaux#[0,1,4,11,17,25]# 5,10,25,60,90,130 cm, id de voxel dans le sol
    nomprof = ['5', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70', '75', '80', '85',
               '90', '95', '100', '105', '110', '115', '120', '125', '130', '135', '140', '145', '150', '155', '160',
               '165', '170', '175', '180', '185', '190', '195', '200']  # a adapter selon dz_sol..

    out_HR = [['var', 'DOY'] + nomprof[0:ncouches_sol]]  # [['DOY','HP5', 'HP10', 'HP25', 'HP60', 'HP90', 'HP130']]

    return outvar, id_out, out_HR


def init_invar(ParamP, nbplantes, epsilon):
    """ Initialise the dictionary 'invar' storing daily state variable of individual plants

    :param ParamP:
    :type ParamP: list
    :param nbplantes:
    :type nbplantes: int
    :param epsilon:
    :type epsilon: float
    :return: invar
    :rtype: dict

    Keys of the 'invar' dictionary (daily state variables):

        * :alive (int): .... - (unit: unitless)
        * :aliveB (int): .... - (unit: unitless)
        * :alivePiv (int): .... - (unit: unitless)
        * :ChangeRoot (float): .... - (unit: ...)
        * :ChangeDepth (float): .... - (unit: ...)
        * :ConcNmoy (float): Average mineral concentration in root voxels - (unit: µmole.L-1)
        * :countSh (int): .... - (unit: ...)
        * :countShExp (int): .... - (unit: ...)
        * :countGelD (int): .... - (unit: ...)
        * :CreservPiv (float): .... - (unit: ...)
        * :cumtranspi (float): .... - (unit: ...)
        * :DemandN_Feuil (float): .... - (unit: ...)
        * :DemandN_Pet (float): .... - (unit: ...)
        * :DemandN_Stem (float): .... - (unit: ...)
        * :DemandN_Tot (float): .... - (unit: ...)
        * :DemandN_TotAer (float): .... - (unit: ...)
        * :DemCp (float): .... - (unit: ...)
        * :DemCp_lf (float): .... - (unit: ...)
        * :DemCp_in (float): .... - (unit: ...)
        * :DemCp_pt (float): .... - (unit: ...)
        * :DiampivMax (float): .... - (unit: ...)
        * :Dplante (float): .... - (unit: ...)
        * :dMSgraine (float): .... - (unit: ...)
        * :dMSenRoot (float): .... - (unit: ...)
        * :dMSenFeuil (float): .... - (unit: ...)
        * :dMSenTige (float): .... - (unit: ...)
        * :dMSenNonRec (float): .... - (unit: ...)
        * :dMSenPiv (float): .... - (unit: ...)
        * :dMSmortGel_aer (float): .... - (unit: ...)
        * :dMSmortPlant_aer (float): .... - (unit: ...)
        * :dMSmortPlant_pivot (float): .... - (unit: ...)
        * :dMSmortPlant_racfine (float): .... - (unit: ...)
        * :dNgraine (float): .... - (unit: ...)
        * :dNmortGel_aer (float): .... - (unit: ...)
        * :dNmortPlant_aer (float): .... - (unit: ...)
        * :dNmortPlant_pivot (float): .... - (unit: ...)
        * :dNmortPlant_racfine (float): .... - (unit: ...)
        * :dRLen2 (float): .... - (unit: ...)
        * :dRLen3 (float): .... - (unit: ...)
        * :dRLenSentot (float): .... - (unit: ...)
        * :dTT (float): .... - (unit: ...)
        * :dTTsol (float): .... - (unit: ...)
        * :dTTphyllo (float): .... - (unit: ...)
        * :firstleaf (float): .... - (unit: ...)
        * :germination (float): .... - (unit: ...)
        * :graineC (float): .... - (unit: ...)
        * :graineN (float): .... - (unit: ...)
        * :Hplante (float): .... - (unit: ...)
        * :isGelDam (int): .... - (unit: ...)
        * :L_Sp (float): .... - (unit: ...)
        * :lsA (float): .... - (unit: ...)
        * :lsAPrev (float): .... - (unit: ...)
        * :lsApexMort (float): .... - (unit: ...)
        * :MS_aerien (float): .... - (unit: ...)
        * :MS_aerienNonRec (float): .... - (unit: ...)
        * :MS_aerienRec (float): .... - (unit: ...)
        * :MS_aer_cumul (float): .... - (unit: ...)
        * :MS_coty (float): .... - (unit: ...)
        * :MS_feuil (float): .... - (unit: ...)
        * :MS_graine (float): .... - (unit: ...)
        * :MS_rac_fine (float): .... - (unit: ...)
        * :MS_rac_fineNet (float): .... - (unit: ...)
        * :MS_pivot (float): .... - (unit: ...)
        * :MS_tige (float): .... - (unit: ...)
        * :MS_senaerien (float): .... - (unit: ...)
        * :MS_tot (float): .... - (unit: ...)
        * :Maerien (float): .... - (unit: ...)
        * :Mfeuil (float): .... - (unit: ...)
        * :Mrac_fine (float): .... - (unit: ...)
        * :Mpivot (float): .... - (unit: ...)
        * :Mtige (float): .... - (unit: ...)
        * :Msenaerien (float): .... - (unit: ...)
        * :Mtot (float): .... - (unit: ...)
        * :NBapexAct (int): .... - (unit: ...)
        * :NBB (int): .... - (unit: ...)
        * :NBBexp (int): .... - (unit: ...)
        * :NBI (int): .... - (unit: ...)
        * :NBD1 (int): .... - (unit: ...)
        * :NBphyto (int): .... - (unit: ...)
        * :NBsh (int): .... - (unit: ...)
        * :Naerien (float): .... - (unit: ...)
        * :NaerienNonRec (float): .... - (unit: ...)
        * :NaerienRec (float): .... - (unit: ...)
        * :Ncoty (float): .... - (unit: ...)
        * :Ngraine (float): .... - (unit: ...)
        * :Npivot (float): .... - (unit: ...)
        * :NreservPiv (float): .... - (unit: ...)
        * :Nrac_fine (float): .... - (unit: ...)
        * :Npc_aer (float): .... - (unit: %)
        * :Npc_aerNonRec (float): .... - (unit: %)
        * :Npc_piv (float): .... - (unit: %)
        * :Npc_rac_fine (float): .... - (unit: %)
        * :Nuptake_sol (float): .... - (unit: ...)
        * :Ndfa (float): .... - (unit: ...)
        * :NNI (float): Individual plant Nitrogen Nutrition Index - (unit: unitless)
        * :PARaPlante (float): .... - (unit: ...)
        * :PARiPlante (float): .... - (unit: ...)
        * :PARaPlanteU (float): .... - (unit: ...)
        * :parap (float): .... - (unit: ...)
        * :parip (float): .... - (unit: ...)
        * :perteN_aerien (float): .... - (unit: ...)
        * :perteN_NonRec (float): .... - (unit: ...)
        * :perteN_rac_fine (float): .... - (unit: ...)
        * :perteN_Piv (float): .... - (unit: ...)
        * :phmgPet (float): .... - (unit: ...)
        * :phmgEntr (float): .... - (unit: ...)
        * :phmgPet_m (float): .... - (unit: ...)
        * :phmgEntr_m (float): .... - (unit: ...)
        * :Qfix (float): .... - (unit: ...)
        * :R_DemandC_Root (float): .... - (unit: ...)
        * :R_DemandC_Shoot (float): .... - (unit: ...)
        * :RDepth (float): .... - (unit: ...)
        * :remob (float): .... - (unit: ...)
        * :rgeq (): .... - (unit: ...)
        * :RLen1 (float): .... - (unit: ...)
        * :RLen2 (float): .... - (unit: ...)
        * :RLen3 (float): .... - (unit: ...)
        * :RLentot (float): .... - (unit: ...)
        * :RLTot (float): .... - (unit: ...)
        * :RLTotNet (float): .... - (unit: m)
        * :RLentotfromRootMass (float): .... - (unit: m)
        * :RLentotfromDev (float): .... - (unit: m)
        * :RUEactu (float): .... - (unit: ...)
        * :RUEpot (float): .... - (unit: ...)
        * :SRL (float): Average plant Specific Root Length - (unit: m.g-1)
        * :SurfPlante (float): .... - (unit: ...)
        * :Surfcoty (float): .... - (unit: ...)
        * :transpi (float): .... - (unit: ...)
        * :TT (float): Thermal time since plant emergence (air temperature)  - (unit: degree.days)
        * :TTsol (float): Thermal time since plant emergence (soil temperature)  - (unit: degree.days)
        * :TTphyllo (float): Phyllochronic time since plant emergence - (unit: phyllochrons)
        * :TTudev (float): ... - (unit: ...)
        * :Udev (float): .... - (unit: ...)
        * :Udevsol (float): .... - (unit: ...)
        * :Udevstress (float): .... - (unit: ...)
    """

    # 1) INVAR
    # dico des variables interne instantanes par plante  utilises dans differents calculs ou preparant des sorties

    ## SurfPlante: liste (par nump) de liste de surface de feuille verte au temps t
    ## PARaPlante: liste (par nump) de liste de PARa feuille verte au temps t
    ## PARiPlante: liste (par nump) de liste de PARi feuille verte+senescente au temps t
    ## parap (anciens dpar / ls_parap): liste de delta de PARa du jour par plante (feuilles vertes) (utilise pour modulation racines notamment)
    ## parip: liste de delta de PARa du jour par plante (feuilles vertes+senescente)
    ## Hplante : liste de hauteur max des plantes
    ## Dplante: liste de diametre max des plantes
    ## Mrac_fine: Liste (par step/jour) de liste de delta de MSracines fines par plante
    ## Mpivot: Liste (par step/jour) de liste de delta de MSpivot par plante
    ## Maerien: Liste (par step/jour) de liste de delta de MSaerien par plante
    ## MS_rac_fine: liste des MSracines_fines cumule au temps t par plante
    ## MS_pivot: liste des MSpivot cumule au temps t par plante
    ## MS_aerien: liste des MSaerien cumule au temps t par plante
    ## MS_aer_cumul: liste des MSaerien cumule, SANS REMISE A ZERO A LA COUPE, pour calcul d'allocation aux racines.
    ## RLTot: liste de total fine root length par plante(m)
    ## Rdepth: liste de profondeur max du pivot par plante (cm)
    ## DiampivMax: liste de diametre max des pivot par plante
    ## countSh: Liste de nb tiges I cumule emis par plante depuis levee (tout A emis de B) - Compteur de Tige I
    ## countShExp: Liste de nb tiges I cumule emis par plante depuis levee par voie B->BDA - Compteur de bourgeons au niveau de la courrone
    ## NBsh : Liste par plante de nb tiges (I ou II) avec nb phytomeres>50% du max
    ## NBI : Liste par plante de nombre de phytomere maxi sur tiges I ayant >75% du max
    ## NBD1 : Liste par plante de nombre de bourgeons dormants D()
    ## NBB : Liste par plante de nombre de bourgeons actifs B()
    ## NBBexp : Liste par plante de nombre de bourgeons actifs B() de statut exp = generateur de nouveaux axes
    ## lsA : Liste par plante de numero de nsh (B() et A()) actifs ou non
    ## lsAPrev : Liste par plante de numero de nsh (B() et A()) actifs ou non au step n-1 : utilise pour maintenir en dormance les D()
    ## lsApexMort : Liste par plante de numero de nsh (A()) mort
    ## DemandN_Feuil : Liste par plante de quantite d'N des feuille pour satisfaire 1 INN = 1 (g.plant-1)
    ## DemandN_Pet : Liste par plante de quantite d'N des petioles pour satisfaire 1 INN = 1 (g.plant-1)
    ## DemandN_Stem : Liste par plante de quantite d'N des tiges pour satisfaire 1 INN = 1 (g.plant-1)
    ## DemandN_Tot : Liste par plante de quantite d'N totales pour satisfaire 1 INN = 1 (g.plant-1)
    ## R_DemandC_Root : Liste par plante de ratio offre/demande C pour croissance des racines
    ## RLen1, Len2,RLen3,RLen4,RLentot : Liste par plante de longueur cumulee de racine d'ordre 1, 2, 3 et total (m)
    ## SRL: Liste par plante de specific root length (m.g-1)
    ## phmgPet: Liste par plante d'effet maximum d'allongement du a la photomorphogenese (petioles)
    ## phmgEntr: Liste par plante d'effet maximum d'allongement du a la photomorphogenese (entrenoeuds)
    ## phmgPet_m: Liste par plante d'effet maximum de reduction de croissance du a la photomorphogenese (petioles)
    ## phmgEntr_m: Liste par plante d'effet maximum de reduction de croissance du a la photomorphogenese (entrenoeuds)
    ## dMSenFeuil: Liste par plante des delta de biomasse de feuille senescent (g.plant-1)
    ## dMSenTige: Liste par plante des delta de biomasse de tige (entre-noeud et petiole) senescent (g.plant-1)
    ## R_DemandC_Shoot: Liste par plante de ratio offre/demande C pour croissance minimu des tiges
    ## NBphyto: Liste par plante de nbr de phytomeres compte sur la nase des entre-noeuds presents
    ## germinaltion: liste of a logical value: 0=no germination, 1=germination, 2=1feuille visible
    ## TTphyllo: liste of phyllochronic time by plant
    # ! RLentot va remplacer RLTot!

    #global epsilon
    #epsilon = 10e-10  # 10e-15

    invar = {'SurfPlante': [], 'PARaPlante': [], 'PARiPlante': [], 'PARaPlanteU': [], 'Hplante': [], 'Dplante': [],
             'RLTot': [], 'RDepth': [], 'parap': [], 'parip': [], 'Mrac_fine': [], 'Mpivot': [], 'Maerien': [],
             'Mfeuil': [], 'MS_coty': [], 'MS_rac_fine': [], 'MS_pivot': [], 'MS_aerien': [], 'MS_feuil': [],
             'MS_aer_cumul': [], 'Mtot': [], 'MS_tot': [], 'DiampivMax': [], 'countSh': [], 'NBsh': [], 'NBI': [],
             'DemandN_Feuil': [], 'DemandN_Pet': [], 'DemandN_Stem': [], 'DemandN_Tot': [], 'DemandN_TotAer': [],
             'NBD1': [], 'NBB': [], 'countShExp': [], 'lsA': [], 'lsAPrev': [], 'lsApexMort': [], 'NBBexp': [],
             'R_DemandC_Root': [], 'RLen1': [], 'RLen2': [], 'RLen3': [], 'RLentot': [], 'SRL': [], 'phmgPet': [],
             'phmgEntr': [], 'phmgPet_m': [], 'phmgEntr_m': [], 'firstleaf': [], 'Naerien': [], 'Npc_aer': [],
             'Npc_piv': [], 'Npc_rac_fine': [], 'Nuptake_sol': [], 'NNI': [], 'Ndfa': [], 'Qfix': [], 'TT': [],
             'TTsol': [], 'dTT': [], 'dTTsol': [], 'dMSenFeuil': [], 'dMSenTige': [], 'R_DemandC_Shoot': [],
             'RUEactu': [], 'DemCp': [], 'DemCp_lf': [], 'DemCp_in': [], 'DemCp_pt': [], 'L_Sp': [], 'remob': [],
             'Nrac_fine': [], 'Npivot': [], 'dRLen2': [], 'dRLen3': [], 'dMSenRoot': [], 'dRLenSentot': [],
             'RLTotNet': [], 'MS_rac_fineNet': [], 'perteN_rac_fine': [], 'Surfcoty': [], 'NBphyto': [],
             'germination': [], 'MS_graine': [], 'Ngraine': [], 'dMSgraine': [], 'dNgraine': [], 'NBapexAct': [],
             'NreservPiv': [], 'rgeq': [], 'transpi': [], 'cumtranspi': [], 'aliveB': [], 'isGelDam': [],
             'dMSmortGel_aer': [], 'dNmortGel_aer': [], 'dTTphyllo': [], 'TTphyllo': [], 'Udev': [], 'Udevstress': [],
             'Udevsol': [], 'TTudev': [], 'RUEpot': [], 'countGelD': [], 'graineC': [], 'graineN': [], 'Ncoty': [],
             'Mtige': [], 'MS_tige': [], 'MS_aerienNonRec': [], 'MS_aerienRec': [], 'NaerienNonRec': [],
             'NaerienRec': [], 'dMSenPiv': [], 'dMSenNonRec': [], 'CreservPiv': [], 'perteN_NonRec': [],
             'perteN_Piv': [], 'perteN_aerien': [], 'Npc_aerNonRec': [], 'Msenaerien': [], 'MS_senaerien': [],
             'dMSmortPlant_aer': [], 'dMSmortPlant_pivot': [], 'dMSmortPlant_racfine': [], 'dNmortPlant_aer': [],
             'dNmortPlant_pivot': [], 'dNmortPlant_racfine': [], 'alivePiv': [], 'alive': [], 'ChangeRoot': [],
             'ChangeDepth': [], 'RLentotfromRootMass': [], 'RLentotfromDev':[], 'ConcNmoy':[]}

    #initialise with empty list for each plant
    for i in range(nbplantes):
        invar['SurfPlante'].append([])
        invar['PARaPlante'].append([])
        invar['PARiPlante'].append([])
        invar['lsA'].append([])
        invar['lsAPrev'].append([])
        invar['lsApexMort'].append([])

    #initialise with null or identical constant value for each plant
    for i in range(nbplantes):
        invar['Hplante'].append(0.)
        invar['Dplante'].append(0.)
        invar['RLTot'].append(0.)
        invar['RDepth'].append(0.)
        invar['parap'].append(0.)
        invar['parip'].append(0.)
        invar['PARaPlanteU'].append(0.)
        invar['DiampivMax'].append(0.1)
        invar['countSh'].append(0)
        invar['NBsh'].append(0.)
        invar['NBI'].append(0.)
        invar['DemandN_Feuil'].append(0.)
        invar['DemandN_Pet'].append(0.)
        invar['DemandN_Stem'].append(0.)
        invar['DemandN_Tot'].append(0.)
        invar['NBD1'].append(0)
        invar['NBB'].append(0)
        invar['countShExp'].append(0)
        invar['NBBexp'].append(0)
        invar['R_DemandC_Root'].append(0)
        invar['RLen1'].append(0)
        invar['RLen2'].append(0)
        invar['RLen3'].append(0)
        invar['RLentot'].append(0)
        invar['SRL'].append(100.)
        invar['phmgPet'].append([1.])
        invar['phmgEntr'].append([1.])
        invar['phmgPet_m'].append([1.])
        invar['phmgEntr_m'].append([1.])
        invar['firstleaf'].append(float("inf"))
        invar['MS_aer_cumul'].append(0.)
        invar['NNI'].append(1.)
        invar['Ndfa'].append(1.)
        invar['TT'].append(0.)
        invar['TTsol'].append(0.)
        invar['dTT'].append(0.)
        invar['dTTsol'].append(0.)
        invar['dMSenFeuil'].append(0.)
        invar['dMSenTige'].append(0.)
        invar['R_DemandC_Shoot'].append(1.)
        invar['DemCp'].append(0.)
        invar['DemCp_lf'].append(0.)
        invar['DemCp_in'].append(0.)
        invar['DemCp_pt'].append(0.)
        invar['MS_pivot'].append(0.)
        invar['remob'].append(0.)
        invar['RLTotNet'].append(0)
        invar['MS_rac_fineNet'].append(0)
        invar['perteN_rac_fine'].append(0)
        invar['NBphyto'].append(0)
        invar['germination'].append(0)
        invar['dMSgraine'].append(0.)
        invar['dNgraine'].append(0.)
        invar['NBapexAct'].append(0)
        invar['NreservPiv'].append(0)
        invar['rgeq'].append(0)
        invar['cumtranspi'].append(0.)
        invar['aliveB'].append(0.)
        invar['isGelDam'].append(0)
        invar['dMSmortGel_aer'].append(0.)
        invar['dNmortGel_aer'].append(0.)
        invar['dTTphyllo'].append(0.)
        invar['TTphyllo'].append(0.)
        invar['Udev'].append(0.)
        invar['Udevstress'].append(0.)
        invar['Udevsol'].append(0.)
        invar['TTudev'].append(0.)
        invar['countGelD'].append(10.)
        invar['graineC'].append(0.)
        invar['graineN'].append(0.)
        invar['MS_aerienRec'].append(0.)
        invar['NaerienRec'].append(0.)
        invar['dMSenPiv'].append(0.)
        invar['dMSenNonRec'].append(0.)
        invar['CreservPiv'].append(0)
        invar['perteN_NonRec'].append(0)
        invar['perteN_Piv'].append(0)
        invar['perteN_aerien'].append(0)
        invar['dMSmortPlant_aer'].append(0)
        invar['dMSmortPlant_pivot'].append(0)
        invar['dMSmortPlant_racfine'].append(0)
        invar['dNmortPlant_aer'].append(0)
        invar['dNmortPlant_pivot'].append(0)
        invar['dNmortPlant_racfine'].append(0)
        invar['alivePiv'].append(0)
        invar['alive'].append(0)
        invar['ChangeRoot'].append(0)
        invar['ChangeDepth'].append(0)
        invar['RLentotfromRootMass'].append(0)
        invar['RLentotfromDev'].append(0)
        invar['ConcNmoy'].append(0)

    # initialise compactments from indidual plant parameters and vectors
    # initialisation des compartiments avec epsilon pour pas que ca bug
    PG = np.array(IOxls.get_lsparami(ParamP, 'PMG')) / 1000.
    invar['MS_graine'] = PG.tolist()
    invar['Mtot'].append(PG.tolist())
    invar['Ngraine'] = np.array(PG) * np.array(IOxls.get_lsparami(ParamP, 'Npc_ini')) / 100.

    frac_coty_ini = np.array(IOxls.get_lsparami(ParamP, 'frac_coty_ini'))  # 0.5 ##a passer en parametre?
    PG = np.array(PG) * epsilon

    invar['Maerien'].append(np.array(PG) * frac_coty_ini)  # 4/5 va aerien
    invar['MS_coty'] = np.array(PG) * frac_coty_ini  # dans les cotyledons
    invar['Mfeuil'].append(np.array(PG) * frac_coty_ini * 0.98)  # tout aerien dans feuil
    invar['Mtige'].append(np.array(PG) * frac_coty_ini * 0.01)  #
    invar['Mrac_fine'].append(np.array(PG) * (1. - frac_coty_ini) * np.array(IOxls.get_lsparami(ParamP, 'frac_rac_fine')))
    invar['Mpivot'].append(np.array(PG) * (1. - frac_coty_ini) * (1. - np.array(IOxls.get_lsparami(ParamP, 'frac_rac_fine'))))
    invar['MS_aerienNonRec'] = np.array(PG) * frac_coty_ini * 0.01
    invar['Msenaerien'].append(np.array(PG) * 0.)

    invar['Naerien'] = np.array(PG) * frac_coty_ini * np.array(IOxls.get_lsparami(ParamP, 'Npc_ini')) / 100.  # meme teneur racine et shoot
    invar['Ncoty'] = np.array(PG) * frac_coty_ini * np.array(IOxls.get_lsparami(ParamP, 'Npc_ini')) / 100.  # meme teneur racine et shoot
    invar['Nrac_fine'] = np.array(PG) * (1. - frac_coty_ini) * np.array(IOxls.get_lsparami(ParamP, 'frac_rac_fine')) * np.array(IOxls.get_lsparami(ParamP, 'Npc_ini')) / 100.  # meme teneur racine et shoot
    invar['Npivot'] = np.array(PG) * (1. - frac_coty_ini) * (1. - np.array(IOxls.get_lsparami(ParamP, 'frac_rac_fine'))) * np.array(IOxls.get_lsparami(ParamP, 'Npc_ini')) / 100.  # meme teneur racine et shoot
    invar['NaerienNonRec'] = np.array(PG) * frac_coty_ini * 0.01 * np.array(IOxls.get_lsparami(ParamP, 'Npc_ini')) / 100.

    return invar


def init_outvar(ParamP, nbplantes, surfsolref):
    """ Initialise the dictionary 'outvar' storing the dynamics of daily state variable of individual plants

    :param ParamP:
    :param nbplantes:
    :return:
    """

    # dico de sorties
    outvar = {'colnames': [], 'pattern': [], 'TT': [], 'TTsol': [], 'SurfPlante': [], 'PARaPlante': [],
              'PARiPlante': [], 'epsi': [], 'Hplante': [], 'Dplante': [], 'dMSaer': [], 'RLTot': [], 'RDepth': [],
              'MS_aerien': [], 'MS_feuil': [], 'MS_tot': [], 'countSh': [], 'demandC': [], 'Leaf_Stem': [], 'NBsh': [],
              'NBI': [], 'time': [], 'FTSW': [], 'Etransp': [], 'DemandN_Feuil': [], 'DemandN_Pet': [],
              'DemandN_Stem': [], 'DemandN_Tot': [], 'Npc': [], 'NBD1': [], 'NBB': [], 'countShExp': [], 'NBBexp': [],
              'R_DemandC_Root': [], 'SRL': [], 'phmgPet': [], 'phmgEntr': [], 'phmgPet_m': [], 'phmgEntr_m': [],
              'Naerien': [], 'Npc_aer': [], 'DemandN_Tot_Aer': [], 'Nuptake_sol': [], 'NNI': [], 'Ndfa': [], 'Qfix': [],
              'dMSenFeuil': [], 'dMSenTige': [], 'MS_pivot': [], 'MS_rac_fine': [], 'R_DemandC_Shoot': [], 'RUE': [],
              'BilanC_PARa': [], 'BilanC_RUE': [], 'BilanCdMStot': [], 'BilanCdMrac_fine': [], 'BilanCdMpivot': [],
              'BilanCdMaer': [], 'BilanCdMSenFeuil': [], 'BilanCdMSenTige': [], 'DemCp': [], 'remob': [], 'Npc_piv': [],
              'Npc_rac_fine': [], 'dRLenSentot': [], 'dMSenRoot': [], 'RLTotNet': [], 'MS_rac_fineNet': [],
              'perteN_rac_fine': [], 'NBphyto': [], 'cutNB': [], 'NBapexAct': [], 'transpi': [], 'cumtranspi': [],
              'aliveB': [], 'dMSmortGel': [], 'dNmortGel': [], 'TTphyllo': [],
              'dTT': [], 'Udevstress': [], 'Udev': [], 'TTudev': [], 'RUEpot': [], 'MS_aerienNonRec': [],
              'MS_aerienRec': [], 'NaerienNonRec': [], 'NaerienRec': [], 'Ncoty': [], 'MS_tige': [], 'graineC': [],
              'graineN': [], 'dMSenPiv': [], 'dMSenNonRec': [], 'BilanCdMSenNonRec': [], 'BilanCdMSenPiv': [],
              'CreservPiv': [], 'NreservPiv': [], 'perteN_NonRec': [], 'perteN_Piv': [], 'perteN_aerien': [],
              'Npc_aerNonRec': [], 'MS_senaerien': [], 'dMSmortPlant_aer': [], 'dMSmortPlant_pivot': [],
              'dMSmortPlant_racfine': [], 'dNmortPlant_aer': [], 'dNmortPlant_pivot': [], 'dNmortPlant_racfine': [],
              'alivePiv': [], 'alive': [], 'ChangeRoot': [], 'RLentotfromRootMass': [], 'RLentotfromDev': [],
              'ConcNmoy': []}

    outvar['pattern'].append(['pattern', 0] + [surfsolref] * nbplantes)
    outvar['colnames'].append(['V1', 'steps'] + IOxls.get_lsparami(ParamP, 'name'))  # ajout des noms d'omglet en 1ere ligne

    return outvar



def default_paramp():
    """ Creates a default parameter dictionnary 'paramp' for defining a set of plant parameters (Fix2)

        Keys of the dictionnary are  parameters for a given individual plant :

        Potential shoot morphogenesis parameters:
            * 'Tdev' :       Parameters for temperature response - Base temperature
            * 'Tmin' :       Tmin beta (Graux)
            * 'Tmax' :       Tmax beta (Graux)
            * 'q' :       Shape parameter beta (Graux)
            * 'phyllochron' :       phyllochron of the primary axis
            * 'phyllochronII' :       phyllochron of the secondary axes
            * 'delai_deb' :       delay of axillary bud budburst on an isolated shoot
            * 'nshoots' :       Maximal number of primary shoots of an isolated plant
            * 'debTallage' :       Stage at with tillering (/primarry branching) occurs in the taproot zone - number of primary leaves
            * 'RvitTallage' :       ratio phyllochone tallage : phyllochrone seminal stem
            * 'delaiMaturBud' :       Delay of bud maturation to reach an equivalent stage of 1 leaf
            * 'nfol' :       Maximal number of leaflets
            * 'Len' :       Maximal Internode length
            * 'Lfeuille' :       Maximal Leaflet length
            * 'Largfeuille' :       Maximal Leaflet Width (negative = not taken into account)
            * 'Lpet' :       Maximal Petiole length
            * 'Lstip' :       Maximal stipule length
            * 'LRS' :       Root segment length
            * 'LenRhiz' :       Rhizome segment length
            * 'ratioM' :       rapport de longueur entre axe I pousse initiale (graine) et autres primaire
            * 'ratioII' :       rapport de longueur entre axe/feuille secondaire et primaire
            * 'HeightTreshAdvRoots' :
            * 'ProbaMaxAdvRoots' :
            * 'delai_AdvRoots' :
            * 'fenetre_AdvRoots' :
            * 'DistLRhizn' :       Binomlial law rhizomes - n parameter
            * 'DistLRhizp' :       Binomlial law rhizomes -p parameter
            * 'aF' :
            * 'aE' :
            * 'aP' :
            * 'aS' :
            * 'delaiF' :
            * 'delaiE' :
            * 'delaiP' :
            * 'delaiS' :
            * 'seuilexpF' :       leaf age 95% of final size
            * 'profilLeafI_Rlens1' :       relative leaf length with respect to node rank - slope 1
            * 'profilLeafI_Rleni1' :       relative leaf length with respect to node rank - ordo 1
            * 'profilLeafI_Rlens2' :       relative leaf length with respect to node rank - slope 2
            * 'profilLeafI_Rleni2' :       relative leaf length with respect to node rank - ordo 2
            * 'profilLeafI_Rlargs1' :       width:length ratio of leaflets with respect to node rank -slope 1
            * 'profilLeafI_Rlargi1' :       width:length ratio of leaflets with respect to node rank - ordo 1
            * 'profilLeafI_Rlargs2' :       width:length ratio of leaflets with respect to node rank -solpe 2
            * 'profilLeafI_Rlargi2' :       width:length ratio of leaflets with respect to node rank - ordo 2
            * 'profilNodeIs1' :       relative node length with respect to node rank - slope 1
            * 'profilNodeIi1' :       relative node length with respect to node rank - ordo 1
            * 'profilNodeIs2' :       relative node length with respect to node rank - slope 2
            * 'profilNodeIi2' :       relative node length with respect to node rank - ordo 2
            * 'profilStpI_ls1' :       relative stipule length with respect to node rank - slope 1
            * 'profilStpI_li1' :       relative stipule length with respect to node rank - ordo 1
            * 'profilStpI_ls2' :       relative stipule length with respect to node rank - slope 2
            * 'profilStpI_li2' :       relative stipule length with respect to node rank - ordo 2
            * 'profilStpI_Rlargs1' :       width/length ratio of stipules with respect to node rank - slope 1
            * 'profilStpI_Rlargi1' :       width/length ratio of stipules with respect to node rank - intercept 1
            * 'profilStpI_Rlargs2' :       width/length ratio of stipules with respect to node rank - slope 2
            * 'profilStpI_Rlargi2' :       width/length ratio of stipules with respect to node rank - intercept 2
            * 'profilPetIs1' :       relative petiole length with respect to node rank - slope 1
            * 'profilPetIi1' :       relative petiole length with respect to node rank - ordo 1
            * 'profilPetIs2' :       relative petiole length with respect to node rank - slope 2
            * 'profilPetIi2' :       relative petiole length with respect to node rank - ordo 2
            * 'profilLeafI_Rnfols' :       relative number of leaflets with respect to node rank - slope
            * 'profilLeafI_Rnfoli' :       relative number of leaflets with respect to node rank - intercept

        Potential root morphogenesis parameters:
            * 'Dmin' :       Minimum root apex diameter
            * 'Dmax' :       Maximum root apex diameter
            * 'DIDm' :       slope of the relationship between mother and daughter root diamter
            * 'IBD' :       Inter-Branch average distance
            * 'ELmax' :       Elongation rate of roots with maximal apex diameter
            * 'DistRA' :       Distance from root apex without branching
            * 'FRD' :       Fine Root Density (RTD chez pages et al)
            * 'GDs' :       Growth duration for root in function of their diameter
            * 'nbnodales' :       number of nodale roots emitted per node
            * 'alloc_rootB' :       Allocation to roots - parameter beta (ratio of total DM production)
            * 'alloc_rootA' :       Allocation to roots - parameter alpha (power of total DM production)
            * 'frac_rac_fine' :       Root DM fraction allocated to fine roots
            * 'frac_remob' :       maximum fraction of taproot biomass in C reserve daily remolisable

        Potential N uptake and allocation parameters:
            * 'PMG' :       Poids de 1000 grains
            * 'DurGraine' :       Duration of seed C and N provision - no stress during this period
            * 'Npc_ini' :       intial %N of germinating seedlings
            * 'frac_coty_ini' :       Fraction of seed mass allocated to the leaf cotyledon
            * 'RUE' :       Radiation use efficiency - whole plant & PAR
            * 'NODcost' :       relative reduction of RUE at 100% fixation
            * 'SLAmin' :       minimum specific leaf area - to compute C demand of new tissues
            * 'SNLmin' :       minimum specific node length - to compute C demand of new tissues
            * 'SPLmin' :       minimum specific petiole length - to compute C demand of new tissues
            * 'SRLmin' :       minimal specific root length
            * 'Frac_piv_sem' :
            * 'Frac_piv_loc' :
            * 'fraction_NonRec' :       fraction of Msaerien non recolte
            * 'ADIL' :       N content of the canopy at 1T.ha-1
            * 'BDILi' :       Exponent of the critical dilution curve for isolated plants
            * 'BDIL' :       Exponent of the critical dilution curve for dense stands
            * 'NoptPiv' :       Optimal N concentration of root organs (taproot)
            * 'NoptFR' :       Optimal N concentration of root organs (fine roots)
            * 'NminPiv' :       Structural N content of the tap root
            * 'Na0' :       Specific leaf nitrogen of top sunny leaves at INN=1
            * 'NL0Sh' :       Lineic stem nitrogen of stems bearing top sunny leaves at INN=1
            * 'NL0Pet' :       Lineic petiole nitrogen of petioles bearing top sunny leaves at INN=1
            * 'Vmax1' :       HATS Vmax
            * 'Kmax1' :       HATS Kmax
            * 'Vmax2' :       LATS Vmax
            * 'Kmax2' :       LATS Kmax
            * 'treshEffRootsN' :       Treshold for maximal total root length density per voxel which is effective for N uptake (cm.cm-3) - STICS:0.5 cm.cm-3
            * 'treshminN' :       lower treshold plant N status for feedback response on root nitrogen uptake (either in NNI unit or %N in roots)
            * 'treshmaxN' :       higher treshold plant N status for feedback response on root nitrogen uptake (either in NNI unit or %N in roots)
            * 'MaxFix' :       Maximal Fixation rate
            * 'DurDevFix' :       Duration after germination to reach Maxial Fixation Capacity

        Parameters for plastic reponses to the environment:
            * 'par_tresh' :       I0 fraction allowing bud outgrowth and shoot organogenesis  (axes I and II)
            * 'par_tresh_til' :       I0 fraction allowing branching and bud outgrowth  (axes I and II)
            * 'photomorphPAR_petini' :       ratio of the unconstrained maximal length for PAR=0
            * 'photomorphPAR_petM' :       maximum ratio of the unconstrained maximal length
            * 'photomorphPAR_pett1' :       PAR for which the maximum ratio is reached
            * 'photomorphPAR_pett2' :       PAR above which no photomorphogenetic effect occurs
            * 'photomorphRFR_pets' :       slope of the relationship between R:FR ratio and ratio of unconstrained maximal length
            * 'photomorphRFR_peti' :       intercept of the relationship between R:FR ratio and ratio of unconstrained maximal length
            * 'photomorphPAR_intini' :       ratio of the unconstrained maximal length for PAR=0
            * 'photomorphPAR_intM' :       maximum ratio of the unconstrained maximal length
            * 'photomorphPAR_intt1' :       PAR for which the maximum ratio is reached
            * 'photomorphPAR_intt2' :       PAR above which no photomorphogenetic effect occurs
            * 'photomorphRFR_ints' :       slope of the relationship between R:FR ratio and ratio of unconstrained maximal length
            * 'photomorphRFR_inti' :       intercept of the relationship between R:FR ratio and ratio of unconstrained maximal length
            * 'MaxSurvOmbr' :       duree max de survie d'un apex secondaire A2 a l'ombre
            * 'WaterTreshExpSurfs' :       slope
            * 'WaterTreshExpSurfd' :       FTSW50
            * 'WaterTreshDevIIs' :       slope
            * 'WaterTreshDevIId' :       FTSW50
            * 'WaterTreshDevIs' :       slope
            * 'WaterTreshDevId' :       FTSW50
            * 'WaterTreshElRootss' :       slope
            * 'WaterTreshElRootsd' :       FTSW50
            * 'WaterTreshAdvRoots' :
            * 'WaterTreshFixs' :       slope
            * 'WaterTreshFixd' :       FTSW50
            * 'WaterTreshRUEs' :       slope
            * 'WaterTreshRUEd' :       FTSW50
            * 'WaterTreshGs' :       FTSW treshold for the onset of transpiration reduction (unit: 0-1 fraction)
            * 'NTreshExpSurfs' :       slope
            * 'NTreshExpSurfd' :       NNI50
            * 'NTreshDevs' :       slope
            * 'NTreshDevd' :       NNI50
            * 'NTreshRUEs' :       slope
            * 'NTreshRUEd' :       NNI50
            * 'NTreshDevIIs' :       slope
            * 'NTreshDevIId' :       NNI50
            * 'limStressTalN' :       limit of NNI stress inhibiting primary shoot production
            * 'limStressTalW' :       limit of FTSW stress inhibiting primary shoot production
            * 'TempTreshRUEb' :       treshold of T response = 0
            * 'TempTreshRUEh' :       treshold of T response = 1
            * 'Tgel' :       frost damage treshold
            * 'PPtreshb' :       treshold of PP response = 0 - base photoperiod
            * 'PPtreshh' :       treshold of PP response = 1 - treshold of response
            * 'leafAlbedo' :       Leaf Albedo

        Parameters for plant senescence:
            * 'spanSen' :       leaf live span until onset of senescence in isolated plants
            * 'spanMrt' :       leaf live span until leaf fall in isolated plants
            * 'LDs' :       Root Life span in function of their diameter
            * 'ombF_Ttresh' :       mini duration treshold inducing onset of leaf senescence for shaded leaves
            * 'ombF_Ltresh' :       mini PAR (??) treshold inducing onset of leaf senescence of shaded leaves
            * 'delai_senperenne' :       delay before onset of turnover for perennial tissues
            * 'TOrate_nonrec' :       turnover rate of perennial aerien non rec
            * 'TOrate_piv' :       turnover rate of perennial raproot

        Geometric parameters:
            * 'phyllotaxy' :       phyllotaxy angle
            * 'leafshape' :       coeff d'allometrie entre (longueur*largeur) et surface d'un foliole
            * 'stipshape' :       coeff d'allometrie entre (longueur*largeur) et surface d'un stipule
            * 'gammaFeuil' :       Average Leaf elevation angle
            * 'gammaFeuilSD' :       standard deviation of leaf  elevation  angle
            * 'IncPet' :       Average Petiole elevation angle
            * 'elv0b' :       parametre ditrib uniforme (borne bass)
            * 'elv0h' :       parametre ditrib uniforme (borne haute)
            * 'elvtresh' :       elv0 treshold sensible to change in initial elevation angle with age ; function of elasticity
            * 'Lmaxeffet' :       longueur de tige ou changement de elv0 atteint son max
            * 'IncRoot0' :       Initial secondary root elevation
            * 'elasticity' :
            * 'g_root' :       root gravitropism
            * 'DPivot2_coeff' :       Slope of the relationship between Mstot taproot (g) and poxer 2 of Dmax taproot (cm2)
            * 'ZPivot_min' :       Soil depth under which  taproot remains at Dmin (just for visualisation)
            * 'offset_diamP' :       distance parameter for plant diameter and bud outgrowth around taproot

        General informations
            * 'name' :       to indicate species name
            * 'type' :       to indicate legume (1) or grass (2)
            * 'ActiveBranch' :       To Activate/desactivate secondary branching of shoots
            * 'gotStip' :       To Activate/desactivate stipules on shoots
            * 'gammagroup' :       group of leaf angle / to regroup entities with similar leaf angle distributions during light interception calculation
            * 'groupe_resid' :       group of residues (default: zero=legume; 1 non-legume)

        Plant residue parameters
            * 'CNRESlf' :       C/N plant residue type1: leaves
            * 'CNRESst' :       C/N plant residue type2: stems
            * 'CNRESr' :       C/N plant residue type3: roots
            * 'CNRESpiv' :       C/N plant residue type4: taproot
            * 'WClf' :       Water Content fraction of fresh plant residue type1: leaves
            * 'WCst' :       Water Content fraction of fresh plant residue type2: stems
            * 'WCr' :       Water Content fraction of fresh plant residue type3: roots
            * 'WCpiv' :       Water Content fraction of fresh plant residue type4: taproot
            * 'CClf' :       Carbon Content fraction of dry plant residue type1: leaves
            * 'CCst' :       Carbon Content fraction of dry  plant residue type2: stems
            * 'CCr' :       Carbon Content fraction of dry  plant residue type3: roots
            * 'CCpiv' :       Carbon Content fraction of dry  plant residue type4: taproot
            * 'Nmireslf' :       Mineral N Content fraction of fresh plant residue type1: leaves
            * 'Nmiresst' :       Mineral N Content fraction of fresh plant residue type2: stems
            * 'Nmiresr' :       Mineral N Content fraction of fresh plant residue type3: roots
            * 'Nmirespiv' :       Mineral N Content fraction of fresh plant residue type4: taproot
            

    :return: Default 'paramp' parameter dictionnary

    .. code-block:: python

        paramp = {
                    'Tdev' : 	5,
                    'Tmin' : 	-7.6,
                    'Tmax' : 	39.4,
                    'q' : 	3.22,
                    'phyllochron' : 	30.58,
                    'phyllochronII' : 	35.59,
                    'delai_deb' : 	65.80816,
                    'nshoots' : 	3,
                    'debTallage' : 	7,
                    'RvitTallage' : 	1,
                    'delaiMaturBud' : 	14,
                    'nfol' : 	3,
                    'Len' : 	3.4875,
                    'Lfeuille' : 	4,
                    'Largfeuille' : 	-1,
                    'Lpet' : 	5.1,
                    'Lstip' : 	0,
                    'LRS' : 	2.50,
                    'LenRhiz' : 	3,
                    'ratioM' : 	1,
                    'ratioII' : 	0.866025404,
                    'HeightTreshAdvRoots' : 	-1,
                    'ProbaMaxAdvRoots' : 	1,
                    'delai_AdvRoots' : 	96,
                    'fenetre_AdvRoots' : 	60,
                    'DistLRhizn' : 	10,
                    'DistLRhizp' : 	0,
                    'aF' : 	0.05,
                    'aE' : 	0.05,
                    'aP' : 	0.05,
                    'aS' : 	0.05,
                    'delaiF' : 	50,
                    'delaiE' : 	130,
                    'delaiP' : 	90,
                    'delaiS' : 	50,
                    'seuilexpF' : 	75.58,
                    'profilLeafI_Rlens1' : 	0.151,
                    'profilLeafI_Rleni1' : 	0.005,
                    'profilLeafI_Rlens2' : 	-0.027,
                    'profilLeafI_Rleni2' : 	1.25,
                    'profilLeafI_Rlargs1' : 	-0.067,
                    'profilLeafI_Rlargi1' : 	1.19,
                    'profilLeafI_Rlargs2' : 	-0.027,
                    'profilLeafI_Rlargi2' : 	0.89,
                    'profilNodeIs1' : 	0.386,
                    'profilNodeIi1' : 	-0.453,
                    'profilNodeIs2' : 	-0.00961,
                    'profilNodeIi2' : 	1.213,
                    'profilStpI_ls1' : 	0.093,
                    'profilStpI_li1' : 	0.0741,
                    'profilStpI_ls2' : 	-0.03,
                    'profilStpI_li2' : 	1.3,
                    'profilStpI_Rlargs1' : 	0.0587,
                    'profilStpI_Rlargi1' : 	0.4408,
                    'profilStpI_Rlargs2' : 	0,
                    'profilStpI_Rlargi2' : 	1,
                    'profilPetIs1' : 	0.376,
                    'profilPetIi1' : 	0.225,
                    'profilPetIs2' : 	-0.0497,
                    'profilPetIi2' : 	1.16,
                    'profilLeafI_Rnfols' : 	0,
                    'profilLeafI_Rnfoli' : 	1,
                    'Dmin' : 	0.019,
                    'Dmax' : 	0.106,
                    'DIDm' : 	0.242,
                    'IBD' : 	0.34,
                    'ELmax' : 	0.109,
                    'DistRA' : 	7.6,
                    'FRD' : 	0.1,
                    'GDs' : 	2000,
                    'nbnodales' : 	-1,
                    'alloc_rootB' : 	0.38,
                    'alloc_rootA' : 	0.87,
                    'frac_rac_fine' : 	0.22,
                    'frac_remob' : 	0.1,
                    'PMG' : 	2.29,
                    'DurGraine' : 	120,
                    'Npc_ini' : 	7,
                    'frac_coty_ini' : 	0.8,
                    'RUE' : 	2.3,
                    'NODcost' : 	0.05,
                    'SLAmin' : 	700,
                    'SNLmin' : 	4.5,
                    'SPLmin' : 	31.5,
                    'SRLmin' : 	250,
                    'Frac_piv_sem' : 	0.1,
                    'Frac_piv_loc' : 	0.7,
                    'fraction_NonRec' : 	0.001,
                    'ADIL' : 	4.8,
                    'BDILi' : 	-0.1,
                    'BDIL' : 	-0.33,
                    'NoptPiv' : 	2,
                    'NoptFR' : 	3 ,
                    'NminPiv' : 	1.4,
                    'Na0' : 	2.13,
                    'NL0Sh' : 	0.025,
                    'NL0Pet' : 	0.003,
                    'Vmax1' : 	0.0018,
                    'Kmax1' : 	50,
                    'Vmax2' : 	0.05,
                    'Kmax2' : 	25000,
                    'treshEffRootsN' : 	1000000,
                    'treshminN' : 	0.8,
                    'treshmaxN' : 	1,
                    'MaxFix' : 	0,
                    'DurDevFix' : 	10,
                    'par_tresh' : 	0.333,
                    'par_tresh_til' : 	0.33,
                    'photomorphPAR_petini' : 	0.74,
                    'photomorphPAR_petM' : 	1.17,
                    'photomorphPAR_pett1' : 	140,
                    'photomorphPAR_pett2' : 	183,
                    'photomorphRFR_pets' : 	-0.68,
                    'photomorphRFR_peti' : 	1.57,
                    'photomorphPAR_intini' : 	0.46,
                    'photomorphPAR_intM' : 	1.44,
                    'photomorphPAR_intt1' : 	138,
                    'photomorphPAR_intt2' : 	183,
                    'photomorphRFR_ints' : 	0,
                    'photomorphRFR_inti' : 	1,
                    'MaxSurvOmbr' : 	1000,
                    'WaterTreshExpSurfs' : 	10,
                    'WaterTreshExpSurfd' : 	0.45,
                    'WaterTreshDevIIs' : 	10,
                    'WaterTreshDevIId' : 	0.4,
                    'WaterTreshDevIs' : 	18,
                    'WaterTreshDevId' : 	0.2,
                    'WaterTreshElRootss' : 	18,
                    'WaterTreshElRootsd' : 	0.2,
                    'WaterTreshAdvRoots' : 	0.3,
                    'WaterTreshFixs' : 	12,
                    'WaterTreshFixd' : 	0.25,
                    'WaterTreshRUEs' : 	12,
                    'WaterTreshRUEd' : 	0.1,
                    'WaterTreshGs' : 	0.4,
                    'NTreshExpSurfs' : 	7,
                    'NTreshExpSurfd' : 0.5,
                    'NTreshDevs' : 	7,
                    'NTreshDevd' : 	0.51,
                    'NTreshRUEs' : 	7,
                    'NTreshRUEd' : 	0.42,
                    'NTreshDevIIs' : 	7,
                    'NTreshDevIId' : 	0.51,
                    'limStressTalN' : 	0,
                    'limStressTalW' : 	0,
                    'TempTreshRUEb' : 	0,
                    'TempTreshRUEh' : 	15,
                    'Tgel' : 	-1,
                    'PPtreshb' : 	5,
                    'PPtreshh' : 	6,
                    'leafAlbedo' : 	0.15,
                    'spanSen' : 	352,
                    'spanMrt' : 	480,
                    'LDs' : 	600000,
                    'ombF_Ttresh' : 	96,
                    'ombF_Ltresh' : 	10,
                    'delai_senperenne' : 	500,
                    'TOrate_nonrec' : 	0.005,
                    'TOrate_piv' : 	0.0005,
                    'phyllotaxy' : 	90,
                    'leafshape' : 	0.71,
                    'stipshape' : 	0.65,
                    'gammaFeuil' : 	-30,
                    'gammaFeuilSD' : 	20,
                    'IncPet' : 	45,
                    'elv0b' : 	0,
                    'elv0h' : 	85,
                    'elvtresh' : 	0,
                    'Lmaxeffet' : 	30,
                    'IncRoot0' : 	70,
                    'elasticity' : 	0.02,
                    'g_root' : 	0.000104,
                    'DPivot2_coeff' : 	0.393,
                    'ZPivot_min' : 	100,
                    'offset_diamP' : 	0.5,
                    'name' : 	'Fix2',
                    'type' : 	1,
                    'ActiveBranch' : 	1,
                    'gotStip' : 	0,
                    'gammagroup' : 	1,
                    'groupe_resid' : 	0,
                    'CNRESlf' : 	16,
                    'CNRESst' : 	16,
                    'CNRESr' : 	16,
                    'CNRESpiv' : 	16,
                    'WClf' : 	0.7,
                    'WCst' : 	0.7,
                    'WCr' : 	0.7,
                    'WCpiv' : 	0.7,
                    'CClf' : 	0.42,
                    'CCst' : 	0.42,
                    'CCr' : 	0.42,
                    'CCpiv' : 	0.42,
                    'Nmireslf' : 	0.00197,
                    'Nmiresst' : 	0.00197,
                    'Nmiresr' : 	0.00197,
                    'Nmirespiv' : 	0.00197

                    }

    """

    paramp = {}

    paramp['Tdev'] = 5 #
    paramp['Tmin'] = -7.6 #
    paramp['Tmax'] = 39.4 #
    paramp['q'] = 3.22 #
    paramp['phyllochron'] = 30.58 #
    paramp['phyllochronII'] = 35.59 #
    paramp['delai_deb'] = 65.80816 #
    paramp['nshoots'] = 3 #
    paramp['debTallage'] = 7 #
    paramp['RvitTallage'] = 1 #
    paramp['delaiMaturBud'] = 14 #
    paramp['nfol'] = 3 #
    paramp['Len'] = 3.4875 #
    paramp['Lfeuille'] = 4 #
    paramp['Largfeuille'] = -1 #
    paramp['Lpet'] = 5.1 #
    paramp['Lstip'] = 0 #
    paramp['LRS'] = 2.50 #
    paramp['LenRhiz'] = 3 #
    paramp['ratioM'] = 1 #
    paramp['ratioII'] = 0.866025404 #
    paramp['HeightTreshAdvRoots'] = -1 #
    paramp['ProbaMaxAdvRoots'] = 1 #
    paramp['delai_AdvRoots'] = 96 #
    paramp['fenetre_AdvRoots'] = 60 #
    paramp['DistLRhizn'] = 10 #
    paramp['DistLRhizp'] = 0 #
    paramp['aF'] = 0.05 #
    paramp['aE'] = 0.05 #
    paramp['aP'] = 0.05 #
    paramp['aS'] = 0.05 #
    paramp['delaiF'] = 50 #
    paramp['delaiE'] = 130 #
    paramp['delaiP'] = 90 #
    paramp['delaiS'] = 50 #
    paramp['seuilexpF'] = 75.58 #
    paramp['profilLeafI_Rlens1'] = 0.151 #
    paramp['profilLeafI_Rleni1'] = 0.005 #
    paramp['profilLeafI_Rlens2'] = -0.027 #
    paramp['profilLeafI_Rleni2'] = 1.25 #
    paramp['profilLeafI_Rlargs1'] = -0.067 #
    paramp['profilLeafI_Rlargi1'] = 1.19 #
    paramp['profilLeafI_Rlargs2'] = -0.027 #
    paramp['profilLeafI_Rlargi2'] = 0.89 #
    paramp['profilNodeIs1'] = 0.386 #
    paramp['profilNodeIi1'] = -0.453 #
    paramp['profilNodeIs2'] = -0.00961 #
    paramp['profilNodeIi2'] = 1.213 #
    paramp['profilStpI_ls1'] = 0.093 #
    paramp['profilStpI_li1'] = 0.0741 #
    paramp['profilStpI_ls2'] = -0.03 #
    paramp['profilStpI_li2'] = 1.3 #
    paramp['profilStpI_Rlargs1'] = 0.0587 #
    paramp['profilStpI_Rlargi1'] = 0.4408 #
    paramp['profilStpI_Rlargs2'] = 0 #
    paramp['profilStpI_Rlargi2'] = 1 #
    paramp['profilPetIs1'] = 0.376 #
    paramp['profilPetIi1'] = 0.225 #
    paramp['profilPetIs2'] = -0.0497 #
    paramp['profilPetIi2'] = 1.16 #
    paramp['profilLeafI_Rnfols'] = 0 #
    paramp['profilLeafI_Rnfoli'] = 1 #

    paramp['Dmin'] = 0.019 #
    paramp['Dmax'] = 0.106 #
    paramp['DIDm'] = 0.242 #
    paramp['IBD'] = 0.34 #
    paramp['ELmax'] = 0.109 #
    paramp['DistRA'] = 7.6 #
    paramp['FRD'] = 0.1 #
    paramp['GDs'] = 2000 #
    paramp['nbnodales'] = -1 #
    paramp['alloc_rootB'] = 0.38 #
    paramp['alloc_rootA'] = 0.87 #
    paramp['frac_rac_fine'] = 0.22 #
    paramp['frac_remob'] = 0.1 #

    paramp['PMG'] = 2.29 #
    paramp['DurGraine'] = 120 #
    paramp['Npc_ini'] = 7 #
    paramp['frac_coty_ini'] = 0.8 #
    paramp['RUE'] = 2.3 #
    paramp['NODcost'] = 0.05 #
    paramp['SLAmin'] = 700 #
    paramp['SNLmin'] = 4.5 #
    paramp['SPLmin'] = 31.5 #
    paramp['SRLmin'] = 250 #
    paramp['Frac_piv_sem'] = 0.1 #
    paramp['Frac_piv_loc'] = 0.7 #
    paramp['fraction_NonRec'] = 0.001 #
    paramp['ADIL'] = 4.8 #
    paramp['BDILi'] = -0.1 #
    paramp['BDIL'] = -0.33 #
    paramp['NoptPiv'] = 2 #
    paramp['NoptFR'] = 3 #
    paramp['NminPiv'] = 1.4 #
    paramp['Na0'] = 2.13 #
    paramp['NL0Sh'] = 0.025 #
    paramp['NL0Pet'] = 0.003 #
    paramp['Vmax1'] = 0.0018 #
    paramp['Kmax1'] = 50 #
    paramp['Vmax2'] = 0.05 #
    paramp['Kmax2'] = 25000 #
    paramp['treshEffRootsN'] = 1000000 #
    paramp['treshminN'] = 0.8 #
    paramp['treshmaxN'] = 1 #
    paramp['MaxFix'] = 0 #
    paramp['DurDevFix'] = 10 #

    paramp['par_tresh'] = 0.333 #
    paramp['par_tresh_til'] = 0.33 #
    paramp['photomorphPAR_petini'] = 0.74 #
    paramp['photomorphPAR_petM'] = 1.17 #
    paramp['photomorphPAR_pett1'] = 140 #
    paramp['photomorphPAR_pett2'] = 183 #
    paramp['photomorphRFR_pets'] = -0.68 #
    paramp['photomorphRFR_peti'] = 1.57 #
    paramp['photomorphPAR_intini'] = 0.46 #
    paramp['photomorphPAR_intM'] = 1.44 #
    paramp['photomorphPAR_intt1'] = 138 #
    paramp['photomorphPAR_intt2'] = 183 #
    paramp['photomorphRFR_ints'] = 0 #
    paramp['photomorphRFR_inti'] = 1 #
    paramp['MaxSurvOmbr'] = 1000 #
    paramp['WaterTreshExpSurfs'] = 10 #
    paramp['WaterTreshExpSurfd'] = 0.45 #
    paramp['WaterTreshDevIIs'] = 10 #
    paramp['WaterTreshDevIId'] = 0.4 #
    paramp['WaterTreshDevIs'] = 18 #
    paramp['WaterTreshDevId'] = 0.2 #
    paramp['WaterTreshElRootss'] = 18 #
    paramp['WaterTreshElRootsd'] = 0.2 #
    paramp['WaterTreshAdvRoots'] = 0.3 #
    paramp['WaterTreshFixs'] = 12 #
    paramp['WaterTreshFixd'] = 0.25 #
    paramp['WaterTreshRUEs'] = 12 #
    paramp['WaterTreshRUEd'] = 0.1 #
    paramp['WaterTreshGs'] = 0.4 #
    paramp['NTreshExpSurfs'] = 7 #
    paramp['NTreshExpSurfd'] = 0.5 #
    paramp['NTreshDevs'] = 7 #
    paramp['NTreshDevd'] = 0.51 #
    paramp['NTreshRUEs'] = 7 #
    paramp['NTreshRUEd'] = 0.42 #
    paramp['NTreshDevIIs'] = 7 #
    paramp['NTreshDevIId'] = 0.51 #
    paramp['limStressTalN'] = 0 #
    paramp['limStressTalW'] = 0 #
    paramp['TempTreshRUEb'] = 0 #
    paramp['TempTreshRUEh'] = 15 #
    paramp['Tgel'] = -1 #
    paramp['PPtreshb'] = 5 #
    paramp['PPtreshh'] = 6 #
    paramp['leafAlbedo'] = 0.15 #

    paramp['spanSen'] = 352 #
    paramp['spanMrt'] = 480 #
    paramp['LDs'] = 600000 #
    paramp['ombF_Ttresh'] = 96 #
    paramp['ombF_Ltresh'] = 10 #
    paramp['delai_senperenne'] = 500 #
    paramp['TOrate_nonrec'] = 0.005 #
    paramp['TOrate_piv'] = 0.0005 #
    paramp['phyllotaxy'] = 90 #
    paramp['leafshape'] = 0.71 #
    paramp['stipshape'] = 0.65 #
    paramp['gammaFeuil'] = -30 #
    paramp['gammaFeuilSD'] = 20 #
    paramp['IncPet'] = 45 #
    paramp['elv0b'] = 0 #
    paramp['elv0h'] = 85 #
    paramp['elvtresh'] = 0 #
    paramp['Lmaxeffet'] = 30 #
    paramp['IncRoot0'] = 70 #
    paramp['elasticity'] = 0.02 #
    paramp['g_root'] = 0.000104 #
    paramp['DPivot2_coeff'] = 0.393 #
    paramp['ZPivot_min'] = 100 #
    paramp['offset_diamP'] = 0.5 #

    paramp['name'] = 'Fix2'  #
    paramp['type'] = 1 #
    paramp['ActiveBranch'] = 1 #
    paramp['gotStip'] = 0 #
    paramp['gammagroup'] = 1 #
    paramp['groupe_resid'] = 0 #

    paramp['CNRESlf'] = 16 #
    paramp['CNRESst'] = 16 #
    paramp['CNRESr'] = 16 #
    paramp['CNRESpiv'] = 16 #
    paramp['WClf'] = 0.7 #
    paramp['WCst'] = 0.7 #
    paramp['WCr'] = 0.7 #
    paramp['WCpiv'] = 0.7 #
    paramp['CClf'] = 0.42 #
    paramp['CCst'] = 0.42 #
    paramp['CCr'] = 0.42 #
    paramp['CCpiv'] = 0.42 #
    paramp['Nmireslf'] = 0.00197 #
    paramp['Nmiresst'] = 0.00197 #
    paramp['Nmiresr'] = 0.00197 #
    paramp['Nmirespiv'] = 0.00197 #

    return paramp

