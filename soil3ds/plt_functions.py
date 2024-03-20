'''

 3DS-soil-model, a 3D Soil model adapted from STICS : plant functions
 ******************************
 Authors: G. Louarn 
 

 

'''


import numpy as np
from numpy import array, multiply, divide, sum
from copy import deepcopy
from soil3ds.miscel_functions import * #soil3ds miscellaneous soil functions


########## diverse plant fonctions - soil water balance


def Transpi_NC(Et0, ls_epsi, ls_FTSW, paramp = [{"leafAlbedo":0.15, "WaterTreshGs":0.4}]):
    """
    
    """
    ls_transp = []
    for i in range(len(ls_epsi)):

        ##Riou
        #coeffTV = ls_epsi[i]/ (1-leafAlbedo)   
        #potentialVineTranspiration = coeffTV * Et0

        ##adaptation eq 7.1 p127 STICS et 7.8 p131 book
        k_eq = 0.7 #LAI equivalent pour un k=0.7; pause k car les deux k et LAI sont inconnu
        LAIeq = -np.log(1-ls_epsi[i])/k_eq
        potentialTranspiration = Et0 * (1 - np.exp(-(k_eq-0.2)*LAIeq))#complement a potentialSoilEvaporation
        #rq: la bonne variable a recuperer a la place de espi serait le taux de couverture ou le transmis vertical pour lequel le correctif 0.2 n'est pas necessaire (Eq. 7.2)
        # faut introduire crop coefficient KMAXp cf eq. 7.7 -> pour le moment limite a Et0 (reference gazon)

      
        if (ls_FTSW[i] > paramp[0]["WaterTreshGs"]):#previousTSW/TTSW
            # la transpiration de la plante est faite a son potentielle
            Ks = 1.
        else:
            # regulation
            # la quantite d'eau presente dans le sol n'est pas suffisante pour
            # que la transpiration de plante se fasse a son maximum   
            Ks = ls_FTSW[i]/paramp[0]["WaterTreshGs"]

        ls_transp.append(Ks*potentialTranspiration)

    return ls_transp
    # leafAlbedo pas necessaire dans cette version
    # parametres pris dans la premiere plante et suppose identiques pour toutes les plantes

#Transpi_NC(2., [0.4], [0.8], leafAlbedo=0.15, FTSWThreshold=0.4)







########## diverse plant fonctions - root distribution



## gestion des racines
##rq: density above 0.5 cm.cm-3 are not taken into account for absorption of water and N in STICS (p90)
##rq2: valeur de 0.1 cm.cm-3 pour definir profondeur d'enracinement dans L2SU(p138)

def vert_roots(dxyz, lvopt):
    """
    
    """
    #""" pour generer un profil vertical de densite racinaire a partir d'une liste - (homogene en x,y)"""
    
    m_roots = []
    for z in range(len(dxyz[2])):
        v = []
        for x in range(len(dxyz[0])):
            vv = []
            for y in range(len(dxyz[1])):
                vv.append(lvopt[z])

            v.append(vv)
        m_roots.append(v)

    return array(m_roots)

#R1 = vert_roots(dxyz, [0.5]*nstrate[2])
#R2 = vert_roots(dxyz, [0.95, 0.85, 0.75, 0.65, 0.55, 0.45, 0.35, 0.25, 0.15, 0.05])


def root_density(ls_roots, S, unit=1):
    """ Compute a list of root length density distribution per plant (unit: cm.cm-3)

    :param ls_roots:
    :param S: Soil object
    :param unit: root length array unit (1: m; 100:cm))
    :return: ls_root_density
    """
    ls_root_density = []
    cor_ = 100. / unit
    for rt in ls_roots:
        plt_density = (rt * cor_) / (S.m_soil_vol * 10 ** 6)  # m->cm / m3 -> cm3
        ls_root_density.append(plt_density)  # cm/cm-3

    return ls_root_density


def effective_root_lengths(ls_roots, tresh = 0.5):
    """
    
    """
    #""" for multiple root systems in competition
    #H0 : perfect symetry locally (in a voxel) for ressource aquisition
    #treshold of effective density similar to single species
    #fraction of effective density to each species/plant = proportion of total root length density"""
    
    m = deepcopy(ls_roots)
    nb_r = len(ls_roots)
    tot_root_dens = sum_mat(ls_roots) #sum bug
    for rt in range(nb_r):
        #frac_root = divide(ls_roots[rt], tot_root_dens)## rq: gere bien les voxels vides
        m[rt] = np.where(tot_root_dens>tresh, tresh*ls_roots[rt]/tot_root_dens, ls_roots[rt])
        
    return m
    #! fournir seuil tresh en densite ou en longueur absolue selon ce qu'on attend


#def effective_root_length(m_root, tresh = 0.5):#faire avec une liste de root systems ls_root_syst pour competition
#    """ for a single root system """
#    m = deepcopy(m_root)
#    for z in range (len(m_root)):
#        for x in range (len(m_root[z])):
#            for y in range (len(m_root[z][x])):
#                if m_root[z][x][y]>tresh:
#                    m[z][x][y] = tresh
#                else:
#                    m[z][x][y] = m_root[z][x][y]
#    return m

#effective_root_length(R2, 0.5)

def build_ls_roots(RLprofil, S):
    """
    
    """
    #""" version 1 root system :a modifier pour prendre en compte une liste de RLprofil"""
    
    idz = list(RLprofil.keys())
    idz.sort()
    RLprofil_ls = []
    for i in idz: RLprofil_ls.append(RLprofil[i])

    R1 = vert_roots(S.dxyz, RLprofil_ls)*100. /  (S.m_soil_vol*100.*100*100.)#profil densite racinaire en cm.cm-3
    ls_roots = [R1]#[R1, R2]#[R3]#
    return ls_roots

#RLprofil = {0: 0.12450386886407872, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0}#sortie l-system
#RLprofil = [{0: 0.12450386886407872, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0}, {0: 0.145, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0}...]#sortie l-system

def build_ls_roots_mult(RLprofil, S):
    """
    
    """
    #""" version 2: prends en compte une liste de RLprofil"""
    
    nbp = len(RLprofil)
    idz = list(RLprofil[0].keys())
    idz.sort()
    ls_roots = []
    for p in range(nbp):
        RLprofil_ls = []
        for i in idz: RLprofil_ls.append(RLprofil[p][i])
        R = vert_roots(S.dxyz, RLprofil_ls)*100. /  (S.m_soil_vol*100.*100*100.)#profil densite racinaire en cm.cm-3
        ls_roots.append(R)# = [R1]#[R1, R2]#[R3]#

    return ls_roots


def RLprof_t(t, ncouches_sol):
    """
    
    """
    #""" pour generer un profil de densite qui varie au cours du temps / phase de test """
    
    RL_profil_max = [0.64]*31+[0.64, 0.6, 0.5, 0.32, 0.16, 0.08, 0.04, 0.02, 0.01]
    n0 = ncouches_sol-t
    if t<ncouches_sol:
        res  = RL_profil_max[-1-t:-1]+[0]*n0
    else:
        res = RL_profil_max
    return res

#for i in range(50):print RLprof_t(i, ncouches_sol)
    






########## diverse plant fonctions - soil nitrogen balance



##### fonctions uptake plantes (avec objet sol + paramp) #####

def HATSvox(paramp, CONCN):
    """

    """
    # """ Specific absorbtion capacity by roots (micromole N.h-1.cm-1 root) eq. 8.36 p 161 """
    # """ paramp = plant parameters """
    HATS = paramp['Vmax1'] * CONCN / (paramp['Kmax1'] + CONCN)
    return HATS

def LATSvox(paramp, CONCN):
    """

    """
    # """ Specific absorbtion capacity by roots (micromole N.h-1.cm-1 root) eq. 8.36 p 161 """
    # """ paramp = plant parameters """
    LATS = paramp['Vmax2']*CONCN / (paramp['Kmax2'] + CONCN)
    return LATS


def VABSvox(paramp, CONCN):
    """ (unit: micromole N.h-1.cm-1 root)
        
    """
    #""" Specific absorbtion capacity by roots (micromole N.h-1.cm-1 root) eq. 8.36 p 161 """
    #""" paramp = plant parameters """
    HATS = HATSvox(paramp, CONCN)
    LATS = LATSvox(paramp, CONCN)
    return HATS+LATS

    #rq: separer HATS et LATS dans des fonctions avec des options?



def FLUXRACs(SN, paramp, ls_lrac):
    """ (unit: kg N.jour-1 per voxel)
        
    """
    #""" list of potential active uptake rate (kg N. voxel-1, day-1) per root system per voxel - fluxTot=total uptake potential by all roots ; 
    #ls_frac_fluxrac = fraction de demande totale par voxel par systeme racinaire - eq. 8.37 p 161 """
    #ls_lrac = ls_ roots (en m) -> a convertir en cm (*100)

    MMA = 71.428 #142.85 #mole d'N.kg-1 (14 g.mole-1)
    ls_flux_rac = []
    fluxTot = SN.m_1*0.
    for i in range(len(ls_lrac)):
        VABS = VABSvox(paramp[i], SN.ConcN())
        #flux = VABS * ls_lrac[i] * 24. / (MMA * 10 ** 6 * 1000) #24h / mole-µmole / g-kg
        flux = VABS * ls_lrac[i] * 100. * 24. / (MMA * 10 ** 6 )  # m-cm/ 24h / mole-µmole

        ls_flux_rac.append(flux)
        fluxTot = fluxTot + flux

    #calcul des fractions par racine
    ls_frac_fluxrac = []
    for i in range(len(ls_flux_rac)):
        ls_frac_fluxrac.append(ls_flux_rac[i] / (fluxTot + 10**-15) )#10**-15 pour eviter division par zero

    return fluxTot, ls_frac_fluxrac
    #renvoie demande flux total et fraction de demande par voxel par systeme racinaire?
    #passer les densite de racines (i.e. ls_roots) plutot que les longueurs (comme STICS)?
    #seuiller lrac  a une longuer efficace max? LVOPTg -> parametre



def FLUXRACs_old(paramp, SN, ls_lrac):
    """

    """
    # """ list of potential active uptake rate (kg N. voxel-1, day-1) per root system per voxel - fluxTot=total uptake potential by all roots ;
    # ls_frac_fluxrac = fraction de demande totale par voxel par systeme racinaire - eq. 8.37 p 161 """
    # ls_lrac = ls_ roots (en m) -> a convertir en cm (*100)

    MMA = 142.85  # mole d'N.kg-1 (g.mole-1) #MMA erreur
    ls_flux_rac = []
    fluxTot = SN.m_1 * 0.
    for i in range(len(ls_lrac)):
        VABS = VABSvox(paramp[i], SN.ConcN_old()) *100. #facteur 100. erreur
        flux = VABS * ls_lrac[i] * 100. * 24. / (MMA * 10 ** 6)  # calcule car retombe pas sur coeff 33.6 indique p 161
        ls_flux_rac.append(flux)
        fluxTot = fluxTot + flux

    # calcul des fractions par racine
    ls_frac_fluxrac = []
    for i in range(len(ls_flux_rac)):
        ls_frac_fluxrac.append(ls_flux_rac[i] / (fluxTot + 10 ** -15))  # 10**-15 pour eviter division par zero

    return fluxTot, ls_frac_fluxrac
    # renvoie demande flux total et fraction de demande par voxel par systeme racinaire?
    # passer les densite de racines (i.e. ls_roots) plutot que les longueurs?



def Potential_NuptakeTot(SN, parSN, paramp, ls_lrac_, ls_mWaterUptakePlt):
    """
        
    """
    #""" Eq. 8.38 """
    # philosophie: 
    # assimilation de l'azote est active, via equation de Devienne (HATS, LATS)-> FLUXRACs( qui est mal nomme!)
    # transport actif peut etre limitee par le flux d'azote a la racine et mobilite de l'N ds le sol (convectif + diffusif = soilNsupply) -> prend le min des 2
    # in fine, n'est rellement preleve que l'N necessaire a la croissance (actual_N uptake) -> borne
    epsilon = 10e-10

    #seuille pour effective root length
    dens_coeff = 100. / (SN.m_soil_vol[0, 0, 0] * 10 ** 6)
    dens_rac = list(map(np.multiply, ls_lrac_, [dens_coeff]*len(ls_lrac_)))  # cm.cm-3
    dens_rac_seuille = effective_root_lengths(dens_rac, tresh=paramp[0]['treshEffRootsN'])  # cm.cm-3 - suppose meme seuil pour tous!
    ls_lrac = list(map(np.multiply, dens_rac_seuille, [1./dens_coeff]*len(ls_lrac_)))  # m - effective root length

    Ntot = (SN.m_NO3 + SN.m_NH4)*SN.m_obstarac
    supNtot, ls_supNtot =  Nsoil_supply(SN, parSN, ls_lrac, ls_mWaterUptakePlt)
    fluxTot, ls_frac_fluxrac = FLUXRACs(SN, paramp, ls_lrac)
    ActUptakeN = SN.m_1*0.
    idmin = SN.m_1*0 #code des facteur limitant uptake N (0=Ntot available, 2=passive soil supply, 1=Active root uptake, 3=No uptake)

    for z in range(len(SN.dxyz[2])):
        for x in range(len(SN.dxyz[0])):
            for y in range(len(SN.dxyz[1])):
                lslim = [Ntot[z][x][y], fluxTot[z][x][y], supNtot[z][x][y]] #mini des 3 sources d'N
                minN = min(lslim)
                id_min = lslim.index(minN)
                if minN > epsilon:
                    ActUptakeN[z][x][y] = minN
                    idmin[z][x][y] = id_min
                else:#no uptake below epsilon treshold
                    ActUptakeN[z][x][y] = 0.
                    idmin[z][x][y] = 3 #none of the three others

    return ActUptakeN, idmin, ls_frac_fluxrac
    # -> ajuste par voxel selon offre du sol et demande des plantes


def Potential_NuptakeTot_Bis(SN, paramp, ls_lrac_):
    """
        
    """
    #""" Eq. 8.38 """
    # test pour compatibilite romain
    # slmt determine par assimilation de l'azote est active, via equation de Devienne (HATS, LATS)-> FLUXRACs( qui est mal nomme!)

    epsilon = 10e-10

    # seuille pour effective root length - suppose identique pour ttes les plantes!
    dens_coeff = 100. / (SN.m_soil_vol[0, 0, 0] * 10 ** 6)
    dens_rac = list(map(np.multiply, ls_lrac_, [dens_coeff] * len(ls_lrac_)))  # cm.cm-3
    dens_rac_seuille = effective_root_lengths(dens_rac, tresh=paramp[0]['treshEffRootsN'])  # cm.cm-3 - suppose meme seuil pour tous!
    ls_lrac = list(map(np.multiply, dens_rac_seuille, [1. / dens_coeff] * len(ls_lrac_)))  # m - effective root length

    Ntot = (SN.m_NO3 + SN.m_NH4)*SN.m_obstarac
    #supNtot, ls_supNtot =  Nsoil_supply(SN, parSN, ls_lrac, ls_mWaterUptakePlt)
    fluxTot, ls_frac_fluxrac = FLUXRACs(SN, paramp, ls_lrac)
    ActUptakeN = SN.m_1*0.
    idmin = SN.m_1*0 #code des facteur limitant uptake N (0=Ntot available, 1=passive soil supply, 2=Active root uptake, 3=No uptake)

    for z in range(len(SN.dxyz[2])):
        for x in range(len(SN.dxyz[0])):
            for y in range(len(SN.dxyz[1])):
                lslim = [Ntot[z][x][y], fluxTot[z][x][y]] #mini des 2 sources d'N
                minN = min(lslim)
                id_min = lslim.index(minN) # (0=Ntot available, 1=Active root uptake)
                if minN > epsilon:
                    ActUptakeN[z][x][y] = minN
                    idmin[z][x][y] = id_min
                else:#no uptake below epsilon treshold
                    ActUptakeN[z][x][y] = 0.
                    idmin[z][x][y] = 3 #none of the three others

    return ActUptakeN, idmin, ls_frac_fluxrac


def Potential_NuptakeTot_old(SN, parSN, paramp, ls_lrac_, ls_mWaterUptakePlt):
    """

    """
    # """ Eq. 8.38 """
    # philosophie:
    # assimilation de l'azote est active, via equation de Devienne (HATS, LATS)-> FLUXRACs( qui est mal nomme!)
    # transport actif peut etre limitee par le flux d'azote a la racine et mobilite de l'N ds le sol (convectif + diffusif = soilNsupply) -> prend le min des 2
    # in fine, n'est rellement preleve que l'N necessaire a la croissance (actual_N uptake) -> borne
    epsilon = 10e-10

    # seuille pour effective root length - suppose identique pour ttes les plantes!
    dens_coeff = 100. / (SN.m_soil_vol[0, 0, 0] * 10 ** 6)
    dens_rac = list(map(np.multiply, ls_lrac_, [dens_coeff] * len(ls_lrac_)))  # cm.cm-3
    dens_rac_seuille = effective_root_lengths(dens_rac, tresh=paramp[0]['treshEffRootsN'])  # cm.cm-3 - suppose meme seuil pour tous!
    ls_lrac = list(map(np.multiply, dens_rac_seuille, [1. / dens_coeff] * len(ls_lrac_)))  # m - effective root length

    Ntot = (SN.m_NO3 + SN.m_NH4) * SN.m_obstarac
    supNtot, ls_supNtot = Nsoil_supply(SN, parSN, ls_lrac, ls_mWaterUptakePlt)
    fluxTot, ls_frac_fluxrac = FLUXRACs_old(paramp, SN, ls_lrac)
    ActUptakeN = SN.m_1 * 0.
    idmin = SN.m_1 * 0  # code des facteur limitant uptake N (0=Ntot available, 1=passive soil supply, 2=Active root uptake, 3=No uptake)

    for z in range(len(SN.dxyz[2])):
        for x in range(len(SN.dxyz[0])):
            for y in range(len(SN.dxyz[1])):
                lslim = [Ntot[z][x][y], supNtot[z][x][y], fluxTot[z][x][y]]  # mini des 3 sources d'N
                minN = min(lslim)
                id_min = lslim.index(minN)
                if minN > epsilon:
                    ActUptakeN[z][x][y] = minN
                    idmin[z][x][y] = id_min
                else:  # no uptake below epsilon treshold
                    ActUptakeN[z][x][y] = 0.
                    idmin[z][x][y] = 3  # none of the three others

    return ActUptakeN, idmin, ls_frac_fluxrac
    # -> ajuste par voxel selon offre du sol et demande des plantes


def Distrib_Potential_Nuptake_Plt(SN, parSN, paramp_, ls_lrac, ls_mWaterUptakePlt):#ls_lrac, ls_mWaterUptakePlt):
    """
        
    """
    #""" Eq. 8.38 - distibue uptake entre NO3/NH4 et entre les differentes plantes"""
    epsilon = 10e-12
    PotUpNtot, idmin, ls_frac_fluxrac = Potential_NuptakeTot(SN, parSN, paramp_, ls_lrac, ls_mWaterUptakePlt)

    #distribution de l'uptakeN par plante
    ##ls_rac_tot = SN.m_1*0.
    ##for rt in ls_lrac: 
    ##   ls_rac_tot += rt

    ls_Pot_Nuptake_plt = []
    ##for rt in ls_lrac: 
    ##    frac_rac_tot = rt / (ls_rac_tot +epsilon)
    ##    ls_Pot_Nuptake_plt.append( frac_rac_tot * PotUpNtot)    

    for rt in ls_frac_fluxrac:
        ls_Pot_Nuptake_plt.append( rt * PotUpNtot)    


    return PotUpNtot, ls_Pot_Nuptake_plt, idmin
    #rq: distribution entre plante se fait en fonction de densite relative de longueur uniquementif relatifs?
    #tenir compte des parametre de transport act: ls_frac_fluxrac de FLUXRACs?? -> OK =fait
    # !! Manque confrontation a demande totale des plantes (8.39): uptake n'excede pas demande des pantes integree sur tout le profil!! -> a reprendre
    # ls_demandeN a fournir!
    #doner l'option de renvoyer zero si paramp, ls_lrac, ls_mWaterUptakePlt sont a None (pour faire tourner en sol nu facilement)



def Distrib_Potential_Nuptake_Plt_Bis(SN, paramp_, ls_lrac):
    """
        
    """
    #""" Eq. 8.38 - distibue uptake entre NO3/NH4 et entre les differentes plantes"""
    #test pour compatibilite romain

    epsilon = 10e-12
    PotUpNtot, idmin, ls_frac_fluxrac = Potential_NuptakeTot_Bis(SN, paramp_, ls_lrac)

    ls_Pot_Nuptake_plt = []
    for rt in ls_frac_fluxrac:
        ls_Pot_Nuptake_plt.append( rt * PotUpNtot)


    return PotUpNtot, ls_Pot_Nuptake_plt, idmin


def Distrib_Potential_Nuptake_Plt_old(SN, parSN, paramp_, ls_lrac, ls_mWaterUptakePlt):  # ls_lrac, ls_mWaterUptakePlt):
    """

    """
    # """ Eq. 8.38 - distibue uptake entre NO3/NH4 et entre les differentes plantes"""
    epsilon = 10e-12
    PotUpNtot, idmin, ls_frac_fluxrac = Potential_NuptakeTot_old(SN, parSN, paramp_, ls_lrac, ls_mWaterUptakePlt)

    # distribution de l'uptakeN par plante
    ##ls_rac_tot = SN.m_1*0.
    ##for rt in ls_lrac:
    ##   ls_rac_tot += rt

    ls_Pot_Nuptake_plt = []
    ##for rt in ls_lrac:
    ##    frac_rac_tot = rt / (ls_rac_tot +epsilon)
    ##    ls_Pot_Nuptake_plt.append( frac_rac_tot * PotUpNtot)

    for rt in ls_frac_fluxrac:
        ls_Pot_Nuptake_plt.append(rt * PotUpNtot)

    return PotUpNtot, ls_Pot_Nuptake_plt, idmin
    # rq: distribution entre plante se fait en fonction de densite relative de longueur uniquementif relatifs?
    # tenir compte des parametre de transport act: ls_frac_fluxrac de FLUXRACs?? -> OK =fait
    # !! Manque confrontation a demande totale des plantes (8.39): uptake n'excede pas demande des pantes integree sur tout le profil!! -> a reprendre
    # ls_demandeN a fournir!
    # doner l'option de renvoyer zero si paramp, ls_lrac, ls_mWaterUptakePlt sont a None (pour faire tourner en sol nu facilement)


def Actual_Nuptake_plt(SN, ls_Pot_Nuptake_plt, ls_demandeN):
    """
        
    """
    #""" Eq. 8.39, 8.40 p162 """
    #caclule de demande sur offre par plante! -> ls_DQN
    ls_DQ_N = []
    for i in range(len(ls_demandeN)):
        DQ = min(ls_demandeN[i] / sum(ls_Pot_Nuptake_plt[i]), 1.) #plafonne a ratio a 1: preleve uniquement a hauteur de demande de la plante!
        ls_DQ_N.append(DQ)

    #calcul d'un actual uptake par plante et recalcul total
    ls_Act_Nuptake_plt = []
    ActUpNtot = SN.m_1*0.
    for i in range(len(ls_Pot_Nuptake_plt)):
        Act_Nuptake_i = ls_Pot_Nuptake_plt[i] * ls_DQ_N[i]
        ls_Act_Nuptake_plt.append(Act_Nuptake_i)
        ActUpNtot = ActUpNtot + Act_Nuptake_i

    ### retire les nitrates et ammomium rellement preleves du sol
    ##frac_NO3 =  SN.m_NO3 / (SN.m_NO3 + SN.m_NH4 + 10**-15)
    ##SN.m_NO3 = SN.m_NO3 - frac_NO3*ActUpNtot
    ##SN.m_NH4 = SN.m_NH4 - (1. - frac_NO3)*ActUpNtot
    ###bilan
    ##SN.bilanN['cumUptakePlt'].append(ActUpNtot/SN.soilSurface *10000)

    return ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N


def Actual_Nuptake_plt_Bis(SN, ls_Pot_Nuptake_plt, ls_PltN, paramp= [{"treshmaxN":1.0, "treshminN":0.8}]):
    """
        
    """
    # calcul en fonction de feedback statut plante = frein a l'uptake si N suffisant
    # ls_PltN = ls_NNI si optNNI=True (sinon, concentration en N des racines)

    ls_frein_N = []
    for i in range(len(ls_PltN)):

        # A passer en parametre (valeurs de NNI ou de Npc racine min et max)
        treshmaxN = paramp[0]['treshmaxN']# 1.0 #NNI  et au dessus -> frein=0 / peut aussi etre teneur en N des racines
        treshminN = paramp[0]['treshminN']#0.8 #NNI et en dessous -> frein=1 / peut aussi etre teneur en N des racines

        if ls_PltN[i] > treshmaxN:
            frein = 0.
        elif ls_PltN[i] < treshminN:
            frein = 1.
        else:
            frein = (ls_PltN[i] - treshminN) / (treshmaxN - treshminN)

        ls_frein_N.append(frein)

    #calcul d'un actual uptake par plante et recalcul total
    ls_Act_Nuptake_plt = []
    ActUpNtot = SN.m_1*0.
    for i in range(len(ls_Pot_Nuptake_plt)):
        Act_Nuptake_i = ls_Pot_Nuptake_plt[i] * ls_frein_N[i]
        ls_Act_Nuptake_plt.append(Act_Nuptake_i)
        ActUpNtot = ActUpNtot + Act_Nuptake_i

    #print('frein', ls_frein_N, ls_PltN)
    ### retire les nitrates et ammomium rellement preleves du sol
    ##frac_NO3 =  SN.m_NO3 / (SN.m_NO3 + SN.m_NH4 + 10**-15)
    ##SN.m_NO3 = SN.m_NO3 - frac_NO3*ActUpNtot
    ##SN.m_NH4 = SN.m_NH4 - (1. - frac_NO3)*ActUpNtot
    ###bilan
    ##SN.bilanN['cumUptakePlt'].append(ActUpNtot/SN.soilSurface *10000)

    return ActUpNtot, ls_Act_Nuptake_plt, ls_frein_N
    # parametres pris dans la premiere plante et suppose identiques pour toutes les plantes


def Actual_Nuptake_plt_old(SN, ls_Pot_Nuptake_plt, ls_demandeN):
    """

    """
    # """ Eq. 8.39, 8.40 p162 """
    # caclule de demande sur offre par plante! -> ls_DQN
    ls_DQ_N = []
    for i in range(len(ls_demandeN)):
        DQ = min(ls_demandeN[i] / sum(ls_Pot_Nuptake_plt[i]),
                 1.)  # plafonne a ratio a 1: preleve uniquement a hauteur de demande de la plante!
        ls_DQ_N.append(DQ)

    # calcul d'un actual uptake par plante et recalcul total
    ls_Act_Nuptake_plt = []
    ActUpNtot = SN.m_1 * 0.
    for i in range(len(ls_Pot_Nuptake_plt)):
        Act_Nuptake_i = ls_Pot_Nuptake_plt[i] * ls_DQ_N[i]
        ls_Act_Nuptake_plt.append(Act_Nuptake_i)
        ActUpNtot = ActUpNtot + Act_Nuptake_i

    ### retire les nitrates et ammomium rellement preleves du sol
    ##frac_NO3 =  SN.m_NO3 / (SN.m_NO3 + SN.m_NH4 + 10**-15)
    ##SN.m_NO3 = SN.m_NO3 - frac_NO3*ActUpNtot
    ##SN.m_NH4 = SN.m_NH4 - (1. - frac_NO3)*ActUpNtot
    ###bilan
    ##SN.bilanN['cumUptakePlt'].append(ActUpNtot/SN.soilSurface *10000)

    return ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N






######

def critN (MS, a=4.8, b=-0.33):
    """ courbe critique de dilution de l'N - marche aussi pour array"""
    # MS = array od MS values (T.ha-1)
    vals = a*MS**b #en %
    if vals.size>1:
        for i in range(len(vals)): vals[i]=min(a, vals[i])
    else:
        vals = min(a, vals)
    return vals 


def demandeNdefaut(MSp,dMSp,Npc, surfsolref, a=4.8, b=-0.33):
    """ demande N pour parties aerienne - suppose meme courbe critique pour tout le monde - base sur N crit de la biomasse totale """
    #MSp = array des MSp (g.plant-1)
    #dMSp = array des dMSp (g.plant-1)
    #Npc = array des Npc plante (%)
    #surfsol sur laquelle sont les plantes #m2

    QN = MSp * Npc/100. #gN.plant-1
    MStot = array(sum(MSp+dMSp))/(surfsolref*100.)#MS new (T.ha-1)
    NcritTot = critN (MStot, a, b)#N crit de MS new
    PotN = (MSp + dMSp) * NcritTot/100. #gN.plant-1
    ls_demandeN = PotN-QN
    ls_demandeN[ls_demandeN<0.]=0.#gN.plant-1
    return ls_demandeN


#surfsolref = 0.05
#MSp = array([1.,1.2, 2.])
#dMSp = array([0.1,0.15,0.2])
#Npc = array([4., 3., 2.])
#demandeNdefaut(MSp,dMSp,Npc, surfsolref)

def demandeNdefaut2(MSp,dMSp,Npc, surfsolref, a=4.8, b1=-0.1 ,b2=-0.33):
    """ demande N pour parties aerienne - suppose meme courbe critique pour tout le monde - base sur N crit de la biomasse totale """
    #MSp = array des MSp (g.plant-1)
    #dMSp = array des dMSp (g.plant-1)
    #Npc = array des Npc plante (%)
    #surfsol sur laquelle sont les plantes #m2

    QN = MSp * Npc/100. #gN.plant-1
    MStot = array(sum(MSp+dMSp))/(surfsolref*100.)#MS new (T.ha-1)
    if MStot>=1.:
        NcritTot = a*MStot**b2#critN (MStot, a, b2)#N crit de MS new dense
    else:
        NcritTot = a*MStot**b1#critN (MStot, a, b1)#N crit de MS new isole

    #filtre NcritTot
    NcritTot[NcritTot>9.]=9.#gN.plant-1

    PotN = (MSp + dMSp) * NcritTot/100. #gN.plant-1
    ls_demandeN = PotN-QN
    ls_demandeN[ls_demandeN<0.]=0.#gN.plant-1
    return ls_demandeN, NcritTot, MStot #renvoie aussi Ncrit et MStot


def demandeNroot(MSpiv,dMSpiv,Npcpiv, surfsolref, Noptpiv):
    """ demande N pour parties racinaire - suppose N critique constant - s'applique aux racines et aux pivots """
    #MSp = array des MSp (g.plant-1)
    #dMSp = array des dMSp (g.plant-1)
    #Npc = array des Npc plante (%)
    #surfsol sur laquelle sont les plantes #m2

    QNpiv = MSpiv * Npcpiv/100. #gN.plant-1
    PotNpiv = (MSpiv + dMSpiv) * Noptpiv/100. #gN.plant-1
    ls_demandeN = PotNpiv-QNpiv
    ls_demandeN[ls_demandeN<0.]=0.#gN.plant-1
    return ls_demandeN




