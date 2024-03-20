'''

 3DS-soil-model, a 3D Soil model adapted from STICS : miscellaneous soil functions
 ******************************
 Authors: G. Louarn 
 

 

'''



import numpy as np
from copy import deepcopy


########## fonctions diverses matrice sol 3D
def sum3(mat):
    """
    
    """
    return mat.sum()#sum(sum(sum(mat)))
    ## pour gerer corectement les sommes totale un peu plus compliquee en muldim...


def sum_mat(ls_m):
    """
    
    """
    #pour gerer correctement les somme de matrices similaires dans une liste (avant sum le faisait bien, mais derniere version bug???)
    s = ls_m[0]
    if len(ls_m)>1:
        for i in range(1,len(ls_m)):
            s=s+ls_m[i]

    return s


def mask(mat, tresh=0.):
    """
    
    """
    #""" retourne matrice masque de 0  et 1
    #utile pour root/no roots; asw/no asw """
    #m = deepcopy(mat) #pas utile->mat cree nouvelle matrice
    return np.where(mat > tresh, 1, 0)
    #m = m*1. #pour convertir en float

#mask(R2)
#mask(asw_t)


def slice_mask(S, id_layer, axis=1):
    """ define a 0-1 mask for a slice of a [z,x,y] np.array

    :param S: Soil object
    :param id_layer: layer id along, x,y or z
    :param axis: slice selection (0: x-y plane; 1: y-z plane; 2:x-z plane)
    :return: [z,x,y] np.array
    """

    mask_ = S.m_1 * 0
    if axis==1: #y-z
        mask_[:, id_layer, :] = S.m_1[:, id_layer, :]
    elif axis==0: #x-y
        mask_[ id_layer, :, :] = S.m_1[ id_layer, :, :]
    elif axis==2: #x-z
        mask_[:, :, id_layer] = S.m_1[:, :, id_layer]

    return mask_

########## diverses fonction parametrage sol

def tetavol_pF_curve(par_sol, pF_):
    """
    
    """
    #""" empirical relationship between volumic water content and pF (=ln(abs(h))where h is soil moisture succion )
    #Driessen, 1986 cited by Penning de Cries et al., 1989
    #WCST = volumic water content at saturation
    #gamma = empirical parameter to account for bulk density and texture"""
    
    teta = float(par_sol['WCST'])*np.exp(-float(par_sol['gamma_theo'])*pF_*pF_)
    return teta
    #tetavol_pF_curve(par_sol, 4.2)


def default_tetaref(par_sol, fc=2., wp=4.2, ad=7.):
    """
    
    """
    par_sol['teta_sat'] = tetavol_pF_curve(par_sol, 0.)
    par_sol['teta_fc'] = tetavol_pF_curve(par_sol, fc)
    par_sol['teta_wp'] = tetavol_pF_curve(par_sol, wp)
    par_sol['teta_ad'] = tetavol_pF_curve(par_sol, ad)
    return par_sol
    #p = default_tetaref(par_sol)


def pF(h):
    """
    
    """
    return np.log10(np.abs(h))




########## diverses fonction - soil water balance



## soil evaporation

def bEV(ACLIMc, ARGIs, HXs):
    """
    
    """
    #""" calcul de  b: coeff empirique sans dimension pour l'evaporation du sol a partir de (STICS manual p 128) : 
    #    ACLIMc : wheather parameter qui depend de la vitesse moyenne du vent (1m-1-> 20; 2m.s-1-> 14)
    #    ARGIs : clay content de la couche superficielle (%)
    #    HXs: volumetric moisture content at fc of the surface layer """
    
    HA = ARGIs/1500.
    return 0.5*ACLIMc*(HXs-HA)*np.power(0.63-HA, 5./3.)

#bEV(ACLIMc=20, ARGIs=20, HXs=0.4)



def soil_EV_1C(Et0, Precip, epsi, previous_state = [0.,0.,0.], leafAlbedo=0.15, U=5., b=0.63):
    """
    
    """
    #""" U: reservoir superficiel = quantite d'eau dans une couche superieure (varie generallement entre 0 et 30 mm
    #    b: coeef empirique sans dimension pour l'evaporation du sol (cf fonction bEV)
    #    meme formulation pour Lebon et al 2003 et STICS"""

    #recup previous state
    previousSES1 = previous_state[0]
    previousSES2 = previous_state[1]
    previousSEP = previous_state[2]

    coeffTV = epsi/ (1.-leafAlbedo)
    potentialSoilEvaporation = (1. - coeffTV) * Et0
    EP = (1. - coeffTV) * Et0
        
    SES1=max(0., previousSES1)
    SES2=max(0.,previousSES2)
    SEP=max(0.,previousSEP)
    
    P = Precip
    P1 = Precip 
    
    if (SES1<U):
        # phase 1
        if (P1>SES1):
            SES1 = 0.
            SES1 = SES1 + EP
        else:
            SES1 = SES1 - P
            SES1 = SES1 + EP
        
        if (SES1>U):
            ER = EP - 0.4*(SES1-U)
            SES2 = 0.6*(SES1-U)
        else:
            ER = EP
    else:
        # phase 2
        if (P1>SES2):
            SEP = 0
            P1 = P1 - SES2
            SES1 = U - P1
            
            if (P1>U):
                SES1 = EP
            else:
                SES1 = SES1 + EP
            
            if (SES1>U):
                ER = EP - 0.4*(SES1-U)
                SES2 = 0.6*(SES1-U)
            else:
                ER = EP
            
        else:
            SES21 = ((2*b*SEP + b**2)**0.5) - b
            SEP = SEP + EP
            SES20 = ((2*b*SEP + b**2)**0.5) - b
            ER = SES20 - SES21
            
            if (P>0):
                ESx = 0.8*P
                
                if (ESx<ER):
                    ESx = ER + P
                
                if (ESx>EP):
                    ESx = EP
                
                ER = ESx
                SES2 = SES2 + ER - P
            else:
                if (ER>EP):
                    ER = EP
                    SES2 = SES2 + ER - P
                else:
                    SES2 = SES2 + ER - P

    state = [SES1, SES2, SEP]
    soilEvaporation = ER 
    return soilEvaporation, state
    #modele Eric #pb effet albedo pour partitionning vegetation

# rq: pourrait gerer par voxel de surface plutot que globalement (ne pas assumer homogeneite horizontale)


def soil_EV_STICS(Et0, Precip, epsi, previous_state=[0., 0., 0.], leafAlbedo=0.15, U=5., b=0.63):
    """
    
    """
    #""" U: reservoir superficiel = quantite d'eau dans une couche superieure (varie generallement entre 0 et 30 mm
    #    b: coeef empirique sans dimension pour l'evaporation du sol (cf fonction bEV)
    #    meme formulation pour Lebon et al 2003 et STICS"""

    # recup previous state
    previousSES1 = previous_state[0]
    previousSES2 = previous_state[1]
    previousSEP = previous_state[2]

    ##adaptation eq 7.1 p127 STICS book
    k_eq = 0.7  # LAI equivalent pour un k=0.7; pause k car les deux k et LAI sont inconnu
    LAIeq = -np.log(1 - epsi) / k_eq
    EP = Et0 * np.exp(-(k_eq - 0.2) * LAIeq)  # potentialSoilEvaporation Eq. 7.1, p127
    # rq: la bonne variable a recuperer a la place de espi serait le taux de couverture ou le transmis vertical pour lequel le correctif 0.2 n'est pas necessaire (Eq. 7.2)

    SES1 = max(0., previousSES1)
    SES2 = max(0., previousSES2)
    SEP = max(0., previousSEP)

    P = Precip
    P1 = Precip

    if (SES1 < U):
        # phase 1
        if (P1 > SES1):
            SES1 = 0.
            SES1 = SES1 + EP
        else:
            SES1 = SES1 - P
            SES1 = SES1 + EP

        if (SES1 > U):
            ER = EP - 0.4 * (SES1 - U)
            SES2 = 0.6 * (SES1 - U)
        else:
            ER = EP
    else:
        # phase 2
        if (P1 > SES2):
            SEP = 0
            P1 = P1 - SES2
            SES1 = U - P1

            if (P1 > U):
                SES1 = EP
            else:
                SES1 = SES1 + EP

            if (SES1 > U):
                ER = EP - 0.4 * (SES1 - U)
                SES2 = 0.6 * (SES1 - U)
            else:
                ER = EP

        else:
            SES21 = ((2 * b * SEP + b ** 2) ** 0.5) - b
            SEP = SEP + EP
            SES20 = ((2 * b * SEP + b ** 2) ** 0.5) - b
            ER = SES20 - SES21

            if (P > 0):
                ESx = 0.8 * P

                if (ESx < ER):
                    ESx = ER + P

                if (ESx > EP):
                    ESx = EP

                ER = ESx
                SES2 = SES2 + ER - P
            else:
                if (ER > EP):
                    ER = EP
                    SES2 = SES2 + ER - P
                else:
                    SES2 = SES2 + ER - P

    state = [SES1, SES2, SEP]
    soilEvaporation = ER
    return soilEvaporation, state
    #combine modele eric et calcul EP avec LAI


# def soil_EV_STICS_old(Et0, Precip, epsi, previous_state=[0., 0., 0.], leafAlbedo=0.15, U=5., b=0.63):
#     """
# 
#     """
#     # """ U: reservoir superficiel = quantite d'eau dans une couche superieure (varie generallement entre 0 et 30 mm
#     #    b: coeef empirique sans dimension pour l'evaporation du sol (cf fonction bEV)
#     #    Et0: Evapotranspiration potenrielle de reference (gazon) en mm
#     #    meme formulation pour Lebon et al 2003 et STICS
#     #    SES1: EP cumule depuis la derniere pluie (previous_state[0])
#     #    SES2: EP cumule depuis passage a phase 2  (previous_state[1])
#     #    SEP : cumul des evap depuis debut phase 2
#     #    ER SoilEvaporation Reel"""
# 
#     ##Riou -> donne des valeur negatives pour forts epsi
#     # coeffTV = epsi/ (1.-leafAlbedo)
#     # EP = (1. - coeffTV) * Et0#potentialSoilEvaporation
# 
#     ##adaptation eq 7.1 p127 STICS book
#     k_eq = 0.7  # LAI equivalent pour un k=0.7; pause k car les deux k et LAI sont inconnu
#     LAIeq = -log(1 - epsi) / k_eq
#     EP = Et0 * exp(-(k_eq - 0.2) * LAIeq)  # potentialSoilEvaporation Eq. 7.1, p127
#     # rq: la bonne variable a recuperer a la place de espi serait le taux de couverture ou le transmis vertical pour lequel le correctif 0.2 n'est pas necessaire (Eq. 7.2)
# 
#     if previous_state[0] + EP < U:  # phase 1
#         SES1 = previous_state[0] + EP
#         SES2 = 0.
#         ER = EP
#         SEP = 0.
#         # mise a jour state!
#     elif previous_state[0] < U and previous_state[0] + EP >= U:  # transition phase1->2
#         fracEP2 = previous_state[0] + EP - U  # precipitation evaporee en phase 2
#         fracEP1 = EP - fracEP2
#         SES1 = previous_state[0] + EP
#         SES2 = fracEP2
#         ER2 = ((2 * b * SES2 + b ** 2) ** 0.5) - b
#         ER = fracEP1 + (fracEP2 / EP) * ER2  # verif
#         SEP = ER2
#     else:  # phase 2
#         SES1 = previous_state[0] + EP
#         SES2 = previous_state[1] + EP
#         ER = ((2 * b * SES2 + b ** 2) ** 0.5) - b - previous_state[
#             2]  # fonction cumul, donc retire valeur de la date precedente
#         SEP = ((2 * b * SES2 + b ** 2) ** 0.5) - b
# 
#     if Precip > 0.:  # comme si pluie vient en fin de journee
#         SES1 = 0.
#         SES2 = 0.
#         SEP = 0.
# 
#     return ER, [SES1, SES2, SEP]  # 0., [SES1, SES2, 0.]#
#     # test par GL
#     # pb de reponse au faiblee pluies: remet tout a zero ; agit fort alors que devrait pas ('effet pompe')


# def soil_EV_STICS_new1(Et0, Precip, epsi, previous_state=[0., 0., 0.], leafAlbedo=0.15, U=5., b=0.63):
#     """
# 
#     """
#     # """ U: reservoir superficiel = quantite d'eau dans une couche superieure (varie generallement entre 0 et 30 mm
#     #    b: coeef empirique sans dimension pour l'evaporation du sol (cf fonction bEV)
#     #    Et0: Evapotranspiration potenrielle de reference (gazon) en mm
#     #    meme formulation pour Lebon et al 2003 et STICS
#     #    SES1: EP cumule depuis la derniere pluie (previous_state[0])
#     #    SES2: EP cumule depuis passage a phase 2  (previous_state[1])
#     #    SEP : cumul des evap depuis debut phase 2
#     #    ER SoilEvaporation Reel"""
# 
#     ##Riou -> donne des valeur negatives pour forts epsi
#     # coeffTV = epsi/ (1.-leafAlbedo)
#     # EP = (1. - coeffTV) * Et0#potentialSoilEvaporation
# 
#     ##adaptation eq 7.1 p127 STICS book
#     k_eq = 0.7  # LAI equivalent pour un k=0.7; pause k car les deux k et LAI sont inconnu
#     LAIeq = -log(1 - epsi) / k_eq
#     EP = Et0 * exp(-(k_eq - 0.2) * LAIeq)  # potentialSoilEvaporation Eq. 7.1, p127
#     # rq: la bonne variable a recuperer a la place de espi serait le taux de couverture ou le transmis vertical pour lequel le correctif 0.2 n'est pas necessaire (Eq. 7.2)
# 
#     if previous_state[0] + EP < U:  # phase 1
#         SES1 = previous_state[0] + EP
#         SES2 = 0.
#         ER = EP
#         SEP = 0.
#         # mise a jour state!
#     elif previous_state[0] < U and previous_state[0] + EP >= U:  # transition phase1->2
#         fracEP2 = previous_state[0] + EP - U  # precipitation evaporee en phase 2
#         fracEP1 = EP - fracEP2
#         SES1 = previous_state[0] + EP
#         SES2 = fracEP2
#         ER2 = ((2 * b * SES2 + b ** 2) ** 0.5) - b
#         ER = fracEP1 + (fracEP2 / EP) * ER2  # verif
#         SEP = ER2
#     else:  # phase 2
#         SES1 = previous_state[0] + EP
#         SES2 = previous_state[1] + EP
#         ER = ((2 * b * SES2 + b ** 2) ** 0.5) - b - previous_state[
#             2]  # fonction cumul, donc retire valeur de la date precedente
#         SEP = ((2 * b * SES2 + b ** 2) ** 0.5) - b
# 
#     # MAJ memoire si pluie - comme si pluie vient en fin de journee, apres calcul ER
#     seuil = U  # 0.7*U
#     if Precip > 0. and Precip >= seuil and SES1 > U:  # si en phase 2 et forte precipitation
#         SES1 = 0.
#         SES2 = 0.
#         SEP = 0.
#     elif Precip > 0. and Precip < seuil and SES1 > U:  # si en phase 2 et faible precipitation
#         # remet en debut de phase 2
#         SES1 = SES1 - SES2
#         SES2 = 0.
#         SEP = 0.
#     elif Precip > 0. and SES1 <= U:  # si en phase 1
#         SES1 = 0.
#         SES2 = 0.
#         SEP = 0.
# 
#     return ER, [SES1, SES2, SEP]  # 0., [SES1, SES2, 0.]#
#     # a tester pour limiter effet peites pluie - effet trop fort! + seuil arbitraire avec tres fort impact


#def evap_frac_zone(dz = [0., 0.3], ZESX=0.3):
#    """ distribution en profondeur du prelevement d'eau evapore sous un voxel de surface
#        integrale de y = -ZESX*ZESX*Z/2 + ZESX -> relation lineaire de la figure 7.4b STICS manual (p130)
#        correspond a equation 7.6 pour CFES=1 (lineaire)"""
#    if dz[0]>=0. and dz[0]<=ZESX: #1ere borne OK
#        if dz[1]>0. and dz[1]<=ZESX: #2e borne OK
#            frac_interval = (-(1./ZESX**2)*(dz[1]**2)+(2./ZESX)*dz[1]) - (-(1./ZESX**2)*(dz[0]**2)+(2./ZESX)*dz[0])
#        elif dz[1]>0. and dz[1]>ZESX: #2e superieure-> remplace par profondeur max ZESX
#            frac_interval = (-(1./ZESX**2)*(ZESX**2)+(2./ZESX)*ZESX) - (-(1./ZESX**2)*(dz[0]**2)+(2./ZESX)*dz[0])
#        else:
#            frac_interval = 0.
#    else:
#        frac_interval = 0.
#
#    return frac_interval

#evap_frac_zone(dz = [0., 0.2], ZESX=0.3)


#def check_soil_evap(evapo_tot, m_frac_evap, tsw_t, m_QH20min, dxyz, ZESX = 0.3):
#    """ verifie que evaporation du sol n'entraine pas de tsw sous la tsw mini d'un sol sec (air) ; sinon ajuste distribution en profondeur(report sur la couche inferieure), puis evapo_tot"""
#    m_evap = deepcopy(m_frac_evap)
#
#    #recherche des id couches sur lesquelles evaporation joue
#    nz, cumz = 0, []
#    for i in range(len(dxyz[2])):
#        cumz.append(dxyz[2][i])
#        nz = nz+1
#        if sum(cumz)>=ZESX:
#            break 
#
#    for x in range(len(dxyz[0])):
#        for y in range(len(dxyz[1])):
#            for z in range(nz-1):#report par colone de sol
#                #res = tsw_t[z][x][y] - m_evap[z][x][y]#m_frac_evap[z][x][y] #- report
#                res = tsw_t[z][x][y] - m_QH20min[z][x][y] - m_evap[z][x][y]#m_frac_evap[z][x][y] #- report
#                if res<=0.:
#                    #m_evap[z+1][x][y] = m_evap[z+1][x][y] + m_evap[z][x][y] - tsw_t[z][x][y]
#                    #m_evap[z][x][y] = tsw_t[z][x][y]
#                    m_evap[z+1][x][y] = m_evap[z+1][x][y]-m_QH20min[z+1][x][y] + m_evap[z][x][y]-m_QH20min[z][x][y] - tsw_t[z][x][y]#??
#                    m_evap[z][x][y] = tsw_t[z][x][y]- m_QH20min[z][x][y]#??
#
#            #derniere couche affectee
#            z = nz-1
#            res = tsw_t[z][x][y] - m_QH20min[z][x][y] - m_evap[z][x][y]#m_frac_evap[z][x][y] #- report
#            if res<=0.:
#                m_evap[z][x][y] = tsw_t[z][x][y] - m_QH20min[z][x][y]
#
#    return sum3(m_evap), m_evap

#test = deepcopy(S.tsw_t)
#test[0][0][0]=0.
#test[1][0][0]=0.
#test[2][0][0]=0.
#test[3][0][0]=0.
#test[4][0][0]=0.#002
#test[5][0][0]=0.#002
#test[6][0][0]=0.#002
#test[7][0][0]=0.002
#evapo_tot,  m_evap = check_soil_evap(sum3(m_frac_evap), m_frac_evap, test, S.dxyz, ZESX = 0.3)


#def distrib_evapNC(m_soil_vox, m_soil_vol, map_Evap0, ZESX=0.3):
#    """ map_Evap0 = liste des Evap des voxel de la couche superieure (matrice x, y) """
#    m = deepcopy(m_soil_vol)
#    if len(m)>1:
#        m.fill(0.)
#        for z in range(len(m)-1):
#            dz = [abs(m_soil_vox[z][0][0][2]), abs(m_soil_vox[z+1][0][0][2])]
#            frac = evap_frac_zone(dz , ZESX)
#            for x in range(len(m[z])):
#                for y in range(len(m[z][x])):
#                    m[z][x][y] = frac*map_Evap0[x][y]
#    else:
#        m.fill(1.)
#
#    return m

#map_Evap0 = array([[1.5]]) ## faire une fonction qui le genere initialement
#m_frac_evap = distrib_evapNC(m_soil_vox, m_soil_vol, map_Evap0, ZESX=0.3)





########## diverses fonction - soil nitrogen balance



## soil N supply


def Convective_Nflux(SN, ls_mWaterUptakePlt):
    """
        
    """
    #""" Eq. 8_34 p 160 - potential passive NO3 uptake with water flow"""
    ls_conv = []
    convtot = SN.m_1*0.
    for i in range(len(ls_mWaterUptakePlt)):
        flux_c = ls_mWaterUptakePlt[i] * SN.ConcNO3()
        ls_conv.append(flux_c)
        convtot = convtot+flux_c

    return convtot, ls_conv
    #renvoyer liste de fraction?

def Diffusive_Nflux(SN, parSN, ls_lrac):
    """
        
    """
    #""" Eq. 8.35 p 161 """
    #ls_lrac = ls_ roots (en m) -> convert en cm (*100)

    ls_diff = []
    difftot = SN.m_1*0.

    FH = (SN.tsw_t - SN.m_QH20wp) / (SN.m_QH20fc - SN.m_QH20wp)
    filtre = FH>0
    filtre = filtre * 1. #remplace valeur negatives par zero
    DIFE = parSN['DIFNg'] * FH * filtre
    coeff = 4*np.pi**0.5
    for i in range(len(ls_lrac)):
        draci = (ls_lrac[i]*100.) / (SN.m_soil_vol*1000000.) #root density (cm.cm-3)
        flux_d = coeff * DIFE * (SN.m_NO3 + SN.m_NH4) *np.sqrt(draci)
        ls_diff.append(flux_d)
        difftot = difftot + flux_d

    return difftot, ls_diff
    #a priori ls_diff pas utilise


def Nsoil_supply(SN, parSN, ls_lrac, ls_mWaterUptakePlt):
    """
        
    """
    #""" Eq. 8.33 p 160 """
    convtot, ls_conv = Convective_Nflux(SN, ls_mWaterUptakePlt)
    difftot, ls_diff = Diffusive_Nflux(SN, parSN, ls_lrac)
    ls_suptot = []
    for i in range(len(ls_conv)):
        ls_suptot.append(ls_conv[i] + ls_diff[i])

    return convtot+difftot, ls_suptot
    #a priori ls_suptot pas utilise

