#from scipy import *
import IOtable
import IOxls
import string
from copy import deepcopy
try:
    from riri5 import RIRI5 as riri #import de la version develop si module soil3ds est installe
except:
    import riri5.RIRI5 as riri

import numpy as np
import pandas as pd


#Temperature response funtions
def betaT(Tmin, Tmax, q, T):
    """ beta de Graux (2011)"""
    Tref = 20
    if T < 0.:  # zero sous zero degres
        fT = 0.
    else:
        fT = ((Tmax - T) / (Tmax - Tref)) * ((T - Tmin) / (Tref - Tmin)) ** q

    return max(fT, 0.)

def TempHoraire(h, Tmin, Tmax):
    """ interpolation cos des temperatures horaires"""
    """ Eq S33a - Evers et al 2010"""
    """ Temperature moyenne journaliere = (Tmin+Tmax/2)"""
    Ta = 0.5*((Tmax+Tmin)+(Tmax+Tmin)*np.cos(np.pi*(h+8.)/12.))
    return Ta


def dTT(vT, p, optT=0):
    """ fonction de cumul du temp thermique; integre reponse non lineaire"""
    """ 3 options: 0=betaD journalier ; 1=betaH Horaire; 2=Tbase lineaire ; vT liste des temperatures"""
    Tref = 20. #reference temperature = 20 degreC
    if optT==0: #betaD
        T = np.mean(vT)
        return max((Tref - p[0]) * betaT(p[1], p[2], p[3], T), 0.)
    elif optT==1: #betaH
        ls_betaH = []
        for T in vT:
            ls_betaH.append(betaT(p[1], p[2], p[3], T))

        return max((Tref - p[0]) * np.mean(ls_betaH), 0.)
    elif optT==2:#Tbase lineaire
        T = np.mean(vT)
        return max((T - p[0]), 0.)
    elif optT==3:#olg bug
        #just for debuging and testing conformity with oler model verions!!! do not use
        T = np.mean(vT)
        return max((T - p[0]) * betaT(p[1], p[2], p[3], T), 0.)


def Calc_Daily_vT(meteo_j, opt_optT):
    # extract dily vector of temperature depending on TT calculation option
    if opt_optT ==0 or opt_optT==2 or opt_optT==3: #betaD ou lineaire ou old debug
        vT = [meteo_j['TmoyDay']]
        vTsol = [meteo_j['Tsol']]
    elif opt_optT ==1: #betaH
        #calcul vecteur des temperatures horaires air
        vT = []
        for h in range(24):
            vT.append(TempHoraire(h, Tmin=meteo_j['Tmin'], Tmax=meteo_j['Tmax']))
        vTsol = [meteo_j['Tsol']]

    return vT, vTsol
    #vT, vTsol = Calc_Daily_vT(meteo_j, opt_optT)

#Light response functions
def DecliSun(DOY):
    """ Declinaison (rad) du soleil en fonction du jour de l'annee """
    alpha = 2. * 3.14 * (DOY - 1) / 365.
    return (0.006918 - 0.399912 * np.cos(alpha) + 0.070257 * np.sin(alpha))

def DayLength(latitude, decli):
    """ photoperiode en fonction de latitude (degre) et declinaison du soleil (rad) """
    lat = np.radians(latitude)
    d = np.arccos(-np.tan(decli) * np.tan(lat))
    if d < 0:
        d = d + 3.14

    # d en 'hour angle'
    return 24. * d / 3.14  # conversion en heures

def trilineaire(x, ratio0, ratiomax, parmaxeff, parnoeff):
    # photomorphogenese
    pentor = (ratiomax - ratio0) / parmaxeff
    pentfin = (1 - ratiomax) / (parnoeff - parmaxeff)
    orlindesc = -(parnoeff * pentfin) + 1
    return min(pentor * x + ratio0, max(1, pentfin * x + orlindesc))

def monomoleculaire(x, Amax, k):
    #reponse saturante type mononmoleculaire
    rate = Amax * (1 - np.exp(-k*x))
    return rate

#general growth functions
def expansion(t, a, delai):
    "croissance sigmoidale"
    return 1/(1 + np.exp(-a*(t-delai)))

def sigmo_stress(v,delai,x):
    "reponse sigmo a stress - FTSW ou INN"
    return 1-1/(1 + np.exp(v*(x-delai)))

def linear_stress(tresh, x):
    "linear response between 0 and tresh - 1 above - FTSW ou INN"
    if x>=tresh:
        resp = 1.
    else:
        resp = (1/tresh)*x
    return resp
    #e.g. response of tranpsiration below FTSW=0.6

def linear_stress2(x, par, opt=0):#treshbas, treshhaut,
    "linear response between 0 and 1, from tresh bas to treshhaut; opt 0: monte; opt 1: baisse"
    treshbas, treshhaut = par[0], par[1]
    if x<=treshbas:
        resp = 0.
    elif x>=treshhaut:
        resp = 1.
    else:
        resp = (1/(treshhaut - treshbas ))* (x-treshbas)

    if opt==1: #reponse descendante
        resp = 1. - resp

    return resp
    #e.g. response to photoperiod
    #linear_stress2(8, 13, 10, opt=0)
    #linear_stress2( 10,[8, 13], opt=0)


# N response functions
def Na_N0(I_I0):
    """ teneur en azote relative en fonction de fraction d'eclairement - Louarn et al. 2014 """
    return I_I0**0.247

def N0(INN):
    """ teneur en azote des feuilles eclairees (g N.m-2) en fonction de INN - Louarn et al. 2014 """
    return 2.08*INN+0.05

def Nl_Nl0(I_I0):
    """ teneur en N lineique relative des tiges """
    Nresidu = 0.17
    return Nresidu + (1-Nresidu)*I_I0**0.51

def NNI_resp(NNI, par):
    # """Belanger, Gastal and lemaire 1992 - RUE/RUEmax = f(Nitrogen Nutrition Index)"""
    # RUE_Eff= max(0,1.05*(1-2.18*exp(-3.13*NNI)))

    res = 1.
    if NNI < 0.95:
        res = sigmo_stress(par[0], par[1], NNI)

    return res  # RUE_Eff

def Ndfa_max(ageTT, DurDevFix, Delfix=100.):
    """ Ndfa (prop d'N issu de fixation) max depend de stade - demare a 100 degree.days (Delfix -> a remonter en parametere)"""
    Delfix = 0.
    slope = 1. / (np.array(DurDevFix) - Delfix)
    ordoOr = -Delfix * slope
    val = slope * np.array(ageTT) + ordoOr
    val[val > 1.] = 1.
    val[val < 0.] = 0.
    return val
    # ages = [100., 200.]
    # DurDevFix = [600., 100.]
    # Ndfa_max(ages, DurDevFix)

def ActualFix(ls_demand, Nuptakes, MaxFix):
    """ calcul fixation a partir de demande, prelev mineral (prioritaire) et MaxFix """
    demande_residuelle = ls_demand - Nuptakes
    fix = []
    for i in range(len(demande_residuelle)):
        fix.append(min(MaxFix[i], demande_residuelle[i]))

    fix = np.array(fix)
    fix[fix < 0] = 0.
    return fix
    # ls_demand = array([1.,2.,3.])
    # Nuptakes = array([1.,1.,1.])
    # MaxFix = array([3.,3.,3.])
    # ActualFix(ls_demand, Nuptakes, MaxFix)


# water stress functions
def FTSW_resp(FTSW, par):  # =[0.4, 0.]):
    # """ facteur de reduction en reponse a FTSW - lineaire entre 0 et 0.4"""
    """ reponse sigmo """
    res = 1.
    # if FTSW <= par[1]:
    #    res=0.
    # elif FTSW > par[1] and FTSW < par[0]:
    #    res=FTSW/(par[0]-par[1])
    if FTSW <= 0.95:
        res = sigmo_stress(par[0], par[1], FTSW)

    return res


#calcul des longueurs et surfaces d'organes
def calc_surF(ParamP, rank, rankp, ordre, l, type=1):
    """ calcul de surface d'une feuille (m2) """
    cor_ordre = ParamP['ratioII'] if ordre == 2 else 1.
    if int(type)==1 or int(type)==2: #legume leaves
        rk = rank + rankp if ordre == 2 else rank
    elif int(type)==3: #graminee
        rk = rank + rankp #rankp utilise pour calcule un rgeq qui tient compte de tallage et coupe

    rk = min(rk, len(ParamP['profilLeafI_l']) - 1)  # au cas ou rank depasse le profil
    nf = ParamP['profilLeafI_nfol'][rk]
    Long = ParamP['profilLeafI_l'][rk] * l * cor_ordre
    larg = ParamP['profilLeafI_larg'][rk] * l * cor_ordre
    if int(type)==1 or int(type)==2: #feuille legumineuse -> losange
        surF = nf * 0.5 * Long * larg / 10000.  # nf*0.5*Long*larg/10000.  #m2
    elif int(type)==3: #graminee -> rectangle
        surF =  Long * larg / 10000.  # nf*0.5*Long*larg/10000.  #m2
    return surF
    # ParamP en parametre

def calc_surS(ParamP, rank, rankp, ordre, l):
    """ calcul de surface d'une stipule (m2) """
    cor_ordre = ParamP['ratioII'] if ordre == 2 else 1.
    rk = rank + rankp if ordre == 2 else rank
    rk = min(rk, len(ParamP['profilStipI_l']) - 1)  # au cas ou rank depasse le profil
    Long = ParamP['profilStipI_l'][rk] * l * cor_ordre
    larg = ParamP['profilStipI_larg'][rk] * l * cor_ordre
    surF = 2 * Long * larg * 0.5 / 10000.  # m2
    return surF
    # ParamP en parametre

def calc_surfcoty(Mcoty, age, DurGraine, carto, ParamP, n_gamagroup, origin_grid, na, dxyz, SLAcoty=600.):
    """ distribution de surface coty dans grille 3D - depend de masse de coty et SLAcoty et graine; et zero apres DurGraine"""
    # valeur de 600 tiree essai RGR2015
    # peut passer SLAcoty en parametre et variable par plante
    mcot = np.zeros([n_gamagroup, na[2], na[1], na[0]])

    for nump in range(len(carto)):
        vox = riri.WhichVoxel(np.array(carto[nump]), origin_grid, na, dxyz)
        if age[nump] <= DurGraine[nump]:  # cotyledons actifs pendant DurGraine
            surfcot = Mcoty[nump] * SLAcoty / 10000.  # m2
        else:
            surfcot = 0.

        mcot[ParamP[nump]['id_grid']][vox[2]][vox[1]][vox[0]] += surfcot

    return mcot

def calc_parapcoty(invar, m_lais, res_abs_i, Mcoty, age, DurGraine, carto, ParamP, n_gamagroup, origin_grid, na, dxyz, SLAcoty=100.):
    """ ajout des PARa des cotyledon a invar['PARiPlante']"""
    for nump in range(len(carto)):
        vox = riri.WhichVoxel(np.array(carto[nump]), origin_grid, na, dxyz)
        sVOX = m_lais[ParamP[nump]['id_grid']][vox[2]][vox[1]][vox[0]]
        if age[nump] <= DurGraine[nump]:  # cotyledons actifs pendant DurGraine
            surfcot = Mcoty[nump] * SLAcoty / 10000.  # m2
        else:
            surfcot = 0.

        if sVOX > 0.:
            PARaF = res_abs_i[ParamP[nump]['id_grid']][vox[2]][vox[1]][vox[0]] * surfcot / sVOX * 3600. * 24 / 1000000.
        else:
            PARaF = 0.

        invar['PARiPlante'][nump].append(PARaF)
        #m_laiPlt[nump][vox[2]][vox[1]][vox[0]] += surfcot
        #lsFeuilBilanR.append([nump, 0, 0, 0, 'coty', surfcot, ParamP[nump]['id_grid'], carto[nump][0], carto[nump][1], carto[nump][2], vox[2], vox[1], vox[0], 0, 0])

        # Plus utilise...

def add_surfcoty(invar, m_lais, m_laiPlt, lsFeuilBilanR, carto, ParamP, origin_grid, na, dxyz, SLAcoty=100.):
    """ ajout de surf cotyledon a m_lais, m_laiPlt, lsFeuilBilanR """
    age = invar['TT']
    Mcoty = invar['MS_coty']
    DurGraine = IOxls.get_lsparami(ParamP, 'DurGraine')

    for nump in range(len(carto)):
        vox = riri.WhichVoxel(np.array(carto[nump]), origin_grid, na, dxyz)
        if age[nump] <= DurGraine[nump]:  # cotyledons actifs pendant DurGraine
            surfcot = Mcoty[nump] * SLAcoty / 10000.  # m2
        else:
            surfcot = 0.

        #m_lais[ParamP[nump]['id_grid']][vox[2]][vox[1]][vox[0]] += surfcot
        #m_laiPlt[nump][vox[2]][vox[1]][vox[0]] += surfcot
        if surfcot > 0.:
            m_lais[ParamP[nump]['id_grid']][vox[2]][vox[1]][vox[0]] += surfcot
            m_laiPlt[nump][vox[2]][vox[1]][vox[0]] += surfcot
            lsFeuilBilanR.append([nump, 0, 0, 0, 'coty', surfcot, ParamP[nump]['id_grid'], carto[nump][0], carto[nump][1], carto[nump][2], vox[2], vox[1], vox[0], 0, 0])


def calc_Lpet(ParamP, rank, rankp, ordre, l, type=1):
    """ calcule de longueur potentielle d'un petiole (m)"""
    cor_ordre = ParamP['ratioII'] if ordre==2 else 1.
    if int(type)==1 or int(type)==2: #legume leaves
        rk = rank + rankp if ordre == 2 else rank
    elif int(type)==3: #graminee
        rk = rank + rankp #rankp utilise pour calcule un rgeq qui tient compte de tallage et coupe

    rk = min(rk, len(ParamP['profilPetI_l'])-1) #au cas ou rank depasse le profil
    lpet = l*ParamP['profilPetI_l'][rk]*cor_ordre/100. #m
    return lpet

def calc_Lent(ParamP, rank, nsh, ordre, l):
    """ calcule de longueur poteielle d'un entre noeud (m)"""
    cor_M = ParamP['ratioM'] if nsh==0 and ordre==1 else 1. #correction tige seminale
    cor_ordre = ParamP['ratioII'] if ordre==2 else 1.
    rank = min(rank, len(ParamP['profilNodeI_l'])-1) #au cas ou rank depasse le profil
    lent = l * ParamP['profilNodeI_l'][rank] * cor_ordre * cor_M/100. #m
    return lent

def calcSurfScale(ParamP, tab, scale):
    """ calcul le cumul de surface foliaire a l'echelle indiquee """
    dp = {}  # dictionnaire a l'echelle choisie: plante/shoot/axe
    for i in range(len(tab['nump'])):
        if scale == 'plt':
            idp = str(tab['nump'][i])
        elif scale == 'sh':
            idp = str(tab['nump'][i]) + '_' + str(tab['nsh'][i])
        elif scale == 'ax':
            idp = str(tab['nump'][i]) + '_' + str(tab['nsh'][i]) + '_' + str(tab['rank'][i])

        age = float(tab['age'][i])
        nump = int(tab['nump'][i])
        ordre = int(tab['ordre'][i])
        rank = int(tab['rank'][i])
        rankp = int(tab['rankp'][i])
        nsh = int(tab['nsh'][i])
        l = float(tab['l'][i])

        surf = 0.
        if tab['organ'][i] == 'Lf':
            surf = calc_surF(ParamP[nump], rank, rankp, ordre, l)  # m2

        if tab['organ'][i] == 'Stp':
            surf = calc_surS(ParamP[nump], rank, rankp, ordre, l)  # m2

        # ajoute une cle pour chaque organe (meme si pas feuille)
        try:
            dp[idp].append(surf)
        except:
            dp[idp] = [surf]

    for k in list(dp.keys()):
        dp[k] = sum(dp[k])

    return dp
    # lsSurfSh = calcSurfScale(IOtable.conv_dataframe(IOtable.t_list(lsOrgans)), 'sh')
    # plus utile car fait dans calcSurfLightScales


def calcSurfLightScales(tab, ParamP):
    """ calcul le cumul de surface foliaire, surface foliaire verte et PARa au echelles plante/shoot_axe - passe la table organe une seulle fois en revue """
    # from dic lsFeuilBilanR : ['nump', 'nsh', 'rank', 'rankp','status', 'surf', 'id_grid', 'X','Y','Z','Vox2','Vox1','Vox0','sVox','paraF']

    dpS, dpSV, dpPARaF, dshS, dshSV, dshPARaF, daxS, daxSV, daxPARaF, daxPARaFsurf = {}, {}, {}, {}, {}, {}, {}, {}, {}, {}
    for i in range(len(tab['nump'])):

        idp = str(tab['nump'][i])
        idsh = str(tab['nump'][i]) + '_' + str(tab['nsh'][i])
        if int(ParamP[tab['nump'][i]]['type']) == 3 :
            #grass -> talle 3 F
            idax = str(tab['nump'][i]) + '_' + str(tab['nsh'][i]) + '_' + str(0)
        else:
            #legume
            idax = str(tab['nump'][i]) + '_' + str(tab['nsh'][i]) + '_' + str(tab['rank'][i])

        PARaF = float(tab['paraF'][i])
        surf = max(float(tab['surf'][i]), 10e-15)  # m2 #max(calc_surF(ParamP[nump], rank, rankp, ordre, l), 10e-15)  # m2
        PARaFsurf = PARaF / (surf * 3600. * 24 / 1000000.)  #
        if tab['status'][i] != 'sen':
            surfV = surf
        else:
            surfV = 0.


        # ajoute une cle pour chaque organe (meme si pas feuille)
        IOxls.append_dic(dpS, idp, surf)
        IOxls.append_dic(dshS, idsh, surf)
        IOxls.append_dic(daxS, idax, surf)
        IOxls.append_dic(dpSV, idp, surfV)
        IOxls.append_dic(dshSV, idsh, surfV)
        IOxls.append_dic(daxSV, idax, surfV)
        IOxls.append_dic(dpPARaF, idp, PARaF)
        IOxls.append_dic(dshPARaF, idsh, PARaF)
        IOxls.append_dic(daxPARaF, idax, PARaF)
        IOxls.append_dic(daxPARaFsurf, idax, PARaFsurf)

    IOxls.sum_ls_dic(dpS)  # surface par plante
    IOxls.sum_ls_dic(dpSV)  # surface par shoot
    IOxls.sum_ls_dic(dshS)  # surface par axe
    IOxls.sum_ls_dic(dshSV)  # surface verte par plante
    IOxls.sum_ls_dic(daxS)  # surface verte par shoot
    IOxls.sum_ls_dic(daxSV)  # surface verte par axe
    IOxls.sum_ls_dic(dpPARaF)  # PARaF par plante
    IOxls.sum_ls_dic(dshPARaF)  # PARaF par shoot
    IOxls.sum_ls_dic(daxPARaF)  # PARaF par axe
    for k in list(daxPARaFsurf.keys()):
        daxPARaFsurf[k] = max(daxPARaFsurf[k])

    return dpS, dpSV, dshS, dshSV, daxS, daxSV, dpPARaF, dshPARaF, daxPARaF, daxPARaFsurf



def AgePivScales(tab, ParamP):
    """ calculage pivot/RacI pour """
    # from dic lsOrgans: ['TT','organ','nump', 'nsh', 'rank', 'rankp', 'strate', 'surf', 'PARaF','statut','age','ordre','l','Long','DOY','cutNB', 'Larg']

    daxAgePiv = {}  # dico des axes qui portent des pivots et de leur age
    for i in range(len(tab['nump'])):

        idp = str(tab['nump'][i])
        idsh = str(tab['nump'][i]) + '_' + str(tab['nsh'][i])
        if int(ParamP[tab['nump'][i]]['type']) == 3 :
            #grass -> talle 3 F
            idax = str(tab['nump'][i]) + '_' + str(tab['nsh'][i]) + '_' + str(0)
        else:
            #legume
            idax = str(tab['nump'][i]) + '_' + str(tab['nsh'][i]) + '_' + str(tab['rank'][i])

        age = float(tab['age'][i])

        if tab['organ'][i] == 'Piv':
            daxAgePiv[idax] = age

    return daxAgePiv
    #separe de calcSurfLightScales



#germination / iitialisation funcions
def germinate(invar, ParamP, nump):
    # mis a jour pour chaque graine a germination
    # creation des cotyledons + reste des reserves

    frac_coty_ini = ParamP['frac_coty_ini']
    invar['MS_coty'][nump] = invar['MS_graine'][nump] * frac_coty_ini
    invar['Ncoty'][nump] = invar['MS_graine'][nump] * frac_coty_ini * ParamP['Npc_ini'] / 100.
    #tout le reste initialement dans reserve

    # met a jour graine qui a germe et defini pools de reserve (ce qui reste dans MSgraine et N graine = reserve pour soutenir croissance ini)
    invar['MS_graine'][nump] -= invar['MS_coty'][nump]
    invar['Ngraine'][nump] -= invar['MS_coty'][nump] * ParamP['Npc_ini'] / 100.

    #rq: de cette facon une grosse partie des reserve N n'est pas utilisee (immobilisee ds coty)
    #je me demande si c pas la le pb de 0???

    #deux constates reutilisee ensuite
    invar['dMSgraine'][nump] = invar['MS_graine'][nump] / ParamP['DurGraine']  # delta MS fourni par degrejour par graine pendant DurGraine
    invar['dNgraine'][nump] = invar['Ngraine'][nump] / ParamP['DurGraine']  # delta QN fourni par degrejour par graine pendant DurGraine

    # cotyledons meurent quand DurGraine atteint -> cf calc_surfcoty
    # en toute logique pas besoin de mettre a jour Mfeuil?

def reserves_graine(invar, ParamP):
    """ calcul des flux venant des reserves de graine et met a jour MS_graine et Ngraine """

    graineC, graineN = [], []
    for nump in range(len(ParamP)):
        dTT = invar['Udev'][nump]  #invar['dTT'][nump]
        if invar['TTudev'][nump] < ParamP[nump]['DurGraine'] and invar['TTudev'][nump] > 0.:
            # suppose consommation reguliere pendant DurGraine
            dMSgraine = invar['dMSgraine'][nump] * dTT
            dNgraine = invar['dNgraine'][nump] * dTT
        else:
            dMSgraine = 0.
            dNgraine = 0.

        #pour eviter negatif
        if dMSgraine > invar['MS_graine'][nump]:
            dMSgraine = invar['MS_graine'][nump]

        if dNgraine > invar['Ngraine'][nump]:
            dNgraine = invar['Ngraine'][nump]

        graineC.append(dMSgraine)
        graineN.append(dNgraine)

    invar['MS_graine'] = np.array(invar['MS_graine']) - np.array(graineC)
    invar['Ngraine'] = np.array(invar['Ngraine']) - np.array(graineN)

    return np.array(graineC), np.array(graineN)


def calc_paraF(dicFeuilBilanR, m_lais, res_abs_i, force_id_grid=None):
    """ update paraF and sVox calculation in lsFeuilBilanR + conv to dico """
    # force_id_grid=n is for external calculation for a single species in a mixture
    # ['nump', 'nsh', 'rank', 'rankp','status', 'surf','id_grid', 'X','Y','Z','Vox2','Vox1','Vox0','sVox','paraF']
    nblignes = len(dicFeuilBilanR['nump']) #- 1
    df = dicFeuilBilanR #IOtable.conv_dataframe(IOtable.t_list(lsFeuilBilanR))
    if nblignes > 0:
        for i in range(nblignes):
            # update sVox
            if (force_id_grid is None):
                id_grid = df['id_grid'][i]
            else:
                id_grid = int(force_id_grid)

            surf = df['surf'][i]
            vox = [df['Vox0'][i], df['Vox1'][i], df['Vox2'][i]]
            sVOX = m_lais[id_grid][vox[2]][vox[1]][vox[0]]
            paraF = res_abs_i[id_grid][vox[2]][vox[1]][vox[0]] * surf / sVOX * 3600. * 24 / 1000000.

            # mise a jour
            df['sVox'][i] = sVOX
            df['paraF'][i] = paraF

    return df


def calc_para_Plt(invar, lsFeuilBilanR):
    """ update invar plant light interception + conv to pd data.frame """
    # conversion data.frame
    lsFeuilBilanR = pd.DataFrame(lsFeuilBilanR)
    nbplt = len(invar['TT'])

    # split pour PARi
    dfPARa = lsFeuilBilanR.groupby(lsFeuilBilanR['nump'])
    groupOK = list(dfPARa.groups.keys())
    pari = []
    for nump in range(nbplt):
        if nump in groupOK:
            x = dfPARa.get_group(nump)
            pari.append(np.sum(x['paraF']))
        else:
            pari.append(0.)

    # split pour PARa (sans feuilles sen)
    lsFeuilBilanR_sen = lsFeuilBilanR[lsFeuilBilanR["status"] != 'sen']
    dfPARa = lsFeuilBilanR_sen.groupby(lsFeuilBilanR_sen['nump'])
    groupOK = list(dfPARa.groups.keys())
    para = []
    for nump in range(nbplt):
        if nump in groupOK:
            x = dfPARa.get_group(nump)
            para.append(np.sum(x['paraF']))
        else:
            para.append(0.)

    # MAJ de invar
    invar['parip'] = np.array(pari)
    invar['parap'] = np.array(para)

    # return lsFeuilBilanR



def Turnover_compart_Perenne(invar, ParamP):
    """ calcul des flux lie a turnover des tissus d'un comaprtiment perennes """

    dMSenNonRec, dMSenPiv = [], []
    perteN_NonRec, perteN_Piv = [], []
    minPiv = 0.001
    for nump in range(len(ParamP)):
        dTT = invar['Udev'][nump]  # invar['dTT'][nump]
        delai_senperenne = ParamP[nump]['delai_senperenne']  # 500. #parametre -> rq: jouer sur ce parametre pour faire mourrir piv TV?
        if invar['TTudev'][nump] > delai_senperenne and invar['aliveB'][nump] == 0: #delai passe et vivante
            TOrate_nonrec = ParamP[nump]['TOrate_nonrec'] #0.005 # a passer en parametre
            TOrate_piv = ParamP[nump]['TOrate_piv'] #0.005  # a passer en parametre
            dMS_aerienNonRec = invar['MS_aerienNonRec'][nump] * TOrate_nonrec * dTT
            if invar['MS_pivot'][nump] > minPiv:
                #TOpivot au dessus d'une taille mini
                dMS_piv = invar['MS_pivot'][nump] * TOrate_piv * dTT
            else:
                dMS_piv = 0.

            perteN_NonRec_i = invar['NaerienNonRec'][nump] * TOrate_nonrec * dTT
            perteN_Piv_i = invar['Npivot'][nump] * TOrate_piv * dTT
        else:
            dMS_aerienNonRec = 0.
            dMS_piv = 0.
            perteN_NonRec_i = 0.
            perteN_Piv_i = 0.

        dMSenNonRec.append(dMS_aerienNonRec)
        dMSenPiv.append(dMS_piv)
        perteN_NonRec.append(perteN_NonRec_i)
        perteN_Piv.append(perteN_Piv_i)

    #MAJ des compartiments
    invar['MS_aerienNonRec'] -= np.array(dMSenNonRec)
    invar['MS_pivot'] = np.array(invar['MS_pivot'])
    invar['MS_pivot'] -= np.array(dMSenPiv) #faire un MSpiv_net?
    invar['MS_pivot'] = invar['MS_pivot'].tolist()

    invar['NaerienNonRec'] -= np.array(perteN_NonRec)
    invar['Npivot'] -= np.array(perteN_Piv)

    #renvoie les flux
    return np.array(dMSenNonRec), np.array(dMSenPiv), np.array(perteN_NonRec), np.array(perteN_Piv)




# Carbon allocation
def rootalloc(parB, parA, SB):
    """ calcule fraction d'alloc racine/shoots en fonction du cumule shoots - Eq. 8 draft article V Migault"""
    """ puis concerti en fraction d'allocation de biomasse totale produite au racine dRB/dMStot a partir du ratio SRB/dSB"""
    #params = [parB, parA]
    nbplantes = len(SB)
    res = [0] * nbplantes
    for nump in range(nbplantes):
        bet = parB[nump]
        alph = parA[nump]
        res[nump] = min(bet, bet * alph * max(SB[nump], 0.00000000001) ** (alph - 1))  # epsilon evitant de calculer un 0 avec puissance negative (cause erreur). Le maximum possible est pour alpha=1, donc beta*1*SB**(1-1) = beta * 1 * 1 = beta.

    dRB_dSB = np.array(res)
    return dRB_dSB / (1 + dRB_dSB)

def calcOffreC(ParamP, tab, scale):
    """ calcul de cumul de para par plante ('plt')/tige('sh')/axe('ax') * RUE-> offreC """
    """ tab = attend dico avec cle ('nump', 'nsh', 'rank', 'PARaF','statut','age','ordre') """
    dp = {}  # dictionnaire a l'echelle choisie: plnate/shoot/axe
    for i in range(len(tab['nump'])):
        if scale == 'plt':
            idp = str(tab['nump'][i])
        elif scale == 'sh':
            idp = str(tab['nump'][i]) + '_' + str(tab['nsh'][i])
        elif scale == 'ax':
            idp = str(tab['nump'][i]) + '_' + str(tab['nsh'][i]) + '_' + str(tab['rank'][i])

        if tab['organ'][i] == 'Lf' and tab['statut'][i] != 'sen':
            nump = int(tab['nump'][i])
            try:
                dp[idp].append(float(tab['PARaF'][i]) * ParamP[nump]['RUE'])
            except:
                dp[idp] = [float(tab['PARaF'][i]) * ParamP[nump]['RUE']]

    for k in list(dp.keys()):
        dp[k] = sum(dp[k])

    return dp
    # pas utilise dans version actuelle
    # approche RUE limite a echelle feuille! ; garder aussi 'sen'?

def calcDemandeC(ParamP, tab, scale, udev, ls_ftswStress, ls_NNIStress):
    """ calcul de demande pour assurer croissance potentielle minimale des Lf(), In() et Pet()en phase d'expansion """
    # distingue les Lf et Stp! car pas meme calcul de surface!
    dp = {}  # dictionnaire a l'echelle choisie: plante/shoot/axe #-> demande tot
    dplf = {}  # demande des feuilles
    dpin = {} #demande des En
    dppt = {} #demande pet
    dTT = udev #temps themique efficace (avec PP)

    for i in range(len(tab['nump'])):
        if scale == 'plt':
            idp = str(tab['nump'][i])
        elif scale == 'sh':
            idp = str(tab['nump'][i]) + '_' + str(tab['nsh'][i])
        elif scale == 'ax':
            idp = str(tab['nump'][i]) + '_' + str(tab['nsh'][i]) + '_' + str(tab['rank'][i])

        age = float(tab['age'][i])
        nump = int(tab['nump'][i])
        ordre = int(tab['ordre'][i])
        rank = int(tab['rank'][i])
        rankp = int(tab['rankp'][i])
        nsh = int(tab['nsh'][i])
        l = float(tab['l'][i])

        if tab['organ'][i] == 'Lf' and tab['statut'][i] == 'exp':
            pot = expansion(age + dTT[nump], ParamP[nump]['aF'], ParamP[nump]['delaiF']) - expansion(age,
                                                                                                     ParamP[nump]['aF'],
                                                                                                     ParamP[nump][
                                                                                                         'delaiF'])
            dl = pot * ls_ftswStress['WaterTreshExpSurf'][nump] * ls_NNIStress['NTreshExpSurf'][nump]
            dSpot = calc_surF(ParamP[nump], rank, rankp, ordre, l + dl) - calc_surF(ParamP[nump], rank, rankp, ordre,
                                                                                    l)  # m2, delta surf potentiel (sans limitation C mais avec stress hydrique)
            dMin = 10000. * dSpot / ParamP[nump]['SLAmin']  # delta masse min feuille
            IOxls.append_dic(dp, idp, dMin)
            IOxls.append_dic(dplf, idp, dMin)

        if tab['organ'][i] == 'Stp' and tab['statut'][i] == 'exp':
            pot = expansion(age + dTT[nump], ParamP[nump]['aS'], ParamP[nump]['delaiS']) - expansion(age,
                                                                                                     ParamP[nump]['aS'],
                                                                                                     ParamP[nump][
                                                                                                         'delaiS'])
            dl = pot * ls_ftswStress['WaterTreshExpSurf'][nump] * ls_NNIStress['NTreshExpSurf'][nump]
            dSpot = calc_surS(ParamP[nump], rank, rankp, ordre, l + dl) - calc_surS(ParamP[nump], rank, rankp, ordre,
                                                                                    l)  # m2, delta surf potentiel (sans limitation C mais avec stress hydrique)
            dMin = 10000. * dSpot / ParamP[nump]['SLAmin']  # delta masse min feuille
            IOxls.append_dic(dp, idp, dMin)
            IOxls.append_dic(dplf, idp, dMin)

        if tab['organ'][i] == 'In' and tab['statut'][i] == 'exp':
            pot = expansion(age + dTT[nump], ParamP[nump]['aE'], ParamP[nump]['delaiE']) - expansion(age,
                                                                                                     ParamP[nump]['aE'],
                                                                                                     ParamP[nump][
                                                                                                         'delaiE'])
            dl = pot * ls_ftswStress['WaterTreshExpSurf'][nump] * ls_NNIStress['NTreshExpSurf'][nump]
            dLpot = calc_Lent(ParamP[nump], rank, nsh, ordre,
                              dl)  # m, delta longueur potentiel (sans limitation C mais avec stress hydrique)
            dMin = dLpot / ParamP[nump]['SNLmin']  # delta masse min En
            IOxls.append_dic(dp, idp, dMin)
            IOxls.append_dic(dpin, idp, dMin)

        if tab['organ'][i] == 'Pet' and tab['statut'][i] == 'exp':
            pot = expansion(age + dTT[nump], ParamP[nump]['aP'], ParamP[nump]['delaiP']) - expansion(age,
                                                                                                     ParamP[nump]['aP'],
                                                                                                     ParamP[nump][
                                                                                                         'delaiP'])
            dl = pot * ls_ftswStress['WaterTreshExpSurf'][nump] * ls_NNIStress['NTreshExpSurf'][nump]
            dLpot = calc_Lpet(ParamP[nump], rank, rankp, ordre, dl)  # m
            dMin = dLpot / ParamP[nump]['SPLmin']  # delta masse min En
            IOxls.append_dic(dp, idp, dMin)
            IOxls.append_dic(dppt, idp, dMin)

    IOxls.sum_ls_dic(dp)
    IOxls.sum_ls_dic(dplf)
    IOxls.sum_ls_dic(dpin)
    IOxls.sum_ls_dic(dppt)

    return dp, dplf, dpin, dppt
    # pourrait decliner en demande pot sans et avec stress?

def Cremob(DemCp, R_DemandC_Shoot, MSPiv, frac_remob=0.1):
    """ remobilisation of C from the taproot to the shoot to ensure minimal growth """
    # frac_remob : fraction remobilisable du pivot par jour (a passer en parametre?)
    ratio_seuil = np.array(deepcopy(R_DemandC_Shoot))
    ratio_seuil[ratio_seuil > 1.] = 1.  # borne ratio demande a 1
    dem_non_couv = np.array(DemCp) * (1 - ratio_seuil)
    dem_non_couv_dispo = frac_remob * np.array(MSPiv) - dem_non_couv  # depend d'un fraction remobilisable du pivot par jour
    dem_non_couv_dispo[dem_non_couv_dispo < 0] = dem_non_couv[dem_non_couv_dispo < 0] + dem_non_couv_dispo[dem_non_couv_dispo < 0]  # borne remobilisation a poids du pivot
    remob = deepcopy(dem_non_couv_dispo)
    remob[dem_non_couv <= 0.] = 0.  # met a zero si couvert
    for i in range(len(remob)):
        remob[i] = min(remob[i], dem_non_couv[i])

    return remob
    # frac_remob dans ParamP -> fait


#calcul variables internes / intermediaire
def calcNB_NI(tab, nbplantes, seuilcountTige=0.5, seuilNItige=0.75):
    """ tab = tableau des apex actifs I et II =lsApex ; seuilNItige: faut au moins cette fraction du max pour etre compter en NI; seuilcountTige: pareil pour etre compter dans le nb tige 'significatives' """
    resall, resI, resNI, resNB = [], [], [], []
    for i in range(nbplantes):
        resall.append([0]);
        resI.append([0]);
        resNI.append([]);
        resNB.append([])

    for i in range(len(tab)):
        idp = int(tab[i][0])
        resall[idp].append(tab[i][2])
        if tab[i][3] == 1:
            resI[idp].append(tab[i][2])

    for i in range(nbplantes):
        resall[i] = max(resall[i])
        resI[i] = max(resI[i])

    for i in range(len(tab)):
        idp = int(tab[i][0])
        if float(tab[i][2]) > seuilcountTige * float(resall[idp]):
            resNB[idp].append(tab[i][2])

        if tab[i][3] == 1 and float(tab[i][2]) > seuilNItige * float(resI[idp]):
            resNI[idp].append(float(tab[i][2]))

    for i in range(nbplantes):
        resNB[i] = len(resNB[i])
        resNI[i] = np.mean(resNI[i])

    return resNB, resI

def cumul_lenIN(tab, tabL, I_I0profilInPlant_, deltaI_I0, nbI_I0):
    """ tab = tableau des apex actifs I et II =lsApex ;tabL = lsOrgans converti en dico"""
    # ajoute un id tige a tab lsApex en derniere colone
    for i in range(len(tab)): tab[i].append(str(tab[i][0]) + '_' + str(tab[i][1]))
    # liste d'id tiges unique
    tab = IOtable.t_list(tab)
    id_sh = list(set(tab[-1]))  # id tige dans derniere colone
    tab = IOtable.t_list(tab)

    # dico des I_I0 max par tige
    res = dict.fromkeys(id_sh, 0)
    for i in range(len(tab)):
        id = tab[i][-1]
        I_I0 = tab[i][4]
        if I_I0 > res[id]:
            res[id] = I_I0

    # dico des longueur cumul par tige
    resL = {}  # dict.fromkeys(id_sh, 0)
    for i in range(len(tabL['organ'])):
        if tabL['organ'][i] == 'In':
            id = str(tabL['nump'][i]) + '_' + str(tabL['nsh'][i])
            try:
                resL[id] += tabL['Long'][i] / 100.  # m
            except:  # si pas dans les cles, la cree
                resL[id] = tabL['Long'][i] / 100.  # pass

    # mise a jour des longueur de tige par classe d'eclairement
    for id in list(res.keys()):
        nump = list(map(int, id.split('_')))[0]#list(map(int, string.split(id, '_')))[0]
        I_I0 = res[id]
        classI_I0 = min(int(I_I0 / deltaI_I0), nbI_I0 - 1)  # pour gerer cas du I_I0=1.
        try:
            cumulL = resL[id]
        except:  # pas encore de In
            cumulL = 0.

        I_I0profilInPlant_[nump][classI_I0] += cumulL

    return I_I0profilInPlant_  # resL
    # bizarre -> donne meme longueur pour tous les axes??
    # faire dico tige pour toutes les tiges et croiser apres avec activeAxes?
    # a continuer avec un profil de longueur par I_I0


#initialiation plante
def MaturBud(delaiMaturBud, NIparent, delta=4):
    """ genere ecart de stade des B() (en phyllocrones) selon stade de developpement tige parente et delaiMaturBud
    - delta= borne min/max d'ecart par defaut +-4phyllo"""
    ecart = NIparent - delaiMaturBud
    if ecart >= 0:
        ecart = min(delta, ecart)
    else:
        ecart = max(-delta, ecart)

    return ecart
    # MaturBud(delaiMaturBud=12, NIparent=15, delta=2)
    # delta : passer ds fichier d'initialisation?


################
# Planter - initialisation scene

def planter_coordinates(type, cote, nbcote, orientRow='X'):
    """

    :param type:
    :param cote:
    :param nbcote:
    :param orientRow:
    :return: carto, list of nd.arrays for the 3D coordiantes of each plant in the scene
    """
    distplantes = cote / nbcote  # 1. #cm

    # pour grand rhizotron
    # yyy = [-12.25, 4.4, 21.05]
    # xxx = [-10.75, -8.35, -5.95, -3.55, -1.15, 1.25, 3.65, 6.05, 8.45, 10.85]

    if type == 'row4' or type=='row4_sp1' or type=='row4_sp2':  # pour carre 4 rangs heterogenes
        Param_, carto = row4([1, 2], Lrow=cote, nbprow=nbcote, orientRow=orientRow)
    elif type == 'row4_Nsp' :  # pour carre 4 rangs heterogenes
        Param_, carto = row4_Nsp([1, 2], Lrow=cote, nbprow=nbcote, orientRow=orientRow)
    elif type == 'random8' or type=="random9" or type=="random8_4" or type=="random10":
        # pour tirage random
        carto = random_planter(nbcote * nbcote, cote, cote)
    elif type=="ilot7":
        carto = Ilot7(distplantes)
    elif type=="damier8" or type=="damier16" or type=="damier9" or type=="homogeneous" or type=="damier8_4" or type=="damier10" :
        # pour carre distance homogene
        carto = regular_square(nbcote, distplantes)
    elif type=="damier8_sp1" or type=="damier16_sp1" or type=="damier8_sp2" or type=="damier16_sp2":
        # pour carre distance homogene
        carto = regular_square(nbcote, distplantes)
        #doit etre reduit apres avec reduce
    else:
        print("unknown type for planter coordinates")
        carto = []

    return carto



def Ilot7(distplantes):
    ## coord pour ilot 7 plantes
    carto = [np.array([0.,0.,0.]), np.array([distplantes,0.,0.]),np.array([-distplantes,0.,0.]),np.array([0.5*distplantes,0.866*distplantes,0.]),np.array([-0.5*distplantes,0.866*distplantes,0.]), np.array([0.5*distplantes,-0.866*distplantes,0.]), np.array([-0.5*distplantes,-0.866*distplantes,0.])]#, array([-10.,0.,0.]), array([0.,7.,0.])] #liste des localisations (1pt par plante) -> a lire en fichier #LF - cos (pi/3) = 0.5   sin (pi/3) = 0.866
    return carto


def regular_square(nbcote, distplantes):
    "regular square planter - used with damier8 / homogeneous "
    yyy = [distplantes / 2.]
    for i in range(1, nbcote): yyy.append(yyy[-1] + distplantes)
    xxx = yyy

    carto = []
    for i in range(len(xxx)):
        for j in range(len(yyy)):
            carto.append(np.array([xxx[i], yyy[j], 0.]))

    return carto


def planter_order_ParamP(ls_g, type, nbcote, opt, shuffle=0):
    """

    :param ls_g:
    :param type:
    :param nbcote:
    :param opt:
    :param speby_row:
    :param shuffle:
    :return: ParamP, list of plant parameter dictionnaries for each plant in the scene
    """

    if type == 'homogeneous':  # cas d'un couvert monospe homogene: prend premier parametrage
        ParamP = [ls_g[0]] * nbcote * nbcote
    elif type == 'damier8' or type == 'random8' or type == 'damier8_sp1' or type == 'damier8_sp2':  # damier binaire 64 plantes
        if nbcote == 8:
            ParamP = damier8(ls_g[0], ls_g[1], opt=opt)
        else:
            # if opt_verbose==1:
            print('Error! :' + type + ' option is for a 64 plant design')
    elif type == 'damier16' or type == 'random16' or type == 'damier16_sp1' or type == 'damier16_sp2':  # damier binaire 256 plantes
        if nbcote == 16:
            ParamP = damier16(ls_g[0], ls_g[1], opt=opt)
        else:
            # if opt_verbose==1:
            print('Error! :' + type + ' option is for a 256 plant design')
    elif type == 'damier9' or type == 'random9':  # damier 3sp 81 plantes
        if nbcote == 9:
            ParamP = damier9_3sp(ls_g, opt=opt)
        else:
            # if opt_verbose==1:
            print('Error! :' + type + ' option is for a 81 plant design')
    elif type == 'damier10' or type == 'random10':  # damier 3sp 81 plantes
        if nbcote == 10:
            ParamP = damier10_5sp(ls_g, opt=opt)
        else:
            # if opt_verbose==1:
            print('Error! :' + type + ' option is for a 100 plant design')
    elif type == 'damier8_4' or type == 'random8_4':  # damier 3sp 81 plantes
        if nbcote == 8:
            ParamP = damier8_4sp(ls_g, opt=opt)
        else:
            # if opt_verbose==1:
            print('Error! :' + type + ' option is for a 64 plant design')
    elif type == 'row4' or type == 'row4_sp1' or type == 'row4_sp2':  # 4 rangs - 2 spe
        ParamP, cart_ = row4(ls_g, nbprow=nbcote, opt=opt)
    elif type == 'row4_Nsp' :  # 4 rangs - N spe
        ParamP, cart_  = row4_Nsp(ls_g, nbprow=nbcote, opt=opt)
    else:
        # defautl= force nb plante comme nbcote pour premier id des parametres
        ParamP = [ls_g[0]] * nbcote

    if shuffle==1:
        # shuffle list in place
        np.random.shuffle(ParamP)

    return ParamP





def random_planter(nbplt, cotex, cotey):
    """ random plant position with scene defined by cotex and cotey"""
    carto=[]
    for i in range(nbplt):
        carto.append(np.array([np.random.uniform(0., cotex), np.random.uniform(0., cotey), 0.]))

    return carto


def damier8(p, vois, opt=4):
    # cree liste de plante ordonnee pour un melange binaire homogene de 64 plantes avec differentes options de proportions
    if opt == 4:  # 50/50
        motif = [p, vois, p, vois, p, vois, p, vois]
    elif opt == 0:  # 0/100
        motif = [vois, vois, vois, vois, vois, vois, vois, vois]
    elif opt == 8:  # 100/0
        motif = [p, p, p, p, p, p, p, p]
    elif opt == 2:  # 25/75
        motif = [p, vois, vois, vois, p, vois, vois, vois]
    elif opt == 6:  # 75/25
        motif = [vois, p, p, p, vois, p, p, p]
    elif opt == 1:  # 1/8
        motif = [p, vois, vois, vois, vois, vois, vois, vois]
    elif opt == 7:  # 7/8
        motif = [vois, p, p, p, p, p, p, p]
    elif opt == 3:  # 5/8
        motif = [vois, p, p, vois, p, vois, p, p]
    elif opt == 5:  # 3/8
        motif = [p, vois, vois, p, vois, p, vois,vois]

    res = []
    for i in range(8):
        res = res + motif[i:8] + motif[0:i]

    return res


def damier16(p, vois, opt=4):
    # cree un melange binaire homogene de 256 plantes avec differentes options de proportions
    if opt == 4:  # 50/50
        motif = [p, vois, p, vois, p, vois, p, vois]+[p, vois, p, vois, p, vois, p, vois]
    elif opt == 0:  # 0/100
        motif = [vois, vois, vois, vois, vois, vois, vois, vois]+[vois, vois, vois, vois, vois, vois, vois, vois]
    elif opt == 8:  # 100/0
        motif = [p, p, p, p, p, p, p, p]+[p, p, p, p, p, p, p, p]
    elif opt == 2:  # 25/75
        motif = [p, vois, vois, vois, p, vois, vois, vois]+[p, vois, vois, vois, p, vois, vois, vois]
    elif opt == 6:  # 75/25
        motif = [vois, p, p, p, vois, p, p, p]+[vois, p, p, p, vois, p, p, p]
    elif opt == 1:  # 1/8
        motif = [p, vois, vois, vois, vois, vois, vois, vois]+[p, vois, vois, vois, vois, vois, vois, vois]
    elif opt == 7:  # 7/8
        motif = [vois, p, p, p, p, p, p, p]+[vois, p, p, p, p, p, p, p]
    elif opt == 3:  # 5/8
        motif = [vois, p, p, vois, p, vois, p, p]+[vois, p, p, vois, p, vois, p, p]
    elif opt == 5:  # 3/8
        motif = [p, vois, vois, p, vois, p, vois,vois]+[p, vois, vois, p, vois, p, vois,vois]

    res = []
    for i in range(16):
        res = res + motif[i:16] + motif[0:i]

    return res


def damier9_3sp(ls_g, opt=4):
    # cree liste de plante de 3 sp. ordonnee pour un melange binaire homogene de 81 plantes avec differentes options de proportions
    if opt == 4:  # 50/50
        motif = [ls_g[0], ls_g[1], ls_g[2], ls_g[0], ls_g[1], ls_g[2], ls_g[0], ls_g[1], ls_g[2]]

    res = []
    for i in range(9):
        res = res + motif[i:9] + motif[0:i]

    return res

def damier10_5sp(ls_g, opt=4):
    # cree liste de plante de 5 sp. ordonnee pour un melange binaire homogene de 100 plantes avec differentes options de proportions
    if opt == 4:  # 50/50
        motif = [ls_g[0], ls_g[1], ls_g[2], ls_g[3], ls_g[4], ls_g[0], ls_g[1], ls_g[2], ls_g[3], ls_g[4]]

    res = []
    for i in range(10):
        res = res + motif[i:10] + motif[0:i]

    return res

def damier8_4sp(ls_g, opt=4):
    # cree liste de plante de 4 sp. ordonnee pour un melange 4 Sp homogene de 64 plantes avec differentes options de proportions
    if opt == 4:  # 50/50
        motif = [ls_g[0], ls_g[1], ls_g[0], ls_g[1], ls_g[0], ls_g[1], ls_g[0], ls_g[1]]
        motif2 = [ls_g[2], ls_g[3], ls_g[2], ls_g[3], ls_g[2], ls_g[3], ls_g[2], ls_g[3]]

    res = []
    for i in range(int(8/2)):
        res = res + motif + motif2

    return res


def row4(ls_g, Lrow=50., nbprow=125,  opt=0, orientRow='X'):
    """ cree un melange 50/50 alterne ou pur sur 4 rangs distance interow chanmp"""
    if opt == 2:  # 50/50
        motif = [ls_g[0], ls_g[1], ls_g[0], ls_g[1]]
    elif opt == 0:  # 0/100
        motif = [ls_g[1], ls_g[1], ls_g[1], ls_g[1]]
    elif opt == 4:  # 100/0
        motif = [ls_g[0], ls_g[0], ls_g[0], ls_g[0]]

    res = []
    for i in range(nbprow):
        res = res + motif

    inter = Lrow / 4.
    onrow = Lrow / nbprow
    xxx = np.arange(0., Lrow, onrow)+onrow/2.
    yyy = [0. * inter + inter/2.] + [1. * inter + inter/2.] + [2. * inter + inter/2.] + [3. * inter + inter/2.]  # 4 rangs
    carto = []
    for i in range(len(xxx)):
        for j in range(len(yyy)):
            if orientRow=='X': #row parralel to X axis
                carto.append(np.array([xxx[i], yyy[j], 0.]))  # +origin
            elif orientRow=='Y': #row parralel to Y axis
                carto.append(np.array([yyy[j], xxx[i], 0.]))  # +origin
            else:
                print("unknown orientRow")

    return res ,carto
    #pourrait renvoyer carto aussi ds homogeneous et damier8....
    #res ,carto=row4(1, 2, Lrow=50., nbprow=125,  opt=0)
    # dans un fichier d'initialiation?
    #prevoir nbprow different par esp... et melange on row...



def row4_Nsp(ls_g, Lrow=50., nbprow=125,  opt=0, orientRow='X'):
    """ """
    if opt==0: # equipropotion on all rows
        speby_row = []
    elif opt==23: #5 spe ; two first on first row, 3 others on second row
        speby_row = [[0,1], [2,3,4], [0,1], [2,3,4]]
    elif opt==21: #3 spe ; two first on first row, 1 others on second row
        speby_row = [[0,1], [2], [0,1], [2]]
    else: #all other cases
        speby_row = []

    #!! pas prop equiprobable entre espece si nb esp different par rang!
    #

    if speby_row == []:
        ls_g1 = ls_g
        ls_g2 = ls_g
        ls_g3 = ls_g
        ls_g4 = ls_g
    else : #list of species by row speby_row provided as 4 lists of IDs in ls_g
        ls_g1 = [ls_g[i] for i in speby_row[0]]
        ls_g2 = [ls_g[i] for i in speby_row[1]]
        ls_g3 = [ls_g[i] for i in speby_row[2]]
        ls_g4 = [ls_g[i] for i in speby_row[3]]

    row1, row2, row3, row4 = [], [], [], []
    for i in range(nbprow):
        row1 += ls_g1
        row2 += ls_g2
        row3 += ls_g3
        row4 += ls_g4

    print('row1', ls_g1, row1)

    #equi prop / opt pas utilise
    res = []
    for i in range(nbprow):
        motif = [row1[i], row2[i], row3[i], row4[i]]
        res = res + motif


    inter = Lrow / 4.
    onrow = Lrow / nbprow
    xxx = np.arange(0., Lrow, onrow)+onrow/2.
    yyy = [0. * inter + inter/2.] + [1. * inter + inter/2.] + [2. * inter + inter/2.] + [3. * inter + inter/2.]  # 4 rangs
    carto = []
    for i in range(len(xxx)):
        for j in range(len(yyy)):
            if orientRow=='X': #row parralel to X axis
                carto.append(np.array([xxx[i], yyy[j], 0.]))  # +origin
            elif orientRow=='Y': #row parralel to Y axis
                carto.append(np.array([yyy[j], xxx[i], 0.]))  # +origin
            else:
                print("unknown orientRow")


    return res ,carto






def reduce_ParamP(ParamP, nom):
    """ pour extraire une espece du ParamP """
    ls_name = IOxls.get_lsparami(ParamP, 'name')
    newParamP, lsid_reduce = [], []
    for i in range(len(ls_name)):
        if ParamP[i]['name']==nom:
            newParamP.append(ParamP[i])
            lsid_reduce.append(i)

    return newParamP, lsid_reduce

def reduce_carto(carto, lsid_reduce):
    """ pour extraire une espece de carto a partir des id_reduce du ParamP """
    newcarto = []
    for i in range(len(carto)):
        if i in lsid_reduce:
            newcarto.append(carto[i])

    return newcarto




def ls_idvois_ordre1(n, cote, nblignes):
    """ pour une plante n, dans un dispocitif regulier arrange en colonnes croissantes de cote indiv"""
    nbindiv = cote * nblignes
    ls_defaut = [n - (cote + 1), n - cote, n - (cote - 1), n - 1, n + 1, n + (cote - 1), n + cote, n + (cote + 1)]

    if n % cote == 0:  # bord haut
        ls_defaut[0] = ls_defaut[0] + cote
        ls_defaut[3] = ls_defaut[3] + cote
        ls_defaut[5] = ls_defaut[5] + cote

    if (n + 1) % cote == 0:  # bord bas
        ls_defaut[2] = ls_defaut[2] - cote
        ls_defaut[4] = ls_defaut[4] - cote
        ls_defaut[7] = ls_defaut[7] - cote

    for i in range(len(ls_defaut)):
        if ls_defaut[i] < 0:  # bord gauche
            ls_defaut[i] = ls_defaut[i] + nbindiv

        if ls_defaut[i] >= nbindiv:  # bord droit
            ls_defaut[i] = ls_defaut[i] - nbindiv

    return ls_defaut
    # pour les voisin de 2e ordre = suffit de faire les voisin de 1er ordre de chasue voisin + unique!
    # a faire dans R #11 et 23 plante / 0 marche
    # print("id vois", ls_idvois_ordre1(4,6,4)) #OK!
#pour pour nump python, pas pour id R!








def updateLargProfile(Lmax, Largmax, profilLeafI_Rlen, profilLeafI_Rlarg):
    """ A function to modify The relative LArg profile given Largmax for the Longest leaf - old parametrization way"""
    # change the intercept of the second line b2
    if Largmax > 0.: #negative values are not condidered for update of the intercept
        ratioM = Largmax/Lmax
        peak = (profilLeafI_Rlen[3] - profilLeafI_Rlen[1]) / (profilLeafI_Rlen[0] - profilLeafI_Rlen[2])
        newb2 = ratioM - (peak * profilLeafI_Rlarg[2])#same slope but pass by Largmax at peak
        newb1 = newb2*(profilLeafI_Rlarg[1]/profilLeafI_Rlarg[3])#same ratio of intercept
        return [profilLeafI_Rlarg[0], newb1, profilLeafI_Rlarg[2], newb2]
    else:
        return profilLeafI_Rlarg #unchanged


def update_shoot_params(ParamP, rankmax=51):
    """ add readable rank profiles in paramP """
    for nump in range(len(ParamP)):

        if int(ParamP[nump]['type']) == 1 or int(ParamP[nump]['type']) == 2:  # feuille legumineuse
            cor_lF = np.sqrt(ParamP[nump]['leafshape'] / 0.5)  # pour afficher feuille avec surface reelle et corriger effet losange
        elif int(ParamP[nump]['type']) == 3:  # graminee
            cor_lF = ParamP[nump]['leafshape']  # pour afficher feuille avec surface reelle et corriger effet recangle (pas sqrt car appliquer que a largeur

        cor_lstp = np.sqrt(ParamP[nump]['stipshape'] / 0.5)  # pour afficher feuille avec surface reelle et corriger effet losange

        ParamP[nump]['profilLeafI_l'] = []
        ParamP[nump]['profilLeafI_larg'] = []
        ParamP[nump]['profilNodeI_l'] = []
        ParamP[nump]['profilPetI_l'] = []
        ParamP[nump]['profilStipI_l'] = []
        ParamP[nump]['profilStipI_larg'] = []
        ParamP[nump]['profilPetI_l'] = []
        ParamP[nump]['profilLeafI_nfol'] = []
        if int(ParamP[nump]['type']) == 1 or int(ParamP[nump]['type']) == 2:  # legumineuse
            ParamP[nump]['k_teta_distf'] = riri.disttetaf(abs(ParamP[nump]['gammaFeuil']), ParamP[nump]['gammaFeuilSD'])  # proportion de feuille par    #    classe d'incli pour calcul des k_teta
        elif int(ParamP[nump]['type']) == 3:  # graminee avec courbure
            nfol = 8
            courbure = -5
            angle_moyen = 0.5 * (ParamP[nump]['gammaFeuil'] + (nfol * courbure) + ParamP[nump]['gammaFeuil'])  # marche pour nb fixe de rectangles (??coupe)
            ParamP[nump]['k_teta_distf'] = riri.disttetaf(abs(angle_moyen), abs(2 * courbure))
            # passer courbure en parametre?? dans 'gammaFeuilSD'??

        ParamP[nump]['profilLeafI_Rlen'] = [ParamP[nump]['profilLeafI_Rlens1'], ParamP[nump]['profilLeafI_Rleni1'], ParamP[nump]['profilLeafI_Rlens2'] , ParamP[nump]['profilLeafI_Rleni2']]
        ParamP[nump]['profilLeafI_Rlarg'] = [ParamP[nump]['profilLeafI_Rlargs1'] , ParamP[nump]['profilLeafI_Rlargi1'], ParamP[nump]['profilLeafI_Rlargs2'] , ParamP[nump]['profilLeafI_Rlargi2']]
        ParamP[nump]['profilLeafI_Rlarg'] = updateLargProfile(ParamP[nump]['Lfeuille'], ParamP[nump]['Largfeuille'], ParamP[nump]['profilLeafI_Rlen'], ParamP[nump]['profilLeafI_Rlarg'])# to consider Largmax as a direct input forcing the relative profile
        #print('profilLeafI_Rlarg', ParamP[nump]['profilLeafI_Rlarg'])
        for rank in range(1, rankmax):  # !limite a 50 noeuds! (rankmax)
            Norml_leaf = min(ParamP[nump]['profilLeafI_Rlens1'] * rank + ParamP[nump]['profilLeafI_Rleni1'], ParamP[nump]['profilLeafI_Rlens2'] * rank + ParamP[nump]['profilLeafI_Rleni2'])
            Normlarg_leaf = max(ParamP[nump]['profilLeafI_Rlarg'][0] * rank + ParamP[nump]['profilLeafI_Rlarg'][1], ParamP[nump]['profilLeafI_Rlarg'][2] * rank + ParamP[nump]['profilLeafI_Rlarg'][3])
            Norml_In = min(ParamP[nump]['profilNodeIs1'] * rank + ParamP[nump]['profilNodeIi1'],ParamP[nump]['profilNodeIs2'] * rank + ParamP[nump]['profilNodeIi2'])
            Norm_pet = min(ParamP[nump]['profilPetIs1'] * rank + ParamP[nump]['profilPetIi1'],ParamP[nump]['profilPetIs2'] * rank + ParamP[nump]['profilPetIi2'])
            Norml_Stp = min(ParamP[nump]['profilStpI_ls1'] * rank + ParamP[nump]['profilStpI_li1'],ParamP[nump]['profilStpI_ls2'] * rank + ParamP[nump]['profilStpI_li2'])
            Normlarg_Stp = min(ParamP[nump]['profilStpI_Rlargs1'] * rank + ParamP[nump]['profilStpI_Rlargi1'],ParamP[nump]['profilStpI_Rlargs2'] * rank + ParamP[nump]['profilStpI_Rlargi2'])
            Normnfol = min(ParamP[nump]['profilLeafI_Rnfols'] * rank + ParamP[nump]['profilLeafI_Rnfoli'],1.)  # nfol est le maximum number of folioles

            if int(ParamP[nump]['type']) == 1 or int(ParamP[nump]['type']) == 2:  # feuille legumineuse
                ParamP[nump]['profilLeafI_l'].append(max(0.001, Norml_leaf * cor_lF * ParamP[nump]['Lfeuille']))
                ParamP[nump]['profilLeafI_larg'].append(max(0.001, Normlarg_leaf * Norml_leaf * cor_lF * ParamP[nump]['Lfeuille']))
            elif int(ParamP[nump]['type']) == 3:  # graminee (applique slmt largeur)
                ParamP[nump]['profilLeafI_l'].append(max(0.001, Norml_leaf * 1. * ParamP[nump]['Lfeuille']))
                ParamP[nump]['profilLeafI_larg'].append(max(0.001, Normlarg_leaf * Norml_leaf * cor_lF * ParamP[nump]['Lfeuille']))

            ParamP[nump]['profilNodeI_l'].append(max(0.001, Norml_In * ParamP[nump]['Len']))
            ParamP[nump]['profilPetI_l'].append(max(0.001, Norm_pet * ParamP[nump]['Lpet']))
            ParamP[nump]['profilStipI_l'].append(max(0.001, Norml_Stp * cor_lstp * ParamP[nump]['Lstip']))
            ParamP[nump]['profilStipI_larg'].append(max(0.001, Normlarg_Stp * Norml_Stp * cor_lstp * ParamP[nump]['Lstip']))
            ParamP[nump]['profilLeafI_nfol'].append(int(max(1, Normnfol * ParamP[nump]['nfol'])))  # max 1 pour interdire les feuilles sans folioles.

        #MAJ param residues (old format)
        ParamP[nump]['CC'] = [ParamP[nump]['CClf'], ParamP[nump]['CCst'],ParamP[nump]['CCr'],ParamP[nump]['CCpiv']]
        ParamP[nump]['WC'] = [ParamP[nump]['WClf'], ParamP[nump]['WCst'], ParamP[nump]['WCr'], ParamP[nump]['WCpiv']]
        ParamP[nump]['Nmires'] = [ParamP[nump]['Nmireslf'], ParamP[nump]['Nmiresst'], ParamP[nump]['Nmiresr'], ParamP[nump]['Nmirespiv']]

    return ParamP




#old - non utilise
# a retirer
def calcLeafStemRatio(ParamP, tab, lsapexI):
    """ calcul de rapport feuille/tige base sur les tigeI actives"""
    # etablit une liste de tige avec A() actif (lsa)
    lsa = []
    for i in range(len(lsapexI)):
        newk = str(lsapexI[i][0]) + '_' + str(lsapexI[i][1])  # id = 'nump_nsh'
        lsa.append(newk)

    # recupere les masses mini des feuilles et tiges des lsa
    dp, dp2 = {}, {}
    for i in range(len(tab['nump'])):
        idp = str(tab['nump'][i])
        idt = str(tab['nump'][i]) + '_' + str(tab['nsh'][i])

        age = float(tab['age'][i])
        nump = int(tab['nump'][i])
        ordre = int(tab['ordre'][i])
        rank = int(tab['rank'][i])
        rankp = int(tab['rankp'][i])
        l = float(tab['l'][i])

        if tab['organ'][i] == 'Lf' and idt in lsa:
            surf = calc_surF(ParamP[nump], rank, rankp, ordre, l)  # m2
            MLf = 10000. * surf / ParamP[nump]['SLAmin']  # masse min feuille
            IOxls.append_dic(dp, idp, MLf)

        if tab['organ'][i] == 'In' and idt in lsa:
            cor_ordre = ParamP[nump]['ratioII'] if ordre == 2 else 1.
            rank = min(rank, len(ParamP[nump]['profilNodeI_l']) - 1)  # au cas ou profil trop long
            Long = l * ParamP[nump]['profilNodeI_l'][
                rank] * cor_ordre / 100.  # m #delta de longueur potentiel (sans limitation C)
            MIn = Long / ParamP[nump]['SNLmin']
            IOxls.append_dic(dp2, idp, MIn)

    # fait somme par plante de Lf et In, puis ratio
    for k in list(dp.keys()):
        leafM = sum(dp[k])
        try:
            inM = sum(dp2[k]) + 0.00000000001
        except:
            inM = 0.00000000001

        dp[k] = leafM / inM

    return dp

def PhylloPot_Grass(rgeq, phylloI, k=0.3):
    #calcule changement de phyllochrone selon le rang
    if rgeq>15:
        phyllo = phylloI
    elif rgeq<1:
        phyllo = monomoleculaire(1., Amax=phylloI, k=k)
    else:
        phyllo = monomoleculaire(rgeq, Amax=phylloI, k=k)

    return phyllo
