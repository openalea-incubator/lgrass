import sys
import os

try:
    import legume
    path_ = os.path.dirname(os.path.abspath(legume.__file__))#local absolute path of L-egume
    #path_ = r'C:\devel\l-egume\legume'#r'C:\devel\grassland'#r'H:\devel\grassland\grassland\luzerne' #r'C:\devel\grassland'#
except:
    path_ = r'C:\devel\l-egume\legume'#r'C:\devel\grassland'#r'H:\devel\grassland\grassland\luzerne' #r'C:\devel\grassland'#

path_leg = os.path.join(path_, 'input')#r'C:\devel\l-egume\l-egume\input'#r'C:\devel\grassland'#r'H:\devel\grassland\grassland\L-gume' #r'C:\devel\grassland'
path_out = os.path.join(path_, 'output')#r'C:\devel\grassland'#r'H:\devel\grassland\grassland\L-gume' #r'C:\devel\grassland'

sys.path.insert(0, path_)
sys.path.insert(0, path_leg)

import IOxls
import ShootMorpho as sh
#from scipy import array
import numpy as np

#meteo
#meteo_path = os.path.join(path_leg,'meteo_exemple.xls')#'meteo_exemple_debugL_gl.xls')##r'H:\devel\grassland\grassland\L-gume\meteo_exemple2.xls'
#ongletM = 'Lusignan30'#'DigitLuz10'#'Lusignan30'#'Lusignan302ans'#'DivLeg15'#'morpholeg15'#'combileg15'#'combileg16'#'Avignon30'#'exemple'#'morpholeg15'#'testJLD'#'competiluz'#
#meteo = IOxls.read_met_file(meteo_path, ongletM)

## lecture fichier management auto
#mn_path = os.path.join(path_leg,'management_auto_exemple.xls')#'management_exemple3_debugL_gl.xls')#r'H:\devel\grassland\grassland\L-gume\management_exemple.xls'
#ongletMn = 'exemple'
#mng_auto = IOxls.read_plant_param(mn_path, ongletMn)

## plante
#path_plante = os.path.join(path_leg,'Parametres_plante_exemple.xls')#'Parametres_plante_v5cLucas.xls')#'Parametres_plante_v18.xls')#'Parametres_plante_v9Lucas_debugL.xls')#r'H:\devel\grassland\grassland\L-gume\Parametres_plante_v5cLucas.xls'
#ongletP = 'Fix2'#'Orca'#'Fet'#'giga'#'solnu'#'Fix1'#'Fix'#'timbale'#'formica'#'canto'##'alfalfa'#'geno_test'#'G3'#'C1'#'8_2'#'kayanne'#'leo'#'G1'#'timb
#g4 = IOxls.read_plant_param(path_plante, ongletP)


## DOYdeb pour les calcul de TT
#DOYdeb = 60
#opt_optT = 0 #option de calcul du cumul de temperature (0=betaD; 1=betaH; 2=lineaireD) # a bien reporter!

def Build_mng_auto(meteo, mng, path_plante, ongletP, DOYdeb, opt_optT):
    """ recalul en entree le management a partir de la meteo et du fichier palnte """
    #prend valeur par defaut d'ongletP pour les calcul de TT
    paramP = IOxls.read_plant_param(path_plante, ongletP)

    mng_auto = {'Date': meteo['Date'], 'year':meteo['year'], 'month':meteo['month'], 'day':meteo['day'], 'DOY':meteo['DOY']}

    vdTT, vTTcum, TTcum = [], [], 0.
    for DOY  in mng_auto['DOY']:
        meteo_j = IOxls.extract_dataframe(meteo, ['TmoyDay','RG','Et0','Precip','Tmin','Tmax','Tsol'], 'DOY', val=DOY)
        vT, vTsol = sh.Calc_Daily_vT(meteo_j, opt_optT) # daily temperature vector
        valdTT = sh.dTT(vT, [paramP['Tdev'], paramP['Tmin'], paramP['Tmax'], paramP['q']], optT=opt_optT)
        vdTT.append(valdTT)
        if DOY>=DOYdeb:
            TTcum += valdTT

        vTTcum.append(TTcum)

    mng_auto['TT'] = vTTcum


    #en fonction d'un vecteur de valeurs
    tt = np.array(mng_auto['TT'])

    #calul des dates de coupes
    mng_auto['Coupe'] = [0]*len(mng_auto['TT'])
    mng_auto['Hcut'] = [0]*len(mng_auto['TT'])

    if int(mng['opt_auto_cut']) == 1:  # vecteur en TT
        for i in range(len(mng['auto_cut_date'])):
            #i = 0
            v1 = tt >= mng['auto_cut_date'][i]
            v1 = v1.tolist()
            try:
                id = v1.index(True)
            except:
                id=-1

            if id>0:
                #print(id)
                mng_auto['Coupe'][id] = 1
                mng_auto['Hcut'][id] = mng['auto_cut_height'][i]

    if int(mng['opt_auto_cut']) == 2:  # regulier en TT
        IntCut = mng['auto_cut_date']
        nbcut = max(tt) / IntCut
        NextCut = mng['auto_cut_date']
        for i in range(nbcut):
            #i = 0
            v1 = tt >= NextCut
            v1 = v1.tolist()
            try:
                id = v1.index(True)
            except:
                id=-1

            if id>0:
                #print(id)
                mng_auto['Coupe'][id] = 1
                mng_auto['Hcut'][id] = mng['auto_cut_height']

            NextCut += mng['auto_cut_date']


    #calul des apports d'N
    mng_auto['FertNO3'] = [0]*len(mng_auto['TT'])
    mng_auto['FertNH4'] = [0]*len(mng_auto['TT'])

    #ammonium
    if int(mng['opt_auto_NH4'])==1: #vecteur en TT
        for i in range(len(mng['auto_NH4_date'])):
            #i = 0
            v1 = tt >= mng['auto_NH4_date'][i]
            v1 = v1.tolist()
            try:
                id = v1.index(True)
            except:
                id=-1

            if id>0:
                #print(id)
                mng_auto['FertNH4'][id] = mng['auto_NH4_ammount'][i]

    if int(mng['opt_auto_NH4'])==2: #regulier en TT
        IntCut = mng['auto_NH4_date']
        nbcut = max(tt) / IntCut
        NextCut = mng['auto_NH4_date']
        for i in range(nbcut):
            #i = 0
            v1 = tt >= NextCut
            v1 = v1.tolist()
            try:
                id = v1.index(True)
            except:
                id=-1

            if id>0:
                #print(id)
                mng_auto['FertNH4'][id] = mng['auto_NH4_ammount']

    #calul des apports d'eau
    mng_auto['Irrig'] = [0]*len(mng_auto['TT'])
    # pas d'irrig
    # a completer au besoin - to do!

    #fonction avec return du mng_auto
    return mng_auto


# mng = Build_mng_auto(meteo, mng_auto, path_plante, ongletP, DOYdeb=60, opt_optT=0)
# mng.keys()
# mng['FertNH4']

#bien s'assurer des appels
#ajouter une option generale qui dit quel type de fichier management est fourni en entree et si doit recalculer le management


#tres sommaire encore!!
# ajouter option 2 reguliere en TT (au lieu de vecteur!)
# pourrait aussi passer vecteur en fonction des DOY
# declenchement Irrig / N en fonction de seuils?
# pas de coupe en fonction de seuil?

#to be added
# a function that modify the mng dict for the next day depending on a condition (e.g. height / water status...) to avoid calculaltion of all in advance




