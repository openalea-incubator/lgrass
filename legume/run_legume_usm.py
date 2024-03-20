
#import the modules necessary to initiate the L-systems
import openalea.lpy as lpy
#from openalea.lpy import *

import os
import sys

try:
    import legume
    path_ = os.path.dirname(os.path.abspath(legume.__file__))#local absolute path of L-egume
except:
    path_ = r'C:\devel\l-egume\legume'#r'C:\devel\grassland'

#print(('path', path_))

sys.path.insert(0, path_)
import IOxls
import IOtable
import pandas as pd
import getopt
#import zipfile




#def lsystemInputOutput_usm(path_, fxls_usm, i=0, foldin = 'input', ongletBatch = 'exemple'):
def lsystemInputOutput_usm(fxls_usm, foldin = 'input', ongletBatch = 'exemple', i=0, path_OUT='output', update_usm_parameters=None):
    """" cree et update l-system en fonction du fichier usm """

    # lecture de la liste des usm
    # path_ = r'H:\devel\grassland\grassland\L-gume'
    usm_path = os.path.join(foldin, fxls_usm)#(path_, foldin, fxls_usm)
    usms = IOxls.xlrd.open_workbook(usm_path)
    ls_usms = IOtable.conv_dataframe(IOxls.get_xls_col(usms.sheet_by_name(ongletBatch)))
    #foldin = pour cas ou fichier d'usm dans un sous dossier different de input / tous les autres sont dans input

    if update_usm_parameters is not None:
        for key, value in update_usm_parameters.items():
            ls_usms[key][i] = value

    #nom fichier en dur (pas en entree de la fonction) + onglet determine par geno
    fscenar = 'liste_scenarios.xls' #'liste_scenarios_exemple.xls'
    fsd = 'exemple_sd.xls' #nom mis a jour mais pas table variance_geno
    fsdx = 'exemple_corr_matrix.xls'
    fopt = 'mod_susm.xls'## fichier option
    fsta = 'stations_exemple.xls' ##fichier station
    ongletSta = 'Lusignan'  # 'exemple'

    #options
    path_opt = os.path.join(foldin, fopt)#(path_, 'input', fopt)
    dic_opt = IOxls.read_plant_param(path_opt, "options")
    #print(dic_opt)

    testsim = {} #dico sortie avec un nom d'usm
    name = str(int(ls_usms['ID_usm'][i])) + '_' + str(ls_usms['l_system'][i])[0:-4]
    seednb = int(ls_usms['seed'][i])
    #names.append(name)
    path_ = os.path.dirname(os.path.abspath(legume.__file__))  # local absolute path of L-egume
    path_lsys = os.path.join(path_, str(ls_usms['l_system'][i]))
    testsim[name] = lpy.Lsystem(path_lsys)  # objet l-system

    # testsim[name].ongletM = str(ls_usms['ongletM'][i])
    meteo_path_ = os.path.join(foldin, str(ls_usms['meteo'][i]))#(path_, 'input', str(ls_usms['meteo'][i]))
    ongletM_ = str(ls_usms['ongletM'][i])
    testsim[name].meteo = IOxls.read_met_file(meteo_path_, ongletM_)

    # testsim[name].ongletMn = str(ls_usms['ongletMn'][i])
    mn_path_ = os.path.join(foldin, str(ls_usms['mng'][i]))#(path_, 'input', str(ls_usms['mng'][i]))
    ongletMn_ = str(ls_usms['ongletMn'][i])
    testsim[name].mng = IOxls.read_met_file(mn_path_, ongletMn_)

    ini_path_ = os.path.join(foldin, str(ls_usms['inis'][i]))#(path_, 'input', str(ls_usms['inis'][i]))
    ongletIni_ = str(ls_usms['ongletIn'][i])
    testsim[name].inis = IOxls.read_plant_param(ini_path_, ongletIni_)

    # testsim[name].ongletP = str(ls_usms['ongletP'][i])
    path_plante = os.path.join(foldin, str(ls_usms['plante'][i]))#(path_, 'input', str(ls_usms['plante'][i]))
    testsim[name].path_plante = path_plante
    path_lsplt = os.path.join(foldin, str(ls_usms['lsplt'][i]))
    mixID = str(ls_usms['mixID'][i])
    tabSpe = pd.read_excel(path_lsplt, sheet_name=mixID)
    ls_Spe = tabSpe["ongletP"].tolist()
    #garde onglet pour les noms
    ongletP = ls_Spe[0] #str(ls_usms['ongletP'][i])
    ongletPvois = ls_Spe[1] #str(ls_usms['ongletVoisin'][i])]


    path_scenar = os.path.join(foldin, fscenar)#(path_, 'input', fscenar)
    testsim[name].mn_sc = path_scenar

    path_variance_geno = os.path.join(foldin, fsd)#(path_, 'input', fsd)
    testsim[name].path_variance_geno = path_variance_geno

    path_variance_matrix = os.path.join(foldin, fsdx)
    testsim[name].path_variance_matrix = path_variance_matrix

    # la, lire scenario et changer parametres
    idscenar1 = int(ls_usms['scenario1'][i])
    idscenar2 = int(ls_usms['scenario2'][i])
    idscenar3 = int(ls_usms['scenario3'][i])
    idscenar4 = int(ls_usms['scenario4'][i])
    idscenar5 = int(ls_usms['scenario5'][i])
    idscenar6 = int(ls_usms['scenario6'][i])

    idscenar1_sd = int(ls_usms['scenario1_sd'][i])
    idscenar2_sd = int(ls_usms['scenario2_sd'][i])
    idscenar3_sd = int(ls_usms['scenario3_sd'][i])
    idscenar4_sd = int(ls_usms['scenario4_sd'][i])
    idscenar5_sd = int(ls_usms['scenario5_sd'][i])
    idscenar6_sd = int(ls_usms['scenario6_sd'][i])

    #ongletScenar2 = ongletPvois  # fait porter les changements sur fichier parametre voisin
    #ongletScenar1 = ongletP

    # sol
    path_sol = os.path.join(foldin, str(ls_usms['sol'][i]))#(path_, 'input', str(ls_usms['sol'][i]))
    ongletS = str(ls_usms['ongletS'][i])
    par_SN, par_sol = IOxls.read_sol_param(path_sol, ongletS)
    par_SN['concrr'] = 0.  # force eau de pluie dans ls test (a retirer)
    # testsim[name].ongletS = str(ls_usms['ongletS'][i])
    testsim[name].par_SN = par_SN
    testsim[name].par_sol = par_sol

    #station
    path_station = os.path.join(foldin, fsta)
    testsim[name].path_station = path_station
    testsim[name].ongletSta = ongletSta


    optdamier = int(ls_usms['optdamier'][i])
    nbcote = int(ls_usms['nbcote'][i])

    if str(ls_usms['typearrangement'][i]) == 'damier8':
        arrang = 'damier' + str(optdamier)
    if str(ls_usms['typearrangement'][i]) == 'damier16':
        arrang = 'damidouble' + str(optdamier)
    elif str(ls_usms['typearrangement'][i]) == 'row4':
        arrang = 'row' + str(optdamier)
    else:
        arrang = str(ls_usms['typearrangement'][i]) + str(optdamier)

    # nommix reste prevu pour melange a 2 !! -> reprendre avec Ls_Spe
    nommix = '_' + ongletP + '-' + ongletPvois + '_' + arrang + '_scenario' + str(idscenar2) + '-' + str(idscenar1)

    #testsim[name].ongletP = ongletP
    #testsim[name].ongletPvois = ongletPvois
    testsim[name].ls_Spe = ls_Spe
    testsim[name].nbcote = nbcote
    testsim[name].opt_sd = int(ls_usms['opt_sd'][i])  # option lue dans ls_usm
    testsim[name].opt_scenar = int(ls_usms['opt_scenar'][i])  # option lue dans ls_usm
    testsim[name].cote = float(ls_usms['cote'][i])
    testsim[name].deltalevmoy = float(ls_usms['deltalevmoy'][i])
    testsim[name].deltalevsd = float(ls_usms['deltalevsd'][i])
    testsim[name].typearrangement = str(ls_usms['typearrangement'][i])
    testsim[name].optdamier = optdamier
    testsim[name].ls_idscenar = [idscenar1, idscenar2, idscenar3, idscenar4, idscenar5, idscenar6]
    #testsim[name].idscenar1 = idscenar1
    #testsim[name].idscenar2 = idscenar2
    testsim[name].ls_idscenar_sd = [idscenar1_sd, idscenar2_sd, idscenar3_sd, idscenar4_sd, idscenar5_sd, idscenar6_sd]
    #testsim[name].ongletScenar2 = ongletScenar2
    #testsim[name].ongletScenar1 = ongletScenar1
    testsim[name].idscenar1_sd = idscenar1_sd
    testsim[name].idscenar2_sd = idscenar2_sd
    testsim[name].Rseed = seednb
    testsim[name].DOYdeb = int(ls_usms['DOYdeb'][i])
    testsim[name].DOYend = int(ls_usms['DOYend'][i])


    #mise a jour des options de simulation
    testsim[name].opt_residu = int(dic_opt['opt_residu'])  # si 0, pas activation de mineralisation
    testsim[name].opt_sd = int(dic_opt['opt_sd'])  # 1 #genere distribution des valeurs de parametres
    testsim[name].opt_covar = int(dic_opt['opt_covar'])  #definie matrice de cavariance a lire dans path_variance_matrix (0 opt_sd generere tirages independants)
    testsim[name].opt_shuffle = int(dic_opt['opt_shuffle']) # 1: for random order of plant species in ParamP ; 0: reular order
    testsim[name].opt_stressN = int(dic_opt['opt_stressN'])  # Active stress N; 1 = stress NNI actif (0= calcule, mais pas applique)
    testsim[name].opt_stressW = int(dic_opt['opt_stressW'])  # Active stressW; 1 = stress FTSW actif (0= calcule, mais pas applique)
    testsim[name].opt_ReadstressN = int(dic_opt['opt_ReadstressN'])  # Force stress N to read input values - for debugging/calibration
    testsim[name].opt_ReadstressW = int(dic_opt['opt_ReadstressW'])  # Force stress FTSW to read input values - for debugging/calibration
    testsim[name].opt_photomorph = int(dic_opt['opt_photomorph'])  # 1 #Activate photomorphogenetic effects on organ growth; 1 Actif (0= calcule, mais pas applique)
    testsim[name].opt_optT = int(dic_opt['opt_optT']) # option de calcul du cumul de temperature (0=betaD; 1=betaH; 2=lineaireD)
    testsim[name].opt_stressGel = int(dic_opt['opt_stressGel']) #Active gel stress option below Tgel
    testsim[name].opt_PP = int(dic_opt['opt_PP']) # Active photoperiodic effects (1 active; 0 inactive)
    testsim[name].opt_Nuptake = int(dic_opt['opt_Nuptake']) #options for calculating plant N uptake - 0:'STICS'  #1:'LocalTransporter'  #2:'old'
    testsim[name].opt_Mng = int(dic_opt['opt_Mng'])  # type of management file to be read: 0: default observed file ; 1: automatic management file #must be consistent with the management file!
    testsim[name].opt_ReadPP = int(dic_opt['opt_ReadPP'])  # Force photoperiod to read input values in management - for indoor experiment
    testsim[name].visu_root = int(dic_opt['visu_root'])  # 1# pour visualisation/interpretation root
    testsim[name].visu_shoot = int(dic_opt['visu_shoot'])  # 1# pour visualisation/interpretation shoot
    testsim[name].visu_leaf = int(dic_opt['visu_leaf'])  # 1# pour visualisation/interpretation feuilles slmt
    testsim[name].visu_sol = int(dic_opt['visu_sol'])  # 1# pour visualisation/interpretation sol
    testsim[name].visu_solsurf = int(dic_opt['visu_solsurf'])  # 0 pour visualisation du pattern
    testsim[name].frDisplay = int(dic_opt['frDisplay'])  # 1 #sauvegarde de la derniere vue
    testsim[name].movDisplay = int(dic_opt['movDisplay'])  # #sauvegarde toutes les vues pour faire un film
    testsim[name].opt_zip = int(dic_opt['opt_zip'])  # if 1, zip and delete the output csv files
    testsim[name].opt_verbose = int(dic_opt['opt_verbose'])  # 0, remove print in the console


    # mise a jour derivartionLength & axiom
    testsim[name].derivationLength = int(ls_usms['DOYend'][i]) - int(ls_usms['DOYdeb'][i])  # derivationLength variable predefinie dans L-py
    arr = str(ls_usms['typearrangement'][i])
    if arr == 'row4':  # carre rang heterogene
        nbplantes = nbcote * 4
    elif arr == 'row4_sp1' or arr == 'row4_sp2' :
        nbplantes = nbcote * 2
    elif arr == 'damier8' or arr == 'damier16' or arr == 'homogeneous' or arr == 'random8' or arr == 'damier9' or arr == 'damier10' or arr == 'damier8_4':  # carre homogene
        nbplantes = nbcote * nbcote
    elif arr == 'damier8_sp1' or arr == 'damier8_sp2' or arr == 'damier16_sp1' or arr == 'damier16_sp2':
        #prends pour le moment que le cas 4 (50/50)
        nbplantes = int(nbcote * nbcote / 2)
    else:
        print('unknown arrangement and nbplant')

    a = lpy.AxialTree()
    a.append(testsim[name].attente(1))
    for j in range(0, nbplantes):
        a.append(testsim[name].Sd(j))

    testsim[name].axiom = a  # passe un axial tree, pas de chaine de caractere

    #pareil : nom prevu pour 2 sp -> reprendre
    if int(ls_usms['opt_sd'][i]) == 1 or int(ls_usms['opt_sd'][i]) == 2:
        sdname = '_SD' + str(idscenar2_sd) + '-' + str(idscenar1_sd)
    else:
        sdname = '_-'

    # path fichiers de sortie
    testsim[name].path_out = path_OUT #os.path.join(path_, str(ls_usms['folder_out'][i]))
    testsim[name].outvarfile = 'toto_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.csv'
    testsim[name].lsorgfile = 'lsAxes_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.csv'
    testsim[name].outHRfile = 'outHR_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.csv'
    testsim[name].resrootfile = 'resroot_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.csv'
    testsim[name].outBilanNfile = 'BilanN_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.csv'
    testsim[name].outimagefile = 'scene_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + sdname + '_' + '.bmp'  # 'scene.bmp'
    testsim[name].outsdfile = 'paramSD_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + '_' + sdname + '_' + '.csv'
    testsim[name].outMngfile = 'MngAuto_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + '_' + sdname + '_' + '.csv'

    # plante si dossier out pas cree
    # pourrait faire la lecture les ls_usm directement dans le l-system pour faciliter...+
    return testsim #dico avec {nom:lsystem}


def runlsystem(lsys, name, clear=1):
    """ run the lsystem from a dict with its name"""
    try:
        lsys[name].derive()
        if clear==1:
            lsys[name].clear()

        print((''.join((name, " - done"))))
    except Exception as e:
        print(e)


def animatelsystem(lsys, name):
    lsys[name].animate()
    lsys[name].clear()
    print((''.join((name," - done"))))

def main():
    # fonction main pour pouvoir l'appeler en python pour l'exemple
    path_input = path_
    usm_file = 'liste_usms_exemple.xls' #ex fxls
    IDusm = 0

    #default exemple simulation
    foldinputs = os.path.join(path_, 'input')
    foldoutputs = os.path.join(path_, 'output')
    ongletB = 'exemple'
    #print(foldoutputs)


    #definition d'arguments avec getopt
    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:i:b:u:o:d", ["usm_file=", "inputs=", "onglet_batch=","usm_scenario=",  "outputs=",])#cf l-grass #"f:o:s:d" #!:d necessaire a la fin sinon lit pas dernier argument
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    #print("opts", opts)
    for opt, arg in opts:
        if opt in ("-f", "--file"):#usm file
            usm_file = arg
        elif opt in ("-i", "--inputs"):#inputs folder
            foldinputs = arg
        elif opt in ("-b", "--batch"):#onglet_batch
            ongletB = arg
        elif opt in ("-u", "--usm"):#usm ID
            IDusm = int(arg)
        elif opt in ("-o", "--outputs"):#outputs folder
            foldoutputs = arg

        # pour le moment output folde fourni dans usm -> a changer

    #mylsys = lsystemInputOutput_usm(path_input, usm_file, IDusm, foldin='multisim', ongletBatch='exemple')
    #print('foldoutputs', foldoutputs)
    mylsys = lsystemInputOutput_usm(usm_file, foldin=foldinputs, ongletBatch=ongletB, i=IDusm, path_OUT=foldoutputs)
    keyname = list(mylsys)[0]
    runlsystem(mylsys, keyname)
    #animatelsystem(mylsys, keyname)


if __name__ == '__main__':
    main()


#finir rendre accessible en externe le fichier de gestion des sorties!
# -> input and output folders -> fait
#rendre output accessible hors usm? -> fait!

#juste id usm: "python run_legume_usm.py -u 1" -> OK!
#test ttes les options: "python run_legume_usm.py -f liste_usms_exemple.xls -i C:\inputs -b exemple -u 1 -o C:\outputs" ->  OK!
#test sur le serveur: "python run_legume_usm.py -f liste_usms_exemple.xls -i /nfs/work/inra_ea/P3F/VGL/inputs -b exemple -u 1 -o /nfs/work/inra_ea/P3F/VGL/outputs"

#a tester en remplacement dans le batch! -> a faire
#prevoir de tout mettre dans le meme dossier en entree -> fichiers a migrer


