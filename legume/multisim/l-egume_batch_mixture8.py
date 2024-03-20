#batch pour L-egume 23/12/17
#version GL utilisee pour melanges binaires (8*8)

#import the modules necessary to initiate the L-systems
from openalea.lpy import *
import multiprocessing

import os
import sys

try:
    import legume
    path_ = os.path.dirname(os.path.abspath(legume.__file__))#local absolute path of L-egume
except:
    path_ = r'C:\devel\l-egume\legume'#r'C:\devel\grassland'

print(('path', path_))

sys.path.insert(0, path_)
import IOxls
import IOtable


global foldin, fxls, ongletBatch, fscenar
#to define if used for multisimulation or non-regression tests
opttest = 4#2#5#1#4#2#'autre'#'exemple'#0#13#'exemple_BA'#
if opttest == 1 or opttest==2 or opttest==3 or opttest==4 or opttest==5: #si multisim des test de non regression (1 or 2)
    #global foldin, fxls, ongletBatch, fscenar
    foldin = 'test\inputs'
    fxls = 'liste_usms_nonregression.xls'
    fscenar = 'liste_scenarios.xls'
    fsd = 'exemple_sd.xls'
    if opttest == 1:#non regression
        ongletBatch = 'test'
    elif opttest == 2:#obssim
        ongletBatch = 'valid'
    elif opttest == 3:#solnu
        ongletBatch = 'solnu'#
    elif opttest == 4:#champ
        ongletBatch = 'test_champ'#
    elif opttest == 5:#pour bea
        ongletBatch = 'test_beajul'#
elif opttest == 'exemple':
    #global foldin, fxls, ongletBatch, fscenar
    foldin = 'multisim'
    fxls = 'liste_usms_exemple.xls'
    ongletBatch = 'exemple'
    fscenar = 'liste_scenarios_exemple.xls'
    fsd = 'exemple_sd.xls'
elif opttest == 'exemple_BA':
    # global foldin, fxls, ongletBatch, fscenar
    foldin = 'multisim'
    fxls = 'liste_usms_exemple_BA.xls'
    ongletBatch = 'exemple'
    # fscenar = 'liste_scenarios_exemple.xls'
    fsd = 'exemple_sd.xls'
else: #other multisimulation to be defined (0)
    #global foldin, fxls, ongletBatch, fscenar
    #to be manually updated
    foldin = 'input'#'multisim'
    fxls = 'liste_usms_essais.xls'#'liste_usms_mix.xls'
    ongletBatch = 'Champs'#'SimTest'
    fscenar = 'liste_scenarios.xls'
    fsd = 'exemple_sd.xls'



#lecture de la liste des usm
#path_ = r'H:\devel\grassland\grassland\L-gume'
mn_path = os.path.join(path_, foldin, fxls)#(path_,'test\inputs','liste_usms_nonregression.xls')#(path_,'multisim','liste_usms_mix.xls')#(path_,'liste_usms_mix_these lucas.xls')#
#ongletBatch = 'test'#'SimTest'#'complement'#'Feuil1'#'Sensi'#
usms = IOxls.xlrd.open_workbook(mn_path)
ls_usms = IOtable.conv_dataframe(IOxls.get_xls_col(usms.sheet_by_name(ongletBatch)))


#lecture scenario de changement des parametres (passer dans le script -> pourrait ne le faire que dans le batch)
mn_sc = os.path.join(path_,'liste_scenarios.xls')#(path_,'liste_scenarios_these_lucas.xls')#
#ongletBatch = 'Feuil1'
#usc = IOxls.xlrd.open_workbook(mn_sc)
#ls_sc = IOtable.conv_dataframe(IOxls.get_xls_col(usc.sheet_by_name(ongletBatch)))



#modif script pour pouvoir lancer en 1 ligne de commande python avec bon noms de fichiers...
# pb exit... ecriture fichier ou pb e dossier sortie?
#+deja dossier test a la racine

#cree la liste de L-systems et liste des noms
testsim={}
names = []
for i in range(len(ls_usms['ID_usm'])):
    if int(ls_usms['torun'][i]) == 1:#si 1 dans la colonne 'torun' l'ajoute a la liste
        name = str(int(ls_usms['ID_usm'][i]))+'_'+str(ls_usms['l_system'][i])[0:-4]
        seednb = int(ls_usms['seed'][i])
        names.append(name)
        path_lsys = os.path.join(path_, str(ls_usms['l_system'][i]))
        testsim[name]=Lsystem(path_lsys)  #objet l-system
        
        
        #testsim[name].ongletM = str(ls_usms['ongletM'][i])
        meteo_path_ =  os.path.join(path_, 'input',str(ls_usms['meteo'][i]))
        ongletM_ = str(ls_usms['ongletM'][i])
        testsim[name].meteo = IOxls.read_met_file(meteo_path_, ongletM_)

        #testsim[name].ongletMn = str(ls_usms['ongletMn'][i])
        mn_path_ = os.path.join(path_, 'input',str(ls_usms['mng'][i]))
        ongletMn_ = str(ls_usms['ongletMn'][i])
        testsim[name].mng = IOxls.read_met_file(mn_path_, ongletMn_)

        ini_path_ = os.path.join(path_, 'input',str(ls_usms['inis'][i]))
        ongletIni_ = str(ls_usms['ongletIn'][i])
        testsim[name].inis = IOxls.read_plant_param(ini_path_, ongletIni_)

        #testsim[name].ongletP = str(ls_usms['ongletP'][i])
        path_plante = os.path.join(path_, 'input',str(ls_usms['plante'][i]))
        ongletP = str(ls_usms['ongletP'][i])
        ongletPvois = str(ls_usms['ongletVoisin'][i])
        testsim[name].path_plante = path_plante
        
        path_variance_geno = os.path.join(path_, 'input', fsd)
        testsim[name].path_variance_geno = path_variance_geno
        #la, lire scenario et changer parametres
        idscenar1 = int(ls_usms['scenario1'][i])
        idscenar2 = int(ls_usms['scenario2'][i])
        idscenar1_sd = int(ls_usms['scenario1_sd'][i])
        idscenar2_sd = int(ls_usms['scenario2_sd'][i])
        ongletScenar2 = ongletPvois #fait porter les changements sur fichier parametre voisin
        ongletScenar1 = ongletP

        #sol
        path_sol = os.path.join(path_, 'input',str(ls_usms['sol'][i]))
        ongletS = str(ls_usms['ongletS'][i])
        par_SN, par_sol = IOxls.read_sol_param(path_sol, ongletS)
        par_SN['concrr'] = 0.  # force eau de pluie dans ls test (a retirer)
        #testsim[name].ongletS = str(ls_usms['ongletS'][i])
        testsim[name].par_SN = par_SN
        testsim[name].par_sol = par_sol

        #nbcote=7 # a passer ext
        #testsim[name].ParamP = [g4]*int(ls_usms['nbplt'][i])#nbcote*nbcote
        optdamier = int(ls_usms['damier'][i])
        nbcote = int(ls_usms['nbcote'][i])
        ### testsim[name].ParamP = damier8(g4,g5,opt=optdamier)
        if str(ls_usms['arrangement'][i]) == 'damier8':
            arrang = 'damier'+str(optdamier)
        elif str(ls_usms['arrangement'][i]) == 'row4':
            arrang = 'row'+str(optdamier)
        else:
            arrang = str(ls_usms['arrangement'][i])+str(optdamier)

        nommix = '_'+ongletP+'-'+ongletPvois+'_'+arrang+'_scenario'+str(idscenar2)+'-'+str(idscenar1)
        

        testsim[name].ongletP = ongletP
        testsim[name].ongletPvois = ongletPvois
        testsim[name].nbcote = nbcote
        testsim[name].opt_sd = int(ls_usms['opt_sd'][i])#1
        testsim[name].cote = float(ls_usms['cote'][i])
        testsim[name].deltalevmoy = float(ls_usms['retard'][i])
        testsim[name].deltalevsd = float(ls_usms['sd_retard'][i])
        testsim[name].typearrangement = str(ls_usms['arrangement'][i])
        testsim[name].optdamier = optdamier
        testsim[name].idscenar1 = idscenar1
        testsim[name].idscenar2 = idscenar2
        testsim[name].ongletScenar2 = ongletScenar2
        testsim[name].ongletScenar1 = ongletScenar1
        testsim[name].idscenar1_sd = idscenar1_sd
        testsim[name].idscenar2_sd = idscenar2_sd
        testsim[name].Rseed = seednb
        testsim[name].DOYdeb = int(ls_usms['DOYdeb'][i])
        testsim[name].DOYend = int(ls_usms['DOYend'][i])
        
        #mise a jour derivartionLength & axiom
        testsim[name].derivationLength = int(ls_usms['DOYend'][i]) - int(ls_usms['DOYdeb'][i])#derivationLength variable predefinie dans L-py
        if str(ls_usms['arrangement'][i]) == 'row4':#carre rang heterogene
            nbplantes = nbcote * 4
        else: #carre homogene
            nbplantes = nbcote*nbcote

        a=AxialTree()
        a.append(testsim[name].attente(1))
        for j in range(0,nbplantes):
           a.append(testsim[name].Sd(j))

        testsim[name].axiom = a#passe un axial tree, pas de chaine de caractere

        if int(ls_usms['opt_sd'][i])==1:
            sdname = '_SD'+str(idscenar2_sd)+'-'+str(idscenar1_sd)
        else:
            sdname = '_-'

        #path fichiers de sortie
        testsim[name].path_out = os.path.join(path_, str(ls_usms['folder_out'][i]))
        testsim[name].outvarfile = 'toto_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'_'+str(ls_usms['ongletM'][i])+sdname+'_'+'.csv'
        testsim[name].lsorgfile = 'lsAxes_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'_'+str(ls_usms['ongletM'][i])+sdname+'_'+'.csv'
        testsim[name].outHRfile = 'outHR_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'_'+str(ls_usms['ongletM'][i])+sdname+'_'+'.csv'
        testsim[name].resrootfile = 'resroot_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'_'+str(ls_usms['ongletM'][i])+sdname+'_'+'.csv'
        testsim[name].outBilanNfile = 'BilanN_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'_'+str(ls_usms['ongletM'][i])+sdname+'_'+'.csv'
        testsim[name].outimagefile = 'scene_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'_'+str(ls_usms['ongletM'][i])+sdname+'_'+'.bmp'#'scene.bmp'
        testsim[name].outsdfile = 'paramSD_' + name + nommix + '_' + str(ls_usms['ongletMn'][i]) + '_' + str(seednb) + '_' + str(ls_usms['ongletM'][i]) + '_'+sdname+'_'+'.csv'

        #plante si dossier out pas cree
        #pourrait faire la lecture les ls_usm directement dans le l-system pour faciliter...+

nb_usms =len(names)#len(ls_usms['ID_usm'])#len(names)#

#print nb_usms, names


#function to run an L-system from the 'testsim' dictionnary
def runlsystem(n):
    testsim[names[n]].derive()
    testsim[names[n]].clear()
    print((''.join((names[n]," - done"))))

def animatelsystem(n):
    testsim[names[n]].animate()
    testsim[names[n]].clear()
    print((''.join((names[n]," - done"))))


#run the L-systems

if __name__ == '__main__':
    multiprocessing.freeze_support()
    CPUnb=multiprocessing.cpu_count()-1 #nombre de processeurs, moins un par prudence. (et pour pouvoir faire d'autres choses en meme temps)
    print('nb CPU: '+str(CPUnb))
    pool = multiprocessing.Pool(processes=CPUnb)
    for i in range(int(nb_usms)):
        pool.apply_async(runlsystem, args=(i,)) #Lance CPUnb simulations en meme temps, lorsqu'une simulation se termine elle est immediatement remplacee par la suivante
        #runlsystem(i) #pour debug hors multisim (messages d'ereur visible)
        #animatelsystem(i)  # pour debug hors multisim (messages d'ereur + sortie archi visible)
    pool.close()
    pool.join()
