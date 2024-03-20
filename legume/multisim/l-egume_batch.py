#modif legume batch to reuse run_l-egume_usm

# import the modules necessary to initiate the L-systems
from openalea.lpy import *
import multiprocessing

import os
import sys

try:
    import legume

    path_ = os.path.dirname(os.path.abspath(legume.__file__))  # local absolute path of L-egume
except:
    path_ = r'C:\devel\l-egume\legume'  # r'C:\devel\grassland'

print(('path', path_))

sys.path.insert(0, path_)
import IOxls
import IOtable
import run_legume_usm as runl


global foldin, fxls, ongletBatch
# to define if used for multisimulation or non-regression tests
opttest = 'exemple'#4#'exemple_BA'#'mayssa'#'sdBea'#'OATbea'#1#2#'Histor'#1#4 ##2#1#  5#2#'autre'#0#13#
if opttest == 1 or opttest == 2 or opttest == 3 or opttest == 4 or opttest == 5:  # si multisim des test de non regression (1 or 2)
    # global foldin, fxls, ongletBatch, fscenar
    foldin =  os.path.join(path_, 'input')#'test\inputs'
    fxls = 'liste_usms_nonregression.xls'
    if opttest == 1:  # non regression
        ongletBatch = 'test'
        foldout = os.path.join(path_, 'test\lastcheck')
    elif opttest == 2:  # obssim
        ongletBatch = 'valid'
        foldout = os.path.join(path_, 'test\lastvalidBis')
    elif opttest == 3:  # solnu
        ongletBatch = 'solnu'  #
        foldout = os.path.join(path_, 'output')
    elif opttest == 4:  # champ
        ongletBatch = 'test_champ'  #
        foldout = os.path.join(path_, 'test', 'test_champ')
    elif opttest == 5:  # pour bea
        ongletBatch = 'test_beajul'  #
        foldout = os.path.join(path_, 'test', 'test2')
elif opttest == 'exemple':
    # global foldin, fxls, ongletBatch, fscenar
    foldin =  os.path.join(path_, 'input')#'multisim'
    foldout = os.path.join(path_, 'output')
    fxls = 'liste_usms_exemple.xls'
    ongletBatch = 'exemple'#'exemple_old' #
    #fscenar = 'liste_scenarios_exemple.xls'
elif opttest == 'exemple_BA':
    # global foldin, fxls, ongletBatch, fscenar
    foldin =  os.path.join(path_, 'input')#'multisim'
    fxls = 'newUSM.xls'#'liste_usms_exemple_BA.xls'#
    ongletBatch = 'pur'#'exemple'
    foldout = os.path.join(path_, 'test', 'test2')
elif opttest == 'mayssa':
    # global foldin, fxls, ongletBatch, fscenar
    # to be manually updated
    foldin =  'C:\inputs\inputs mayssa\DIGITLUZ'#r'C:\inputs\inputs test variance BLW'#os.path.join(path_, 'input')#'input'  # 'multisim'
    fxls = 'liste_usms_test.xls'#'liste_usms_eval.xls'#'liste_usms_exemple.xls'#'liste_usms_essais.xls'  # 'liste_usms_mix.xls'
    ongletBatch = 'Test'#'valid'#'Param1GL'#'OATbea'#'Histor'#'Champs'  # 'SimTest'#
    foldout =  'C:\inputs\inputs mayssa\output'#os.path.join(path_, 'output')
else:  # to personalize - other multisimulation to be defined (0)
    # global foldin, fxls, ongletBatch, fscenar
    # to be manually updated
    foldin =  r'C:\inputs\inputs test variance BLW'#os.path.join(path_, 'input')#'input'  # 'multisim'
    fxls = 'liste_usms_exemple.xls'#'liste_usms_essais.xls'  # 'liste_usms_mix.xls'
    ongletBatch = 'Param1GL'#'OATbea'#'Histor'#'Champs'  # 'SimTest'#
    foldout = os.path.join(path_, 'output')

#scenar et fsd definis en dur dans la fonction lsystemInputOutput_usm
#fscenar = 'liste_scenarios.xls'
#fsd = 'exemple_sd.xls'



#lecture de la liste des usm
#path_ = r'H:\devel\grassland\grassland\L-gume'
mn_path = os.path.join(path_, foldin, fxls)#(path_,'test\inputs','liste_usms_nonregression.xls')#(path_,'multisim','liste_usms_mix.xls')#(path_,'liste_usms_mix_these lucas.xls')#
#ongletBatch = 'test'#'SimTest'#'complement'#'Feuil1'#'Sensi'#
usms = IOxls.xlrd.open_workbook(mn_path)
ls_usms = IOtable.conv_dataframe(IOxls.get_xls_col(usms.sheet_by_name(ongletBatch)))



# cree la liste de L-systems et liste des noms
testsim = {}
names = []
for i in range(len(ls_usms['ID_usm'])):
    if int(ls_usms['torun'][i]) == 1:  # si 1 dans la colonne 'torun' l'ajoute a la liste
        #mylsys = runl.lsystemInputOutput_usm(path_, fxls, i, foldin=foldin, ongletBatch=ongletBatch)
        mylsys = runl.lsystemInputOutput_usm(fxls, foldin=foldin, ongletBatch=ongletBatch, i=i, path_OUT=foldout)
        name = list(mylsys)[0]
        names.append(name)
        testsim[name] = mylsys[name]


nb_usms = len(names)  # len(ls_usms['ID_usm'])#len(names)#

# print (nb_usms, names)


#function to run an L-system from the 'testsim' dictionnary
def runlsystem(n):
    #testsim : global
    testsim[names[n]].derive()
    testsim[names[n]].clear()
    print((''.join((names[n]," - done"))))

#def runlsystem_name(name):
#    #testsim : global
#    testsim[name].derive()
#    testsim[name].clear()
#    print((''.join((name," - done"))))

#def runlsystem_obj(lsys):
#    """ run the lsystem from lsys object"""
#    lsys.derive()
#    lsys.clear()
#    #print((''.join((name, " - done"))))
#    #marche pas en multiprocessing



#run the L-systems

if __name__ == '__main__':
    multiprocessing.freeze_support()
    CPUnb=multiprocessing.cpu_count()-1 #nombre de processeurs, moins un par prudence. (et pour pouvoir faire d'autres choses en meme temps)
    print('nb CPU: '+str(CPUnb))
    pool = multiprocessing.Pool(processes=CPUnb)
    for i in range(int(nb_usms)):
        pool.apply_async(runlsystem, args=(i,)) #marche
        #runlsystem(i) #pour debug hors multisim (messages d'ereur visible)

        #pool.apply_async(runlsystem_name, args=(names[i],)) #marche aussi
        #pool.apply_async(runlsystem_obj, args=(testsim[names[i]],)) #marche pas!
        # pool.apply_async(runlsystem, args=(testsim[names[i]],)) #Lance CPUnb simulations en meme temps, lorsqu'une simulation se termine elle est immediatement remplacee par la suivante
        # pool.apply_async(runl.runlsystem, kwds = {'lsys':testsim ,'name': names[i]} ) #marche pas
        #runl.runlsystem(testsim, names[i]) #pour debug hors multisim (messages d'ereur visible)
        #runl.animatelsystem(testsim, names[i])  # pour debug hors multisim (messages d'ereur + sortie archi visible)
    pool.close()
    pool.join()

#marche en single run
#marche avec fonction locale a 1 argument si nom ou id..
#marche pas en multiple run avec runl.runlsystem et argument mutiples???
#marche pas avec deux arguments, meme si fonction locale et sans gestion d'exception -> vraiment la gestion des deux arguments??
#fonction a un argument marche si passe numero ou nom en args, pas si passe directement l'objet lsys? comme deepcopy - #RuntimeError: Pickling of "openalea.lpy.__lpy_kernel__.Lsystem" instances is not enabled
#car pickling pas autorise : https://stackoverflow.com/questions/27868395/python-multiprocessing-object-passed-by-value
#d'ou peut pas utiliser meme fonction que run_legume_usm et besoin d'un obljet liste de lsystem en variable globale