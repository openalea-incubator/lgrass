#single un pour L-egume permettant le lancement en ligne de commande avec options 15/08/18
#lecture des entrees dans dossier input / ecrture des sorties dans les fichiers par defaut du output


#import the modules necessary to initiate the L-systems
from openalea.lpy import *

import os
import sys
import argparse #pour parser les argument en ligne de commande

try:
    import legume
    path_ = os.path.dirname(os.path.abspath(legume.__file__))#local absolute path of L-egume
except:
    path_ = r'C:\devel\l-egume\legume'#r'C:\devel\grassland'


sys.path.insert(0, path_)
import IOxls
import IOtable

#print(sys.argv)
#options = sys.argv #recupere une liste d'arguments passes en ligne de commandes (non structure, sans aide)

parser = argparse.ArgumentParser() #permet de personnalise les options et le fomat attendu des entrees + aide
group = parser.add_mutually_exclusive_group() #deux options eclusive : usm/detail modes (par defaut: lance en usm avec les les valeurs par defaut)
group.add_argument("-u", "--usm", nargs=3, action="store", help="launch l-egume in usm mode", metavar=('usmsXLS', 'ongletXLS', 'usmNB'), default=('liste_usms_mix.xls', 'SimTest', '1'))#, nargs="*", action="store_true")
group.add_argument("-d", "--detail", help="launch l-egume in detail input mode", action="store_true")
# liste des option du mode detail
parser.add_argument("-lsys", help="l-system file (detail mode)", default="l-egume.lpy")
parser.add_argument("-met", help="meteo file (detail mode)", nargs=2, default=("meteo_exemple_debugL_gl.xls", "Lusignan30"),  metavar=('file', 'onglet'))
parser.add_argument("-mng", help="management file (detail mode)", nargs=2, default=("management_exemple3_debugL_gl.xls", "Lusignan30IrrN"),  metavar=('file', 'onglet'))
parser.add_argument("-ini", help="initialsation file (detail mode)", nargs=2, default=("Initialisation_sol_exemple.xls", "Lusignan30_5x5"),  metavar=('file', 'onglet'))
parser.add_argument("-sol", help="soil file (detail mode)", nargs=2, default=("Parametres_sol_exemple2_debugL_glbis.xls", "lusignan99"),  metavar=('file', 'onglet'))
parser.add_argument("-plt", help="plant file (detail mode)", nargs=3, default=("Parametres_plante_v5cLucas.xls", "Fix1", "nonFix1"),  metavar=('file', 'ongletP', 'ongletVois'))
parser.add_argument("-doy", help="Start / End DOYs (detail mode)", nargs=2, default=("60", "335"),  metavar=('DOYdeb', 'DOYend'))
parser.add_argument("-scn", help="Scene options (detail mode)", nargs=4, default=("damier8", "4", "40.", "8"),  metavar=('typearrangement', 'optdamier', 'cote', 'nbcote'))
parser.add_argument("-sd", help="Seed options (detail mode)", default="0",  metavar=('seed'))
parser.add_argument("-rd", help="Random ini options (detail mode)", nargs=2, default=("30", "15"),  metavar=('deltalevmoy', 'deltalevsd'))
#...
args = parser.parse_args()


print ('args')#, args, args.usm, args.detail)
#print('u1', args.usm[1])
#print('options', options)


#create the l-system object
if not args.detail:#mode usm; pas en mode detail #len(options)>1:
    # idee: passer un nom d'usm (a aller lire dans un fichier xls ou tout est decrit) y compris le lsystem
    usmsXLS = args.usm[0]#options[1]
    ongletXLS = args.usm[1]#options[2]
    usmNB = args.usm[2]#options[3]

    # lecture de la liste des usm
    # path_ = r'H:\devel\grassland\grassland\L-gume'
    mn_path = os.path.join(path_, 'input', usmsXLS)  #
    ongletBatch = ongletXLS  # 'SimTest'#'Feuil1'#'Sensi'#
    usms = IOxls.xlrd.open_workbook(mn_path)
    ls_usms = IOtable.conv_dataframe(IOxls.get_xls_col(usms.sheet_by_name(ongletBatch)))

    #trouve le numero de ligne correspondant a usmNB
    id = -1
    for i in range(len(ls_usms['ID_usm'])):
        if str(int(ls_usms['ID_usm'][i])) == str(usmNB):
            id = i
            break
    # moins efficace que multi ca doit chercher la bonne ligne pour chaque simul( demander le umero de ligne a la place pour eviter la boucle?)

    #cree l'objet l-system
    # name = 'l-egume-1.0.37hquad_mixtures_initc.lpy'
    if id != -1:
        #name = str(ls_usms['l_system'][id])
        path_lsys = os.path.join(path_, str(ls_usms['l_system'][id]))
        testsim = Lsystem(path_lsys)  # objet l-system
    else:
        print ('erreur: usm inconnue')
else:# a faire
    # idee: ouvrir lsystem et specifie et ajuster au cas par cas les entree par defaut
    #name = args.lsys
    id = -1
    path_lsys = os.path.join(path_, args.lsys)
    testsim = Lsystem(path_lsys)

    #pass
    #reuperer les differentes valers d'entrees par defaut



#get the input variables values
if not args.detail and id !=-1 :#mode usm et id OK
    name = str(int(ls_usms['ID_usm'][id])) + '_' + str(ls_usms['l_system'][id])[0:-4]
    #path_lsys = os.path.join(path_, str(ls_usms['l_system'][id]))

    # meteo / mng / ini / plante input files
    meteo_path_ = os.path.join(path_, 'input', str(ls_usms['meteo'][id]))
    ongletM_ = str(ls_usms['ongletM'][id])

    mn_path_ = os.path.join(path_, 'input', str(ls_usms['mng'][id]))
    ongletMn_ = str(ls_usms['ongletMn'][id])

    ini_path_ = os.path.join(path_, 'input', str(ls_usms['inis'][id]))
    ongletIni_ = str(ls_usms['ongletIn'][id])

    path_plante = os.path.join(path_, 'input', str(ls_usms['plante'][id]))
    ongletP = str(ls_usms['ongletP'][id])
    ongletPvois = str(ls_usms['ongletVoisin'][id])

    sol_path_ = os.path.join(path_, 'input', str(ls_usms['sol'][id]))
    ongletS = str(ls_usms['ongletS'][id])

    #pas acces aux options de scenario
    ##  lire scenario et changer parametres
    #idscenar1 = int(ls_usms['scenario1'][id])
    #idscenar2 = int(ls_usms['scenario2'][id])
    #ongletScenar2 = ongletPvois  # fait porter les changements sur fichier parametre voisin

    # options de la scene
    typearrangement = str(ls_usms['arrangement'][id])
    optdamier = int(ls_usms['damier'][id])
    cote = float(ls_usms['cote'][id])
    nbcote = int(ls_usms['nbcote'][id])

    #options de simul
    seednb = int(ls_usms['seed'][id])
    DOYdeb = int(ls_usms['DOYdeb'][id])
    DOYend = int(ls_usms['DOYend'][id])
    derivationLength = int(ls_usms['DOYend'][id]) - int(ls_usms['DOYdeb'][id])
    deltalevmoy = float(ls_usms['retard'][id])
    deltalevsd = float(ls_usms['sd_retard'][id])

    #sorties
    #nommix = '_' + ongletP + '-' + ongletPvois + '_' + 'damier' + str(optdamier) + '_scenario' + str(idscenar2) + '-' + str(idscenar1)

    #en faire un dico (en faits'est est deja un qui est a mettre a jour) et passer le dico?
else:#mode detail
    name = str(args.lsys)[0:-4]

    # meteo / mng / ini / plante input files
    meteo_path_ = os.path.join(path_, 'input', args.met[0])
    ongletM_ = args.met[1]

    mn_path_ = os.path.join(path_, 'input', args.mng[0])
    ongletMn_ = args.mng[1]

    ini_path_ = os.path.join(path_, 'input', args.ini[0])
    ongletIni_ = args.ini[1]

    path_plante = os.path.join(path_, 'input', args.plt[0])
    ongletP = args.plt[1]
    ongletPvois = args.plt[2]

    sol_path_ = os.path.join(path_, 'input', args.sol[0])
    ongletS = args.sol[1]

    # options de la scene
    typearrangement = args.scn[0]#str(ls_usms['arrangement'][id])
    optdamier = int(args.scn[1])#int(ls_usms['damier'][id])
    cote = float(args.scn[2])#float(ls_usms['cote'][id])
    nbcote = int(args.scn[3])#int(ls_usms['nbcote'][id])
    #pas mal d'autres options de scene qui peuvent etre passees!

    # options de simul
    seednb = int(args.sd)
    DOYdeb = int(args.doy[0])
    DOYend = int(args.doy[1])
    derivationLength = int(args.doy[1]) - int(args.doy[0])
    deltalevmoy = float(args.rd[0])
    deltalevsd = float(args.rd[1])

    #pass
    #pas mal d'autres variables d'entrees encore a endre accessibles...


#update the default l-system object with optinal input values
if (not args.detail and id != -1) or args.detail: #tout sauf mode usm avec pb d'usm
    testsim.meteo = IOxls.read_met_file(meteo_path_, ongletM_)
    testsim.inis = IOxls.read_plant_param(ini_path_, ongletIni_)
    testsim.mng = IOxls.read_met_file(mn_path_, ongletMn_)
    testsim.par_SN, testsim.par_sol =  IOxls.read_sol_param(sol_path_, ongletS)
    testsim.par_SN['concrr'] = 0.  # force eau de pluie(serait a enlever)

    testsim.ongletS = ongletS#sert a rien: ce qu'il faut c'est lire fchier sol...
    testsim.ongletP = ongletP
    testsim.ongletPvois = ongletPvois
    testsim.nbcote = nbcote
    testsim.cote = cote
    testsim.deltalevmoy = deltalevmoy
    testsim.deltalevsd = deltalevsd
    testsim.typearrangement = typearrangement
    testsim.optdamier = optdamier
    #testsim.idscenar1 = idscenar1
    #testsim.idscenar2 = idscenar2
    #testsim.ongletScenar2 = ongletScenar2
    #testsim.ongletScenar1 = ongletScenar1
    testsim.Rseed = seednb
    testsim.DOYdeb = DOYdeb
    testsim.DOYend = DOYend
    testsim.derivationLength = derivationLength

    # mise a jour axiom: as active car nbplante pas bien calcule pour toutes les options?
    nbplantes = nbcote * nbcote
    a = AxialTree()
    a.append(testsim.attente(1))
    for j in range(0, int(nbplantes)):
        a.append(testsim.Sd(j))

    testsim.axiom = a  # passe un axial tree, pas de chaine de caractere

    #pas mise a jour des noms des sorties
    # testsim.path_out = os.path.join(path_, str(ls_usms['folder_out'][id]))
    # testsim.outvarfile = 'toto_' + name + nommix + '_' + str(ls_usms['ongletMn'][id]) + '_' + str(seednb) + '.csv'
    # testsim.lsorgfile = 'lsAxes_' + name + nommix + '_' + str(ls_usms['ongletMn'][id]) + '_' + str(seednb) + '.csv'
    # testsim.outHRfile = 'outHR_' + name + nommix + '_' + str(ls_usms['ongletMn'][id]) + '_' + str(seednb) + '.csv'
    # testsim.resrootfile = 'resroot_' + name + nommix + '_' + str(ls_usms['ongletMn'][id]) + '_' + str(seednb) + '.csv'
    # testsim.outBilanNfile = 'BilanN_' + name + nommix + '_' + str(ls_usms['ongletMn'][id]) + '_' + str(seednb) + '.csv'
    # testsim.outimagefile = 'scene_'+name+nommix+'_'+str(ls_usms['ongletMn'][i])+'_'+str(seednb)+'.bmp'#'scene.bmp'


#print ('name:', name)


#faire passer les autres valeurs de l'usm dans le lsystem!!
# en faire une fonction eventuellement reutilisable...



#function to run an L-system from the 'testsim' dictionnary
def runlsystem(testsim):
    testsim.derive()
    testsim.clear()
    print((''.join((name," - done"))))



if __name__ == '__main__':
    runlsystem(testsim)
    #pass


#en ligne de commande (option usm): python l-egume_run.py -u liste_usms.xls Onglet idUSM
#en ligne de commande (option detail): python l-egume_run.py -d -lsys l-egume.lpy
#option detail est encore.. a detailler


# https://docs.python.org/3/howto/argparse.html
#pourrait facilement rendre le batch acceesible en ligne de commande en passant les arguments des fichiers et onglets xls en options (l-egume_run_multi?)
#faire deux options: run -u (usm) ou -d (detail) qui prendraient des arguments differents?

#serait interessant d'ajouter des options -vparamNames, -vparamVals qui seraient a passer pour mettre a jour le fichier plante et utilisees pour faire de l'optimisation ciblee sur ces parametes
#avec un nags='*'
#solutin eventuelle pour lancer depuis intenet? installer serveur local sur ordi modelisation et lancer un scrip python/php simple depuis internet qui appelle le run.py?