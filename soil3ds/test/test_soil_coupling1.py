## test of a 1D soil with homogeneous root distribution requiring only soil3ds (no dependency)

import os
import sys


import soil3ds

path_ = os.path.dirname(os.path.abspath(soil3ds.__file__))  # local absolute path of soil3ds
print(('path', path_))
sys.path.insert(0, path_)

from soil3ds import soil_moduleN as solN


import numpy as np
import IOxls
import IOtable





#################
## exemple 1: initialisation sol personalisee en n couches + fonctionnement VGL 'Local Transporter' avec racines homogenes
#################

# lecture des parametre sol directement a partir d'un fichier sol
# personalisation des initialisations (sol personnalise / creation solN)
# loop de 100 jours avec option 'Local Transporter' pour plant uptake

foldin =  os.path.join(path_, 'test','inputs')
fxls_sol = 'Parametres_sol_exemple.xls'
ongletS = 'lusignan99'

path_sol = os.path.join(path_,foldin,fxls_sol)
par_SN, par_sol = IOxls.read_sol_param(path_sol, ongletS)

# initialisation taille scene / discretisation (1D - homogene pour ttes les couches)
cote = 100 #cm
pattern8 = [[0, 0], [cote, cote]]
Lsol = max((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
largsol = min((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
#surfsolref = Lsol * largsol  # m2
dz_sol = 5. #cm
ncouches_sol = 20
discret_solXY = [1,1] #nombre de voxel selon X et Y
dxyz=[[Lsol / discret_solXY[0]] * discret_solXY[0], [largsol / discret_solXY[1]] * discret_solXY[1], [dz_sol / 100.]* ncouches_sol]
#lims_sol = rtd.lims_soil(pattern8, dxyz=[[Lsol / discret_solXY[0]] * discret_solXY[0], [largsol / discret_solXY[1]] * discret_solXY[1], [dz_sol / 100.] * ncouches_sol])


# vecteurs pour initialisation des propietes des couches de sol
vsoilnumbers = [1]*ncouches_sol #numeros de sol du profil -> mesures acsyd11
vDA = [1.31]*ncouches_sol #densite apparente de sol (mesure pesees initial aschyd11)
vMO = [par_SN['MO0_30']]*ncouches_sol
vCN = [par_SN['CN0_30']]*ncouches_sol
vARGIs = [par_SN['ARGIs']]*ncouches_sol
vCALCs = [par_SN['CALCs']]*ncouches_sol
vNH4 = [2.]*ncouches_sol # #!! kg d'N.ha-1 (entree de STICS)
vNO3 = [2.]*ncouches_sol # kg d'N.ha-1 (entree de STICS)
HRpinit = [] #default


# lecture meteo / mng journaliere = fixe
meteo_j = {'TmoyDay': 11., 'RG': 846.7, 'Et0': 2.0, 'Precip': 0., 'Tmin': '', 'Tmax': '', 'Tsol': 10., 'I0': 47.04, 'durjour': 10.88}
mng_j = {'Coupe': 0.0, 'Irrig': 0.0, 'FertNO3': 0.0, 'FertNH4': 0.0, 'Hcut': 3.0, 'ForceNNI': 1.0}
Tsol = meteo_j['Tsol']


# options pour simul sol
opt_residu = 0 #sans residus a mineraliser
opt_Nuptake = 1 #Local Transporters


#############
#init de l'objet sol
S = solN.SoilN(par_sol, par_SN, soil_number=vsoilnumbers,
                   dxyz=dxyz, vDA=vDA, vCN=vCN, vMO=vMO, vARGIs=vARGIs, vNO3=vNO3,
                   vNH4=vNH4, vCALCs=vCALCs, Tsol=Tsol, obstarac=None, pattern8=pattern8)



###### Creation variables plante entree pour simul plante-sol: ls_epsi / ls_roots / concentration N racine ou NNI = fixes
nb_plt = 200 #nombre de plantes dans pattern8
ls_epsi = [0.5/200]*nb_plt #fration de espsi.plant-1 (equivalent a transmis de 50%)
MSrac_plt = np.array([10./200]*nb_plt) #g plt-1 (equivalent a 1 T.ha-1)
SRL = 250 #m.g-1
LENrac_plt = MSrac_plt*SRL #m.plt-1
ls_N = np.array([0.75]*nb_plt)#invar['NNI']

#calcul de ls_roots adapte format sol S avec distribution homogene dans tout le sol
ls_roots = [] # cm !!
for i in range(nb_plt):
    #longueur de racine
    rootLen_i = S.m_1 * LENrac_plt[i]/ncouches_sol*100
    ls_roots.append(rootLen_i)

#np.sum(ls_roots[0])
#LENrac_plt[0]

#lecture parametre plante ParamP utilise dans les calculs du sol
path_plante = os.path.join(path_,foldin,'Parametres_plante_exemple.xls')#'Parametres_plante_v5cLucas.xls')#'Parametres_plante_v18.xls')#'Parametres_plante_v9Lucas_debugL.xls')#r'H:\devel\grassland\grassland\L-gume\Parametres_plante_v5cLucas.xls'
ongletP = 'Fix2'
g4 = IOxls.read_plant_param(path_plante, ongletP)
ParamP = [g4]*nb_plt
#utilise pourquoi / quel param precisement utilise dans le sol? -> revoir pour rendre explicite


######### loop pour n_jour
n_jour = 100
for j in range(n_jour):
    meteo_j = {'TmoyDay': 11., 'RG': 846.7, 'Et0': 2.0, 'Precip': 0., 'Tmin': '', 'Tmax': '', 'Tsol': 10., 'I0': 47.04, 'durjour': 10.88}
    mng_j = {'Coupe': 0.0, 'Irrig': 0.0, 'FertNO3': 0.0, 'FertNH4': 0.0, 'Hcut': 3.0, 'ForceNNI': 1.0}

    # Step Sol avec les inputs prevues dans VGL
    tag_inputs_soil_step = [S, par_SN, meteo_j, mng_j, ParamP, ls_epsi, ls_roots, ls_N, opt_residu, opt_Nuptake]  # input tag
    S, stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, temps_sol = solN.step_bilanWN_solVGL(*tag_inputs_soil_step)


    #print de sorties N
    kgNO3solHa = S.m_NO3.sum() / S.surfsolref*10000
    kgNH4solHa = S.m_NH4.sum() / S.surfsolref * 10000
    kgNsolHa = kgNO3solHa + kgNH4solHa
    uptNO3PltHa = ls_Act_Nuptake_plt[0].sum()*nb_plt/S.surfsolref*10000
    lix = S.lixiNO3  # /S.surfsolref*10000
    uptNO3Pltmole = ls_Act_Nuptake_plt[0].sum()*1000000*1000/(14*24) #en µmole d'N / plante / h

    cumMinN_j = S.bilanN['cumMinN'][-1]
    Lix_j = S.bilanN['cumLix'][-1]
    UptakePlt_j = S.bilanN['cumUptakePlt'][-1].sum()
    azomes = S.bilanN['azomes'][-1]
    MinN = S.bilanN['cumMinN'][-1]

    #print('N', j, azomes, kgNO3solHa, kgNH4solHa, Lix_j, MinN, UptakePlt_j)
    #print('N', j, azomes,kgNsolHa, kgNO3solHa,kgNH4solHa, lix, Lix_j, uptNO3PltHa, UptakePlt_j)
    #print('N', j, uptNO3PltHa, UptakePlt_j)
    #print('N plt', j, np.array(list(map(np.sum, ls_Act_Nuptake_plt)))) # en kg N/plant-1.d-1
    print('N mole', j, uptNO3Pltmole)  # en µmole d'N / plante / h

    # print de sorties W
    transp = sum(ls_transp)
    tsw = S.tsw_t.sum()

    #print('Water', j, tsw,transp, stateEV)







