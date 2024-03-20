#############
# GL - test sol 3D, avec racines 3D
# utilise racines du L-system
#############

import os
import soil3ds
from soil3ds import soil_moduleN as solN
from soil3ds import soil_wrapper as soil_interface

import openalea.lpy as lpy

path_ = os.path.dirname(os.path.abspath(soil3ds.__file__))  # path ou trouver les inputs
path_leg = os.path.join(path_, 'test')



# test external coupling archisimple eau-N
# Run de Archisimple avec un couplage interne a soil3DS dans le EndEach (eau et N: step_bilanWN_solVGL)
# couplage pas a tous les steps: seulement 1 fois par jour!


nomlsystem = 'ArchiSimple GL2_coupledWN.lpy'
lsysR = lpy.Lsystem(os.path.join(path_leg, nomlsystem))
# print(lsysR)

# simulation avec couplage sol
# active le calcule interne dans EndEach / desactive external coupling
lsysR.opt_external_coupling = 1
lsysR.opt_verbose = 0  # retire les print

# initialisation avec fonction du lsystem
lstring = lsysR.axiom  # initialise lstring avex l'axiom du l-system
nb_iter = lsysR.derivationLength
S, stateEV, intsoil, par_SN = lsysR.initiatisation_soil_default(lsysR.pattern8, lsysR.dz, lsysR.size, lsysR.stateEV, lsysR.properties_3ds)

# boucle pour un calcul step by step
for i in range(nb_iter + 1):
    # print('iter ',i)
    lstring = lsysR.derive(lstring, i, 1)
    lscene = lsysR.sceneInterpretation(lstring)

    ############
    # step couplage sol
    ############

    var_soil_coupling = lsysR.var_soil_coupling
    TPS, IDj, TTdays, S_1, intsoil_1, stateEV_1, Et0, Tsol, Rain, Irrig, FertNO3, FertNH4, epsi, ParamPN, ls_N, opt_residu, opt_Nuptake = var_soil_coupling

    if TPS in TTdays:  # fait step wtaer balance pas tous les degres jours: uniquement 1 fois par jour (TT indiques dans TTdays)

        j = IDj  # id du jour dans la liste meteo /mng

        # calcul grille densite sur intsoil_1 remplie dans le L-system
        ls_roots = [soil_interface.soil3Dw2s3DSprop(intsoil_1, S_1, 'root_length')]

        # update daily variables
        ls_epsi = [epsi[j]]
        meteo_j = {'Et0': Et0[j], 'Precip': Rain[j], 'Tsol': Tsol[j]}
        mng_j = {'Irrig': Irrig[j], 'FertNO3': FertNO3[j], 'FertNH4': FertNH4[j]}

        # step water et N balance avec 1 seule entite (ls_roots et epsi)
        tag_inputs_soil_step = [S, par_SN, meteo_j, mng_j, ParamPN, ls_epsi, ls_roots, ls_N, opt_residu, opt_Nuptake]  # input tag
        res_soil_step = solN.step_bilanWN_solVGL(*tag_inputs_soil_step)
        S, stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, temps_sol = res_soil_step

        print('tp water-N balance:', TPS)

        # step de 1 jour
        IDj = IDj + 1

    # mise a jour des variables du lsystem
    lsysR.S = S
    # lsysR.stateEV = stateEV

    print('ftsw', ls_ftsw[0])


##termes du bilan hydrique global
S.CloseWbalance()  # -> equilibre
S.CloseCbalance()  # -> equilibre
S.CloseNbalance()  # -> equilibre






