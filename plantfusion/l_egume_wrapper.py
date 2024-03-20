import os
from copy import deepcopy

import numpy
import pandas
from plantfusion.indexer import Indexer
import scipy

from openalea.lpy import *

import legume.IOxls as IOxls
import legume.IOtable as IOtable
import legume.run_legume_usm as runl
import legume.ShootMorpho as sh
import legume.daily_loop as loop
from legume.initialisation import init_plant_residues_fromParamP

import riri5.RIRI5 as riri

from plantfusion.utils import create_child_folder
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.indexer import Indexer


class L_egume_wrapper(object):
    """Wrapper for l-egume model

    construction creates the lsystem

    1 instance = 1 usm = 1 plant specy

    """

    def __init__(
        self,
        name="legume",
        indexer=Indexer(),
        in_folder="",
        out_folder=None,
        nameconfigfile="liste_usms_exemple.xls",
        ongletconfigfile="exemple",
        IDusm=None,
        planter=None,
        caribu_scene=False
    ) -> None:
        if out_folder is not None:
            try:
                os.mkdir(os.path.normpath(out_folder))
                print("Directory ", out_folder, " Created ")
            except FileExistsError:
                pass

            # output folder for l-egume
            out_folder = os.path.normpath(out_folder)
            self.out_folder = os.path.join(out_folder, name)
            try:
                os.mkdir(os.path.normpath(self.out_folder))
                print("Directory ", self.out_folder, " Created ")
            except FileExistsError:
                pass
            create_child_folder(self.out_folder, "brut")
            create_child_folder(self.out_folder, "graphs")

        else:
            self.out_folder = ""

        self.name = name
        self.indexer = indexer
        self.global_index = indexer.global_order.index(name)
        self.legume_index = indexer.legume_names.index(name)

        # read l-egume configuration files
        mn_path = os.path.join(in_folder, nameconfigfile)
        usms = IOxls.xlrd.open_workbook(mn_path)
        ls_usms = IOtable.conv_dataframe(IOxls.get_xls_col(usms.sheet_by_name(ongletconfigfile)))

        # create list of lsystem from config file
        if planter.generation_type == "default":
            planter = None

            if all([i == 1 for i in ls_usms["torun"]]):
                i = ls_usms["torun"].index(1)
                self.__load_lsystem(
                    nameconfigfile, in_folder, ongletconfigfile, i, os.path.join(self.out_folder, "brut"), planter
                )
        else:
            self.__load_lsystem(
                nameconfigfile,
                in_folder,
                ongletconfigfile,
                ls_usms["ID_usm"].index(IDusm),
                os.path.join(self.out_folder, "brut"),
                planter,
            )
        
        self.lstring = self.lsystem.axiom
        self.lsystem.opt_external_coupling = 1
        self.lsystem.opt_Nuptake = 0

        if caribu_scene:
            self.lsystem.visu_leaf = 1
            self.lsystem.visu_root = 0
            self.lsystem.visu_shoot = 0
            self.lsystem.visu_sol = 0
            self.lsystem.visu_solsurf = 0

        # in order to compute tag_inputs_loop
        lstrings_temp = self.lsystem.derive(self.lstring, 0, 1)

        self.number_of_species = len(self.lsystem.tag_loop_inputs[17])
        self.number_of_plants = len(self.lsystem.tag_loop_inputs[3])

        if self.number_of_species > 1 :
            if self.indexer.legume_number_of_species[self.legume_index] != self.number_of_species:
                self.indexer.legume_number_of_species[self.legume_index] = self.number_of_species
                self.indexer.update_legume_several_species(self.name)
            self.global_index = [index for index, item in enumerate(self.indexer.global_order) if item == self.name]
            self.legume_index = [index for index, item in enumerate(self.indexer.global_order) if item == self.name]


        if planter is not None:
            if isinstance(self.global_index, list):
             self.index_in_global_plants = [sum(planter.number_of_plants[: self.global_index[0]]),
                                                    sum(planter.number_of_plants[: self.global_index[-1] + 1])]
            else:
                self.index_in_global_plants = [
                    sum(planter.number_of_plants[: self.global_index]),
                    sum(planter.number_of_plants[: self.global_index + 1]),
                ]

        self.res_trans = None
        self.res_abs_i = None
        self.invar: dict = {}
        self.domain = None

    def __load_lsystem(self, nameconfigfile, in_folder, ongletconfigfile, i, path_OUT, planter=None):
        update_parameters = {}
        if planter is not None:
            update_parameters["typearrangement"] = planter.legume_typearrangement
            update_parameters["nbcote"] = planter.legume_nbcote[self.legume_index]
            update_parameters["cote"] = planter.legume_cote
            update_parameters["optdamier"] = planter.legume_optdamier

        mylsys = runl.lsystemInputOutput_usm(
            nameconfigfile,
            foldin=in_folder,
            ongletBatch=ongletconfigfile,
            i=i,
            path_OUT=path_OUT,
            update_usm_parameters=update_parameters,
        )
        name = list(mylsys)[0]
        self.simulation_name = list(mylsys)[0]
        self.lsystem = mylsys[name]

    def derive(self, t):
        self.lstring = self.lsystem.derive(self.lstring, t, 1)
        self.invar = self.lsystem.tag_loop_inputs[0]

    def light_inputs(self, elements="triangles"):
        if elements == "triangles":
            return self.lsystem.sceneInterpretation(self.lstring)

        elif elements == "voxels":
            leaf_area = self.lsystem.tag_loop_inputs[13]
            angle_distrib = self.lsystem.tag_loop_inputs[17]
            legume_grid = {"LA": leaf_area, "distrib": angle_distrib}
            return legume_grid

        else:
            print("Unknown light model")
            raise

    def light_results(self, energy, lighting: Light_wrapper, selective_global_index=None) -> None:
        if selective_global_index is not None:
            saved_global_index = self.global_index
            self.global_index = selective_global_index

        if lighting.lightmodel == "caribu":
            self.res_trans = self.transfer_caribu_legume(
                energy=energy,
                nb0=lighting.nb_empty_z_layers(),
                elements_outputs=lighting.results_organs(),
                sensors_outputs=lighting.results_sensors(),
            )

        elif lighting.lightmodel == "ratp":
            self.res_abs_i, self.res_trans = self.transfer_ratp_legume(
                energy, lighting.results_voxels(), lighting.nb_empty_z_layers()
            )

        elif lighting.lightmodel == "riri5":
            self.res_trans = lighting.res_trans()
            self.res_abs_i = lighting.res_abs_i()

        self.compute_plants_interception(lighting.results_organs(), energy, lighting.soil_energy())
        self.compute_potential_plant_growth()

        if selective_global_index is not None:
            self.global_index = saved_global_index

    def soil_inputs(self):
        ls_roots = self.lsystem.tag_loop_inputs[21]
        ParamP = self.lsystem.tag_loop_inputs[3]
        par_SN = self.lsystem.tag_loop_inputs[19]
        meteo_j = self.lsystem.tag_loop_inputs[6]
        mng_j = self.lsystem.tag_loop_inputs[7]
        opt_residu = self.lsystem.tag_loop_inputs[-2]
        opt_Nuptake = self.lsystem.opt_Nuptake
        soil = self.lsystem.tag_loop_inputs[18]

        ls_N = []
        if opt_Nuptake == 0 or opt_Nuptake == 2:  # 'STICS' or 'old':
            ls_N = self.demandeN

        elif opt_Nuptake == 1:
            ls_N = numpy.array(self.invar["NNI"])

        # gere l'aggregation des entrees par plante
        ls_epsi = self.epsi.tolist()
        ls_N = ls_N.tolist()

        # step soil en commun
        self.inputs_soil = [
            soil,
            par_SN,
            meteo_j,
            mng_j,
            ParamP,
            ls_epsi,
            ls_roots,
            ls_N,
            opt_residu,
            opt_Nuptake,
        ]
        N_content_roots_per_plant = ls_N
        roots_length_per_plant_per_soil_layer = ls_roots
        plants_soil_parameters = ParamP
        plants_light_interception = ls_epsi

        return (
            N_content_roots_per_plant,
            roots_length_per_plant_per_soil_layer,
            plants_soil_parameters,
            plants_light_interception,
        )

    def soil_results(self, results_soil, planter=None, selective_global_index=None) -> None:
        if selective_global_index is not None:
            saved_global_index = self.global_index
            self.global_index = selective_global_index
            if isinstance(self.global_index, list):
                self.index_in_global_plants = [sum(planter.number_of_plants[: self.global_index[0]]),
                                                    sum(planter.number_of_plants[: self.global_index[-1] + 1])]
            else:
                self.index_in_global_plants = [
                    sum(planter.number_of_plants[: self.global_index]),
                    sum(planter.number_of_plants[: self.global_index + 1]),
                ]
            self.global_index = saved_global_index

        if planter is not None:
            if isinstance(self.global_index, list):
                self.index_in_global_plants = [sum(planter.number_of_plants[: self.global_index[0]]),
                                                        sum(planter.number_of_plants[: self.global_index[-1] + 1])]
            else:
                self.index_in_global_plants = [
                    sum(planter.number_of_plants[: self.global_index]),
                    sum(planter.number_of_plants[: self.global_index + 1]),
                ]
            
        soil, stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, temps_sol = results_soil

        # if 
        ls_ftsw = ls_ftsw[self.index_in_global_plants[0] : self.index_in_global_plants[1]]
        ls_transp = ls_transp[self.index_in_global_plants[0] : self.index_in_global_plants[1]]
        ls_Act_Nuptake_plt_leg = ls_Act_Nuptake_plt[self.index_in_global_plants[0] : self.index_in_global_plants[1]]
        temps_sol = temps_sol[self.index_in_global_plants[0] : self.index_in_global_plants[1]]

        self.results_soil = [soil, stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt_leg, temps_sol]

    def compute_plants_interception(self, organs_results=None, energy=1, pari_soil_in=-1):
        surf_refVOX = self.lsystem.tag_loop_inputs[15]
        dicFeuilBilanR = self.lsystem.tag_loop_inputs[14]
        leaf_area = self.lsystem.tag_loop_inputs[13]

        if self.domain is None:
            surfsolref = self.lsystem.tag_loop_inputs[12]
        else:
            surfsolref = (self.domain[1][0] - self.domain[0][0]) * (self.domain[1][1] - self.domain[0][1])

        # R_FR voxel (calcul de zeta)
        tag_light_inputs2 = [self.res_trans / (energy * surf_refVOX)]  # input tag
        self.rfr = riri.rfr_calc_relatif(*tag_light_inputs2)

        # interception au sol
        if pari_soil_in < 0:
            transmi_sol = numpy.sum(self.res_trans[-1][:][:]) / (energy * surfsolref)
            pari_soil = max(1.0 - transmi_sol, 1e-15)
        else:
            pari_soil = 1 - pari_soil_in

        # interception des plantes
        # res_abs_i existe donc on est passé soit par ratp soit riri, il faut màj invar['parip']
        if self.res_abs_i is not None:
            dicFeuilBilanR = sh.calc_paraF(dicFeuilBilanR, leaf_area, self.res_abs_i)
            sh.calc_para_Plt(self.invar, dicFeuilBilanR)

        pari_canopy = numpy.sum(self.invar["parip"])
        
        # on ajoute les rayonnements des autres fspm de la scène
        if organs_results is not None and not organs_results.empty:
            if isinstance(self.global_index, list):
                species_not_legume = [i for i in organs_results["VegetationType"].unique() if i not in self.global_index]
            else:
                species_not_legume = [i for i in organs_results["VegetationType"].unique() if i != self.global_index]
            filtered_data = organs_results[(organs_results.VegetationType.isin(species_not_legume))]
            if not filtered_data.empty:
                pari_canopy += numpy.sum(filtered_data["par Ei"]) * energy

        ratio_pari_plante = self.invar["parip"] / (pari_canopy + 10e-15)
        self.epsi = pari_soil * ratio_pari_plante

        print(self.simulation_name, "epsi = ", sum(self.epsi))

    def compute_potential_plant_growth(self):
        outvar = self.lsystem.tag_loop_inputs[1]
        ParamP = self.lsystem.tag_loop_inputs[3]
        meteo_j = self.lsystem.tag_loop_inputs[6]
        mng_j = self.lsystem.tag_loop_inputs[7]
        nbplantes = self.lsystem.tag_loop_inputs[11]
        ls_ftswStress = self.lsystem.tag_loop_inputs[24]
        ls_NNIStress = self.lsystem.tag_loop_inputs[25]
        ls_TStress = self.lsystem.tag_loop_inputs[26]
        lsApex = self.lsystem.tag_loop_inputs[27]
        lsApexAll = self.lsystem.tag_loop_inputs[28]
        opt_stressW = self.lsystem.tag_loop_inputs[38]
        opt_stressN = self.lsystem.tag_loop_inputs[39]
        opt_stressGel = self.lsystem.tag_loop_inputs[40]

        if self.domain is None:
            surfsolref = self.lsystem.tag_loop_inputs[12]
        else:
            surfsolref = (self.domain[1][0] - self.domain[0][0]) * (self.domain[1][1] - self.domain[0][1])

        self.demandeN = []
        self.temps = []
        self.invar, outvar, ls_demandeN_bis, temps = loop.daily_growth_loop(
            ParamP,
            self.invar,
            outvar,
            self.epsi,
            meteo_j,
            mng_j,
            nbplantes,
            surfsolref,
            ls_ftswStress,
            ls_NNIStress,
            ls_TStress,
            lsApex,
            lsApexAll,
            opt_stressW,
            opt_stressN,
            opt_stressGel,
        )

        self.demandeN = ls_demandeN_bis
        self.temps = temps

    def run(self):
        # pour les variables communes
        (
            invar0,
            outvar,
            invar_sc,
            ParamP,
            station,
            carto,
            meteo_j,
            mng_j,
            DOY,
            cutNB,
            start_time,
            nbplantes,
            surfsolref0,
            m_lais,
            dicFeuilBilanR,
            surf_refVOX,
            triplets,
            ls_dif,
            S0,
            par_SN,
            lims_sol,
            ls_roots,
            ls_mat_res,
            vCC,
            ls_ftswStress,
            ls_NNIStress,
            ls_TStress,
            lsApex,
            lsApexAll,
            dicOrgans,
            deltaI_I0,
            nbI_I0,
            I_I0profilLfPlant,
            I_I0profilPetPlant,
            I_I0profilInPlant,
            NlClasses,
            NaClasses,
            NlinClasses,
            opt_stressW,
            opt_stressN,
            opt_stressGel,
            opt_residu,
            dxyz,
        ) = self.lsystem.tag_loop_inputs

        if self.domain is None:
            surfsolref = surfsolref0
        else:
            surfsolref = (self.domain[1][0] - self.domain[0][0]) * (self.domain[1][1] - self.domain[0][1])

        [soil, stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, temps_sol] = self.results_soil

        ##########
        # setp update plant stress variables
        ##########
        tag_inputs_stress = [
            ParamP,
            self.invar,
            invar_sc,
            self.temps,
            DOY,
            nbplantes,
            surfsolref,
            self.epsi,
            ls_ftsw,
            ls_transp,
            ls_Act_Nuptake_plt,
            self.demandeN,
            ls_ftswStress,
            ls_TStress,
            dicOrgans,
            dicFeuilBilanR,
            lsApex,
            start_time,
            cutNB,
            deltaI_I0,
            nbI_I0,
            I_I0profilLfPlant,
            I_I0profilPetPlant,
            I_I0profilInPlant,
            NlClasses,
            NaClasses,
            NlinClasses,
            outvar,
        ]

        (
            self.invar,
            invar_sc,
            outvar,
            I_I0profilInPlant,
            ls_ftswStress,
            ls_NNIStress,
            ls_TStress,
        ) = loop.Update_stress_loop(*tag_inputs_stress)

        ##########
        # step update soil residues senescence
        ##########

        # refait initialisation des residues au step 1 avec ensemble des plante (ParamP commun)
        current_iter = self.doy() - self.lsystem.DOYdeb
        if current_iter == 0 and opt_residu == 1:
            CC = init_plant_residues_fromParamP(soil, opt_residu, ParamP, par_SN)

        if opt_residu == 1:  # option residu activee: mise a jour des cres
            invar_merge = self.invar

            tag_inputs_residue_updt = [
                ls_mat_res,
                soil,
                ls_roots,
                par_SN["PROFHUMs"],
                ParamP,
                invar_merge,
                opt_stressGel,
            ]
            ls_mat_res = loop.distrib_residue_mat_frominvar(
                *tag_inputs_residue_updt
            )  # update la matrice des residus (propre a l-egume/VGL)
            soil = loop.merge_residue_mat(ls_mat_res, vCC, soil)  # update du sol

        #########
        # reinjecte les sorties midiee dans le lsystem
        #########
        self.lsystem.invar = self.invar
        self.lsystem.outvar = outvar
        self.lsystem.invar_sc = invar_sc

        self.lsystem.S = soil
        self.lsystem.stateEV = stateEV
        self.lsystem.ls_mat_res = ls_mat_res

        self.lsystem.res_trans = self.res_trans
        if self.res_abs_i is not None:
            self.lsystem.res_abs_i = numpy.array([self.res_abs_i])
        self.lsystem.res_rfr = self.rfr

        self.lsystem.ls_ftswStress = ls_ftswStress
        self.lsystem.ls_NNIStress = ls_NNIStress
        self.lsystem.ls_TStress = ls_TStress
        self.lsystem.I_I0profilInPlant = I_I0profilInPlant

    def end(self):
        if self.lsystem.tag_loop_inputs[8] < self.lsystem.DOYend :
            #fermeture des bilans sol -> dicout
            dicout = loop.sol_dicout_endsim(self.lsystem.S, self.lsystem.outvar, self.lsystem.DOYdeb, self.lsystem.DOYend, self.lsystem.opt_residu)
            
            #rasemble noms et objets pour ecritures sorties
            ls_outf_names = [self.lsystem.outvarfile, self.lsystem.outBilanNfile, self.lsystem.outHRfile, self.lsystem.resrootfile, self.lsystem.lsorgfile, self.lsystem.outMngfile, self.lsystem.outsdfile] #noms des fichiers de sorties
            ls_objw = [self.lsystem.outvar, dicout, self.lsystem.out_HR, self.lsystem.res_root, self.lsystem.savelsOrgans, self.lsystem.mng, self.lsystem.res_sd] #objets de donnees
            
            # liste de cle a verifier pour ecriture des sorties journaliere
            ls_keyvar_pot = ['colnames','pattern','TT','time','cutNB','SurfPlante', 'PARaPlante', 'PARiPlante', 'epsi', 'dMSaer', 'Hplante', 'Dplante','RLTot','RDepth','MS_aerien','MS_feuil','MS_tot','countSh','countShExp','demandC','Leaf_Stem','NBsh','NBI','NBD1','NBB','FTSW','Etransp','DemandN_Feuil','DemandN_Pet', 'DemandN_Stem','DemandN_Tot', 'DemandN_Tot_Aer', 'Npc', 'Npc_aer', 'NNI','Ndfa', 'Qfix','Naerien','Nuptake_sol','R_DemandC_Root', 'SRL','dMSenFeuil','dMSenTige', 'MS_pivot', 'MS_rac_fine','R_DemandC_Shoot','RUEpot','RUE','Npc_piv','Npc_rac_fine','dRLenSentot','dMSenRoot','RLTotNet','MS_rac_fineNet','perteN_rac_fine','NBphyto','NBapexAct','transpi','cumtranspi','aliveB','dMSmortGel','dNmortGel','TTphyllo','DemCp','dTT','Udev','Udevstress','TTudev','MS_aerienNonRec', 'MS_aerienRec', 'NaerienNonRec','NaerienRec','Ncoty','MS_tige','graineC','graineN','CreservPiv','NreservPiv','dMSenPiv','dMSenNonRec', 'perteN_NonRec', 'perteN_Piv', 'perteN_aerien', 'Npc_aerNonRec','MS_senaerien','dMSmortPlant_aer','dMSmortPlant_pivot','dMSmortPlant_racfine','dNmortPlant_aer','dNmortPlant_pivot','dNmortPlant_racfine','alivePiv','alive','ChangeRoot','RLentotfromRootMass','RLentotfromDev','ConcNmoy']
            
            # ecriture + liste des fichiers ecrits en sortie et affichee en fin de simul
            ls_fileOUT = loop.write_vgl_outf(self.lsystem.outf, self.lsystem.path_out, ls_outf_names, ls_objw, ls_keyvar_pot, self.lsystem.outfvar)
    
    
            #zippe les sorties si option opt_zip
            if self.lsystem.opt_zip==1:
                if self.lsystem.opt_verbose == 1:
                    print(ls_fileOUT)
        
            nomzip = self.lsystem.outvarfile[0:-4]
            IOtable.Outzip(self.lsystem.path_out, nomzip+'.zip', ls_fileOUT)
            IOtable.Outdel(ls_fileOUT)

        print(("".join((self.simulation_name, " - done"))))
        # désallocation des lsystem
        self.lsystem.clear()

    def voxels_size(self):
        return self.lsystem.tag_loop_inputs[-1]

    def number_of_voxels(self):
        leaf_area_per_voxels = self.lsystem.tag_loop_inputs[13]
        return leaf_area_per_voxels.shape

    def transfer_ratp_legume(self, energy, voxels_outputs, nb0, epsilon=1e-8):
        m_lais = self.lsystem.tag_loop_inputs[13]
        # initialize absorb energy array
        res_abs_i = numpy.zeros((m_lais.shape[0], m_lais.shape[1], m_lais.shape[2], m_lais.shape[3]))

        ratpzlayers = max(voxels_outputs["Nz"])

        # voxel top side area
        dS = self.lsystem.tag_loop_inputs[15]
        res_trans = numpy.ones((m_lais.shape[1], m_lais.shape[2], m_lais.shape[3]))
        # maximum transmitted energy is total incoming energy per area
        res_trans = res_trans * (energy * dS)

        for ix in range(m_lais.shape[3]):
            for iy in range(m_lais.shape[2]):
                for iz in range(ratpzlayers):
                    legume_iz = iz + nb0

                    condition_x = voxels_outputs.Nx == m_lais.shape[2] - iy
                    vox_data = voxels_outputs[
                        condition_x & (voxels_outputs.Ny == ix + 1) & (voxels_outputs.Nz == iz + 1)
                    ]
                    if not vox_data.empty:
                        a = min(sum(vox_data["Transmitted"]), dS)
                        res_trans[legume_iz, iy, ix] = energy * a

                    s_entity = 0
                    for k in range(m_lais.shape[0]):
                        s_entity += m_lais[k][legume_iz][iy][ix]

                    if s_entity > 0.0:
                        if isinstance(self.global_index, list):
                            indices = self.global_index
                        elif isinstance(self.global_index, int):
                            indices = [self.global_index]
                        for id_legume,id_global in enumerate(indices):
                            if len(vox_data) > 0:
                                v_dat = vox_data[vox_data.VegetationType == id_global + 1]
                                v = v_dat["Intercepted"].values[0]
                                if v > epsilon:
                                    res_abs_i[id_legume, legume_iz, iy, ix] = energy * v

                                # if a voxel has leaf area > 0, it must have a minimum intercepted energy value
                                else:
                                    res_abs_i[id_legume, legume_iz, iy, ix] = epsilon

        return res_abs_i, res_trans

    def transfer_caribu_legume(
        self,
        energy,
        nb0,
        elements_outputs,
        sensors_outputs,
        epsilon=1e-8,
    ):
        # initialize absorb energy
        nplantes = len(self.invar["Hplante"])
        self.invar["parap"] = scipy.array([0.0] * nplantes)
        self.invar["parip"] = scipy.array([0.0] * nplantes)

        ent_organs_outputs = pandas.DataFrame({})
        if isinstance(self.global_index, list):
            filter = elements_outputs["VegetationType"].isin(self.global_index)
        else:
            filter = elements_outputs.VegetationType == self.global_index
        ent_organs_outputs = elements_outputs[filter]

        # non empty scene
        for i in range(len(ent_organs_outputs)):
            organe_id = int(ent_organs_outputs.iloc[i]["Organ"])

            # PAR in W/m²
            par_intercept = ent_organs_outputs.iloc[i]["par Ei"] * energy
            S_leaf = ent_organs_outputs.iloc[i]["Area"]

            id_plante = self.lstring[organe_id][0]
            p_s = par_intercept * S_leaf
            a = float(self.invar["parip"][id_plante])
            self.invar["parip"][id_plante] = a + p_s

            # we remove senescent leaves
            if self.lstring[organe_id][9] != "sen":
                a = float(self.invar["parap"][id_plante])
                self.invar["parap"][id_plante] = a + p_s

        # all non empty plant must have a minimum intercepted energy
        plants_surface = self.lsystem.tag_loop_inputs[14]["surf"]
        if plants_surface != []:
            if len(self.invar["parip"]) == len(plants_surface):
                for p in range(len(self.invar["parip"])):
                    if self.invar["parip"][p] == 0.0 and plants_surface[p] > 0.0:
                        self.invar["parip"][p] = epsilon

        # conversion
        c = (3600 * 24) / 1000000
        self.invar["parap"] *= c
        self.invar["parip"] *= c

        ## Transmitted radiations throughout a grid of voxels
        m_lais = self.lsystem.tag_loop_inputs[13]
        res_trans = numpy.ones((m_lais.shape[1], m_lais.shape[2], m_lais.shape[3]))

        # if non empty scene
        if not elements_outputs.empty:
            if [x for x in self.indexer.legume_names if x != self.name] != []:
                sensors_specy_id = self.global_index
            else:
                sensors_specy_id = 0
            sensors = sensors_outputs[sensors_outputs.VegetationType==sensors_specy_id]
            ID_capt = 0
            for ix in range(m_lais.shape[3]):
                for iy in range(m_lais.shape[2]):
                    for iz in range(m_lais.shape[1] - nb0):
                        a = min(sensors.iloc[ID_capt]["PAR"], 1.0)
                        res_trans[((m_lais.shape[1] - 1)) - iz][iy][ix] = a
                        ID_capt += 1

        # surface d'une face d'un voxel
        dS = self.lsystem.tag_loop_inputs[15]
        # gives maximum transmitted energy
        res_trans = res_trans * energy * dS

        return res_trans

    def energy(self):
        meteo_j = self.lsystem.tag_loop_inputs[6]
        energy = 0.48 * meteo_j["RG"] * 10000 / (3600 * 24)
        return energy

    def doy(self):
        DOY = self.lsystem.tag_loop_inputs[8]
        return DOY

    def set_domain(self, domain):
        self.domain = domain

    @staticmethod
    def fake_scene():
        epsilon = 1e-14
        return {0: [[(0.0, 0.0, 0.0), (0.0, epsilon, 0.0), (0.0, epsilon, epsilon)]]}


def passive_lighting(data, energy, DOY, scene, legume_wrapper, lighting_wrapper):
    invar_saved = deepcopy(legume_wrapper.invar)
    lighting_wrapper.run(scenes=[scene], energy=energy, day=DOY, parunit="RG")

    legume_wrapper.light_results(energy, lighting_wrapper)
    legume_wrapper.invar = deepcopy(invar_saved)

    data["epsi"].extend(legume_wrapper.epsi)
    data["parip"].extend(legume_wrapper.invar["parip"])
    data["t"].extend([DOY] * len(legume_wrapper.epsi))
