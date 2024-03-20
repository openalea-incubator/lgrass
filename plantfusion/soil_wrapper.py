import numpy
import pandas
import os
import math

from soil3ds import soil_moduleN as solN
import soil3ds.IOxls as IOxls
import soil3ds.IOtable as IOtable
import legume.initialisation as initial

from plantfusion.planter import Planter
from plantfusion.utils import create_child_folder


class Soil_wrapper(object):
    def __init__(
        self,
        in_folder="",
        out_folder="",
        IDusm=1,
        nameconfigfile="liste_usms_exemple.xls",
        ongletconfigfile="exemple",
        opt_residu=0,
        opt_Nuptake=1,
        planter=Planter(),
        legume_pattern=False,
        legume_wrapper=None,
        only_water_balance=False,
        save_results=False,
    ) -> None:
        self.only_water_balance = only_water_balance

        if not legume_pattern:  
            self.planter = planter

            # lecture des parametre sol directement a partir d'un fichier sol
            # initialisation taille scene / discretisation (1D - homogene pour ttes les couches)
            # initialisation d'une première scène pour avoir les côtés du domaine simulé

            # lecture meteo / mng journaliere = fixe
            usms_path = os.path.join(in_folder, nameconfigfile)
            usms = IOxls.xlrd.open_workbook(usms_path)
            ls_usms = IOtable.conv_dataframe(IOxls.get_xls_col(usms.sheet_by_name(ongletconfigfile)))
            id = ls_usms["ID_usm"].index(float(IDusm))

            fxls_sol = ls_usms["sol"][id]
            ongletS = ls_usms["ongletS"][id]
            path_sol = os.path.join(in_folder, fxls_sol)
            par_SN, par_sol = IOxls.read_sol_param(path_sol, ongletS)

            fxls_mng = ls_usms["mng"][id]
            ongletMn = ls_usms["ongletMn"][id]
            fxls_Met = ls_usms["meteo"][id]
            ongletMet = ls_usms["ongletM"][id]
            path_met = os.path.join(in_folder, fxls_Met)
            path_mng = os.path.join(in_folder, fxls_mng)
            fxls_ini = ls_usms["inis"][id]
            ongletIn = ls_usms["ongletIn"][id]
            path_inis = os.path.join(in_folder, fxls_ini)
            inis = IOxls.read_plant_param(path_inis, ongletIn)
            met = IOxls.read_met_file(path_met, ongletMet)
            mng = IOxls.read_met_file(path_mng, ongletMn)

            dz_sol = inis["dz_sol"]
            pattern8 = [[v * 100 for v in x] for x in planter.domain]  # conversion m en cm
            discret_solXY = list(map(int, inis["discret_solXY"]))  # [10,10]# nb de discretisation du sol en X et en Y

            # meteo / mng journalier
            meteo_j = IOxls.extract_dataframe(
                met, ["TmoyDay", "RG", "Et0", "Precip", "Tmin", "Tmax", "Tsol"], "DOY", val=ls_usms["DOYdeb"][id]
            )
            for k in list(meteo_j.keys()):
                meteo_j[k] = meteo_j[k][0]

            # initialisation du sol avec fonction pour Lpy (sans residus)
            soil, Tsol = initial.init_sol_fromLpy(
                inis, meteo_j, par_sol, par_SN, discret_solXY, dz_sol, pattern8, opt_residu=0
            )

            self.parameters_SN = par_SN
            self.meteo = met
            self.management = mng
            self.soil = soil
            self.option_residu = opt_residu
            self.option_Nuptake = opt_Nuptake

            plant_path = os.path.join(in_folder, ls_usms["plante"][id])
            self.soil_plants_parameters_bare_soil = IOxls.read_plant_param(plant_path, "solnu")

            if legume_wrapper is not None:
                legume_wrapper.lsystem.S = soil

        else:
            if only_water_balance :
                raise AttributeError("l-egume soil computes N balance")
            self.soil = legume_wrapper.lsystem.tag_loop_inputs[18]
            self.option_residu = legume_wrapper.lsystem.tag_loop_inputs[-2]
            self.option_Nuptake = legume_wrapper.lsystem.opt_Nuptake
            self.meteo = legume_wrapper.lsystem.meteo
            self.management = legume_wrapper.lsystem.mng
            self.parameters_SN = legume_wrapper.lsystem.par_SN
            self.soil_plants_parameters_bare_soil = IOxls.read_plant_param(legume_wrapper.lsystem.path_plante, "solnu")

        self.soil_dimensions = [len(self.soil.dxyz[i]) for i in [2, 0, 1]]

        self.save_results = save_results
        if save_results:
            create_child_folder(out_folder, "soil")
            self.out_folder = os.path.join(out_folder, "soil")
            self.data: dict = {
                # clé ajouté
                "DOY": [],
                "kgNO3solHa": [],
                "kgNH4solHa": [],
                "kgNsolHa": [],
                "uptNO3PltHa": [],
                "lix": [],
                "cumMinN_j": [],
                "Lix_j": [],
                "cumUptakePlt_j": [],
                "azomes": [],
                "transp": [],
                "tsw": [],
                # clé dans soil_moduleN.CloseNbalance (soil3ds)
                "FinalInertN": [],
                "FinalActiveN": [],
                "finalNZygo": [],
                "finalNres": [],
                "ResidueMinNtot": [],
                "NminfromNresCum": [],
                "HumusMinNtot": [],
                "MinNtot": [],
                "InputNtot": [],
                "OutputNtot": [],
                "TotNRain": [],
                "TotNIrrig": [],
                "TotFertNO3": [],
                "TotFertNH4": [],
                "InputNmintot": [],
                "FinalNO3": [],
                "FinalNH4": [],
                "Lixtot": [],
                "N2Otot": [],
                "TotUptPlt": [],
                "OutputNmintot": [],
            }

    def run(
        self,
        day=1,
        N_content_roots_per_plant=[],
        roots_length_per_plant_per_soil_layer=[],
        soil_plants_parameters=[],
        plants_light_interception=[],
    ):
        meteo_j = IOxls.extract_dataframe(
            self.meteo, ["TmoyDay", "RG", "Et0", "Precip", "Tmin", "Tmax", "Tsol"], "DOY", val=day
        )
        mng_j = IOxls.extract_dataframe(
            self.management, ["Coupe", "Irrig", "FertNO3", "FertNH4", "Hcut"], "DOY", val=day
        )
        for k in list(meteo_j.keys()):
            meteo_j[k] = meteo_j[k][0]
        for k in list(mng_j.keys()):
            mng_j[k] = mng_j[k][0]

        if not self.only_water_balance :
            self.inputs = []
            self.inputs.append(self.soil)
            self.inputs.append(self.parameters_SN)
            self.inputs.append(meteo_j)
            self.inputs.append(mng_j)
            v = []
            for t in soil_plants_parameters:
                v.extend(t)
            self.inputs.append(v)
            v = []
            for t in plants_light_interception:
                v.extend(t)
            self.inputs.append(v)
            v = []
            for t in roots_length_per_plant_per_soil_layer:
                v.extend(t)
            self.inputs.append(v)
            v = []
            for t in N_content_roots_per_plant:
                v.extend(t)
            self.inputs.append(numpy.array(v))
            self.inputs.append(self.option_residu)
            self.inputs.append(self.option_Nuptake)

            self.nb_plants = len(self.inputs[4])
            self.results = solN.step_bilanWN_solVGL(*self.inputs)

            if self.save_results:
                self.append_results(day)
        else:
            self.inputs = []
            self.inputs.append(meteo_j['Et0'] * self.soil.surfsolref)

            v = []
            for t in roots_length_per_plant_per_soil_layer:
                v.extend(t)
            self.inputs.append(v)
            v = []
            for t in plants_light_interception:
                v.extend(t)
            self.inputs.append(v)

            self.inputs.append(meteo_j['Precip'] * self.soil.surfsolref)
            self.inputs.append(mng_j['Irrig'] * self.soil.surfsolref)
            self.inputs.append(self.soil.stateEV)
            self.nb_plants = len(self.inputs[1])
            self.results = self.soil.stepWBmc(*self.inputs)

    def bare_soil_inputs(self, epsilon=1e-10):
        R1 = self.soil.m_1 * epsilon  # pas zero sinon bug FTSW
        ls_roots = [R1]
        ls_epsi = [0.0]
        ls_N = [1.0]
        ls_paramP = [self.soil_plants_parameters_bare_soil]

        return ls_N, ls_roots, ls_paramP, ls_epsi

    def append_results(self, day):
        Nuptake_1plant = self.results[4][0].sum()
        soil = self.inputs[0]

        # calcul le bilan
        soil.CloseNbalance()

        # append les valeurs copiées dans le bilan soil3ds
        for key in soil.bilanN.keys():
            if key in self.data:
                if isinstance(soil.bilanN[key], list):
                    self.data[key].append(soil.bilanN[key][-1])
                else:
                    self.data[key].append(soil.bilanN[key])

        # calcul des valeurs custom
        kgNO3solHa = soil.m_NO3.sum() / soil.surfsolref * 10000
        kgNH4solHa = soil.m_NH4.sum() / soil.surfsolref * 10000
        kgNsolHa = kgNO3solHa + kgNH4solHa
        uptNO3PltHa = Nuptake_1plant * self.nb_plants / soil.surfsolref * 10000
        lix = soil.lixiNO3

        cumMinN_j = soil.bilanN["cumMinN"][-1]
        Lix_j = soil.bilanN["cumLix"][-1]
        cumUptakePlt_j = soil.bilanN["cumUptakePlt"][-1].sum()

        transp = sum(self.results[3])
        tsw = soil.tsw_t.sum()

        # append ces valeurs
        self.data["DOY"].append(day)
        self.data["kgNO3solHa"].append(kgNO3solHa)
        self.data["kgNH4solHa"].append(kgNH4solHa)
        self.data["kgNsolHa"].append(kgNsolHa)
        self.data["uptNO3PltHa"].append(uptNO3PltHa)
        self.data["lix"].append(lix)
        self.data["cumMinN_j"].append(cumMinN_j)
        self.data["Lix_j"].append(Lix_j)
        self.data["cumUptakePlt_j"].append(cumUptakePlt_j)
        self.data["transp"].append(transp)
        self.data["tsw"].append(tsw)

    def end(self):
        try:
            pandas.DataFrame(self.data).to_csv(os.path.join(self.out_folder, "outputs_soil3ds.csv"))
        except AttributeError:
            print("Soil save results not activated")
        finally:
            pass

    def whichvoxel_xy(self, p):
        """
        dxyz : [dx, dy, dz]
        """
        ix = math.floor((p[0] - self.soil.origin[0]) / self.soil.dxyz[0][0])
        iy = math.floor((p[1] - self.soil.origin[1]) / self.soil.dxyz[1][0])

        return ix, iy
