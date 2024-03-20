import os
import numpy
import pandas
import math

from alinea.adel.adel_dynamic import AdelDyn
from alinea.adel.echap_leaf import echap_leaves

from lightvegemanager.stems import extract_stems_from_MTG

from fspmwheat import caribu_facade
from fspmwheat import cnwheat_facade
from fspmwheat import elongwheat_facade
from fspmwheat import farquharwheat_facade
from fspmwheat import growthwheat_facade
from fspmwheat import senescwheat_facade
from fspmwheat import fspmwheat_facade

from cnwheat import (
    simulation as cnwheat_simulation,
    model as cnwheat_model,
    parameters as cnwheat_parameters,
    tools as cnwheat_tools,
)

from soil3ds.IOxls import read_plant_param

from plantfusion.utils import create_child_folder, save_df_to_csv
from plantfusion.planter import Planter
from plantfusion.indexer import Indexer


class Wheat_wrapper(object):
    HOUR_TO_SECOND_CONVERSION_FACTOR = 3600
    AXES_INDEX_COLUMNS = ["t", "plant", "axis"]
    ELEMENTS_INDEX_COLUMNS = ["t", "plant", "axis", "metamer", "organ", "element"]
    HIDDENZONES_INDEX_COLUMNS = ["t", "plant", "axis", "metamer"]
    ORGANS_INDEX_COLUMNS = ["t", "plant", "axis", "organ"]
    SOILS_INDEX_COLUMNS = ["t", "plant", "axis"]

    def __init__(
        self,
        name="wheat",
        in_folder="",
        out_folder=None,
        N_fertilizations={},
        tillers_replications={},
        planter=Planter(),
        indexer=Indexer(),
        run_from_outputs=False,
        external_soil_model=False,
        nitrates_uptake_forced=False,
        initialize_nitrates_uptake=0.25,
        update_parameters_all_models=None,
        stored_times=None,
        option_static=False,
        LIGHT_TIMESTEP=4,
        SENESCWHEAT_TIMESTEP=1,
        FARQUHARWHEAT_TIMESTEP=1,
        ELONGWHEAT_TIMESTEP=1,
        GROWTHWHEAT_TIMESTEP=1,
        CNWHEAT_TIMESTEP=1,
        AXES_INITIAL_STATE_FILENAME="axes_initial_state.csv",
        ORGANS_INITIAL_STATE_FILENAME="organs_initial_state.csv",
        HIDDENZONES_INITIAL_STATE_FILENAME="hiddenzones_initial_state.csv",
        ELEMENTS_INITIAL_STATE_FILENAME="elements_initial_state.csv",
        SOILS_INITIAL_STATE_FILENAME="soils_initial_state.csv",
        METEO_FILENAME="meteo_Ljutovac2002.csv",
        NITRATES_UPTAKE_FORCINGS_FILENAME="nitrates_uptake_forcings.csv",
        SOIL_PARAMETERS_FILENAME="",
        SOIL_PARAMETERS_SHEETNAME="cnwheat",
        AXES_OUTPUTS_FILENAME="axes_outputs.csv",
        ORGANS_OUTPUTS_FILENAME="organs_outputs.csv",
        HIDDENZONES_OUTPUTS_FILENAME="hiddenzones_outputs.csv",
        ELEMENTS_OUTPUTS_FILENAME="elements_outputs.csv",
        SOILS_OUTPUTS_FILENAME="soils_outputs.csv",
        AXES_POSTPROCESSING_FILENAME="axes_postprocessing.csv",
        ORGANS_POSTPROCESSING_FILENAME="organs_postprocessing.csv",
        HIDDENZONES_POSTPROCESSING_FILENAME="hiddenzones_postprocessing.csv",
        ELEMENTS_POSTPROCESSING_FILENAME="elements_postprocessing.csv",
        SOILS_POSTPROCESSING_FILENAME="soils_postprocessing.csv",
    ) -> None:
        self.N_fertilizations = N_fertilizations
        self.tillers_replications = tillers_replications

        self.external_soil_model = external_soil_model
        self.nitrates_uptake_forced = nitrates_uptake_forced
        self.option_static = option_static

        self.last_year_doy = 0

        self.name = name
        self.indexer = indexer
        self.global_index = indexer.global_order.index(name)
        self.wheat_index = indexer.wheat_names.index(name)

        self.nb_plants = planter.number_of_plants[self.global_index]
        if name in planter.plant_density:
            self.plant_density = {1: planter.plant_density[name]}
        else:
            self.plant_density = planter.plant_density
        self.generation_type = planter.generation_type

        self.LIGHT_TIMESTEP = LIGHT_TIMESTEP
        self.SENESCWHEAT_TIMESTEP = SENESCWHEAT_TIMESTEP
        self.FARQUHARWHEAT_TIMESTEP = FARQUHARWHEAT_TIMESTEP
        self.ELONGWHEAT_TIMESTEP = ELONGWHEAT_TIMESTEP
        self.GROWTHWHEAT_TIMESTEP = GROWTHWHEAT_TIMESTEP
        self.CNWHEAT_TIMESTEP = CNWHEAT_TIMESTEP

        self.AXES_OUTPUTS_FILENAME = AXES_OUTPUTS_FILENAME
        self.ORGANS_OUTPUTS_FILENAME = ORGANS_OUTPUTS_FILENAME
        self.HIDDENZONES_OUTPUTS_FILENAME = HIDDENZONES_OUTPUTS_FILENAME
        self.ELEMENTS_OUTPUTS_FILENAME = ELEMENTS_OUTPUTS_FILENAME
        self.SOILS_OUTPUTS_FILENAME = SOILS_OUTPUTS_FILENAME

        in_folder = os.path.normpath(in_folder)

        self.AXES_POSTPROCESSING_FILENAME = AXES_POSTPROCESSING_FILENAME
        self.ORGANS_POSTPROCESSING_FILENAME = ORGANS_POSTPROCESSING_FILENAME
        self.HIDDENZONES_POSTPROCESSING_FILENAME = HIDDENZONES_POSTPROCESSING_FILENAME
        self.ELEMENTS_POSTPROCESSING_FILENAME = ELEMENTS_POSTPROCESSING_FILENAME
        self.SOILS_POSTPROCESSING_FILENAME = SOILS_POSTPROCESSING_FILENAME

        if out_folder is not None:
            self.out_folder = os.path.join(os.path.normpath(out_folder), name)
            try:
                os.mkdir(os.path.normpath(self.out_folder))
                print("Directory ", self.out_folder, " Created ")
            except FileExistsError:
                pass

            create_child_folder(self.out_folder, "brut")
            create_child_folder(self.out_folder, "postprocessing")
            create_child_folder(self.out_folder, "graphs")

        if external_soil_model:
            SOILS_INITIAL_STATE_FILENAME = None

            if nitrates_uptake_forced:
                nitrates_uptake_data_filepath = os.path.join(in_folder, NITRATES_UPTAKE_FORCINGS_FILENAME)
                nitrates_uptake_data_df = pandas.read_csv(nitrates_uptake_data_filepath)
                self.nitrates_uptake_data_grouped = nitrates_uptake_data_df.groupby(
                    cnwheat_simulation.Simulation.ORGANS_T_INDEXES
                )

            else:
                wheat_paramp = read_plant_param(os.path.normpath(SOIL_PARAMETERS_FILENAME), SOIL_PARAMETERS_SHEETNAME)
                self.soil_parameters = [wheat_paramp] * self.nb_plants

        self.inputs_dataframes = {}
        new_start_time = -1
        if run_from_outputs:
            previous_outputs_dataframes = {}

            for initial_state_filename, outputs_filename, index_columns in (
                (AXES_INITIAL_STATE_FILENAME, AXES_OUTPUTS_FILENAME, self.AXES_INDEX_COLUMNS),
                (ORGANS_INITIAL_STATE_FILENAME, ORGANS_OUTPUTS_FILENAME, self.ORGANS_INDEX_COLUMNS),
                (HIDDENZONES_INITIAL_STATE_FILENAME, HIDDENZONES_OUTPUTS_FILENAME, self.HIDDENZONES_INDEX_COLUMNS),
                (ELEMENTS_INITIAL_STATE_FILENAME, ELEMENTS_OUTPUTS_FILENAME, self.ELEMENTS_INDEX_COLUMNS),
                (SOILS_INITIAL_STATE_FILENAME, SOILS_OUTPUTS_FILENAME, self.SOILS_INDEX_COLUMNS),
            ):
                if not initial_state_filename:
                    continue
                previous_outputs_dataframe = pandas.read_csv(os.path.join(out_folder, outputs_filename))
                # Convert NaN to None
                previous_outputs_dataframes[outputs_filename] = previous_outputs_dataframe.replace({numpy.nan: None})

                assert "t" in previous_outputs_dataframes[outputs_filename].columns

                last_t_step = max(previous_outputs_dataframes[outputs_filename]["t"])
                new_start_time = last_t_step + 1

                if initial_state_filename == ELEMENTS_INITIAL_STATE_FILENAME:
                    elements_previous_outputs = previous_outputs_dataframes[outputs_filename]
                    new_initial_state = elements_previous_outputs[~elements_previous_outputs.is_over.isnull()]
                else:
                    new_initial_state = previous_outputs_dataframes[outputs_filename]
                idx = (
                    new_initial_state.groupby([col for col in index_columns if col != "t"])["t"].transform(max)
                    == new_initial_state["t"]
                )
                self.inputs_dataframes[initial_state_filename] = new_initial_state[idx].drop(["t"], axis=1)

            # Make sure boolean columns have either type bool or float
            bool_columns = [
                "is_over",
                "is_growing",
                "leaf_is_emerged",
                "internode_is_visible",
                "leaf_is_growing",
                "internode_is_growing",
                "leaf_is_remobilizing",
                "internode_is_remobilizing",
            ]
            for df in [
                self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME],
                self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME],
            ]:
                for cln in bool_columns:
                    if cln in df.keys():
                        df[cln].replace(to_replace="False", value=0.0, inplace=True)
                        df[cln].replace(to_replace="True", value=1.0, inplace=True)
                        df[cln] = pandas.to_numeric(df[cln])
        else:
            for inputs_filename in (
                AXES_INITIAL_STATE_FILENAME,
                ORGANS_INITIAL_STATE_FILENAME,
                HIDDENZONES_INITIAL_STATE_FILENAME,
                ELEMENTS_INITIAL_STATE_FILENAME,
                SOILS_INITIAL_STATE_FILENAME,
            ):
                if not inputs_filename:
                    continue
                inputs_dataframe = pandas.read_csv(os.path.join(in_folder, inputs_filename))
                self.inputs_dataframes[inputs_filename] = inputs_dataframe.replace({numpy.nan: None})

        # Start time of the simulation
        self.start_time = max(0, new_start_time)

        # Name of the CSV files which contains the meteo data
        self.meteo = pandas.read_csv(os.path.join(in_folder, METEO_FILENAME), index_col="t")

        # -- OUTPUTS CONFIGURATION --

        # Save the outputs with a full scan of the MTG at each time step (or at selected time steps)
        UPDATE_SHARED_DF = False
        if stored_times is None:
            stored_times = "all"
        if not (stored_times == "all" or isinstance(stored_times, list)):
            print("stored_times should be either 'all', a list or an empty list.")
            raise
        self.stored_times = stored_times

        # create empty dataframes to shared data between the models
        self.shared_axes_inputs_outputs_df = pandas.DataFrame()
        self.shared_organs_inputs_outputs_df = pandas.DataFrame()
        self.shared_hiddenzones_inputs_outputs_df = pandas.DataFrame()
        self.shared_elements_inputs_outputs_df = pandas.DataFrame()
        self.shared_soils_inputs_outputs_df = pandas.DataFrame()

        # define lists of dataframes to store the inputs and the outputs of the models at each step.
        self.axes_all_data_list = []
        self.organs_all_data_list = []
        self.hiddenzones_all_data_list = []
        self.elements_all_data_list = []
        self.soils_all_data_list = []

        self.all_simulation_steps = []

        # -- ADEL and MTG CONFIGURATION --

        # read adelwheat inputs at t0
        self.adel_wheat = AdelDyn(seed=1, scene_unit="m", leaves=echap_leaves(xy_model="Soissons_byleafclass"))
        self.g = self.adel_wheat.load(dir=in_folder)

        # ---------------------------------------------
        # ----- CONFIGURATION OF THE wrapperS -------
        # ---------------------------------------------

        # -- ELONGWHEAT (created first because it is the only wrapper to add new metamers) --
        # Initial states
        elongwheat_hiddenzones_initial_state = self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
            elongwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS
            + [
                i
                for i in elongwheat_facade.simulation.HIDDENZONE_INPUTS
                if i in self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()
        elongwheat_elements_initial_state = self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            elongwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS
            + [
                i
                for i in elongwheat_facade.simulation.ELEMENT_INPUTS
                if i in self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()
        elongwheat_axes_initial_state = self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            elongwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS
            + [
                i
                for i in elongwheat_facade.simulation.AXIS_INPUTS
                if i in self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        phytoT_df = pandas.read_csv(os.path.join(in_folder, "phytoT.csv"))

        # update parameters if specified
        if update_parameters_all_models and "elongwheat" in update_parameters_all_models:
            update_parameters_elongwheat = update_parameters_all_models["elongwheat"]
        else:
            update_parameters_elongwheat = None

        # wrapper initialisation
        self.elongwheat_facade_ = elongwheat_facade.ElongWheatFacade(
            self.g,
            ELONGWHEAT_TIMESTEP * self.HOUR_TO_SECOND_CONVERSION_FACTOR,
            elongwheat_axes_initial_state,
            elongwheat_hiddenzones_initial_state,
            elongwheat_elements_initial_state,
            self.shared_axes_inputs_outputs_df,
            self.shared_hiddenzones_inputs_outputs_df,
            self.shared_elements_inputs_outputs_df,
            self.adel_wheat,
            phytoT_df,
            update_parameters_elongwheat,
            update_shared_df=UPDATE_SHARED_DF,
        )

        # -- CARIBU --
        self.caribu_facade_ = caribu_facade.CaribuFacade(
            self.g, self.shared_elements_inputs_outputs_df, self.adel_wheat, update_shared_df=UPDATE_SHARED_DF
        )

        # -- SENESCWHEAT --
        # Initial states
        senescwheat_roots_initial_state = (
            self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]
            .loc[self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]["organ"] == "roots"][
                senescwheat_facade.converter.ROOTS_TOPOLOGY_COLUMNS
                + [
                    i
                    for i in senescwheat_facade.converter.SENESCWHEAT_ROOTS_INPUTS
                    if i in self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns
                ]
            ]
            .copy()
        )

        senescwheat_elements_initial_state = self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            senescwheat_facade.converter.ELEMENTS_TOPOLOGY_COLUMNS
            + [
                i
                for i in senescwheat_facade.converter.SENESCWHEAT_ELEMENTS_INPUTS
                if i in self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        senescwheat_axes_initial_state = self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            senescwheat_facade.converter.AXES_TOPOLOGY_COLUMNS
            + [
                i
                for i in senescwheat_facade.converter.SENESCWHEAT_AXES_INPUTS
                if i in self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        # update parameters if specified
        if update_parameters_all_models and "senescwheat" in update_parameters_all_models:
            update_parameters_senescwheat = update_parameters_all_models["senescwheat"]
        else:
            update_parameters_senescwheat = None

        # wrapper initialisation
        self.senescwheat_facade_ = senescwheat_facade.SenescWheatFacade(
            self.g,
            SENESCWHEAT_TIMESTEP * self.HOUR_TO_SECOND_CONVERSION_FACTOR,
            senescwheat_roots_initial_state,
            senescwheat_axes_initial_state,
            senescwheat_elements_initial_state,
            self.shared_organs_inputs_outputs_df,
            self.shared_axes_inputs_outputs_df,
            self.shared_elements_inputs_outputs_df,
            update_parameters_senescwheat,
            update_shared_df=UPDATE_SHARED_DF,
        )

        # -- FARQUHARWHEAT --
        # Initial states
        farquharwheat_elements_initial_state = self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            farquharwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS
            + [
                i
                for i in farquharwheat_facade.converter.FARQUHARWHEAT_ELEMENTS_INPUTS
                if i in self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        farquharwheat_axes_initial_state = self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            farquharwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS
            + [
                i
                for i in farquharwheat_facade.converter.FARQUHARWHEAT_AXES_INPUTS
                if i in self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        # Use the initial version of the photosynthesis sub-model (as in Barillot et al. 2016, and in Gauthier et al. 2020)
        update_parameters_farquharwheat = {"SurfacicProteins": False, "NSC_Retroinhibition": False}

        # wrapper initialisation
        self.farquharwheat_facade_ = farquharwheat_facade.FarquharWheatFacade(
            self.g,
            farquharwheat_elements_initial_state,
            farquharwheat_axes_initial_state,
            self.shared_elements_inputs_outputs_df,
            update_parameters_farquharwheat,
            update_shared_df=UPDATE_SHARED_DF,
        )

        # -- GROWTHWHEAT --
        # Initial states
        growthwheat_hiddenzones_initial_state = self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
            growthwheat_facade.converter.HIDDENZONE_TOPOLOGY_COLUMNS
            + [
                i
                for i in growthwheat_facade.simulation.HIDDENZONE_INPUTS
                if i in self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        growthwheat_elements_initial_state = self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            growthwheat_facade.converter.ELEMENT_TOPOLOGY_COLUMNS
            + [
                i
                for i in growthwheat_facade.simulation.ELEMENT_INPUTS
                if i in self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        growthwheat_root_initial_state = (
            self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]
            .loc[self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME]["organ"] == "roots"][
                growthwheat_facade.converter.ROOT_TOPOLOGY_COLUMNS
                + [
                    i
                    for i in growthwheat_facade.simulation.ROOT_INPUTS
                    if i in self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns
                ]
            ]
            .copy()
        )

        growthwheat_axes_initial_state = self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME][
            growthwheat_facade.converter.AXIS_TOPOLOGY_COLUMNS
            + [
                i
                for i in growthwheat_facade.simulation.AXIS_INPUTS
                if i in self.inputs_dataframes[AXES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        # update parameters if specified
        if update_parameters_all_models and "growthwheat" in update_parameters_all_models:
            update_parameters_growthwheat = update_parameters_all_models["growthwheat"]
        else:
            update_parameters_growthwheat = None

        # wrapper initialisation
        self.growthwheat_facade_ = growthwheat_facade.GrowthWheatFacade(
            self.g,
            GROWTHWHEAT_TIMESTEP * self.HOUR_TO_SECOND_CONVERSION_FACTOR,
            growthwheat_hiddenzones_initial_state,
            growthwheat_elements_initial_state,
            growthwheat_root_initial_state,
            growthwheat_axes_initial_state,
            self.shared_organs_inputs_outputs_df,
            self.shared_hiddenzones_inputs_outputs_df,
            self.shared_elements_inputs_outputs_df,
            self.shared_axes_inputs_outputs_df,
            update_parameters_growthwheat,
            update_shared_df=UPDATE_SHARED_DF,
        )

        # -- CNWHEAT --
        # Initial states
        cnwheat_organs_initial_state = self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME][
            [
                i
                for i in cnwheat_facade.cnwheat_converter.ORGANS_VARIABLES
                if i in self.inputs_dataframes[ORGANS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        cnwheat_hiddenzones_initial_state = self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME][
            [
                i
                for i in cnwheat_facade.cnwheat_converter.HIDDENZONE_VARIABLES
                if i in self.inputs_dataframes[HIDDENZONES_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        cnwheat_elements_initial_state = self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME][
            [
                i
                for i in cnwheat_facade.cnwheat_converter.ELEMENTS_VARIABLES
                if i in self.inputs_dataframes[ELEMENTS_INITIAL_STATE_FILENAME].columns
            ]
        ].copy()

        if not external_soil_model:
            cnwheat_soils_initial_state = self.inputs_dataframes[SOILS_INITIAL_STATE_FILENAME][
                [
                    i
                    for i in cnwheat_facade.cnwheat_converter.SOILS_VARIABLES
                    if i in self.inputs_dataframes[SOILS_INITIAL_STATE_FILENAME].columns
                ]
            ].copy()
        else:
            cnwheat_soils_initial_state = None

        # update parameters if specified
        if update_parameters_all_models and "cnwheat" in update_parameters_all_models:
            update_parameters_cnwheat = update_parameters_all_models["cnwheat"]
        else:
            update_parameters_cnwheat = {}

        # wrapper initialisation
        self.cnwheat_facade_ = cnwheat_facade.CNWheatFacade(
            self.g,
            CNWHEAT_TIMESTEP * self.HOUR_TO_SECOND_CONVERSION_FACTOR,
            self.plant_density,
            update_parameters_cnwheat,
            cnwheat_organs_initial_state,
            cnwheat_hiddenzones_initial_state,
            cnwheat_elements_initial_state,
            cnwheat_soils_initial_state,
            self.shared_axes_inputs_outputs_df,
            self.shared_organs_inputs_outputs_df,
            self.shared_hiddenzones_inputs_outputs_df,
            self.shared_elements_inputs_outputs_df,
            self.shared_soils_inputs_outputs_df,
            update_shared_df=UPDATE_SHARED_DF,
            external_soil_model=external_soil_model,
        )

        # Force the senescence and photosynthesis of the population
        if external_soil_model:
            if self.nitrates_uptake_forced:
                t = 0
                self.force_nitrates_uptake(t)
            else:
                self.uptake_nitrate_hour = initialize_nitrates_uptake
                self.update_Nitrates_cnwheat_mtg()

        # -- FSPMWHEAT --
        # wrapper initialisation
        self.fspmwheat_facade_ = fspmwheat_facade.FSPMWheatFacade(self.g)

        # update geometry
        self.adel_wheat.update_geometry(self.g)

    def light_inputs(self, planter):
        if self.generation_type == "default":
            scene_wheat = planter.create_heterogeneous_canopy(
                self.adel_wheat,
                mtg=self.g,
                stem_name="StemElement",
                leaf_name="LeafElement1",
                indice_wheat_instance=self.wheat_index,
            )

        elif self.generation_type == "random":
            scene_wheat = planter.generate_random_wheat(
                self.adel_wheat,
                mtg=self.g,
                indice_wheat_instance=self.wheat_index,
                stem_name="StemElement",
                leaf_name="LeafElement1",
            )

        elif self.generation_type == "row":
            scene_wheat = planter.generate_row_wheat(
                self.adel_wheat, self.g, self.wheat_index, stem_name="StemElement", leaf_name="LeafElement1"
            )

        else:
            print("can't recognize positions generation type, choose between default, random and row")
            raise

        stems = extract_stems_from_MTG(self.g, self.global_index)

        return scene_wheat, stems

    def light_results(self, energy, lighting, selective_global_index=None):
        if selective_global_index is not None:
            saved_global_index = self.global_index
            self.global_index = selective_global_index

        results = lighting.results_organs()
        lightmodel = lighting.lightmodel

        # crée un tableau comme dans caribu_facade de fspm-wheat
        dico_par = {}
        para_dic = {}
        erel_dic = {}

        df_outputs_esp = results[results.VegetationType == self.global_index]
        for s in df_outputs_esp["Organ"]:
            d = df_outputs_esp[df_outputs_esp.Organ == s]

            if lightmodel == "caribu":
                para_dic[s] = d["par Eabs"].values[0] * energy
                erel_dic[s] = d["par Eabs"].values[0]  # lighting ran with energy = 1.

            elif lightmodel == "ratp":
                para_dic[s] = d["PARa"].values[0] * energy
                erel_dic[s] = d["Intercepted"].values[0]

        dico_par["PARa"] = para_dic
        dico_par["Erel"] = erel_dic

        for param in dico_par:
            if param not in self.g.properties():
                self.g.add_property(param)
            # update the self.g
            self.g.property(param).update(dico_par[param])

        if selective_global_index is not None:
            self.global_index = saved_global_index

    def soil_inputs(self, soil, planter, lighting):
        # ls_N
        N_content_roots = self.compute_N_content_roots()
        N_content_roots_per_plant = [N_content_roots] * self.nb_plants

        # ls_roots
        roots_length_per_plant_per_soil_layer = self.compute_roots_length(soil, planter)

        # ls_epsi
        organs_results = lighting.results_organs()
        filtered_data = organs_results[organs_results.VegetationType.isin([self.global_index])]
        plant_leaf_area = numpy.sum(filtered_data["Area"].values) / self.nb_plants
        plants_light_interception = self.compute_plants_light_interception(plant_leaf_area, lighting.soil_energy())

        return (
            N_content_roots_per_plant,
            roots_length_per_plant_per_soil_layer,
            self.soil_parameters,
            plants_light_interception,
        )

    def soil_results(self, uptakeN_per_plant, planter=None, selective_global_index=None):
        if selective_global_index is not None:
            saved_global_index = self.global_index
            self.global_index = selective_global_index

        if planter is not None:
            index_in_global_plants = [
                sum(planter.number_of_plants[: self.global_index]),
                sum(planter.number_of_plants[: self.global_index + 1]),
            ]
        else:
            index_in_global_plants = [0, self.nb_plants]

        # considère que tous les blé sont identiques: moyenne
        uptake_nitrate = 0.0
        for i in range(*index_in_global_plants):
            uptake_nitrate += numpy.sum(uptakeN_per_plant[i])
        uptake_nitrate *= 1 / self.nb_plants

        self.uptake_nitrate_hour = self.convert_uptake_nitrate(uptake_nitrate)

        if selective_global_index is not None:
            self.global_index = saved_global_index

    def run(self, t_light):
        if not ((t_light % self.LIGHT_TIMESTEP == 0) and (self.PARi_next_hours(t_light) > 0)):
            Erel = self.g.property("Erel")
            PARa_output = {k: v * self.energy(t_light) for k, v in Erel.items()}
            outputs = {}
            outputs.update({"PARa": PARa_output})
            for param in outputs.keys():
                if param not in self.g.properties():
                    self.g.add_property(param)
                # update the MTG
                self.g.property(param).update(outputs[param])

        # suite de la simu
        for t_senescwheat in range(t_light, t_light + self.SENESCWHEAT_TIMESTEP, self.SENESCWHEAT_TIMESTEP):
            # run SenescWheat
            self.senescwheat_facade_.run()

            # Test for dead plant # TODO: adapt in case of multiple plants
            if (
                not self.shared_elements_inputs_outputs_df.empty
                and numpy.nansum(
                    self.shared_elements_inputs_outputs_df.loc[
                        self.shared_elements_inputs_outputs_df["element"].isin(["StemElement", "LeafElement1"]),
                        "green_area",
                    ]
                )
                == 0
            ):
                # append the inputs and outputs at current step to global lists
                self.all_simulation_steps.append(t_senescwheat)
                self.axes_all_data_list.append(self.shared_axes_inputs_outputs_df.copy())
                self.organs_all_data_list.append(self.shared_organs_inputs_outputs_df.copy())
                self.hiddenzones_all_data_list.append(self.shared_hiddenzones_inputs_outputs_df.copy())
                self.elements_all_data_list.append(self.shared_elements_inputs_outputs_df.copy())
                self.soils_all_data_list.append(self.shared_soils_inputs_outputs_df.copy())
                break

            # Run the rest of the model if the plant is alive
            for t_farquharwheat in range(
                t_senescwheat, t_senescwheat + self.SENESCWHEAT_TIMESTEP, self.FARQUHARWHEAT_TIMESTEP
            ):
                # get the meteo of the current step
                Ta, ambient_CO2, RH, Ur = self.meteo.loc[
                    t_farquharwheat, ["air_temperature", "ambient_CO2", "humidity", "Wind"]
                ]

                # run FarquharWheat
                self.farquharwheat_facade_.run(Ta, ambient_CO2, RH, Ur)

                for t_elongwheat in range(
                    t_farquharwheat, t_farquharwheat + self.FARQUHARWHEAT_TIMESTEP, self.ELONGWHEAT_TIMESTEP
                ):
                    # run ElongWheat
                    Tair, Tsoil = self.meteo.loc[t_elongwheat, ["air_temperature", "soil_temperature"]]
                    self.elongwheat_facade_.run(Tair, Tsoil, option_static=self.option_static)

                    # Update geometry
                    self.adel_wheat.update_geometry(self.g)

                    for t_growthwheat in range(
                        t_elongwheat, t_elongwheat + self.ELONGWHEAT_TIMESTEP, self.GROWTHWHEAT_TIMESTEP
                    ):
                        # run GrowthWheat
                        self.growthwheat_facade_.run()

                        for t_cnwheat in range(
                            t_growthwheat, t_growthwheat + self.GROWTHWHEAT_TIMESTEP, self.CNWHEAT_TIMESTEP
                        ):
                            print("t cnwheat is {}".format(t_cnwheat))

                            # N fertilization if any
                            if (
                                not self.external_soil_model
                                and self.N_fertilizations is not None
                                and len(self.N_fertilizations) > 0
                            ):
                                if t_cnwheat in self.N_fertilizations.keys():
                                    self.cnwheat_facade_.soils[(1, "MS")].nitrates += self.N_fertilizations[t_cnwheat]

                            if t_cnwheat > 0:
                                if self.external_soil_model:
                                    if self.nitrates_uptake_forced:
                                        self.force_nitrates_uptake(t_cnwheat)
                                    else:
                                        self.update_Nitrates_cnwheat_mtg()

                                # run CNWheat
                                Tair = self.meteo.loc[t_elongwheat, "air_temperature"]
                                Tsoil = self.meteo.loc[t_elongwheat, "soil_temperature"]
                                self.cnwheat_facade_.run(Tair, Tsoil, self.tillers_replications)

                            # append outputs at current step to global lists
                            if (self.stored_times == "all") or (t_cnwheat in self.stored_times):
                                (
                                    axes_outputs,
                                    elements_outputs,
                                    hiddenzones_outputs,
                                    organs_outputs,
                                    soils_outputs,
                                ) = self.fspmwheat_facade_.build_outputs_df_from_MTG()

                                self.all_simulation_steps.append(t_cnwheat)
                                self.axes_all_data_list.append(axes_outputs)
                                self.organs_all_data_list.append(organs_outputs)
                                self.hiddenzones_all_data_list.append(hiddenzones_outputs)
                                self.elements_all_data_list.append(elements_outputs)
                                self.soils_all_data_list.append(soils_outputs)

    def end(self, run_postprocessing=False):
        PRECISION = 4

        outputs_brut = os.path.join(self.out_folder, "brut")

        # save des données brutes
        outputs_df_dict = {}
        for outputs_df_list, outputs_filename, index_columns in (
            (self.axes_all_data_list, self.AXES_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.AXES_T_INDEXES),
            (self.organs_all_data_list, self.ORGANS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.ORGANS_T_INDEXES),
            (
                self.hiddenzones_all_data_list,
                self.HIDDENZONES_OUTPUTS_FILENAME,
                cnwheat_simulation.Simulation.HIDDENZONE_T_INDEXES,
            ),
            (
                self.elements_all_data_list,
                self.ELEMENTS_OUTPUTS_FILENAME,
                cnwheat_simulation.Simulation.ELEMENTS_T_INDEXES,
            ),
            (self.soils_all_data_list, self.SOILS_OUTPUTS_FILENAME, cnwheat_simulation.Simulation.SOILS_T_INDEXES),
        ):
            data_filepath = os.path.join(outputs_brut, outputs_filename)
            outputs_df = pandas.concat(outputs_df_list, keys=self.all_simulation_steps, sort=False)
            outputs_df.reset_index(0, inplace=True)
            outputs_df.rename({"level_0": "t"}, axis=1, inplace=True)
            outputs_df = outputs_df.reindex(
                index_columns + outputs_df.columns.difference(index_columns).tolist(), axis=1, copy=False
            )
            outputs_df.fillna(value=numpy.nan, inplace=True)  # Convert back None to NaN
            save_df_to_csv(outputs_df, data_filepath, PRECISION)
            outputs_file_basename = outputs_filename.split(".")[0]
            outputs_df_dict[outputs_file_basename] = outputs_df.reset_index()

        # construit un premier postprocessing
        if run_postprocessing:
            POSTPROCESSING_DIRPATH = os.path.join(self.out_folder, "postprocessing")

            time_grid = list(outputs_df_dict.values())[0].t
            delta_t = (time_grid.loc[1] - time_grid.loc[0]) * self.HOUR_TO_SECOND_CONVERSION_FACTOR

            axes_postprocessing_file_basename = self.AXES_POSTPROCESSING_FILENAME.split(".")[0]
            hiddenzones_postprocessing_file_basename = self.HIDDENZONES_POSTPROCESSING_FILENAME.split(".")[0]
            organs_postprocessing_file_basename = self.ORGANS_POSTPROCESSING_FILENAME.split(".")[0]
            elements_postprocessing_file_basename = self.ELEMENTS_POSTPROCESSING_FILENAME.split(".")[0]
            soils_postprocessing_file_basename = self.SOILS_POSTPROCESSING_FILENAME.split(".")[0]

            postprocessing_df_dict = {}
            (
                postprocessing_df_dict[axes_postprocessing_file_basename],
                postprocessing_df_dict[hiddenzones_postprocessing_file_basename],
                postprocessing_df_dict[organs_postprocessing_file_basename],
                postprocessing_df_dict[elements_postprocessing_file_basename],
                postprocessing_df_dict[soils_postprocessing_file_basename],
            ) = cnwheat_facade.CNWheatfacade.postprocessing(
                axes_outputs_df=outputs_df_dict[self.AXES_OUTPUTS_FILENAME.split(".")[0]],
                hiddenzone_outputs_df=outputs_df_dict[self.HIDDENZONES_OUTPUTS_FILENAME.split(".")[0]],
                organs_outputs_df=outputs_df_dict[self.ORGANS_OUTPUTS_FILENAME.split(".")[0]],
                elements_outputs_df=outputs_df_dict[self.ELEMENTS_OUTPUTS_FILENAME.split(".")[0]],
                soils_outputs_df=outputs_df_dict[self.SOILS_OUTPUTS_FILENAME.split(".")[0]],
                delta_t=delta_t,
            )

            for postprocessing_file_basename, postprocessing_filename in (
                (axes_postprocessing_file_basename, self.AXES_POSTPROCESSING_FILENAME),
                (hiddenzones_postprocessing_file_basename, self.HIDDENZONES_POSTPROCESSING_FILENAME),
                (organs_postprocessing_file_basename, self.ORGANS_POSTPROCESSING_FILENAME),
                (elements_postprocessing_file_basename, self.ELEMENTS_POSTPROCESSING_FILENAME),
                (soils_postprocessing_file_basename, self.SOILS_POSTPROCESSING_FILENAME),
            ):
                postprocessing_filepath = os.path.join(POSTPROCESSING_DIRPATH, postprocessing_filename)
                postprocessing_df_dict[postprocessing_file_basename].to_csv(
                    postprocessing_filepath, na_rep="NA", index=False, float_format="%.{}f".format(PRECISION)
                )

        self.outputs_df_dict = outputs_df_dict

    def update_Nitrates_cnwheat_mtg(self):
        mtg_plants_iterator = self.g.components_iter(self.g.root)
        for plant in self.cnwheat_facade_.population.plants:
            cnwheat_plant_index = plant.index
            while True:
                mtg_plant_vid = next(mtg_plants_iterator)
                if int(self.g.index(mtg_plant_vid)) == cnwheat_plant_index:
                    break
            mtg_axes_iterator = self.g.components_iter(mtg_plant_vid)
            for axis in plant.axes:
                # update in cnwheat_facade_
                axis.roots.__dict__["Uptake_Nitrates"] = self.uptake_nitrate_hour

                # Update Nitrates uptake in MTG
                cnwheat_axis_label = axis.label
                while True:
                    mtg_axis_vid = next(mtg_axes_iterator)
                    if self.g.label(mtg_axis_vid) == cnwheat_axis_label:
                        break
                mtg_roots_properties = self.g.get_vertex_property(mtg_axis_vid)["roots"]
                mtg_roots_properties["Uptake_Nitrates"] = self.uptake_nitrate_hour

    def compute_N_content_roots(self):
        organs_df = self.organs_all_data_list[-1]
        hiddenzones_df = self.hiddenzones_all_data_list[-1]
        elements_df = self.elements_all_data_list[-1]

        # constantes
        C_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.C_MOLAR_MASS
        N_MOLAR_MASS = cnwheat_model.EcophysiologicalConstants.N_MOLAR_MASS
        HEXOSE_MOLAR_MASS_C_RATIO = cnwheat_model.EcophysiologicalConstants.HEXOSE_MOLAR_MASS_C_RATIO
        NITRATES_MOLAR_MASS_N_RATIO = cnwheat_model.EcophysiologicalConstants.NITRATES_MOLAR_MASS_N_RATIO
        AMINO_ACIDS_MOLAR_MASS_N_RATIO = cnwheat_model.EcophysiologicalConstants.AMINO_ACIDS_MOLAR_MASS_N_RATIO

        # récupère grandeurs
        structure = organs_df.fillna(0)["structure"]
        sucrose = organs_df.fillna(0)["sucrose"]
        starch = organs_df.fillna(0)["starch"]
        nitrates = organs_df.fillna(0)["nitrates"]
        amino_acids = organs_df.fillna(0)["amino_acids"]
        proteins = organs_df.fillna(0)["proteins"]
        mstruct = organs_df.fillna(0)["mstruct"]
        Nstruct = organs_df.fillna(0)["Nstruct"]

        ## DRY MASS
        organs_df["sum_dry_mass"] = (
            ((structure + starch) * 1e-6 * C_MOLAR_MASS / HEXOSE_MOLAR_MASS_C_RATIO)
            + (sucrose * 1e-6 * C_MOLAR_MASS) / HEXOSE_MOLAR_MASS_C_RATIO
            + (starch * 1e-6 * C_MOLAR_MASS) / HEXOSE_MOLAR_MASS_C_RATIO
            + (nitrates * 1e-6 * N_MOLAR_MASS) / NITRATES_MOLAR_MASS_N_RATIO
            + (amino_acids * 1e-6 * N_MOLAR_MASS) / AMINO_ACIDS_MOLAR_MASS_N_RATIO
            + (proteins * 1e-6 * N_MOLAR_MASS) / AMINO_ACIDS_MOLAR_MASS_N_RATIO
            + mstruct
        )
        dry_mass_roots = (
            organs_df[(organs_df["organ"] == "roots")].groupby(["plant", "axis"])["sum_dry_mass"].agg("sum")
        )

        # calcule de la proportion partie roots du phloem
        hz_df_MS = hiddenzones_df[hiddenzones_df["axis"] == "MS"].copy()
        elt_df_MS = elements_df[elements_df["axis"] == "MS"].copy()
        hz_df_MS["mstruct_tillers"] = hz_df_MS["mstruct"] * hz_df_MS["nb_replications"]
        elt_df_MS["mstruct_tillers"] = elt_df_MS["mstruct"] * elt_df_MS["nb_replications"]
        sum_mstruct_shoot = hz_df_MS.groupby(["plant", "axis"])["mstruct_tillers"].agg("sum") + elt_df_MS.groupby(
            ["plant", "axis"]
        )["mstruct_tillers"].agg("sum")
        sum_mstruct_shoot.fillna(0, inplace=True)
        sum_mstruct_roots = organs_df[(organs_df["organ"] == "roots")].groupby(["plant", "axis"])["mstruct"].agg("sum")
        shoot_roots_mstruct_ratio = sum_mstruct_shoot / sum_mstruct_roots
        phloem_shoot_root = 1 / (1 + 1 / shoot_roots_mstruct_ratio)

        # dry mass de la partie roots du phloem
        sum_dry_mass_phloem = (
            organs_df[(organs_df["organ"] == "phloem")].groupby(["plant", "axis"])["sum_dry_mass"].agg("sum")
        )
        sum_dry_mass_phloem_roots = sum_dry_mass_phloem * (1 - phloem_shoot_root)

        # dry total des roots
        sum_dry_mass_roots = sum_dry_mass_phloem_roots + dry_mass_roots

        ## N MASS
        organs_df["N_g"] = (
            (nitrates * 1e-6 * N_MOLAR_MASS)
            + (amino_acids * 1e-6 * N_MOLAR_MASS)
            + (proteins * 1e-6 * N_MOLAR_MASS)
            + Nstruct
        )

        # N mass du phloem roots
        sum_N_g_phloem = organs_df[(organs_df["organ"] == "phloem")].groupby(["plant", "axis"])["N_g"].agg("sum")
        sum_N_g_phloem_shoot = sum_N_g_phloem * (1 - phloem_shoot_root)

        # N mass du reste des roots
        sum_N_g_roots = organs_df[(organs_df["organ"] == "roots")].groupby(["plant", "axis"])["N_g"].agg("sum")

        # N mass totale des roots
        sum_N_g_roots = sum_N_g_roots + sum_N_g_phloem_shoot

        # finalement pourcentage de N dans les racines
        N_content_roots = sum_N_g_roots / sum_dry_mass_roots * 100

        return N_content_roots.values[0]

    def compute_roots_length(self, soil_wrapper, planter):
        """
        soil_dimensions : [z, x, y]
        """
        # une seule plante dans cnwheat
        roots_mass = [0]
        for i, plant in enumerate(self.cnwheat_facade_.population.plants):
            for axis in plant.axes:
                roots_mass[i] += axis.roots.mstruct  # masse en g

        positions = planter.wheat_positions
        self.compute_SRL_wheat(roots_mass[0])

        # longueur spécifique x masse en gramme/nbplantes
        ls_roots = []
        for p in positions:
            # on répartit de manière homogène les racines à travers les couches du sol
            # convertit m en cm # --> peut etre en metre finalement
            ix, iy = soil_wrapper.whichvoxel_xy(p)
            roots_length_per_voxel = self.rootsdistribution(roots_mass[0], ix, iy, soil_wrapper)
            ls_roots.append(roots_length_per_voxel)

        return ls_roots

    def rootsdistribution(self, roots_mass, ix, iy, soil_wrapper, distribtype="homogeneous"):
        roots_length_per_voxel = numpy.zeros(soil_wrapper.soil_dimensions)
        if distribtype == "homogeneous":
            roots_length_per_voxel[:, ix, iy] = (roots_mass * self.SRL) / soil_wrapper.soil_dimensions[0]
        
        return roots_length_per_voxel

    def compute_SRL_wheat(self, mass_roots):
        if mass_roots < 0.6:
            a = 334
            self.SRL =  a * mass_roots + 60
            # a = 8.247933150630281 # math.log(141)/0.6
            # return np.exp(mass_roots * a) + 59
        else:
            self.SRL = 200

    def compute_plants_light_interception(self, plant_leaf_area, soil_energy):
        # conversion
        c = (3600 * 24) / 1000000
        # portion au sol pour chaque plante
        return [(1 - soil_energy) * ((plant_leaf_area * c) / self.nb_plants)] * self.nb_plants

    @staticmethod
    def convert_uptake_nitrate(uptake_nitrate):
        # masse atomic de l'azote 14.0067 g.mol-1
        atomic_mass_N = 14.0067

        # conversion kg en µmol d'azote
        uptake_nitrate_g = uptake_nitrate * 10**3
        uptake_mol = uptake_nitrate_g / atomic_mass_N
        uptake_micromol = uptake_mol * 10**6

        # daily_uptake_nitrate = 10**9 * uptake_nitrate / atomic_mass_N

        # on répartie la grandeur sur chaque heure
        daily_uptake_nitrate = uptake_micromol / 24

        return daily_uptake_nitrate

    def force_nitrates_uptake(self, t):
        mtg_plants_iterator = self.g.components_iter(self.g.root)
        for plant in self.cnwheat_facade_.population.plants:
            cnwheat_plant_index = plant.index
            while True:
                mtg_plant_vid = next(mtg_plants_iterator)
                if int(self.g.index(mtg_plant_vid)) == cnwheat_plant_index:
                    break
            mtg_axes_iterator = self.g.components_iter(mtg_plant_vid)
            for axis in plant.axes:
                # Update Nitrates uptake in dataframe
                group = self.nitrates_uptake_data_grouped.get_group((t, plant.index, axis.label, "roots"))
                nitrates_uptake_data_to_use = (
                    group.loc[group.first_valid_index(), group.columns.intersection(["Uptake_Nitrates"])]
                    .dropna()
                    .to_dict()
                )
                axis.roots.__dict__.update(nitrates_uptake_data_to_use)

                # Update Nitrates uptake in MTG
                cnwheat_axis_label = axis.label
                while True:
                    mtg_axis_vid = next(mtg_axes_iterator)
                    if self.g.label(mtg_axis_vid) == cnwheat_axis_label:
                        break
                mtg_roots_properties = self.g.get_vertex_property(mtg_axis_vid)["roots"]
                mtg_roots_properties.update(nitrates_uptake_data_to_use)

    def energy(self, t):
        return self.meteo.loc[t, ["PARi"]].iloc[0]

    def doy(self, t, soil3ds=False):
        if soil3ds:
            if t > 0:
                if self.meteo.loc[t - 1, ["DOY"]].iloc[0] > self.meteo.loc[t, ["DOY"]].iloc[0]:
                    self.last_year_doy += self.meteo.loc[t - 1, ["DOY"]].iloc[0]
            return self.meteo.loc[t, ["DOY"]].iloc[0] + self.last_year_doy

        else:
            return self.meteo.loc[t, ["DOY"]].iloc[0]

    def hour(self, t):
        return self.meteo.loc[t, ["hour"]].iloc[0]

    def PARi_next_hours(self, t):
        return self.meteo.loc[range(t, t + self.LIGHT_TIMESTEP), ["PARi"]].sum().values[0]

    def next_day_next_hour(self, t):
        return self.meteo.loc[t + self.SENESCWHEAT_TIMESTEP, ["DOY"]].iloc[0]

    @staticmethod
    def fake_scene():
        epsilon = 1e-14
        return {
            19: [[(0.0, 0.0, 0.0), (0.0, epsilon, 0.0), (0.0, epsilon, epsilon)]],
            34: [[(0.0, 0.0, 0.0), (epsilon, 0.0, 0.0), (epsilon, 0.0, epsilon)]],
        }


def passive_lighting(data, t, DOY, scene, lighting_wrapper, stems=None):
    lighting_wrapper.run(scenes=[scene], day=DOY, parunit="micromol.m-2.s-1", stems=stems)

    results = lighting_wrapper.results_organs()

    para = results["Organ"] * results["Area"]
    para *= 1 / results["Area"].sum()

    data["PARa"].append(para)
    data["t"].append(t)
