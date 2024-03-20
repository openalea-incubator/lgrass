"""
    LightVegeManager
    ****************

    Main class of the tool. Calls all the other modules in ``src``.

    3 inputs dict for setting all parameters:

    .. code-block:: python
    
        geometry = {
                    "scenes" : [scene0, scene1, scene2, ...] ,
                    "domain" : ((xmin, ymin), (xmax, ymax)),
                    "stems id" : [(id_element, id_scene), ...],
                    "transformations" : {
                                            "scenes unit" : kwarg ,
                                            "rescale" : kwarg ,
                                            "translate" : kwarg ,
                                            "xyz orientation" : kwarg
                                            }
                        }

    .. code-block:: python

        environment = {
                        "coordinates" : [latitude, longitude, timezone] ,
                        
                        "sky" : "turtle46" ,
                        "sky" : ["file", filepath] ,
                        "sky" : [nb_azimut, nb_zenith, "soc" or "uoc"] ,

                        "direct" : bool, # sun radiations
                        "diffus" : bool, # sky radiations
                        "reflected" : bool, # reflected radiation in the canopy
                        "infinite" : bool, # infinitisation of the scene
                        }

    Currently LightVegeManager handles the light models RATP and CARIBU:

    .. code-block:: python

        caribu_args = {
                        "sun algo" : "ratp",
                        "sun algo" : "caribu",

                            "caribu opt" : {
                                            band0 = (reflectance, transmittance),
                                            band1 = (reflectance, transmittance),
                                            ...
                                            },
                            "debug" : bool,
                            "soil mesh" : bool,
                            "sensors" : ["grid", dxyz, nxyz, orig]
                        }

    .. code-block:: python

        ratp_args = {
                        # Grid specifications
                        "voxel size" : [dx, dy, dz],
                        "voxel size" : "dynamic",
                        
                        "origin" : [xorigin, yorigin, zorigin],
                        "origin" : [xorigin, yorigin],

                        "number voxels" : [nx, ny, nz],
                        "grid slicing" : "ground = 0."
                        "tesselation level" : int

                        # Leaf angle distribution
                        "angle distrib algo" : "compute global",
                        "angle distrib algo" : "compute voxel",
                        "angle distrib algo" : "file",

                        "nb angle classes" : int,
                        "angle distrib file" : filepath,

                        # Vegetation type
                        "soil reflectance" : [reflectance_band0, reflectance_band1, ...],
                        "reflectance coefficients" : [reflectance_band0, reflectance_band1, ...],
                        "mu" : [mu_scene0, mu_scene1, ...]
                    }

    .. seealso:: For more details :ref:`Inputs description <inputs>`
"""
import time
import os
import subprocess
import numpy
import math

from lightvegemanager.trianglesmesh import (
    isatriangle,
    chain_triangulations,
    apply_transformations,
    compute_area_max,
    compute_minmax_coord,
    compute_trilenght_max,
)
from lightvegemanager.VTK import VTKtriangles
from lightvegemanager.defaultvalues import default_LightVegeManager_inputs


class LightVegeManager(object):
    """Main class for the tool LightVegeManager

    Common simulation order:

    ``input geometries -> build and prepare data -> call light model -> transfer results to plant models``

    It includes:

        Main methods:

        * ``__init__``: initializes and builds static object for the rest of simulation
        * :meth:`build`: builds and prepare all geometric meshes
        * :meth:`run`: calls a light model and manages its inputs and outputs

        Transfer methods:

        * :meth:`to_MTG`: transfers ligthing results to a MTG table
        * :meth:`to_l_egume`: transfers ligthing results to l-egume by creating two arrays as inputs for the plant model

        Analysis tools: analyses a set of triangles and formats them as a turbid medium inputs

        * :meth:`s5`: fortran tool
        * :meth:`s2v`: c++ tool

        Visualisation tools:

        * :meth:`plantGL_nolight`: return a plantGL scene from the geometry
        * :meth:`plantGL_light`: return a plantGL scene from the geometry with lighting results
        * :meth:`plantGL_sensors`: return a plantGL scene from virtual sensors if created

        * :meth:`VTK_nolight`: write VTK file with only geometric informations
        * :meth:`VTK_light`: write VTK file with geometric informations and associated light results
        * :meth:`VTK_sun`: write VTK file representing the sun as a line

        Getters: to use light results with external routines

        * :meth:`riri5_transmitted_light`
        * :meth:`riri5_intercepted_light`
        * :meth:`elements_outputs`
        * :meth:`triangles_outputs`
        * :meth:`voxels_outputs`
        * :meth:`sensors_outputs`
        * :meth:`sun`
        * :meth:`soilenergy`
        * :meth:`maxtrianglearea`
        * :meth:`legume_empty_layers`
        * :meth:`tesselationtime`
        * :meth:`modelruntime`
        * :meth:`leafangledistribution`

    :param environment: Environment parameters, defaults to {}
    :type environment: dict, optional
    :param lightmodel: either ``"ratp"`` or ``"caribu"``, defaults to ""
    :type lightmodel: str, optional
    :param lightmodel_parameters: light model parameters, defaults to {}
    :type lightmodel_parameters: dict, optional
    :param main_unit: measure unit for the global scene where the light will be computed, defaults to "m"
    :type main_unit: str, optional
    :raises ValueError: lightmodel entry not valid, either ``'ratp'`` or ``'caribu'``
    """

    ## MAIN METHODS ##
    def __init__(self, environment={}, lightmodel="", lightmodel_parameters={}, main_unit="m"):
        # gets default variables from LightVegeManager_defaultvalues.py
        (
            default_environnement,
            default_ratp_parameters,
            default_caribu_parameters,
            default_riri5_parameters,
        ) = default_LightVegeManager_inputs()

        # check if choosen light model is known by the tool
        if lightmodel != "ratp" and lightmodel != "caribu" and lightmodel != "riri5":
            raise ValueError("Unknown lightmodel: can be either 'ratp' or 'caribu' or 'riri5' ")

        # save inputs in the instance
        self.__environment = environment
        self.__lightmodel = lightmodel
        self.__lightmodel_parameters = lightmodel_parameters
        self.__main_unit = main_unit

        # copy of default values in input parameters on not initialized keys
        for key, value in default_environnement.items():
            if key not in self.__environment:
                self.__environment[key] = value

        if lightmodel == "caribu":
            lghtdict = default_caribu_parameters
        elif lightmodel == "ratp" or lightmodel == "riri5":
            lghtdict = default_ratp_parameters
            if lightmodel == "riri5":
                lghtdict.update(default_riri5_parameters)

        for key, value in lghtdict.items():
            if key not in self.__lightmodel_parameters:
                self.__lightmodel_parameters[key] = value

        # sky building
        skytype = self.__environment["sky"]
        if lightmodel == "caribu":
            if self.__environment["diffus"]:
                from lightvegemanager.sky import CARIBUsky

                self.__sky = CARIBUsky(skytype)

        elif lightmodel == "ratp":
            from lightvegemanager.sky import RATPsky

            self.__sky = RATPsky(skytype)

        elif lightmodel == "riri5":
            if skytype == "turtle46":
                self.__sky = ["soc", "VXpXmYpYm"]
            else:
                self.__sky = skytype

    def build(self, geometry={}, global_scene_tesselate_level=0):
        """Builds a mesh of the simulation scene in the right light model format

        :param geometry: geometric parameters, contains geometric scenes, defaults to {}
        :type geometry: dict, optional
        :param global_scene_tesselate_level: option to subdivide all triangles of the mesh a certain number of times (to fine tuning the mesh), defaults to 0
        :type global_scene_tesselate_level: int, optional
        :raises ValueError: Currently, converting voxels mesh to triangles mesh is not possible
        """
        # pre-check of scenes input, if it has only one triangle or one list of triangles
        if isatriangle(geometry) or all(isatriangle(s) for s in geometry):
            geometry = {"scenes": geometry}

        self.__geometry = geometry

        # First process of the scenes list, it gathers all triangulations
        self.__complete_trimesh, self.__matching_ids, legume_grid, id_legume_scene = chain_triangulations(
            self.__geometry["scenes"]
        )

        # applies geometric transformations on inputs scenes if precised
        if "transformations" in self.__geometry:
            apply_transformations(
                self.__complete_trimesh, self.__matching_ids, self.__geometry["transformations"], self.__main_unit
            )

        self.__areamax = compute_area_max(self.__complete_trimesh)

        # global tesselation of triangulation
        if self.__matching_ids and global_scene_tesselate_level > 0:
            from lightvegemanager.tesselator import iterate_triangles

            new_trimesh = {}
            for id_ele, triangles in self.__complete_trimesh.items():
                new_tr_scene = []
                for t in triangles:
                    level = 0
                    iterate_triangles(t, level, global_scene_tesselate_level, new_tr_scene)
                new_trimesh[id_ele] = new_tr_scene
            self.__complete_trimesh = new_trimesh

        self.__pmin, self.__pmax = compute_minmax_coord(self.__complete_trimesh)

        self.__triangleLmax = compute_trilenght_max(self.__complete_trimesh)

        if self.__lightmodel == "caribu":
            if legume_grid:
                raise ValueError(
                    "Conversion from voxels grid to triangles \
                                is not possible yet"
                )
            if "stems id" not in self.__geometry:
                self.__geometry["stems id"] = None

            # build sensors
            self.__issensors = "sensors" in self.__lightmodel_parameters

            if self.__issensors and (
                self.__lightmodel_parameters["sensors"][0] == "grid"
                or self.__lightmodel_parameters["sensors"][0][0] == "grid"
            ):
                from lightvegemanager.CARIBUinputs import create_caribu_legume_sensors
                from lightvegemanager.voxelsmesh import reduce_layers_from_trimesh

                try:
                    import openalea.plantgl.all as pgl
                except ImportError:
                    pass

                if isinstance(self.__lightmodel_parameters["sensors"], list):
                    dxyz = self.__lightmodel_parameters["sensors"][1]
                    nxyz = self.__lightmodel_parameters["sensors"][2]
                    orig = self.__lightmodel_parameters["sensors"][3]
                    arg = (dxyz, nxyz, orig, self.__pmax, self.__complete_trimesh, self.__matching_ids, None, True)
                    sensors_caribu, sensors_plantgl, Pmax_capt = create_caribu_legume_sensors(*arg)

                    id = [-1]
                    self.__sensors_plantgl = sensors_plantgl
                    self.__sensors_caribu = sensors_caribu
                    self.__nb0 = reduce_layers_from_trimesh(
                        self.__complete_trimesh,
                        self.__pmax,
                        self.__lightmodel_parameters["sensors"][1],
                        self.__lightmodel_parameters["sensors"][2],
                        self.__matching_ids,
                        id,
                    )
                elif isinstance(self.__lightmodel_parameters["sensors"], dict):
                    start_id = 0
                    sensors_plantgl = pgl.Scene()
                    sensors_caribu = {}
                    nb0 = 9999999
                    id = [-1]
                    for specy_indice, sensors_parameters in self.__lightmodel_parameters["sensors"].items():
                        dxyz = sensors_parameters[1]
                        nxyz = sensors_parameters[2]
                        orig = sensors_parameters[3]
                        arg = (
                            dxyz,
                            nxyz,
                            orig,
                            self.__pmax,
                            self.__complete_trimesh,
                            self.__matching_ids,
                            None,
                            True,
                            start_id,
                        )
                        _caribu, _plantgl, Pmax_capt = create_caribu_legume_sensors(*arg)
                        sensors_caribu.update(_caribu)
                        sensors_plantgl += _plantgl
                        start_id = max(sensors_caribu.keys())

                        _nb0 = reduce_layers_from_trimesh(
                            self.__complete_trimesh,
                            self.__pmax,
                            sensors_parameters[1],
                            sensors_parameters[2],
                            self.__matching_ids,
                            id,
                        )
                        nb0 = min(nb0, _nb0)
                    self.__nb0 = nb0
                    self.__sensors_plantgl = sensors_plantgl
                    self.__sensors_caribu = sensors_caribu

        # Builds voxels grid from input geometry
        elif self.__lightmodel == "ratp" or self.__lightmodel == "riri5":
            # number of input species
            numberofentities = 0

            # initialize number of empty layers
            self.__nb0 = 0

            if legume_grid:
                numberofentities = sum([self.__geometry["scenes"][i]["LA"].shape[0] for i in id_legume_scene])

            # triangles in the inputs
            if self.__matching_ids:
                from lightvegemanager.buildRATPscene import build_RATPscene_from_trimesh

                # separates stem elements in a new specy
                if "stems id" in self.__geometry:
                    from lightvegemanager.stems import manage_stems_for_ratp

                    manage_stems_for_ratp(
                        self.__geometry["stems id"], self.__matching_ids, self.__lightmodel_parameters
                    )
                else:
                    self.__geometry["stems id"] = None

                # search for number of species after stems processing (it adds one separate specy)
                numberofentities += max([v[1] for v in self.__matching_ids.values()]) + 1

                # init mu and reflectance  values if not precised
                self.__lightmodel_parameters["mu"] = [1.0 for x in range(numberofentities)]
                self.__lightmodel_parameters["reflectance coefficients"] = [[0.0, 0.0] for x in range(numberofentities)]

                arg = (
                    self.__complete_trimesh,
                    [self.__pmin, self.__pmax],
                    self.__triangleLmax,
                    self.__matching_ids,
                    self.__lightmodel_parameters,
                    self.__environment["coordinates"],
                    self.__environment["reflected"],
                    self.__environment["infinite"],
                    self.__geometry["stems id"],
                    len(self.__geometry["scenes"]),
                    self.__lightmodel_parameters["full grid"],
                )

                (
                    self.__complete_voxmesh,
                    self.__matching_tri_vox,
                    self.__angle_distrib,
                    self.__complete_trimesh,
                ) = build_RATPscene_from_trimesh(*arg)

            # creates an empty RATP grid of voxels if geometric inputs are empty
            else:
                from lightvegemanager.buildRATPscene import build_RATPscene_empty

                arg = (
                    self.__lightmodel_parameters,
                    [self.__pmin, self.__pmax],
                    self.__environment["coordinates"],
                    self.__environment["infinite"],
                )

                self.__complete_voxmesh, self.__angle_distrib = build_RATPscene_empty(*arg)
                self.__matching_tri_vox = {}

            # if there is a grid of voxels in the inputs, we converts it in a RATP grid
            if legume_grid and not self.__matching_ids and self.__lightmodel == "ratp":
                from lightvegemanager.buildRATPscene import legumescene_to_RATPscene, concatene_legumescenes

                if len(id_legume_scene) > 1:
                    concatened_legumescene = concatene_legumescenes(
                        [self.__geometry["scenes"][i] for i in id_legume_scene]
                    )
                else:
                    concatened_legumescene = self.__geometry["scenes"][id_legume_scene[0]]

                arg = (
                    concatened_legumescene,
                    self.__lightmodel_parameters,
                    self.__environment["coordinates"],
                    self.__environment["reflected"],
                    self.__environment["infinite"],
                )

                self.__complete_voxmesh, self.__angle_distrib, self.__nb0 = legumescene_to_RATPscene(*arg)

        if self.__lightmodel == "riri5":
            dS = self.__lightmodel_parameters["voxel size"][0] * self.__lightmodel_parameters["voxel size"][1]
            if legume_grid:
                self.__LA_riri5 = self.__geometry["scenes"][id_legume_scene[0]]["LA"] / dS
                self.__distrib_riri5 = self.__geometry["scenes"][id_legume_scene[0]]["distrib"]
            if not legume_grid:
                from lightvegemanager.RiRi5inputs import ratpgrid_to_riri5

                self.__LA_riri5 = ratpgrid_to_riri5(self.__complete_voxmesh)
                self.__distrib_riri5 = self.__angle_distrib["global"]

    def run(self, energy=0.0, day=0, hour=0, parunit="micromol.m-2.s-1", truesolartime=False, id_sensors=None):
        """Calls the light model and formats lighting results

        :param energy: input radiation energy, defaults to 0
        :type energy: float, optional
        :param day: simulation day, defaults to 0
        :type day: int, optional
        :param hour: simulation hour, defaults to 0
        :type hour: int, optional
        :param parunit: input energy unit, light models manages radiations in different unit, you can precise input unit and LightVegeManager will convert it in the right unit, defaults to "micromol.m-2.s-1"
        :type parunit: str, optional
        :param truesolartime: simulation hour is a true solar time or local time (depending on simulation coordinates), defaults to False
        :type truesolartime: bool, optional
        :param id_sensors: if you use CARIBU with a grid of virtual sensors, you have to precise which input scenes the grid must match, defaults to None
        :type id_sensors: list, optional
        :raises ValueError: with CARIBU you can precise the sun algorithm to calculate sun position, can be either ``"caribu"`` or ``"ratp"``
        :raises ValueError: valid radiations are ``"direct"``, ``"diffuse"``, ``"reflected"``
        """
        self.__energy = energy

        ## RATP ##
        if self.__lightmodel == "ratp":
            from lightvegemanager.RATPinputs import RATP_vegetation, RATP_meteo

            vegetation = RATP_vegetation(
                self.__lightmodel_parameters, self.__angle_distrib, self.__environment["reflected"]
            )
            meteo = RATP_meteo(
                energy,
                day,
                hour,
                self.__environment["coordinates"],
                parunit,
                truesolartime,
                self.__environment["direct"],
                self.__environment["diffus"],
            )

            if self.__complete_voxmesh.nveg > 0:
                from alinea.pyratp.runratp import runRATP
                from lightvegemanager.outputs import out_ratp_voxels

                # Run of RATP
                start = time.time()
                res = runRATP.DoIrradiation(self.__complete_voxmesh, vegetation, self.__sky, meteo)
                self.__time_runmodel = time.time() - start

                # output management
                self.__voxels_outputs = out_ratp_voxels(self.__complete_voxmesh, res, parunit)

                # if there are triangulations in the inputs
                if self.__matching_ids:
                    from lightvegemanager.outputs import out_ratp_triangles, out_ratp_elements

                    arg = (self.__complete_trimesh, self.__matching_ids, self.__matching_tri_vox, self.__voxels_outputs)
                    self.__triangles_outputs = out_ratp_triangles(*arg)

                    arg = (
                        self.__matching_ids,
                        self.__environment["reflected"],
                        self.__lightmodel_parameters["reflectance coefficients"],
                        self.__triangles_outputs,
                    )
                    self.__elements_outputs = out_ratp_elements(*arg)

            else:
                from lightvegemanager.outputs import out_ratp_empty_grid

                print("--- Empty RATP grid")
                self.__voxels_outputs, self.__triangles_outputs, self.__elements_outputs = out_ratp_empty_grid(
                    day, hour
                )

        ## CARIBU ##
        elif self.__lightmodel == "caribu":
            from lightvegemanager.outputs import out_caribu_elements, out_caribu_triangles

            sun_up = False
            if self.__environment["direct"]:
                # computes sun position
                arg = (day, hour, self.__environment["coordinates"], truesolartime)
                if self.__lightmodel_parameters["sun algo"] == "ratp":
                    from lightvegemanager.sun import ratp_sun

                    self.__sun = ratp_sun(*arg)
                elif self.__lightmodel_parameters["sun algo"] == "caribu":
                    from lightvegemanager.sun import caribu_sun

                    self.__sun = caribu_sun(*arg)
                else:
                    raise ValueError("sun algo not recognize")

                # check if sun is up
                # criteria is set to elevation > 2°, from alinea shortwave_balance.F90 -> DirectBeam_Interception, l.68
                sun_up = math.asin(-self.__sun[2]) * (180 / math.pi) > 2.0

            compute = (sun_up and self.__environment["direct"]) or self.__environment["diffus"]
            if compute:
                from lightvegemanager.CARIBUinputs import Prepare_CARIBU, run_caribu
                from alinea.caribu.CaribuScene import CaribuScene
                from alinea.caribu.sky_tools.spitters_horaire import RdRsH
                from lightvegemanager.outputs import out_caribu_mix, out_caribu_nomix, out_caribu_sensors

                # CARIBU preparations
                arg = (
                    self.__complete_trimesh,
                    self.__geometry,
                    self.__matching_ids,
                    (self.__pmin, self.__pmax),
                    self.__lightmodel_parameters,
                    self.__environment["infinite"],
                    id_sensors,
                )
                opt, sensors_caribu, debug, matching_sensors_species = Prepare_CARIBU(*arg)
                issoilmesh = self.__lightmodel_parameters["soil mesh"] != -1
                self.__domain = self.__geometry["domain"]

                if self.__matching_ids:
                    # Initialize a CaribuScene
                    if self.__environment["diffus"] and self.__environment["direct"]:
                        c_scene_sky = CaribuScene(
                            scene=self.__complete_trimesh,
                            light=self.__sky,
                            opt=opt,
                            scene_unit=self.__main_unit,
                            pattern=self.__geometry["domain"],
                            soil_mesh=self.__lightmodel_parameters["soil mesh"],
                            debug=debug,
                        )
                        sun_sky_option = "mix"
                        light = [tuple((1.0, self.__sun))]

                    else:
                        if self.__environment["diffus"]:
                            light = self.__sky
                            sun_sky_option = "sky"

                        elif self.__environment["direct"]:
                            light = [tuple((1.0, self.__sun))]
                            sun_sky_option = "sun"

                        else:
                            raise ValueError("Error with radiative inputs")

                    c_scene = CaribuScene(
                        scene=self.__complete_trimesh,
                        light=light,
                        opt=opt,
                        scene_unit=self.__main_unit,
                        pattern=self.__geometry["domain"],
                        soil_mesh=self.__lightmodel_parameters["soil mesh"],
                        debug=debug,
                    )

                    # Runs CARIBU
                    arg = [
                        c_scene,
                        not self.__environment["reflected"],
                        self.__environment["infinite"],
                        sensors_caribu,
                        self.__energy,
                    ]
                    if sun_sky_option == "mix":
                        start = time.time()
                        raw_sun, aggregated_sun = run_caribu(*arg)
                        arg[0] = c_scene_sky
                        raw_sky, aggregated_sky = run_caribu(*arg)
                        self.__time_runmodel = time.time() - start

                        #: Spitters's model estimating for the diffuse:direct ratio
                        # % de sky dans la valeur d'énergie finale
                        Rg = energy / 2.02  #: Global Radiation (W.m-2)
                        #: Diffuse fraction of the global irradiance
                        rdrs = RdRsH(Rg=Rg, DOY=day, heureTU=hour, latitude=self.__environment["coordinates"][0])

                        raw, aggregated, self.__sensors_outputs, self.__soilenergy = out_caribu_mix(
                            rdrs,
                            c_scene,
                            c_scene_sky,
                            raw_sun,
                            aggregated_sun,
                            raw_sky,
                            aggregated_sky,
                            self.__issensors,
                            issoilmesh,
                        )
                    else:
                        start = time.time()
                        raw, aggregated = run_caribu(*arg)
                        self.__time_runmodel = time.time() - start

                        self.__sensors_outputs, self.__soilenergy = out_caribu_nomix(
                            c_scene, aggregated, self.__issensors, issoilmesh
                        )
                else:
                    aggregated, raw = {}, {}
                    self.__sensors_outputs = {"par": {}}
                    for id, triangles in sensors_caribu.items():
                        self.__sensors_outputs["par"][id] = 1.0

                if self.__issensors:
                    self.__sensors_outputs_df = out_caribu_sensors(
                        day, hour, self.__sensors_outputs, matching_sensors_species
                    )

            # Outputs management
            arg = [day, hour, self.__complete_trimesh, self.__matching_ids, raw, compute]
            self.__triangles_outputs = out_caribu_triangles(*arg)
            arg[4] = aggregated
            arg.append(self.__triangles_outputs)
            self.__elements_outputs = out_caribu_elements(*arg)

            # if "sensors" in self.__lightmodel_parameters and \
            #     ((isinstance(self.__lightmodel_parameters["sensors"], list) and self.__lightmodel_parameters["sensors"][-1] == "vtk") \
            #     or (isinstance(self.__lightmodel_parameters["sensors"], dict) and self.__lightmodel_parameters["sensors"][0][-1] == "vtk")):
            #     # create list with radiative value for each sensor triangle (2 triangles per sensor)
            #     var = []
            #     if isinstance(self.__lightmodel_parameters["sensors"], list) :
            #         sensor_path = os.path.join(self.__lightmodel_parameters["sensors"][-2], sensor_name)
            #     elif isinstance(self.__lightmodel_parameters["sensors"], dict) :
            #         sensor_path = os.path.join(self.__lightmodel_parameters["sensors"][0][-2], sensor_name)

            #     for id, triangles in sensors_caribu.items():
            #         for t in triangles:
            #             var.append(self.__sensors_outputs["par"][id])

            #     sensor_name = ("sensors_h"      +
            #                     str(int(hour))  +
            #                     "_d"            +
            #                     str(int(day))   +
            #                     ".vtk")
            #     VTKtriangles(sensors_caribu, [var], ["intercepted"], sensor_path)

        ## RiRi (l-egume) ##
        elif self.__lightmodel == "riri5":
            from riri5.RIRI5 import calc_extinc_allray_multi_reduced

            dS = self.__lightmodel_parameters["voxel size"][0] * self.__lightmodel_parameters["voxel size"][1]

            tag_light_inputs = [
                self.__LA_riri5,
                None,
                self.__distrib_riri5,
                energy * dS,
            ]

            # mise a jour de res_trans, res_abs_i, res_rfr, ls_epsi
            start = time.time()
            self.__riri5_transmitted_light, self.__riri5_intercepted_light = calc_extinc_allray_multi_reduced(
                *tag_light_inputs, optsky=self.__sky[0], opt=self.__sky[1]
            )
            self.__time_runmodel = time.time() - start

    ## TRANSFER METHODS ##
    def to_MTG(self, energy=1.0, mtg=None, id=None):
        """Transfers lighting results to a MTG table.

        .. warning:: The :meth:`run` must have been called before to have results dataframes.

        Results are ``pandas.Dataframe`` stored in the ``self``

        :param energy: input energy, defaults to 1
        :type energy: float, optional
        :param mtg: MTG table with ``"PARa"`` and ``"Erel"`` entries in its properties, defaults to None
        :type mtg: MTG, optional
        :param id: you can precise to which input scenes the MTG table corresponds, defaults to None
        :type id: list or tuple, optional
        :raises AttributeError: you need to call :func:run first
        """
        if not hasattr(self, "_LightVegeManager__elements_outputs"):
            raise AttributeError("No results yet, run a light modeling first")

        # crée un tableau comme dans caribu_facade de fspm-wheat
        dico_par = {}
        para_dic = {}
        erel_dic = {}
        if id is None:
            for s in self.__elements_outputs["Organ"]:
                d = self.__elements_outputs[self.__elements_outputs.Organ == s]

                if self.__lightmodel == "caribu":
                    para_dic[s] = d["par Eabs"].values[0] * energy
                    erel_dic[s] = d["par Eabs"].values[0] / self.__energy

                elif self.__lightmodel == "ratp":
                    para_dic[s] = d["PARa"].values[0]
                    erel_dic[s] = d["Intercepted"].values[0]

        elif type(id) == list or type(id) == tuple:
            for esp in id:
                df_outputs_esp = self.__elements_outputs[self.__elements_outputs.VegetationType == esp]
                for s in df_outputs_esp["Organ"]:
                    d = df_outputs_esp[df_outputs_esp.Organ == s]

                    if self.__lightmodel == "caribu":
                        para_dic[s] = d["par Eabs"].values[0] * energy
                        erel_dic[s] = d["par Eabs"].values[0] / self.__energy

                    elif self.__lightmodel == "ratp":
                        para_dic[s] = d["PARa"].values[0]
                        erel_dic[s] = d["Intercepted"].values[0]

        dico_par["PARa"] = para_dic
        dico_par["Erel"] = erel_dic

        for param in dico_par:
            if param not in mtg.properties():
                mtg.add_property(param)
            # update the MTG
            mtg.property(param).update(dico_par[param])

    def to_l_egume(self, energy=1.0, m_lais=[], list_lstring=[], list_dicFeuilBilanR=[], list_invar=[], id=None):
        """Transfers lighting results to l-egume

        .. warning:: The :meth:`run` must have been called before to have results dataframes.

        .. warning:: l-egume needs transmitted energy informations located in a grid of voxels. You need to have the same dimensions in the lighting results.

            * With RATP, RATP grid must have the same dimensions as l-egume intern grid.

            * With CARIBU, you must create a grid of virtual sensors in the same dimensions as l-egume intern grid.

        Results are ``pandas.Dataframe`` stored in the ``self``

        :param energy: input energy, defaults to 1
        :type energy: float, optional
        :param m_lais: leaf area represented in a numpy.array of dimension [number of species, number of z layers, number of y layers, number of x layers], defaults to []
        :type m_lais: numpy.array, optional
        :param list_lstring: from l-egume, each element corresponds to an input specy of l-egume. Each element is a dict lstring stores the l-system of each plant, defaults to []
        :type list_lstring: list of dict, optional
        :param list_dicFeuilBilanR: from l-egume, each element corresponds to an input specy of l-egume. Each element is a dict dicFeuiBilanR stores correspondances between voxels grid and each plant, defaults to []
        :type list_dicFeuilBilanR: list of dict, optional
        :param list_invar: from l-egume, each element corresponds to an input specy of l-egume. Each element is a dict invar stores instant intern variables of l-egume., defaults to []
        :type list_invar: list of dict, optional
        :param id: list of indices from input scenes which corresponds to current l-egume instance. If you have several plantmodel among the input scenes, you need to precise which results you want to transfer to which instance of l-egume. defaults to None
        :type id: list, optional
        :raises AttributeError: you need to call :meth:`run` first
        :raises ValueError: unknown light model

        :return:

            if light model is RATP :func:`transfer_ratp_legume` :

                * ``res_abs_i``: absorbed energy in each voxels of a grid matching the dimensions of l-egume intern grid of voxels. One value for each specy.

                * ``res_trans``: transmitted energy in each voxels of a grid matching the dimensions of l-egume intern grid of voxels

            if light model is CARIBU :func:transfer_caribu_legume :

                * update of list_invar: updates the keys ``"parap"`` and ``"parip"`` for each specy. Cumulatative energy per plant

                * ``res_trans``: transmitted energy in each voxels of a grid matching the dimensions of l-egume intern grid of voxels

        :rtype: numpy.array
        """
        if not hasattr(self, "_LightVegeManager__elements_outputs"):
            raise AttributeError("No results yet, run a light modeling first")

        epsilon = 1e-14

        if self.__lightmodel == "ratp":
            from lightvegemanager.transfer import transfer_ratp_legume

            return transfer_ratp_legume(
                m_lais, energy, self.__complete_voxmesh, self.__voxels_outputs, self.__nb0, epsilon
            )

        elif self.__lightmodel == "caribu":
            from lightvegemanager.voxelsmesh import reduce_layers_from_trimesh
            from lightvegemanager.transfer import transfer_caribu_legume

            skylayer = reduce_layers_from_trimesh(
                self.__complete_trimesh,
                self.__pmax,
                self.__lightmodel_parameters["sensors"][1],
                self.__lightmodel_parameters["sensors"][2],
                self.__matching_ids,
                id,
            )

            return transfer_caribu_legume(
                energy,
                skylayer,
                id,
                self.__elements_outputs,
                self.__sensors_outputs,
                self.__lightmodel_parameters["sensors"][1],
                self.__lightmodel_parameters["sensors"][2],
                m_lais,
                list_invar,
                list_lstring,
                list_dicFeuilBilanR,
                self.__environment["infinite"],
                epsilon,
            )

        else:
            raise ValueError("Unknown light model (ratp or caribu)")

    ## EXTERN TOOLS ##
    def s5(self):
        """Creates inputs files for s5 and runs it
        s5 is an external tool made to analyse a set of triangles in order to use a grid of voxels.
        It also computes leaf angle distribution from the triangulation

        .. note:: All files are located in ``s5`` folder

        **Input files created**

            * fort.51: contains triangulation. Possibility to precise stem elements it is registered in the instance of LightVegeManager

            * s5.par: stores grid of voxels informations and number of entities

        **Output files created**

            * fort.60

                dimensions xy of the grid

                stats by specy

                    - total leaf area
                    - leaf area index
                    - global zenital leaf angle distribution
                    - global azimutal leaf angle distribution

                stats by voxels

                    - #specy | #ix | #iy | #iz (coordinate xyz of the voxel) | eaf area density
                    - local zenital leaf angle distribution
                    - local azimutal leaf angle distribution

            * leafarea: for each specy, for each voxel

                ix | iy | iz | #specy | LAD | zenital-distribution | azimutal-distribution

        **example**

        >>> myscene # a plantgl Scene
        >>> testofs5 = LightVegeManager() # create a instance
        >>> testofs5.build( geometry={ "scenes" : [myscene] } ) # build the geometry
        >>> testofs5.s5() # run of s5, creates input and output files

        """
        currentfolder = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        s5folder = os.path.join(currentfolder, os.path.normpath("s5"))
        fort51 = os.path.join(s5folder, os.path.normpath("fort.51"))
        s5par = os.path.join(s5folder, os.path.normpath("s5.par"))

        # écriture du fichier fort.51 contenant la triangulation
        f = open(fort51, "w")
        c_tr = 1
        for id, triangles in self.__complete_trimesh.items():
            for t in triangles:
                if self.__geometry["stems id"] is not None:
                    if tuple(self.__matching_ids[id]) in self.__geometry["stems id"]:
                        stem = "000"
                else:
                    stem = "001"
                label = (
                    str(self.__matching_ids[id][1] + 1) + str("%05i" % (self.__matching_ids[id][1] + 1)) + stem + "000"
                )

                f.write("p\t1\t%s\t3\t" % (label))
                for i in range(3):
                    f.write("%f\t%f\t%f\t" % (t[i][0], t[i][1], t[i][2]))
                f.write("\n")
                c_tr += 1
        f.close()

        # écriture du fichier s5.par contenant les informations de la grille
        f = open(s5par, "w")
        f.write("%i\t%i\t%i\t%i\t0\n" % (self.__complete_voxmesh.nent, 9, 9, self.__complete_voxmesh.njz))
        for i in range(self.__complete_voxmesh.njz):
            f.write("%f\t" % (self.__complete_voxmesh.dz[i]))
        f.write("\n")

        njx = self.__complete_voxmesh.njx
        njy = self.__complete_voxmesh.njy
        dx = self.__complete_voxmesh.dx
        dy = self.__complete_voxmesh.dy
        f.write("%f\t%i\t%f\t%f\t%i\t%f\n" % (njx * dx, njx, dx, njy * dy, njy, dy))
        f.close()

        # exécution de s5 dans un sous process
        subprocess.call(".\s5.exe", shell=True, cwd=s5folder)

        print("\n" + "--- Fin de s5.f")

    def s2v(self):
        """Creates inputs files for s2v and runs it
        s5 is an external tool made to analyse a set of triangles in order to use a grid of voxels.
        It also computes leaf angle distribution from the triangulation

        .. note:: All files are located in ``s2v`` folder

        **Input files created**

            * fort.51: stores triangulation. Possibility to precise stem elements it is registered in the instance of LightVegeManager

            * s2v.par: stores grid of voxels informations and number of entities

        **Output files created**

            * s2v.log

                logs about processing

                global statistics

                    - total leaf area per specy
                    - total leaf area
                    - leaf area index
                    - global zenital leaf angle distribution

            * s2v.can

                z layer where each triangle is located (and copy its vertices)

            * s2v.area

                triangle id | z layer | triangle area

            * out.dang: SAIL file

                line 1: global leaf area index for specy 1
                line 2: global zenital leaf angle distribution for specy 1

            * leafarea: SAIL file

                each line: 0 | idz | leaf area index for each slope class in current z layer | 0 | 0 | leaf area density on the layer

        **example**

        >>> myscene # a plantgl Scene
        >>> testofs2v = LightVegeManager() # create a instance
        >>> testofs2v.build( geometry={ "scenes" : [myscene] } ) # build the geometry
        >>> testofs2v.s2v() # run of s2v, creates input and output files

        """
        currentfolder = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        s2vfolder = os.path.join(currentfolder, os.path.normpath("s2v"))
        fort51 = os.path.join(s2vfolder, os.path.normpath("fort.51"))
        s2vpar = os.path.join(s2vfolder, os.path.normpath("s2v.par"))

        # écriture du fichier fort.51 contenant la triangulation
        f = open(fort51, "w")
        c_tr = 1
        for id, triangles in self.__complete_trimesh.items():
            for t in triangles:
                if self.__geometry["stems id"] is not None:
                    if tuple(self.__matching_ids[id]) in self.__geometry["stems id"]:
                        stem = "000"
                else:
                    stem = "001"
                label = (
                    str(self.__matching_ids[id][1] + 1) + str("%05i" % (self.__matching_ids[id][1] + 1)) + stem + "000"
                )

                f.write("p\t1\t%s\t3\t" % (label))
                for i in range(3):
                    f.write("%f\t%f\t%f\t" % (t[i][0], t[i][1], t[i][2]))
                f.write("\n")
                c_tr += 1
        f.close()

        f = open(s2vpar, "w")
        # ligne 1
        f.write("9 9 %i\n" % (self.__complete_voxmesh.njz))

        # ligne 2
        for i in range(self.__complete_voxmesh.njz):
            f.write("%f " % (self.__complete_voxmesh.dz[i]))
        f.write("\n")

        # ligne 3
        njx = self.__complete_voxmesh.njx
        njy = self.__complete_voxmesh.njy
        dx = self.__complete_voxmesh.dx
        dy = self.__complete_voxmesh.dy
        f.write("%f %i %f %f %i %f %i\n" % (njx * dx, njx, dx, njy * dy, njy, dy, self.__complete_voxmesh.nent))

        # ligne 4
        f.write(
            "%f %f %f\n" % (self.__complete_voxmesh.xorig, self.__complete_voxmesh.yorig, self.__complete_voxmesh.zorig)
        )

        f.close()

        # exécution de s2v dans un sous process
        subprocess.call(".\s2v++.exe", shell=True, cwd=s2vfolder)

        print("--- Fin de s2v.cpp")

    def to_VTK(
        self,
        lighting=False,
        path="",
        i=0,
        printtriangles=True,
        printvoxels=True,
        virtual_sensors=False,
        sun=False,
        sun_scale=2,
        sun_origin=(0, 0, 0),
        sun_center=True,
    ):
        """Writes a VTK from mesh(es) in ``self``

        .. warning:: If ``lighting=True``, the :meth:`run` must have been called before to have results dataframes.

        .. note:: if ``sun_center`` is False, ``sun_origin`` is the starting point of the line and it ends at ``sun_origin + sun_scale*sun.position``, if ``sun_center`` is True, ``sun_origin`` is the middle point of the line.

        :param lighting: for writing lighting information associated with each element, defaults to False
        :type lighting: bool, optional
        :param path: file and path names for the VTK file(s)
        :type path: string, optional
        :param i: time iteration of the current VTk file(s)
        :type i: int, optional
        :param printtriangles: write triangulation if one has been created in :func:build, defaults to True
        :type printtriangles: bool, optional
        :param printvoxels: write grid of voxels if one has been created in :func:build, defaults to True
        :type printvoxels: bool, optional
        :param virtual_sensors: write grid of virtual sensors if one has been created in :func:build, defaults to False
        :type virtual_sensors: bool, optional
        :param sun: write a line representing the sun position, defaults to False
        :type sun: bool, optional
        :param sun_scale: size of the line, defaults to 2
        :type sun_scale: int, optional
        :param sun_origin: starting or middle point of the line, defaults to (0,0,0)
        :type sun_origin: tuple, optional
        :param sun_center: if True orig is the middle of the line, otherwise it is the starting point, defaults to True
        :type sun_center: bool, optional
        :raises AttributeError: it needs to have lighting informations
        :raises AttributeError: it needs to have a grid of virtual sensors
        """
        if lighting:
            # deactivate lighting if if no lighting results
            if self.__lightmodel == "ratp" and (not hasattr(self, "_LightVegeManager__voxels_outputs")):
                print("--- VTK:  No light data, run the simulation")
                lighting = False
            elif self.__lightmodel == "caribu" and (not hasattr(self, "_LightVegeManager__triangles_outputs")):
                print("--- VTK:  No light data, run the simulation")
                lighting = False

            if sun and not hasattr(self, "_LightVegeManager__sun"):
                raise AttributeError("No results yet, run a light modeling first")

        if self.__lightmodel == "ratp" and printvoxels:
            from lightvegemanager.VTK import ratp_prepareVTK

            if lighting:
                filepath = path + "_voxels_light" + "_" + str(i) + ".vtk"
                datanames = ["ShadedPAR", "SunlitPAR", "ShadedArea", "SunlitArea", "PARa", "Intercepted", "Transmitted"]
                ratp_prepareVTK(self.__complete_voxmesh, filepath, datanames, self.__voxels_outputs)

            else:
                filepath = path + "_voxels_nolight" + "_" + str(i) + ".vtk"
                ratp_prepareVTK(self.__complete_voxmesh, filepath)

        if self.__matching_ids and printtriangles:
            data = []
            datanames = []
            if lighting:
                filepath = path + "_triangles_light" + "_" + str(i) + ".vtk"

                if self.__lightmodel == "ratp":
                    datanames = [
                        "ShadedPAR",
                        "SunlitPAR",
                        "ShadedArea",
                        "SunlitArea",
                        "PARa",
                        "Intercepted",
                        "Transmitted",
                    ]
                elif self.__lightmodel == "caribu":
                    datanames = []
                    for band in self.__lightmodel_parameters["caribu opt"].keys():
                        datanames.append(band + " Eabs")
                        datanames.append(band + " Ei")
                data = [list(self.__triangles_outputs[name].fillna(0.0)) for name in datanames]

            else:
                filepath = path + "_triangles_nolight" + "_" + str(i) + ".vtk"

            VTKtriangles(self.__complete_trimesh, data, datanames, filepath)

        if sun:
            from lightvegemanager.VTK import VTKline

            filepath = path + "_sun" + "_" + str(i) + ".vtk"
            if sun_center:
                start = tuple([a + b * sun_scale for a, b in zip(sun_origin, self.__sun)])
                end = tuple([a + b * -sun_scale for a, b in zip(sun_origin, self.__sun)])
            else:
                start = sun_origin
                end = tuple([a + b * sun_scale for a, b in zip(sun_origin, self.__sun)])

            VTKline(start, end, filepath)

        if virtual_sensors:
            if not self.__issensors:
                raise AttributeError("No virtual sensors in the inputs")
            else:
                filepath = path + "_virtualsensors" + "_" + str(i) + ".vtk"
                data = []
                datanames = []
                if lighting:
                    var = []

                    for id, triangles in self.__sensors_caribu.items():
                        for t in triangles:
                            var.append(self.__sensors_outputs["par"][id])
                    data.append(var)
                    datanames.append("intercepted")
                indice = []
                for id, triangles in self.__sensors_caribu.items():
                    for t in triangles:
                        indice.append(
                            self.__sensors_outputs_df[self.__sensors_outputs_df.Sensor == id]["VegetationType"].values[
                                0
                            ]
                        )
                data.append(indice)
                datanames.append("VegetationType")

                VTKtriangles(self.__sensors_caribu, data, datanames, filepath)

    def to_plantGL(self, lighting=False, printtriangles=True, printvoxels=False, virtual_sensors=False):
        """Return a plantGL Scene from mesh(es) in ``self``

        .. warning:: If ``lighting=True``, the :meth:`run` must have been called before to have results dataframes.

        :param lighting: for writing lighting information associated with each element, defaults to False
        :type lighting: bool, optional
        :param printtriangles: write triangulation if one has been created in :meth:`build`, defaults to True
        :type printtriangles: bool, optional
        :param printvoxels: write grid of voxels if one has been created in :meth:`build`, defaults to False
        :type printvoxels: bool, optional
        :param virtual_sensors: write grid of virtual sensors if one has been created in :func:build, defaults to False
        :type virtual_sensors: bool, optional
        :raises AttributeError:  AttributeError
        :return: Returns plantGL Scene with triangles and/or voxels from the final mesh
        :rtype: plantgl.Scene
        :return: Returns plantGL Scene with a grid of virtual sensors
        :rtype: plantgl.Scene
        """
        import openalea.plantgl.all as pgl

        if lighting:
            # deactivate lighting if if no lighting results
            if self.__lightmodel == "ratp" and (not hasattr(self, "_LightVegeManager__voxels_outputs")):
                print("--- VTK:  No light data, run the simulation")
                lighting = False
            elif self.__lightmodel == "caribu" and (not hasattr(self, "_LightVegeManager__triangles_outputs")):
                print("--- VTK:  No light data, run the simulation")
                lighting = False

        plantgl_voxscene = pgl.Scene()
        plantgl_triscene = pgl.Scene()
        if self.__lightmodel == "ratp" and printvoxels:
            from lightvegemanager.plantGL import ratpgrid_to_plantGLScene

            if lighting:
                if printtriangles:
                    transparency = 0.35
                    plantgl_voxscene = ratpgrid_to_plantGLScene(
                        self.__complete_voxmesh, transparency=transparency, plt_cmap="Greens"
                    )
                else:
                    transparency = 0.0
                    plantgl_voxscene = ratpgrid_to_plantGLScene(
                        self.__complete_voxmesh,
                        outputs=self.__voxels_outputs,
                        plt_cmap="seismic",
                        transparency=transparency,
                    )
            else:
                transparency = 0.0
                if printtriangles:
                    transparency = 0.35
                plantgl_voxscene = ratpgrid_to_plantGLScene(
                    self.__complete_voxmesh, transparency=transparency, plt_cmap="Greens"
                )
        if self.__matching_ids and printtriangles:
            if lighting:
                from lightvegemanager.plantGL import cscene_to_plantGLScene_light

                if self.__lightmodel == "caribu":
                    column_name = "par Ei"
                elif self.__lightmodel == "ratp":
                    column_name = "PARa"
                plantgl_triscene = cscene_to_plantGLScene_light(
                    self.__complete_trimesh, outputs=self.__triangles_outputs, column_name=column_name
                )
            else:
                from lightvegemanager.plantGL import cscene_to_plantGLScene_stems

                plantgl_triscene = cscene_to_plantGLScene_stems(
                    self.__complete_trimesh, stems_id=self.__geometry["stems id"], matching_ids=self.__matching_ids
                )

        plantgl_scene = pgl.Scene()
        for s in plantgl_triscene:
            plantgl_scene.add(s)
        for s in plantgl_voxscene:
            plantgl_scene.add(s)

        if virtual_sensors:
            if not self.__issensors:
                raise AttributeError("No virtual sensors in the inputs")
            else:
                if lighting:
                    var = [v for v in self.__sensors_outputs["par"].values()]

                    plt_cmap = "seismic"
                    minvalue = numpy.min(var)
                    maxvalue = numpy.max(var)
                    colormap = pgl.PglMaterialMap(minvalue, maxvalue, plt_cmap)

                    for i, s in enumerate(self.__sensors_plantgl):
                        s.appearance = colormap(var[i])
                return plantgl_scene, self.__sensors_plantgl

        return plantgl_scene

    ## GETTERS ##
    @property
    def riri5_transmitted_light(self):
        """transmitted results if light model is RiRi light

        :return: transmitted energy following a grid of voxels
        :rtype: numpy.array
        """
        return self.__riri5_transmitted_light

    @property
    def riri5_intercepted_light(self):
        """Intercepted results if light model is RiRi light,

        :return: intercepted energy following a grid of voxels, for each specy
        :rtype: numpy.array
        """
        return self.__riri5_intercepted_light

    @property
    def elements_outputs(self):
        """Lighting results aggregate by element

        :return: Column names can change depending on the light model. Commonly, there is element indice, its area and intercepted energy
        :rtype: pandas.Dataframe
        """
        return self.__elements_outputs

    @property
    def triangles_outputs(self):
        """Lighting results aggregate by triangle, if it has a triangulation in its inputs

        :return: .. seealso:: :mod:`outputs` for column names
        :rtype: pandas.Dataframe
        """
        return self.__triangles_outputs

    @property
    def voxels_outputs(self):
        """Lighting results aggregate by voxels, only with RATP as the selected light model

        :return: .. seealso:: :mod:`outputs` for column names
        :rtype: pandas.Dataframe
        """
        return self.__voxels_outputs

    def sensors_outputs(self, dataframe=False):
        """Lighting results aggregate by sensors.
        Only with CARIBU if you activated the virtual sensors option

        :return: The output format is the same as CARIBU output format. Dict with ``Eabs`` key for absorbed energy and ``Ei`` for incoming energy
        :rtype: dict
        """
        try:
            if dataframe:
                return self.__sensors_outputs_df
            else:
                return self.__sensors_outputs
        except AttributeError:
            return None

    @property
    def sun(self):
        """Return sun position of the last :func:run call

        :return: vector (x, y, z)
        :rtype: tuple
        """
        return self.__sun

    @property
    def soilenergy(self):
        """Return soil energy, only with CARIBU if soilmesh option is activated

        :return: The output format is the same as CARIBU output format. Dict with ``Eabs`` key for absorbed energy and ``Ei`` for incoming energy
        :rtype: dict
        """
        return self.__soilenergy

    @property
    def maxtrianglearea(self):
        """Returns the largest triangle of triangles mesh
        Computed in :meth:`build`

        :return: area the triangle
        :rtype: float
        """
        try:
            return self.__areamax
        except AttributeError:
            return 0.0

    @property
    def legume_empty_layers(self):
        """Returns number of empty layers between top of the canopy and number of z layers expected by l-egume

        :return: result of :func:`fill_ratpgrid_from_legumescene`
        :rtype: int
        """
        return self.__nb0

    @property
    def tesselationtime(self):
        """_summary_

        :return: _description_
        :rtype: _type_
        """
        try:
            return self.__tess_time
        except AttributeError:
            return 0.0

    @property
    def modelruntime(self):
        """Returns running time of the light model

        :return: time in s
        :rtype: float
        """
        try:
            return self.__time_runmodel
        except AttributeError:
            return 0.0

    @property
    def leafangledistribution(self):
        try:
            return self.__angle_distrib
        except AttributeError:
            return {}

    @property
    def domain(self):
        """Returns the largest triangle of triangles mesh
        Computed in :meth:`build`

        :return: area the triangle
        :rtype: float
        """
        try:
            return self.__domain
        except AttributeError:
            return 0.0
