import os
import math

from lightvegemanager.LVM import LightVegeManager
from plantfusion.utils import create_child_folder
from plantfusion.planter import Planter
from plantfusion.indexer import Indexer


class Light_wrapper(object):
    def __init__(
        self,
        planter=Planter(),
        indexer=Indexer(),
        lightmodel="",
        sky="turtle46",
        direct=False,
        diffuse=True,
        reflected=False,
        coordinates=[46.4, 0.0, 1.0],
        infinite=True,
        legume_wrapper=None,
        caribu_opt={"par": (0.10, 0.07)},
        voxels_size=[1.0, 1.0, 1.0],
        angle_distrib_algo="compute global",
        nb_angle_class=9,
        mu=1.0,
        writegeo=False,
        out_folder="",
    ):
        self.transformations = planter.transformations
        self.indexer = indexer
        self.writegeo = writegeo
        self.compute_sensors = False
        self.direct=direct
        if writegeo:
            create_child_folder(os.path.normpath(out_folder), "light")
            self.out_folder = os.path.join(os.path.normpath(out_folder), "light")
            create_child_folder(self.out_folder, "vtk")
            create_child_folder(self.out_folder, "plantgl")

        self.lightmodel = lightmodel

        self.type_domain = planter.type_domain
        self.domain = planter.domain

        self.environment = {
            "coordinates": coordinates,
            "sky": sky,
            "direct": direct,
            "diffus": diffuse,
            "reflected": reflected,
            "infinite": infinite,
        }

        # calcul du nombre d'espèce
        self.number_of_species = len(indexer.global_order)

        # les instances de l-egume doivent suivre la même grille
        dxyz_legume = [0.0] * 3
        nxyz_legume = [0] * 3
        if legume_wrapper is not None:
            if isinstance(legume_wrapper, list):
                nxyz_legume = []
                for wrap in legume_wrapper:
                    nxyz_legume.append([
                        wrap.number_of_voxels()[3],
                        wrap.number_of_voxels()[2],
                        wrap.number_of_voxels()[1],
                    ])
                nxyz_legume = [max(*x) for x in zip(*nxyz_legume)]
                dxyz_legume = [x * 0.01 for x in legume_wrapper[0].voxels_size()]  # conversion de cm à m
            else:
                nxyz_legume = [
                    legume_wrapper.number_of_voxels()[3],
                    legume_wrapper.number_of_voxels()[2],
                    legume_wrapper.number_of_voxels()[1],
                ]
                dxyz_legume = [x * 0.01 for x in legume_wrapper.voxels_size()]  # conversion de cm à m

        lightmodel_parameters = {}

        if lightmodel == "caribu":
            lightmodel_parameters["caribu opt"] = caribu_opt
            lightmodel_parameters["sun algo"] = "caribu"
            if legume_wrapper is not None:
                self.compute_sensors = True
                x_trans, y_trans = 0, 0
                if isinstance(legume_wrapper, list):
                    lightmodel_parameters["sensors"] = {}
                    for wrap in legume_wrapper:
                        if "translate" in planter.transformations:
                            if wrap.global_index in planter.transformations["translate"]:
                                x_trans, y_trans, z = planter.transformations["translate"][wrap.global_index]
                        orig = [self.domain[0][0]+x_trans, self.domain[0][1]+y_trans, 0.0]
                        nxyz_legume = [
                            wrap.number_of_voxels()[3],
                            wrap.number_of_voxels()[2],
                            wrap.number_of_voxels()[1],
                        ]
                        dxyz_legume = [x * 0.01 for x in wrap.voxels_size()]  # conversion de cm à m
                        lightmodel_parameters["sensors"][wrap.global_index] = ["grid", dxyz_legume, nxyz_legume, orig]
                else:                
                    if "translate" in planter.transformations:
                        if isinstance(legume_wrapper.global_index, list):
                            if any([i in planter.transformations["translate"] for i in legume_wrapper.global_index]) :
                                x_trans, y_trans, z = planter.transformations["translate"][legume_wrapper.global_index[0]]
                        else:
                            if legume_wrapper.global_index in planter.transformations["translate"] :
                                x_trans, y_trans, z = planter.transformations["translate"][legume_wrapper.global_index]
                    orig = [self.domain[0][0]+x_trans, self.domain[0][1]+y_trans, 0.0]   
                    lightmodel_parameters["sensors"] = ["grid", dxyz_legume, nxyz_legume, orig]

            lightmodel_parameters["debug"] = False
            lightmodel_parameters["soil mesh"] = 1

        elif lightmodel == "ratp" or lightmodel == "riri5":
            if legume_wrapper is not None:
                lightmodel_parameters["voxel size"] = dxyz_legume
                lightmodel_parameters["origin"] = [0.0, 0.0, 0.0]
                lightmodel_parameters["full grid"] = True
            else:
                lightmodel_parameters["voxel size"] = voxels_size
                lightmodel_parameters["full grid"] = False

            lightmodel_parameters["mu"] = [mu] * self.number_of_species
            lightmodel_parameters["reflectance coefficients"] = [[0.0, 0.0]] * self.number_of_species

            if "/" in angle_distrib_algo or "\\" in angle_distrib_algo:
                lightmodel_parameters["angle distrib algo"] = "file"
                lightmodel_parameters["angle distrib file"] = angle_distrib_algo
            else:
                lightmodel_parameters["angle distrib algo"] = angle_distrib_algo
                lightmodel_parameters["nb angle classes"] = nb_angle_class
                lightmodel_parameters["soil reflectance"] = [0.0, 0.0]

        else:
            print("lightmodel not recognize")
            raise

        self.light = LightVegeManager(
            environment=self.environment,
            lightmodel=lightmodel,
            lightmodel_parameters=lightmodel_parameters,
            main_unit="m",
        )

        self.i_vtk = 0

    def run(self, energy=1.0, scenes=[], day=1, hour=12, parunit="RG", stems=None):
        geometry = {"scenes": scenes, "domain": self.domain, "transformations": self.transformations, "stems id": stems}

        self.light.build(geometry)
        self.light.run(
            energy=energy,
            day=day,
            hour=hour,
            truesolartime=True,
            parunit=parunit,
            id_sensors=self.indexer.legume_index,
        )

        if self.writegeo:
            file_project_name = os.path.join(os.path.normpath(self.out_folder), "vtk", "plantfusion_")
            
            if self.lightmodel == "ratp":
                printvoxels = True
            elif self.lightmodel == "caribu":
                printvoxels = True

            self.light.to_VTK(lighting=True, 
                                path=file_project_name, 
                                i=self.i_vtk, 
                                printtriangles=True, 
                                printvoxels=printvoxels, 
                                virtual_sensors=self.compute_sensors, 
                                sun=self.direct)
            scene_plantgl= self.light.to_plantGL(lighting=True, 
                                                    printtriangles=True, 
                                                    printvoxels=printvoxels, 
                                                    virtual_sensors=self.compute_sensors)
            
            if self.compute_sensors:
                scene_plantgl[0].save(
                    os.path.join(self.out_folder, "plantgl", "scene_light_plantgl_" + str(self.i_vtk)) + ".bgeom"
                )
                scene_plantgl[1].save(
                    os.path.join(self.out_folder, "plantgl", "sensors_plantgl_" + str(self.i_vtk)) + ".bgeom"
                )
            else:
                scene_plantgl.save(
                    os.path.join(self.out_folder, "plantgl", "scene_light_plantgl_" + str(self.i_vtk)) + ".bgeom"
                )

            self.i_vtk += 1

    def results_organs(self):
        try:
            return self.light.elements_outputs
        except AttributeError:
            return None

    def results_voxels(self):
        return self.light.voxels_outputs

    def results_triangles(self):
        return self.light.triangles_outputs

    def results_sensors(self):
        return self.light.sensors_outputs(dataframe=True)

    def res_trans(self):
        return self.light.riri5_transmitted_light

    def res_abs_i(self):
        return self.light.riri5_intercepted_light

    def soil_energy(self):
        try:
            return self.light.soilenergy["Qi"]
        except AttributeError:
            return -1

    def nb_empty_z_layers(self):
        return self.light.legume_empty_layers

    def xydomain_lightvegemanager(self):
        return self.light.domain
    
    def plantgl(self, lighting=False, printtriangles=True, printvoxels=False):
        return self.light.to_plantGL(lighting=lighting, 
                                        printtriangles=printtriangles, 
                                        printvoxels=printvoxels, 
                                        virtual_sensors=self.compute_sensors) 
