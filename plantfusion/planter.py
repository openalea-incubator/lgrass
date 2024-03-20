import numpy
import math
import pandas

from alinea.adel.Stand import AgronomicStand
import openalea.plantgl.all as plantgl
import random

from plantfusion.indexer import Indexer


class Planter:
    def __init__(
        self,
        generation_type="default",
        indexer=Indexer(),
        legume_cote={},
        legume_number_of_plants={},
        inter_rows=0.15,
        plant_density={1: 250},
        xy_plane=None,
        xy_square_length=0.5,
        translate=None,
        noise_plant_positions=0.0,
        save_wheat_positions=False,
        seed=None,
    ) -> None:
        self.generation_type = generation_type
        self.plant_density = plant_density
        self.save_wheat_positions = save_wheat_positions
        self.noise_plant_positions = noise_plant_positions
        self.indexer = indexer
        self.number_of_plants: list = [0 for i in indexer.global_order]

        # les lsystem l-egume sont par défaut en cm et le reste en m
        self.transformations: dict = {"scenes unit": {}}
        for i in range(len(self.indexer.global_order)):
            self.transformations["scenes unit"][i] = "m"
        for i in self.indexer.legume_index:
            self.transformations["scenes unit"][i] = "cm"

        if generation_type == "default":
            self.__default_preconfigured(legume_cote, inter_rows, plant_density, xy_plane, translate, seed)
            for name, nb_plt in legume_number_of_plants.items():
                self.number_of_plants[indexer.global_order.index(name)] = nb_plt

        elif generation_type == "random":
            self.__random(plant_density, xy_square_length)
            self.type_domain = "mix"

        elif generation_type == "row":
            self.__row(plant_density, inter_rows)
            self.type_domain = "mix"

    def __random(self, plant_density, xy_square_length):
        # fait un carré à partir de (0,0)
        # xy_plane = longueur d'un cote du carré du sol
        self.legume_nbcote = []
        self.wheat_positions = []
        self.other_positions = []

        self.domain = ((0.0, 0.0), (xy_square_length, xy_square_length))

        for name, density in plant_density.items():
            if name in self.indexer.legume_names:
                # on réajuste le domaine pour avoir 64 plantes
                xy_square_length = math.sqrt(64 / density)
                self.legume_nbcote.append(8)
                self.domain = ((0.0, 0.0), (xy_square_length, xy_square_length))

                self.legume_typearrangement = "random8"
                # conversion m en cm
                self.legume_cote = xy_square_length * 100

                self.legume_optdamier = 8

                self.number_of_plants[self.indexer.global_order.index(name)] = 64

            else:
                self.number_of_plants[self.indexer.global_order.index(name)] = int(
                    xy_square_length * xy_square_length * density
                )

    def __row(self, plant_density, inter_rows):
        self.total_n_rows = 2 * len(self.indexer.global_order)

        xy_square_length = inter_rows * self.total_n_rows
        self.domain = ((0.0, 0.0), (xy_square_length, xy_square_length))

        self.inter_rows = inter_rows
        self.legume_nbcote = []
        self.wheat_positions = []
        self.other_positions = []

        for name, density in plant_density.items():
            if name in self.indexer.legume_names:
                if self.indexer.legume_number_of_species[self.indexer.legume_names.index(name)] == 1:
                    self.legume_typearrangement = "row4_sp1"
                elif self.indexer.legume_number_of_species[self.indexer.legume_names.index(name)] == 2:
                    self.legume_typearrangement = "row4"

                # conversion m en cm
                self.legume_cote = inter_rows * self.total_n_rows * 100
                self.legume_nbcote.append(int(xy_square_length * xy_square_length * density / 2))
                self.legume_optdamier = 2

                # if self.indexer.legume_number_of_species[self.indexer.legume_names.index(name)] > 1:
                indices = [index for index, item in enumerate(self.indexer.global_order) if item == name]
                for i in indices:
                    self.number_of_plants[i] = int(xy_square_length * xy_square_length * density)

            else:
                self.number_of_plants[self.indexer.global_order.index(name)] = int(
                    xy_square_length * xy_square_length * density
                )

        self.transformations["translate"] = {}
        if self.total_n_rows > 4:
            for i in self.indexer.legume_index:
                self.transformations["translate"][i] = (0.0, (i) * inter_rows, 0.0)
            for i in self.indexer.wheat_index:
                self.transformations["translate"][i] = (0.0, (i - 0.5) * inter_rows, 0.0)
            for i in self.indexer.other_index:
                self.transformations["translate"][i] = (0.0, (i - 0.5) * inter_rows, 0.0)

        # il y a que deux espèces
        else:
            # 2 wheats ou combinaison avec 1 ou 2 autres FSPM
            if len(self.indexer.wheat_names) > 1  or len(self.indexer.other_names) > 1:
                self.transformations["translate"][0] = (0.0, -inter_rows, 0.0)

            # 2 legume
            elif len(self.indexer.legume_names) > 1:
                self.transformations["translate"][1] = (0.0, inter_rows, 0.0)

    def __default_preconfigured(
        self, legume_cote={}, inter_rows=0.15, plant_density={1: 250}, xy_plane=None, translate=None, seed=None
    ):
        self.plant_density = plant_density
        self.inter_rows = inter_rows

        for i in self.indexer.wheat_index:
            self.number_of_plants[i] = 50
        if self.indexer.wheat_active:
            self.type_domain = "create_heterogeneous_canopy"

            # on vient récupérer le domain de AgronomicStand
            if seed is not None:
                random.seed(seed)
                numpy.random.seed(seed)

            stand = AgronomicStand(
                sowing_density=self.plant_density[1],
                plant_density=self.plant_density[1],
                inter_row=self.inter_rows,
                noise=self.noise_plant_positions,
            )
            _, domain, _, _ = stand.smart_stand(nplants=50, at=self.inter_rows, convunit=1)
            self.domain = domain

        # temporaire : on applique les translations que aux instances l-egume
        # xy_translate : dict {"name isntance" : (x,y,z)}
        if translate is not None:
            self.transformations["translate"] = {}
            for name, vector in translate.items():
                id = self.indexer.global_order.index(name)
                self.transformations["translate"][id] = vector

        # gestion du domain xy de définition
        if xy_plane is None:
            # si recalcul le domain via create_heterogeneous_canopy
            if not self.indexer.legume_active:
                self.type_domain = "create_heterogeneous_canopy"
            elif not self.indexer.wheat_active:
                self.type_domain = "l-egume"
                # convertit domain cm en m
                domains = []
                for name, cote in legume_cote.items():
                    vector = (0, 0, 0)
                    if translate is not None and name in translate :
                        vector = translate[name]
                    domains.append([
                        (0.0 + vector[0], 0.0  + vector[1]),
                        (cote * 0.01  + vector[0], cote * 0.01  + vector[1]),
                    ])
                self.domain = ((min([x[0][0] for x in domains]), min([x[0][1] for x in domains])), 
                            ((max([x[1][0] for x in domains]), max([x[1][1] for x in domains]))))

            else:
                self.type_domain = "mix"
                if self.indexer.legume_active:
                    domains = []
                    for name, cote in legume_cote.items():
                        vector = (0, 0, 0)
                        if translate is not None and name in translate :
                            vector = translate[name]
                        domains.append([
                            (0.0 + vector[0], 0.0  + vector[1]),
                            (cote * 0.01  + vector[0], cote * 0.01  + vector[1]),
                        ])
                    legume_domain = ((min([x[0][0] for x in domains]), min([x[0][1] for x in domains])), 
                                ((max([x[1][0] for x in domains]), max([x[1][1] for x in domains]))))

                # a été calculé au-dessus à  l'appel de create_heterogeneous_canopy
                wheat_domain = self.domain
                if translate is not None:
                    if self.indexer.wheat_names[0] in translate :
                        vector = translate[self.indexer.wheat_names[0]]
                        wheat_domain = (
                            (wheat_domain[0][0] + vector[0], wheat_domain[0][1] + vector[1]),
                            (wheat_domain[1][0] + vector[0], wheat_domain[1][1] + vector[1]),
                        )
                self.domain = (
                    (min(legume_domain[0][0], wheat_domain[0][0]), min(legume_domain[0][1], wheat_domain[0][1])),
                    (max(legume_domain[1][0], wheat_domain[1][0]), max(legume_domain[1][1], wheat_domain[1][1])),
                )
        else:
            self.type_domain = "input"
            self.domain = xy_plane


    def generate_random_other(self, indice_instance=0, seed=None):
        if seed is not None:
            s = seed
        else:
            s = 1234
        random.seed(s)
        numpy.random.seed(s)

        # tirage des positions
        # list de 3-tuple des positions
        if self.other_positions != [] :
            positions = self.other_positions
        else:
            positions = []
            for i in range(self.number_of_plants[self.indexer.other_index[indice_instance]]):
                positions.append(
                    (numpy.random.uniform(0.0, self.domain[1][0]), numpy.random.uniform(0.0, self.domain[1][0]), 0.0)
                )

        self.other_positions = positions

        return positions
    
    def generate_random_wheat(
        self, adel_wheat, mtg, indice_wheat_instance=0, stem_name="stem", leaf_name="leaf", seed=None
    ):
        var_leaf_inclination = 0.157
        var_leaf_azimut = 1.57
        var_stem_azimut = 0.157

        if seed is not None:
            s = seed
        else:
            s = 1234
        random.seed(s)
        numpy.random.seed(s)

        initial_scene = adel_wheat.scene(mtg)

        # tirage des positions
        # list de 3-tuple des positions
        if self.wheat_positions != [] and self.save_wheat_positions:
            positions = self.wheat_positions
        else:
            positions = []
            for i in range(self.number_of_plants[self.indexer.wheat_index[indice_wheat_instance]]):
                positions.append(
                    (numpy.random.uniform(0.0, self.domain[1][0]), numpy.random.uniform(0.0, self.domain[1][0]), 0.0)
                )

        self.wheat_positions = positions

        generated_scene = self.__generate_wheat_from_positions(
            initial_scene, mtg, positions, var_leaf_inclination, var_leaf_azimut, var_stem_azimut, stem_name, leaf_name
        )

        return generated_scene

    def generate_row_wheat(self, indice_instance=0, seed=None):

        if seed is not None:
            s = seed
        else:
            s = 1234
        random.seed(s)
        numpy.random.seed(s)

        if self.other_positions != []:
            positions = self.other_positions
        else:
            positions = []

            inter_plants = (
                2 * self.domain[1][1] / self.number_of_plants[self.indexer.wheat_index[indice_instance]]
            )
            nrows = 2
            self.total_n_rows

            # first row on left 1/2 interrow, then 1 out of 2 row is wheat
            rows_y = [self.inter_rows * 1.5, ((self.total_n_rows / nrows) + 1.5) * self.inter_rows]
            for y in rows_y:
                for ix in range(int(self.number_of_plants[self.indexer.wheat_index[indice_instance]] / nrows)):
                    x = inter_plants * (0.5 + ix)
                    p = (
                        random.uniform(x - self.noise_plant_positions, x + self.noise_plant_positions),
                        random.uniform(y - self.noise_plant_positions, y + self.noise_plant_positions),
                        0.0,
                    )
                    positions.append(p)

        self.other_positions = positions

        return positions

    def generate_row_wheat(
        self, adel_wheat, mtg, indice_wheat_instance=0, stem_name="stem", leaf_name="leaf", seed=None
    ):
        var_leaf_inclination = 0.157
        var_leaf_azimut = 1.57
        var_stem_azimut = 0.157

        if seed is not None:
            s = seed
        else:
            s = 1234
        random.seed(s)
        numpy.random.seed(s)

        initial_scene = adel_wheat.scene(mtg)

        if self.wheat_positions != [] and self.save_wheat_positions:
            positions = self.wheat_positions
        else:
            positions = []

            inter_plants = (
                2 * self.domain[1][1] / self.number_of_plants[self.indexer.wheat_index[indice_wheat_instance]]
            )
            nrows = 2
            self.total_n_rows

            # first row on left 1/2 interrow, then 1 out of 2 row is wheat
            rows_y = [self.inter_rows * 1.5, ((self.total_n_rows / nrows) + 1.5) * self.inter_rows]
            for y in rows_y:
                for ix in range(int(self.number_of_plants[self.indexer.wheat_index[indice_wheat_instance]] / nrows)):
                    x = inter_plants * (0.5 + ix)
                    p = (
                        random.uniform(x - self.noise_plant_positions, x + self.noise_plant_positions),
                        random.uniform(y - self.noise_plant_positions, y + self.noise_plant_positions),
                        0.0,
                    )
                    positions.append(p)

        self.wheat_positions = positions

        generated_scene = self.__generate_wheat_from_positions(
            initial_scene,
            mtg,
            positions,
            var_leaf_inclination,
            var_leaf_azimut,
            var_stem_azimut,
            stem_name=stem_name,
            leaf_name=leaf_name,
        )

        return generated_scene

    def create_heterogeneous_canopy(
        self,
        geometrical_model,
        mtg=None,
        var_leaf_inclination=0.157,
        var_leaf_azimut=1.57,
        var_stem_azimut=0.157,
        stem_name="stem",
        leaf_name="leaf",
        indice_wheat_instance=0,
        seed=None,
    ):
        """
        Duplicate a plant in order to obtain a heterogeneous canopy.

        :param int nplants: the desired number of duplicated plants
        :param float var_plant_position: variability for plant position (m)
        :param float var_leaf_inclination: variability for leaf inclination (rad)
        :param float var_leaf_azimut: variability for leaf azimut (rad)
        :param float var_stem_azimut: variability for stem azimut (rad)
        :param string id_type: precise how to set the shape id of the elements : None, plant or organ

        :return: duplicated heterogenous scene and its domain
        :rtype: openalea.plantgl.all.Scene, (float)
        """
        if seed is not None:
            random.seed(seed)
            numpy.random.seed(seed)

        # Planter
        stand = AgronomicStand(
            sowing_density=self.plant_density[1],
            plant_density=self.plant_density[1],
            inter_row=self.inter_rows,
            noise=self.noise_plant_positions,
        )
        _, domain, positions, _ = stand.smart_stand(
            nplants=self.number_of_plants[self.indexer.wheat_index[indice_wheat_instance]],
            at=self.inter_rows,
            convunit=1,
        )
        self.wheat_positions = positions

        random.seed(1234)

        generated_scene = self.__generate_wheat_from_positions(
            geometrical_model,
            mtg,
            positions,
            var_leaf_inclination,
            var_leaf_azimut,
            var_stem_azimut,
            stem_name=stem_name,
            leaf_name=leaf_name,
        )

        if self.type_domain == "create_heterogeneous_canopy":
            self.domain = domain

        return generated_scene

    def __generate_wheat_from_positions(
        self,
        geometrical_model,
        mtg=None,
        positions=[(0.0, 0.0, 0.0)],
        var_leaf_inclination=0.157,
        var_leaf_azimut=1.57,
        var_stem_azimut=0.157,
        stem_name="stem",
        leaf_name="leaf",
    ):
        # Load scene
        if not isinstance(geometrical_model, plantgl.Scene):
            initial_scene = geometrical_model.scene(mtg)
        else:
            initial_scene = geometrical_model

        alea_canopy = pandas.DataFrame()

        # Built alea table if does not exist yet
        if alea_canopy.empty and mtg is not None:
            elements_vid_list = []
            for mtg_plant_vid in mtg.components_iter(mtg.root):
                for mtg_axis_vid in mtg.components_iter(mtg_plant_vid):
                    for mtg_metamer_vid in mtg.components_iter(mtg_axis_vid):
                        for mtg_organ_vid in mtg.components_iter(mtg_metamer_vid):
                            for mtg_element_vid in mtg.components_iter(mtg_organ_vid):
                                if mtg.label(mtg_element_vid) == leaf_name:
                                    elements_vid_list.append(mtg_element_vid)

            elements_vid_df = pandas.DataFrame({"vid": elements_vid_list, "tmp": 1})
            positions_df = pandas.DataFrame(
                {"pos": range(len(positions)), "tmp": 1, "azimut_leaf": 0, "inclination_leaf": 0}
            )
            alea = pandas.merge(elements_vid_df, positions_df, on=["tmp"])
            alea = alea.drop("tmp", axis=1)
            for vid in elements_vid_list:
                numpy.random.seed(vid)
                alea.loc[alea["vid"] == vid, "azimut_leaf"] = numpy.random.uniform(
                    -var_leaf_azimut, var_leaf_azimut, size=len(positions)
                )
                alea.loc[alea["vid"] == vid, "inclination_leaf"] = numpy.random.uniform(
                    -var_leaf_inclination, var_leaf_inclination, size=len(positions)
                )
            alea_canopy = alea

        # Duplication and heterogeneity
        duplicated_scene = plantgl.Scene()
        position_number = 0
        for pos in positions:
            azimut_stem = random.uniform(-var_stem_azimut, var_stem_azimut)
            for shp in initial_scene:
                if mtg.label(shp.id) == stem_name:
                    rotated_geometry = plantgl.EulerRotated(azimut_stem, 0, 0, shp.geometry)
                    translated_geometry = plantgl.Translated(plantgl.Vector3(pos), rotated_geometry)
                    new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=shp.id)
                    duplicated_scene += new_shape
                elif mtg.label(shp.id) == leaf_name:
                    # Add shp.id in alea_canopy if not in yet:
                    if shp.id not in list(alea_canopy["vid"]):
                        new_vid_df = pandas.DataFrame({"vid": shp.id, "pos": range(len(positions))})
                        numpy.random.seed(shp.id)
                        new_vid_df["azimut_leaf"] = numpy.random.uniform(
                            -var_leaf_azimut, var_leaf_azimut, size=len(positions)
                        )
                        new_vid_df["inclination_leaf"] = numpy.random.uniform(
                            -var_leaf_inclination, var_leaf_inclination, size=len(positions)
                        )
                        alea_canopy = alea_canopy.copy().append(new_vid_df, sort=False)
                    # Translation to origin
                    anchor_point = mtg.get_vertex_property(shp.id)["anchor_point"]
                    trans_to_origin = plantgl.Translated(-anchor_point, shp.geometry)
                    # Rotation variability
                    azimut = alea_canopy.loc[
                        (alea_canopy.pos == position_number) & (alea_canopy.vid == shp.id), "azimut_leaf"
                    ].values[0]
                    inclination = alea_canopy.loc[
                        (alea_canopy.pos == position_number) & (alea_canopy.vid == shp.id), "inclination_leaf"
                    ].values[0]
                    rotated_geometry = plantgl.EulerRotated(azimut, inclination, 0, trans_to_origin)
                    # Restore leaf base at initial anchor point
                    translated_geometry = plantgl.Translated(anchor_point, rotated_geometry)
                    # Translate leaf to new plant position
                    translated_geometry = plantgl.Translated(pos, translated_geometry)
                    new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=shp.id)
                    duplicated_scene += new_shape

            position_number += 1

        return duplicated_scene
