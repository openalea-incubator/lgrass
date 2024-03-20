"""
    trianglesmesh
    **************

    It builds and handles triangulation mesh.

    The main triangulation format in LightVegeManager is the CARIBU format:
    scene (dict): a `{primitive_id: [triangle,]}` dict. A triangle is a list of 3-tuples points coordinates.
    example:
    
    .. code-block:: python
    
        scene = {   # element 0
                    0 : [ 
                            [(1,1,1), (0,0,0), (1,2,2)],  # triangle 0 of element 0
                            [(10,10,10), (0,0,0), (10,20,20)] # triangle 1 of element 0
                        ] ,
                    
                    # element 1
                    1 : [
                            [(12,12,12), (10,20,30), (21,22,22)], 
                            [(10,10,10), (10,20,30), (20,20,20)]
                        ]
                }

    .. seealso:: For more details :ref:`Inputs description <inputs>`

"""
import itertools
import numpy
import pandas
import bisect
import random
import math

from lightvegemanager.basicgeometry import triangle_area, rescale, translate, zrotate

def triangles_entity(cscene, entity_id, matching_ids):
    """Return a list of triangles belonging to the specy ``entity_id``

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :param entity_id: specy indice in matching_ids
    :type entity_id: int
    :param matching_ids:
        dict that matches new element indices in cscene with specy indice and
        input element indice,
        .. code:: matching_ids = { new_element_id : (input_element_id, specy_id)}

        this dict allows us to look how species there is the inputs geometric data

    :type matching_ids: dict
    :return: list of triangles belonging to ``entity_id`` in ``cscene``, a triangle is ``[(x1, y1, z1), (x2, y2, z2), (x3, y3, z3)]``
    :rtype: list
    """
    select_triangles = [t for k, t in cscene.items() if matching_ids[k][1] == entity_id]
    return list(itertools.chain(*select_triangles))


def globalid_to_elementid(cscene, triangleid):
    """Return the element index in cscene from a global triangle indice,
    as if triangles were in a list without a sorting by element.

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :param triangleid: indice from 0 to total number of triangles in cscene
    :type triangleid: int
    :raises IndexError: if triangleid > total number of triangles in cscene
    :return: element indice where triangleid belongs in cscene
    :rtype: int
    """
    cumul_triangles_per_ele = numpy.cumsum([len(v) for v in cscene.values()])
    if triangleid >= cumul_triangles_per_ele[-1]:
        raise IndexError("id %i > number of triangles: %i" % (triangleid, cumul_triangles_per_ele[-1]))

    # returns i as if cumul_triangles_per_ele[i-1] <= triangleid < cumul_triangles_per_ele[i]
    return bisect.bisect(cumul_triangles_per_ele.tolist(), triangleid)


def globalid_to_triangle(cscene, triangleid):
    """Return the corresponding triangle in cscene from a global triangle indice,
    as if triangles were in a list without a sorting by element.

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :param triangleid: indice from 0 to total number of triangles in cscene
    :type triangleid: int
    :raises IndexError: if triangleid > total number of triangles in cscene
    :return: Corresponding triangle in cscene
    :rtype: list of tuple
    """
    cumul_triangles_per_ele = numpy.cumsum([len(v) for v in cscene.values()])
    if triangleid >= cumul_triangles_per_ele[-1]:
        raise IndexError("id %i > number of triangles: %i" % (triangleid, cumul_triangles_per_ele[-1]))
    ele = globalid_to_elementid(cscene, triangleid)
    return cscene[ele][triangleid - sum(cumul_triangles_per_ele[0:ele])]


def compute_area_max(cscene):
    """Maximum triangle area in a triangulation

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :return: maximum triangle area in ``cscene``
    :rtype: float
    """
    amax = -999999
    for t in itertools.chain(*[v for v in cscene.values()]):
        a = triangle_area(t)
        if a > amax:
            amax = a
    return amax


def compute_minmax_coord(cscene):
    """Maximum and minimum point in all triangulation mesh

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :return: minimum point, maximum point
    :rtype: list, list
    """
    (xmax, ymax, zmax) = (-999999 for i in range(3))
    (xmin, ymin, zmin) = (999999 for i in range(3))
    for t in itertools.chain(*[v for v in cscene.values()]):
        for p in t:
            if p[0] > xmax:
                xmax = p[0]
            if p[0] < xmin:
                xmin = p[0]
            if p[1] > ymax:
                ymax = p[1]
            if p[1] < ymin:
                ymin = p[1]
            if p[2] > zmax:
                zmax = p[2]
            if p[2] < zmin:
                zmin = p[2]

    return [xmin, ymin, zmin], [xmax, ymax, zmax]


def compute_trilenght_max(cscene):
    """Search for maximum side of an axis oriented voxel where all triangles could fit in

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :return: maximum side of an axis oriented voxel where all triangles in cscene could fit in
    :rtype: float
    """
    lenmax = -999999
    # iterate on all triangles
    for t in itertools.chain(*[v for v in cscene.values()]):
        (xmax, ymax, zmax) = (-999999 for i in range(3))
        (xmin, ymin, zmin) = (999999 for i in range(3))
        for p in t:
            if p[0] > xmax:
                xmax = p[0]
            if p[0] < xmin:
                xmin = p[0]
            if p[1] > ymax:
                ymax = p[1]
            if p[1] < ymin:
                ymin = p[1]
            if p[2] > zmax:
                zmax = p[2]
            if p[2] < zmin:
                zmin = p[2]
        m = max(xmax - xmin, ymax - ymin, zmax - zmin)
        if m > lenmax:
            lenmax = m

    return lenmax

def isatriangle(t):
    if not isinstance(t, list) or len(t) != 3:
        return False
    
    for item in t:
        if not isinstance(item, (list, tuple)) or len(item) != 3:
            return False
        
        if not all(isinstance(val, float) for val in item):
            return False
    
    return True

def chain_triangulations(scenes):
    """Aggregates all input geometric scenes in one global triangulation

        Current known formats:

            * plantGL scene
            * VGX file
            * CARIBU triangulation format (which is used in LightVegeManager for the global scene)
            * MTG table with ``"geometry"`` identifier
            * l-egume grid: dict of two entries, leaf area through a voxels grid and leaf angle distribution, for each specy

        Only l-egume grid can stores multiple species, otherwise each entry scene must represent one specy

    :param scenes: ``geometric["scenes"]`` from LightVegeManager inputs, list of geometric scenes. Scenes can be in different format in the list.
    :type scenes: list
    :return:
        it returns 4 objects

            * ``complete_trimesh``: global triangulation which aggregates all input scenes. It is in CARIBU format
            * ``matching_ids``: dict which stores correspondances between new element indice, element indice in the input and specy indice
            * ``legume_grid``: boolean specifing if at least one l-egume grid is among the input scenes
            * ``id_legume_scene``: indice in the input scenes of a l-egume grid

    :rtype: dict of list, dict of list, bool, int
    """
    try:
        import openalea.plantgl.all as pgl
    except ImportError:
        pass
    try:
        from alinea.caribu import plantgl_adaptor
    except ImportError:
        pass
    try:
        from openalea.mtg.mtg import MTG
    except ImportError:
        pass

    # pre-check of scenes input, if it has only one triangle or one list of triangles
    if isatriangle(scenes) :
        scenes = [scenes]
    elif all(isatriangle(s) for s in scenes) and scenes != []:
        scenes = [scenes]
    
    complete_trimesh = {}
    matching_ids = {}
    legume_grid = False
    id_legume_scene = []

    element_count = 0
    for entity, scene in enumerate(scenes):
        cscene = {}

        # single triangle such as: scene = [(x1, y1, z1), (x2, y2, z2), (x3, y3, z3)]
        if isinstance(scene, list) :
            if isatriangle(scene) :
                cscene = { 0 : [ scene ] }
            elif isatriangle(scene[0]) :
                cscene = {0 : scene}

        # scene planteGL
        elif isinstance(scene, pgl.Scene):
            # keys/element indice are shape.id
            cscene = plantgl_adaptor.scene_to_cscene(scene)

        # fichier VGX
        elif isinstance(scene, str) and scene.split(".")[-1] == "vgx":
            # only one element id equal to element_count
            cscene = vgx_to_caribu(scene, element_count)

        # scene already in format CARIBU
        elif isinstance(scene, dict) and isinstance(list(scene.values())[0], list):
            cscene = scene

        # MTG table with "geometry" identifier
        elif isinstance(scene, MTG):
            # mtg["geometry"] is plantGL scene
            cscene = plantgl_adaptor.mtg_to_cscene(scene)

        # voxels grid in l-egume format
        elif isinstance(scene, dict) and len(scene) == 2:
            legume_grid = True
            id_legume_scene.append(entity)

        # creates a new numerotation for elements id adds them to the final triangulation
        for key, val in cscene.items():
            matching_ids[element_count] = [key, entity]
            val = [list(a) for a in val]  # converts triangles in list
            complete_trimesh[element_count] = val
            element_count += 1
    
    if not id_legume_scene :
        id_legume_scene = None
    
    return complete_trimesh, matching_ids, legume_grid, id_legume_scene


def vgx_to_caribu(file, id_element):
    """Reads VGX file and converts them in CARIBU scene format
    line is a leaf triangle if R column != 42
    :param file: path name of the file
    :type file: string
    :param id_element: element indice where the triangulation will be stored
    :type id_element: int
    :return: triangulation in CARIBU format ``{id_element : [triangle,]}``
    :rtype: dict of list
    """
    f = open(file, "r")
    lines = f.readlines()
    scene = {}
    triangles = []
    for l in lines[1:]:
        l_list = l.split("\t")

        # we consider as leaf elements line with R (RGB) != 42
        if l_list[10] != "42":
            tr = [
                (float(l_list[13]), float(l_list[14]), float(l_list[15])),
                (float(l_list[16]), float(l_list[17]), float(l_list[18])),
                (float(l_list[19]), float(l_list[20]), float(l_list[21])),
            ]

            if triangle_area(tr) > 0:
                triangles.append(tr)

    f.close()

    scene[id_element] = triangles
    return scene


def apply_transformations(cscene, matching_ids, transformations, cscene_unit):
    """Applies geometric transformations on some part of the triangulation set in input datas
    Each transformation applies on all triangles from the same input scene.

    :param cscene: LightVegeManager triangulations mesh
    :type cscene: dict of list
    :param matching_ids:
        dict that matches new element indices in cscene with specy indice and
        input element indice,

        .. code-block::

            matching_ids = { new_element_id : (input_element_id, specy_id)}

        this dict allows us to look how species there is the inputs geometric data

    :type matching_ids: dict
    :param transformations: dict containing which geometric to apply on which element

        Possible transformations:

        * ``["scenes unit] = {id_scene: unit}``
        * ``["rescale"] = {id_scene: h}``
        * ``["translate"] = {id_scene: vector}``
        * ``["xyz orientation"] = {id_scene: orientation}``

    ``id_scene`` is the index in the input geometric scenes list

    :type transformations: dict of dict
    :param cscene_unit: measure unit of cscene
    :type cscene_unit: string
    :raises ValueError: unit is not in the list ``['mm','cm','dm','m','dam','hm','km']``
    :raises ValueError: unit is not in the list ``['mm','cm','dm','m','dam','hm','km']``
    :raises ValueError: xyz orientation is not one of ``"x+ = S", "x+ = E", "x+ = W", "x+ = N"``
    :return: ``cscene`` with transformations applied to inputs scene parts in the global mesh
    :rtype: dict of list
    """
    # Rescale based on measure unit
    if "scenes unit" in transformations:
        units = {"mm": 0.001, "cm": 0.01, "dm": 0.1, "m": 1, "dam": 10, "hm": 100, "km": 1000}

        if cscene_unit not in units:
            raise ValueError(
                "Unknown final scene unit: select one in this \
                                    list ['mm','cm','dm','m','dam','hm','km']"
            )

        for key, scene_unit in transformations["scenes unit"].items():
            if scene_unit not in units:
                raise ValueError(
                    "Unknown scene %i unit: select one in this \
                            list ['mm','cm','dm','m','dam','hm','km']"
                    % (key)
                )

            if scene_unit != cscene_unit:
                h = units[scene_unit] / units[cscene_unit]
                for id in [k for k, v in matching_ids.items() if v[1] == key]:
                    cscene[id] = rescale(cscene[id], h)

    if "rescale" in transformations:
        for key, value in transformations["rescale"].items():
            for id in [k for k, v in matching_ids.items() if v[1] == key]:
                cscene[id] = rescale(cscene[id], value)

    if "translate" in transformations:
        for key, value in transformations["translate"].items():
            for id in [k for k, v in matching_ids.items() if v[1] == key]:
                cscene[id] = translate(cscene[id], value)

    # on ramène la scène à la convention x+ = N
    if "xyz orientation" in transformations:
        for key, value in transformations["xyz orientation"].items():
            if value == "x+ = S":
                rot = 180
            elif value == "x+ = W":
                rot = 90
            elif value == "x+ = E":
                rot = -90
            elif value == "x+ = N":
                rot = 0
            else:
                raise ValueError("Unknown xyz orientation in scene %i" % (key))
            for id in [k for k, v in matching_ids.items() if v[1] == key]:
                cscene[id] = zrotate(cscene[id], rot)


def create_heterogeneous_canopy(
    geometrical_model,
    mtg=None,
    nplants=50,
    var_plant_position=0.03,
    var_leaf_inclination=0.157,
    var_leaf_azimut=1.57,
    var_stem_azimut=0.157,
    plant_density=250,
    inter_row=0.15,
    id_type=None,
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
    from alinea.adel.Stand import AgronomicStand
    import openalea.plantgl.all as plantgl
    import random

    if seed is not None:
        random.seed(seed)
        numpy.random.seed(seed)

    # Load scene
    if not isinstance(geometrical_model, plantgl.Scene):
        initial_scene = geometrical_model.scene(mtg)
    else:
        initial_scene = geometrical_model

    alea_canopy = pandas.DataFrame()

    # Planter
    stand = AgronomicStand(
        sowing_density=plant_density, plant_density=plant_density, inter_row=inter_row, noise=var_plant_position
    )
    _, domain, positions, _ = stand.smart_stand(nplants=nplants, at=inter_row, convunit=1)

    random.seed(1234)

    # Built alea table if does not exist yet
    if alea_canopy.empty and mtg is not None:
        elements_vid_list = []
        for mtg_plant_vid in mtg.components_iter(mtg.root):
            for mtg_axis_vid in mtg.components_iter(mtg_plant_vid):
                for mtg_metamer_vid in mtg.components_iter(mtg_axis_vid):
                    for mtg_organ_vid in mtg.components_iter(mtg_metamer_vid):
                        for mtg_element_vid in mtg.components_iter(mtg_organ_vid):
                            if mtg.label(mtg_element_vid) == "LeafElement1":
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
    new_id = 0
    if id_type == "organ":
        new_id = max([sh.id for sh in initial_scene]) + 1

    for pos in positions:
        azimut_stem = random.uniform(-var_stem_azimut, var_stem_azimut)
        for shp in initial_scene:
            if mtg is not None:
                if mtg.label(shp.id) == "StemElement":
                    rotated_geometry = plantgl.EulerRotated(azimut_stem, 0, 0, shp.geometry)
                    translated_geometry = plantgl.Translated(plantgl.Vector3(pos), rotated_geometry)
                    if id_type == "organ":
                        new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=new_id)
                    elif id_type == "plant":
                        new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=position_number)
                    else:
                        new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=shp.id)
                    duplicated_scene += new_shape
                elif mtg.label(shp.id) == "LeafElement1":
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

                    if id_type == "organ":
                        new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=new_id)
                    elif id_type == "plant":
                        new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=position_number)
                    else:
                        new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=shp.id)
                    duplicated_scene += new_shape

            else:
                rotated_geometry = plantgl.EulerRotated(azimut_stem, 0, 0, shp.geometry)
                translated_geometry = plantgl.Translated(plantgl.Vector3(pos), rotated_geometry)
                if id_type == "organ":
                    new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=new_id)
                elif id_type == "plant":
                    new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=position_number)
                else:
                    new_shape = plantgl.Shape(translated_geometry, appearance=shp.appearance, id=shp.id)
                duplicated_scene += new_shape
            new_id += 1
        position_number += 1

    return duplicated_scene, domain

def random_triangle_generator(worldsize=(0,100), 
                                spheresize=(1.,1.), 
                                sigma_angle=(math.pi, math.pi), 
                                theta_angle=(math.pi/4, math.pi/5)):
    """Generate a random based on parameters
    
    Vertices are generated on a surface of a sphere

    Args:
        worldsize (tuple, optional): min and max where sphere center can be generated. Defaults to (0,100).
        spheresize (tuple, optional): mean and std of the sphere size. Defaults to (1.,1.).
        sigma_angle (tuple, optional): mean and std of the spherical angle on xy plane. Defaults to (math.pi, math.pi).
        theta_angle (tuple, optional): mean and std of the zenithal angle. Defaults to (math.pi/4, math.pi/5).

    Returns:
        list of 3 3-tuples: triangles defined by 3 xyz points
    """    
    r = random.gauss(spheresize[0], spheresize[1])
    x0, y0, z0 = [random.uniform(worldsize[0], worldsize[1]) for i in range(3)]
    triangle = []
    for i in range(3) :
        s = random.gauss(sigma_angle[0], sigma_angle[1])
        t = random.gauss(theta_angle[0], theta_angle[1])
        x = x0 + r * math.cos(s) * math.sin(t)
        y = y0 + r * math.sin(s) * math.sin(t)
        z = z0 + r * math.cos(t)
        triangle.append((x,y,z))
    return triangle