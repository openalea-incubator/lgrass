'''
    CARIBUinputs
    *************

    Manages and prepares input information for CARIBU.

    In this module we use two of the LightVegeManager's inputs dict.

    *  ``parameters`` corresponding to CARIBU parameters and contains

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
    
    * ``geometry`` corresponding to the geometric information with the scenes inputs

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
    
    .. seealso:: For more details :ref:`Inputs description <inputs>`

'''

def Prepare_CARIBU(trimesh, 
                    geometry, 
                    matching_ids, 
                    minmax, 
                    parameters,
                    infinite, 
                    idsensors) :
    """Format optical parameters and prepare virtual sensors and debug if activated

    :param trimesh: triangles mesh aggregated by indice elements
        
        .. code-block:: 
        
            { id : [triangle1, triangle2, ...]}

    :type trimesh: dict
    :param geometry: geometric parameters from inputs of LightVegeManager
    :type geometry: dict
    :param matching_ids: 
        dict that matches new element indices in trimesh with specy indice and
        input element indice, 
        
        .. code-block:: 
        
            matching_ids = { new_element_id : (input_element_id, specy_id)}

    :type matching_ids: dict
    :param minmax: list of mininuml point and maximum point in a triangles mesh
    :type minmax: [3-tuple, 3-tuple]
    :param parameters: RATP parameters from inputs of LightVegeManager
    :type parameters: dict
    :param infinite: if the user wishes to activate infinitisation of the grid
    :type infinite: bool
    :param idsensors: 
        Sets which input scenes the grid of virtual sensors will match. 
        Elements of this list refers to indices in the list ``geometry["scenes"]``.
        If no input indices are given, the grid will match the global scene ``trimesh``. 
        Otherwise, you can match the grid to a specific set of input scenes
    :type idsensors: list of int
    
    :return:
        
        *  ``opt``: optical parameters formatted for CARIBU. It takes the form of dict, where each 
            key is a bandlength (PAR, NIR etc...) and values are transmittance and reflectance
            (or only reflectance for stems elements)
        
        *  ``sensors_caribu``: geometric data of the virtual sensors in CARIBU scene format. Sensors
            are horizontal square made of two triangles
            
            .. code-block:: 
                
                sensors_caribu = { sensor_id : [triangle1, triangle2], ...}

        * ``debug``: boolean if user wants to activate th edebug option in CARIBU

        * ``matching_sensors_species``: dict where each key is a sensor ID and value, its matching input vegetation type

    :rtype: dict, dict, bool, dict
    """                    
    # manage stems element and creates optical parameters
    stems = []
    if 'stems id' in geometry : stems = geometry["stems id"]
    opt = CARIBU_opticals(matching_ids, 
                            parameters, 
                            stems)

    # debug, sensors, domain
    sensors_caribu = None
    matching_sensors_species = {}
    if "sensors" in parameters and isinstance(parameters["sensors"], list):
        if "sensors" in parameters and parameters["sensors"][0] == "grid" :
            dxyz = parameters["sensors"][1]
            nxyz = parameters["sensors"][2]
            orig = parameters["sensors"][3]
            arg = (dxyz,
                    nxyz,
                    orig,
                    minmax[1],
                    trimesh,
                    matching_ids,
                    idsensors,
                    infinite)
            sensors_caribu, \
            sensors_plantgl, \
            Pmax_capt = create_caribu_legume_sensors(*arg)
            for i in sensors_caribu.keys():
                matching_sensors_species[i] = 0
    elif "sensors" in parameters and isinstance(parameters["sensors"], dict):
        sensors_caribu = {}
        matching_sensors_species = {}
        start_id = 0
        for specy_indice, sensors_parameters in parameters["sensors"].items():
            dxyz = sensors_parameters[1]
            nxyz = sensors_parameters[2]
            orig = sensors_parameters[3]
            arg = (dxyz,
                    nxyz,
                    orig,
                    minmax[1],
                    trimesh,
                    matching_ids,
                    idsensors,
                    infinite,
                    start_id)
            sensors_dict, \
            sensors_plantgl, \
            Pmax_capt = create_caribu_legume_sensors(*arg)
            sensors_caribu.update(sensors_dict)
            start_id = max(sensors_dict.keys()) + 1

            for key in sensors_dict.keys():
                matching_sensors_species[key] = specy_indice


    # infinite scene pattern if not precised in the inputs
    if "domain" not in geometry :
        # special case for l-egume and virtual sensors
        if "sensors" in parameters and parameters["sensors"][0] == "grid" :
            reduction_domain = 0
            d = ((orig[0] + reduction_domain, 
                    orig[1] + reduction_domain), 
                        (nxyz[0] * dxyz[0] - reduction_domain, 
                            nxyz[1] * dxyz[1] - reduction_domain))
        else:
            d = ((minmax[0][0], minmax[0][1]), 
                    (minmax[1][0], minmax[1][1]))

        geometry["domain"] = d

    debug = False
    if "debug" in parameters and parameters["debug"] :  debug = True
    
    return opt, sensors_caribu, debug, matching_sensors_species
    

def CARIBU_opticals(matching_ids, parameters, stems_id=None) :
    """Sets optical parameters for CARIBU from LightVegeManager inputs
    It removes transmittance for stems elements if precised

    :param matching_ids: 
        dict that matches new element indices in trimesh with specy indice and
        input element indice, :code:`matching_ids = { new_element_id : (input_element_id, specy_id)}`
    :type matching_ids: dict
    :param parameters: RATP parameters from inputs of LightVegeManager
    :type parameters: dict
    :param stems_id: list of potential stems element in the input scenes, defaults to None
    :type stems_id: list of 2-tuple, optional
    :return: 
        optical parameters formatted for CARIBU. It takes the form of dict, where each 
        key is a bandlength (PAR, NIR etc...) and values are transmittance and reflectance
        (or only reflectance for stems elements)
    :rtype: dict
    """    
    opt = {}
    for band, coef in parameters["caribu opt"].items() :
        opt[band] = {}
    for id, val in matching_ids.items():
        for band, coef in parameters["caribu opt"].items() :
            # id is element, stems_id is a list of (element id, specy id)
            # if stem then no transmittance, the object is opaque
            if stems_id is not None and (val[0], val[1]) in stems_id:
                #: (reflectance,) of the stems ou organe opaque
                opt[band][id] = (coef[0],)  
            else:
                # A 2-tuple encode a symmetric translucent material defined
                # by a reflectance and a transmittance
                # A 4-tuple encode an asymmetric translucent material defined
                # the reflectance and transmittance
                # of the upper and lower side respectively 
                #: (reflectance, transmittance) of the adaxial side of the 
                # leaves, élément translucide symétrique
                opt[band][id] = coef

    return opt

def create_caribu_legume_sensors(dxyz, 
                                nxyz, 
                                orig, 
                                pmax, 
                                trimesh, 
                                matching_ids, 
                                id_sensors, 
                                infinite,
                                start_id=0) :
    """Creates a set of virtual sensors following a voxels grid
    each sensor is a square made by two triangles and takes place on the bottom face of a voxel
    The grid follow the xyz axis (and so the voxels)

    :param dxyz: [dx, dy, dz] size of a voxel in each direction
    :type dxyz: list
    :param nxyz: [nx, ny, nz] number of voxels on each axis
    :type nxyz: list
    :param orig: [x_origin, y_origin, z_origin], origin of the grid, cartesian coordinates
    :type orig: list (or tuple)
    :param pmax: [x, y, z] maximum point of trimesh in the xyz space
    :type pmax: lsit (or tuple)
    :param trimesh: triangles mesh aggregated by indice elements
        
        .. code-block:: 
        
            { id : [triangle1, triangle2, ...]}

    :type trimesh: dict
    :param matching_ids: 
        dict that matches new element indices in trimesh with specy indice and
        input element indice, :code:`matching_ids = { new_element_id : (input_element_id, specy_id)}`
    :type matching_ids: dict
    :param idsensors: 
        Sets which input scenes the grid of virtual sensors will match. 
        Elements of this list refers to indices in the list ``geometry["scenes"]``.
        If no input indices are given, the grid will match the global scene ``trimesh``. 
        Otherwise, you can match the grid to a specific set of input scenes
    :type idsensors: list of int
    :param infinite: if the user wishes to activate infinitisation of the grid
    :type infinite: bool
    :return: 
        it returns 3 objects:

            * ``sensors_caribu``, geometric data of the virtual sensors in CARIBU scene format. Sensors
                are horizontal square made of two triangles :code:`sensors_caribu = { sensor_id : [triangle1, triangle2], ...}`

            * ``s_capt``, PlantGL scene of the virtual sensors

            * ``sensors_maxcenter``, [x, y, z] representing point on the highest sensors layer and middle of
                xy plane in the grid

    :rtype: dict, PlantGL scene, list
    """
    import openalea.plantgl.all as pgl
    from lightvegemanager.voxelsmesh import reduce_layers_from_trimesh

    # if number of filled layers is less than number of expected layers (nxyz[2])
    skylayer = reduce_layers_from_trimesh(trimesh, 
                                            pmax,
                                            dxyz,
                                            nxyz,
                                            matching_ids, 
                                            id_sensors)
            
    # initialize the PlantGL scene
    s_capt = pgl.Scene()
    
    # point to output
    sensors_maxcenter = [ orig[0] + 0.5*dxyz[0]*nxyz[0] , 
                            orig[1] + 0.5*dxyz[1]*nxyz[1] ,  
                            orig[2] + dxyz[2] * (nxyz[2] - skylayer -1)]

    # orientation changes depending if the scene is infinite or not
    if infinite : 
        # square looks upward
        points = [orig, 
                    (orig[0] + dxyz[0], orig[1],  orig[2]), 
                    (orig[0] + dxyz[0], orig[1] + dxyz[1],  orig[2]), 
                    (orig[0], orig[1] + dxyz[1],  orig[2])]
    else :
        # square looks downwards
        points = [orig,  
                    (orig[0], dxyz[1],  orig[2]), 
                    (dxyz[0], dxyz[1],  orig[2]), 
                    (dxyz[0], orig[1],  orig[2])]  
    normals = [(0, 0, 1) for i in range(4)]
    indices = [(0, 1, 2, 3)]

    # plantGL square
    square_pgl = pgl.QuadSet(points, indices, normals, indices)
    
    # generate the sensors in plantGL format among the grid
    ID_capt = start_id
    dico_translat = {}
    for ix in range(nxyz[0]):
        for iy in range(nxyz[1]):
            for iz in range(nxyz[2] - skylayer):
                # translation vector
                tx = ix * dxyz[0]
                ty = iy * dxyz[1]
                tz = iz * dxyz[2]

                # save the vector
                dico_translat[ID_capt] = [tx, ty, tz]

                # adds a squared sensor to the plantGL scene
                sensor = pgl.Translated(geometry=square_pgl, translation=(tx, ty, tz))
                s_capt.add(pgl.Shape(geometry=sensor, id=ID_capt))
                ID_capt += 1

    # sensors in CARIBU scene format with triangles
    caribu_sensors = {}
    for c in s_capt:
        # Preparation
        pt_lst = c.geometry.geometry.pointList
        idx_lst = c.geometry.geometry.indexList

        # two triangles for each squarred sensor
        for i in range(0, len(idx_lst)):
            triangles = []
            
            # upper triangle
            x = [pt_lst[idx_lst[i][j]][0] + dico_translat[c.id][0] for j in range(3)]
            y = [pt_lst[idx_lst[i][j]][1] + dico_translat[c.id][1] for j in range(3)]
            z = [pt_lst[idx_lst[i][j]][2] + dico_translat[c.id][2] for j in range(3)]
            triangles.append((list(zip(x, y, z))))

            # lower triangle
            ids = [0, 2, 3]
            x = [pt_lst[idx_lst[i][j]][0] + dico_translat[c.id][0] for j in ids]
            y = [pt_lst[idx_lst[i][j]][1] + dico_translat[c.id][1] for j in ids]
            z = [pt_lst[idx_lst[i][j]][2] + dico_translat[c.id][2] for j in ids]
            triangles.append((list(zip(x, y, z))))
            
            # save the triangles
            caribu_sensors[c.id] = triangles
    
    return caribu_sensors, s_capt, sensors_maxcenter

def run_caribu(c_scene, direct_active, infinite, sensors, energy=1.) :
    """runs caribu depending on input options

    :param c_scene: instance of CaribuScene containing geometry, light source(s), opt etc...
    :type c_scene: CaribuScene
    :param direct_active: Whether only first order interception is to be computed, default is True (no rediffusions)
    :type direct_active: bool
    :param infinite: if the user wishes to activate infinitisation of the grid
    :type infinite: bool
    :param sensors: geometric data of the virtual sensors in CARIBU scene format. Sensors are horizontal square made of two triangles
        
        .. code:: python
        
            sensors_caribu = { sensor_id : [triangle1, triangle2], ...}
    
    :type sensors: dict
    :param energy: input energy value
    :type energy: float, optional default is 1.
    :return: results are stored in two dict

        - raw (dict of dict) a {band_name: {result_name: property}} dict of dict.
            Except for result_name='sensors', each property is a {primitive_id: [values,]} dict containing results
            for individual triangles of the primitive
        - aggregated (dict of dict) : a {band_name: {result_name: property}}
            Except for result_name='sensors', each property is a {primitive_id: value} dict containing aggregated
            results for each primitive
        
        result_name are :

                    - area (float): the individual areas (m2)
                    - Eabs (float): the surfacic density of energy absorbed (m-2)
                    - Ei (float): the surfacic density of energy incoming  (m-2) additionally, if split_face is True
                    - Ei_inf (float): the surfacic density of energy incoming on the inferior face (m-2)
                    - Ei_sup (float): the surfacic density of energy incoming on the superior face (m-2)
                    - sensors (dict of dict): area, surfacic density of incoming direct energy and surfacic density of incoming total energy of sensors grouped by id, if any
    
    :rtype: dict of dict, dict of dict
    """    
    if sensors is None :
        raw, aggregated = c_scene.run(direct=direct_active, infinite=infinite)
    else :
        raw, aggregated = c_scene.run(direct=direct_active, infinite=infinite, 
                                                            sensors=sensors)

    # adjust relative results with energy input value
    for band in raw :
        for radtype in ("Eabs", "Ei") :
            for k in raw[band][radtype].keys(): raw[band][radtype][k] = [x * energy for x in raw[band][radtype][k]]
            for k in aggregated[band][radtype].keys(): aggregated[band][radtype][k] *= energy
                
    return raw, aggregated