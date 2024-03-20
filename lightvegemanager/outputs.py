"""
   
    outputs
    *******

    Manages and reformat output results from the light models in pandas.Dataframe with similar columns
    names. 
    
    Light models managed in this module: 

        - RATP
        - CARIBU

    Column names for voxels:

        * ``'VegetationType'``: id of the specy
        * ``'Iteration'``: RATP can handle multiple iterations in inputs, but it is not properly exploited in LightVegeManager yet
        * ``'Day'``: input day if LightVegeManager
        * ``'Hour'``: input hour if LightVegeManager
        * ``'Voxel'``: index of current voxel
        * ``'Nx'``: index in grid following x axis, corresponds to pyrapt.grid.numx
        * ``'Ny'``: index in grid following y axis, corresponds to pyrapt.grid.numy
        * ``'Nz'``: index in grid following z axis, corresponds to pyrapt.grid.numz
        * ``'ShadedPAR'``: shaded energy in the current voxel, 0 if no direct light in the input
        * ``'SunlitPAR'``: sunlit energy
        * ``'ShadedArea'``: leaf area in voxel which is shaded
        * ``'SunlitArea'``: leaf area in voxel which is sunlit
        * ``'Area'``: total leaf area in the voxel
        * ``'PARa'``: ( ShadedPAR * ShadedArea + SunlitPAR * SunlitArea ) / ( ShadedArea + SunlitArea )
        * ``'Intercepted'``: relative portion of input rays intercepted by the voxel
        * ``'Transmitted'``: relative portion of input rays leaving out the voxel

    Column names for triangles:

        * ``'VegetationType'``: id of the specy
        * ``'Day'``: input day if LightVegeManager
        * ``'Hour'``: input hour if LightVegeManager
        * ``'Triangle'``: indice of the triangle
        * ``'Organ'``: element where the triangle belongs

        
        if ligh model is CARIBU:

            * ``'Area'``: area of the triangles
        
            for each bandwidth in the inputs you have two entries:

            * ``band + " Eabs"``: energy absorbed by the triangle
            * ``band + " Ei"``: energy intercepted by the triangle
        
        If light model is RATP:

            * ``'Area'``: area of total leaf area in the voxel where the triangle is located
            * ``'primitive_area'``: area of the triangle
            * ``'Voxel'``: index of the voxel where the triangle is located
            * ``'Nx'``: index in grid following x axis, corresponds to pyrapt.grid.numx
            * ``'Ny'``: index in grid following y axis, corresponds to pyrapt.grid.numy
            * ``'Nz'``: index in grid following z axis, corresponds to pyrapt.grid.numz
            * ``'ShadedPAR'``: shaded energy in the current voxel
            * ``'SunlitPAR'``: sunlit energy, 0 if no direct light in the input
            * ``'ShadedArea'``: leaf area in voxel which is shaded
            * ``'SunlitArea'``: leaf area in voxel which is sunlit
            * ``'PARa'``: ( ShadedPAR * ShadedArea + SunlitPAR * SunlitArea ) / ( ShadedArea + SunlitArea )
            * ``'Intercepted'``: relative portion of input rays intercepted by the voxel
            * ``'Transmitted'``: relative portion of input rays leaving out the voxel

    Column names for elements:

        * ``'Day'``: input day if LightVegeManager
        * ``'Hour'``: input hour if LightVegeManager
        * ``'Organ'``: id of current element
        * ``"VegetationType"``: id of the specy
        * ``"Area"``: area of the element (sum of all triangles belonging to current element)

        if light model is CARIBU, for each bandwidth in the optical inputs you have two entries:

            * ``band + " Eabs"``: energy absorbed by the element, integrated on all triangles belonging to current element
            * ``band + " Ei"``: energy intercepted by the triangle, integrated on all triangles belonging to current element


        if light model is RATP, the following entries are integrated on all triangles belonging to current element 
        with E, current element, t triangles in E and V an entry: :math:`\\frac{\\sum_{t \\in E} area_t * V_t}{area_E}`

            * ``"PARa"``
            * ``"Intercepted"``
            * ``"Transmitted"``
            * ``"SunlitPAR"``
            * ``"SunlitArea"``
            * ``"ShadedPAR"``
            * ``"ShadedArea"``
"""

import pandas
import itertools

from lightvegemanager.basicgeometry import triangle_area

# from lightvegemanager.trianglesmesh import *

def out_ratp_empty_grid(day, hour) :
    """Returns an empty dataframe results for RATP

    :param day: day of simulation
    :type day: int
    :param hour: hour of simulation
    :type hour: int
    :return: returns a dataframe with all the keys expected with a RATP simulation but results are empty
        used if input geometry is empty 
    :rtype: pandas.DataFrame
    """    
    df_voxels =  pandas.DataFrame({'VegetationType':[0],
                                'Iteration':[1],
                                'Day':day,
                                'Hour':hour,
                                'Voxel':[0],
                                'Nx':[0],
                                'Ny':[0],
                                'Nz':[0],
                                'ShadedPAR':[0],
                                'SunlitPAR':[0],
                                'ShadedArea':[0],
                                'SunlitArea': [0],
                                'Area': [0],
                                'PARa': [0],
                                'Intercepted': [0],
                                'Transmitted': [0]
                            })
    
    df_triangles =  pandas.DataFrame({'VegetationType':[0],
                            'Iteration':[1],
                            'Day':day,
                            'Hour':hour,
                            'Voxel':[0],
                            'Nx':[0],
                            'Ny':[0],
                            'Nz':[0],
                            'Triangle': [0],
                            'Organ': [0], 
                            'primitive_area': [0],
                            'ShadedPAR':[0],
                            'SunlitPAR':[0],
                            'ShadedArea':[0],
                            'SunlitArea': [0],
                            'Area': [0],
                            'PARa': [0],
                            'Intercepted': [0],
                            'Transmitted': [0]
                        })
    df_organs =  pandas.DataFrame({'VegetationType':[0],
                        'Day':day,
                        'Hour':hour,
                        'Organ': [0], 
                        'ShadedPAR':[0],
                        'SunlitPAR':[0],
                        'ShadedArea':[0],
                        'SunlitArea': [0],
                        'Area': [0],
                        'PARa': [0],
                        'Intercepted': [0],
                        'Transmitted': [0]
                    })
    
    return df_voxels, df_triangles, df_organs


def out_ratp_voxels(ratpgrid, res, parunit) :
    """Converts RATP results to pandas dataframe compared to the voxels

    :param ratpgrid: RATP grid with voxels informations
    :type ratpgrid: pyratp.grid
    :param res: output table of RATP
    :type res: pyratp.ratp.out_rayt
    :param parunit: energy unit of input

        RATP expects W.m-2 in input and return results in µmol.s-1.m-2.
        if ``parunit="W.m-2"`` the outputs is converted to W.m-2, otherwise results are in µmol.s-1.m-2

    :type parunit: string
    :return: res.T converted in a pandas Dataframe with voxels relative columns. The soil is not part of the results
    :rtype: pandas.DataFrame
    """    
    # decompress all outputs of RATP
    # each output is a numpy.array of one dimension of size number of voxels x number of iterations
    VegetationType, \
    Iteration,      \
    day,            \
    hour,           \
    VoxelId,        \
    ShadedPAR,      \
    SunlitPAR,      \
    ShadedArea,     \
    SunlitArea,     \
    xintav,         \
    Ptransmitted = res.T

    if parunit == "W.m-2" :
        # ('PAR' is expected in  Watt.m-2 in RATP input, whereas output is in 
        # micromol => convert back to W.m2 (cf shortwavebalance, line 306))
        ShadedPAR = ShadedPAR / 4.6
        SunlitPAR = SunlitPAR / 4.6
    
    para_list=[]
    for i in range(len(ShadedPAR)):
        if (ShadedArea[i] + SunlitArea[i]) > 0 :
            shaded = ShadedPAR[i] * ShadedArea[i]
            sunlit = SunlitPAR[i] * SunlitArea[i]
            area = ShadedArea[i] + SunlitArea[i]
            para_list.append((shaded + sunlit )/area)
        else:
            para_list.append(0.)
    
    # check if we don't have negative erel
    erel_list=[]
    for i in range(len(xintav)):
        if xintav[i] >= 1e-6 :
            erel_list.append(xintav[i])
        else:
            erel_list.append(0.)

    # voxel indices list read in RATP grid
    numx=[]
    numy=[]
    numz=[]
    for v in VoxelId :
        numx.append(ratpgrid.numx[int(v)-1])
        numy.append(ratpgrid.numy[int(v)-1])
        numz.append(ratpgrid.numz[int(v)-1])

    dfvox =  pandas.DataFrame({'VegetationType':VegetationType,
                                    'Day':day,
                                    'Hour':hour,
                                    'Voxel':VoxelId,
                                    'Nx':numx,
                                    'Ny':numy,
                                    'Nz':numz,
                                    'ShadedPAR':ShadedPAR,
                                    'SunlitPAR':SunlitPAR,
                                    'ShadedArea':ShadedArea,
                                    'SunlitArea': SunlitArea,
                                    'Area': ShadedArea + SunlitArea,
                                    'PARa': para_list,
                                    'Intercepted': erel_list, 
                                    'Transmitted': Ptransmitted
                                })
                
    # we don't return the soil
    return dfvox[dfvox['VegetationType'] > 0]

def out_ratp_triangles(trimesh,
                        matching_ele_ent,
                        matching_tr_vox,
                        voxels_outputs) :
    """Converts RATP results to pandas dataframe compared to the triangles (and voxels)

    :param trimesh: triangles mesh aggregated by indice elements
        .. code-block:: { id : [triangle1, triangle2, ...]}
    :type trimesh: dict
    :param matching_ele_ent: 
        dict that matches new element indices in trimesh with specy indice and
        input element indice, 
        .. code:: matching_ids = { new_element_id : (input_element_id, specy_id)}

    :type matching_ele_ent: dict
    :param matching_tr_vox: dict where key is a triangle indice and value the matching voxel indice where the 
        barycenter of the triangle is located   
    :type matching_tr_vox: dict
    :param voxels_outputs: output of :func:out_ratp_voxels
    :type voxels_outputs: pandas.Dataframe
    :return: output of :func:out_ratp_voxels merged with associated triangles
    :rtype: pandas.Dataframe
    """   
    from lightvegemanager.trianglesmesh import globalid_to_elementid

    # creation of triangle table
    entity = {}
    for id, match in matching_ele_ent.items():
        entity[id] = match[1] + 1
    
    index = range(len(matching_tr_vox))
    vox_id = [matching_tr_vox[str(i)] + 1 for i in index]
    triangles = list(itertools.chain(*trimesh.values()))
    sh_id = [globalid_to_elementid(trimesh, i) for i in range(len(triangles))]
    s = [triangle_area(t) for t in triangles]

    # new dataframe with triangles indexed
    dftriangles = pandas.DataFrame({'Triangle': index,
                                'Organ': sh_id, 
                                'Voxel':vox_id, 
                                'VegetationType':[entity[id] for id in sh_id], 
                                'primitive_area':s})

    # common columns for merging : VegetationType, VoxelID
    trianglesoutputs = pandas.merge(dftriangles, voxels_outputs)  
    
    # sort lines by triangle indices
    return  trianglesoutputs.sort_values('Triangle')

def out_ratp_elements(matching_ele_ent, 
                        reflected, 
                        reflectance_coef, 
                        trianglesoutputs) :
    """If trimesh is not empty, aggregates RATP results by element (keys in dict trimesh)

    :param matching_ele_ent: 
        dict that matches new element indices in trimesh with specy indice and
        input element indice, 
        .. code:: matching_ids = { new_element_id : (input_element_id, specy_id)}

    :type matching_ele_ent: dict
    :param reflected: if the user wishes to activate reflected radiations
    :type reflected: bool
    :param reflectance_coef: coefficient for each specy and each input bandwidth 
    :type reflectance_coef: list of list
    :param trianglesoutputs: output of :func:out_ratp_triangles
    :type trianglesoutputs: pandas.Dataframe
    :return: dataframe integrated on element informations
    :rtype: pandas.Dataframe
    """                        
    # save values by element and specy
    nshapes = len(matching_ele_ent)
    s_shapes = []
    s_area=[]
    s_para=[]
    s_pari=[]
    s_intercepted=[]
    s_transmis=[]
    s_day=[]
    s_hour=[]
    s_ent=[]
    s_parsun=[]
    s_parsha=[]
    s_areasun=[]
    s_areasha=[]
    for id in range(nshapes):
        # iterations starts at 1 (fortran)
        nent = matching_ele_ent[id][1]
        dffil = trianglesoutputs[(trianglesoutputs.Organ == id)]
        
        t_areas = dffil["primitive_area"]
        sum_area = sum(t_areas)
        
        s_hour.append(dffil["Hour"].values[0])
        s_day.append(dffil["Day"].values[0])
        s_area.append(sum_area)

        # if reflected rays has been computed in RATP
        if reflected :
            para = sum(t_areas * dffil['PARa']) / sum_area
            s_para.append(para)
        
        else:
            # incident PAR
            pari = sum(t_areas * dffil['PARa']) / sum_area
            s_pari.append(pari)
            s_para.append(pari - (pari * reflectance_coef[nent][0]))
        
        s_parsun.append(sum(t_areas * dffil['SunlitPAR']) / sum_area)
        s_parsha.append(sum(t_areas * dffil['ShadedPAR']) / sum_area)
        s_areasun.append(sum(t_areas * dffil['SunlitArea']) / sum_area)
        s_areasha.append(sum(t_areas * dffil['ShadedArea']) / sum_area)
        s_intercepted.append(sum(t_areas * dffil['Intercepted']) /sum_area)
        s_transmis.append(sum(t_areas * dffil['Transmitted']) /sum_area)
        s_ent.append(dffil["VegetationType"].values[0])
        s_shapes.append(matching_ele_ent[id][0])
    
    return pandas.DataFrame({
                                "Day" : s_day,
                                "Hour" : s_hour,
                                "Organ" : s_shapes,
                                "VegetationType" : s_ent,
                                "Area" : s_area,
                                "PARa" : s_para,
                                "Intercepted" : s_intercepted,
                                "Transmitted" : s_intercepted,
                                "SunlitPAR" : s_parsun,
                                "SunlitArea" : s_areasun,
                                "ShadedPAR" : s_parsha,
                                "ShadedArea" : s_areasha
                            })

def out_caribu_mix( rdrs,
                    c_scene_sun, c_scene_sky,
                    raw_sun, aggregated_sun, 
                    raw_sky, aggregated_sky,
                    issensors,
                    issoilmesh) :
    """Mix diffuse and direct rays results according to a ratio rdrs

    :param rdrs: Spitters's model estimating for the diffuse:direct ratio
    :type rdrs: float
    :param c_scene_sun: CaribuScene with direct lighting computed (from a sun)
    :type c_scene_sun: CaribuScene
    :param c_scene_sky: CaribuScene with diffuse lighting computed (from a sky)
    :type c_scene_sky: CaribuScene
    :param raw_sun: results of c_scene_sun for triangles
    :type raw_sun: dict of dict
    :param aggregated_sun: results of c_scene_sun aggregated by element
    :type aggregated_sun: dict of dict
    :param raw_sky: results of c_scene_sky for triangles
    :type raw_sky: dict of dict
    :param aggregated_sky:  results of c_scene_sky aggregated by element
    :type aggregated_sky: dict of dict
    :param issensors: if option virtual sensors is activated
    :type issensors: bool
    :param issoilmesh: if option soil mesh is activated
    :type issoilmesh: bool
    :return:
        if sensors is activated, it extracts its results from aggregated_sun and aggregated_sky

        if soilmesh is activated, it extracts its results from c_scene_sun and c_scene_sky
        
        then it mixes raw, aggregated, virtual sensors and soil energy results such as 
        
        ``result = rdrs * diffuse_result + (1 - rdrs) * direct_result``

    :rtype: dict of dict, dict of dict, dict of list, dict of float
    """
    raw, aggregated, sensors, soil_energy = raw_sun, aggregated_sun, {}, {}
    for band in raw.keys() : 
        for ray in ['Eabs', 'Ei'] :
            for ide, v_sun in raw[band][ray].items() :
                v_sky = raw_sky[band][ray][ide]
                raw[band][ray][ide] = [rdrs * sky + (1 - rdrs) * sun for sun,sky in zip(v_sun, v_sky)]

            for ide, v_sun in aggregated[band][ray].items() :
                v_sky = aggregated_sky[band][ray][ide]
                aggregated[band][ray][ide] = rdrs * v_sky + (1 - rdrs) * v_sun

        if issensors :
            sensors[band] = {}
            for id, v_sun in aggregated_sun[band]['sensors']['Ei'].items() :
                v_sky = aggregated_sky[band]['sensors']['Ei'][id]
                sensors[band][id] = rdrs * v_sky + (1 - rdrs) * v_sun

    if issoilmesh : 
        q1, e1 = c_scene_sky.getSoilEnergy()
        q2, e2 = c_scene_sun.getSoilEnergy()
        q = rdrs * q1 + (1 - rdrs) * q2
        e = rdrs * e1 + (1 - rdrs) * e2

        soil_energy['Qi'] = q
        soil_energy['Einc'] = e

    return raw, aggregated, sensors, soil_energy

def out_caribu_nomix(c_scene, aggregated, issensors, issoilmesh) :
    """Extracts only sensors and soilmesh results if activated

    :param c_scene: instance of CaribuScene containing geometry, light source(s), opt etc...
    :type c_scene: CaribuScene
    :param aggregated: result of CaribuScene.run 
    :type aggregated: dict of dict
    :param issensors: if option virtual sensors is activated
    :type issensors: bool
    :param issoilmesh: if option soil mesh is activated
    :type issoilmesh: bool
    :return: extracts sensors results from aggregated and soil energy from c_scene
    :rtype: dict, dict
    """    
    sensors, soil_energy = {}, {}
    
    if issensors :
        for band in aggregated.keys() :
            sensors[band] =  aggregated[band]['sensors']['Ei']

    if issoilmesh :
        q, e = c_scene.getSoilEnergy()
        
        soil_energy['Qi'] = q
        soil_energy['Einc'] = e

    return sensors, soil_energy

def out_caribu_elements(day, hour, trimesh, matching_ids, aggregated, sun_up, caribu_triangles={}) :
    """Converts aggregated in a pandas.Dataframe following indices in LightVegeManager

    :param day: day of simulation
    :type day: int
    :param hour: hour of simulation
    :type hour: int
    :param trimesh: triangles mesh aggregated by indice elements
        .. code-block:: { id : [triangle1, triangle2, ...]}
    :type trimesh: dict
    :param matching_ele_ent: 
        dict that matches new element indices in trimesh with specy indice and
        input element indice, 
        .. code:: matching_ids = { new_element_id : (input_element_id, specy_id)}

    :type matching_ele_ent: dict
    :param aggregated: result of CaribuScene.run 
    :type aggregated: dict of dict
    :param sun_up: if sun elevation is > 2° and direct rays are > 0 W.m²
    :type sun_up: bool
    :return: aggregated results rearrange in a Dataframe with element correspondance in LightVegeManager
    :rtype: pandas.Dataframe
    """
    (s_shapes, 
    s_area, 
    s_day, 
    s_hour, 
    s_ent) = ([0] * len(matching_ids) for i in range(5))
    
    for key,val in matching_ids.items():
        s_shapes[key] = val[0]
        s_area[key] = sum([triangle_area(t) for t in trimesh[key]])
        s_day[key] = day
        s_hour[key] = hour
        s_ent[key] = matching_ids[key][1]

    dico_shape = {
        "Day" : s_day,
        "Hour" : s_hour,
        "Organ" : s_shapes,
        "VegetationType" : s_ent,
        "Area" : s_area
    }
    
    s_Eabs = {}
    s_Ei = {}
    for band, dict_val in aggregated.items() : 
        s_Eabs[band] = [0]*len(matching_ids)
        s_Ei[band] = [0]*len(matching_ids)
        for key,val in matching_ids.items():
            if sun_up: s_Eabs[band][key] = dict_val['Eabs'][key]
            if sun_up: s_Ei[band][key] = dict_val['Ei'][key]

        if sum(s_Eabs[band]) <= 0. :
            v,w = [], []
            for i in range(max(caribu_triangles["Organ"]) + 1):
                data = caribu_triangles[caribu_triangles.Organ == i]
                v.append(sum([x*y for x,y in zip(data["Area"], data[band + " Eabs"])]) / dico_shape["Area"][i])
                w.append(sum([x*y for x,y in zip(data["Area"], data[band + " Ei"])]) / dico_shape["Area"][i])
            
            dico_shape[band + " Eabs"] = v
            dico_shape[band + " Ei"] = w
        else :
            dico_shape[band + " Eabs"] = s_Eabs[band]
            dico_shape[band + " Ei"] = s_Ei[band]

    return pandas.DataFrame(dico_shape)

def out_caribu_triangles(day, hour, trimesh, matching_ids, raw, sun_up) :
    """Converts raw in a pandas.Dataframe following simulation datas

    :param day: day of simulation
    :type day: int
    :param hour: hour of simulation
    :type hour: int
    :param trimesh: triangles mesh aggregated by indice elements
        .. code-block:: { id : [triangle1, triangle2, ...]}
    :type trimesh: dict
    :param matching_ele_ent: 
        dict that matches new element indices in trimesh with specy indice and
        input element indice, 
        .. code:: matching_ids = { new_element_id : (input_element_id, specy_id)}

    :type matching_ele_ent: dict
    :param raw: triangles results from CaribuScene.run
    :type raw: dict of dict
    :param sun_up: if sun elevation is > 2° and direct rays are > 0 W.m²
    :type sun_up: bool
    :return: triangle results rearrange in a Dataframe
    :rtype: pandas.Dataframe
    """    
    nb_triangles = sum([len(triangles) for triangles in trimesh.values()])
    
    (s_shapes, 
    s_tr,
    s_area, 
    s_day, 
    s_hour, 
    s_ent) = ([0] * nb_triangles for i in range(6))

    id_tr = 0
    for ele, triangles in trimesh.items():
        for t in triangles :
            s_shapes[id_tr] = ele
            s_tr[id_tr] = id_tr
            s_area[id_tr] = triangle_area(t)
            s_day[id_tr] = day
            s_hour[id_tr] = hour
            s_ent[id_tr] = matching_ids[ele][1]
            id_tr += 1

    dico_tr = {
        "Day" : s_day,
        "Hour" : s_hour,
        "Triangle" : s_tr,
        "Organ" : s_shapes,
        "VegetationType" : s_ent,
        "Area" : s_area
    }
    
    s_Eabs = {}
    s_Ei = {}
    for band, dict_val in raw.items() : 
        if sun_up: s_Eabs[band] = \
            list(itertools.chain(*[dict_val['Eabs'][ele] for ele in trimesh.keys()]))
        if sun_up: s_Ei[band] =  \
            list(itertools.chain(*[dict_val['Ei'][ele] for ele in trimesh.keys()]))
    
        dico_tr[band + " Eabs"] = s_Eabs[band]
        dico_tr[band + " Ei"] = s_Ei[band]

    return pandas.DataFrame(dico_tr)

def out_caribu_sensors(day, hour, sensors, matching_sensors_species) :
    """Converts raw in a pandas.Dataframe following simulation datas

    :param day: day of simulation
    :type day: int
    :param hour: hour of simulation
    :type hour: int
    :param sensors: lighting results for each virtual sensors in caribu
    :type sensors: dict
    :param matching_sensors_species: dict that match each sensor with a plant specy
    :type matching_sensors_species: dict
    :return: sensors results rearrange in a Dataframe
    :rtype: pandas.Dataframe
    """       
    s_par = [value for value in sensors["par"].values()]
    s_sen = [key for key in sensors["par"].keys()]
    s_day = [day] * len(sensors["par"])
    s_hour = [hour] * len(sensors["par"])
    s_ent = [value for value in matching_sensors_species.values()]

    dico_sensors = {
        "Day" : s_day,
        "Hour" : s_hour,
        "Sensor" : s_sen,
        "PAR" : s_par,
        "VegetationType" : s_ent,
    }
    
    return pandas.DataFrame(dico_sensors)