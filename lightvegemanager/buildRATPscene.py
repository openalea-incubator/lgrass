"""
    buildRATPscene
    **************

    Creation of a PyRATP.grid from inputs. The following functions create and initialize a grid from
    either a triangles mesh, a l-egume grid or an empty geometric input

    The argument ``parameters`` refers to one the three inputs dict of LightVegeManager. It is 
    structured as so:

    .. code:: python

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
import numpy

def extract_grid_origin(parameters, minmax):
    """Create grid origin from either input parameters or from the scene

    :param parameters: RATP parameters from inputs of LightVegeManager
    :type parameters: dict
    :param minmax: list of mininuml point and maximum point in a triangles mesh
    :type minmax: [3-tuple, 3-tuple]
    :return: origin of the grid depending of input parameters and size of the mesh
    :rtype: 3-tuple
    """
    # option 1: origin not specified in the input, it is then the minimum point
    if "origin" not in parameters:
        [xorig, yorig, zorig] = minmax[0]
        zorig = -zorig
    else:
        # option 2: only x and y are specified, z = - zmin
        if len(parameters["origin"]) == 2:
            xorig, yorig, zorig = parameters["origin"][0], parameters["origin"][1], -minmax[0][2]
        # option 3: origin is specified
        elif len(parameters["origin"]) == 3:
            [xorig, yorig, zorig] = parameters["origin"]

    return xorig, yorig, zorig


def build_RATPscene_from_trimesh(
    trimesh,
    minmax,
    triLmax,
    matching_ids,
    parameters,
    coordinates,
    reflected,
    infinite,
    stems_id=None,
    nb_input_scenes=0,
    fullgrid=False
):
    """Build a RATP grid from a triangles mesh.

    :param trimesh: triangles mesh aggregated by indice elements
        .. code-block:: { id : [triangle1, triangle2, ...]}

    :type trimesh: dict
    :param minmax: list of mininuml point and maximum point in a triangles mesh
    :type minmax: [3-tuple, 3-tuple]
    :param triLmax: longest side of all the triangles in trimesh
    :type triLmax: float
    :param matching_ids:
        dict that matches new element indices in trimesh with specy indice and
        input element indice
        .. code:: matching_ids = { new_element_id : (input_element_id, specy_id)}
    :type matching_ids: dict
    :param parameters: RATP parameters from inputs of LightVegeManager
    :type parameters: dict
    :param coordinates: [latitude, longitude, timezone]
    :type coordinates: list
    :param reflected: if the user wishes to activate reflected radiations
    :type reflected: bool
    :param infinite: if the user wishes to activate infinitisation of the grid
    :type infinite: bool
    :param stems_id: list of potential stems element in the input scenes, defaults to None
    :type stems_id: list of 2-tuple, optional
    :param nb_input_scenes: number of input scenes in the geometry dict, defaults to 0
    :type nb_input_scenes: int, optional

    :return:
        It returns 3 objects:

        * ``ratpgrid``
            pyratp.grid filled with areas of triangles in trimesh and input parameters

        * ``matching_tri_vox``
            dict where key is a triangle indice and value the matching voxel indice where the
            barycenter of the triangle is located

        * ``distrib``
            dict with a ``"global"`` key and if ``["angle distrib algo"] = "compute voxel"``
            a second entry where key is ``"voxel"``

            value for ``"voxel"`` is a ``numpy.array`` of dimension
            ``(numberofvoxels,numberofentities,numberofclasses)``

            value for ``"global"`` is a list of ``numberofentities`` list of ``numberofclasses``
            elements

    :rtype: pyratp.grid, dict, dict
    """
    from alinea.pyratp import grid
    from lightvegemanager.voxelsmesh import fill_ratpgrid_from_trimesh

    # computes number of species from matching_ids
    # indices in fortran starts from 1
    numberofentities = max([v[1] for v in matching_ids.values()]) + 1

    ## Building Leaf Angle Distribution (only global distribution) ##
    distrib = {}
    distrib_algo = parameters["angle distrib algo"]

    # computes global distribution
    if distrib_algo == "compute global" or distrib_algo == "compute voxel":
        from lightvegemanager.leafangles import compute_distrib_globale
        distrib["global"] = compute_distrib_globale(trimesh, matching_ids, parameters["nb angle classes"])
    # read a file with leaf angle distribution
    elif distrib_algo == "file":
        from lightvegemanager.leafangles import read_distrib_file
        distrib["global"] = read_distrib_file(parameters["angle distrib file"], numberofentities)

    ## Initialize Grid ##
    # voxels size
    if parameters["voxel size"] == "dynamic":
        dv = [3 * triLmax for i in range(3)]
    else:
        dv = parameters["voxel size"]

    # grid origin
    xorig, yorig, zorig = extract_grid_origin(parameters, minmax)

    # number of voxels
    if "number voxels" in parameters:
        [nx, ny, nz] = parameters["number voxels"]

    else:
        from lightvegemanager.voxelsmesh import compute_grid_size_from_trimesh
        grid_slicing = None
        if "grid slicing" in parameters:
            grid_slicing = parameters["grid slicing"]

        nx, ny, nz = compute_grid_size_from_trimesh(minmax[0], minmax[1], dv, grid_slicing)

    if reflected:
        soil_reflectance = parameters["soil reflectance"]
    else:
        soil_reflectance = [0.0, 0.0]
    ratpgrid = grid.Grid.initialise(
        nx,
        ny,
        nz,
        dv[0],
        dv[1],
        dv[2],
        xorig,
        yorig,
        zorig,
        coordinates[0],
        coordinates[1],
        coordinates[2],
        numberofentities,
        soil_reflectance,
        toric=infinite,
    )

    ## Filling the Grid ##
    # tesselation of triangles on the grid
    if parameters["tesselation level"] > 0:
        from lightvegemanager.voxelsmesh import tesselate_trimesh_on_grid

        trimesh = tesselate_trimesh_on_grid(trimesh, ratpgrid, parameters["tesselation level"])

    # filling
    ratpgrid, matching_tri_vox = fill_ratpgrid_from_trimesh(trimesh, matching_ids, ratpgrid, stems_id, nb_input_scenes)

    # leaf anfle distribution locally by voxels
    if distrib_algo == "compute voxel":
        from lightvegemanager.leafangles import compute_distrib_voxel

        distrib["voxel"] = compute_distrib_voxel(
            trimesh, matching_ids, parameters["nb angle classes"], int(ratpgrid.nveg), matching_tri_vox
        )

    # if you want to activate all voxels in the grid
    if fullgrid:
        dvolume = ratpgrid.dx * ratpgrid.dy * ratpgrid.dz[0]
        k = ratpgrid.nveg
        n_vox_per_nent = []
        for ne in range(ratpgrid.nent):
            k = ratpgrid.nveg
            for ix in range(ratpgrid.njx):
                for iy in range(ratpgrid.njy):
                    for iz in range(ratpgrid.njz):                        
                        if ratpgrid.kxyz[ix, iy, iz] == 0 :
                            S_voxel = 1e-14

                            # ajouter 1 pour utilisation f90
                            ratpgrid.kxyz[ix, iy, iz] = k + 1
                            ratpgrid.numx[k] = ix + 1
                            ratpgrid.numy[k] = iy + 1
                            ratpgrid.numz[k] = iz + 1
                            ratpgrid.nume[ne, k] = ne + 1
                            ratpgrid.nje[k] = max(ne + 1, ratpgrid.nje[k])
                            ratpgrid.nemax = max(ratpgrid.nemax, ratpgrid.nje[k])

                            ratpgrid.leafareadensity[ne, k] += S_voxel / dvolume
                            ratpgrid.s_vt_vx[ne, k] += S_voxel
                            ratpgrid.s_vx[k] += S_voxel
                            ratpgrid.s_vt[ne] += S_voxel
                            ratpgrid.s_canopy += S_voxel

                            k += 1
            n_vox_per_nent.append(k)
        ratpgrid.nveg = max(n_vox_per_nent)


    return ratpgrid, matching_tri_vox, distrib, trimesh


def build_RATPscene_empty(parameters, minmax, coordinates, infinite):
    """Build a RATP grid from an empty geometric input

    :param parameters: RATP parameters from inputs of LightVegeManager
    :type parameters: dict
    :param minmax: list of mininuml point and maximum point in a triangles mesh
    :type minmax: [3-tuple, 3-tuple]
    :param coordinates: [latitude, longitude, timezone]
    :type coordinates: list
    :param infinite: if the user wishes to activate infinitisation of the grid
    :type infinite: bool

    :return:
        * ``ratpgrid``: pyratp.grid sets with input parameters

        * ``distrib``: dict with one entry ``"global" : [1.]``

    :rtype: pyratp.grid, dict
    """
    from alinea.pyratp import grid

    if "number voxels" in parameters:
        [nx, ny, nz] = parameters["number voxels"]
    else:
        nx, ny, nz = 1, 1, 1

    dv = parameters["voxel size"]

    xorig, yorig, zorig = extract_grid_origin(parameters, minmax)

    ratpgrid = grid.Grid.initialise(
        nx,
        ny,
        nz,
        dv[0],
        dv[1],
        dv[2],
        xorig,
        yorig,
        zorig,
        coordinates[0],
        coordinates[1],
        coordinates[2],
        1,
        [0.0, 0.0],
        toric=infinite,
    )
    distrib = {"global": [[1.0]]}

    return ratpgrid, distrib


def legumescene_to_RATPscene(legumescene, parameters, coordinates, reflected, infinite):
    """Creates a RATP grid of voxels from a l-egume grid format

    :param legumescene:

        l-egume grid represented by a dict with two entries:

        *  ``"LA"``: equivalent to m_lais in l-egume, a numpy.array of dimension ``(nent, nz, ny, nx)`` which represents leaf area in each voxel for each specy

        *  ``"distrib"``: equivalent to ls_dif in l-egume, a numpy.array of dimension ``(nent, nclasses)`` which represents the global leaf angle distribution for each input specy

        .. note:: legumescene is the only input geometric scene which can handle several species

    :type legumescene: dict
    :param parameters: RATP parameters from inputs of LightVegeManager
    :type parameters: dict
    :param coordinates: [latitude, longitude, timezone]
    :type coordinates: list
    :param reflected: if the user wishes to activate reflected radiations
    :type reflected: bool
    :param infinite: if the user wishes to activate infinitisation of the grid
    :type infinite: bool

    :return:

        It returns 3 objects

        * ``ratpgrid`` : pyratp.grid filled with areas of triangles in trimesh and input parameters

        * ``distrib`` : dict with only a ``"global"`` key, value for ``"global"`` is a list of ``numberofentities`` list of ``numberofclasses`` elements

        * ``nb0`` : number of empty layers between top canopy and maximum layer in legumescene

    :rtype: pyratp.grid, dict, int
    """
    from alinea.pyratp import grid
    from lightvegemanager.voxelsmesh import fill_ratpgrid_from_legumescene

    distrib = {"global": legumescene["distrib"]}

    # extract grid parameters
    numberofentities = legumescene["LA"].shape[0]
    dv = parameters["voxel size"]
    xorig, yorig, zorig = 0, 0, 0
    nx = legumescene["LA"].shape[3]
    ny = legumescene["LA"].shape[2]
    nz = legumescene["LA"].shape[1]

    # count all the empty layers between top canopy and maximum layers in legumescene
    nb0 = 0
    laicum = numpy.sum(legumescene["LA"], axis=0)
    laicumvert = numpy.sum(laicum, axis=(1, 2))
    for i in range(len(laicumvert)):
        if laicumvert[i] == 0.0:
            nb0 += 1
        else:
            break
    nz = nz - nb0

    # grid initialization
    if reflected:
        soil_reflectance = parameters["soil reflectance"]
    else:
        soil_reflectance = [0.0, 0.0]
    ratpgrid = grid.Grid.initialise(
        nx,
        ny,
        nz,
        dv[0],
        dv[1],
        dv[2],
        xorig,
        yorig,
        zorig,
        coordinates[0],
        coordinates[1],
        coordinates[2],
        numberofentities,
        soil_reflectance,
        toric=infinite,
    )

    # grid filling
    fill_ratpgrid_from_legumescene(legumescene, ratpgrid, nb0)

    return ratpgrid, distrib, nb0

def concatene_legumescenes(scenes):
    return_scene = scenes[0]
    for grid in scenes[1:]:
        return_scene["LA"] = numpy.append(return_scene["LA"], grid["LA"], axis=0)
        return_scene["distrib"].append(grid["distrib"][0])
    return return_scene