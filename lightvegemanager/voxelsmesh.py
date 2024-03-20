"""
    voxelsmesh
    **********

    Builds and handles axis oriented voxels mesh
"""

import itertools

from lightvegemanager.basicgeometry import triangle_barycenter, triangle_area

def compute_grid_size_from_trimesh(pmin, pmax, dv, grid_slicing=None):
    """Dynamically compute number of voxels for each axis in the grid
    The grid is ajusted to be the smallest box containing another mesh

    :param pmin: Minimum point of a mesh ``[x, y, z]``
    :type pmin: list
    :param pmax: Maximum point of a mesh ``[x, y, z]``
    :type pmax: list
    :param dv: size a voxel in each direction ``[dx, dy, dz]``
    :type dv: list
    :param grid_slicing: possibility to force the ground to be at z=0, defaults to None
    :type grid_slicing: string, optional
    :return: number of voxels in each direction x, y and z
    :rtype: int, int, int
    """
    [dx, dy, dz] = dv

    # if pmax == pmin we ajust to have one layer of voxels
    for i in range(3):
        if pmin[i] == pmax[i]:
            pmax[i] += dv[i]
            pmin[i] -= dv[i]

    nx = int((pmax[0] - pmin[0]) // dx)
    ny = int((pmax[1] - pmin[1]) // dy)
    nz = 1
    if grid_slicing is not None:
        if grid_slicing == "ground = 0.":
            nz = int((pmax[2] - 0.0) // dz)
    else:
        nz = int((pmax[2] - pmin[2]) // dz)
    if (pmax[0] - pmin[0]) % dx > 0:
        nx += 1
    if (pmax[1] - pmin[1]) % dy > 0:
        ny += 1
    if (pmax[2] - pmin[2]) % dz > 0:
        nz += 1

    return nx, ny, nz


def tesselate_trimesh_on_grid(trimesh, ratpgrid, levelmax):
    """Loop on all triangles of a triangulation to tesselate them for better matching a grid of voxels
    Triangles will subdivide on sides of voxels respecting a certain maximum level

    :param trimesh: triangles mesh aggregated by indice elements

        .. code-block::

            { id : [triangle1, triangle2, ...]}

    :type trimesh: dict of list
    :param ratpgrid: RATP grid of voxels
    :type ratpgrid: pyratp.grid
    :param levelmax: maximum level for subdividing triangles
    :type levelmax: int
    :return: a copy of trimesh with subdivided triangles
    :rtype: dict of list
    """
    from lightvegemanager.tesselator import iterate_trianglesingrid

    new_trimesh = {}
    for id_ele, triangles in trimesh.items():
        new_tr_scene = []
        for t in triangles:
            level = 0
            iterate_trianglesingrid(t, ratpgrid, level, levelmax, new_tr_scene)
        new_trimesh[id_ele] = new_tr_scene

    return new_trimesh


def fill_ratpgrid_from_trimesh(trimesh, matching_ids, ratpgrid, stems_id=None, nb_input_scenes=0):
    """Fills a RATP grid from a triangulation.
    It gives barycenters and areas of triangles to update the leaf area density in the corresponding
    voxel.

    :param trimesh: triangles mesh aggregated by indice elements

        .. code-block::

            { id : [triangle1, triangle2, ...]}

    :type trimesh: dict of list
    :param matching_ids:
        dict that matches new element indices in cscene with specy indice and
        input element indice,

        .. code-block::

            matching_ids = { new_element_id : (input_element_id, specy_id)}

        this dict allows us to look how species there is the inputs geometric data

    :type matching_ids: dict
    :param ratpgrid: RATP grid of voxels
    :type ratpgrid: pyratp.grid
    :param stems_id: list of potential stems element in the input scenes, defaults to None
    :type stems_id: list of 2-tuple, optional
    :param nb_input_scenes: number of input geometrical scenes. It can be different with number of species if there is a l-egume grid in the input with several species in it, defaults to 0
    :type nb_input_scenes: int, optional
    :return: copy of ``ratpgrid`` with leaf area density values from barycenters and areas of the input triangles
    :rtype: pyratp.grid
    """
    from alinea.pyratp import grid

    entity, barx, bary, barz, a, n = [], [], [], [], [], []
    for id, ts in trimesh.items():
        for t in ts:
            bar = triangle_barycenter(t)
            barx.append(bar[0])
            bary.append(bar[1])
            barz.append(bar[2])

            # if element is a stem, i.e. an opaque organ, we divide its area by 2 (only one face receive lighting)
            if (isinstance(stems_id, list) or isinstance(stems_id, tuple)) and \
            ((matching_ids[id][0], matching_ids[id][1] - nb_input_scenes) in stems_id or \
            (matching_ids[id][0],matching_ids[id][1]) in stems_id):
                a.append(triangle_area(t) / 2)
            else:
                a.append(triangle_area(t))

            n.append(0.0)
            entity.append(matching_ids[id][1])

    return grid.Grid.fill_1(entity, barx, bary, barz, a, n, ratpgrid)


def fill_ratpgrid_from_legumescene(legumescene, ratpgrid, nb0):
    """Fills a RATP grid from a l-egume grid
    It updates ``ratpgrid``.

    :param legumescene:
        l-egume grid represented by a dict with two entries:

        *  ``"LA"``: equivalent to m_lais in l-egume, a numpy.array of dimension ``(nent, nz, ny, nx)`` which represents leaf area in each voxel for each specy. nz=0 is the top layer

        *  ``"distrib"``: equivalent to ls_dif in l-egume, a numpy.array of dimension ``(nent, nclasses)`` which represents the global leaf angle distribution for each input specy

        .. note:: legumescene is the only input geometric scene which can handle several species
    :type legumescene: dict
    :param ratpgrid: RATP grid of voxels. nz=0 is the top layer
    :type ratpgrid: pyratp.grid
    :param nb0: number of empty layers from top of the canopy and maximum z layers in m_lais
    :type nb0: int
    """
    dvolume = ratpgrid.dx * ratpgrid.dy * ratpgrid.dz[0]
    # remplissage des voxels avec la grille
    n_vox_per_nent = []
    for ne in range(ratpgrid.nent):
        k = 0
        for ix in range(ratpgrid.njx):
            for iy in range(ratpgrid.njy):
                for iz in range(ratpgrid.njz):
                    legume_iz = iz + nb0

                    # on force les voxels vides dans les couches non vides à
                    # être interprétés par RATP
                    S_voxel = max(1e-14, legumescene["LA"][ne][legume_iz][iy][ix])

                    # ajouter 1 pour utilisation f90
                    # sens xy sont inversés entre l-egume et RATP
                    ratpgrid.kxyz[ratpgrid.njy - (iy + 1), ix, iz] = k + 1
                    ratpgrid.numx[k] = (ratpgrid.njy - (iy + 1)) + 1
                    ratpgrid.numy[k] = ix + 1
                    ratpgrid.numz[k] = iz + 1
                    ratpgrid.nume[ne, k] = ne + 1
                    ratpgrid.nje[k] = max(ne + 1, ratpgrid.nje[k])
                    ratpgrid.nemax = max(ratpgrid.nemax, ratpgrid.nje[k])

                    ratpgrid.leafareadensity[ne, k] += S_voxel / dvolume
                    ratpgrid.s_vt_vx[ne, k] += S_voxel
                    ratpgrid.s_vx[k] += S_voxel
                    ratpgrid.s_vt[ne] += S_voxel
                    ratpgrid.s_canopy += S_voxel

                    k = k + 1
        n_vox_per_nent.append(k)

    ratpgrid.nveg = max(n_vox_per_nent)
    ratpgrid.nsol = ratpgrid.njx * ratpgrid.njy  # Numbering soil surface areas
    for jx in range(ratpgrid.njx):
        for jy in range(ratpgrid.njy):
            ratpgrid.kxyz[jx, jy, ratpgrid.njz] = ratpgrid.njy * jx + jy + 1

    for k in range(ratpgrid.nveg):
        for je in range(ratpgrid.nje[k]):
            if je == 0:
                # !!! on considère l-egume avec une hauteur fixe de voxel
                # Incrementing total canopy volume
                ratpgrid.volume_canopy[ratpgrid.nent] = ratpgrid.volume_canopy[ratpgrid.nent] + dvolume
            if ratpgrid.s_vt_vx[je, k] > 0.0:
                ratpgrid.volume_canopy[ratpgrid.nume[je, k] - 1] = (
                    ratpgrid.volume_canopy[ratpgrid.nume[je, k] - 1] + dvolume
                )
                ratpgrid.voxel_canopy[ratpgrid.nume[je, k] - 1] = ratpgrid.voxel_canopy[ratpgrid.nume[je, k] - 1] + 1


def reduce_layers_from_trimesh(trimesh, pmax, dxyz, nxyz, matching_ids, ids=None):
    """Number of empty layers in a grid of voxels, from the top of the canopy to the last expected layer.
    It computes this number from a triangles mesh

    :param trimesh: triangles mesh aggregated by indice elements

        .. code-block::

            { id : [triangle1, triangle2, ...]}

    :type trimesh: dict of list
    :param pmax: Maximum point of a mesh ``[x, y, z]``
    :type pmax: list
    :param dxyz: size of sides of a voxel ``[dx, dy, dz]``
    :type dxyz: list
    :param nxyz: number of voxels in each direction ``[nx, ny, nz]``, nz is the expected number of layers
    :type nxyz: list
    :param matching_ids:
        dict that matches new element indices in cscene with specy indice and
        input element indice,

        .. code-block::

            matching_ids = { new_element_id : (input_element_id, specy_id)}

        this dict allows us to look how species there is the inputs geometric data

    :type matching_ids: dict
    :param ids: list of specy indices considered not empty, defaults to None
    :type ids: list of int, optional
    :return: number of empty layers between top of the canopy (represented by ``trimesh``) and ``nxyz[2]``
    :rtype: int
    """
    from lightvegemanager.trianglesmesh import triangles_entity
    
    # empty scene
    if not matching_ids :
        skylayer = int(nxyz[2] - 2)

    else:
        # reduction si number of filled layers < expected number of layers
        if ids is None:
            skylayer = (pmax[2]) // dxyz[2]
            if skylayer < nxyz[2] and pmax[2] > 0:
                skylayer = int(nxyz[2] - 1 - skylayer)

            # otherwise we keep the initial number of layers
            else:
                skylayer = 0

        # we compute the maximum z of trimesh by omitting some species listed in ids
        else:
            # we firstly check if some input species are not empty
            i = 0
            specie_empty = matching_ids[i][1] not in ids
            while (matching_ids[i][1] not in ids) and (i + 1 < len(matching_ids)):
                i += 1
                if matching_ids[i][1] in ids:
                    specie_empty = False

            # geometry is empty
            if specie_empty:
                skylayer = nxyz[2] - 1

            # we compute empty layers from the non empty species
            else:
                zmax = -999999
                if trimesh:
                    for i in ids:
                        triangles = triangles_entity(trimesh, i, matching_ids)
                        for z in list(itertools.chain(*[[p[2] for p in t] for t in triangles])):
                            if z > zmax:
                                zmax = z

                skylayer = zmax // dxyz[2]
                if skylayer < nxyz[2] and zmax > 0:
                    skylayer = int(nxyz[2] - 1 - skylayer)
    
    return skylayer
