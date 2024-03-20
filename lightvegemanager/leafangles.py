"""
    leafangles
    **********

    Handles leaf angle distribution, both in its dynamic computing or reading a file

    An angle distribution is a list of n elements, where n is the number of angle class between 0 and
    90°. Each element of this list is a percentage for leaves to be oriented from 0 to 90/n degree.
    The sum of all the list entries must equal 1.

"""
import numpy
import bisect

from lightvegemanager.basicgeometry import triangle_area, triangle_elevation

def read_distrib_file(path, numberofentities):
    """Reads global leaf angle distribution in a file
    the file must matches the format : one specy distribution per line
    each percentage separated by a coma ','

    **example**

    .. code-block:: bash

        0.1382,0.1664,0.1972,0.1925,0.1507,0.0903,0.0425,0.0172,0.005
        0.3,0.1,0.15,0.2,0.05,0.2

    :param path: path to the file
    :type path: string
    :param numberofentities: number of species, aka numbre of lines in the file
    :type numberofentities: int
    :return: distribution for each specy, number of entries on one line
    :rtype: list of list
    """
    f_angle = open(path, "r")
    distrib = []
    for i in range(numberofentities):
        line = f_angle.readline()
        distrib.append([float(x) for x in line.split(",")[:]])
    f_angle.close()

    return distrib


def compute_distrib_globale(trimesh, matching_ids, numberofclasses):
    """Calculation of a global leaf angle distribution from a triangle mesh

    :param trimesh: triangles mesh aggregated by indice elements :code:`{ id : [triangle1, triangle2, ...]}`
    :type trimesh: dict
    :param matching_ids:
        dict that matches new element indices in trimesh with specy indice and
        input element indice, :code:`matching_ids = { new_element_id : (input_element_id, specy_id)}`

        this dict allows us to look how species there is the inputs geometric data

    :type matching_ids: dict
    :param numberofclasses: number angle class wanted between 0 and 90 degree
    :type numberofclasses: int

    :return: a leaf angle distribution for each specy. Each distribution is a length numberofclasses the distribution is computed on all trimesh
    :rtype: list of list
    """
    from lightvegemanager.trianglesmesh import triangles_entity

    # dimensions [specy][angle class]
    distrib = []

    angles = list(numpy.linspace(90 / numberofclasses, 90, numberofclasses))

    # number of species
    numberofentities = max([v[1] for v in matching_ids.values()]) + 1

    for e in range(numberofentities):
        tri_entity = triangles_entity(trimesh, e, matching_ids)
        area_entity = sum([triangle_area(t) for t in tri_entity])

        classes = [0] * numberofclasses
        if area_entity > 0:
            for t in tri_entity:
                # id as angles[id-1] < triangle_elevation(t) <= angles[id]
                id_classe = bisect.bisect_left(angles, triangle_elevation(t))
                classes[id_classe] += triangle_area(t)

            classes = [a / area_entity for a in classes]
        distrib.append(classes)

    return distrib


def compute_distrib_voxel(trimesh, matching_ids, numberofclasses, numberofvoxels, matching_tri_vox):
    """Calculation of a local leaf angle distribution from a triangle mesh on each voxel of a grid

    :param trimesh: triangles mesh aggregated by indice elements :code:`{ id : [triangle1, triangle2, ...]}`
    :type trimesh: dict
    :param matching_ids:
        dict that matches new element indices in trimesh with specy indice and
        input element indice, :code:`matching_ids = { new_element_id : (input_element_id, specy_id)}`
        this dict allows us to look how species there is the inputs geometric data

    :type matching_ids: dict
    :param numberofclasses: number angle class wanted between 0 and 90 degree
    :type numberofclasses: int
    :param numberofclasses: number of non empty voxels in the grid
    :type numberofclasses: int
    :param matching_tri_vox: dict where key is a triangle indice and value the matching voxel indice where the
        barycenter of the triangle is located
    :type matching_tri_vox: dict
    :return: array of dimension [number of voxels][number of species][number of angle classes]
        i.e. it returns for each voxel a leaf angle distribution for each specy
    :rtype: numpy.array
    """
    from lightvegemanager.trianglesmesh import globalid_to_elementid, globalid_to_triangle

    angles = list(numpy.linspace(90 / numberofclasses, 90, numberofclasses))

    numberofentities = max([v[1] for v in matching_ids.values()]) + 1

    # dimensions [#voxel][#specy][#angle class]
    distrib = numpy.zeros([numberofvoxels, numberofentities, numberofclasses])

    # loops follow distrib dimensions
    for k in range(numberofvoxels):
        # triangles list in the current voxel
        tri_in_vox = [int(id) for id, v in matching_tri_vox.items() if int(v) == k]
        # sorting by specy
        for n in range(numberofentities):
            shape_ent = [ke for ke, v in matching_ids.items() if v[1] == n]
            # select triangles belonging to specy n
            tri_in_vox_in_ent = [i for i in tri_in_vox if globalid_to_elementid(trimesh, i) in shape_ent]

            # if current specy is present in the current voxel
            if tri_in_vox_in_ent:
                area_vox_ent = 0
                for id in tri_in_vox_in_ent:
                    t = globalid_to_triangle(trimesh, id)
                    id_classe = bisect.bisect_left(angles, triangle_elevation(t))

                    at = triangle_area(t)
                    distrib[k][n][id_classe] += at
                    area_vox_ent += at

                distrib[k][n] = distrib[k][n] / area_vox_ent
    return distrib
