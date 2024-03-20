"""
    tesselator
    ****************************

    Handles tesselation of a triangulation, the subdivision of triangulation.
    The tesselation was initially made to have a better matching of a triangle mesh in a grid of voxels
    and avoid triangle to be located in several voxels

    Tesselation operates as a recursive function until a certain level is reached
"""
from lightvegemanager.basicgeometry import middle

## NOT USED, PyRATP.pyratp.grid.grid_index is preferred
def whichvoxel(p, mygrid):
    """Returns in which voxel p is located

    :param p: 3D vector located in the grid
    :type p: list or tuple
    :param mygrid:  RATP grid
    :type mygrid: pyratp.grid
    :return: returns the voxels where p is located in a list [kx, ky, kz] where ki is the indice
        voxel on each axis, [0, 0, 0] is origin voxel
    :rtype: list
    """
    vox = [(int)(p[0] / mygrid.dx)]
    vox.append((int)(p[1] / mygrid.dy))

    # recherche de le couche en z
    jz = -10
    arret = 0
    j = 1
    while (j < mygrid.njz) and (arret == 0):
        if p[2] > mygrid.dz[j]:
            arret = 1
        j += 1
    if arret == 1:
        jz = j - 1
    else:
        jz = mygrid.njz - 1
    vox.append(jz)
    return vox

## NOT USED, using impleted python functions is faster
def samevoxel(voxels):
    """Check if all voxels in list voxels are the same

    :param voxels: list of voxels in format [kx, ky, kz]
    :type voxels: list
    :return: if all elements of voxels are the same
    :rtype: bool
    """
    # test si tous les sommets sont dans le même triangle
    test = 0
    for i in range(len(voxels)):
        for j in range(3):
            test += abs((float)(voxels[i][j] - voxels[(i + 1) % 3][j]))

    # triangle dans voxel
    return test == 0


def tesselate(mygrid, triangle):
    """Check if the triangle is strictly inside a voxel or between several voxels. If so, it cuts out
    triangle in 4 smaller triangles.

    :param mygrid: RATP voxels grid
    :type mygrid: pyratp.grid
    :param triangle: triangle represented by its 3 vertices :code:`triangle = [(x1, y1, z1), (x2, y2, z2), (x3, y3, z3)]`
    :type triangle: list of tuple
    :return:

        - if the triangle is inside a voxel, it returns the input triangle in a 1-list
        - else if the triangle is between several voxels, it cuts out triangle in four smaller triangles from its vertices and barycenter
    :rtype: list of triangle
    """
    from alinea.pyratp import grid

    # Get voxels indices where the triangle vertices are located
    wh = []
    for i in range(3):
        # wh.append(whichvoxel(triangle[i], mygrid))
        Jx, Jy, Jz = grid.grid_index([triangle[i][0]], [triangle[i][1]], [triangle[i][2]], mygrid, toric=False)
        wh.append([Jx[0], Jy[0], Jz[0]])

    # tests if all the vertices are located in the same voxel
    if wh.count(wh[0]) == len(wh):
        return [triangle]

    # we subdivide in 4 triangles
    else:
        smalltriangles = []

        # compute vertices middles
        middles = [middle(triangle[i], triangle[(i + 1) % 3]) for i in range(3)]

        # triangle 1
        smalltriangles.append([triangle[0], middles[0], middles[2]])

        # triangle 2
        smalltriangles.append([middles[0], triangle[1], middles[1]])

        # triangle 3
        smalltriangles.append([middles[1], triangle[2], middles[2]])

        # triangle 4
        smalltriangles.append([middles[0], middles[1], middles[2]])

        return smalltriangles


def tesselate2(triangle):
    """Same as :func:tesselate but without the grid matching

    :param triangle: triangle represented by its 3 vertices :code:`triangle = [(x1, y1, z1), (x2, y2, z2), (x3, y3, z3)]`
    :type triangle: list of tuple
    :return: it cuts out triangle in four smaller triangles from its vertices and barycenter
    :rtype: list of triangle
    """
    smalltriangles = []

    # compute vertices middles
    middles = [middle(triangle[i], triangle[(i + 1) % 3]) for i in range(3)]

    # triangle 1
    smalltriangles.append([triangle[0], middles[0], middles[2]])

    # triangle 2
    smalltriangles.append([middles[0], triangle[1], middles[1]])

    # triangle 3
    smalltriangles.append([middles[1], triangle[2], middles[2]])

    # triangle 4
    smalltriangles.append([middles[0], middles[1], middles[2]])

    return smalltriangles


def iterate_trianglesingrid(triangle, mygrid, level, levelmax, triangles_shape):
    """Recursion on a triangle, while its subdivision doesn't match an RATP grid or exceed a certain level of tesselation

    :param triangle: triangle represented by its 3 vertices :code:`triangle = [(x1, y1, z1), (x2, y2, z2), (x3, y3, z3)]`
    :type triangle: list of tuple
    :param mygrid: RATP voxels grid
    :type mygrid: pyratp.grid
    :param level: current level of tesselation, i.e. how many times we cut out the initial triangle
    :type level: int
    :param levelmax: final level of tesselation to not be exceeded
    :type levelmax: int
    :param triangles_shape: list of triangles used to stock new triangles or input ``triangle`` if no tesselation is computed
    :type triangles_shape: list of triangles
    :return: ``triangles_shape`` updates with input ``triangle`` or new triangles if tesselation was activated
    :rtype: list of triangles
    """
    level += 1
    ltriangle = tesselate(mygrid, triangle)
    if len(ltriangle) == 1:
        triangles_shape.append(ltriangle[0])
    else:
        if level == levelmax:
            triangles_shape.append(triangle)

        else:
            for subt in ltriangle:
                iterate_trianglesingrid(subt, mygrid, level, levelmax, triangles_shape)


def iterate_triangles(triangle, level, levelmax, triangles_shape):
    """Subdivide a triangle until its tesselation reaches levelmax. It stocks the new triangles in ``triangle_shape``
    Recursive function

    :param triangle: triangle represented by its 3 vertices :code:`triangle = [(x1, y1, z1), (x2, y2, z2), (x3, y3, z3)]`
    :type triangle: list of tuple
    :param level: current level of tesselation, i.e. how many times we cut out the initial triangle
    :type level: int
    :param levelmax: final level of tesselation to not be exceeded
    :type levelmax: int
    :param triangles_shape: list of triangles used to stock new triangles or input ``triangle`` if no tesselation is computed
    :type triangles_shape: list of triangles
    :return: ``triangles_shape`` updates with new triangles from tesselation
    :rtype: list of triangles
    """
    level += 1
    ltriangle = tesselate2(triangle)
    if level == levelmax:
        triangles_shape.append(triangle)

    else:
        for subt in ltriangle:
            iterate_triangles(subt, level, levelmax, triangles_shape)
