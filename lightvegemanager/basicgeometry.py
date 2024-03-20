"""
    basicgeometry
    *************

    Provides basic geometric operations in 3D used by LightVegeManager. 
    Vectors are represented by their cartesian coordinates in a 3-tuple of floats such as v = (x, y, z)
    Triangles are a list of 3 vectors representing its vertices

"""

import numpy
import math


def crossproduct(v1, v2) :
    """Crossproduct between two vectors v1^v2

    :param v1: vector 1
    :type v1: 3-tuple
    :param v2: vector 2
    :type v2: 3-tuple
    :return: crossproduct
    :rtype: 3-tuple
    """    
    x = v1[1] * v2[2] - v2[1] * v1[2]
    y = v1[2] * v2[0] - v2[2] * v1[0]
    z = v1[0] * v2[1] - v2[0] * v1[1]
    return (x,y,z)

def middle(v1, v2) :
    """Middle point between v1 and v2

    :param v1: vector 1
    :type v1: 3-tuple
    :param v2: vector 2
    :type v2: 3-tuple
    :return: middle point :math:`p = (v1+v2)/2`
    :rtype: 3-tuple
    """    
    p = tuple([a+b for a,b in zip(v1, v2)])
    return (p[0]/2, p[1]/2, p[2]/2)

def triangle_normal(triangle) :
    """Computes normal of an oriented triangle (3D)

    .. note:: a triangle is ``[(x1,y1,z1),(x2,y2,z2),(x3,y3,z3)]``

    :param triangle: triangle represented by its 3 vertices
    :type triangle: list
    :return: normalized vector
    :rtype: 3-tuple
    """    
    side1 = [x-y for x, y in zip(triangle[1], triangle[0])]
    side2 = [x-y for x, y in zip(triangle[2], triangle[1])]
    side3 = [x-y for x, y in zip(triangle[0], triangle[2])]

    v12 = crossproduct(side1, side2)
    v23 = crossproduct(side2, side3)
    v31 = crossproduct(side3, side1)
    n = [x+y+z for x,y,z in zip (v12, v23, v31)]
    norm = math.sqrt(sum([c**2 for c in n]))
    return (n[0]/norm, n[1]/norm, n[2]/norm)

def triangle_elevation(triangle) :
    """Computes normal elevation of triangle

    .. note:: a triangle is ``[(x1,y1,z1),(x2,y2,z2),(x3,y3,z3)]``

    :param triangle: triangle represented by its 3 vertices
    :type triangle: list
    :return: 
        angle in degree 
        elevation starts from ground to normal vector
        must be between 0 and 90°
    
    :rtype: float
    """    
    n = triangle_normal(triangle)
    # compute normal elevation in radian
    e = math.acos(abs(n[2]))
    # converts in degree
    e *= 180/math.pi
    # elevation must be between 0 and 90
    if e > 90 and e < 180 : e = 90-(e%90)
    else : e = e%90
    
    return e

def triangle_area(triangle) :
    """triangle area

    .. note:: a triangle is ``[(x1,y1,z1),(x2,y2,z2),(x3,y3,z3)]``

    :param triangle: triangle represented by its 3 vertices
    :type triangle: list
    :return: area

        .. note:: algorithm is a copy of ``_surf`` in ``alinea.caribu.caributriangleset``
    
    :rtype: float
    """
    a, b, c = tuple(map(numpy.array, triangle))
    x, y, z = numpy.cross(b - a, c - a).tolist()
    return numpy.sqrt(x ** 2 + y ** 2 + z ** 2) / 2.0

def triangle_barycenter(triangle) :
    """triangle barycenter
    
    .. note:: a triangle is ``[(x1,y1,z1),(x2,y2,z2),(x3,y3,z3)]``
    
    :param triangle: triangle represented by its 3 vertices
    :type triangle: list
    :return: isobarycenter :math:`(s1 + s2 + s3)/3`
    :rtype: float
    """    
    return tuple([s/3 for s in [x+y+z for x,y,z in zip(*triangle)]])

def rescale(triangles, h) :
    """Multiplies all triangle vertices by a ratio h

    .. note:: a triangle is ``[(x1,y1,z1),(x2,y2,z2),(x3,y3,z3)]``

    :param triangles: list of triangles in the same format as below
    :type triangles: list of list
    :param h: ratio
    :type h: float
    :return: list of triangles with each triangle from input rescaled by h
    :rtype: list of list

    example
    -------

    >>> triangles = [[(0,0,0), (0,1,0), (0,1,1)], [(0,0,1), (0,0,0), (1,0,1)]]
    >>> rescale(triangles, 3)
    [[(0,0,0), (0,3,0), (0,3,3)], [(0,0,3), (0,0,0), (3,0,3)]]
    """
    return [[tuple([x*h for x in p]) for p in t] for t in triangles] 

def translate(triangles, tvec) :
    """Moves a list of triangles with a translation of vector tvec

    .. note:: a triangle is ``[(x1,y1,z1),(x2,y2,z2),(x3,y3,z3)]``

    :param triangles: list of triangles in the same format as below
    :type triangles: list of list
    :param tvec: vector
    :type tvec: 3-tuple
    :return: list of triangles with each triangle from input translated by tvec
    :rtype: list of list

    example
    -------

    >>> triangles = [[(0,0,0), (0,1,0), (0,1,1)], [(0,0,1), (0,0,0), (1,0,1)]]
    >>> tvec = (3,0,1)
    >>> translate(triangles, tvec )
    [[(3,0,1), (3,1,1), (3,1,2)], [(3,0,2), (3,0,1), (4,0,2)]]
    
    """
    return [[tuple([x+y for x,y in zip(p,tvec)]) for p in t] \
                                                    for t in triangles]

def zrotate(triangles, omegadeg) :
    """Rotates a list of triangles in the xy plane from an angle
    
    .. note:: a triangle is ``[(x1,y1,z1),(x2,y2,z2),(x3,y3,z3)]``
    
    :param triangles: list of triangles in the same format as below
    :type triangles: list of list
    :param omegadeg: angle in degree
    :type omegadeg: float
    :return: list of triangles with each triangle from input rotated around z axis by omegadeg
    :rtype: list of list

    example
    -------

    >>> triangles = [[(0,0,0), (0,1,0), (0,1,1)], [(0,0,1), (0,0,0), (1,0,1)]]
    >>> zrotate(triangles, 90 )
    [[(0,0,0), (-1,0,0), (-1,0,1)], [(0,0,1), (0,0,0), (0,1,1)]]

    """
    omegarad = omegadeg * math.pi/180
    newtriangles = []
    for t in triangles :
        newt = []
        for p in t :
            x = math.cos(omegarad)*p[0] - math.sin(omegarad)*p[1]
            y = math.sin(omegarad)*p[0] + math.cos(omegarad)*p[1]
            newt.append(tuple([x,y,p[2]]))
        newtriangles.append(newt)
    return newtriangles

