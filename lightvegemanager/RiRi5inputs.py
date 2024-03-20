"""
    RiRi5inputs
    **********

    Translate RATP grid in RiRi5 format
"""
import numpy

def ratpgrid_to_riri5(ratpgrid):
    """Transform a RATP grid of voxels into a RiRi5 grid format

    :param ratpgrid: grid of voxels from RATP
    :type ratpgrid: pyratp.grid
    :return: grid of voxels containing leaf area in each voxel for each specy
    :rtype: numpy.array
    """    
    la = numpy.zeros([ratpgrid.nent, ratpgrid.njz + 1, ratpgrid.njy, ratpgrid.njx])
    for ne in range(ratpgrid.nent):
        for ix in range(ratpgrid.njx):
            for iy in range(ratpgrid.njy):
                for iz in range(ratpgrid.njz):
                    k = ratpgrid.kxyz[ix, iy, iz]
                    if k > 0 :
                        la[ne][iz+1][ratpgrid.njy - (iy+1)][ix] = ratpgrid.leafareadensity[ne, k]

    for ne in range(ratpgrid.nent):
        zeros_array = numpy.zeros([ratpgrid.njy, ratpgrid.njx])
        la[ne][0][:][:] = zeros_array

    return la
