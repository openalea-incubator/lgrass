"""
    transfer
    ********

    Handles transfer of LightVegeManager results to plant Models.

    Currently, it manages CN-Wheat and l-egume. Only l-egume needs additionnal processes to
    converts Dataframe results in a compatible format.

    l-egume expects the absorb PAR either per plant or locally following a grid of voxels, and the transmitted PAR
    locally following a grid of voxels. Here, we will convert and transform results for CARIBU and RATP to those l-egume format.
"""
import numpy
import scipy

def transfer_ratp_legume(m_lais, energy, ratp_grid, voxels_outputs, nb0, epsilon=1e-8):
    """Transfers LightVegeManager outputs from RATP to l-egume
    Absorb and transmitted PAR will follow a RATP grid of voxels, matching the dimensions of intern l-egume grid of voxels.

    :param m_lais: leaf area represented in a numpy.array of dimension
        [number of species, number of z layers, number of y layers, number of x layers]
    :type m_lais: numpy.array
    :param energy: input energy
    :type energy: float
    :param ratp_grid: RATP grid of voxels
    :type ratp_grid: pyratp.grid
    :param voxels_outputs: results from LightVegeManager
    :type voxels_outputs: pandas.Dataframe
    :param nb0: number of empty layers from top of the canopy and maximum z layers in m_lais
    :type nb0: int
    :param epsilon: criteria of minimum intercepted portion of PAR in a non empty voxel
    :type epsilon: float
    :return:

        two array with lighting informations

            * ``res_abs_i``: absorb PAR in each voxel for RATP grid of voxels. It has the same dimensions as ``m_lais``

            * ``res_trans``: transmitted PAR in each voxel, i.e. the energy leaving the voxels from input rays. This value is not dependent on specy

        dimensions are (number of z layers, number of y layers, number of x layers)

    :rtype: numpy.array, numpy.array
    """
    # initialize absorb energy array
    res_abs_i = numpy.zeros((m_lais.shape[0], m_lais.shape[1], m_lais.shape[2], m_lais.shape[3]))

    # voxel top side area
    dS = ratp_grid.dx * ratp_grid.dy
    res_trans = numpy.ones((m_lais.shape[1], m_lais.shape[2], m_lais.shape[3]))
    # maximum transmitted energy is total incoming energy per area
    res_trans = res_trans * (energy * dS)

    for ix in range(m_lais.shape[3]):
        for iy in range(m_lais.shape[2]):
            for iz in range(ratp_grid.njz):
                legume_iz = iz + nb0

                condition_x = voxels_outputs.Nx == m_lais.shape[2] - iy
                vox_data = voxels_outputs[condition_x & (voxels_outputs.Ny == ix + 1) & (voxels_outputs.Nz == iz + 1)]
                if not vox_data.empty :
                    a = min(sum(vox_data["Transmitted"]), dS)
                    res_trans[legume_iz, iy, ix] = energy * a

                s_entity = 0
                for k in range(m_lais.shape[0]):
                    s_entity += m_lais[k][legume_iz][iy][ix]

                if s_entity > 0.0:
                    for ie in range(m_lais.shape[0]):
                        if len(vox_data) > 0:
                            v_dat = vox_data[vox_data.VegetationType == ie + 1]
                            v = v_dat["Intercepted"].values[0]
                            if v > epsilon:
                                res_abs_i[ie, legume_iz, iy, ix] = energy * v

                            # if a voxel has leaf area > 0, it must have a minimum intercepted energy value
                            else:
                                res_abs_i[ie, legume_iz, iy, ix] = epsilon

    return res_abs_i, res_trans


def transfer_caribu_legume(
    energy,
    skylayer,
    id,
    elements_outputs,
    sensors_outputs,
    sensors_dxyz,
    sensors_nxyz,
    m_lais,
    list_invar,
    list_lstring,
    list_dicFeuilBilanR,
    infinite,
    epsilon,
):
    """Transfers LightVegeManager outputs from CARIBU to l-egume
    We will update list_invar which stores the total intercepted energy for each plant, and return
    an array storing transmitted energy following the intern grid of voxels in l-egume. To do so, we
    used virtual sensors in CARIBU to get incoming radiations in selected locations.

    :param energy: input energy
    :type energy: float
    :param skylayer: number of empty layers from top of the canopy and maximum z layers in m_lais
    :type skylayer: int
    :param id: list of indices of input scenes associated with l-egume
    :type id: list of int
    :param elements_outputs: Dataframe results of elements formatted by :func:outputs.out_caribu_elements
    :type elements_outputs: pandas.Dataframe
    :param sensors_outputs: lighting results of virtual sensors form CARIBU in the format for each bandwidth computed,
        
        .. code-block::

            sensors_outputs[band+" Eabs"] = {sensor_id : energy}
            sensors_outputs[band+" Ei"] = {sensor_id : energy}

    :type sensors_outputs: dict of list
    :param sensors_dxyz: size of sides of a voxel in the grid of virtual sensors [dx, dy, dz]
    :type sensors_dxyz: list
    :param sensors_nxyz: number of sensors in each direction in the grid [nx, ny, nz]
    :type sensors_nxyz: list
    :param m_lais: leaf area represented in a numpy.array of dimension
        [number of species, number of z layers, number of y layers, number of x layers]
    :type m_lais: numpy.array
    :param list_invar: from l-egume, each element corresponds to an input specy of l-egume. Each element is a dict invar stores instant intern variables of l-egume.
    :type list_invar: list of dict
    :param list_lstring: from l-egume, each element corresponds to an input specy of l-egume. Each element is a dict lstring stores the l-system of each plant
    :type list_lstring: list of dict
    :param list_dicFeuilBilanR: from l-egume, each element corresponds to an input specy of l-egume. Each element is a dict dicFeuiBilanR stores correspondances between voxels grid and each plant
    :type list_dicFeuilBilanR: list of dict
    :param infinite: if the user wishes to activate infinitisation of the grid
    :type infinite: bool
    :param epsilon: criteria of minimum intercepted portion of PAR in a voxel (if res_abs_i is zero, l-egume will crash)
    :type epsilon: float
    :raises ValueError: Virtual sensors and finite scene doesn't work yet with CARIBU
    :return:

        * it updates ``list_invar`` and its key entries ``"parap"`` and ``"parip"``, each element if the scipy.array is the sum of all intercepted energy for each plant. This process is a rewrite of ``calc_paraF`` in ShootMorpho.py module of l-egume, adapted to LightVegeManager numerotation of triangles

        * ``res_trans`` an array of transmitted energy for each voxel in a grid of dimensions ``sensors_dxyz * sensors_nxyz``

    :rtype: numpy.array
    """
    ## Absorb radiations for each plant of each specy
    for k in range(len(list_invar)):
        # initialize absorb energy
        nplantes = len(list_invar[k]["Hplante"])
        list_invar[k]["parap"] = scipy.array([0.0] * nplantes)
        list_invar[k]["parip"] = scipy.array([0.0] * nplantes)

        if id == None:
            filter = elements_outputs.VegetationType == k
            ent_organs_outputs = elements_outputs[filter]
        elif type(id) == list or type(id) == tuple:
            filter = elements_outputs.VegetationType == id[k]
            ent_organs_outputs = elements_outputs[filter]

        # non empty scene
        for i in range(len(ent_organs_outputs)):
            organe_id = int(ent_organs_outputs.iloc[i]["Organ"])

            # PAR in W/m²
            par_intercept = ent_organs_outputs.iloc[i]["par Ei"] * energy
            S_leaf = ent_organs_outputs.iloc[i]["Area"]

            id_plante = list_lstring[k][organe_id][0]
            p_s = par_intercept * S_leaf
            a = float(list_invar[k]["parip"][id_plante])
            list_invar[k]["parip"][id_plante] = a + p_s

            # we remove senescent leaves
            if list_lstring[k][organe_id][9] != "sen":
                a = float(list_invar[k]["parap"][id_plante])
                list_invar[k]["parap"][id_plante] = a + p_s

        # all non empty plant must have a minimum intercepted energy
        if len(list_invar[k]["parip"]) == len(list_dicFeuilBilanR[k]["surf"]) :
            for p in range(len(list_invar[k]["parip"])):
                if list_invar[k]["parip"][p] == 0.0 and list_dicFeuilBilanR[k]["surf"][p] > 0.0:
                    list_invar[k]["parip"][p] = epsilon

        # conversion
        c = (3600 * 24) / 1000000
        list_invar[k]["parap"] *= c
        list_invar[k]["parip"] *= c

    ## Transmitted radiations throughout a grid of voxels
    res_trans = numpy.ones((m_lais.shape[1], m_lais.shape[2], m_lais.shape[3]))

    # if non empty scene
    if not elements_outputs.empty:
        # different treatment if scene is infinite, (issues with virtual sensors and finite scene with CARIBU)
        if infinite:
            nb0 = min(m_lais.shape[1] - sensors_nxyz[2], m_lais.shape[1])
            ID_capt = 0
            for ix in range(sensors_nxyz[0]):
                for iy in range(sensors_nxyz[1]):
                    for iz in range(sensors_nxyz[2] - skylayer):
                        a = min(sensors_outputs["par"][ID_capt], 1.)
                        res_trans[nb0 + ((sensors_nxyz[2] - 1)) - iz][iy][ix] = a
                        ID_capt += 1

        else:
            raise ValueError("CARIBU Sensor + no infinite -> Doesn't work yet")

    # gives maximum transmitted energy
    dS = sensors_dxyz[0] * sensors_dxyz[1]
    res_trans = res_trans * energy * dS

    return res_trans
