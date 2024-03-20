"""
    
    RATPinputs
    **********

    Manage vegetation and meteo input informations for RATP

    The argument `parameters` refers to one the three inputs dict of LightVegeManager. It is 
    structured as so:

    .. code:: python

        ratp_args = {
                # Grid specifications
                "voxel size" : [dx, dy, dz],
                "voxel size" : "dynamic",
                
                "origin" : [xorigin, yorigin, zorigin],
                "origin" : [xorigin, yorigin],

                "number voxels" : [nx, ny, nz],
                "grid slicing" : "ground = 0.",
                "tesselation level" : int,

                "full grid" : bool,

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

    .. seealso:: :ref:`Inputs description <inputs>`

"""
from math import *

def RATP_vegetation(parameters, angle_distrib, reflected):
    """Initialise a RATP Vegetation object from LightVegeManager input datas

    :param parameters: RATP parameters from inputs of LightVegeManager
    :type parameters: dict
    :param angle_distrib: leaf angle distribution
    :type angle_distrib: dict
    :param reflected: if the user wishes to activate reflected radiations
    :type reflected: bool
    :return: Vegetation types contains clumoing effect ratio, leaf angle distribution and reflectance/transmittance of leaves for each specy
    :rtype: PyRATP.pyratp.vegetation.Vegetation
    """
    from alinea.pyratp.vegetation import Vegetation

    entities_param = []
    if parameters["angle distrib algo"] != "compute voxel":
        for id, mu_ent in enumerate(parameters["mu"]):
            if reflected:
                reflectance_coef = parameters["reflectance coefficients"][id]
            else:
                reflectance_coef = [0.0, 0.0]
                parameters["reflectance coefficients"].append(reflectance_coef)

            entities_param.append({"mu": mu_ent, "distinc": angle_distrib["global"][id], "rf": reflectance_coef})

        return Vegetation.initialise(entities_param)
    else:
        for id, mu_ent in enumerate(parameters["mu"]):
            if reflected:
                reflectance_coef = parameters["reflectance coefficients"][id]
            else:
                reflectance_coef = [0.0, 0.0]
            entities_param.append({"mu": mu_ent, "rf": reflectance_coef})

        return Vegetation.initialise(entities_param, pervoxel=True, distribvox=angle_distrib["voxel"])


def RATP_meteo(energy, day, hour, coordinates, parunit, truesolartime, direct, diffus):
    """Initialise a RATP MicroMeteo object from LightVegeManager input parameters

    :param energy: input ray energy
    :type energy: float
    :param day: input day
    :type day: int
    :param hour: input hour
    :type hour: int
    :param coordinates: [latitude, longitude, timezone]
    :type coordinates: list
    :param parunit: unit of energy argument, "micromol.m-2.s-1" or else
    :type parunit: string
    :param truesolartime: activates true solar time or local time to compute sun position
    :type truesolartime: bool
    :param direct: if direct rays are activated
    :type direct: bool
    :param diffus: if diffuse rays are activated
    :type diffus: bool
    :return: input meteorological data at current time step
    :rtype: PyRATP.pyratp.micrometeo.MicroMeteo
    """
    from alinea.pyratp.micrometeo import MicroMeteo

    # RATP expects W.m-2 for input energy
    if parunit == "micromol.m-2.s-1":
        #: Spitters's model estimating for the diffuse:direct ratio
        # coefficient 2.02 : 4.6 (conversion en W.m-2) x 0.439 (PAR -> global)
        RdRs = RdRsH(Rg=energy / 2.02, DOY=day, heureTU=hour, latitude=coordinates[0])

        # coeff 4.6 : https://www.researchgate.net/post/Can-I-convert-PAR-photo-active-radiation-value-of-micro-mole-M2-S-to-Solar-radiation-in-Watt-m2
        energy = energy / 4.6  # W.m-2
    else:
        RdRs = RdRsH(Rg=energy / 0.439, DOY=day, heureTU=hour, latitude=coordinates[0])

    # PAR et Dif en W.m^-2
    if diffus:
        # Direct et diffus
        if direct:
            rdif = energy * RdRs
        # uniquement diffus
        else:
            rdif = energy
    # uniquement direct
    else:
        rdif = 0.0

    return MicroMeteo.initialise(doy=day, hour=hour, Rglob=energy, Rdif=rdif, truesolartime=truesolartime)

############ G Louarn - adaptation de spitters.c EGC grignon
def RdRsH(Rg, DOY, heureTU, latitude):
    """ fraction diffus/Global en fonction du rapport Global(Sgd)/Extraterrestre(Sod)- pas de temps horaire """
    
    def DecliSun(DOY):
        """ Declinaison (rad) du soleil en fonction du jour de l'annee """
        alpha = 2 * pi * (DOY - 1) / 365
        return (0.006918 - 0.399912 * cos(alpha) + 0.070257 * sin(alpha))
    
    hrad = 2 * pi / 24 * (heureTU - 12)
    lat = radians(latitude)
    dec = DecliSun(DOY)
    costheta = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(hrad)
    Io = 1370 * (1 + 0.033 * cos(2 * pi * (DOY - 4) / 366))#eclairement (w/m2) a la limitte de l'atmosphere dans un plan perpendiculaire aux rayons du soleil, fonction du jour
    So = Io * costheta #eclairement dans un plan parallele a la surface du sol
    RsRso = Rg / So
    R = 0.847 - 1.61 * costheta + 1.04 * costheta * costheta
    K = (1.47 - R) / 1.66
    
    if (RsRso <= 0.22) :
        return(1)
    elif (RsRso <= 0.35) :
        return(1 - 6.4 * (RsRso - 0.22)**2)
    elif (RsRso <= K) :
        return(1.47 - 1.66 * RsRso)
    else:
        return(R)