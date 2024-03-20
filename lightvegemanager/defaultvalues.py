"""
    defaultvalues
    *************

    Default for simulation fixed parameters
    These parameters are sets in LightVegeManager if the user did not precise thoses inputs
    
"""


def default_LightVegeManager_inputs():
    """returns default parameters for LightVegeManager constructor ``__init__``

    :return: environment and light model parameters
    :rtype: dict, dict
    """
    default_environnement = {}
    default_environnement["coordinates"] = [46.4, 0.0, 1.0]  # INRAE Lusignan
    default_environnement["sky"] = "turtle46"
    default_environnement["diffus"] = True
    default_environnement["direct"] = True
    default_environnement["reflected"] = False
    default_environnement["infinite"] = False

    default_ratp_parameters = {}
    default_ratp_parameters["voxel size"] = "dynamic"
    default_ratp_parameters["soil reflectance"] = [0.0, 0.0]
    default_ratp_parameters["reflectance coefficients"] = []
    default_ratp_parameters["mu"] = []
    default_ratp_parameters["tesselation level"] = 0
    default_ratp_parameters["angle distrib algo"] = "compute global"
    default_ratp_parameters["nb angle classes"] = 9
    default_ratp_parameters["full grid"] = False

    default_caribu_parameters = {}
    default_caribu_parameters["caribu opt"] = {"par": (0.10, 0.05)}
    default_caribu_parameters["sun algo"] = "caribu"
    default_caribu_parameters["soil mesh"] = -1
    default_caribu_parameters["debug"] = False

    default_riri5_parameters = {}
    default_riri5_parameters["voxel size"] = [1., 1., 1.]

    return default_environnement, default_ratp_parameters, default_caribu_parameters, default_riri5_parameters
