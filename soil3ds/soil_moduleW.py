'''

 3DS-soil-model, a 3D Soil model adapted from STICS : soil water balance module
 ******************************
 Authors: G. Louarn 
 

 

'''



from copy import deepcopy
import numpy as np
from numpy import array, multiply, divide, sum
from soil3ds.miscel_functions import * #soil3ds miscellaneous soil functions
from soil3ds.plt_functions import * #soil3ds miscellaneous plant functions


#from openalea.plantgl.all import *
#import sys
#path_ = r'H:\devel\grassland\grassland\luzerne'
#sys.path.insert(0, path_)
#from Obj3Dutils import *


##changements version 5
#ajout d'un patern en entree pour initialisation + stocke info surfsolref (facilite couplage racine dans impulse)-> pas d'action sur calculs
# ajout d'une entree obstarac (matrice 2D de profondeur d'obtacle en m ; None par defaut) et d'une matrice d'effet sur la disponibilite des ressource (m_obstrac)
# mise a jour d calcul de FTSW pour tenir compte de m_obstrac
# retire des ; dans le calcul de transpiration
# defini une fonction test_uni() pour faire un test uitaire
# retire dependence a opealea: fonctions de visu dans un module supplementaire:soil_modulevisu1.py

#A faire: retirer


## changements version 4
# distinguer les entree irrigation et rain pour les comptabiliser dans le bilan -> OK
# reprendre distrib_water_uptakeNC pour sortir les uptake par plante -> besoin pour les prelv d'N par plante -> OK



#a faire: map_Evap0 devrait etre une entree (pas evapo_tot) -> encore calcule!
#a faire: k_eq utilise pour evaporation: fixe = a sortir en parametre!
# pour eviter effet de marche de ftsw qd racine passe ds un nveau compartiment: considerer fraction du volume comme pour ZESX?
# gestion du drainage: 2 passe des precipitation -> simplifier (1 passe seulement suffit?) pour gere infiltration des nitrates?
# evapo_tot -> renvoie la somme (pas la moyenne)





class Soil(object):
    """ Main class for the soil object of the 'soil3ds' model 
        
    The class :class:`soil3ds.soil_moduleW.Soil` includes descriptions of soil physical properties and methods to compute water balance adpated for a 3D grid from STICS model

    :param par_sol: A dictionnary defining water humidity thresholds for each soil type necessary to build the Soil object
    :type par_sol: dict
    :param soil_number: List of nz soil type IDs along soil depth, default to [1]*10
    :type soil_number: list
    :param dxyz: List of three lists indicating voxel discretization along [x,y,z] axes (unit: m), default to 1 m3 in 10 vertical layers
    :type dxyz: list
    :param vDA: List of nz soil bulk density along soil depth (unit: g.cm-3), default to [1.2]*10
    :type vDA: list
    :param ZESX: STICS parameter ZESX - (unit: m), default to 0.3
    :type ZESX: float
    :param CFES: STICS parameter CFES - (unitless), default to 1.
    :type CFES: float
    :param obstarac: List of nz voxel fractions affected by 'obstarac' along soil depth, default to None (no root obstacles)
    :type obstarac: list
    :param pattern8: list of two [x,y] points defining the soil limits (unit: cm), default to [[0,0], [100,100]]
    :type pattern8: list


    Soil object Attributes:
        * :dxyz (list): Copy of dxyz input (unit: m)
        * :pattern (list): Copy of pattern8 input (unit: cm)
        * :obstarac (list): Copy of obstarac input
        * :size (list): List [nz,nx,ny] indicating 3D grid size
        * :origin (array): Origin [x,y,z] of the 3D grid (unit: m)
        * :surfsolref (float): Surface area of the soil object(unit: m2)
        * :corners (list): List of voxel corners along axes [x,y,z] (unit: m)
        
        * :m_1 (nd.array): Unitary array of size [nz,nx,ny]
        * :m_soil_vox (nd.array): Array of size [nz,nx,ny] storing voxel (x,y,z) ID
        * :m_soil_vol (nd.array): Array of size [nz,nx,ny] storing voxel volume (unit: m3)
        * :m_vox_surf (nd.array): Array of size [nz,nx,ny] storing voxel area (unit: m2)
        * :m_soil_number (nd.array):  Array of size [nz,nx,ny] storing voxel soil type ID
        * :m_DA (nd.array):  Array of size [nz,nx,ny] storing voxel soil apparent density (unit: g.cm-3)
        
        * :m_teta_fc (nd.array):  Array of size [nz,nx,ny] storing voxel soil volumic humidity at field capacity (0-1 fraction)
        * :m_teta_ad (nd.array):  Array of size [nz,nx,ny] storing voxel soil volumic humidity at minimal water content (0-1 fraction)
        * :m_teta_wp (nd.array):  Array of size [nz,nx,ny] storing voxel soil volumic humidity at wilting point (0-1 fraction)
        * :m_teta_sat (nd.array):  Array of size [nz,nx,ny] storing voxel soil volumic humidity at saturation (0-1 fraction)
        * :m_QH20fc (nd.array):  Array of size [nz,nx,ny] storing voxel water content at field capacity (unit: mm)
        * :m_QH20max (nd.array):  Array of size [nz,nx,ny] storing voxel available water content at field capacity (unit: mm)
        * :m_QH20wp (nd.array):  Array of size [nz,nx,ny] storing voxel water content at wilting point (unit: mm)
        * :m_QH20min (nd.array):  Array of size [nz,nx,ny] storing voxel water content at minimal water content (unit: mm)
        * :m_QH20sat (nd.array):  Array of size [nz,nx,ny] storing voxel water content at saturation (unit: mm)
        * :m_frac_evapZ (nd.array):  Array of size [nz,nx,ny] storing voxel relative verticale distribution of soil evaporation at field capacipy (0-1 fraction)
        * :m_corr_vol_evap (nd.array):  Array of size [nz,nx,ny] storing voxel volume fraction available to soil evaporation (0-1 fraction)
        * :m_obstarac (nd.array):  Array of size [nz,nx,ny] storing voxel proportion affected by 'obstarac' (root growth obstacle; 0-1 fraction)
        
        * :asw_t (nd.array):  Array of size [nz,nx,ny] storing voxel daily amount of available soil water (unit: mm)
        * :tsw_t (nd.array):  Array of size [nz,nx,ny] storing voxel daily total amount of soil water (unit: mm)
        * :ftsw_t (nd.array):  Array of size [nz,nx,ny] storing voxel daily fraction of available soil water (0-1 fraction)
        * :bilanW (dict):  Dictionnary storing soil Water balance daily outputs
        
    Main methods:

        * ``__init__``: Initializes and builds static Soil object for the rest of simulation
        * :meth:`build_teta_m`: Build an array [nz, nx, ny] of soil humidity threshold according to par_sol[key]
        * :meth:`compute_teta_lim`: Initialise the different thresholds of volumetric water from 'par_sol' parameters
        * :meth:`init_mat_obstarac`: Initialise 'm_obstarac' attribute with an array of size [nz,nx,ny]
        * :meth:`init_asw`: Initialise asw_t attribute with an array of size [nz,nx,ny]
        * :meth:`HRv`: Compute relative soil humidity from soil attributes
        * :meth:`HRp`: Compute ponderal water content from soil attributes
        * :meth:`HRv_capteur`: Compute volumetric water content from soil attributes
        * :meth:`update_asw`: Update of daily available soil water array 'asw_t'
        * :meth:`update_ftsw`: Update of daily fraction of available soil water array 'ftsw_t'
        * :meth:`get_vox_coordinates`: Get corner coordinates of voxel [z, x, y]
        * :meth:`ls_1storder_vox`: Determines the list of first order neigbhour voxels

    Water evaporation methods:

        * :meth:`Kzi`: Compute coefficients K(Z,I) - Water proportion available for evaporation
        * :meth:`distrib_frac_evap`: Initialise 'm_frac_evapZ' attribute
        * :meth:`corr_vol_voxel_evap`: Initialise 'm_corr_vol_evap' attribute
        * :meth:`distrib_evapSTICS`: Compute distribtion of evaporative demande per voxel
    
    
    Water infiltration methods:

        * :meth:`infil_layer`: Compute water infiltration for soil layer idz
        * :meth:`distrib_PI`: Distribute water infiltration into the Soil grid from rain and/or irrigation map
    
    
    Water uptake methods:

        * :meth:`water_uptakeVox`: Compute water uptake relative distribution for a list of individual plants in voxel [idz,idx,idy]
        * :meth:`distrib_water_uptakeNC`: Compute daily water partitioning between individual plants and the spatial distribution of plant water uptake
    

    Water balance and daily step methods:

        * :meth:`stepWBmc`: Compute daily step Soil water balance
        * :meth:`OpenWbalance`: Initialise bilanW attribute, a dictionary storing water balance information
        * :meth:`UpdateWbalance`: Daily update of bilanW attribute
        * :meth:`CloseWbalance`: Finalise calculation of whole simulation variables in bilanW attribute
        * :meth:`PrintWbalance`: Print a summary table of bilanW attribute

    """
    
    def __init__(self, par_sol, soil_number = [1]*10, dxyz = [[1.], [1.], [0.2]*10], vDA=[1.2]*10, ZESX=0.3, CFES=1., obstarac=None, pattern8=[[0,0],[100.,100.]]):
        """ Initializes and builds static Soil object for the rest of simulation

        """
        # """
        # dxyz : dimensions des voxels selon x, y, z (m)
        # m_soil_vol : volume des voxels (m3)
        # m_DA : densite apparente du sol par voxel (g.cm-3)
        # m_soil_number: index du type de sol par voxel(lien avec dictionnaire des parametres hydraulique du sol)
        # m_vox_surf: surface de sol par voxel (m2)
        # soilSurface: surface de sol (m2)
        #
        # asw_t : quantite d'eau libre pour plante dans le voxel au temps t (mm)
        # tsw_t : quantite d'eau totale dans le voxel au temps t (mm)
        # ftsw_t : fraction d'eau libre disponible pour les plantes (Transpirable) par voxel (0-1, sans dimension)
        # Kzi() : effet hydrique fraction d'eau disponible pour evaporation (Evaporable) par voxel (0-1, sans dimension)
        # HRv() : humidite volumique du sol (%)
        # HRp() : humidite ponderale du sol (%)
        #
        # m_QH20sat : quantite d'eau du sol a saturation par voxel (mm)
        # m_QH20fc : quantite d'eau du sol a la capacite au champs par voxel (mm)
        # m_QH20max : quantite d'eau libre dispo pour la plante  a la capacite au champs par voxel (mm)(field capacipi-wilting point)
        # m_QH20wp : quantite d'eau du sol restante au wilting point (non accessible a la plante) par voxel (mm)
        # m_QH20min : quantite d'eau du sol minimum a areet de l'evaporation (sol s'aseche pas en dessous) (mm)
        #
        # m_teta_sat : humidite volumique du sol a saturation par voxel (fraction 0-1)
        # m_teta_fc : humidite volumique du sol a la capacite au champs par voxel  (fraction 0-1)
        # m_teta_wp : humidite volumique du sol au wilting point par voxel  (fraction 0-1)
        # m_teta_ad : humidite volumique minimum du sol sec  (fraction 0-1)
        #
        # m_frac_evapZ : distribution verticale relative de evaporation a la capacite au champ (qd Kzi=1) (0-1, sans dimension)
        # m_corr_vol_evap : complement a la fraction de volume effectivement accessible a evaporation par voxel (0-1, sans dimension)
        #
        # """

        self.dxyz = dxyz
        size = list(map(len, dxyz))
        self.size = size[2:]+size[0:2] #z,x,y
        self.pattern = pattern8 #cm
        self.origin = array(pattern8[0]+[0.])/100. #m
        Lsol = max((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
        largsol = min((pattern8[1][0] - pattern8[0][0]) / 100., (pattern8[1][1] - pattern8[0][1]) / 100.)  # m
        self.surfsolref = Lsol * largsol  # m2
        self.soilSurface = sum(self.dxyz[0])*sum(self.dxyz[1]) #?? doublon? surfsolref?

        self.m_1 = []
        self.m_soil_vox = []#voxel positions
        self.m_soil_vol = []#voxel volume
        self.m_soil_number = []#strate de type de sol
        self.m_DA = []
        self.m_vox_surf = []
        for z in range(len(dxyz[2])):
            v, vvox, v2, v3, vd, vs = [], [], [], [], [], []
            for x in range(len(dxyz[0])):
                vv, vvvox, vv2, vv3, vvDA, vvs = [], [], [], [], [], []
                for y in range(len(dxyz[1])):
                    vv.append(1)
                    vol = dxyz[0][x] * dxyz[1][y] * dxyz[2][z]
                    dists = [sum(dxyz[0][0:x]), sum(dxyz[1][0:y]), sum(dxyz[2][0:z])]
                    vvvox.append(array(dists))
                    vv2.append(vol)
                    vv3.append(soil_number[z])
                    vvDA.append(vDA[z])
                    vvs.append(dxyz[0][x] * dxyz[1][y])

                v.append(vv)
                vvox.append(vvvox)
                v2.append(vv2)
                v3.append(vv3)
                vd.append(vvDA)
                vs.append(vvs)

            self.m_1.append(v)
            self.m_soil_vox.append(vvox)
            self.m_soil_vol.append(v2)
            self.m_soil_number.append(v3)
            self.m_DA.append(vd)
            self.m_vox_surf.append(vs)

        self.m_1, self.m_soil_vox, self.m_soil_vol, self.m_soil_number, self.m_DA, self.m_vox_surf = array(self.m_1), array(self.m_soil_vox), array(self.m_soil_vol), array(self.m_soil_number), array(self.m_DA), array(self.m_vox_surf)
               
        self.m_frac_evapZ = self.distrib_frac_evap(ZESX, CFES)
        self.m_corr_vol_evap = self.corr_vol_voxel_evap(ZESX)
        self.compute_teta_lim(par_sol)
        self.obstarac = obstarac
        self.init_mat_obstarac()
        self.init_asw(HRp_init=None)  # initialise par defaut a la capacite au champ...
        # self.update_ftsw()

        xs = np.cumsum(np.array(self.dxyz[0])) - self.dxyz[0][0] + self.pattern[0][0] / 100 #S.pattern  in cm #S.dxyz  # m
        ys = np.cumsum(np.array(self.dxyz[1])) - self.dxyz[1][0] + self.pattern[0][1] / 100
        zs = -np.cumsum(np.array(self.dxyz[2])) + self.dxyz[2][0]
        self.corners = [xs,ys,zs]#m

    
    
    #####  Main methods : initialise/ update / examine soil object #####
    

    def get_vox_coordinates(self, z, x, y):
        """ Get corner coordinates of voxel [z, x, y]

        :param z: z value of voxel ID in grid [nz, nx, ny]
        :type z: int
        :param x: x value of voxel ID in grid [nz, nx, ny]
        :type x: int
        :param y: y value of voxel ID in grid [nz, nx, ny]
        :type y: int

        :return: 3D coordinates of the voxel upper corner
        :rtype: numpy.array
        """
        #""" to get the corner coordinates of a voxel from his IDs"""
        return np.array([self.corners[0][x], self.corners[1][y], self.corners[2][z]])


    def build_teta_m(self, par_sol, key='teta_fc'):
        """ Build an array [nz, nx, ny] of soil humidity threshold according to par_sol[key]

        :param par_sol: Parameter dictionnary 'par_sol' defining water humidity thresholds in a Soil object
        :type par_sol: dict
        :param key: A key in par_sol dictionnary, default = 'teta_fc'
        :param key: string

        :return: Array of size [nz,nx,ny] storing voxel volumetric water content for threshold 'key'
        :rtype: numpy.array
        """
        m_teta = []  #m_teta = matrice de teneur en eau volumique
        for z in range(len(self.m_soil_number)):
            v = []
            for x in range(len(self.m_soil_number[z])):
                vv = []
                for y in range(len(self.m_soil_number[z][x])):
                    vv.append(par_sol[str(self.m_soil_number[z][x][y])][key])

                v.append(vv)
            m_teta.append(v)
        return array(m_teta)


    def compute_teta_lim(self, par_sol):
        """ Initialise the different thresholds of volumetric water content of the grid voxels from 'par_sol' parameters

        :param par_sol: Parameter dictionnary 'par_sol' defining water humidity thresholds in a Soil object
        :type par_sol: dict

        """
        self.m_teta_fc = self.build_teta_m(par_sol, key='teta_fc')
        self.m_teta_ad = self.build_teta_m(par_sol, key='teta_ad')
        self.m_teta_wp = self.build_teta_m(par_sol, key='teta_wp')
        self.m_teta_sat = self.build_teta_m(par_sol, key='teta_sat')
        self.m_QH20fc = multiply(self.m_soil_vol, self.m_teta_fc)*1000. #quantite totale max en mm d'eau dans le sol a la capacite au champ
        self.m_QH20max = multiply(self.m_soil_vol, self.m_teta_fc - self.m_teta_wp)*1000. #quantite max en mm dispo pour la plante (field capacipi-wilting point)
        self.m_QH20wp = multiply(self.m_soil_vol, self.m_teta_wp)*1000. #quantite eau restante au wp (non accessible a la plante)
        self.m_QH20min = multiply(self.m_soil_vol, self.m_teta_ad)*1000. #quantite min en mm dispo pour evaporation (sol s'aseche pas en dessous)
        self.m_QH20sat = multiply(self.m_soil_vol, self.m_teta_sat)*1000. #quantite totale max d'eau dans le sol a saturation 
        

    def init_asw(self, HRp_init=None):
        """ Initialise asw_t attribute with an array of size [nz,nx,ny]

        :param HRp_init: List of nz initial soil ponderal humidity by vertical layers, default to None imply firld capacity value
        :type HRp_init: list, optional

        """
        #""" par defaut a la capacite au champ """
        #'HRp_init : a 1D vector/list of initial HRp (%) in the vertical z axis'
        
        if HRp_init==None:
            #avant, faut avoir lance compute_teta_lim(par_sol)
            F = self.m_1*1.
            self.tsw_t = multiply(F, self.m_QH20fc) #cree matrice tsw_t
            self.asw_t = multiply(F, self.m_QH20max) #cree matrice asw_t
            self.update_ftsw()
            self.OpenWbalance() #initialise W balance
        else:#initialisation avec un vecteur d'humidite ponderales
            for z in range(len(self.dxyz[2])):
                self.tsw_t[z,:,:] = self.m_1[z,:,:]*HRp_init[z]*(10*self.m_DA[z,:,:])*self.m_soil_vol[z,:,:]

            self.update_asw()
            self.update_ftsw()
            self.OpenWbalance()#initialise W balance

    def update_asw(self):
        """ Update calculation of daily available soil water array 'asw_t'
    
        """
        self.asw_t = self.tsw_t - self.m_QH20wp

    def update_ftsw(self):
        """ Update calculation of daily fraction of available soil water array 'ftsw_t'
    
        """
        ftsw = self.asw_t / self.m_QH20max
        negs = ftsw >0.
        self.ftsw_t = negs * 1. *ftsw * self.m_obstarac # pour eviter valeurs negatives

    def HRv(self):
        """ Compute relative soil humidity from soil attributes - 100 = field capacity (unit: %)

        :return: Relative soil humidity

        """
        return divide(self.tsw_t, self.m_QH20fc)*100.
        # peut etre 100 pour cent?? -> OK par rapport aux parametetre de reponse de mineralisation a l'humidite du sol (proportion de capciate au champ)
        #ambigu car pas HRv mesure par capteurs sols!

    def HRv_capteur(self):
        """ Compute volumetric water content from soil attributes - 100 = total voxel volume (unit: %)

        :return: Volumetric water content

        """
        return (self.tsw_t / self.m_soil_vol)/10.
        # correspond a mesure par capteurs sols!

    def HRp(self):
        """ Compute ponderal water content from soil attributes - 100 = voxel dry soil weight (unit: g H2O.g dry soil-1)

        :return: Ponderal water content

        """
        #""" compute relative soil humidity - Humidite ponderale (%, 100 * g H2O.g sol-1)"""
        return (self.tsw_t / self.m_soil_vol)/ (self.m_DA * 10.)
        # correspond a mesure par pesees sols!

    def init_mat_obstarac(self):
        """ Initialise 'm_obstarac' attribute with an array of size [z,x,y] - Defines voxel volume proportion affected by 'obstarac' (reduced volume for roots)

        """
        
        if self.obstarac is None or type(self.obstarac)!=type(array([0.])):#pas d'obstarac
            mat_obstarac = 1. * self.m_1
        else: #obstrac = matrice 2D de valeurs negatives (m)
            ls_prof = [0.]
            for i in self.dxyz[2]:
                ls_prof.append(ls_prof[-1] - i)

            mat_obstarac = 1. * self.m_1
            for x in range(len(self.dxyz[0])):
                for y in range(len(self.dxyz[1])):
                    obstXY = self.obstarac[x, y]
                    for z in range(len(self.dxyz[2])):
                        if ls_prof[z] > obstXY >= ls_prof[z + 1]:
                            ratio = (obstXY - ls_prof[z]) / (ls_prof[z + 1] - ls_prof[z])
                            mat_obstarac[z, x, y] = ratio
                        elif obstXY > ls_prof[z + 1]:
                            mat_obstarac[z, x, y] = 0.

        self.m_obstarac  =  mat_obstarac


    def ls_1storder_vox(self, x,y,z, opt=1):
        """ Determines the list of first order neigbhour voxels for voxel [z, x, y] in the same layer z

        :param z: z value of voxel ID in grid [nz, nx, ny]
        :type z: int
        :param x: x value of voxel ID in grid [nz, nx, ny]
        :type x: int
        :param y: y value of voxel ID in grid [nz, nx, ny]
        :type y: int
        :param opt: If opt=1, bypass and returns only initial voxel ID

        :return: A list of voxel IDs

        """
        #""" definit liste de voxels concomittents de 1er ordre / considere sol thorique, opt=1 renvoie uniquement voxel pour ecoulement vertical direct"""
        ## por gerer infiltrations 
        
        dxyz = self.dxyz
        
        lx,ly = len(dxyz[0]), len(dxyz[1])
        if lx<3 or ly<3 or opt==1: #pas au moins 9 voxel ou ecoulement direct (opt=1)
            return [[x,y,z]]
        else: 
            if x==0: #bord gauche
                if y>0 and y<ly-1:#pas coin
                    return [[lx-1, y+1, z], [x, y+1, z], [x+1, y+1, z], [lx-1, y, z], [x, y, z], [x+1, y, z], [lx-1, y-1, z], [x, y-1, z], [x+1, y-1, z]]
                elif y==0:#coin bas
                    return [[lx-1, y+1, z], [x, y+1, z], [x+1, y+1, z], [lx-1, y, z], [x, y, z], [x+1, y, z], [lx-1, ly-1, z], [x, ly-1, z], [x+1, ly-1, z]]
                elif y==ly-1:#coin haut
                    return [[lx-1, 0, z], [x, 0, z], [x+1, 0, z], [lx-1, y, z], [x, y, z], [x+1, y, z], [lx-1, y-1, z], [x, y-1, z], [x+1, y-1, z]]
            elif x == lx-1:#bord droit
                if y>0 and y<ly-1:#pas coin
                    return [[x-1, y+1, z], [x, y+1, z], [0, y+1, z], [x-1, y, z], [x, y, z], [0, y, z], [x-1, y-1, z], [x, y-1, z], [0, y-1, z]]
                elif y==0:#coin bas
                    return [[x-1, y+1, z], [x, y+1, z], [0, y+1, z], [x-1, y, z], [x, y, z], [0, y, z], [x-1, ly-1, z], [x, ly-1, z], [0, ly-1, z]]
                elif y==ly-1:#coin haut
                    return [[x-1, 0, z], [x, 0, z], [0, 0, z], [x-1, y, z], [x, y, z], [0, y, z], [x-1, y-1, z], [x, y-1, z], [0, y-1, z]]
            elif y==0:#cote bas sans les coins
                return [[x-1, y+1, z], [x, y+1, z], [x+1, y+1, z], [x-1, y, z], [x, y, z], [x+1, y, z], [x-1, ly-1, z], [x, ly-1, z], [x+1, ly-1, z]]
            elif y==ly-1:#cote haut sans les coins
                return [[x-1, 0, z], [x, 0, z], [x+1, 0, z], [x-1, y, z], [x, y, z], [x+1, y, z], [x-1, y-1, z], [x, y-1, z], [x+1, y-1, z]]
            else: #au centre
                return [[x-1, y+1, z], [x, y+1, z], [x+1, y+1, z], [x-1, y, z], [x, y, z], [x+1, y, z], [x-1, y-1, z], [x, y-1, z], [x+1, y-1, z]]

        #ls_1storder_vox(S.dxyz, 0,0,0, opt='1')
    
    
    #####  soil evaporation functions #####

    def Kzi(self):
        """ Compute coefficients K(Z,I) - Water proportion available for evaporation -
        (STICS book, Eq. 7.6, p 129)

        :return: Array of size [nz,nx,ny] with K(Z,I) coefficients

        """
        # seuilmax = (self.m_QH20max - self.m_QH20min) * (1 - self.m_frac_evapZ) + self.m_QH20min
        # return (self.tsw_t - seuilmax) / (self.m_QH20fc - seuilmax)
        return (self.tsw_t - self.m_QH20min) / (self.m_QH20fc - self.m_QH20min)


    def distrib_frac_evap(self, ZESX=0.3, CFES=1.):
        """ Initialise 'm_frac_evapZ' attribute with relative vertical distribution of soil evaporation at field capacity (i.e. Kzi=1) -
        (STICS book, Eq. 7.6 p 129) - integration with  dz = 1cm

        :param ZESX: STICS parameter ZESX - (unit: m)
        :type ZESX: float
        :param CFES: STICS parameter CFES - (unitless)
        :type CFES: float

        :return: Array of size [nz,nx,ny] with relative distribution of soil evaporation for Kzi=1

        """
        #""" distribution verticale des fractions de l'evaporation en fonction de ZESX et CFES pour Kzi=1 p- integration sur dz=1cm - Eq. 7.6 p 129"""
        
        #calcul des ponderations par cm
        nbstratesZ = len(self.dxyz[2])
        v = [self.dxyz[2][0]]
        for i in range(1, nbstratesZ):
            v.append(v[-1]+self.dxyz[2][i])

        nbcm = int(v[-1]*100)
        calc = []
        for i in range(nbcm):
            if float(i)/100. <= ZESX:
                calc.append((ZESX-float(i)/100.)**CFES)
            else:
                calc.append(0.)

        calc = array(calc)/sum(calc)
        array(list(range(nbcm)))

        #repartition dans les n strates du sol selon dxyz
        dics = {}
        for i in range(nbstratesZ):
            dics[i] = array(list(range(nbcm)))/100./v[i] >= 1.
            dics[i] = dics[i]*1

        strate = array([0]*len(calc))
        for k in list(dics.keys()):
            strate += dics[k]

        vsum = [0.]*nbstratesZ
        for i in range(len(calc)):
            vsum[strate[i]] += calc[i]

        #creation de la matrice
        m = deepcopy(self.m_1)*1.
        for z in range(nbstratesZ):
            m[z,:,:] = m[z,:,:]*vsum[z]

        return m

    def corr_vol_voxel_evap(self, ZESX=0.3):
        """ Initialise 'm_corr_vol_evap' attribute to account for the transition voxel

        :param ZESX: STICS parameter ZESX - (unit: m)
        :type ZESX: float

        :return: Array of size [nz,nx,ny] storing voxel volume fraction available to soil evaporation

        """
        #""" correction pour tenir compe du fait que voxel de transition = toute l'eau n'est pas dispo (sous ZESX) -> a stocker """
        
        v = [self.dxyz[2][0]]
        for i in range(1, len(self.dxyz[2])):
            v.append(v[-1]+self.dxyz[2][i])

        epsilon = 0.000001 #pour gerer pb d'arrondi au limite des voxels
        test = array(v)/ZESX -1.  <= epsilon
        test = test *1.
        restes = (array(v)*test)%ZESX
        test = test.tolist()
        id_transition = test.index(0.)
        test = array(test)*1.
        test[id_transition] = (self.dxyz[2][id_transition] - restes[id_transition]) / self.dxyz[2][id_transition]
        m_cor_vol = deepcopy(self.m_1)*1.
        for z in range(m_cor_vol.shape[0]): m_cor_vol[z,:,:] = m_cor_vol[z,:,:]*test[z]

        comp = (1.-m_cor_vol)
        vide = comp <1.
        comp = comp * vide*1.

        return comp #m_cor_vol #renvoie le complement finalement, car utilise dans les calculs

    def distrib_evapSTICS(self, map_Evap0):
        """ Compute distribution of evaporative demande per voxel

        :param map_Evap0: array [x, y] of evapavorative demand at soil surface (unit: mm)
        :type map_Evap0: nd.array

        :return: map_Evap, demandeEvap - representing corrected [x,y] map of evaporation and [z,x,y] array of evaporative demand, respectively - (unit: mm)

        """
        #""" map_Evap0 = liste des Evap par colone de voxel (matrice x, y) """
        
        # remplace sans bug check_soil_evap et distrib_evapNC
        coeffs = self.Kzi()*self.m_frac_evapZ
        coeffs1 = deepcopy(self.m_1)*1.
        for z in range(coeffs1.shape[0]): coeffs1[z,:,:] = sum(coeffs, axis=0) #etend en 3D la somme

        coeffsOK = coeffs/(coeffs1 + 0.0000000001) #coeffs equation 7.6, livre STICS

        demande = deepcopy(self.m_1)*1.
        for z in range(demande.shape[0]): demande[z,:,:] = map_Evap0 #etend en 3D la demande evapo

        demandeEvap = demande * coeffsOK
        dispo =  self.tsw_t - self.m_QH20min - self.m_corr_vol_evap*(self.m_QH20fc - self.m_QH20min)
        #dispo = self.tsw_t - ((self.m_QH20fc - self.m_QH20min)*(1-self.m_frac_evapZ) + self.m_QH20min)
        #dispo = self.tsw_t - ((self.m_QH20fc - self.m_QH20min) * (self.m_frac_evapZ) + self.m_QH20min)
        #dispo = self.tsw_t - ((self.m_QH20fc - self.m_QH20min) * 0.5 + self.m_QH20min) #forcage - marche bien sur AssosolNu

        diff = dispo - demandeEvap
        negs = diff<=0.
        demandeEvapOK = demandeEvap + negs*diff
        map_EvapOK = sum(demandeEvapOK,axis=0)

        return map_EvapOK, demandeEvapOK
        #S.tsw_t - demandeEvapOK


    #####  water infiltration #####

    def infil_layer(self, tsw_temp, in_ , idz, opt_infil=1):#in_=map_PI,
        """ Compute water infiltration for soil layer idz

        :param tsw_temp: Temporary [nz,nx,ny] array storing voxel total amount of soil water (unit: mm)
        :type tsw_temp: nd.array
        :param in_: Input infiltration water array coming from (idz-1) layer (unit: mm)
        :type in_: nd.array
        :param idz: z value of voxel/layer ID in grid [nz, nx, ny]
        :type idz: int
        :param opt_infil: Option indicating if infiltration is vertical (=1) or distributed between first order voxels of the layer below (=2), default=1
        :type opt_infil: int

        :return: new_tsw, out - indicating new temporary tsw array after infiltration and output water towards (idz+1) layer, respectively (unit: mm)

        """
        tsw_t = tsw_temp
        m_QH20max = self.m_QH20fc
        dxyz = self.dxyz
        
        new = tsw_t[idz]+in_
        out_ = deepcopy(in_)
        out_.fill(0.)
        for x in range(len(tsw_t[idz])):
            for y in range(len(tsw_t[idz][x])):
                if new[x][y] > m_QH20max[idz][x][y]:# si superieur a fc infiltration
                    q_out = new[x][y]-m_QH20max[idz][x][y]
                    ls_v = self.ls_1storder_vox(x,y,idz, opt_infil)#distribution entre les 1st order ; mettre opt=1 si veut forcer verticalement / 2 si
                    if len(ls_v)>1:
                        ponder = [0.0416666, 0.0416666, 0.0416666, 0.0416666, 2/3., 0.0416666, 0.0416666, 0.0416666, 0.0416666]# 2/3 en dessous 1/3 au premier ordre
                    else:
                        ponder = [1.]

                    for i in range(len(ls_v)):#distribution dans les voxels du dessous
                        nx, ny = ls_v[i][0], ls_v[i][1]
                        out_[nx][ny] = out_[nx][ny]+ponder[i]*q_out #bug: out_[nx][ny] = ponder[i]*q_out Resolu! 
                     
                    new[x][y] = m_QH20max[idz][x][y]# ce qui reste dans le voxel = fc

        return new, out_
        #infil_layer(asw_t, m_QH20max, array([[5.]]), dxyz, 0)


    def distrib_PI(self, tsw_temp, map_PI, opt_infil=1):
        """ Distribute water infiltration into the soil grid from rain and/or irrigation distribution map at the soil surface

        :param tsw_temp: Temporary [nz,nx,ny] array storing voxel total amount of soil water (unit: mm)
        :type tsw_temp: nd.array
        :param map_PI: Input infiltration water map coming from rain and/or irrigation in the top soil layer (unit: mm)
        :type map_PI: nd.array
        :param opt_infil: Option indicating if infiltration is vertical (=1) or distributed between first order voxels of the layer below (=2), default=1
        :type opt_infil: int

        :return: tsw_t, out - indicating new tsw_t array after infiltration and output water drainage out of the last soil layer, respectively (unit: mm)

        """
        #""" distribution des precipitations+irrigation dans le sol 
        #ecoulement vertical uniquement 
        #tout ce qui est au dessus de la capacite au champ est transfere en 1 jour - drainage = sortie de la couche la plus profonde """
        
        tsw_t = tsw_temp
        m_QH20max = self.m_QH20fc
        dxyz = self.dxyz
        
        in_ = map_PI
        ls_out = []
        for z in range(len(tsw_t)):
            new, out_ = self.infil_layer(tsw_t, in_, z, opt_infil)
            ls_out.append(out_)
            in_ = out_
            tsw_t[z] = new

        return tsw_t, ls_out #tsw_t, out_
        # prevoir ditribution normale en fonction de distance au centre du voxel pour "homogeneiser" en profondeur et eviter chemin preferentiel vers drainage et limiter sensibilite a la taille des voxels
        # peut aussi ajouter stockage jusqua saturation puis relargage selon vitesse d'infiltration (STICS p224)
        #map_PI = array([[5.]])
        #distrib_PI(asw_t, m_QH20max, map_PI)


    #####  plant water uptake #####

    def water_uptakeVox(self, ls_masked_asw, ls_roots_eff, ls_transp, idx, idy, idz, opt_Wvox=1):
        """ Compute water uptake relative distribution for a list of individual plants in voxel [idz,idx,idy] -
        (STICS book, p 140)

        :param ls_masked_asw: Array of size [nbplt,z,x,y] storing available soil water (i.e. if root length>0 and teta>wp; else=0) to each individual plant (unit: mm)
        :type ls_masked_asw:  nd.array
        :param ls_roots_eff: List of [z,x,y] arrays for efficient root length distribution per individual plant (unit: m)
        :type ls_roots_eff: list
        :param ls_transp: List of individual plant daily transpiration (unit: mm)
        :type ls_transp: list
        :param idx: x value of voxel ID in grid [nz, nx, ny]
        :type idx: int
        :param idy: y value of voxel ID in grid [nz, nx, ny]
        :type idy: int
        :param idz: z value of voxel ID in grid [nz, nx, ny]
        :type idz: int
        :param opt_Wvox: option for computing water uptake partitioning between plants and voxels (default: 2=weighted asw and root length; 1 = STICS book eq 7.24; 2: root length weighted by available water)
        :type opt_Wvox: int

        :return: ls_uptakeVox, list of individual plant water uptake in voxel [idz,idx,idy] (unit: mm)

        """


        ls_uptake = []

        if opt_Wvox==1:  #""" STICS manual (p 140 7.24) adapte pour liste de root systems"""
            for p in range(len(ls_roots_eff)):
                if ls_masked_asw[p][idz][idx][idy]>0.:#y a de l'eau dispo
                    frac_root = ls_roots_eff[p][idz][idx][idy]/sum3(ls_roots_eff[p]) # fraction de racine de la plante p dans le voxel
                    frac_asw = ls_masked_asw[p][idz][idx][idy]/sum3(ls_masked_asw[p]) # fraction d'asw de la plante p dans le voxel
                    uptake_vox = (ls_transp[p] * (frac_root+frac_asw))/2.
                else:
                    uptake_vox = 0.

                ls_uptake.append(uptake_vox)
        elif opt_Wvox==2: # new option with weighted relative length (p231 new book 2023)
            for p in range(len(ls_roots_eff)):
                if ls_masked_asw[p][idz][idx][idy]>0.:#y a de l'eau dispo
                    frac_weigh = (ls_roots_eff[p][idz][idx][idy] * ls_masked_asw[p][idz][idx][idy]) / sum3(ls_roots_eff[p] * ls_masked_asw[p])
                    uptake_vox = ls_transp[p] * frac_weigh
                    # pas fraction de longueu de la plante en relatif des autres plantes dans le voxel
                    # gere distribution de la demande de plante p entre voxel.
                    # partitionnement de la demande se fait sur la lumiere -> ls_trans[p]
                    # et apres si pas assez??? -> gere avec condition asw > 0
                else:
                    uptake_vox = 0.

                ls_uptake.append(uptake_vox)
        else:
            print("unknown opt_Wvox option")

        return ls_uptake
        #tester si prend plus... ajuster?
        #water_uptakeVox(ls_masked_asw, ls_roots_eff, ls_Et0, 0, 0, 0)


    def distrib_water_uptakeNC(self, ls_masked_asw, ls_roots_eff, ls_transp):
        """ Compute daily water partitioning between individual plants and the spatial distribution of plant water uptake (unit: mm)

        :param ls_masked_asw: Array of size [nbplt,z,x,y] storing available soil water (i.e. if root length>0 and teta>wp; else=0) to each individual plant (unit: mm)
        :type ls_masked_asw:  nd.array
        :param ls_roots_eff: List of [z,x,y] arrays for efficient root length distribution per individual plant (unit: m)
        :type ls_roots_eff: list
        :param ls_transp: List of individual plant daily transpiration (unit: mm)
        :type ls_transp: list

        :return: uptake, ls_uptake - corresponding to a [z,x,y] array for total water uptake by plants and a [nbplt,z,x,y] array storing individual plant daily water uptakes, respectively

        """
        asw_t = deepcopy(self.asw_t)
        upt = self.m_1*0.

        ls_upt = []
        for i in range(len(ls_roots_eff)):
            ls_upt.append(self.m_1*0.)

        for z in range(len(asw_t)):
            for x in range(len(asw_t[z])):
                for y in range(len(asw_t[z][x])):
                    upvx = self.water_uptakeVox(ls_masked_asw, ls_roots_eff, ls_transp, x, y, z)
                    upt[z][x][y] = sum(upvx) #gere negatif avec masque en amont
                    for plt in range(len(ls_upt)):
                        ls_upt[plt][z][x][y] = upvx[plt]

        return upt, ls_upt
        #renvoyer info par plante?ici distribue slmt entre voxel -> oui!! besoin pour les flux d'N

        #upt = distrib_water_uptakeNC(asw_t, ls_masked_asw, ls_roots_eff, ls_Et0)
        #sum(distrib_water_uptakeNC(asw_t, ls_masked_asw, ls_roots_eff, ls_Et0)), sum(ls_Et0)


    #####  Daily loop #####

    def stepWBmc(self, Et0, ls_roots, ls_epsi, Rain, Irrig, previous_state, ZESX=0.3, leafAlbedo=0.15, U=5., b=0.63, FTSWThreshold=0.4, treshEffRoots=0.5, opt=1):
        """ Compute daily step Soil water balance -
        Apply a sequence: Plant transipration and soil evaporation -> Rain and Irrigation application -> Water infiltration and Drainage

        :param Et0: Daily potential evapotranspiration (unit: mm)
        :type Et0: float
        :param ls_roots: List of [z,x,y] arrays for root length distribution per individual plant (unit: m)
        :type ls_roots: list
        :param ls_epsi: List of individual plant daily fraction of incoming global radiation interception (0-1 fraction)
        :type ls_epsi: list
        :param Rain: Daily rain (unit: mm)
        :type Rain: float
        :param Irrig: Daily irrigation (unit: mm)
        :type Irrig: float
        :param previous_state: List of memory variables storing cumulative previous day evaporation used for computing daily soil evaporation
        :type previous_state: list
        :param ZESX: STICS parameter ZESX - (unit: m), default to 0.3
        :type ZESX: float
        :param leafAlbedo: STICS parameter leafAlbedo - (0-1 fraction), default to 0.15
        :type leafAlbedo: float
        :param U: Total topsoil water reservoir accessible to evaporation (unit: mm), default to 5.
        :type U: float
        :param b: Empirical coefficient for soil evaporation (unitless), default to 0.63
        :type b: float
        :param FTSWThreshold: STICS FTSW threshold for reducing plant transpiration (unit: 0-1 fraction), default to 0.4
        :type FTSWThreshold: float
        :param treshEffRoots: STICS maximal threshold for root length density efficient in water uptake (unit: cm.cm-3), default to 0.5
        :type treshEffRoots: float
        :param opt_infil: Option indicating if infiltration is vertical (=1) or distributed between first order voxels of the layer below (=2), default=1
        :type opt_infil: int

        :return:
            * 'ls_transp': List of daily transpiration by individual plants (unit: mm)
            * 'evapo_tot': Daily water evaporation (unit: mm)
            * 'ls_drainage': Daily water drainage out of the last soil layer (unit: mm)
            * 'state': List of current day memory states for computing soil evaporation
            * 'ls_m_transpi': A [nbplt,z,x,y] nd.array storing individual plant daily water uptake distributions (unit: mm)
            * 'm_evap': A [z,x,y] nd.array storing daily water evaporation distribution (unit: mm)
            * 'ls_ftsw': List of daily ftsw by individual plants (0-1 fraction)

        """

        #""" calcul du bilan hydrique journalier """
        # SORTIES:
        # evapo_tot : quantite d'eau evaporee totale au temps t (mm)
        # m_evap : quantite d'eau evaporee au temps t par voxel (mm)
        # m_transpi : quantite d'eau transpiree au temps t par voxel (mm) - cumul des n plantes par voxel
        # ls_drainage : quantite d'eau drainee au temps t par voxel (mm) - liste de matrice par horizon
        # ls_transp : quantite d'eau transpiree au temps t par plante (mm) - liste de cumul par plante
        # ls_ftsw : liste des fraction d'eau disponible pour la transpiration par plante (fraction 0-1) - tient compte des distributions des racines

        Precip = Rain + Irrig

        ## transpi 
        ls_roots_eff = effective_root_lengths(ls_roots, tresh = treshEffRoots)
        ls_mask = list(map(mask, ls_roots_eff))#liste de masques racine
        masked_aswt = multiply(mask(self.asw_t), self.asw_t) #pour sol, met a zero si teta sous wp
        ls_masked_asw = list(map(multiply, ls_mask, [masked_aswt]*len(ls_mask)))#pour chaque root system, asw(mm) si racine et teta> wp; zero sinon
        ls_masked_tswt = list(map(multiply, ls_mask, [self.m_QH20max]*len(ls_mask)))
        ls_ftsw = list(map(divide, list(map(float, list(map(sum3, ls_masked_asw)))), list(map(float, list(map(sum3, ls_masked_tswt))))))
        paramp = [{"leafAlbedo":leafAlbedo, "WaterTreshGs":FTSWThreshold}] # a faire passer en amont
        ls_transp = Transpi_NC(Et0, ls_epsi, ls_ftsw, paramp)
        #m_transpi = distrib_water_uptakeNC(self.asw_t, ls_masked_asw, ls_roots_eff, ls_transp)
        m_transpi, ls_m_transpi = self.distrib_water_uptakeNC(ls_masked_asw, ls_roots_eff, ls_transp)
        ls_transp = list(map(sum, ls_m_transpi)) #correct ls_trans with actual transpiration

        ## evaporation
        #evapo_tot, state = soil_EV_1C(Et0, Precip, sum(ls_epsi), previous_state , leafAlbedo, U, b)
        #evapo_tot, state = soil_EV_STICS(Et0 , Precip , sum(ls_epsi),previous_state, leafAlbedo, U, b)
        #map_Evap0 = self.m_1[0] * evapo_tot  / float(len(self.dxyz[0]) * len(self.dxyz[1]))  # map_Evap0 -> devrait etre une entree? (interaction avec bilan radiatif?)

        evapo_tot, state = soil_EV_STICS(Et0/self.soilSurface, Precip/self.soilSurface, sum(ls_epsi), previous_state , leafAlbedo, U, b)
        map_Evap0 = self.m_1[0]*evapo_tot*self.soilSurface/float(len(self.dxyz[0])*len(self.dxyz[1])) #map_Evap0 -> devrait etre une entree? (interaction avec bilan radiatif?)

        ## bilan WB
        map_PI = self.m_1[0]*Precip/float(len(self.dxyz[0])*len(self.dxyz[1]))
        # 1) applique pluie et calcule D0 (drainage sans transpi pour saturer)
        ##m, D0 = distrib_PI(self.asw_t, self.m_QH20max, map_PI, self.dxyz,opt)
        m, ls_D0 = self.distrib_PI(self.tsw_t, map_PI, opt)
        D0 = ls_D0[-1]
        # 2) retire evapo et transpi
        self.tsw_t = m 
        map_EvapOK, m_evap = self.distrib_evapSTICS(map_Evap0) # evap tient compte des pluies
        evapo_tot = sum(map_EvapOK)#mean(map_EvapOK) #corrige pour eviter valeur sous mini
        m = m - m_evap - m_transpi

        # 3) remet ce qui a drainer et recalcule drainage residuel apres transpi (=transpire en priorite ce qui est tombe pdt la journee)
        ##m, D1 = distrib_PI(m, self.m_QH20max, D0, self.dxyz)
        m, ls_drainage = self.distrib_PI(m, D0)
        #D1 = ls_drainage[-1] #drainage profond

        self.tsw_t = m
        self.update_asw()
        self.update_ftsw()
        self.UpdateWbalance(Rain, Irrig, evapo_tot, ls_transp, ls_drainage) #Irrigj pas distingue en entree pour le moment!-> a faire

        return ls_transp, evapo_tot, ls_drainage, state,  ls_m_transpi, m_evap, ls_ftsw


    #def plot_soil_properties (self, vals, MaScene=Scene(), col_scale=5):#dxyz, m_soil_vox, asw_t pourraient etre remplacee par objet sol
    #    """ vals = matrice de valeur entre 0 et 1 de propriete de sol a visualiser / e.g S.ftsw_t """
    #    bx = Box(Vector3(1.,1.,1.))
    #    for z in range(len(self.dxyz[2])):
    #        for x in range(len(self.dxyz[0])):
    #            for y in range(len(self.dxyz[1])):
    #                dims = [self.dxyz[0][x]/2., self.dxyz[1][y]/2., self.dxyz[2][z]/2.]
    #                p_ini = self.m_soil_vox[z][x][y]
    #                col = couleur (col_scale, max(0.,vals[z][x][y]))
    #                b = transformation(bx, dims[0], dims[1], dims[2], 0,0,0, p_ini[0], p_ini[1], -p_ini[2]-self.dxyz[2][z]/2.)
    #                MaScene.add(Shape(b, Material(Color3(col[0],col[1],col[2]), transparency=0.65)))
    #
    #    return MaScene
    ## mis dans soil_modulevisu1.py


    def OpenWbalance(self):
        """ Initialise bilanW attribute, a dictionary storing water balance information

            * Keys for simulation inputs:
                - 'intialWC': Initial soil water (unit: mm)
                - 'Irrigtot': Total simulation irrigation (unit: mm)
                - 'PPtot': Total simulation precipitations (unit: mm)
            * Keys for simulation outputs:
                - 'FinalWC': Final soil water (unit: mm)
                - 'EVtot': Total simulation soil evaporation (unit: mm)
                - 'Tranptot': Total simulation plant transpiration (unit: mm)
                - 'Drainagetot': Total simulation water drainage (unit: mm)
            * Keys for simulation totals:
                - 'InputWtot': Total simulation water inputs (unit: mm)
                - 'OutputWtot': Total simulation water outputs (unit: mm)
            * Keys for daily inputs/outputs :
                - 'cumPP': List of daily precipitations (unit: mm)
                - 'cumIrrig': List of daily irrigation (unit: mm)
                - 'cumD': List of daily soil water drainage (unit: mm)
                - 'cumEV': List of daily soil evaporation (unit: mm)
                - 'cumTransp': List of daily total plant transpiration (unit: mm)


        """
        #""" Dictionnary for soil Water balance (mm)
        #Keys for Daily outputs: 'cumPP', 'cumIrrig', 'cumD', 'cumEV', 'cumTransp'
        #Keys for Total Input: 'intialWC', 'Irrigtot','PPtot'
        #Keys for Total outputs: 'FinalWC', 'EVtot', 'Tranptot','Drainagetot'
        #Keys for totals: 'InputWtot', 'OutputWtot'
        #"""
        
        self.bilanW = {}
        self.bilanW['intialWC'] = sum3(self.tsw_t) / self.soilSurface
        self.bilanW['cumEV'], self.bilanW['cumPP'], self.bilanW['cumIrrig'], self.bilanW['cumD'], self.bilanW['cumTransp'] = [],[],[], [], []
        self.bilanW['TSWt'] = []
        #cumET0 = []
        
    def UpdateWbalance(self, PPj, Irrigj, evapo_tot, ls_transp, Drainage):
        """ Daily update of bilanW attribute

        """
        surfsolref = self.soilSurface
        self.bilanW['cumPP'].append(PPj / surfsolref)
        self.bilanW['cumIrrig'].append(Irrigj / surfsolref)
        self.bilanW['cumTransp'].append(sum(ls_transp) / surfsolref)
        self.bilanW['cumEV'].append(evapo_tot/ surfsolref)
        self.bilanW['cumD'].append(sum(Drainage[-1]) / surfsolref)
        self.bilanW['TSWt'].append(sum3(self.tsw_t) / self.soilSurface)


    def CloseWbalance(self, print_=1):
        """ Finalise calculation of whole simulation variables in bilanW attribute
    
        """
        #input
        self.bilanW['PPtot'] = sum(self.bilanW['cumPP'])
        self.bilanW['Irrigtot'] = sum(self.bilanW['cumIrrig'])
        self.bilanW['InputWtot'] = self.bilanW['intialWC'] + self.bilanW['PPtot'] + self.bilanW['Irrigtot']
        #output
        self.bilanW['FinalWC'] = sum3(self.tsw_t) / self.soilSurface
        self.bilanW['EVtot'] = sum(self.bilanW['cumEV'])
        self.bilanW['Tranptot'] = sum(self.bilanW['cumTransp'])
        self.bilanW['Drainagetot'] = sum(self.bilanW['cumD'])
        self.bilanW['OutputWtot'] = self.bilanW['FinalWC'] + self.bilanW['EVtot'] + self.bilanW['Tranptot'] + self.bilanW['Drainagetot']

        #print self.bilanW['OutputWtot'], self.bilanW['InputWtot']
        if print_==1:
            self.PrintWbalance()
        #pourrait le diriger vers un fichier de sortie texte?

    def PrintWbalance(self):
        """ Print a summary table of bilanW attribute
    
        """
        bilanW = self.bilanW
        print ("")
        print ("Water Balance Input (mm)\t\t\t\t\t Water Balance Output (mm)")
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Initial Soil Water:\t {0:8.1f}\t\t\t\t Final Soil Water:\t {1:8.1f}".format(bilanW['intialWC'], bilanW['FinalWC'])))
        print(("Precipitation:\t\t {0:8.1f}\t\t\t\t Transpiration:\t\t {1:8.1f}".format(bilanW['PPtot'], bilanW['Tranptot'])))
        print(("Irrigation:\t\t\t {0:8.1f}\t\t\t\t Evaporation:\t\t {1:8.1f}".format(bilanW['Irrigtot'],bilanW['EVtot'])))
        print(("                            \t\t\t\t Deep infiltration:\t {0:8.1f}".format(bilanW['Drainagetot'])))
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Total:\t\t\t\t {0:8.1f}\t\t\t\t Total:\t\t\t\t {1:8.1f}".format(bilanW['InputWtot'], bilanW['OutputWtot'])))
        print ("")




def default_par_sol():
    """ Creates a default parameter dictionnary 'par_sol' for defining the different soil types found in Soil object
    and used for computing water balance - Example with 4 possible soil layers / types

        Keys of the dictionnary are  parameters for a given soil layer:
            * 'soil number': soil type layer ID
            * 'soil type': soil type layer name
            * 'teta_sat': soil volumic humidity at water saturation (0-1 fraction)
            * 'teta_fc': soil volumic humidity at field capacity (0-1 fraction)
            * 'teta_wp': soil volumic humidity at wilting point (0-1 fraction)
            * 'teta_ad': soil volumic humidity at minimal water content (0-1 fraction)
            * 'DA': soil bulk density (unit: g.cm-3)

    :return: Default 'par_sol' parameter dictionnary

    .. code-block:: python

        par_sol = {
                    '1': {'soil number': '1', 'soil type': 'H5cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.479,
                    'teta_wp': 0.380, 'teta_ad': 0.340, 'DA': 1.36},

                    '2': {'soil number': '2', 'soil type': 'H10cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.347,
                    'teta_wp': 0.275, 'teta_ad': 0.248, 'DA': 1.36}
                    }


    """
    # parameters 'WCST' and 'gamma_theo' are not used / necessary in current model version

    par_sol = {}
    # add soil type 1
    par_sol['1'] = {'soil number': '1', 'soil type': 'H5cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.479,
                    'teta_wp': 0.380, 'teta_ad': 0.340, 'DA': 1.36, 'WCST': '0.35', 'gamma_theo': '0.07'}
    # add soil type 2
    par_sol['2'] = {'soil number': '2', 'soil type': 'H10cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.347,
                    'teta_wp': 0.275, 'teta_ad': 0.248, 'DA': 1.36, 'WCST': '0.35', 'gamma_theo': '0.07'}
    # add soil type 3
    par_sol['3'] = {'soil number': '3', 'soil type': 'H25cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.377,
                    'teta_wp': 0.301, 'teta_ad': 0.248, 'DA': 1.36, 'WCST': '0.35', 'gamma_theo': '0.07'}
    # add soil type 4
    par_sol['4'] = {'soil number': '4', 'soil type': 'H90cm - ASCHYD11', 'teta_sat': 0.55, 'teta_fc': 0.426,
                    'teta_wp': 0.312, 'teta_ad': 0.248, 'DA': 1.36, 'WCST': '0.35', 'gamma_theo': '0.07'}

    return par_sol



#def plot_soilWC (dxyz, m_soil_vox, ftsw_t, MaScene=Scene(), col_scale=5):#dxyz, m_soil_vox, asw_t pourraient etre remplacee par objet sol
#    """ """
#    for z in range(len(dxyz[2])):
#        for x in range(len(dxyz[0])):
#            for y in range(len(dxyz[1])):
#                dims = [dxyz[0][x]/2., dxyz[1][y]/2., dxyz[2][z]/2.]
#                p_ini = m_soil_vox[z][x][y] 
#                col = couleur (col_scale, max(0.,ftsw_t[z][x][y]))
#                b = transformation(bx, dims[0], dims[1], dims[2], 0,0,0, p_ini[0], p_ini[1], -p_ini[2]-dxyz[2][z]/2.)
#                #b.setName(luz_label(1, 0, 0, 0, 0))#ajout d'un label canestra a ajuster
#                MaScene.add(Shape(b, Material(Color3(col[0],col[1],col[2]), transparency=0.65)))
#
#    return MaScene
## a passer comme methode classe







##  puis faire bilan
#def TSW_NC(previousTSW, map_PI, m_frac_evap, m_frac_transpi, m_QH20max, dxyz, opt=1):
#    """ """
#    # 1) applique pluie et calcule D0 (drainage sans transpi pour saturer)
#    m, D0 = distrib_PI(previousTSW, m_QH20max, map_PI, dxyz,opt)
#    # 2) retire evapo et transpi
#    m = m - m_frac_evap - m_frac_transpi
#    # 3) remet ce qui a drainer et recalcule drainage residuel apres transpi (=transpire en priorite ce qui est tombe pdt la journee)
#    m, D1 = distrib_PI(m, m_QH20max, D0, dxyz)
#    return m, D1



