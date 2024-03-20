'''

 3DS-soil-model, a 3D Soil model adapted from STICS : soil Nitrogen balance module
 ******************************
 Authors: G. Louarn 
 

 

'''



#from scipy import *
#from rpy import r
#import sys
#path_ = r'H:\devel\grassland\grassland\luzerne'
#path2_ = r'H:\devel\grassland\grassland\L-gume'
#sys.path.insert(0, path_)
#sys.path.insert(0, path2_)
#from soil_module import * #soil_module5
from soil3ds.soil_moduleW import * #soil3ds installe comme module


class SoilN(Soil):
    """ Main class for the soilN object of the 'soil3ds' model

    The class :class:`soil3ds.soil_moduleN.SoilN` includes descriptions of soil properties and methods to compute nitrogen balance adpated for a 3D grid from STICS model.
    The class inherits from 'Soil' class which includes descriptions of soil physical properties and methods to compute water balance

    :param par_sol: A dictionnary defining water humidity thresholds for each soil type necessary to build the Soil object
    :type par_sol: dict
    :param parSN: A dictionnary defining general soil parameters derived from STICS
    :type parSN: dict
    :param soil_number: List of nz soil type IDs along soil depth
    :type soil_number: list
    :param dxyz: List of three lists indicating voxel discretization along [x,y,z] axes (unit: m)
    :type dxyz: list
    :param vDA: List of nz soil bulk density along soil depth (unit: g.cm-3)
    :type vDA: list
    :param vCN: List of nz soil C:N ratio along soil depth (unitless)
    :type vCN: list
    :param vMO: List of nz soil soil organic matter concentration along soil depth (unit: g C.kg soil-1)
    :type vMO: list
    :param vARGIs: List of nz soil soil clay content along soil depth (unit: %)
    :type vARGIs: list
    :param vNO3: List of nz soil soil nitrate content along soil depth to initialize 'm_NO3' (unit: kg N-NO3.ha-1)
    :type vNO3: list
    :param vNH4: List of nz soil soil ammonium content along soil depth to initialize 'm_NH4' (unit: kg N-NH4.ha-1)
    :type vNH4: list
    :param vCALCs: List of nz soil calcacerous content along soil depth (unit: %)
    :type vCALCs: list
    :param Tsol: Daily average topsoil temperature (unit: degree Celsius)
    :type Tsol: float
    :param obstarac: List of nz voxel fractions affected by 'obstarac' along soil depth, default to None (no root obstacles)
    :type obstarac: list
    :param pattern8: list of two [x,y] points defining the soil limits (unit: cm), default to [[0,0], [100,100]]
    :type pattern8: list


    SoilN object Attributes:
        * :CALCs (float): Calcareous content in topsoil layer - value for non calcareous soils (p220 STICS-book) - (unit: %)
        * :pHeau (float): Soil water pH
        * :ARG (float): Clay content in topsoil layer

        * :stateEV (list): List of memory variables storing cumulative previous day evaporation used for computing daily soil evaporation
        * :Uval (float): Total topsoil water reservoir accessible to evaporation (unit: mm)
        * :b (float): Empirical coefficient for soil evaporation (unitless)

        * :Corg (nd.array): Array of size [nz,nx,ny] storing voxel soil organic Carbon content (unit: kg C per voxel)
        * :Norg (nd.array): Array of size [nz,nx,ny] storing voxel soil organic Nitrogen content (unit: kg N per voxel)
        * :InertCorg (nd.array): Array of size [nz,nx,ny] storing voxel soil inert organic Carbon content (unit: kg C per voxel)
        * :InertNorg (nd.array): Array of size [nz,nx,ny] storing voxel soil inert organic Nitrogen content (unit: kg N per voxel)
        * :K2HUM (nd.array): Array of size [nz,nx,ny] storing voxel Potential rate of SOM mineralisation (unit: ...)
        * :m_MO (nd.array): Array of size [nz,nx,ny] storing voxel soil organic matter (SOM) concentration (unit: g C.kg soil-1)
        * :m_CNHUM (nd.array): Array of size [nz,nx,ny] storing voxel soil C:N ratio (unitless)
        * :m_NH4 (nd.array): Array of size [nz,nx,ny] storing voxel content in mineral N-NH4  (unit: kg N-NH4 per voxel)
        * :m_NO3 (nd.array): Array of size [nz,nx,ny] storing voxel content in mineral N-N03  (unit: kg N-NO3 per voxel)
        * :m_Tsol (nd.array): Array of size [nz,nx,ny] storing voxel daily average temperature (unit: degree Celsius)

        * :ls_CRES (list): List storing [nz,nx,ny] 'm_CRES' arrays by type of residue (unit: kg C per voxel)
        * :ls_CBio (list): List storing [nz,nx,ny] 'm_CBio' arrays by type of residue (unit: kg C per voxel)
        * :parResi (dict): A dictionnary storing the mineralisation parameters of all soil organic residues

        * :CO2respSoil (float): Cumulative amount of C-CO2 respired from all soil residues (unit: kg C)
        * :N2ONitrif (float): Cumulative amount of N-N2O emitted from nitification reactions (unit: kg N)
        * :N2ODenitrif (float): Cumulative amount of N-N2O emitted from denitification reactions (unit: kg N)
        * :lixiNO3 (float): Cumulative amount of N-NO3 tranferred below bottom soil layer (unit: kg N)

        * :bilanC (dict): Dictionnary storing soil Carbon balance daily outputs
        * :bilanN (dict): Dictionnary storing soil Nitrogen balance daily outputs


    Main methods:

        * ``__init__``: initializes and builds static object for the rest of simulation
        * :meth:`ConcNO3`: Calculation of the nitrate concentration (unit: kg N.mm-1)
        * :meth:`moleN`: Calculation of the mineral N molar content (unit: micromole N per voxel)
        * :meth:`ConcN`: Calculation of the molar concentration of mineral nitrogen (unit: micromole N.L-1)
        * :meth:`ConcN_roots`: Calculation of the average molar concentration of mineral nitrogen in each plant root zone (unit: micromole N.L-1)
        * :meth:`ConcN_old`: Calculation of the molar concentration of mineral nitrogen - former calulation to account for 1 Ha
        * :meth:`init_memory_EV`: Inititalise attributes corresponding to memory variables and parameters to compute soil evaporation
        * :meth:`update_memory_EV`: Update 'stateEV' attribute with new list of memory variables
        * :meth:`updateTsol`: Update daily soil temperature array 'm_Tsol'

    Soil mineralisation methods:

        * :meth:`mask_PROFUM`: Compute a [nz,nx,ny] mask array to apply 'PROFHUMs' mineralisation parameter
        * :meth:`Pot_rate_SOMMin`: Compute Potential rate of SOM mineralisation (unit: kg N.day-1)
        * :meth:`SOMMin_RespT`: Compute response of mineralisation rate to soil temperature (0-1 fraction)
        * :meth:`SOMMin_RespHum`: Compute response of mineralisation rate to soil relative humidity 'HRv'
        * :meth:`Act_rate_SOMMin`: Compute Actual rate of SOM mineralisation (unit: kg N.day-1)
        * :meth:`stepNB`: Compute daily step for SOM mineralisation and update SoilN object

    Residue mineralisation methods:

        * :meth:`init_residues`: Initialise 'parResi' attribute storing the properties of soil residues
        * :meth:`addResPAR`: Add a new set of residue parameters to input dictionary of type 'parResi'
        * :meth:`VdistribResidues`:  Distribute a new residue in a [nz,nx,ny] 'm_CRES' array assuming horizontal homogeneity (unit: kg C per Voxel)
        * :meth:`addResMat`: Update 'ls_CRES', 'bilanC' and 'bilanN' attributes when adding a new residue
        * :meth:`Pot_rate_ResidueMin`: Compute the rate of residue mineralisation consirering soil water and temperature effects (unit: kg C.day-1)
        * :meth:`Pot_Ndemand_microbialBio`: Compute total microbial N demand for all residues (unit: kg de N per voxel)
        * :meth:`FN_factor`: Compute 'FN' reduction factor related to Nmin availability to support residue mineralisation (0-1 fraction)
        * :meth:`FBIO_factor`: Compute 'FBIO' reduction factor (0-1 fraction)
        * :meth:`ls_NRES`: Compute a list of [nz,nx,ny] arrays storing residue N distribution per type of residue (unit: kg N per voxel)
        * :meth:`ls_NBio`: Compute a list of [nz,nx,ny] arrays storing microbial N distribution per type of residue (unit: kg N per voxel)
        * :meth:`mixResMat`: Mix a new input residue [nz,nx,ny] array with an existing residue type in 'parResi' and 'ls_CRES' attributes
        * :meth:`stepResidueMin`: Compute daily Carbon and Nitrogen fluxes associated with mineralisation of all residues
        * :meth:`stepMicrobioMin`: Compute daily Carbon and Nitrogen fluxes associated with microbial biomass of all residues

    Soil nitrification methods:

        * :meth:`Nitrif_RespHum`: Compute response of nitrification rate to soil Humidity 'HRv' (0-1 fraction)
        * :meth:`Nitrif_RespPH`: Compute response of nitrification rate to soil pH (0-1 fraction)
        * :meth:`Nitrif_RespT`: Compute response of nitrification rate to soil temperature (0-1 fraction)
        * :meth:`stepNitrif`: Compute daily N fluxes associated to NH4+ nitrification

    Nitrates infiltration methods:

        * :meth:`infil_layerNO3`: Compute nitrate infiltration for soil layer idz
        * :meth:`distrib_NO3`: Distribute nitrate infiltration into the soil grid from mineralisation and/or fertilizer distribution map at the soil surface
        * :meth:`stepNINFILT`: Compute daily N fluxes associated to nitrate infiltration and fertilizer application

    Plant N uptake methods:

        * :meth:`stepNuptakePlt`: Compute daily N fluxes associated to plant mineral N uptake

    Nitrogen and Carbon balance :

        * :meth:`OpenCbalance`: Initialise bilanC attribute, a dictionary storing carbon balance information
        * :meth:`CloseCbalance`: Finalise calculation of whole simulation variables in bilanC attribute
        * :meth:`PrintCbalance`: Print a summary table of bilanC attribute
        * :meth:`OpenNbalance`: Initialise bilanN attribute, a dictionary storing nitrogen balance information
        * :meth:`CloseNbalance`: Finalise calculation of whole simulation variables in bilanN attribute
        * :meth:`PrintNbalance`: Print a summary table of bilanN attribute

    """
    
    def __init__(self, par_sol, parSN, soil_number , dxyz, vDA,  vCN, vMO, vARGIs, vNO3, vNH4,  vCALCs, Tsol, obstarac=None, pattern8=[[0,0],[100.,100.]]):
        """ Initializes and builds static SoilN object for the rest of simulation

        """
        # """
        # par_sol SN contient en plus pour l'azote:
        #     'FMIN1G'        #(day-1) (p145)
        #     'FMIN2G'        #(%ARGIS-1) para pot rate min a ARGIs (p145)
        #     'FMIN3G'        #(% CALC-1)para pot rate min a CALCss (p145)
        #     'FINERTG'       #0.65 = default fraction of N pool inactive for minearlisation (p145) (could be smaller in grassland & forest)
        #     'PROFHUMs'      # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220
        #     'HMinMg'        #Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
        #     'HoptMg'        #Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)
        #     'TRefg'         #reference temperature (degreC)
        #     'FTEMHAg'       #asymptotic value of FTH (seen as a logistic response)
        #     'FTEMHg'        #(K-1)
        #     'FTEMHB'
        #     'FNXg'          #(day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
        #     'PHMinNITg'     #pH min de nitrification (prop of Field Capacity) #value p149
        #     'PHMaxNITg'     #pH max de nitrification (prop of Field Capacity) #value p149
        #     'HMinNg'        #Humidite min de nitrification #value p149
        #     'HoptNg'        #Humidite opt de nitrification  #value p149
        #     'TNITMINg'      #Temperature min de nitrification #  (degreC)value p151
        #     'TNITOPTg'      #Temperature opt de nitrification #  (degreC)value p151
        #     'TNITMAXg'      #Temperature max de nitrification #  (degreC)value p151
        #     'RATIONITs'     #proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
        #     'DIFNg'         #N diffusion coefficient at field capacity (cm2.day-1, p 161)
        #
        # CALCs
        # pHeau
        # ARGIs
        # m_Tsol
        #
        # Corg: (kg C dans le voxel)
        # Norg: (kg N dans le voxel)
        # InertCorg: (kg C dans le voxel)
        # InertNorg: (kg N dans le voxel)
        # m_NH4: (kg N NH4 dans le voxel)
        # m_NO3: (kg N NO3 dans le voxel)
        # m_MO
        # m_CNHUM
        #
        # K2HUM
        # N2ONitrif
        # N2ODenitrif
        # lixiNO3
        #
        # bilanC et bilanN: dictionnaires contenant les variables dynamiques et les cumuls necessaire a l'etablissement des bilans C et N (kg C/N.ha-1)
        # """

        #initialisation sol et teneur en eau
        Soil.__init__(self,par_sol, soil_number, dxyz, vDA, parSN['ZESX'], parSN['CFES'], obstarac, pattern8)
        self.compute_teta_lim(par_sol)
        self.init_asw()
        # init variables memoire evaporation
        self.init_memory_EV(parSN)

        self.CALCs = vCALCs[0] #(%) value for non calcareous soils (p220 STICS) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! A mesurer/mettre a jour (rq: mesure CaO anterieure parcelle C2: 1.58 g.kg-1 dans 0-30)
        self.pHeau = parSN['pH']#pH
        self.ARG = vARGIs[0]
        #a faire: NORGs matrices de NORGs, CORGs, K2HUM...
        
        #initialisation matrice Temp et et compartiments azote sol
        self.m_MO, self.m_CNHUM, self.K2HUM, self.m_Tsol, self.m_NH4, self.m_NO3 = [], [], [], [], [], []
        for z in range(len(dxyz[2])):
            v, v1, v2, v3, v4, v5 = [], [], [], [], [], []
            for x in range(len(dxyz[0])):
                vv, vv1, vv2, vv3, vv4, vv5 = [], [], [], [], [], []
                for y in range(len(dxyz[1])):
                    surf= self.dxyz[0][x] * self.dxyz[1][y]
                    vv.append(vMO[z])
                    vv1.append(vCN[z])
                    vv2.append(self.Pot_rate_SOMMin(vCALCs[z], vARGIs[z], parSN))
                    vv3.append(Tsol) #suppose cste et egale a Tair -> a adapter
                    vv4.append(vNH4[z]/10000.*surf) #en kg d'N par voxel
                    vv5.append(vNO3[z]/10000.*surf) #en kg d'N par voxel

                v.append(vv)
                v1.append(vv1)
                v2.append(vv2)
                v3.append(vv3)
                v4.append(vv4)
                v5.append(vv5)

            self.m_MO.append(v)
            self.m_CNHUM.append(v1)
            self.K2HUM.append(v2)
            self.m_Tsol.append(v3)
            self.m_NH4.append(v4)
            self.m_NO3.append(v5)

        self.m_MO, self.m_CNHUM,  self.m_Tsol, self.m_NH4, self.m_NO3  = array(self.m_MO), array(self.m_CNHUM),  array(self.m_Tsol), array(self.m_NH4), array(self.m_NO3)
        self.K2HUM = array(self.K2HUM) 
        self.Corg = multiply(self.m_soil_vol , self.m_MO)* self.mask_PROFUM(parSN) * self.m_DA / 1.72 #en kg de C (#aplique sur Profhum)
        self.Norg = divide(self.Corg, self.m_CNHUM) #en kg de N
        self.InertCorg = self.Corg * parSN['FINERTG']
        self.InertNorg = self.Norg * parSN['FINERTG']
        self.OpenCbalance()
        self.OpenNbalance()
        
        self.N2ONitrif, self.N2ODenitrif = 0., 0. #kg N for the whole soil volume
        self.lixiNO3 = 0. #kg N
        #reprendre les Temperature avec un modele plus elabore!! ici = Tair!


    #####  Main methods : initialise/ update / examine soilN object #####
    
    def ConcNO3(self):
        """ Calculation of the nitrate concentration (unit: kg N.mm-1)
        (STICS book, Eq. 8.34, p 160)
        
        """
        #mm d'eau libre? (eau liee retiree -  pas dans les eaux de drainage)?
        #non - rq: dans devienne-baret (2000): utilise toute l'eau du sol pour calcul de concentration
        return self.m_NO3 / self.tsw_t#(S.tsw_t - S.m_QH20min + 0.00000001)

    def moleN(self):
        """ Calculation of the mineral N molar content (unit: µmole N per voxel)

        """
        # """ calculation of the molar concentration of mineral nitrogen (micromole N.L-1) by voxel - 8.36 p 161 """
        # L d'eau libre (eau liee retiree - car pas dans les eaux de drainage)
        MMA = 71.428  # 142.85 #mole d'N.kg-1 (14 g.mole-1)
        moleN_ = (self.m_NO3 + self.m_NH4) * (MMA) * 10 ** 6  # micromole d'N
        return moleN_

    def ConcN(self):
        """ Calculation of the molar concentration of mineral nitrogen (unit: micromole N.L-1)
        (STICS book, Eq. 8.36, p 161)
        
        """
        return self.moleN() / self.tsw_t  #micromole d'N.L-1


    def ConcN_old(self):
        """ Calculation of the molar concentration of mineral nitrogen - former calulation to account for 1 Ha
        
        """
        #""" calculation of the molar concentration of mineral nitrogen (micromole N.L-1) by voxel - 8.36 p 161 """
        # L d'eau libre (eau liee retiree - car pas dans les eaux de drainage)
        MMA = 142.85  # mole d'N.kg-1 (g.mole-1)
        moleN_ = (self.m_NO3 + self.m_NH4) / (MMA) * 10 ** 6  # micromole d'N

        return moleN_ / (self.tsw_t * self.m_vox_surf) * 10000  # remis pour conc sur 1ha pour coller au parametrage de sTICS + erreur voxsurf

    def ConcN_roots(self, ls_roots):
        """ Calculation of the average molar concentration of mineral nitrogen in each plant root zone (unit: micromole N.L-1)

        :param ls_roots:
        :type ls_roots: list

        :return: ls_concN:
        """
        ls_mask_r = list(map(mask, ls_roots))  # mask de racines
        ls_masked_moleN = list(map(np.multiply, ls_mask_r, [self.moleN()] * len(ls_mask_r))) #micromole N
        ls_vol_eau = list(map(np.multiply, ls_mask_r, [self.tsw_t] * len(ls_mask_r))) #L
        ls_concN = list(map(np.divide, list(map(np.sum, ls_masked_moleN)), list(map(np.sum, ls_vol_eau))))

        return ls_concN
        # rq: np.sum(S.moleN()) / (np.sum(S.tsw_t * S.m_vox_surf)) # whole soil average N molar concentration (micromole N .L-1)

    def init_memory_EV(self, parSN):
        """ Inititalise attributes corresponding to memory variables and parameters to compute soil evaporation

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        """
        # """  """
        self.Uval = parSN['q0']
        HXs = self.m_teta_fc[0, 0, 0]  # par_sol[str(vsoilnumbers[0])]['teta_fc']  # humidite a la capacite au champ de l'horizon de surface
        self.b_ = bEV(parSN['ACLIMc'], parSN['ARGIs'], HXs)
        self.stateEV = [0., 0., 0.]  # pour le calcul de l'evaporation du sol (memoire du cumul evapore depuis derniere PI)
        # marche seulement pour solN (car faut parSN)


    def update_memory_EV(self, new_vals):
        """ Update 'stateEV' attribute with new list of memory variables

        """
        self.stateEV = new_vals


    def updateTsol(self, Tair, optTsol=1):
        """ Update daily soil temperature array 'm_Tsol'

        :param Tair: Daily air temperature
        :type Tair: float
        :param optTsol: Option for computing soil temperature, default to 1 (1:assume Tsol=Tair)
        :type optTsol: int

        """

        if optTsol==1:
            self.m_Tsol = self.m_1 * Tair
        else:
            self.m_Tsol = None
            # A ameliorer avec options bilan d'E! et TCULT!


    #####  soil mineralisation functions #####


    def mask_PROFUM(self, parSN):
        """ Compute a [nz,nx,ny] mask array to apply 'PROFHUMs' mineralisation parameter over soil depth

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        :return: ...
        :rtype: nd.array

        """
        # """ pour creer un mask pour profhum """
        PROFHUM = parSN['PROFHUMs'] / 100.  # in m
        limz = [0.]
        for i in range(len(self.dxyz[2])):
            limz.append(limz[-1] + self.dxyz[2][i])

        limz = array(limz)

        v = limz < PROFHUM
        v = v * 1.
        v = v.tolist()
        idlim = v.index(0)
        limz[idlim - 1]
        v[idlim] = (PROFHUM - limz[idlim - 1]) / (limz[idlim] - limz[idlim - 1])
        v = v[1:]  # mask 1D

        # applique a matrice sol
        res = self.m_1 * 1.
        for i in range(len(v)):
            res[i, :, :] = res[i, :, :] * v[i]

        return res


    def Pot_rate_SOMMin(self, CALCs, ARGIs, parSN):
        """ Compute Potential rate of SOM mineralisation (unit: kg N.day-1)
        (STICS book, Eq. 8.5, p 145)

        :param CALCs: Voxel calcacerous content (unit: %)
        :type CALCs: float
        :param ARGIs: Voxel clay content (unit: %)
        :type ARGIs: float
        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        :return: Voxel potential rate of SOM mineralisation 'K2HUMi' (unit: kg N.day-1)
        :rtype: float
        """

        K2HUMi = parSN['FMIN1G']*np.exp(-parSN['FMIN2G']*ARGIs)/(1+parSN['FMIN3G']*CALCs)
        return K2HUMi
        #? revoir ARGIs et CALCs depuis objet sol -> non laisse la possibilite d'adapter le calcul en fonction de profondeur

    def SOMMin_RespT(self, parSN):
        """ Compute response of mineralisation rate to soil temperature
        (STICS book, Eq. 8.3, p 143) - described as a sigmoid process (FTH)

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        :return: 'FTH', a relative factor  accounting for temperature effect on basal mineralisation (0-1 fraction)
        :rtype: float
        """
        #""" reponse de la mineralisation (ammonification) de la SOM a la temperature - described as a sigmoid process (FTH) - Eq 8.3 p143 et (pour residus FTR p146-147) """
        if self.m_Tsol.min()<=0.:#min(min(min(self.m_Tsol)))<=0.:
            FTH=0.*self.m_1
        else:
            FTH = parSN['FTEMHAg'] / (1 + parSN['FTEMHB']*np.exp(-parSN['FTEMHg']*self.m_Tsol))
        return FTH 

    def SOMMin_RespHum(self, parSN):
        """ Compute response of mineralisation rate to soil relative humidity 'HRv'
        (STICS book, Eq. 8.2, p 143)

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        :return: 'FH', a relative factor  accounting for soil water effect on basal mineralisation (0-1 fraction)
        :rtype: float
        """
        #""" reponse de la mineralisation (ammonification) de la SOM a l'humidite relative (FH) - Eq. 8.2 p 143 - aussi utilise pour residus """
        HR = self.HRv()
        FH = (HR - parSN['HMinMg']*100) / (parSN['HoptMg']*100 - parSN['HMinMg']*100)
        for i in range(len(FH)):
            for j in range(len(FH[i])):
                for k in range(len(FH[i][j])):
                    if FH[i][j][k]<0.:
                        FH[i][j][k]=0.
                    elif FH[i][j][k]>1.:
                        FH[i][j][k]=1.
        return FH

    def Act_rate_SOMMin(self, parSN):
        """ Compute Actual rate of SOM mineralisation (unit: kg N.day-1)
        (STICS book, Eq. 8.1, p 142)

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        :return: 'K2i', basal voxel mineralisation rate (unit: kg N.day-1)
        :rtype: float
        
        """
        return self.K2HUM*self.SOMMin_RespHum(parSN)*self.SOMMin_RespT(parSN)

    def stepNB(self, parSN):
        """ Compute daily step for SOM mineralisation and update SoilN object
        
        """
        #Mineralisation of soil organic matter
        NHUM = self.Norg - self.InertNorg#* (1-par['FINERTG']) #!!! *PROFHUMs : mask pour retirer profondeur ou mineralisation pas significative???
        dN_NH4 = NHUM * self.Act_rate_SOMMin(parSN)
        dC_NH4 = dN_NH4*self.m_CNHUM
        self.Norg = self.Norg - dN_NH4
        self.Corg = self.Corg - dC_NH4 #suppose CN constant/pas affecte par SOM mineralisation!
        self.m_MO = self.Corg*1.72 # pas vraiment utile
        self.m_NH4 = self.m_NH4 + dN_NH4 #verif unite!!! ici en kg d'N dans le voxel de sol ->OK
        #Update bilanCN
        self.bilanN['cumMinN'].append(sum3(dN_NH4) / self.soilSurface *10000)
        self.bilanC['cumMinC'].append(sum3(dC_NH4) / self.soilSurface *10000)


    #####  Residues mineralisation functions #####


    def init_residues(self, parSN, vCNRESt=[], vAmount=[], vProps=[], vWC=[], vCC=[], forced_Cres=None):
        """ Initialise 'parResi' attribute storing the properties of soil residues

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict
        :param vCNRESt: A list of CN ratio for the different types of residues, default to empty list
        :type vCNRESt: list
        :param vAmount: A list of fresh amounts for the different types of residues (unit: T fresh matter.ha-1), default to empty list
        :type vAmount: list
        :param vProps: A list of residue proportion by soil layer for the different types of residues (0-1 fractions), default to empty list
        :type vProps: list
        :param vWC: A list of water content for the different types of residues (0-1 fraction), default to empty list
        :type vWC: list
        :param vCC: A list of carbon content for the different types of residues (0-1 fraction), default to empty list
        :type vCC: list
        :param forced_Cres: A list  of [nz,nx,ny] arrays to use as forced input to 'ls_CRES', default to None
        :type forced_Cres: list


        .. code-block:: python

            self.parResi = {
                        'CNRESt': [],
                        'CNBio': [],
                        'KRES': [],
                        'YRES': [],
                        'KBio': [],
                        'HRES': [],
                        'TRefg': 15.,
                        'FTEMHAg': 12.,
                        'FTEMHg': 0.103,
                        'FTEMHB': 52.,
                        }

        """
        #""" initialisation des compartiments en relation avec gestion des residus """
        #dictionnaire de parametres des residus
        self.parResi = {}
        self.parResi['CNRESt'], self.parResi['CNBio'], self.parResi['KRES'], self.parResi['YRES'], self.parResi['KBio'], self.parResi['HRES'] = [],[],[],[],[],[]

        # same temperature and humidity response for all residues (Recous et al., 1995, p147)
        self.parResi['TRefg'] = parSN['TRefg']#15. #reference temperature
        self.parResi['FTEMHAg'] = 12. #asymptotic value of FTR (seen as a logistic response) p147
        self.parResi['FTEMHg'] = 0.103 #(K-1) p147
        self.parResi['FTEMHB'] = 52. #from eq 8.3

        for val in vCNRESt:
            self.addResPAR(self.parResi, val)

        #liste de matrices pour CRES et CBio
        self.ls_CRES, self.ls_CBio = [], []
        self.bilanN['NminfromNres'] = []  # liste de liste de deltaNmin produits par residus pour les bilans
        for i in range(len(vAmount)):
            self.addResMat(vAmount[i], vProps[i], vWC[i], vCC[i], forced_Cres) #utilise fonction de distribution verticale par default, mais peut forcer une matrice donnee

        #update bilanN
        self.bilanN['initialNres'] = sum(self.ls_NRES()) / self.soilSurface * 10000

        #CO2 resp
        self.CO2respSoil = 0. #separer (kg de C par le volume total de sol) #rq: suppose aucune recapture


    def addResPAR(self, dic, CNRes):
        """ Add a new set of residue parameters to input dictionary  -
        Parameters calculated according to Nicolardot et al. (2001) as functions of residue CN ratio

        :param dic: A dictionnary of type 'parResi'
        :type dic: dict
        :param CNRes: CN ratio of the new residue
        :type CNRes: float

        :return: Updated dictionnary of type 'parResi'
        :rtype: dict
        """
        #""" add a new series of parammeters for a residue according to Nicolardot et al. (2001) in function of its CN ratio"""
        dic['CNRESt'].append(CNRes) #CSURNRESt C/N des residus #liste pour les different residus
        dic['CNBio'].append(max(7.8, 16.1-123./CNRes)) #CN ratio of the zymogeneous biomass -> fontion du CN des residus (eq. 8.6) -> des biomasses microbiennes associee a chaque type de residu?
        dic['KRES'].append(0.07+1.94/CNRes) # decomposition rate constant (day-1- normalised day at 15dC) from organic residue to microbil biomass?  (fig 8.4 p146) -> fonction de CNRESt
        dic['YRES'].append(0.62) # Assimilation yield of residue-C by microbial biomass - partition parameter between CO2 and microbil biomass (fig 8.4 p146) -> constant
        dic['KBio'].append(0.0110) # decomposition rate constant (day-1 - normalised day at 15dC) from microbial biomass to humus? (fig 8.4 p146) -> constant
        dic['HRES'].append(1-(0.69*CNRes)/(11.2+CNRes)) #  partition parameter between CO2 and humus - Humification rate of microbial biomass -(fig 8.4 p146) -> fonction de CNRESt
        return dic


    def VdistribResidues(self, Amount, Vprop, Wcontent=0.8,Ccontent=0.42):
        """ Distribute a new residue in a [nz,nx,ny] 'm_CRES' array assuming horizontal homogeneity (unit: kg C per Voxel)

        :param Amount: Fresh amount of the residue (unit: T fresh matter.ha-1)
        :type Amount: float
        :param Vprop: List of nz residue proportion by vertical soil layer
        :type Vprop: list
        :param Wcontent: Residue water content (0-1 fraction of fresh weight), default to 0.8
        :type Wcontent: float
        :param Ccontent: Residue carbon content (0-1 fraction of dry weight), default to 0.42
        :type Ccontent: float

        :return: 'm_CRES', a [nz,nx,ny] array describing residue C distribution into the soil (unit: kg C per Voxel)
        :rtype: nd.array

        """
        #"""initialise CRES (Amount of decompasble C in the residue) for a given amount/gratient """
        #Amount (T fresh matter.ha-1)
        #Wcontent (water content; proportion)
        #Ccontent (carbon content; propostion)
        #Vprop : list of proportion per horizon (distribution assumed homogeneous for a given horizon) -> same number of elements than dz
        
        Q = (Amount/10000.)*self.soilSurface*1000.  #kg of fresh matter
        QC = Q*(1-Wcontent)*Ccontent #kg of C
        nbcase_layer = len(self.dxyz[0])*len(self.dxyz[1])

        m_CRES = []
        for z in range(len(self.dxyz[2])):
            v = []
            for x in range(len(self.dxyz[0])):
                vv = []
                for y in range(len(self.dxyz[1])):
                    vv.append(Vprop[z]*QC/nbcase_layer)

                v.append(vv)
            m_CRES.append(v)

        m_CRES = array(m_CRES)
        return m_CRES
        

    def addResMat(self, Amount, Vprop, Wcontent=0.8,Ccontent=0.42, forced_Cres=None):
        """ Update 'ls_CRES', 'bilanC' and 'bilanN' attributes when adding a new residue

        :param Amount: Fresh amount of the residue (unit: T fresh matter.ha-1)
        :type Amount: float
        :param Vprop: List of nz residue proportion by vertical soil layer
        :type Vprop: list
        :param Wcontent: Residue water content (0-1 fraction of fresh weight), default to 0.8
        :type Wcontent: float
        :param Ccontent: Residue carbon content (0-1 fraction of dry weight), default to 0.42
        :type Ccontent: float
        :param forced_Cres: A [nz,nx,ny] array to use as forced input to 'ls_CRES', default to None
        :type forced_Cres: nd.array
        
        """
        #""" add matrice associee au C des residus (CRES) et de leur microbial biomass """
        self.ls_CBio.append(0*self.m_1)#matrice for the Amount of C in the microbial biomass (kg C in the voxel) -> initialise a zero
        if forced_Cres == None:
            cres = self.VdistribResidues(Amount, Vprop, Wcontent, Ccontent)
            self.ls_CRES.append(cres)
            #bilan
            self.bilanC['initialCres'] += sum(cres)/ self.soilSurface *10000
            self.bilanN['NminfromNres'].append([])  # ajoute une liste vide pour le residu
        else: #pour gerer le cas ou passera une matrice (3D) deja faite a partir de maquettes
            self.ls_CRES.append(forced_Cres)
            #bilan
            self.bilanC['initialCres'] += sum(forced_Cres)/ self.soilSurface *10000
            self.bilanN['NminfromNres'].append([])  # ajoute une liste vide pour le residu


    def Pot_rate_ResidueMin(self, res_id, parSN):
        """ Compute the potential rate of residue mineralisation consirering soil water and temperature effects (unit: kg C.day-1)
        (STICS book, Eq. 8.7, p 146)

        :param res_id: residue ID in 'parResi' dict
        :type res_id: int
        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        :return: Potential rate of residue mineralisation (unit: kg C.day-1)
        :rtype: float

        """
        #""" changement potentiel en C du residu id (dans conditions donnes d'humidite et T) eq. 8.7 p 146"""

        pDCRES = -self.parResi['KRES'][res_id] * self.ls_CRES[res_id] * self.SOMMin_RespHum(parSN) * self.SOMMin_RespT(self.parResi) #sans *FN (dispo en N) ->  pot microbial growth dans ces condition
        return pDCRES
        #faudrait lire les FH plutot que de les recalculer a chaque fois
        #laisser en valeur positive?

    def Pot_Ndemand_microbialBio(self, parSN):
        """ Compute total microbial N demand for all residues (unit: kg de N per voxel)

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        :return: an array of size [nz,nx,ny] storing total microbial N demand to support residue mineralisation (unit: kg de N per voxel)
        :rtype: nd.array
        
        """
        #"""  demande totale en azote pour atteindre croissance optimale microbio de tous les residus (kg de N par voxel) - somme des demandes pour chaque residu """
        #parSN: pour reponse a humidite
        res = self.m_1*0
        for i in range(len(self.parResi['KRES'])):
            pDCRESi = -self.Pot_rate_ResidueMin(i, parSN)
            CBio_poti = pDCRESi *   self.parResi['YRES'][i]
            deltaNi = CBio_poti/self.parResi['CNBio'][i] - pDCRESi/self.parResi['CNRESt'][i] #ce qu'il manque entre N issu de biomasse decomposee et N requis pour microbio! (>0 de ce qu'il faut en N pour atteindre potentiel)
            res = res + deltaNi

        return res

    def FN_factor(self, parSN):
        """ Compute 'FN' reduction factor related to Nmin availability to support residue mineralisation

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict
        
        """
        #""" calcul du facteur de reduction lie a la disponibilite en Nmin a proximite des residus """
        #!! demande des plantes pas prise en compte -> servie slmt s'il en reste apres microbio!
        MND = self.Pot_Ndemand_microbialBio(parSN)
        Nmin = self.m_NH4 + self.m_NO3
        delta = Nmin - MND #<0 -> manque de N 
        ratio = Nmin/(MND+1e-15)
        FN = deepcopy(delta)
        for i in range(len(FN)):
            for j in range(len(FN[i])):
                for k in range(len(FN[i][j])):
                    if FN[i][j][k]>=0.:#pas de limitation
                        FN[i][j][k] = 1.
                    else:#<0 -> manque de N -> au max prend Nmin d'ou reduction de Nmin/MND
                        FN[i][j][k] = ratio[i][j][k] #rajouter un petit% de securite?
        return FN

    def FBIO_factor(self):
        """  Compute 'FBIO' reduction factor (0-1 fraction), default to 1. -
        (STICS book, Eq. 8.11, p 147)
        
        """
        #""" p147 - a faire et introduire dans stepMicrobioMin"""
        return 1.

    def stepResidueMin(self, parSN):
        """ Compute Daily Carbon and Nitrogen fluxes associated with mineralisation of all residues -

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict
        
        """

        FN = self.FN_factor(parSN) #!! faudrait calculer FN une seule fois pour tous les residus!
        cumNRes1, cumNRes2= [],[]
        for i in range(len(self.ls_CRES)):
            #fux de C
            DCRESi = self.Pot_rate_ResidueMin(i, parSN) * FN
            self.ls_CRES[i] = self.ls_CRES[i] + DCRESi #ou - si garde DCRESi positif
            DCBioi = DCRESi*self.parResi['YRES'][i]
            self.ls_CBio[i] = self.ls_CBio[i] - DCBioi #ou + si garde DCRESi positif
            DCO2 = DCRESi-DCBioi
            self.CO2respSoil = self.CO2respSoil - sum3(DCO2)#ou + si garde DCRESi positif
            #bilan
            self.bilanC['cumCO2Res1'].append(- sum3(DCO2)/ self.soilSurface *10000)

            #flux de N
            #ajouter a NH4 le complement N de CO2 emission?
            DNCO2 = DCO2/self.parResi['CNRESt'][i]
            self.m_NH4 = self.m_NH4 - DNCO2 #ou + si garde DCRESi positif
            Nmin = self.m_NH4 + self.m_NO3 + 0.000000000001 #pour eviter les Nan dans f_NH4 et f_NO3
            Ndemandi = -DCBioi/self.parResi['CNBio'][i] + DCBioi/self.parResi['CNRESt'][i] #verif signes!
            #preleve dans NH4 et NO3 a hauteur de leur contribution a Nmin + verif que aucun devient negatif
            f_NH4 = self.m_NH4/Nmin
            f_NO3 = self.m_NO3/Nmin
            self.m_NH4 = self.m_NH4 - f_NH4*Ndemandi
            self.m_NO3 = self.m_NO3 - f_NO3*Ndemandi
            #bilan
            cumNRes1.append(- sum3(DNCO2) / self.soilSurface *10000)
            cumNRes2.append(- sum3(Ndemandi)/ self.soilSurface *10000)
            self.bilanN['NminfromNres'][i].append(- sum3(DNCO2) / self.soilSurface *10000)

        self.bilanN['cumNRes1'].append(sum(cumNRes1))#pour ajouter uniquement cumul journalier de tous les residus
        self.bilanN['cumNRes2'].append(sum(cumNRes2))


    def stepMicrobioMin(self, parSN):
        """ Compute Daily Carbon and Nitrogen fluxes associated with microbial biomass of all residues -

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        """
        #parSN : pour reponse a humidite
        cumNRes3= []
        for i in range(len(self.ls_CBio)):
            #fux de C
            DCBioi = self.parResi['KBio'][i] * self.ls_CBio[i] * self.SOMMin_RespHum(parSN) * self.SOMMin_RespT(self.parResi)
            self.ls_CBio[i] = self.ls_CBio[i] - DCBioi
            DCHUM = DCBioi * self.parResi['HRES'][i]
            self.Corg = self.Corg + DCHUM
            DCO2 = DCBioi - DCHUM
            self.CO2respSoil = self.CO2respSoil + sum3(DCO2)
            #bilan
            self.bilanC['cumCO2Res2'].append(sum3(DCO2)/ self.soilSurface *10000)

            #flux de N
            #equilibre Norg 
            self.Norg = self.Norg + DCHUM/self.m_CNHUM
            #ajouter a NH4 le complement N de CO2 emission? + ajouter N lie au changement de CN de DCHUM
            Neccesi = DCHUM/self.parResi['CNBio'][i] - DCHUM/self.m_CNHUM #en ecces par rapport a CN de MOS #!!! si <0 pourrait potentiellement faire passer Nmin <0 -> a ajuster avec facteur FBIO!
            self.m_NH4 = self.m_NH4 + DCO2/self.parResi['CNBio'][i] + Neccesi 
            
            #bilan
            cumNRes3.append(sum3(DCO2/self.parResi['CNBio'][i] + Neccesi)/ self.soilSurface *10000)

        self.bilanN['cumNRes3'].append(sum(cumNRes3))

        #facteur FBIO!! = 1 par defaut -> peut augmenter <-> baisse de CN ratio de la biomasse microbienne si Nmin exhauted (pour eviter bilan negatif)
        #a faire! avec approche similaire a FN_factor!


    def ls_NRES(self):
        """ Compute a list of [nz,nx,ny] arrays storing residue N distribution per type of residue (unit: kg N per voxel)

        """
        #""" calcul N dans les residus - liste equivalente de ls_CRES """
        lsNRES = []
        for i in range(len(self.ls_CRES)):
            lsNRES.append(self.ls_CRES[i]/self.parResi['CNRESt'][i])
        return lsNRES

    def ls_NBio(self):
        """ Compute a list of [nz,nx,ny] arrays storing microbial N distribution per type of residue (unit: kg N per voxel)
        
        """
        #""" calcul N dans les biomasse microbienne de residus - liste equivalente de ls_CBio """
        lsNbio = []
        for i in range(len(self.ls_CBio)):
            lsNbio.append(self.ls_CBio[i]/self.parResi['CNBio'][i])
        return lsNbio

    def mixResMat(self, mat_res, idres, Ccontent=0.42):
        """ Mix a new input residue [nz,nx,ny] array with an existing residue type in 'parResi' and 'ls_CRES' attributes

        :param mat_res: An input residue array of size [nz,nx,ny] (unit: g dry matter per voxel)
        :type mat_res: nd.array
        :param idres: residue ID in 'parResi' dict
        :type idres: int
        :param Ccontent: Residue carbon content (0-1 fraction of dry weight), default to 0.42
        :type Ccontent: float

        """
        #""" add matrice associee des residus (en g MS par voxl) au C (CRES) d' un residu avec id deja existant """
        # suppose CsurN comme residu existant

        cres = mat_res * Ccontent / 1000.  # conversion #kg of C per voxel
        self.ls_CRES[idres] += cres
        # bilan
        self.bilanC['initialCres'] += sum(cres) / self.soilSurface * 10000
        self.bilanN['initialNres'] += sum(self.ls_CRES[idres] / self.parResi['CNRESt'][idres]) / self.soilSurface * 10000

    #  a faire (?) cree noueau residu si tres different? ajuster bilan C/N # faire evoluer lrd parametres des reisdus selon C/N vrai?


    #####  soil nitrification functions #####


    def Nitrif_RespHum(self, parSN):
        """ Compute response of nitrification rate to soil Humidity 'HRv' (0-1 fraction) -
        (STICS book, Eq. 8.14, p 150)

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        :return: 'FHN', a relative factor  accounting for soil humidity on nitrification (0-1 fraction)
        :rtype: float
        
        """
        #""" reponse de la nitrification de NH4+ a l'humidite relative (FHN) - Eq. 8.14 p 150 - increasing sigmoid-like curve """
        HR = self.HRv()
        FHN = (HR - parSN['HMinNg']*100) / (parSN['HoptNg']*100 - parSN['HMinNg']*100)
        for i in range(len(FHN)):
            for j in range(len(FHN[i])):
                for k in range(len(FHN[i][j])):
                    if FHN[i][j][k]<0.:
                        FHN[i][j][k]=0.
                    elif FHN[i][j][k]>1.:
                        FHN[i][j][k]=1.
        return FHN

    def Nitrif_RespPH(self, parSN):
        """ Compute response of nitrification rate to soil pH (0-1 fraction) -
        (STICS book, Eq. 8.13, p 150)

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        :return: 'FPHN', a relative factor  accounting for soil pH on nitrification (0-1 fraction)
        :rtype: float
        
        """
        #""" reponse de la nitrification de NH4+ au pH (FH) - Eq. 8.13 p 150 - increasing sigmoid-like curve """
        pHs = self.pHeau
        FPHN = (pHs - parSN['PHMinNITg']) / (parSN['PHMaxNITg'] - parSN['PHMinNITg'])
        if FPHN<0.:
            FPHN=0.
        elif FPHN>1.:
            FPHN=1.
        return FPHN #a priori scalaire a ne calculer qu'une fois (comme pH change pas)

    def Nitrif_RespT(self, parSN):
        """ Compute response of nitrification rate to soil temperature (0-1 fraction) -
        (STICS book, Eq. 8.15, p 151)

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict

        :return: 'FTN', a relative factor  accounting for soil temperature on nitrification (0-1 fraction)
        :rtype: float
        
        """
        #""" reponse de la nitrification de NH4+ a la temperature (FTN) - Eq. 8.15 p 151 - bilinear beta-like curve """
        FTN = deepcopy(self.m_Tsol)
        for i in range(len(FTN)):
            for j in range(len(FTN[i])):
                for k in range(len(FTN[i][j])):
                    Tsol = FTN[i][j][k]
                    if Tsol <= parSN['TNITOPTg']:
                        ratio = (Tsol - parSN['TNITMINg']) / (parSN['TNITOPTg'] - parSN['TNITMINg'])
                        FTN[i][j][k] = max(0., ratio)
                    else:#superieur a Topt
                        ratio = (Tsol - parSN['TNITMAXg']) / (parSN['TNITOPTg'] - parSN['TNITMAXg'])
                        FTN[i][j][k] = max(0., ratio)
        return FTN

    def stepNitrif(self, parSN):
        """ Compute daily N fluxes associated to NH4+ nitrification -
        (STICS book, Eq. 8.12 and 8.16, p 151)

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict
        
        """

        TNITRIF = self.m_NH4 * parSN['FNXg'] * self.Nitrif_RespHum(parSN) * self.Nitrif_RespPH(parSN) * self.Nitrif_RespT(parSN)
        NITRIF = (1 - parSN['RATIONITs']) * TNITRIF
        self.m_NO3 = self.m_NO3 + NITRIF
        self.m_NH4 = self.m_NH4 - TNITRIF
        dN2ONitrif = sum3(TNITRIF-NITRIF)
        self.N2ONitrif = self.N2ONitrif + dN2ONitrif
        self.bilanN['cumN2O'].append(dN2ONitrif / self.soilSurface *10000) #UpdateNminbalance


    #####  Nitrates infiltration functions #####


    def infil_layerNO3(self, in_N, out_Water , idz, opt_infil=1):
        """ Compute nitrate infiltration for soil layer idz

        :param in_N: Input infiltration N-NO3 array of size [nx,ny] coming from (idz-1) layer (unit: kg N per voxel)
        :type in_N: nd.array
        :param out_Water: Output infiltration water array of size [nx,ny] going to (idz+1) layer (unit: mm)
        :type out_Water: nd.array
        :param idz: z value of voxel/layer ID in grid [nz, nx, ny]
        :type idz: int
        :param opt_infil: Option indicating if infiltration is vertical (=1) or distributed between first order voxels of the layer below (=2), default=1
        :type opt_infil: int

        :return: new, out_N

        """
        new = self.m_NO3[idz] + in_N
        propNO3 = out_Water / (self.m_QH20max[idz] + out_Water) #prop de nitrate qui part est fraction du volume d'eau max qui passe (jamais >1)
        #putmask(propNO3, propNO3>1. ,1.)#!!!verif pas superieur a 1 et sinon remplace par 1!!  syntaxe interessante
        out_N = new*propNO3
        new = new-out_N

        if opt_infil==2: #ditribution pas juste verticale
            out_N2 = deepcopy(out_N)
            out_N2.fill(0.)
            for x in range(len(out_N2)):
                for y in range(len(out_N2[x])):
                    q_out = out_N[x][y]
                    ls_v = self.ls_1storder_vox(x, y, idz, opt_infil) #distribution entre les 1st order ; mettre opt=1 si veut forcer verticalement / 2 si
                    if len(ls_v)>1:
                        ponder = [0.0416666, 0.0416666, 0.0416666, 0.0416666, 2/3., 0.0416666, 0.0416666, 0.0416666, 0.0416666]# 2/3 en dessous 1/3 au premier ordre
                    else:
                        ponder = [1.]

                    for i in range(len(ls_v)):#distribution dans les voxels du dessous
                        nx, ny = ls_v[i][0], ls_v[i][1]
                        out_N2[nx][ny] = out_N2[nx][ny]+ponder[i]*q_out  

            out_N = out_N2

        return new, out_N


    def distrib_NO3(self, map_N, ls_outWater, opt_infil=1): #map_N = map application nitrates en surface
        """ Distribute nitrate infiltration into the soil grid from mineralisation and/or fertilizer distribution map at the soil surface

        :param map_N: Input infiltration N-NO3 array at soil surface (unit: kg N per voxel)
        :type map_N: nd.array
        :param ls_outWater: List of nz output infiltration water array of size [nx,ny] by layer (unit: mm)
        :type ls_outWater: list
        :param opt_infil: Option indicating if infiltration is vertical (=1) or distributed between first order voxels of the layer below (=2), default=1
        :type opt_infil: int

        :return: matNO3_t, out_N
        """
        
        in_N = map_N
        matNO3_t = deepcopy(self.m_NO3)
        #ls_out = []
        for z in range(len(matNO3_t)):
            new, out_N = self.infil_layerNO3(in_N, ls_outWater[z], z, opt_infil)
            #ls_out.append(out_)
            in_N = out_N
            matNO3_t[z] = new

        return matNO3_t, out_N
        #A faire: distinguer les entree N pluie, irrigation, fertilisation


    def stepNINFILT(self, mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, ls_outWater, opt_infil=1):#(self, map_N, ls_outWater, opt=1):
        """ Compute daily N fluxes associated to nitrate infiltration and fertilizer application

        :param mapN_Rain: Input N from rain in a array of size [nx,ny] at soil surface (unit: kg N per voxel)
        :type mapN_Rain: nd.array
        :param mapN_Irrig: Input N from irrigation water in a array of size [nx,ny] at soil surface (unit: kg N per voxel)
        :type mapN_Irrig: nd.array
        :param mapN_fertNO3: Input N from NO3- fertilizers in a array of size [nx,ny] at soil surface (unit: kg N per voxel)
        :type mapN_fertNO3: nd.array
        :param mapN_fertNH4: Input N from NH4+ fertilizers in a array of size [nx,ny] at soil surface (unit: kg N per voxel)
        :type mapN_fertNH4: nd.array
        :param ls_outWater: List of nz output infiltration water array of size [nx,ny] by layer (unit: mm)
        :type ls_outWater: list
        :param opt_infil: Option indicating if infiltration is vertical (=1) or distributed between first order voxels of the layer below (=2), default=1
        :type opt_infil: int

        """

        #ajout N NO3 mobile
        map_N = mapN_Rain + mapN_Irrig + mapN_fertNO3 #+ mapN_fertNH4
        matNO3_t, out_NO3 = self.distrib_NO3(map_N, ls_outWater, opt_infil)
        self.m_NO3 = matNO3_t
        Lix = sum(sum(out_NO3))
        self.lixiNO3 = self.lixiNO3 + Lix

        #ajout N NH4 non mobile dans 1ere couche
        self.m_NH4[0,:,:] = self.m_NH4[0,:,:] + mapN_fertNH4

        #bilans
        self.bilanN['cumRain'].append(sum(mapN_Rain)/ self.soilSurface *10000) 
        self.bilanN['cumIrrig'].append(sum(mapN_Irrig)/ self.soilSurface *10000)
        self.bilanN['cumfertNO3'].append(sum(mapN_fertNO3)/ self.soilSurface *10000) 
        self.bilanN['cumfertNH4'].append(sum(mapN_fertNH4)/ self.soilSurface *10000)
        self.bilanN['cumLix'].append(Lix / self.soilSurface *10000) 
     

    #####  Plant N uptake functions #####

    def stepNuptakePlt(self, parSN, paramp=[{}], ls_lrac=None, ls_mWaterUptakePlt=None, ls_demandeN=None, optNuptake=1):
        """ Compute daily N fluxes associated to plant mineral N uptake

        :param parSN: A dictionnary defining general soil parameters derived from STICS
        :type parSN: dict
        :param paramp: List of individual plant parameters
        :type paramp: list
        :param ls_lrac: List of [z,x,y] arrays for root length distribution per individual plant - equivalent to 'ls_roots' (unit: m)
        :type ls_lrac: list
        :param ls_mWaterUptakePlt: A [nbplt,z,x,y] nd.array storing individual plant daily water uptake distributions (unit: mm) - equivalent to 'ls_m_transpi' (unit: mm)
        :type ls_mWaterUptakePlt: nd.array
        :param ls_demandeN: List of plant Nitrogen status per individual plant ; for optNuptake=0, a list of plant N demand (kg N.plt-1); for optNuptake=1, a list of root N concentration (unit: %) or NNI (unitless)
        :type ls_demandeN: list
        :param optNuptake: Options for computing plant mineral N uptake, either 0: 'original STICS' (minimum of plant demand, plant uptake capacity and soil provision) or 1: 'LocalTransporters' (plant uptake from root transporters)
        :type optNuptake: int

        :return: ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N, idmin
        """
        # ls_lrac = ls_roots!
        # """ calculation of actual N uptake by plant - if no plants (baresoil) -> let None in ls_rac,ls_mWaterUptakePlt, ls_demandeN """

        if optNuptake == 0:  # 'STICS':
            # version initiale tiree de STICS avec correction unites concN
            if ls_lrac is None or ls_mWaterUptakePlt is None or ls_demandeN is None:  # si pas de plante (au moins fournir le paramp qui donne un nbre de plante)
                ActUpNtot = self.m_1 * 0.
                ls_Act_Nuptake_plt = [self.m_1 * 0.] * len(paramp)
                ls_DQ_N = [1.] * len(paramp)
                idmin = self.m_1 * -1.
            else:  # si plante
                PotUpNtot, ls_Pot_Nuptake_plt, idmin = Distrib_Potential_Nuptake_Plt(self, parSN, paramp, ls_lrac, ls_mWaterUptakePlt)
                ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N = Actual_Nuptake_plt(self, ls_Pot_Nuptake_plt, ls_demandeN)
        elif optNuptake == 1:  # 'LocalTransporter':
            # version reponse locale racine-trasporter
            if ls_lrac is None or ls_mWaterUptakePlt is None or ls_demandeN is None:  # si pas de plante (au moins fournir le paramp qui donne un nbre de plante)
                ActUpNtot = self.m_1 * 0.
                ls_Act_Nuptake_plt = [self.m_1 * 0.] * len(paramp)
                ls_DQ_N = [1.] * len(paramp)
                idmin = self.m_1 * -1.
            else:  # si plante
                ls_PltN = ls_demandeN  # avec cette option doit etre ls valeur de NNI
                PotUpNtot, ls_Pot_Nuptake_plt, idmin = Distrib_Potential_Nuptake_Plt_Bis(self, paramp, ls_lrac)
                ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N = Actual_Nuptake_plt_Bis(self, ls_Pot_Nuptake_plt, ls_PltN, paramp)
        elif optNuptake == 2:  # 'old':
            # version ancienne (bug concN)
            if ls_lrac is None or ls_mWaterUptakePlt is None or ls_demandeN is None:  # si pas de plante (au moins fournir le paramp qui donne un nbre de plante)
                ActUpNtot = self.m_1 * 0.
                ls_Act_Nuptake_plt = [self.m_1 * 0.] * len(paramp)
                ls_DQ_N = [1.] * len(paramp)
                idmin = self.m_1 * -1.
            else:  # si plante
                ls_PltN = ls_demandeN  # avec cette option doit etre ls valeur de NNI
                PotUpNtot, ls_Pot_Nuptake_plt, idmin = Distrib_Potential_Nuptake_Plt_old(self, parSN, paramp, ls_lrac, ls_mWaterUptakePlt)
                ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N = Actual_Nuptake_plt_old(self, ls_Pot_Nuptake_plt, ls_PltN)

        # retire les nitrates et ammomium rellement preleves du sol
        frac_NO3 = self.m_NO3 / (self.m_NO3 + self.m_NH4 + 10 ** -15)
        self.m_NO3 = self.m_NO3 - frac_NO3 * ActUpNtot
        self.m_NH4 = self.m_NH4 - (1. - frac_NO3) * ActUpNtot

        # bilan
        self.bilanN['cumUptakePlt'].append(sum(ActUpNtot) / self.soilSurface * 10000)
        self.bilanN['azomes'].append((sum(self.m_NO3) + sum(self.m_NH4)) / self.soilSurface * 10000)

        return ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N, idmin


    #####  Nitrogen and Carbon balance #####


    def OpenCbalance(self):
        """ Initialise bilanC attribute, a dictionary storing carbon balance information (unit: kg C.ha-1)

            * Keys for simulation inputs:
                - 'intialInertC': Initial soil organic C in the inert humified pool (unit: kg C.ha-1)
                - 'intialActiveC': Initial soil organic C in the active humified pool (unit: kg C.ha-1)
                - 'initialCres': Initial soil organic C in the organic residue pools (unit: kg C.ha-1)
                - 'initialCZygo': Initial soil organic C in the Zymogeneous biomass pools (unit: kg C.ha-1)
            * Keys for simulation outputs:
                - 'FinalInertC': Final soil organic C in the inert humified pool (unit: kg C.ha-1)
                - 'FinalActiveC': Final soil organic C in the active humified pool (unit: kg C.ha-1)
                - 'finalCres': Final soil organic C in the organic residue pools (unit: kg C.ha-1)
                - 'finalCZygo': Final soil organic C in the Zymogeneous biomass pools (unit: kg C.ha-1)
            * Keys for simulation totals:
                - 'InputCtot':  Total simulation organic C inputs (unit: kg C.ha-1)
                - 'OutputCtot': Total simulation organic C outputs (unit: kg C.ha-1)
                - 'MinCtot': Total simulation net C produced by mineralisation (humus+residues) (unit: kg C.ha-1)
            * Keys for daily inputs/outputs :
                - 'cumMinC': List of daily C-CO2 emitted from humus mineralisation (unit: kg C.ha-1)
                - 'cumCO2Res1': List of daily C-CO2 emitted from direct residue mineralisation (unit: kg C.ha-1)
                - 'cumCO2Res2': List of daily C-CO2 emitted from residue microbial biomass turnover (unit: kg C.ha-1)

        """
        #""" Dictionnary for soil Carbon balance (kg C.ha-1)
        #Keys for Daily outputs: 'cumMinC', 'cumCO2Res1', 'cumCO2Res2'
        #Keys for Total Input: 'intialInertC', 'intialActiveC', 'initialCZygo', 'initialCres'
        #Keys for Total outputs: 'FinalInertC', 'FinalActiveC', 'MinCtot'
        #Keys for totals: 'InputCtot', 'OutputCtot'
        #"""
        
        surfsolref = self.soilSurface
        self.bilanC = {}
        #Humus
        self.bilanC['intialInertC'] = sum3(self.InertCorg) / surfsolref *10000
        self.bilanC['intialActiveC'] = sum3(self.Corg-self.InertCorg) / surfsolref *10000
        self.bilanC['cumMinC'] = []
        #residus
        self.bilanC['cumCO2Res1'], self.bilanC['cumCO2Res2'] = [], []
        try:
            self.bilanC['initialCZygo'] = sum(self.ls_CBio)/ surfsolref *10000
            self.bilanC['initialCres'] = sum(self.ls_CRES)/ surfsolref *10000
        except:#si pas de residus ajoute
            self.bilanC['initialCZygo'] = 0.
            self.bilanC['initialCres'] = 0.

        # autres termes avec microbio et residus a ajouter!

    #def UpdateCbalance(self, dCMin):
    #    surfsolref = self.soilSurface
    #    self.bilanC['cumMinC'].append(sum3(dCMin) / self.soilSurface *10000)
    #    # autres termes avec microbio et residus a ajouter!!

    def CloseCbalance(self, print_=1):
        """ Finalise calculation of whole simulation variables in bilanC attribute
        
        """
        surfsolref = self.soilSurface
        #Humus Mineralisation
        #input  
        self.bilanC['InputCtot'] = self.bilanC['intialInertC'] + self.bilanC['intialActiveC'] + self.bilanC['initialCZygo'] + self.bilanC['initialCres']
        #output
        self.bilanC['FinalInertC'] = sum3(self.InertCorg) / surfsolref *10000
        self.bilanC['FinalActiveC'] = sum3(self.Corg-self.InertCorg) / surfsolref *10000

        #Residus
        try:
            self.bilanC['finalCZygo'] = sum(self.ls_CBio)/ surfsolref *10000
            self.bilanC['finalCres'] = sum(self.ls_CRES)/ surfsolref *10000
        except:#si pas de residus ajoute
            self.bilanC['finalCZygo'] = 0.
            self.bilanC['finalCres'] = 0.

        self.bilanC['MinCtot'] = sum(self.bilanC['cumMinC']) + sum(self.bilanC['cumCO2Res1']) + sum(self.bilanC['cumCO2Res2'])
        self.bilanC['OutputCtot'] = self.bilanC['FinalInertC'] + self.bilanC['FinalActiveC'] + self.bilanC['MinCtot'] + self.bilanC['finalCZygo'] + self.bilanC['finalCres']

        if print_==1:
            self.PrintCbalance()
        #pourrait le diriger vers un fichier de sortie texte?

    def PrintCbalance(self):
        """ Print a summary table of bilanC attribute
        
        """
        bilanC = self.bilanC
        print ("")
        print ("Carbon Balance Input (kg C.ha-1)\t\t\t Carbon Balance Output (kg C.ha-1)")
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Active Humified Pool:\t {0:8.1f}\t\t\t Active Humified Pool:\t {1:8.1f}".format(bilanC['intialActiveC'], bilanC['FinalActiveC'])))
        print(("Inert Humified Pool:\t {0:8.1f}\t\t\t Inert Humified Pool:\t {1:8.1f}".format(bilanC['intialInertC'], bilanC['FinalInertC'])))
        print(("Zymogeneous Bio Pool:\t {0:8.1f}\t\t\t Zymogeneous Bio Pool:\t {1:8.1f}".format(bilanC['initialCZygo'], bilanC['finalCZygo'])))
        print(("Added organic matter:\t {0:8.1f}\t\t\t Added organic matter:\t {1:8.1f}".format(bilanC['initialCres'], bilanC['finalCres'])))
        print(("                            \t\t\t\t Mineralisation:\t\t {0:8.1f}".format(bilanC['MinCtot'])))
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Total:\t\t\t\t\t {0:8.1f}\t\t\t Total:\t\t\t\t\t {1:8.1f}".format(bilanC['InputCtot'], bilanC['OutputCtot'])))
        print ("")

    def OpenNbalance(self):
        """ Initialise bilanN attribute, a dictionary storing nitrogen balance information (unit: kg N.ha-1)

            * Keys for simulation inputs:
                - 'intialInertN': Initial soil organic N in the inert humified pool (unit: kg N.ha-1)
                - 'intialActiveN': Initial soil organic N in the active humified pool (unit: kg N.ha-1)
                - 'intialNO3': Initial soil mineral N-NO3 (unit: kg N.ha-1)
                - 'intialNH4': Initial soil mineral N-NH4 (unit: kg N.ha-1)
                - 'initialNres': Initial soil organic N in the organic residue pools (unit: kg N.ha-1)
                - 'initialNZygo': Initial soil organic N in the Zymogeneous biomass pools (unit: kg N.ha-1)

                - 'TotNRain': Total simulation mineral N inputs with rain (unit: kg N.ha-1)
                - 'TotNIrrig': Total simulation mineral N inputs with irrigation water (unit: kg N.ha-1)
                - 'TotFertNO3': Total simulation N-NO3 fertilizer application (unit: kg N.ha-1)
                - 'TotFertNH4': Total simulation N-NH4 fertilizer application (unit: kg N.ha-1)

            * Keys for simulation outputs:
                - 'FinalInertN': Final soil organic N in the inert humified pool (unit: kg N.ha-1)
                - 'FinalActiveN': Final soil organic N in the active humified pool (unit: kg N.ha-1)
                - 'FinalNO3': Final soil mineral N-NO3 (unit: kg N.ha-1)
                - 'FinalNH4': Final soil mineral N-NH4 (unit: kg N.ha-1)
                - 'finalNres': Final soil organic N in the organic residue pools (unit: kg N.ha-1)
                - 'finalNZygo': Final soil organic N in the Zymogeneous biomass pools (unit: kg N.ha-1)

                - 'Lixtot': Total simulation N-NO3 lixiviation in drainage water (unit: kg N.ha-1)
                - 'N2Otot': Total simulation N-N2O emissions in the atmosphere (unit: kg N.ha-1)
                - 'TotUptPlt': Total simulation mineral N uptake by plants (unit: kg N.ha-1)
                - 'HumusMinNtot': Total simulation mineral N produced by humus mineralisation (unit: kg N.ha-1)
                - 'ResidueMinNtot': Total simulation net mineral N produced by residue mineralisation (Res1+Res2+Res3) (unit: kg N.ha-1)
                - 'MinNtot': Total simulation net mineral N produced by mineralisation (humus+residues) (unit: kg N.ha-1)

            * Keys for simulation totals:
                - 'InputNtot': Total simulation organic N inputs (unit: kg N.ha-1)
                - 'InputNmintot': Total simulation mineral N inputs (unit: kg N.ha-1)
                - 'OutputNtot': Total simulation organic N outputs (unit: kg N.ha-1)
                - 'OutputNmintot': Total simulation mineral N outputs (unit: kg N.ha-1)
            * Keys for daily inputs/outputs :
                - 'cumMinN': List of daily soil humus mineralisation (unit: kg N.ha-1)
                - 'cumRain': List of daily mineral N inputs with rain (unit: kg N.ha-1)
                - 'cumIrrig': List of daily mineral N inputs with irrigation water (unit: kg N.ha-1)
                - 'cumfertNO3': List of daily mineral N-NO3 fertilizer application (unit: kg N.ha-1)
                - 'cumfertNH4': List of daily mineral N-NH4 fertilizer application (unit: kg N.ha-1)
                - 'cumUptakePlt': List of daily total plant mineral N uptake (unit: kg N.ha-1)
                - 'azomes': List of daily total soil mineral N (unit: kg N.ha-1)
                - 'cumLix': List of daily N-NO3 lixiviation in drainage water (unit: kg N.ha-1)
                - 'cumN2O': List of daily N-N2O emmission in the atmosphere (unit: kg N.ha-1)
                - 'cumNRes1': List of daily total mineral N direcly released by residue mineralisation (unit: kg N.ha-1)
                - 'cumNRes2': List of daily total mineral N uptake for microbial biomass growth (unit: kg N.ha-1)
                - 'cumNRes3': List of daily total mineral N released from microbial biomass turnover (unit: kg N.ha-1)
                - 'NminfromNresCum': List of daily residue mineralisation by residue type (unit: kg N.ha-1)


        """
        #""" Dictionnary for soil Organic and Mineral balance (kg N.ha-1)
        #Keys for Daily outputs: 
        #Keys for Total Input: 
        #Keys for Total outputs: 
        #Keys for totals: 'InputNtot', 'OutputNtot', 'InputNmintot', 'OutputNmintot'
        #"""
        
        self.bilanN = {}
        self.bilanN['intialInertN'] = sum3(self.InertNorg) / self.soilSurface *10000
        self.bilanN['intialActiveN'] = sum3(self.Norg-self.InertNorg) / self.soilSurface *10000
        self.bilanN['cumMinN'] = []
        self.bilanN['intialNO3'] = sum3(self.m_NO3) / self.soilSurface *10000
        self.bilanN['intialNH4'] = sum3(self.m_NH4) / self.soilSurface *10000
        self.bilanN['cumLix'], self.bilanN['cumN2O'] = [],[]
        self.bilanN['cumRain'] = []
        self.bilanN['cumIrrig'] = []
        self.bilanN['cumfertNO3'] = []
        self.bilanN['cumfertNH4'] = []
        self.bilanN['cumUptakePlt'] = []
        self.bilanN['azomes'] = [] #equivalent a stics
        #residus
        self.bilanN['cumNRes1'], self.bilanN['cumNRes2'], self.bilanN['cumNRes3'] = [], [], []
        self.bilanN['NminfromNresCum'] = []

        try:
            self.bilanN['initialNZygo'] = sum(self.ls_NBio())/ self.soilSurface *10000
            self.bilanN['initialNres'] = sum(self.ls_NRES())/ self.soilSurface *10000
        except:#si pas de residus ajoute
            self.bilanN['initialNZygo'] = 0.
            self.bilanN['initialNres'] = 0.

        #calcul des Nzigo et Nres prevus nulle part: seulement C qui etait compte (a developper!)


    #def UpdateNorgbalance(self, dNMin):
    #    surfsolref = self.soilSurface
    #    self.bilanN['cumMinN'].append(sum3(dNMin) / surfsolref *10000)
    #    # autres termes avec microbio et residus a ajouter!!

    #distribue dans differentes fonctions
    #def UpdateNminbalance(self, Lix, dN2O):
    #    surfsolref = self.soilSurface
    #    self.bilanN['cumLix'].append(Lix / surfsolref *10000)
    #    self.bilanN['cumN2O'].append(dN2O/ surfsolref *10000)
    #    # autres termes avec plantes...

    def CloseNbalance(self, print_=1):
        """ Finalise calculation of whole simulation variables in bilanC attribute
        
        """
        
        surfsolref = self.soilSurface
        #N org humus
        self.bilanN['FinalInertN'] = sum3(self.InertNorg) / surfsolref *10000
        self.bilanN['FinalActiveN'] = sum3(self.Norg-self.InertNorg) / surfsolref *10000


        #residus
        try:
            self.bilanN['finalNZygo'] = sum(self.ls_NBio())/ surfsolref *10000
            self.bilanN['finalNres'] = sum(self.ls_NRES())/ surfsolref *10000
        except:#si pas de residus ajoute
            self.bilanN['finalNZygo'] = 0.
            self.bilanN['finalNres'] = 0.

        self.bilanN['ResidueMinNtot'] = sum(self.bilanN['cumNRes1']) + sum(self.bilanN['cumNRes2']) + sum(self.bilanN['cumNRes3'])
        try:
            self.bilanN['NminfromNresCum'] = list(map(sum, self.bilanN['NminfromNres']))
        except:
            self.bilanN['NminfromNresCum'] = 0.

        self.bilanN['HumusMinNtot'] = sum(self.bilanN['cumMinN'])
        self.bilanN['MinNtot'] = self.bilanN['ResidueMinNtot'] + self.bilanN['HumusMinNtot'] 
        self.bilanN['InputNtot'] = self.bilanN['intialInertN'] + self.bilanN['intialActiveN'] + self.bilanN['initialNres'] + self.bilanN['initialNZygo']
        self.bilanN['OutputNtot'] = self.bilanN['FinalInertN'] + self.bilanN['FinalActiveN'] + self.bilanN['finalNres'] + self.bilanN['finalNZygo'] + self.bilanN['MinNtot'] 


        #Input Min
        self.bilanN['TotNRain'] = sum(self.bilanN['cumRain'])
        self.bilanN['TotNIrrig'] = sum(self.bilanN['cumIrrig'])
        self.bilanN['TotFertNO3'] = sum(self.bilanN['cumfertNO3'])
        self.bilanN['TotFertNH4'] = sum(self.bilanN['cumfertNH4'])

        self.bilanN['InputNmintot'] = self.bilanN['intialNO3'] + self.bilanN['intialNH4'] + self.bilanN['MinNtot'] + self.bilanN['TotNRain'] + self.bilanN['TotNIrrig'] + self.bilanN['TotFertNO3'] + self.bilanN['TotFertNH4']
        
        #Output Min
        self.bilanN['FinalNO3'] = sum3(self.m_NO3) / surfsolref *10000
        self.bilanN['FinalNH4'] = sum3(self.m_NH4) / surfsolref *10000
        self.bilanN['Lixtot'] = sum(self.bilanN['cumLix'])
        self.bilanN['N2Otot'] = sum(self.bilanN['cumN2O']) #!! manque denitrif!
        self.bilanN['TotUptPlt'] = sum(self.bilanN['cumUptakePlt'])
        self.bilanN['OutputNmintot'] = self.bilanN['FinalNO3'] + self.bilanN['FinalNH4'] + self.bilanN['Lixtot'] + self.bilanN['N2Otot'] + self.bilanN['TotUptPlt']

        if print_==1:
            self.PrintNbalance()


    def PrintNbalance(self):
        """ Print a summary table of bilanN attribute
        
        """
        
        bilanN= self.bilanN
        #Norg
        print ("")
        print ("Organic N Balance Input (kg N.ha-1)\t\t\t Organic N Balance Output (kg N.ha-1)")
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Active Humified Pool:\t {0:8.1f}\t\t\t Active Humified Pool:\t {1:8.1f}".format(bilanN['intialActiveN'], bilanN['FinalActiveN'])))
        print(("Inert Humified Pool:\t {0:8.1f}\t\t\t Inert Humified Pool:\t {1:8.1f}".format(bilanN['intialInertN'], bilanN['FinalInertN'])))
        print(("Zymogeneous Bio Pool:\t {0:8.1f}\t\t\t Zymogeneous Bio Pool:\t {1:8.1f}".format(bilanN['initialNZygo'], bilanN['finalNZygo'])))
        print(("Added organic matter:\t {0:8.1f}\t\t\t Added organic matter:\t {1:8.1f}".format(bilanN['initialNres'], bilanN['finalNres'])))
        print(("                            \t\t\t\t Mineralisation:\t\t {0:8.1f}".format(bilanN['MinNtot'])))
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Total:\t\t\t\t\t {0:8.1f}\t\t\t Total:\t\t\t\t\t {1:8.1f}".format(bilanN['InputNtot'], bilanN['OutputNtot'])))
        print ("")

        #Nmin
        print ("")
        print ("Mineral N Balance Input (kg N.ha-1)\t\t\t Mineral N Balance Output (kg N.ha-1)")
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Initial soil NO3 :\t\t {0:8.1f}\t\t\t Final soil NO3:\t {1:8.1f}".format(bilanN['intialNO3'], bilanN['FinalNO3'])))
        print(("Initial soil NH4:\t\t {0:8.1f}\t\t\t Final soil NH4:\t {1:8.1f}".format(bilanN['intialNH4'], bilanN['FinalNH4'])))
        print(("Humus Mineralisation:\t {0:8.1f}\t\t\t Leaching:\t\t\t {1:8.1f}".format(bilanN['HumusMinNtot'], bilanN['Lixtot'])))
        print(("Resid. Mineralisation:\t {0:8.1f}\t\t\t N2O:\t\t\t\t {1:8.1f}".format(bilanN['ResidueMinNtot'], bilanN['N2Otot'])))
        print(("Rain:\t\t\t\t\t {0:8.1f}\t\t\t Uptake plant:\t\t {1:8.1f}".format(bilanN['TotNRain'], bilanN['TotUptPlt'])))
        print(("Irrigation:\t\t\t\t {0:8.1f}\t\t\t ".format(bilanN['TotNIrrig'])))
        print(("Fertilizers NO3:\t\t {0:8.1f}\t\t\t ".format(bilanN['TotFertNO3'])))
        print(("Fertilizers NH4:\t\t {0:8.1f}\t\t\t ".format(bilanN['TotFertNH4'])))


        #print ("Rain:\t\t\t\t\t {0:8.1f}\t\t\t ".format(bilanN['TotNRain']))
        print ("----------------------------\t\t\t\t ----------------------------")
        print(("Total:\t\t\t\t\t {0:8.1f}\t\t\t Total:\t\t\t\t {1:8.1f}".format(bilanN['InputNmintot'], bilanN['OutputNmintot'])))
        print ("")







def step_bilanWN_solVGL(S, par_SN, meteo_j,  mng_j, ParamP, ls_epsi, ls_roots, ls_demandeN_bis, opt_residu, opt_Nuptake):
    """ Daily step for soil Water and Nitrogen balance from L-egume lsystem inputs, meteo and management

    :param S: Previous day SoilN object
    :type S: `soil3ds.soil_moduleN.SoilN`
    :param par_SN: A dictionnary defining general soil parameters derived from STICS
    :type par_SN: dict
    :param meteo_j: A dictionnary containing daily meterorological data
    :type meteo_j: dict
    :param mng_j: A dictionnary containing daily management of irrigation and N fertilisation
    :type mng_j: dict
    :param ParamP: List of individual plant parameters for the L-egume model
    :type ParamP: list
    :param ls_epsi: List of individual plant daily fraction of incoming global radiation interception (0-1 fraction)
    :type ls_epsi: list
    :param ls_roots: List of [z,x,y] arrays for root length distribution per individual plant (unit: m)
    :type ls_roots: list
    :param ls_demandeN_bis: List of individual plant Nitrogen demand
    :type ls_demandeN_bis: list
    :param opt_residu: Option to activate calculation of residue mineralisation (0=deactivate or 1=activate)
    :type opt_residu: int
    :param opt_Nuptake: Option to define the calculation of plant N uptake (0:'STICS', 1:'LocalTransporter', 2:'old')
    :type opt_Nuptake: int

    :return:
        * [S,  new_stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, temps_sol]
            - S (`soil3ds.soil_moduleN.SoilN`): Updated SoilN object
            - new_stateEV (list): Updated memory variables to compute soil transpiration
            - ls_ftsw (list): List of daily ftsw by individual plants (0-1 fraction)
            - ls_transp (list): List of daily transpiration by individual plants (unit: mm)
            - ls_Act_Nuptake_plt (list):  List of Array of size [nz,nx,ny] storing voxel daily individual plant mineral N uptake (unit: kg N per voxel.plant-1.day-1)
            - temps_sol (list): List of other output variables [evapo_tot, Drainage, ls_m_transpi, m_evap, ActUpNtot, ls_DQ_N, idmin]

    """
    # testRL = updateRootDistrib(invar['RLTot'][0], ls_systrac[0], lims_sol)
    # ls_roots = rtd.build_ls_roots_mult(invar['RLTot'], ls_systrac, lims_sol) #ancien calcul base sur SRL fixe
    #ls_roots = rtd.build_ls_roots_mult(array(invar['RLTotNet']) * 100. + 10e-15, ls_systrac, lims_sol)  # !*100 pour passer en cm et tester absoption d'azote (normalement m) #a passer apres calcul de longuer de racine!
    # esternalise calcul de ls_roots -> wrapper prend grille en entree et plus geom

    surfsolref = S.surfsolref

    # preparation des entrees eau
    Rain = meteo_j['Precip']
    Irrig = mng_j['Irrig']  # ['irrig_Rh1N']#R1N = sol_nu

    # preparation des entrees azote
    mapN_Rain = 1. * S.m_1[0, :, :] * Rain * par_SN['concrr'] * S.m_vox_surf[0,:,:] # Nmin de la pluie kg N par voxel
    mapN_Irrig = 1. * S.m_1[0, :, :] * Irrig * par_SN['concrr'] * S.m_vox_surf[0,:,:] # Nmin de l'eau d'irrigation kg N par voxel
    mapN_fertNO3 = 1. * S.m_1[0, :, :] * mng_j['FertNO3'] * S.m_vox_surf[0, :, :] / 10000.  # kg N par voxel
    mapN_fertNH4 = 1. * S.m_1[0, :, :] * mng_j['FertNH4'] * S.m_vox_surf[0, :, :] / 10000.  # kg N par voxel

    S.updateTsol(meteo_j['Tsol'])  # (meteo_j['TmoyDay'])#(meteo_j['Tsol'])# #Tsol forcee comme dans STICS (Tsol lue!!)

    #############
    # step  sol
    #############
    treshEffRoots_ = 10e10  # valeur pour forcer a prendre densite effective
    stateEV, Uval, b_ = S.stateEV, S.Uval, S.b_
    ls_transp, evapo_tot, Drainage, new_stateEV, ls_m_transpi, m_evap, ls_ftsw = S.stepWBmc(meteo_j['Et0'] * surfsolref,
                                                                                        ls_roots, ls_epsi,
                                                                                        Rain * surfsolref,
                                                                                        Irrig * surfsolref, stateEV,
                                                                                        par_SN['ZESX'], leafAlbedo=0.15,
                                                                                        U=Uval, b=b_, FTSWThreshold=0.4,
                                                                                        treshEffRoots=treshEffRoots_,
                                                                                        opt=1)
    S.update_memory_EV(new_stateEV)
    S.stepNB(par_SN)
    if opt_residu == 1:  # s'ily a des residus
        S.stepResidueMin(par_SN)
        S.stepMicrobioMin(par_SN)
    S.stepNitrif(par_SN)
    #ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N, idmin = S.stepNuptakePlt(par_SN, ParamP, ls_roots, ls_m_transpi,ls_demandeN_bis)
    S.stepNINFILT(mapN_Rain, mapN_Irrig, mapN_fertNO3, mapN_fertNH4, Drainage, opt_infil=1)
    ActUpNtot, ls_Act_Nuptake_plt, ls_DQ_N, idmin = S.stepNuptakePlt(par_SN, ParamP, ls_roots, ls_m_transpi, ls_demandeN_bis, opt_Nuptake)
    #print(amin(S.m_NO3), unique(idmin, return_counts=True),ls_DQ_N)

    temps_sol = [evapo_tot, Drainage, ls_m_transpi, m_evap, ActUpNtot, ls_DQ_N, idmin] #other output variables

    return [S,  new_stateEV, ls_ftsw, ls_transp, ls_Act_Nuptake_plt, temps_sol]
    #lims_sol et surfsolref pourrait pas etre fournie via S.?
    #pourquoi b_ et Uval trainent la? (paramtres sol??)
    #return more output variables?? -> OK temps_sol



def default_parSN():
    """ Creates a default parameter dictionnary 'parSN' for defining the STICS parameters necessary for computing N balance

        Keys of the dictionnary are  parameters for a given soil object:

        General STICS Nitrogen parameters
            * 'FMIN1G' :       (day-1) (p145)
            * 'FMIN2G' :       (%ARGIS-1) para pot rate min a ARGIs (p145)
            * 'FMIN3G' :       (% CALC-1)para pot rate min a CALCss (p145)
            * 'FINERTG' :      0.65 = default fraction of N pool inactive for minearlisation (p145) (could be smaller in grassland & forest)
            * 'HMinMg' :       Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
            * 'HoptMg' :       Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)
            * 'TRefg' :        reference temperature (degreC)
            * 'FTEMHAg' :      asymptotic value of FTH (seen as a logistic response)
            * 'FTEMHg' :       (K-1)
            * 'FTEMHB' :
            * 'FNXg' :         (day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
            * 'PHMinNITg' :    pH min de nitrification (prop of Field Capacity) #value p149
            * 'PHMaxNITg' :    pH max de nitrification (prop of Field Capacity) #value p149
            * 'HMinNg' :       Humidite min de nitrification #value p149
            * 'HoptNg' :       Humidite opt de nitrification  #value p149
            * 'TNITMINg' :     Temperature min de nitrification #  (degreC)value p151
            * 'TNITOPTg' :     Temperature opt de nitrification #  (degreC)value p151
            * 'TNITMAXg' :     Temperature max de nitrification #  (degreC)value p151
            * 'RATIONITs' :    proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
            * 'DIFNg' :        N diffusion coefficient at field capacity (cm2.day-1, p 161)

        soil specific water and nitrogen parameters
            * 'ZESX' :          Maximal depth of soil affected by evaporation (unit: m)
            * 'CFES' :          Shape parameter to define relative contribution of soil depth to evaporation (unitless; 1. = proportional to soil depth)
            * 'q0' :            Cumulative evaporation treshold indicating the end of maximum evaporation rate (unit: mm)

            * 'Norg' :          Soil Organic N content in the topsoil layer (unit: g N.kg-1 soil)
            * 'PROFHUMs' :      Depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220 (unit: cm)

            * 'pH' :            Soil pH (unit: pH unit)
            * 'ARGIs' :     Soil clay content (unit: %)
            * 'ACLIMc' :        Climatic demand of soil evaporation parameter (unit: mm)
            * 'concrr' :        nitrate concentration in rainfall (unit: kg N.m-2.mm-1)

    :return: Default 'parSN' parameter dictionnary

    .. code-block:: python

        parSN = {
                    'FMIN1G' : 0.0006,
                    'FMIN2G' : 0.0272,
                    'FMIN3G' : 0.0167,
                    'FINERTG' : 0.65,
                    'HMinMg' : 0.3,
                    'HoptMg' : 1.,
                    'TRefg' :  15. ,
                    'FTEMHAg' : 25. ,
                    'FTEMHg' :  0.12 ,
                    'FTEMHB' : 145.,
                    'FNXg' :   0.5,
                    'PHMinNITg' :  3. ,
                    'PHMaxNITg' : 5.5 ,
                    'HMinNg' : 0.67,
                    'HoptNg' : 1.,
                    'TNITMINg' : 5.,
                    'TNITOPTg' : 20.,
                    'TNITMAXg' :  45.,
                    'RATIONITs' : 0.,
                    'DIFNg' :  0.018,
                    'ZESX' : 0.30,
                    'CFES' : 1.,
                    'q0' :  9.42,
                    'Norg' : 1.1,
                    'PROFHUMs' : 30.,
                    'pH' : 7.1,
                    'ARGIs': 18.3,
                    'ACLIMc' : 26.,
                    'concrr' : 0.000002
                    }

    """
    # parameters 'WCST' and 'gamma_theo' are not used / necessary in current model version

    # General STICS Nitrogen parameters
    par_SN = {}
    par_SN['FMIN1G'] = 0.0006  # (day-1) (p145)
    par_SN['FMIN2G'] = 0.0272  # (%ARGIS-1) para pot rate min a ARGIs (p145)
    par_SN['FMIN3G'] = 0.0167  # (% CALC-1)para pot rate min a CALCss (p145)

    par_SN['FINERTG'] = 0.65  # 0.65 = default fraction of stable pool (p145) should be smaller in grassland & forest

    par_SN['HMinMg'] = 0.3  # Humidite min de mineralisation (prop of Field Capacity) #value p 142 (cite Rodrigo et al. 1997)
    par_SN['HoptMg'] = 1.  # Humidite opt de mineralisation (prop of Field capacity) #value p 142 (cite Rodrigo et al. 1997)

    par_SN['TRefg'] = 15.  # reference temperature
    par_SN['FTEMHAg'] = 25.  # asymptotic value of FTH (seen as a logistic response)
    par_SN['FTEMHg'] = 0.120  # (K-1)
    par_SN['FTEMHB'] = 145.

    par_SN['FNXg'] = 0.5  # (day-1) maximum fraction of NH4 transformed by nitrification every day in the nitrification layer(default value in STICS parameter excel file; 0.5 p149 for tropical soils)
    par_SN['PHMinNITg'] = 3.  # pH min de nitrification (prop of Field Capacity) #value p149
    par_SN['PHMaxNITg'] = 5.5  # pH max de nitrification (prop of Field Capacity) #value p149
    par_SN['HMinNg'] = 0.67  # Humidite min de nitrification #value p149
    par_SN['HoptNg'] = 1.  # Humidite opt de nitrification  #value p149
    par_SN['TNITMINg'] = 5.  # Temperature min de nitrification #  (degreC)value p151
    par_SN['TNITOPTg'] = 20.  # Temperature opt de nitrification #  (degreC)value p151
    par_SN['TNITMAXg'] = 45.  # Temperature max de nitrification #  (degreC)value p151
    par_SN['RATIONITs'] = 0.  # proportion of nitrtified NH4 converted to N2O (pas trouve de valeur par defaut - p151) #0-> N2O pas active
    par_SN['DIFNg'] = 0.018  # N diffusion coefficient at field capacity (cm2.day-1, p 161)

    # soil specific water and nitrogen balance parameters
    par_SN['ZESX'] = 0.30  #  (m)
    par_SN['CFES'] = 1.  #  (unitless)
    par_SN['q0'] = 9.46  #  (mm)

    par_SN['Norg'] = 1.1  # (g N.kg-1 sol)
    par_SN['PROFHUMs'] = 30.  # (cm) depth of soil contributing to SOM mineralisation #-> peut constituer un masque pour considerer certaines couches ou pas default value p220

    par_SN['pH'] = 7.1  # (pH unit)
    par_SN['ARGIs'] = 18.3  #
    par_SN['ACLIMc'] = 26.  #
    par_SN['concrr'] = 0.000002  #concentration en N de la pluie (kg N.m-2.mm-1 )

    return par_SN
    # ajouter parametre hydriques generaux?




def default_paramp():
    """ Creates a default parameter dictionnary 'paramp' for defining the STICS plant parameters necessary for computing N balance

        Keys of the dictionnary are  parameters for a given soil object:

        General STICS Nitrogen uptake parameters
            * 'Vmax1' : HATS Vmax (unit: micromole.cm-1.h-1)
            * 'Kmax1' : HATS Kmax (unit: micromole.L-1)
            * 'Vmax2' : LATS Vmax (unit: micromole.cm-1.h-1)
            * 'Kmax2' : LATS Kmax (unit: micromole.L-1)
            * 'treshmaxN' : Lower treshold of plant N status feedback effect on root nitrogen uptake (unit: either in NNI unit or %)
            * 'treshminN' : Higher treshold of plant N status feedback effect on root nitrogen uptake (unit: either in NNI unit or %)
            * 'leafAlbedo' : Leaf albedo (unit: 0-1 fraction)
            * 'WaterTreshGs' : FTSW treshold for the onset of transpiration reduction (unit: 0-1 fraction)
            * 'treshEffRootsN' : Treshold for maximal total root length density per voxel which is effective for N uptake (cm.cm-3)

    :return: Default 'paramp' parameter dictionnary

    .. code-block:: python

        paramp =  {
                    'Vmax1': 0.0018,
                     'Kmax1': 50.0,
                     'Vmax2': 0.05,
                     'Kmax2': 25000.0,
                     'treshmaxN': 1.0,
                     'treshminN': 0.8,
                     'leafAlbedo': 0.15,
                     'WaterTreshGs': 0.4,
                     'treshEffRootsN': 10**6

                    }

    """

    paramp = {}
    paramp['Vmax1'] = 0.0018
    paramp['Kmax1'] = 50.
    paramp['Vmax2'] = 0.05
    paramp['Kmax2'] = 25000.0
    paramp['treshmaxN'] = 1.0
    paramp['treshminN'] = 0.8

    paramp['leafAlbedo'] = 0.15
    paramp['WaterTreshGs'] = 0.4
    paramp['treshEffRootsN'] = 10**6

    return paramp

