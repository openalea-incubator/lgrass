import os
from copy import deepcopy

import numpy as np
import pandas as pd
import scipy

import openalea.lpy as lpy

from lgrass import meteo_ephem
from lgrass import param_reproduction_functions as prf
from lgrass import cuts
from lgrass import run_caribu_lgrass
from lgrass import gen_lstring
import lgrass

from plantfusion.utils import create_child_folder
from plantfusion.light_wrapper import Light_wrapper
from plantfusion.indexer import Indexer


class L_grass_wrapper():
    """
    Wrapper for l-grass model

    mandatory methods : init, derive, light_inputs, light_results, soil_inputs, soil_results, run
    optional method : end
    other methods can be added for readability

    """

    def __init__(
            self,
            name="lgrass",
            indexer=Indexer(),
            in_folder="", out_folder=None,
            nameconfigfile="liste_usms_exemple.xls", ongletconfigfile="exemple", # for soil ?
            IDusm=None,
            planter=None,
            caribu_scene=False
            ):
        """
        Initialize FSPM

        - gestion des path (inputs et outputs)
        - gestion du plan (via Planter & plan_simulation.csv)
        - définition des paramètres des plantes (avec les fichiers de lgrass)
        - définition des paramètres de simulation
        - Initialisation du lsystem (via self.load_lsystem())

        """
        # Manage output path



        # Manage Planter
        self.planter = planter
        self.name = name
        self.load_lsystem(lpy_filename='lgrass.lpy')





    def load_lsystem(self, lpy_filename, id_scenario=0, id_gener=1):
        """
        Initialize Lsystem
        lpy_filename : Lsystem file 
        plan_sim : pandas Dataframe (see plan_simulation.csv)
        id_scenario : row index of plan_sim to use
        
        """
                # Fichiers d'entrée
        genet_file = 'ped.r'
        # un fichier de remplacement du modèle génétique qui génère une population de C et détermine le nombre de plantes du couvert
        param_plant_file = 'liste_plantes.csv'

        # Répertoires de lecture/écriture
        INPUTS_DIRPATH = 'inputs'
        OUTPUTS_DIRPATH = 'outputs'
        GENET_DIRPATH = 'modelgenet'

        plan_sim = pd.read_csv('plan_simulation.csv')
        row = plan_sim.iloc[id_scenario]

        self.lsystem = lpy.Lsystem(lpy_filename)
        

        # Choix du fichier de lecture du C en fonction de l'option de reproduction des plantes
        opt_repro = row["option_reproduction"]
        if opt_repro != "False":
            self.in_genet_file = os.path.join(GENET_DIRPATH, genet_file)
        else:
            self.in_genet_file = None
        self.in_param_file = os.path.join(INPUTS_DIRPATH, param_plant_file)


        self.lsystem.name_sim = self.name

        self.lsystem.option_tallage = row["option_tallage"]
        self.lsystem.option_senescence = row["option_senescence"]
        self.lsystem.option_floraison = row["option_floraison"]
        self.lsystem.option_tiller_regression = row["option_tiller_regression"]
        self.lsystem.option_morphogenetic_regulation_by_carbone = row["option_morphogenetic_regulation_by_carbone"]
        self.lsystem.derivationLength = int(row["derivationLength"])
        self.lsystem.sowing_date = row["sowing_date"]
        self.lsystem.site = row["site"]
        self.lsystem.meteo = meteo_ephem.import_meteo_data(row["meteo_path"], row['sowing_date'], row['site'])
        self.lsystem.option_caribu = row["option_caribu"]
        #self.lsystem.output_induction_file_name = name + '_' + 'induction'
        #self.lsystem.output_organ_lengths_file_name = name + '_' + 'organ_lengths'

        # Parametres des plantes
        # ATTENTION : le PLANTER doit overwrite les paramètres nb_plantes, NBlignes, NBcolonnes, posplante 
        self.lsystem.ParamP, self.lsystem.nb_plantes, self.lsystem.NBlignes, self.lsystem.NBcolonnes, self.lsystem.posPlante, self.lsystem.Plantes, self.lsystem.Genotypes, self.lsystem.flowering_model = prf.define_param(
            in_param_file=self.in_param_file, in_genet_file=self.in_genet_file, out_param_file=os.path.join(OUTPUTS_DIRPATH, self.name + '.csv'),
            id_gener=id_gener, opt_repro=opt_repro)
        
        # overwrite posPlantes with planter data
        if self.planter is not None:
            self.indice_instance = self.planter.indexer.other_names.index(self.name) # get the index of the instance in the other_names list (other than legume or wheat)
            self.plant_positions = self.planter.generate_random_other(indice_instance=self.indice_instance) # generate positions depending on instance plant density and domain.
            self.lsystem.nb_plantes = len(self.plant_positions)

        self.lsystem.lstring = self.lsystem.axiom

        #init caribu

        if self.lsystem.option_caribu :
            self.dico_caribu = run_caribu_lgrass.init(meteo=self.lsystem.meteo,
                                                       nb_plantes=self.lsystem.nb_plantes,
                                                       scenario=row) 


        def derive(self, dd=1):
            """
            Derive FSPM : update l-string for scene interpretation

            /!\ degree-day-wise derivation

            NB : if combined with a one-day timestep FSPM, a loop checking for date must be implemented in main file
            For combination of 1-day and 1-hour (l-egume & wheat), see :
            https://github.com/mwoussen/plantfusion/blob/develop/simulations/full_coupling_random.py

            """
            self.lsystem.lstring = self.lsystem.derive(self.lsystem.lstring, dd, 1)


        def light_inputs(self):
            """
            Compute scene used for running light model

            Returns a geometric scene in one of those formats :
            - Plantgl scene
            - MTG adelwheat
            - VGX file path
            - triangle dictionnary, ordered by organ
            - voxel grid
            --> WHICH ONE IS NEEDED FOR L-GRASS ?

            see : https://lightvegemanager.readthedocs.io/en/latest/inputs.html#scenes

            """
            self.lsystem.lscene = self.lsystem.sceneInterpretationtself(self.lsystem.lstring)
            return self.lsystem.lscene


        def light_results(self, lighting):
            """
            Interpretation of updated lighting object (after lighting.run()) for l-grass
            Depends on the compatible light models

            See L-grass documentation for implementation

            for inspiration : 
                https://github.com/mwoussen/plantfusion/blob/develop/plantfusion/l_egume_wrapper.py
                https://github.com/mwoussen/plantfusion/blob/develop/plantfusion/wheat_wrapper.py
            
            """
            # Caribu doit être exécuté une fois par jour

            self.lsystem.BiomProd, dico_caribu['radiation_interception'], dico_caribu['ray'] = run_caribu_lgrass.runcaribu(
                self.lsystem.lstring, 
                self.lsystel.lscene,
                self.lsystem.current_day,
                self.lsystel.tiller_appearance,
                self.lsystem.nbplantes,
                dico_caribu
            )


        def soil_inputs(self):
            """
            Get soil "scene"

            Returns 4 tables :
            - root N proportion for each plant ([0, 1])
            - root lenght (m) for each plant and each soil voxel (?)
            - variety-specific parameters
            - light interception capacity for each plant ([0, 1])

        """
            pass


        def soil_results(self):
            """
            Interpretation of updated soil object (after soil.run()) for l-grass

            """
            pass


        def run(self):
            """
            Update FSPM parameters to be ready for next derivation 
            
            """
            pass


        def end(self):
            """
            optional custom ending

            """
            pass


