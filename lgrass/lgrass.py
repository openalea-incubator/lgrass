# coding: utf8

# script contenant les classes de lgrass.lpy, pour que le module pickle puisse les lire

import numpy as np
#########   Definition Plante   #########


class ParamPlante:
    def __init__(self, id_plante=0, id_geno=0):
        self.id_plante = id_plante
        self.id_geno = id_geno


class ParamPhytomere:
    def __init__(self, id_plante=0, id_talle=0, id_rang=1, topology=(0,), age=0, angletal=0,
                 organ_lengths={}):  # , Carbdescend=0, Carbmonte=0):
        self.id_plante = id_plante
        self.id_talle = id_talle
        self.id_rang = id_rang
        self.topology = topology  # Topological information
        self.age = age
        self.angletal = angletal
        self.organ_lengths = organ_lengths


class ParamFeuille:
    def __init__(self, age=0, Agecroiss=0, Taillefeuille=0, id_plante=0, id_talle=0, id_rang=1, topology=(0,), Ymax=0,
                 Difftps=0, Taillefinalelimbe=0, Taillefinalegaine=0., Taillelimbe=0, Taillegaine=0, Phase='cachee',
                 rapportK=0.2, coupe=0., Cutstatus='intact', angleinsert=0, angletal=0., id_geno=0, surface_limbe=0,
                 surface_gaine=0, biomass=0., Besoinencroiss=0, TailleEmergence=0., R=0.973):
        self.age = age  # Age de la feuille
        self.Agecroiss = Agecroiss  #
        self.Taillefeuille = Taillefeuille  # Taille de la feuille au temps t
        self.id_plante = id_plante  # Identifiant de la plante (compris entre 0 et nb_plantes)
        self.id_talle = id_talle  # Identifiant de la talle
        self.id_rang = id_rang  # Rang de la feuille
        self.topology = topology  # Topological information
        self.Ymax = Ymax  # Taille finale de la feuille
        self.Difftps = Difftps  # Delai entre l'apparition de 2 feuilles
        self.Taillefinalelimbe = Taillefinalelimbe  # Taille finale du limbe
        self.Taillefinalegaine = Taillefinalegaine  # Taille finale de la gaine
        self.Taillelimbe = Taillelimbe  # Taille du limbe au temps t
        self.Taillegaine = Taillegaine  # Taille de la gaine au temps t
        self.Phase = Phase  # Stade de la feuille : cachee, visible, 0 0
        self.rapportK = rapportK  # Rapport gaine/limbe
        self.coupe = coupe  # Longueur de coupe
        self.Cutstatus = Cutstatus  # Quantite de coupe : intact, partiellementcoupee(=feuille coupee et gaine intacte), entierementcoupee(=feuille et gaine coupee)
        self.angleinsert = angleinsert  # Angle limbe
        self.angletal = angletal  # Angle de la talle qui ecarte la feuille
        self.id_geno = id_geno
        self.surface_limbe = surface_limbe  # Surface limbe (mm2)
        self.surface_gaine = surface_gaine  # Surface gaine (mm2)
        self.biomass = biomass
        self.Besoinencroiss = Besoinencroiss
        self.TailleEmergence = TailleEmergence
        self.R = R


class ParamSegFeuille:
    def __init__(self, id_plante, id_talle, id_rang, idLong, hauteur):
        self.id_plante = id_plante
        self.id_talle = id_talle
        self.id_rang = id_rang
        self.idLong = idLong
        self.hauteur = hauteur


class ParamEntrenoeud:
    def __init__(self, id_plante, id_talle, id_rang=1, topology=(0,), internode_type='short_internode', length=0,
                 width=2, final_length=10):
        self.id_plante = id_plante  # Id of the plant (between 0 and number of plants-1)
        self.id_talle = id_talle  # Id of the tiller
        self.id_rang = id_rang  # Rank of the phytomere
        self.topology = topology  # Topological information
        self.internode_type = internode_type  # Type of the internode mean : short_internode or long_internode
        self.length = length  # Current length of the internode when time = t
        self.width = width  # Width of the internode, constant
        self.final_length = final_length  # Final length of the internode (function of the sheath of the previous phytomere)


class ParamApex:
    def __init__(self, id_plante, id_talle, id_rang=1, topology=(0,), retard=0, total_emerged_leaves=0,
                 total_created_leaves=0, day=0, primary_induction_index=0, secondary_induction_index=0,
                 phenological_state='vegetative', primordia_number=3, spikelet_created=0, final_leaf_number=None,
                 final_spikelet_number=None, organ_lengths_dict={}, terminal_spikelet_height=0, heading_date=None,
                 appearance_date=None, id_geno=0):  # PARAMETRE
        self.id_plante = id_plante
        self.id_talle = id_talle
        self.id_rang = id_rang
        self.topology = topology  # Topological information
        self.retard = retard
        self.total_emerged_leaves = total_emerged_leaves
        self.total_created_leaves = total_created_leaves  # number of leaves produced by the apex
        self.day = day
        self.primary_induction_index = primary_induction_index
        self.secondary_induction_index = secondary_induction_index
        self.phenological_state = phenological_state  # Phenological state of the apex. Values : vegetative, reproductive, dead
        self.primordia_number = primordia_number
        self.spikelet_created = spikelet_created
        self.final_leaf_number = final_leaf_number
        self.final_spikelet_number = final_spikelet_number
        self.organ_lengths_dict = organ_lengths_dict
        self.terminal_spikelet_height = terminal_spikelet_height
        self.heading_date = heading_date
        self.appearance_date = appearance_date
        self.id_geno = id_geno


class ParamApexTal:
    def __init__(self, id_plante, id_talle=0, id_rang=1, topology=(0,), retard=0, day=0, primary_induction_index=0,
                 secondary_induction_index=0, phenological_state='vegetative', id_geno=0):
        self.id_plante = id_plante
        self.id_talle = id_talle
        self.id_rang = id_rang
        self.topology = topology  # Topological information
        self.retard = retard
        self.day = day
        self.primary_induction_index = primary_induction_index
        self.secondary_induction_index = secondary_induction_index
        self.phenological_state = phenological_state  # Phenological state of the apex. Values : vegetative, reproductive, dead
        self.id_geno = id_geno


class ParamEntrenoeudEpi:
    def __init__(self, id_plante, id_talle, id_rang, topology=(0,), long_entrenoeud=0):
        self.id_plante = id_plante
        self.id_talle = id_talle
        self.id_rang = id_rang
        self.topology = topology  # Topological information
        self.long_entrenoeud = long_entrenoeud  # Longueur moyenne des entrenoeuds (entre chaque insertion d epillet)


class ParamEpillet:
    def __init__(self, id_plante, id_talle, id_rang, topology=(0,), long_epillet=0):
        self.id_plante = id_plante
        self.id_talle = id_talle
        self.id_rang = id_rang
        self.topology = topology  # Topological information
        self.long_epillet = long_epillet  # Longueur moyenne des epillets


class ParambourgeonRoot:
    def __init__(self, id_plante, id_talle=0, id_rang=1, nb_prod_root=0, id_geno=0):
        self.id_plante = id_plante
        self.id_talle = id_talle
        self.id_rang = id_rang
        self.nb_prod_root = nb_prod_root
        self.id_geno = id_geno


class ParamRacine:
    def __init__(self, age=0, id_plante=0, id_talle=0, id_rang=0, axe_arret=0, lateral_arret=0, feuil_ref=0):
        self.age = age
        self.id_plante = id_plante
        self.id_talle = id_talle
        self.id_rang = id_rang
        self.axe_arret = axe_arret
        self.lateral_arret = lateral_arret
        self.feuil_ref = feuil_ref


class pte:
    def __init__(self, id_pointe=0, id_plante=0, age=0, diametre=0., distPrimInit=0., longueur=0.01, profondeur=0.,
                 dateDerniereCreation=0, posO=np.array([0., 0., 0.]),
                 Tortue=np.array([[0., 0., -1.], [0., 1., 0], [-1., 0., 0.]]), isaxe=False, arretee=False, senile=False,
                 axe_mort=False, segment=0):
        self.id_pointe = id_pointe
        self.id_plante = id_plante
        self.age = age
        self.diametre = diametre
        self.distPrimInit = distPrimInit
        self.longueur = longueur
        self.profondeur = profondeur
        self.dateDerniereCreation = dateDerniereCreation
        self.posO = posO  # Position de la pointe dans l'espace
        self.Tortue = Tortue
        self.isaxe = isaxe
        self.arretee = arretee  #
        self.senile = senile  # Etat de senescence de la pointe
        self.axe_mort = axe_mort
        self.segment = segment


class primord:
    def __init__(self, id_plante, id_primord, age, diametre, Tortue=np.array([[0., 0., -1], [0., 1, 0], [-1., 0., 0.]]),
                 avorte=False, id_pointe_axe=0, posO=np.array([0., 0., 0.])):
        self.id_plante = id_plante
        self.id_primord = id_primord
        self.age = age
        self.diametre = diametre
        self.Tortue = Tortue
        self.avorte = avorte
        self.id_pointe_axe = id_pointe_axe
        self.posO = posO


class seg:
    def __init__(self, jourForm, diametre, longueur=0., id_pointe_axe=0):
        self.jourForm = jourForm
        self.diametre = diametre
        self.longueur = longueur
        self.id_pointe_axe = id_pointe_axe


class Horizon:  # Horizon de sol
    def __init__(self, croiss, ramif, iCMeca, oCMeca):
        self.croiss  # Coefficient de croissance, compris entre 0 et 1
        self.ramif  # Coefficient multiplicateur de distance inter-ramif
        self.iCMeca  # Intensite de la contrainte mecanique
        self.oCMeca  # Orientation de la contrainte mecanique (O iso, ou 1 vert)
