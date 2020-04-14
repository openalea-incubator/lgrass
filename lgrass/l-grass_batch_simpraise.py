# coding=utf-8
# '''
# Created on 19/03/2020
#
# @author: modelisation - TR
# '''

# import the modules necessary to initiate the L-systems
import os
import pandas as pd
import openalea.lpy as opy
import openalea.plantgl as opal
from lgrass import param_reproduction_functions as prf
from lgrass import meteo_ephem
from lgrass import caribu


def runlsystem(plan_sim=None, id_scenario=0, id_gener=1):
    if plan_sim is None:
        raise NameError('Pas de plan de simulation chargé.')

    # Répertoires de lecture/écriture
    INPUTS_DIRPATH = 'inputs'
    OUTPUTS_DIRPATH = 'outputs'
    GENET_DIRPATH = 'modelgenet'

    # Charger le plan de simulation
    row = plan_sim.iloc[id_scenario]
    name = str(row["name"])

    # Charger le lsystem
    lpy_filename = os.path.join('lgrass.lpy')
    lsystem = opy.Lsystem(lpy_filename)

    # Parametres des plantes
    param = prf.get_param()
    opt_repro = row["option_reproduction"]
    if opt_repro:
        in_genet_file = os.path.join(GENET_DIRPATH, 'ped.r')     # nom de fichier en dur
    else:
        in_genet_file = os.path.join(INPUTS_DIRPATH, 'donnees_C.csv')    # nom de fichier en dur
    # lsystem.ParamP, lsystem.nb_plantes, lsystem.NBlignes, lsystem.NBcolonnes, lsystem.posPlante, lsystem.Plantes, lsystem.Genotypes, lsystem.flowering_model = reproduction.define_param()
    lsystem.ParamP, lsystem.nb_plantes, lsystem.NBlignes, lsystem.NBcolonnes, lsystem.posPlante, lsystem.Plantes, lsystem.Genotypes, lsystem.flowering_model = prf.define_param(
        in_genet_file=in_genet_file, out_param_file=os.path.join(INPUTS_DIRPATH, name+'_G'+str(id_gener)+'.csv'), param=param, id_gener=id_gener, opt_repro=opt_repro)

    # Parametres de simulation
    lsystem.option_tallage = row["option_tallage"]
    lsystem.option_senescence = row["option_senescence"]
    lsystem.option_floraison = row["option_floraison"]
    lsystem.option_tiller_regression = row["option_tiller_regression"]
    lsystem.option_mophogenetic_regulation_by_carbone = row["option_mophogenetic_regulation_by_carbone"]
    lsystem.derivationLength = int(row["derivationLength"])
    lsystem.sowing_date = row["sowing_date"]
    lsystem.site = row["site"]
    lsystem.output_induction_file_name = name + '_' + 'induction'
    lsystem.output_organ_lengths_file_name = name + '_' + 'organ_lengths'
    lsystem.meteo = meteo_ephem.import_meteo_data(row["meteo_path"], row['sowing_date'], row['site'])
    lsystem.cutting_dates = [] if pd.isna(row["cutting_dates"]) \
        else [row["cutting_dates"]] if isinstance(row["cutting_dates"], int) \
        else [int(i) for i in row["cutting_dates"].split("_")]

    # Gestion caribu
    opt_caribu = row["option_caribu"]
    if opt_caribu:
        day = 1
        dico_caribu = caribu.init(meteo=lsystem.meteo, nb_plantes=lsystem.nb_plantes, scenario=row)
        lsystem.BiomProd = [0.] * lsystem.nb_plantes

    # Rédaction d'un fichier de sortie
    path_out = os.path.join(OUTPUTS_DIRPATH, row[0] + '.csv')
    output = open(path_out, 'w')
    output.write("GDD;Date;Day;nb_talles;biomasse_aerienne;surface_foliaire" + "\n")

    # Lancement du lsystem
    lstring = lsystem.axiom
    for dd in range(lsystem.derivationLength):
        lstring = lsystem.derive(lstring, dd, 1)
        lscene = lsystem.sceneInterpretation(lstring)
        if opt_caribu:
            try:
                lsystem.BiomProd, dico_caribu['radiation_interception'], dico_caribu[
                    'Ray'] = caribu.runcaribu(lstring, lscene, lsystem.current_day,
                                              lsystem.tiller_appearance,
                                              lsystem.nb_plantes, dico_caribu, day)
            except:
                continue
            day = lsystem.current_day
            if lsystem.current_day > day:
                output.write(";".join(
                    [str(lsystem.TPS), str(lsystem.sowing_date), str(lsystem.current_day), str(lsystem.nb_talle),
                        str(lsystem.BiomProd[0]), str(lsystem.rapportS9_SSol_dict[0])]) + "\n")
        opal.all.Viewer.display(lscene)

    output.close()

    # Matrice de croisement des plantes
    mat = prf.create_seeds(lstring, param, lsystem.nb_plantes)

    lsystem.clear()
    print(''.join((name, " - done")))

    if opt_repro:
        return mat


# Algorithme de reproduction des générations via le modèle génétique
def simpraise(plan_sim=None, id_scenario=0):
    if plan_sim is None:
        raise NameError('Pas de plan de simulation chargé.')
    row = plan_sim.iloc[id_scenario]
    if (not row['option_reproduction']) | (not row['option_floraison']):   #reproduction impossible
        return runlsystem(plan_sim=plan_sim, id_scenario=id_scenario, id_gener=1)

    # Config des fichiers d'entrée
    INPUTS_DIRPATH = 'inputs'
    src = os.path.join(INPUTS_DIRPATH, 'insim.txt')
    dst = 'modelgenet'
    exe = 'simpraise.exe'

    # Génération des fondateurs
    prf.rungenet(src, dst, exe, None)

    # Boucle des générations
    for i in range(1, row['num_gener']+1):
        mat = runlsystem(plan_sim=plan_sim, id_scenario=id_scenario, id_gener=i)
        prf.rungenet(src, dst, exe, mat)


plan = pd.read_csv(os.path.join('inputs', "plan_simulation.csv"))

runlsystem(plan_sim=plan, id_scenario=1)

#simpraise(plan_sim=plan, id_scenario=0)
