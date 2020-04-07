# coding=utf-8
# '''
# Created on 19/03/2020
#
# @author: modelisation - TR
# '''

# import the modules necessary to initiate the L-systems
import sys
import os
import pandas as pd
import openalea.lpy as opy
import openalea.plantgl as opal
from lgrass import param_plantes
from lgrass import reproduction
from lgrass import meteo_ephem
from lgrass import caribu


def runlsystem(id_scenario=0):
    # Charger le plan de simulation
    INPUTS_DIRPATH = 'inputs'
    OUTPUTS_DIRPATH = 'outputs'
    simul_conditions = pd.read_csv(os.path.join(INPUTS_DIRPATH, "plan_simulation.csv"))
    row = simul_conditions.iloc[id_scenario]
    name = str(row["name"])

    # Charger le lsystem
    lpy_filename = os.path.join('lgrass.lpy')
    lsystem = opy.Lsystem(lpy_filename)

    # Parametres des plantes
    param = param_plantes.get_param()
    lsystem.ParamP, lsystem.nb_plantes, lsystem.NBlignes, lsystem.NBcolonnes, lsystem.posPlante, lsystem.Plantes, lsystem.Genotypes, lsystem.flowering_model = param_plantes.define_param(
        out_param_file=os.path.join(INPUTS_DIRPATH,name), param=param)

    # Donnees meteo
    lsystem.meteo = meteo_ephem.import_meteo_data(row["meteo_path"], row['sowing_date'], row['site'])

    # Gestion caribu
    opt_caribu = row["option_caribu"]
    if opt_caribu:
        day = 1
        dico_caribu = caribu.init(meteo=lsystem.meteo, nb_plantes=lsystem.nb_plantes, scenario=row)
        lsystem.BiomProd = [0.] * lsystem.nb_plantes

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
    lsystem.cutting_dates = [] if pd.isna(row["cutting_dates"]) \
        else [row["cutting_dates"]] if isinstance(row["cutting_dates"], int) \
        else [int(i) for i in row["cutting_dates"].split("_")]

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
            if lsystem.current_day>day :
                output.write(";".join(
                    [str(lsystem.TPS), str(lsystem.sowing_date), str(lsystem.current_day), str(lsystem.nb_talle),
                        str(lsystem.BiomProd[0]), str(lsystem.rapportS9_SSol_dict[0])]) + "\n")
        opal.all.Viewer.display(lscene)

    # Gestion de la reproduction
    print(reproduction.create_seeds(lstring, param, lsystem.nb_plantes))

    output.close()
    lsystem.clear()
    print(''.join((name, " - done")))


runlsystem(id_scenario=0)
# runlsystem(id_scenario=1)
# runlsystem(id_scenario=2)
