# '''
# Created on 19/03/2020
#
# @author: modelisation - TR
# '''


# import the modules necessary to initiate the L-systems
import sys
import os
import pandas as pd
import random
from lgrass import param_plantes
from openalea.lpy import Lsystem
from openalea.plantgl.all import *
from lgrass import meteo_ephem
from lgrass import caribu

nb_graines_par_epillet = 3


def runlsystem(id_scenario=0):
    # Charger le plan de simulation
    INPUTS_DIRPATH = 'inputs'
    simul_conditions = pd.read_csv(os.path.join(INPUTS_DIRPATH, "plan_simulation.csv"))
    row = simul_conditions.iloc[id_scenario]
    name = str(row["name"])

    # Charger le lsystem
    lpy_filename = os.path.join('lgrass.lpy')
    lsystem = Lsystem(lpy_filename)

    # Parametres des plantes
    lsystem.ParamP, lsystem.nb_plantes, lsystem.NBlignes, lsystem.NBcolonnes, lsystem.posPlante, lsystem.Plantes, lsystem.Genotypes, lsystem.flowering_model = param_plantes.define_param(out_param_file=name)

    # Donnees meteo
    lsystem.meteo = meteo_ephem.import_meteo_data(row["meteo_path"], row['sowing_date'], row['site'])

    # Gestion caribu
    opt_caribu = row["option_caribu"]
    if opt_caribu:
        dico_caribu = caribu.init_caribu(meteo=lsystem.meteo, nb_plantes=lsystem.nb_plantes, scenario=row)

    # Parametres de simulation
    lsystem.option_tallage = row["option_tallage"]
    lsystem.option_tiller_regression = row["option_tiller_regression"]
    lsystem.option_mophogenetic_regulation_by_carbone = row["option_mophogenetic_regulation_by_carbone"]
    lsystem.derivationLength = int(row["derivationLength"])
    # lsystem.meteo_path = os.path.join()
    lsystem.sowing_date = row["sowing_date"]
    lsystem.site = row["site"]
    # lsystem.flowering_model = flowering_param
    lsystem.output_induction_file_name = name + '_' + 'induction'
    lsystem.output_organ_lengths_file_name = name + '_' + 'organ_lengths'
    lsystem.cutting_dates = [] if pd.isna(row["cutting_dates"]) \
        else [row["cutting_dates"]] if isinstance(row["cutting_dates"], int) \
        else [int(i) for i in row["cutting_dates"].split("_")]

    # Gestion de la creation d une nouvelle generation
    Epillets = []
    path_out = 'outputs/' + row[0] + '.csv'
    output = open(path_out, 'w')
    output.write("GDD;Date;Day;nb_talles;biomasse_aerienne;surface_foliaire" + "\n")

    # Lancement du lsystem
    lsystem.BiomProd = [0.] * lsystem.nb_plantes
    lstring = lsystem.axiom
    for dd in range(lsystem.derivationLength):
        try:
            day = lsystem.current_day
        except:
            day = 1
        lstring = lsystem.derive(lstring, dd, 1)
        lscene = lsystem.sceneInterpretation(lstring)
        if opt_caribu:
            try:
                lsystem.BiomProd, dico_caribu['radiation_interception'], dico_caribu[
                    'Ray'] = caribu.apply_caribu_lgrass(lstring, lscene, lsystem.TPS,
                                                        lsystem.current_day,
                                                        lsystem.tiller_appearance,
                                                        lsystem.nb_plantes, dico_caribu, day)
            except:
                continue
        if lsystem.current_day>day :
            output.write(";".join(
                [str(lsystem.TPS), str(lsystem.sowing_date), str(lsystem.current_day), str(lsystem.nb_talle),
                    str(lsystem.BiomProd[0]), str(lsystem.rapportS9_SSol_dict[0])]) + "\n")

        Viewer.display(lscene)

    # Gestion des graines
    for mod in lstring:
        if mod.name in ('apex',):
            if mod[0].phenological_state == 'reproductive' and mod[0].final_spikelet_number is not None:
                Epillets.append((mod[0].id_plante, mod[0].id_talle, mod[0].final_spikelet_number))
    nb_tot_graines = []
    for tupple in Epillets:
        for i in range(int(tupple[2]) * nb_graines_par_epillet):
            nb_tot_graines.append(tupple)
    Graines = []
    for rand in range(lsystem.nb_plantes):
        if len(nb_tot_graines) != 0:
            Graines.append(nb_tot_graines[random.randint(0, len(nb_tot_graines) - 1)])

    output.close()

    lsystem.clear()
    print(''.join((name, " - done")))


# for x in range(simul_conditions.shape[0]):
runlsystem(id_scenario=0)
runlsystem(id_scenario=1)
