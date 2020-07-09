# -*- coding: latin-1 -*-
"""
    test lgrass
    ~~~~~~~~~~~~~~~

    Test the model L-grass.

    You must first install :mod:`lgrass` (and add it to your PYTHONPATH)
    before running this script with the command `python`.

    :copyright: Copyright TODO, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Verdenal et al, 2008.
"""


import os

import numpy as np
import pandas as pd

from openalea.lpy import Lsystem
import lgrass

INPUTS_DIRPATH = 'inputs'
OUTPUTS_DIRPATH = 'outputs'

NSTEP = 800
PRECISION = 3
RELATIVE_TOLERANCE = 10**-PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE

DESIRED_SERIE_FOLIAIRE_FILENAME = 'desired_sorties_feuilles_finales_50_{}_30.csv'.format(NSTEP)
DESIRED_SORTIE_SURFACE_BIOMASS_FILENAME = 'desired_surface_biomass_50_{}_30.csv'.format(NSTEP)
DESIRED_OUTPUT_INDUCTION_FILENAME = 'desired_output_induction.csv'.format(NSTEP)
DESIRED_OUTPUT_ORGAN_LENGTHS_FILENAME = 'output_organ_lengths.csv'.format(NSTEP)


def compare_actual_to_desired(data_dirpath, actual_data_df, desired_data_filename, actual_data_filename=None, overwrite_desired_data=False):
    # read desired data

    if actual_data_filename is not None:
        actual_data_filepath = actual_data_filename
        actual_data_df.to_csv(actual_data_filepath, na_rep='NA', index=False)

    if overwrite_desired_data:
        desired_data_filepath = os.path.join(data_dirpath, desired_data_filename)
        actual_data_df.to_csv(desired_data_filepath, na_rep='NA', index=False)
    else:
        desired_data_filepath = os.path.join(data_dirpath, desired_data_filename)
        desired_data_df = pd.read_csv(desired_data_filepath)

        # keep only numerical data
        for column in ('topology', 'Phase', 'Date', 'Site', 'Organ'):
            if column in desired_data_df.columns:
                del desired_data_df[column]
                del actual_data_df[column]

        # compare to the desired data
        np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run(overwrite_desired_data=False):
    lpy_filename = os.path.join(lgrass.__path__[0], "lgrass.lpy")
    lsys = Lsystem(lpy_filename)
    axiom = lsys.axiom

    lsys.meteo_path = os.path.join(INPUTS_DIRPATH, 'meteo_file.csv')
    lsys.INPUTS_DIRPATH = INPUTS_DIRPATH
    lsys.OUTPUTS_DIRPATH = OUTPUTS_DIRPATH
    lsys.derivationLength = NSTEP
    lsys.DureeExp = NSTEP
    lsys.option_tallage = True
    lsys.option_senescence = True
    lsys.option_floraison = True
    lsys.option_caribu = 'Off'
    lsys.option_tiller_regression = False
    lsys.option_mophogenetic_regulation_by_carbone = False
    lsys.derive(axiom, NSTEP)

    # convert the outputs to Pandas dataframe
    surface_biomass = pd.read_csv(lsys.chemin_fichier1.name)
    evol = pd.read_csv(lsys.chemin_fichier2.name)
    output_induction_file_path = os.path.join(OUTPUTS_DIRPATH, 'output_induction.csv')
    output_induction = pd.read_csv(output_induction_file_path)
    output_organ_lengths_file_path = os.path.join(OUTPUTS_DIRPATH, 'output_organ_lengths.csv')
    output_organ_lengths = pd.read_csv(output_organ_lengths_file_path)

    # # compare outputs
    compare_actual_to_desired(OUTPUTS_DIRPATH, surface_biomass, DESIRED_SORTIE_SURFACE_BIOMASS_FILENAME, lsys.chemin_fichier1.name, overwrite_desired_data)
    compare_actual_to_desired(OUTPUTS_DIRPATH, evol, DESIRED_SERIE_FOLIAIRE_FILENAME, lsys.chemin_fichier2.name, overwrite_desired_data)
    compare_actual_to_desired(OUTPUTS_DIRPATH, output_induction, DESIRED_OUTPUT_INDUCTION_FILENAME, output_induction_file_path, overwrite_desired_data)
    compare_actual_to_desired(OUTPUTS_DIRPATH, output_organ_lengths, DESIRED_OUTPUT_ORGAN_LENGTHS_FILENAME, output_organ_lengths_file_path, overwrite_desired_data)

    if overwrite_desired_data:
        print("New desired files written")
    else:
        print("Test passed successfully")

if __name__ == '__main__':
    test_run(overwrite_desired_data=False)
