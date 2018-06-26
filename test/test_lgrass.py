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

"""
    Information about this versioned file:
        $LastChangedBy: cchambon $
        $LastChangedDate: 2016-10-24 17:11:32 +0200 (lun., 24 oct. 2016) $
        $LastChangedRevision: 29 $
        $URL: https://subversion.renater.fr/senesc-wheat/trunk/test/test_senescwheat.py $
        $Id: test_senescwheat.py 29 2016-10-24 15:11:32Z cchambon $
"""

import os

import numpy as np
import pandas as pd

from openalea.lpy import Lsystem
import lgrass

INPUTS_DIRPATH = 'inputs'
OUTPUTS_DIRPATH = 'outputs'

DESIRED_EVOL_FILENAME = 'desired_evol.csv'
DESIRED_SORTIE_SIMULVALIDCOUPE_FILENAME = 'desired_sortie simul ValidCoupe.csv'
DESIRED_SORTIE_SURFACE_BIOMASS_FILENAME = 'desired_surface_biomass.csv'

NSTEP = 200
PRECISION = 6
RELATIVE_TOLERANCE = 10**-PRECISION
ABSOLUTE_TOLERANCE = RELATIVE_TOLERANCE


def compare_actual_to_desired(data_dirpath, actual_data_df, desired_data_filename, actual_data_filename=None):
    # read desired data
    desired_data_filepath = os.path.join(data_dirpath, desired_data_filename)
    desired_data_df = pd.read_csv(desired_data_filepath)

    if actual_data_filename is not None:
        actual_data_filepath = actual_data_filename
        actual_data_df.to_csv(actual_data_filepath, na_rep='NA', index=False)

    # keep only numerical data
    for column in ('nb_feuillemergees', 'organ', 'element'):
        if column in desired_data_df.columns:
            del desired_data_df[column]
            del actual_data_df[column]

    # compare to the desired data
    np.testing.assert_allclose(actual_data_df.values, desired_data_df.values, RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE)


def test_run():
    lpy_filename = os.path.join(lgrass.__path__[0], "lgrass.lpy")
    lsys = Lsystem(lpy_filename)
    axiom = lsys.axiom

    lsys.INPUTS = INPUTS_DIRPATH
    lsys.OUTPUTS = OUTPUTS_DIRPATH
    lstring = lsys.derive(axiom, NSTEP)

    # convert the outputs to Pandas dataframe
    surface_biomass = pd.read_csv(lsys.chemin_fichier.name)
    sortie_simul_ValidCoupe = pd.read_csv(lsys.chemin_fichier3.name)
    evol = pd.read_csv(lsys.chemin_fichier5.name)

    # # compare outputs
    compare_actual_to_desired(OUTPUTS_DIRPATH, surface_biomass, DESIRED_SORTIE_SURFACE_BIOMASS_FILENAME, lsys.chemin_fichier.name)
    compare_actual_to_desired(OUTPUTS_DIRPATH, sortie_simul_ValidCoupe, DESIRED_SORTIE_SIMULVALIDCOUPE_FILENAME, lsys.chemin_fichier3.name)
    compare_actual_to_desired(OUTPUTS_DIRPATH, evol, DESIRED_EVOL_FILENAME, lsys.chemin_fichier5.name)

    print "Test passed successfully"

if __name__ == '__main__':
    test_run()