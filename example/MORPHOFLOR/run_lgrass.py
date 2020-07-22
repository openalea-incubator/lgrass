import os
import sys
import getopt

import pandas as pd
from lgrass import flowering_functions
import lgrass
from openalea.lpy import Lsystem


def run_lgrass(scenario_id=1, inputs_dir_path='inputs', outputs_dir_path='outputs'):
    """ Run lgrass.lpy file from external script

    :param str inputs_dir_path: the path to the input files (meteo, list of simulations...)
    :param str outputs_dir_path: the path to the outputs files
    :param int scenario_id: the index of the scenario to be run from the list of simulations CSV file

    """

    # Scenario to be run
    scenarii_df = pd.read_csv(os.path.join('inputs', 'plan_simulation.csv'), index_col='Scenario')
    scenario = scenarii_df.loc[scenario_id].to_dict()
    scenario_name = scenario['name']

    # Update parameters of flowering model
    flowering_model = lgrass.flowering_functions.FloweringFunctions()
    flowering_model.param.__dict__.update((k, scenario[k]) for k in set(scenario).intersection(flowering_model.param.__dict__))

    # Load lsystem
    lpy_filename = os.path.join(lgrass.__path__[0], "lgrass.lpy")
    lsystem = Lsystem(lpy_filename)

    # Update parameters
    lsystem.option_profile = "plateau"
    lsystem.flowering_model = flowering_model
    lsystem.derivationLength = int(scenario['derivationLength'])
    lsystem.option_tallage = scenario['option_tallage']
    lsystem.option_senescence = scenario['option_senescence']
    lsystem.option_floraison = scenario['option_floraison']
    lsystem.meteo_path = os.path.join(inputs_dir_path, scenario['meteo_filename'])
    lsystem.sowing_date = scenario['sowing_date']
    lsystem.site = scenario['site']
    lsystem.output_induction_file_name = '{}_induction'.format(scenario_name)
    lsystem.output_organ_lengths_file_name = '{}_organ_lengths'.format(scenario_name)
    lsystem.cutting_dates = [] if pd.isna(scenario['cutting_dates']) \
        else [scenario['cutting_dates']] if isinstance(scenario['cutting_dates'], int) \
        else [int(i) for i in scenario['cutting_dates'].split("_")]
    lsystem.ParamP[0]['C'] = scenario['value_C']
    lsystem.ParamP[0]['Premiecroiss'] = scenario['Premiecroiss']
    lsystem.ParamP[0]['PS_compensation_point'] = scenario['PS_compensation_point']
    if outputs_dir_path:
        lsystem.OUTPUTS_DIRPATH = outputs_dir_path

    # Run the lsystem
    try:
        lsystem.derive()  # permet le declenchement de la fonction "End" du script lpy
        lsystem.clear()
    except Exception as e:
        print(e)


if __name__ == '__main__':
    inputs = 'inputs'
    outputs = 'outputs'
    scenario = 14

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:s:d", ["inputs=", "outputs=", "scenario="])
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-i", "--inputs"):
            inputs = arg
        elif opt in ("-o", "--outputs"):
            outputs = arg
        elif opt in ("-s", "--scenario"):
            scenario = int(arg)

    run_lgrass(scenario_id=scenario, inputs_dir_path=inputs, outputs_dir_path=outputs)
