from example.MORPHOFLOR import run_lgrass
import multiprocessing as mp
import os
import pandas as pd
import time
from functools import partial
from shutil import copyfile

if __name__ == '__main__':

    simul_name = 'test11'
    outputs_dir_path = os.path.join('outputs', simul_name)

    if not os.path.exists(outputs_dir_path):
        os.mkdir(outputs_dir_path)
    else:
        raise Warning('Outputs_dir_path already exists')

    copyfile(os.path.join('inputs', 'plan_simulation.csv'), os.path.join(outputs_dir_path, 'plan_simulation.csv'))

    scenarii_df = pd.read_csv(os.path.join('inputs', 'plan_simulation.csv'), index_col='Scenario')
    scenarii_df['Scenario'] = scenarii_df.index
    scenarii = scenarii_df.Scenario

    tstart = time.time()
    num_processes = mp.cpu_count() - 1
    p = mp.Pool(num_processes)
    run_lg = partial(run_lgrass.run_lgrass,
                     inputs_dir_path='inputs',
                     outputs_dir_path=outputs_dir_path)
    mp_solutions = p.map(run_lg, list(scenarii))
    p.terminate()
    p.join()

    tend = time.time()
    tmp = (tend - tstart) / 60.

    print("multiprocessing: %8.3f minutes" % tmp)

