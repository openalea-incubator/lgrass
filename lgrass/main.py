# coding: utf8
# '''
# Created on 21/04/2020
#
# @author: modelisation - TR
# '''
import multiprocessing
import time
import pandas as pd
import lgrass_batch_simpraise as batch


if __name__ == '__main__':
    timing = time.time()
    plan = pd.read_csv("inputs/plan_simulation.csv", sep=';')

    # for i in range(9,17):
    #     batch.runlsystem(plan, i, 1)

    ###    Utilisation de Pool multiprocess    ###

    multiprocessing.freeze_support()
    CPUnb = int(multiprocessing.cpu_count()) # nombre de processeurs utilises
    pool = multiprocessing.Pool(processes=3)
    for j in range(1,3):
        pool.apply_async(batch.runlsystem, args=(plan, j, 1))
    pool.close()
    pool.join()
    print('Global execution time : ', time.time() - timing)

    ###   Utilisation de Process multiprocess   ###

    # q = multiprocessing.Queue()
    # p = multiprocessing.Process(target=batch.runlsystem, args=(plan, 16, 1, q))
    # p.start()
    # p.join()