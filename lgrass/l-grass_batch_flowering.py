'''
Created on 11/10/2018

@author: modelisation - SR
'''
#batch pour L-grass with flowering



#import the modules necessary to initiate the L-systems
import time
import math
import re
from openalea.lpy import *
from openalea.plantgl.all import *
from alinea.caribu.CaribuScene import CaribuScene
import multiprocessing
import sys
path_ = r'D:\Simon\Python\lgrass\lgrass'
#sys.path.insert(0, path_)

import os
#from generateScene import *
import numpy as np
import datetime
from lgrass import flowering_functions
import csv
from openalea.lpy import Lsystem
import itertools

OUTPUTS_DIRPATH = 'outputs'
# date_string = datetime.datetime.now() .strftime('%Y_%m_%d_%Hh%M')
# output_batch_path = open(os.path.join(OUTPUTS_DIRPATH, 'output_batch_path', date_string), "wb")
# output_batch = csv.writer(output_batch_path)
# output_batch.writerow(['vai', 'vbee', 'sldl', 'absolute_max_leaf_number', 'absolute_min_leaf_number'])



#cree la liste de L-systems et liste des noms
testsim={}
names = []


nb_simul = 0
daily_vern_rate_list = np.arange(0.001, 0.010, 1)
basic_vern_rate_list = np.arange(0.01, 0.20, 1)
photoperiod_min_list = np.arange(10, 12, 10)
photoperiod_max_list = np.arange(14, 16, 10)


for x in itertools.product(daily_vern_rate_list, basic_vern_rate_list, photoperiod_min_list, photoperiod_max_list):
    print(x)
    nb_simul += 1
    flowering_param = flowering_functions.FloweringFunctions()
    flowering_param.param.daily_vern_rate = x[0]
    flowering_param.param.basic_vern_rate = x[1]
    flowering_param.param.photoperiod_min = x[2]
    flowering_param.param.photoperiod_max = x[3]
    name = str(x[0]) + '_' + str(x[1]) + '_' + str(x[2]) + '_' + str(x[3])
    names.append(name)
    lpy_filename = os.path.join('lgrass.lpy')
    testsim[name] = Lsystem(lpy_filename)

    testsim[name].flowering_model = flowering_param

# function to run an L-system from the 'testsim' dictionnary


def runlsystem(n):
    testsim[names[n]].derive() # permet le declenchement de la fonction "End" du script lpy
    print(testsim[names[n]].output_dict)
    with open(os.path.join(OUTPUTS_DIRPATH, 'sortie_test_path', str(n) + '.csv'), 'wb') as sortie_test_path:
        sortie_test = csv.writer(sortie_test_path)
        sortie_test.writerows(testsim[names[n]].output_dict.items())
    testsim[names[n]].clear()
    print(''.join((names[n]," - done")))

    # print(testsim[names[n]].output_dict)


runlsystem(2)

# #run the L-systems
# if __name__ == '__main__':
#     multiprocessing.freeze_support()
#     CPUnb = multiprocessing.cpu_count()-1 #nombre de processeurs, moins un par prudence. (et pour pouvoir faire d'autres choses en meme temps)
#     print 'nb CPU: ' + str(CPUnb)
#     pool = multiprocessing.Pool(processes=CPUnb)
#     for i in range(int(nb_simul)):
#         pool.apply_async(runlsystem, args=(i,)) #Lance CPUnb simulations en meme temps, lorsqu'une simulation se termine elle est immediatement remplacee par la suivante
#     pool.close()
#     pool.join()