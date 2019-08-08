# coding: utf8
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math


class GraphicOutputs:
    def __init__(self):
        self.color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

    def graph_induction(self, path, csv_name):
        df = pd.read_csv(path + csv_name + '.csv')
        plt.clf()
        number_of_tiller = float(len(set(df.Id_talle.tolist())))
        colors = int(math.ceil(number_of_tiller / 8)) * self.color_list
        for tiller in set(df.Id_talle.tolist()):
            col = colors[int(tiller)]
            plt.xlabel(u'GDD (°Cj)')
            plt.ylabel(u'Induction rate (%)')
            plt.plot(u'GDD', u'Vernalisation_rate', data=df[(df.Id_talle == tiller)], color=col, linestyle='-')
            plt.plot(u'GDD', u'Secondary_induction_rate', data=df[(df.Id_talle == tiller)], color=col, linestyle=':')
        plt.savefig(path + csv_name + '.png')

    def graph_length(self, path, csv_name):
        df = pd.read_csv(path + csv_name + '.csv')
        plt.clf()
        number_of_tiller = len(np.unique(df.Id_talle.tolist()))
        plt.subplots_adjust(wspace=0.5, hspace=0.5)
        i = 0
        for tiller in set(df.Id_talle.tolist()):
            i += 1
            plt.subplot(number_of_tiller, 1, i)
            plt.title(u'Tiller ' + str(tiller))
            plt.xlabel(u'GDD (°Cj)')
            plt.ylabel(u'Length (mm)')
            rank_number = float(len(set(df[df.Id_talle == tiller].Id_rang.tolist())))
            colors = int(math.ceil(rank_number/8))*self.color_list
            for rank in set(df[df.Id_talle == tiller].Id_rang.tolist()):
                col = colors[rank - 1]
                plt.plot(u'GDD', u'Length', data=df[(df.Id_talle == tiller) & (df.Id_rang == rank) &
                                                    (df.Organ == u'internode')], color=col, linestyle='-')
                plt.plot(u'GDD', u'Length', data=df[(df.Id_talle == tiller) & (df.Id_rang == rank) &
                                                    (df.Organ == u'sheath')], color=col, linestyle=':')
                plt.plot(u'GDD', u'Length', data=df[(df.Id_talle == tiller) & (df.Id_rang == rank) &
                                                    (df.Organ == u'limb')], color=col, linestyle='--')
        plt.savefig(path + csv_name + '.png')
