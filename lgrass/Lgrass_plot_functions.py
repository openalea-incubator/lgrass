# coding: utf8
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
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
            plt.xlabel(u'Day from Oct 1rst')
            plt.ylabel(u'Induction rate (%)')
            plt.plot(u'Day', u'Vernalisation_rate', data=df[(df.Id_talle == tiller)], color=col, linestyle='-')
            plt.plot(u'Day', u'Secondary_induction_rate', data=df[(df.Id_talle == tiller)], color=col, linestyle=':')
        plt.savefig(path + csv_name + '.pdf')

    def graph_length(self, path, csv_name):
        pdf_file = matplotlib.backends.backend_pdf.PdfPages(path + csv_name + '.pdf')
        df = pd.read_csv(path + csv_name + '.csv')
        number_of_tiller = len(np.unique(df.Id_talle.tolist()))
        i = 0
        for tiller in set(df.Id_talle.tolist()):
            plt.figure()
            ax = plt.subplot(111)
            ax.set_title(u'Tiller ' + str(tiller))
            ax.set_xlabel(u'Day from Oct 1rst')
            ax.set_ylabel(u'Length (mm)')
            rank_number = float(len(set(df[df.Id_talle == tiller].Id_rang.tolist())))
            colors = int(math.ceil(rank_number/8))*self.color_list
            for rank in set(df[df.Id_talle == tiller].Id_rang.tolist()):
                col = colors[rank - 1]
                ax.plot(u'Day', u'Length', data=df[(df.Id_talle == tiller) & (df.Id_rang == rank) &
                                                    (df.Organ == u'internode')], color=col, linestyle='-')
                ax.plot(u'Day', u'Length', data=df[(df.Id_talle == tiller) & (df.Id_rang == rank) &
                                                    (df.Organ == u'sheath')], color=col, linestyle=':')
                ax.plot(u'Day', u'Length', data=df[(df.Id_talle == tiller) & (df.Id_rang == rank) &
                                                    (df.Organ == u'limb')], color=col, linestyle='--')
            pdf_file.savefig()
        pdf_file.close()


graph = GraphicOutputs()
graph.graph_length(r'D:\Simon\Python\lgrass\lgrass\outputs', '\output_organ_lengths')