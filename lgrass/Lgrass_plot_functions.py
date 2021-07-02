# -*- coding: latin-1 -*-

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import numpy as np
import pandas as pd
import math
import os


class GraphicOutputs:
    def __init__(self, graph_dirpath):
        self.graph_dirpath = graph_dirpath
        self.color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    def graph_induction(self, csv_name):
        df = pd.read_csv(csv_name + '.csv')
        plt.ioff()
        plt.figure()
        number_of_tiller = float(len(set(df.Id_talle.tolist())))
        colors = int(math.ceil(number_of_tiller / len(self.color_list))) * self.color_list
        for tiller in set(df.Id_talle.tolist()):
            col = colors[int(tiller)]
            plt.xlabel(u'Day from sowing')
            plt.ylabel(u'Induction rate (%)')
            plt.plot(u'Day', u'Vernalisation_rate', data=df[(df.Id_talle == tiller)], color=col, linestyle='-')
            plt.plot(u'Day', u'Secondary_induction_rate', data=df[(df.Id_talle == tiller)], color=col, linestyle=':')
        plt.savefig(os.path.join(csv_name + '_' + '.pdf'))
        plt.close()

    def graph_length(self, csv_name):
        pdf_file = matplotlib.backends.backend_pdf.PdfPages(os.path.join(csv_name + '_' +
                                                                         'longueur_feuilles' + '.pdf'))
        df = pd.read_csv(csv_name + '.csv')
        for tiller in set(df.Id_talle.tolist()):
            plt.ioff()
            plt.figure()
            ax = plt.subplot(111)
            ax.set_title(u'Tiller ' + str(tiller))
            ax.set_xlabel(u'Day from sowing')
            ax.set_ylabel(u'Length (mm)')
            rank_number = float(len(set(df[df.Id_talle == tiller].Id_rang.tolist())))
            colors = int(math.ceil(rank_number / len(self.color_list)))*self.color_list
            for rank in set(df[df.Id_talle == tiller].Id_rang.tolist()):
                col = colors[rank - 1]
                ax.plot(u'Day', u'Length', data=df[(df.Id_talle == tiller) & (df.Id_rang == rank) &
                                                   (df.Organ == u'internode')], color=col, linestyle='-')
                ax.plot(u'Day', u'Length', data=df[(df.Id_talle == tiller) & (df.Id_rang == rank) &
                                                   (df.Organ == u'sheath')], color=col, linestyle=':')
                ax.plot(u'Day', u'Length', data=df[(df.Id_talle == tiller) & (df.Id_rang == rank) &
                                                   (df.Organ == u'limb')], color=col, linestyle='--')
            pdf_file.savefig()
            plt.close()
        pdf_file.close()

    def graph_tiller_number(self, csv_name):
        plt.figure()

        # Read csv file and group by plant
        df = pd.read_csv(csv_name)
        df_grouped = df.groupby('id plante')

        # Plot tiller number per plant
        for id_plant, data in df_grouped:
            data.groupby('TPS')['id talle'].count().plot(label=id_plant, legend=True)

        # Legends
        plt.xlabel(u'Temps (heure)')
        plt.ylabel(u'Nombre de talles par n° plante')
        # Save plot
        plt.savefig(os.path.join(self.graph_dirpath, 'Nombre_talles.PNG'))
        plt.close()

    def graph_LAI(self, csv_name, pattern):
        plt.figure()
        plt.ioff()

        # Read csv file
        df = pd.read_csv(csv_name)
        # plot LAI
        surface_pattern = pattern**2
        df['LAI'] = df['Surface_feuilles_emergees'] / surface_pattern
        df.groupby('TPS')['LAI'].sum().plot()

        # Legends
        plt.xlabel(u'Temps (heure)')
        plt.ylabel('LAI')
        # Save plot
        plt.savefig(os.path.join(self.graph_dirpath, 'LAI.PNG'))
        plt.close()

    def graph_leaf_number(self, csv_name):
        pdf_file = matplotlib.backends.backend_pdf.PdfPages(os.path.join(csv_name + '_' +
                                                                         'Nombre_feuilles' + '.pdf'))
        # Read csv file
        df = pd.read_csv(csv_name)
        df_grouped = df.groupby('topology')

        # plot leaf number
        for topo, data in df_grouped:
            plt.figure()
            data.groupby('TPS')['nb_feuille_emergees'].mean().plot(label=topo, legend=True)
            # Legends
            plt.xlabel(u'Temps (heure)')
            plt.ylabel('Leaf number per axis')
            # Save plot
            pdf_file.savefig()
            plt.close()
        pdf_file.close()
