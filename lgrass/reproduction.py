# coding: utf8
import random
import numpy as np


# Créer la matrice de croisement des plantes pour établir une nouvelle génération
def create_seeds(lstring, param, nb_plantes):
    seed_number = int(param.loc[param['name'] == 'seeds_by_spikelet', :]['value'])
    matrix = np.zeros((nb_plantes, nb_plantes))
    seeds = []

    # construction des graines et de leurs parents
    for mod in lstring:
        if mod.name in ('apex',):
            if mod[0].final_spikelet_number is not None:
                for i in range(int(mod[0].final_spikelet_number * seed_number)):
                    id_mere = mod[0].id_plante  # séléction de la mère
                    peres = [j for j in range(nb_plantes)]
                    peres.remove(id_mere)
                    id_pere = random.choice(peres)  # séléction aléatoire du père
                    seeds.append((id_mere, id_pere))

    # construction de la matrice de croisement des plantes
    for i in range(len(seeds)):
        matrix[seeds[i][0], seeds[i][1]] += 1

    return matrix
