# legume : 1 wrapper = 1 usm = 1 espèce de plante


class Indexer:
    def __init__(self, global_order=[], wheat_names=[], legume_names=[], other_names=[]) -> None:
        self.global_order: list = global_order
        self.wheat_names: list = wheat_names
        self.legume_names: list = legume_names
        self.other_names: list = other_names

        self.legume_active: bool = legume_names != []
        self.wheat_active: bool = wheat_names != []

        self.legume_number_of_species: list = [1] * len(legume_names)

        self.wheat_index: list = [global_order.index(name) for name in wheat_names if name in global_order]
        self.legume_index: list = [global_order.index(name) for name in legume_names if name in global_order]
        self.other_index: list = [global_order.index(name) for name in other_names if name in global_order]

    def update_legume_several_species(self, name, number_of_species=None):
        global_id = self.global_order.index(name)
        legume_id = self.legume_names.index(name)

        if number_of_species is not None:
            self.legume_number_of_species[legume_id] = number_of_species

        # insert another instance name
        for i in range(self.legume_number_of_species[legume_id] - 1):
            self.global_order.insert(global_id, name)
            self.legume_names.insert(legume_id, name)

        # update global index of each fspm
        self.wheat_index = [self.global_order.index(name) for name in self.wheat_names if name in self.global_order]
        self.other_index = [self.global_order.index(name) for name in self.other_names if name in self.global_order]
        self.legume_index = []
        n_before = ""
        for n in self.legume_names:
            if n != n_before:
                self.legume_index.extend([index for index, item in enumerate(self.global_order) if item == n])
                n_before = n


    def light_scenes_mgmt(self, scenes_dict):
        scenes = [0.] * len(self.global_order)
        for name, geo in scenes_dict.items():
            scenes[self.global_order.index(name)] = geo
        return scenes
    
    def soil_inputs(self, soil_dict):
        out = [[0.] * len(self.global_order) for i in range(4)]
        for name, tuple_var in soil_dict.items():
            for v, w in zip(out, tuple_var):
                v[self.global_order.index(name)] = w

        return tuple(out)
