# Generate a virtual grassland with 2 FSPMs : L-grass & L-egume
# using PlantFusion
# See example simulations : https://github.com/mwoussen/plantfusion/blob/develop/simulations

from plantfusion.l_egume_wrapper import L_egume_wrapper # L-egume FSPM management
from custom_wrappers.l_grass_wrapper import L_grass_wrapper # L-grass FSPM management : TO BE IMPLEMENTED
# see : https://github.com/mwoussen/plantfusion/blob/develop/simulations/other_fspm.py
from plantfusion.light_wrapper import Light_wrapper # light management using LightVegeManager
from plantfusion.soil_wrapper import Soil_wrapper # soil management
from plantfusion.indexer import Indexer # index management
from plantfusion.planter import Planter # plant position management

import time
import datetime
import os

def simulation(in_folder, out_folder, idusm1, idusm2, writegeo=False):
    try:
        # Create target Directory
        os.mkdir(os.path.normpath(out_folder))
        print("Directory ", os.path.normpath(out_folder), " Created ")
    except FileExistsError:
        print("Directory ", os.path.normpath(out_folder), " already exists")
    

    ######################
    ### INITIALIZATION ###
    ######################
        
    legume_name = "legume1"
    lgrass_name = "lgrass"
    indexer = Indexer(global_order=[legume_name, lgrass_name], legume_names=[legume_name], lgrass_names=[lgrass_name])

    generation_type = "random"
    plant_density = {legume_name : 250, lgrass_name : 450} # plantes.m-2
    planter = Planter(generation_type=generation_type, indexer=indexer, plant_density=plant_density)
    #note : arg "inter_rows" is not needed for random generation (I think), but is set to 0.15 anyway.

    legume = L_egume_wrapper(
        name=legume_name,
        indexer=indexer,
        in_folder=in_folder, out_folder=out_folder,
        planter=planter,
        caribu_scene=True
        
    )

    lgrass = L_grass_wrapper(
        name=lgrass_name,
        indexer=indexer,
        in_folder=in_folder, out_folder=out_folder,
        planter=planter
        #ADD ARGS SPECIFIC TO L-GRASS
    )

    sky = "turtle46"
    lighting = Light_wrapper(
        lightmodel="caribu",
        out_folder=out_folder,
        sky=sky,
        planter=planter, # why is it needed ?
        indexer=indexer,
        writegeo=writegeo
    )

    soil = Soil_wrapper(in_folder=in_folder, 
                        out_folder=out_folder, 
                        IDusm=1711, # unit of simulation
                        planter=planter, 
                        opt_residu=0, 
                        save_results=True)
    

    ##################
    ### SIMULATION ###
    ##################

    current_time_of_the_system = time.time()
    for t in range(legume.lsystem.derivationLenght) :# chose derivation lenght

        legume.derive()

        while lgrass.lsystem.current_day > t:
            lgrass.derive()

        scene_legume = legume.light_inputs(elements="triangles")
        scene_lgrass = lgrass.light_inputs() # ARGS DEPENDS ON LGRASS CLASS IMPLEMENTATION
        scenes = indexer.light_scenes_mgmt({lgrass_name : scene_lgrass, legume_name : scene_legume})

        lighting.run(scenes=scenes, day=legume.doy(), parunit="RG") #what is parunit="RG" ?
        legume.light_results(legume.energy(), lighting)
        lgrass.light_results() # ARGS DEPENDS ON LGRASS CLASS IMPLEMENTATION

        soil_legume_inputs = legume.soil_inputs()
        soil_lgrass_inputs = lgrass.soil_inputs()

        (
            N_content_roots_per_plant,
            roots_length_per_plant_per_soil_layer,
            plants_soil_parameters,
            plants_light_interception, # what are these output's type ?
        ) = indexer.soil_inputs({legume_name : soil_legume_inputs, lgrass_name : soil_lgrass_inputs})

        soil.run( # update soil inputs
            legume.doy(), #doy : day of the year (julian day)
            N_content_roots_per_plant,
            roots_length_per_plant_per_soil_layer,
            plants_soil_parameters,
            plants_light_interception,
        )

        #compute new plant parameters for next derivation. To be implemented in Lgrass wrapper
        legume.run()
        lgrass.run()
    execution_time = int(time.time() - current_time_of_the_system)
    print("\n" "Simulation run in {}".format(str(datetime.timedelta(seconds=execution_time))))

    # end() effect is specific to the model
    legume.end()
    lgrass.end()
    soil.end()

if __name__ == "__main__":
    in_folder = "inputs_soil_legume"
    out_folder = "outputs/two_legume_caribu"
    id1 = 17111
    id2 = 17112
    writegeo = True

    simulation(in_folder, out_folder, id1, id2, writegeo=writegeo)