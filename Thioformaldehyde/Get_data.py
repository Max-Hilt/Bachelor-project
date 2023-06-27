# -*- coding: utf-8 -*-
# Author: Max Hilt
# Importing libraries
import numpy as np
import yaml
import os
import os.path
import argparse
from yaml.loader import BaseLoader
import shutil


def get_path_data():

    # Creating datastructure
    total_energy = np.empty((steps, num_excitations+1))

    for i in range(steps):
        # First we do the ground state energy
        print(f'Collecting data of step {i}')
        with open(f'gradient_path/step_{i}/molgw_gw.yaml') as f:
            data = yaml.load(f, Loader=BaseLoader)

            # Performing a simple check to see if an obvious simulation error has occured
            if os.path.isfile(f'gradient_path/step_{i}/ENERGY_QP'):
                total_energy[i,0] = float(data['scf energy']['total']) * Ha_to_ev

                # Now we add the BSE excitation energy spectra
                for j in range(num_excitations):
                    with open(f'gradient_path/step_{i}/molgw_bse.yaml') as f:
                        data = yaml.load(f, Loader=BaseLoader)
                        total_energy[i,j+1] = float(data['optical spectrum']['excitations']['energies'][f'{j+1}'])
            else:
                print(f'Iteration {i} has some error in the simulation')
                total_energy[i,0] = 99999

                # Now we add the BSE excitation energy spectra
                for j in range(num_excitations):
                    with open(f'gradient_path/step_{i}/molgw_bse.yaml') as f:
                        data = yaml.load(f, Loader=BaseLoader)
                        total_energy[i,j+1] = 99999


        
    # Saving this data
    np.save(f'Path_energies',total_energy)


    return


parser = argparse.ArgumentParser()
parser.add_argument("-n", dest="number_of_iterations", required=True, type=int, help="The total number of iterations that happened")
args = parser.parse_args()



# Setting parameters --This can maybe be automated better--
steps = args.number_of_iterations
num_excitations = 3
Ha_to_ev = 27.2114

get_path_data()