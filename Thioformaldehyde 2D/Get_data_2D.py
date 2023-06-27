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
    total_energy = np.empty((steps,steps, num_excitations+1))

    # Performing two for loops to save all the data
    for i in range(steps):
        for j in range(steps):
            # First we take the groundstate energy of each step
            print(f'Collecting data of step {i*steps + j}')

            with open(f'Simulations/location_{i*steps + j}/molgw_gw.yaml') as f:
                data = yaml.load(f, Loader=BaseLoader)

            # Performing a simple check to see if an obvious simulation error has occured
            if os.path.isfile(f'Simulations/location_{i*steps + j}/ENERGY_QP'):

                # And saving the groundstate energy
                total_energy[i,j,0] = float(data['scf energy']['total']) * Ha_to_ev

                # Now we add the BSE excitation energy spectra
                with open(f'Simulations/location_{i*steps + j}/molgw_bse.yaml') as f:
                    data = yaml.load(f, Loader=BaseLoader)
                    for k in range(num_excitations): 
                        total_energy[i,j,k+1] = float(data['optical spectrum']['excitations']['energies'][f'{k+1}'])
            else:
                print(f'Iteration {i*steps + j} has some error in the simulation')

                # So we set the value of the energy to 99999 for all states
                total_energy[i,j,0] = 99999

                for k in range(num_excitations):
                    total_energy[i,j,k+1] = 0


        
    # Saving this data
    np.save(f'Surface_energy',total_energy)

    return


# Setting parameters 
steps = int(input('How many steps did the simulation have? '))
print(f'The simulation has {steps} steps')
num_excitations = 3
Ha_to_ev = 27.2114

get_path_data()