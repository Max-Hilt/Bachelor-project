# -*- coding: utf-8 -*-
# Author: Max Hilt
# Importing libraries
import numpy as np
import yaml
import os
import os.path
from yaml.loader import BaseLoader
import shutil


def get_gw_data():
    # For GW with BSE we have two seperate yaml files per simulation
    # Creating datastructure
    total_energy = np.empty((steps, num_excitations+1))

    for i in range(steps):

        # First we do the ground state energy
        print(i)
        with open(f'{path}simulations/sim{i}/GW/molgw.yaml') as f:
            data = yaml.load(f, Loader=BaseLoader)

            # Performing a simple check to see if an obvious simulation error has occured
            if os.path.isfile(f'{path}simulations/sim{i}/BSE/ENERGY_QP'):
                total_energy[i,0] = float(data['scf energy']['total']) * Ha_to_ev

                # Now we add the BSE excitation energy spectra
                for j in range(num_excitations):
                    with open(f'{path}simulations/sim{i}/BSE/molgw.yaml') as f:
                        data = yaml.load(f, Loader=BaseLoader)
                        total_energy[i,j+1] = float(data['optical spectrum']['excitations']['energies'][f'{j+1}'])
            else:
                print(f'Iteration {i} has some error in the simulation')
                total_energy[i,0] = 99999

                # Now we add the BSE excitation energy spectra
                for j in range(num_excitations):
                    with open(f'{path}simulations/sim{i}/BSE/molgw.yaml') as f:
                        data = yaml.load(f, Loader=BaseLoader)
                        total_energy[i,j+1] = 99999


        
    # Saving this data
    np.save(f'Energies_{scf}_{basis}_{method}',total_energy)

    return


def get_dft_data():
    # For DFT with BSE we only have one yaml file per simulation

    # Creating datastructure
    total_energy = np.empty((steps, num_excitations+1))

    for i in range(steps):

        # First we do the ground state energy
        print(i)
        with open(f'{path}simulations/sim{i}/DFT_BSE/molgw.yaml') as f:
            data = yaml.load(f, Loader=BaseLoader)

            # Performing a check is the simulation was perforrmed correctly
            if os.path.isfile(f'{path}simulations/sim{i}/DFT_BSE/dynamical_dipole_polarizability.dat'):
                total_energy[i,0] = float(data['scf energy']['total']) * Ha_to_ev

                # Now we add the BSE excitation energy spectra
                for j in range(num_excitations):
                    total_energy[i,j+1] = float(data['optical spectrum']['excitations']['energies'][f'{j+1}'])
            else:
                print(f'Iteration {i} has some error in the simulation')
                total_energy[i,0] = 99999
                # Now we add the BSE excitation energy spectra
                for j in range(num_excitations):
                    total_energy[i,j+1] = 99999
        
    # Saving this data
    np.save(f'Energies_{scf}_{basis}_{method}',total_energy)


    ## Copying the data to the root folder for easy acces
    #src_path = f'Energies_{scf}_{basis}_{method}.npy'
    #dst_path = f'../../Energies_{scf}_{basis}_{method}.npy'
    #shutil.copy(src_path, dst_path)

    return

# Setting constants
Ha_to_ev = 27.2114

# Setting parameters
num_excitations = 40      # How many of the excitations of BSE are included

# Getting simulation settings
with open('Post_processing.yaml') as file:
    data = yaml.load(file, Loader=BaseLoader)
    steps = int(data['Steps'])
    scf = data['Scf']
    basis = data['Basis']
    method = data['Method']
    path = data['Path']
    print('The following simulations settings were found\n')
    print(f'Steps: {steps}, Scf: {scf}, Basis: {basis} and Method: {method}')


if method == 'GW':
    get_gw_data()
else:
    get_dft_data()
