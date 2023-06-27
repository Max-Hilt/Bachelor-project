"""
Created on Sun May 28 

@author: Max Hilt
"""

# Importing libraries
import numpy as np
import utils as util
import argparse
import os
from yaml.loader import BaseLoader
import yaml




def calculate_derivative(left,right,delta_step):
    '''
    INPUTS:
        left : float
            This is the value one delta_step to the left of the point we are calculating the derivative of.
        right : float
            This is the value one delta_step to the right of the point we are calculating the derivative of.
        delta_step : float
            This is half the difference between the values where left and right were evaluated.

    OUTPUTS:
        derivative : float
            The derivative at the point in between the points left and right.

    
    This function calculates the central derivative given the point to the left, the point to the right and the delta distance performed by the step
    '''

    derivative = (right - left)/(2*delta_step)

    return derivative


def get_data(path):
    '''
    INPUTS:
        path : string
            This is the path to the YAML file containing the simulation results

    OUTPUTS:
        energies : numpy array
            The energies of the groundstate plus num_excitations excited states in [eV]

    This function collects the data saved to the YAML files after simulation.
    '''
    energies = np.zeros(num_excitations + 1)

    with open(f'{path}/molgw_gw.yaml') as f:
        data = yaml.load(f, Loader=BaseLoader)
        energies[0] = float(data['scf energy']['total']) * Ha_to_ev

    with open(f'{path}/molgw_bse.yaml') as f:
        data = yaml.load(f, Loader=BaseLoader)
        for j in range(num_excitations):
            # Here we add the excitation energy of the YAML file immediately to the ground state energy.
            energies[j+1] = float(data['optical spectrum']['excitations']['energies'][f'{j+1}']) + energies[0]

    return energies


def write_new_zmat(derivatives):


    # We start with reading the previous zmat
    atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist = util.readzmat(current_folder + f'Thio_{current_index}.zmat')

    # Than we make our adjustments
    print('This is specially written for Thioformaldehyde')
        
    filename1 = f'Thio.zmat'
    filename2 = f'Thio.xyz'
    
    # First we make the relative corrections

    # In the distances
    for i in range(len(rlist)):
        rlist[i] = rlist[i] - derivatives[i]
    
    # In the angles
    for i in range(len(alist)):
        alist[i] = alist[i] - derivatives[i + len(rlist)]
    
    # In the dihedral angles
    for i in range(len(dlist)):
        dlist[i] = dlist[i] - derivatives[i + len(rlist) + len(alist)]

    # Then we write our two new files
    util.write_adjusted_zmat(filename1, atomnames, rlist, alist, dlist)
    util.write_xyz(filename2, atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist)
    
    return

def stop_simulation(derivatives):
    '''
    INPUTS

    OUTPUTS

    This function creates a file called simulation.stop if the derivative has become small enough, once the bash file recognizes that this file has been created the simulation is stopped.

    '''
    # We set the threshold value here
    threshold = 0.05
    filename = "simulation.stop"

    # Finally if the length is under our threshold we are satisfied that we have converged enough.
    if all(abs(derivatives/learning_rate) < threshold) and (abs(derivatives[0]/learning_rate) < threshold):
        with open(filename, 'w') as f: 
            f.write(f'Simulation stopped on iteration {current_index} excited state {excitedstate} was simulated')  
            print(f'The simulation was stopped because all derivatives are small enough. We have done {current_index} iterations')   
    return
        

# We require the calling to also pass through the current folder we are in.
parser = argparse.ArgumentParser()
parser.add_argument("-id", dest="current_index", required=True, type=int, help="The index of the current iteration")
parser.add_argument("-ex", dest="excited_state_num", required=True, type=int, help="What excited state is being optimized")
args = parser.parse_args()

current_index = args.current_index
excitedstate = args.excited_state_num
current_folder = f'Simulations/iteration_{current_index}/'

# Defining parameters & constants
Ha_to_ev = 27.2114
num_excitations = excitedstate + 5
dimensions = 6

# Settings for gradient descent
delta_distance = 0.01
delta_angle = 0.1
delta_dihedral = 0.1
learning_rate = 0.005


# We simply go over all the different parameters that can be changed and calculate their derivatives.
derivatives = np.zeros(dimensions)

for i in range(dimensions):
    left_energies = get_data(f'{current_folder}par_{i}_low')
    right_energies = get_data(f'{current_folder}par_{i}_high')

    # Making sure the correct delta value is used
    if i > (dimensions-3)/3 +1:
        delta_distance = delta_angle

    derivatives[i] = calculate_derivative(left_energies[excitedstate], right_energies[excitedstate], delta_distance)


print('These are the derivatives: \n')
print(derivatives)
print('\n\n')

# Now we want to 'normalize'each direction such that the step length is limited.
derivatives = learning_rate*np.copy(derivatives)

print('These are the changes: \n')
print(derivatives)
print('\n\n')

# Final step is to calculate the new zmat file.
write_new_zmat(derivatives)

# And we check if we have already converged enough to stop the simulation
stop_simulation(derivatives)