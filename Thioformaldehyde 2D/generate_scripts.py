"""
Created on Sun May 28 

@author: Max Hilt
"""

# Importing libraries
import numpy as np
import argparse
import os
import utils as util


def generate_files():
    ''' 
    This function generates the slightly adjusted z mat files such that the derivative can be calculated.
    
    '''

    # We create two lists of the parameters where we create a grid of
    variable_1_list = np.linspace(min_variable1, max_variable1, steps)
    variable_2_list = np.linspace(min_variable2, max_variable2, steps)


    parameters = len(rlist) + len(alist) + len(dlist)
    for i in range(steps):
        for j in range(steps):
            
            # First we define the filenames of the files we will be creating
            filename1 = f'./Simulations/location_{i + steps*j}/structure.zmat'
            filename2 = f'./Simulations/location_{i + steps*j}/structure.xyz'
            
            # Then we actually create the folder where the files will be stored
            os.mkdir(f'./Simulations/location_{i + steps*j}')

            # For each iteration we copy the current lists that have been created since we don't want any permanent changes on the original list.
            # This is inefficient code but still fast enough that it is not yet a priority to improve.
            rlist1 = np.copy(rlist)
            alist1 = np.copy(alist)
            dlist1 = np.copy(dlist)

            # Now we look at what value must be altered and write the corresponding files
            if variable1 < len(rlist):
                rlist1[variable1] = variable_1_list[i]
            elif variable1 < len(rlist) + len(alist):
                alist1[variable1-len(rlist)] = variable_1_list[i]
            else:
                dlist1[variable1 - len(rlist) - len(alist)] = variable_1_list[i]

            if variable2 < len(rlist):
                rlist1[variable2] = variable_2_list[j]
            elif variable2 < len(rlist) + len(alist):
                alist1[variable2-len(rlist)] = variable_2_list[j]
            else:
                dlist1[variable2 - len(rlist) - len(alist)] = variable_2_list[j]     

            # Now all the variables are as we want them to be so we can write the files,

            util.write_adjusted_zmat(filename1, atomnames, rlist1, alist1, dlist1)
            util.write_xyz(filename2, atomnames, rconnect, rlist1, aconnect, alist1, dconnect, dlist1)
            


# First we must let the user select the grid that is used
variable1 = int(input('What is the first variable do you want to change? [0-5]  '))
print(f'Variable {variable1} was selected')
min_variable1 = float(input('What is the min value you want?   '))
print(f'The minimum of variable 1 is {min_variable1}')
max_variable1 = float(input('What is the max value you want?   '))
print(f'The minimum of variable 1 is {max_variable1}')
variable2 = int(input('What is the second variable do you want to change? [0-5]  '))
print(f'Variable {variable2} was selected')
min_variable2 = float(input('What is the min value you want?   '))
print(f'The minimum of variable 1 is {min_variable2}')
max_variable2 = float(input('What is the max value you want?   '))
print(f'The minimum of variable 1 is {max_variable2}')
steps = int(input('How many steps should the grid be divided   '))
print(f'The grid is divided in {steps} steps')


# We read the current zmat file such that we gat the variables that we won't change
filename = 'Thio.zmat'
atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist = util.readzmat(filename)


# Then we get all the derivative files around it
os.mkdir('./Simulations')
generate_files()

# Finally we save a file so we know what the simulation settings were
f = open(f'./Simulation_settings.txt', 'w+')
f.write('##### This simulation had the following settings #####\n\n')
f.write(f'Variable 1 {variable1}\n')
f.write(f'Minlength: {min_variable1}\n')
f.write(f'Maxlength: {max_variable1}\n')
f.write(f'Variable 2 {variable2}\n')
f.write(f'Minlength: {min_variable2}\n')
f.write(f'Maxlength: {max_variable2}\n')
f.write(f'Steps: {steps}\n')


print('Files have been generated')
