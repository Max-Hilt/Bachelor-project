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
import gcutil as gc


def get_co_bondlength(steps):
    '''
    INPUTS:


    OUPUTS:



    '''

    # Creating data_structure
    bond_lengths = np.empty(steps)

    # Collecting the data

    for i in range(steps):
        # Defining the filename of the current iteration
        filename = f'iteration_{i}/CO_{i}.zmat'

        # Gatering the data from that zmat file
        atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist = gc.readzmat(filename)

        # Saving the bondlength to the bond_lengths file that is being created.
        bond_lengths[i] = rlist[0]

        # To give an indication we also print the result
        print(f'Iteration {i} has bondlength {rlist[0]}')


    # Finally when all data is collected we save the results
    np.save(f'CO_bondlengths',bond_lengths)
    return 

parser = argparse.ArgumentParser()
parser.add_argument("-n", dest="number_of_iterations", required=True, type=int, help="The total number of iterations that happened")
args = parser.parse_args()

# Setting parameters 
steps = args.number_of_iterations


# Collecting the bondlength data to one file
get_co_bondlength(steps)