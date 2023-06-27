"""
Created on Sun May 28 

@author: Max Hilt
"""

# Importing libraries
import numpy as np
import argparse
import os
import utils as util

# First we look at out xyz output file
parser = argparse.ArgumentParser()
parser.add_argument("-zmat", dest="zmatfile", required=True, type=str, help="Zmat file we calculate the derivate around")
parser.add_argument("-i", dest="iteration", required=True, type=int, help="What integer of the simulation we are on")
args = parser.parse_args()

filename = args.zmatfile
iteration = args.iteration




def make_derivative_files():
    ''' 
    This function generates the slightly adjusted z mat files such that the derivative can be calculated.
    
    '''


    parameters = len(rlist) + len(alist) + len(dlist)
    for i in range(parameters):
        
        filename1 = f'./Simulations/iteration_{iteration}/par_{i}_low/structure.zmat'
        filename2 = f'./Simulations/iteration_{iteration}/par_{i}_high/structure.zmat'
        filename3 = f'./Simulations/iteration_{iteration}/par_{i}_low/structure.xyz'
        filename4 = f'./Simulations/iteration_{iteration}/par_{i}_high/structure.xyz'
        
        # Making the folder there it goes in
        os.mkdir(f'./Simulations/iteration_{iteration}/par_{i}_low')
        os.mkdir(f'./Simulations/iteration_{iteration}/par_{i}_high')

        if i < len(rlist):
            # We change one of the lengths
            rlist1 = np.copy(rlist)
            rlist1[i] -= delta_distance 
            rlist2 = np.copy(rlist)
            rlist2[i] += delta_distance
            
            util.write_adjusted_zmat(filename1, atomnames, rlist1, alist, dlist)
            util.write_adjusted_zmat(filename2, atomnames, rlist2, alist, dlist)
            util.write_xyz(filename3, atomnames, rconnect, rlist1, aconnect, alist, dconnect, dlist)
            util.write_xyz(filename4, atomnames, rconnect, rlist2, aconnect, alist, dconnect, dlist)
        
        elif i < len(rlist) + len(alist):
            # We change one of the angles
            
            alist1 = np.copy(alist)
            alist[i-len(rlist)] -= delta_angle
            alist2 = np.copy(alist)
            alist[i-len(rlist)] += delta_angle
            
            util.write_adjusted_zmat(filename1, atomnames, rlist, alist1, dlist)
            util.write_adjusted_zmat(filename2, atomnames, rlist, alist2, dlist)
            util.write_xyz(filename3, atomnames, rconnect, rlist, aconnect, alist1, dconnect, dlist)
            util.write_xyz(filename4, atomnames, rconnect, rlist, aconnect, alist2, dconnect, dlist)

        else:
            # We change of the dihedral angles
            
            dlist1 = np.copy(dlist)
            dlist1[i - len(rlist) - len(alist)] -= delta_dihedral
            dlist2 = np.copy(dlist)
            dlist2[i - len(rlist) - len(alist)] += delta_dihedral
                   
            util.write_adjusted_zmat(filename1, atomnames, rlist, alist, dlist1)
            util.write_adjusted_zmat(filename2, atomnames, rlist, alist, dlist2)
            util.write_xyz(filename3, atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist1)
            util.write_xyz(filename4, atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist2)
        



# Settings for gradient descent
delta_distance = 0.01
delta_angle = 0.1
delta_dihedral = 0.1


# We read the current zmat file that we want to look around
atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist = util.readzmat(filename)


# Then we get all the derivative files around it
make_derivative_files()

