# Written by Max Hilt

# Importing libraries
import numpy as np
import gcutil as gc
import argparse
import os
from yaml.loader import BaseLoader
import yaml

def replace_vars(vlist, variables):
    """ Replaces a list of variable names (vlist) with their values
        from a dictionary (variables).
    """
    for i, v in enumerate(vlist):
        if v in variables:
            vlist[i] = variables[v]
        else:
            try:
                # assume the "variable" is a number
                vlist[i] = float(v)
            except:
                print("Problem with entry " + str(v))

def write_xyz(filename, atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist):
    """Prints out an xyz file from a decomposed z-matrix"""
    npart = len(atomnames)

    with open(filename, 'w') as f:
        f.write(f'{npart}'+ '\n')
        f.write('CO' + '\n')
        
        # put the first atom at the origin
        xyzarr = np.zeros([npart, 3])
        if (npart > 1):
            # second atom at [r01, 0, 0]
            xyzarr[1] = [rlist[0], 0.0, 0.0]

        if (npart > 2):
            # third atom in the xy-plane
            # such that the angle a012 is correct 
            i = rconnect[1] - 1
            j = aconnect[0] - 1
            r = rlist[1]
            theta = alist[0] * np.pi / 180.0
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            a_i = xyzarr[i]
            b_ij = xyzarr[j] - xyzarr[i]
            if (b_ij[0] < 0):
                x = a_i[0] - x
                y = a_i[1] - y
            else:
                x = a_i[0] + x
                y = a_i[1] + y
            xyzarr[2] = [x, y, 0.0]

        for n in range(3, npart):
            # back-compute the xyz coordinates
            # from the positions of the last three atoms
            r = rlist[n-1]
            theta = alist[n-2] * np.pi / 180.0
            phi = dlist[n-3] * np.pi / 180.0
            
            sinTheta = np.sin(theta)
            cosTheta = np.cos(theta)
            sinPhi = np.sin(phi)
            cosPhi = np.cos(phi)

            x = r * cosTheta
            y = r * cosPhi * sinTheta
            z = r * sinPhi * sinTheta
            
            i = rconnect[n-1] - 1
            j = aconnect[n-2] - 1
            k = dconnect[n-3] - 1
            a = xyzarr[k]
            b = xyzarr[j]
            c = xyzarr[i]
            
            ab = b - a
            bc = c - b
            bc = bc / np.linalg.norm(bc)
            nv = np.cross(ab, bc)
            nv = nv / np.linalg.norm(nv)
            ncbc = np.cross(nv, bc)
            
            new_x = c[0] - bc[0] * x + ncbc[0] * y + nv[0] * z
            new_y = c[1] - bc[1] * x + ncbc[1] * y + nv[1] * z
            new_z = c[2] - bc[2] * x + ncbc[2] * y + nv[2] * z
            xyzarr[n] = [new_x, new_y, new_z]
                
        # print results
        for i in range(npart):
            f.write('{:<4s}\t{:>11.5f}\t{:>11.5f}\t{:>11.5f}\n'.format(atomnames[i], xyzarr[i][0], xyzarr[i][1], xyzarr[i][2]))

def write_adjusted_zmat(filename, atomnames, rlist, alist, dlist, rvar=False, avar=False, dvar=False):
    """Prints a z-matrix from xyz coordinates, distances, and atomnames,
       optionally with the coordinate values replaced with variables.
    """
    
    with open(filename, 'w') as f:
        npart = len(rlist)+1
        ncoord = 3
        if npart > 0:
            # Write the first atom
            f.write(atomnames[0] + '\n')
            
            if npart > 1:
                # and the second, with distance from first
                n = atomnames[1]
                if (rvar):
                    r = 'R1'
                else:
                    r = '{:>11.5f}'.format(rlist[0])
                f.write('{:<3s} {:>4d}  {:11s}'.format(n, 1, r) + '\n')
                
                if npart > 2:
                    n = atomnames[2]
                    
                    if (rvar):
                        r = 'R2'
                    else:
                        r = '{:>11.5f}'.format(rlist[1])
                    
                    if (avar):
                        t = 'A1'
                    else:
                        t = '{:>11.5f}'.format(alist[0])
    
                    f.write('{:<3s} {:>4d}  {:11s} {:>4d}  {:11s}'.format(n, 1, r, 2, t) + '\n')
                    
                    if npart > 3:
                        for i in range(3, npart):
                            n = atomnames[i]
    

                            if (rvar):
                                r = 'R{:<4d}'.format(i)
                            else:
                                r = '{:>11.5f}'.format(rlist[i-1])
    
                            if (avar):
                                t = 'A{:<4d}'.format(i-1)
                            else:
                                t = '{:>11.5f}'.format(alist[i-2])
                            
                            if (dvar):
                                d = 'D{:<4d}'.format(i-2)
                            else:
                                d = '{:>11.5f}'.format(dlist[i-3])
                            f.write('{:3s} {:>4d}  {:11s} {:>4d}  {:11s} {:>4d}  {:11s}'.format(n, i-2, r, i-1, t, i, d) + '\n')

def readzmat(filename):
    """ Reads in a z-matrix in standard format,
        returning a list of atoms and coordinates.
    """
    zmatf = open(filename, 'r')
    atomnames = []
    rconnect = []  # bond connectivity
    rlist = []     # list of bond length values
    aconnect = []  # angle connectivity
    alist = []     # list of bond angle values
    dconnect = []  # dihedral connectivity
    dlist = []     # list of dihedral values
    variables = {} # dictionary of named variables
    
    if not zmatf.closed:
        for line in zmatf:
            words = line.split()
            eqwords = line.split('=')
            
            if len(eqwords) > 1:
                # named variable found 
                varname = str(eqwords[0]).strip()
                try:
                    varval  = float(eqwords[1])
                    variables[varname] = varval
                except:
                    print("Invalid variable definition: " + line)
            
            else:
                # no variable, just a number
                # valid line has form
                # atomname index1 bond_length index2 bond_angle index3 dihedral
                if len(words) > 0:
                    atomnames.append(words[0])
                if len(words) > 1:
                    rconnect.append(int(words[1]))
                if len(words) > 2:
                    rlist.append(words[2])
                if len(words) > 3:
                    aconnect.append(int(words[3]))
                if len(words) > 4:
                    alist.append(words[4])
                if len(words) > 5:
                    dconnect.append(int(words[5]))
                if len(words) > 6:
                    dlist.append(words[6])
    
    # replace named variables with their values
    replace_vars(rlist, variables)
    replace_vars(alist, variables)
    replace_vars(dlist, variables)
    
    zmatf.close()

    return (atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist) 

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
    atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist = readzmat(current_folder + f'CO_{current_index}.zmat')

    # Than we make our adjustments
    print('This is specially written for CO')
        
    filename1 = f'CO.zmat'
    filename2 = f'CO.xyz'
    
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
    write_adjusted_zmat(filename1, atomnames, rlist, alist, dlist)
    write_xyz(filename2, atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist)
    
    return
        

# We require the calling to also pass through the current folder we are in.
parser = argparse.ArgumentParser()
parser.add_argument("-id", dest="current_index", required=True, type=int, help="The index of the current iteration")
parser.add_argument("-ex", dest="excited_state_num", required=True, type=int, help="What excited state is being optimized")
args = parser.parse_args()

current_index = args.current_index
excitedstate = args.excited_state_num
current_folder = f'iteration_{current_index}/'

# Defining parameters & constants
Ha_to_ev = 27.2114
num_excitations = excitedstate + 5
dimensions = 1

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
    print(f'Excited state is {excitedstate}')
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