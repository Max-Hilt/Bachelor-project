"""
Created on Sun May 28 13:31:16 2023

@author: Max Hilt
"""

# Importing libraries
import numpy as np
import gcutil as gc
import argparse
import os

# First we look at out xyz output file
parser = argparse.ArgumentParser()
parser.add_argument("-zmat", dest="zmatfile", required=True, type=str, help="Zmat file we calculate the derivate around")
parser.add_argument("-i", dest="iteration", required=True, type=int, help="What integer of the simulation we are on")
args = parser.parse_args()

filename = args.zmatfile
iteration = args.iteration

print(f'we have {iteration} is type {type(iteration)}')

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




def make_derivative_files():
    ''' 
    This function generates the slightly adjusted z mat files such that the derivative can be calculated.
    
    '''


    parameters = len(rlist) + len(alist) + len(dlist)
    for i in range(parameters):
        
        filename1 = f'./iteration_{iteration}/par_{i}_low/structure.zmat'
        filename2 = f'./iteration_{iteration}/par_{i}_high/structure.zmat'
        filename3 = f'./iteration_{iteration}/par_{i}_low/structure.xyz'
        filename4 = f'./iteration_{iteration}/par_{i}_high/structure.xyz'
        
        # Making the folder there it goes in
        os.mkdir(f'./iteration_{iteration}/par_{i}_low')
        os.mkdir(f'./iteration_{iteration}/par_{i}_high')

        if i < len(rlist):
            # We change one of the lengths
            rlist1 = np.copy(rlist)
            rlist1[i] -= delta_distance 
            rlist2 = np.copy(rlist)
            rlist2[i] += delta_distance
            
            write_adjusted_zmat(filename1, atomnames, rlist1, alist, dlist)
            write_adjusted_zmat(filename2, atomnames, rlist2, alist, dlist)
            write_xyz(filename3, atomnames, rconnect, rlist1, aconnect, alist, dconnect, dlist)
            write_xyz(filename4, atomnames, rconnect, rlist2, aconnect, alist, dconnect, dlist)
        
        elif i < len(rlist) + len(alist):
            # We change one of the angles
            
            alist1 = np.copy(alist)
            alist[i-len(rlist)] -= delta_angle
            alist2 = np.copy(alist)
            alist[i-len(rlist)] += delta_angle
            
            write_adjusted_zmat(filename1, atomnames, rlist, alist1, dlist)
            write_adjusted_zmat(filename2, atomnames, rlist, alist2, dlist)
            write_xyz(filename3, atomnames, rconnect, rlist, aconnect, alist1, dconnect, dlist)
            write_xyz(filename4, atomnames, rconnect, rlist, aconnect, alist2, dconnect, dlist)

        else:
            # We change of the dihedral angles
            
            dlist1 = np.copy(dlist)
            dlist1[i - len(rlist) - len(alist)] -= delta_dihedral
            dlist2 = np.copy(dlist)
            dlist2[i - len(rlist) - len(alist)] += delta_dihedral
                   
            write_adjusted_zmat(filename1, atomnames, rlist, alist, dlist1)
            write_adjusted_zmat(filename2, atomnames, rlist, alist, dlist2)
            write_xyz(filename3, atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist1)
            write_xyz(filename4, atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist2)
        

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


# Settings for gradient descent
delta_distance = 0.01
delta_angle = 0.1
delta_dihedral = 0.1


# We read the current zmat file that we want to look around
atomnames, rconnect, rlist, aconnect, alist, dconnect, dlist = readzmat(filename)


# Then we get all the derivative files around it
make_derivative_files()


