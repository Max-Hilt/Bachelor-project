# -*- coding: utf-8 -*-
"""
Created on Wed May 10 18:46:17 2023

@author: Max Hilt
"""

# Importing libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,mark_inset)


def zoomed_image_plot(data, angstrom, data_string):
    
    '''
    This function makes a figure with a zoomed in subplot inside of the figure
    of the groundstate energy
    
    INPUTS:
        
        data [2D numpy array]: The sorted dataset with all the excitations
        angstrom [1D numpy array]: The range of distances that were simulated
        data_string [String]: The plot title
        
    OUTPUTS:
        
        None
        
    '''
    ground_state = data[:,0]

    
    # Creating the main figure
    fig, ax1 = plt.subplots()

    ax1.plot(angstrom, ground_state, 'x', c='b', mew=2, alpha=0.8, label='Simulations')
    ax1.set_xticks(list(plt.xticks()[0]) + [1.2], [text_obj for text_obj in plt.xticks()[0]] + ['Measurement'])
    ax1.set_xlabel(r'Bond length [$\AA$]')
    ax1.set_ylabel(r'Groundstate energy [eV]')
    ax1.legend(loc=0)


    ax2 = plt.axes([0,0,1,1])
    # Setting the position of the subplot (x,y,width, heigth)
    ip = InsetPosition(ax1, [0.45,0.4,0.5,0.5])
    ax2.set_axes_locator(ip)
    # Mark the region that is in the second plot. loc1, loc2 are corners with connected lines ec is inverse opacity.
    mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

    # Here we plot the data for the inset figure

    ax2.plot(angstrom[150:250], ground_state[150:250], 'x', c='b', mew=2, alpha=0.8, label = 'Simulations')
    ax2.axvline(1.128, color = 'black', label = 'Measurements')

    # And creating the legend
    ax2.legend(loc=0)

    # Some ad hoc tweaks.1
    ax1.set_ylim(-3085, -3050)
    ax1.set_xlim(0.9, 1.5)
    #ax2.set_yticks(np.arange(0,2,0.4))
    #ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
    #ax2.tick_params(axis='x', which='major', pad=8)

    plt.title(data_string)
    plt.show()
    
    return


def correct_data(data, crossing_treshold=0.01, degenerate_treshold=0.01):
    '''
    This function attempts to correctly identify and label the excited state energies after they have crossed eachother. 
    
    ! Known limitations if one of the two lines is nondegenerate because two lines lie ontop of eachother the code does not work properly
    
    INPUTS:
        
        data [2D numpy array]: The incorrectly labled data
        crossing_treshold [float]: Within what distance we see two lines as having crossed
        degenerate_treshold [float]: Used in an attempt to get rid of the degenerate lines by realizing when two lines are ontop of eachother
        
    OUTPUTS:
        
        sorted_data [2D numpy array]: A dataset that correctly identifies crossings of nondegenerate lines
        
    '''
    # First we calculate the energies of the ground state and excited states  
    corr_data = np.empty((np.shape(data)))
    corr_data[:,0] = data[:,0]
    for i in range(np.shape(data)[1]-1):
        corr_data[:,i+1] = data[:,0] + data[:,i+1] 
        
    # Now we must find the places where two lines may have crossed each other
    crossings_list = []
    
    # We loop over all possible columns combinations
    for j in range(np.shape(corr_data)[1]):
        # We must keep track of the crossings in the current comparison
        crossings = []
        for i in range(np.shape(corr_data)[1]):
            
            # To make sure that we only look at interactions once
            if i > j:
            
                # Seek segments where a crossing is possible
                abs_difference = abs(corr_data[:,i] - corr_data[:,j])
                crossing_mask = abs_difference < crossing_treshold
                
                # Now we must identify the single crossing points
                while np.any(crossing_mask):
                    counter = 0    
                    temp = []
                    
                    # First we identify the next occurance of a true.
                    first_index = np.argmax(crossing_mask)
                    
                    # Now we add all the consequative true indices
                    while crossing_mask[first_index + counter] == True:
                        # We make sure to set all values we visited to false
                        crossing_mask[first_index + counter] = False
                        
                        # And we add the correct index to a new temp array
                        temp.append(first_index + counter)
                        counter += 1
                        
                        if counter + first_index >= len(crossing_mask):
                            break
                
                    # Now to find the minimum distance in this array segment
                    additional_index = np.argmin(abs_difference[temp])
                    crossing_point = first_index + additional_index
                    crossings.append(crossing_point)
        
        crossings_list.append(crossings)
                          
    
    # Now that we know the crossing locations we can start getting non degenerate data
    degenerate_mask = abs(data[0,:-1] - data[0,1:]) < degenerate_treshold
    numbers = np.arange(0, len(degenerate_mask), 1, dtype=int)
    starting_indices = numbers[degenerate_mask]

    # We start with a complete copy of the data
    non_degenerate_data = np.copy(corr_data)
    
    
    # Now we get rid of wrong values
    for i, num in enumerate(starting_indices):
        current_column = num
        previous_cross = 0
         
        whileloop = True
        # Now we repeat the whileloop until the column is all changed.
        while_count = 0
        while whileloop:
            while_count += 1
            if while_count>10000:
                print('Stuck in whileloop 1')
            
            # We must make sure the crossing list is nonempty
            if crossings_list[current_column] == []:
                non_degenerate_data[previous_cross:, current_column] = np.ones(np.shape(non_degenerate_data[previous_cross:, current_column]))
                whileloop = False

            
            for i, cross_index in enumerate(crossings_list[current_column]):
                
                if cross_index > previous_cross:
                    # We set the nondegenerate entries equal to one
                    non_degenerate_data[previous_cross:cross_index,current_column] = np.ones(np.shape(non_degenerate_data[previous_cross:cross_index,current_column]))
                                
                    if current_column == 0:
                        current_column = 1
                        
                    elif current_column == (len(crossings_list)-1):
                        current_column = current_column - 1
                        
                    elif abs(data[cross_index + 1,current_column - 1] - data[cross_index + 1,current_column]) < abs(data[cross_index + 1,current_column + 1] - data[cross_index + 1,current_column]):
                        current_column = current_column - 1
                        
                    else:
                        current_column = current_column + 1
                        
                    previous_cross = cross_index
                    break
                
                # If we have achieved the last value and apperently not switched were done
                if cross_index == crossings_list[current_column][-1]:
                    non_degenerate_data[previous_cross:, current_column] = np.ones(np.shape(non_degenerate_data[previous_cross:, current_column]))
                                    
                    whileloop = False
    

    # Finally from the identifier matrix we can construct our correct data
    non_degenerate_corr_data = np.zeros(np.shape(data))
    for i in range(np.shape(data)[0]):
        counter = 0
        for j in range(np.shape(data)[1]):
            if non_degenerate_data[i,j] != 1:
                non_degenerate_corr_data[i,counter] = non_degenerate_data[i,j]
                counter += 1
    THISLIST = non_degenerate_corr_data
    # Now we have our non degenerate dataset so we redo the intersection calculations with this new and better dataset
    crossings_list = []
    
    # We loop over all possible columns combinations
    for j in range(np.shape(non_degenerate_corr_data)[1]):
        # We must keep track of the crossings in the current comparison
        crossings = []
        for i in range(np.shape(non_degenerate_corr_data)[1]):
            
            # To make sure that we only look at interactions once
            if i > j:
            
                # Seek segments where a crossing is possible
                abs_difference = abs(non_degenerate_corr_data[:,i] - non_degenerate_corr_data[:,j])
                crossing_mask = abs_difference < crossing_treshold
                
                # Now we must identify the single crossing points
                while_count = 0
                while np.any(crossing_mask):
                    
                    counter = 0    
                    temp = []
                    
                    # First we identify the next occurance of a true.
                    first_index = np.argmax(crossing_mask)
                    
                    # Now we add all the consequative true indices
                    while crossing_mask[first_index + counter] == True:
                        while_count += 1
                        if while_count > 10000:
                            print('Stuck in while loop')
                        
                        
                        # We make sure to set all values we visited to false
                        crossing_mask[first_index + counter] = False
                        
                        # And we add the correct index to a new temp array
                        temp.append(first_index + counter)
                        counter += 1
                        
                        if counter + first_index >= len(crossing_mask):
                            break
                
                    # Now to find the minimum distance in this array segment
                    additional_index = np.argmin(abs_difference[temp])
                    crossing_point = first_index + additional_index
                    crossings.append(crossing_point)
        
        crossings_list.append(crossings)

    num_rows = np.shape(non_degenerate_corr_data)[0]
    num_cols = np.shape(non_degenerate_corr_data)[1]
    # Here we create the identifier matrix that can be seen as the map to assign correct values
    identifier_matrix = np.arange(0, num_cols)[np.newaxis, :].repeat(num_rows, axis=0)
    
    
    wheres_what = np.arange(0, np.shape(identifier_matrix)[1], 1, dtype=int)
    
    # Now we must correctly swap the identifier matrix list to account of all the crossings
    # And we must make sure that we also keep the order of the crossings correct
    for k in range(np.shape(non_degenerate_corr_data)[0]):
        
        for i, crossings in enumerate(crossings_list):
            for j in range(len(crossings)):
                # Check if the crossing must be performed now
                if crossings[j] == k:
                    
                    if i < len(crossings_list)-1:   # Dit even doordenken
                        
                        # Performing the swap
                        temp = np.copy(identifier_matrix[k:,wheres_what[i]])
                        identifier_matrix[k:,wheres_what[i]] = np.copy(identifier_matrix[k:,wheres_what[i+1]])
                        identifier_matrix[k:,wheres_what[i+1]] = temp
                        
                        tempel = wheres_what[i] 
                        wheres_what[i] = wheres_what[i+1]
                        wheres_what[i+1] = tempel
                        # And now also the crossingslist must be modified
 
    sorted_data = np.zeros(np.shape(non_degenerate_corr_data))
    
    for i in range(np.shape(sorted_data)[0]):
        for j in range(np.shape(sorted_data)[1]):
            sorted_data[i,j] = non_degenerate_corr_data[i,identifier_matrix[i,j]]
       
    return sorted_data


def one_excitation(data, angstrom, data_string, min_ground = 0, min_excited = 0, excited_num = 1):

    '''
    This function makas a plot of the Ground state and Selected excited state
    and also shows their minima if they are provided.
    
    INPUTS:
        
        data [2D numpy array]: The sorted dataset with all the excitations
        angstrom [1D numpy array]: The range of distances that were simulated
        data_string [String]: The plot title
        excited_num [Int]: What excitation is plotted
        min_ground [Float]: Location minima of the ground_state
        min_excited [Float]: Location minima of the excited state.
        
    OUTPUTS:
        
        None
        
    '''    

    # Selecting the correct groundstate and excited state energies
    ground_state_energy = data[:,0]
    excited_energy = data[:,excited_num]
    
    # Plotting the data
    plt.figure()
    plt.title(data_string)
    plt.plot(angstrom, ground_state_energy,'.', label='Ground state energy')
    plt.plot(angstrom, excited_energy,'.', label=f'Excited state {excited_num} energy')
    plt.plot(angstrom, data[:,2], color = 'black', alpha = 0.4, label='higher order 2')
    plt.plot(angstrom, data[:,3], color = 'black', alpha= 0.4, label='higher order 3')
    
    if min_ground != 0:
        plt.axvline(min_ground,ls = '--', label = 'Minimum ground state')
    if min_excited != 0:
        plt.axvline(min_excited,ls = '--', color='C1', label = 'Minimum excited state')
    plt.legend()
    
    return

def one_excitation_report(data, angstrom, data_string, min_ground = 0, min_excited = 0, excited_num = 1):

    '''
    This function is an almost identical copy of the one_excitation function
    Only slight modifications have been made to make the plot desired to be put
    T
    
    INPUTS:
        
        data [2D numpy array]: The sorted dataset with all the excitations
        angstrom [1D numpy array]: The range of distances that were simulated
        data_string [String]: The plot title
        excited_num [Int]: What excitation is plotted
        min_ground [Float]: Location minima of the ground_state
        min_excited [Float]: Location minima of the excited state.
        
    OUTPUTS:
        
        None
        
    '''    

    # Selecting the correct groundstate and excited state energies
    ground_state_energy = data[:,0]
    excited_energy = data[:,excited_num]
    
    # Plotting the data
    plt.figure()
    plt.title(data_string)
    plt.plot(angstrom, ground_state_energy,'.', label='Ground-state')
    plt.plot(angstrom, excited_energy,'.', label=f'First excited-state')
    plt.plot(angstrom, data[:,2], color = 'black', alpha = 0.4, label='Second excited-state')
    plt.plot(angstrom, data[:,3], color = 'black', alpha= 0.6, label='Third excited-state')
    
    plt.xlabel('Bond length [$\AA$]')
    plt.ylabel('Energy [eV]')
    plt.legend(bbox_to_anchor=(1, 1), loc=1, borderaxespad=0)
    
    if min_ground != 0:
        plt.axvline(min_ground,ls = '--', label = 'Minimum ground state')
    if min_excited != 0:
        plt.axvline(min_excited,ls = '--', color='C1', label = 'Minimum excited state')
    
    return

def five_excitation(data, angstrom):

    '''
    This function plots the first five excitations together with the ground state energy.
    
    INPUTS:
        
        data [2D numpy array]: The sorted dataset with all the excitations
        angstrom [1D numpy array]: The range of distances that were simulated
        
    OUTPUTS:
        
        None
        
    '''      

    ground_state_energy = data[:,0]

    # Creating the figure
    plt.figure()
    plt.title('First five excited states')
    plt.plot(angstrom, ground_state_energy,'.', label='Groundstate energy')
    for i in range(5):
        excited_energy = data[:,i+1]
        plt.plot(angstrom, excited_energy, '.', label=f'Excitedstate {i+1} energy')
    plt.legend()
    
    return

def all_plot_opacity(data, angstrom):

    '''
    This function plots all the excitations with reduced opacity so degeneracies can more easily
    be seen.
    
    INPUTS:
        
        data [2D numpy array]: The unsorted dataset with all the excitations
        angstrom [1D numpy array]: The range of distances that were simulated
        
    OUTPUTS:
        
        None
        
    '''   
    # Creating the figure
    plt.figure()
    plt.title('All excited states')
    for i in range(np.shape(data)[1]):
        if i == 0:
            excited_energy = data[:,i]
        else:
            excited_energy = data[:,i] + data[:,0]
        plt.plot(angstrom, excited_energy, color = 'C0',linewidth = 3, alpha = 0.4, label=f'Excited energy {i}')
    plt.legend()
    
    return

def find_minima(data, angstrom, excited_num = 1):
    ground_state_energy = data[:,0]
    excited_energy = data[:,excited_num]
    
    # Finding the minima in the data
    min_ground = np.argmin(ground_state_energy)
    min_excited = np.argmin(excited_energy)
    
    dr = angstrom[1] - angstrom[0]
    
    force = (excited_energy[min_ground+1]-excited_energy[min_ground-1])/(2*dr)
    
    return angstrom[min_ground], angstrom[min_excited], force




def get_sim_data(method, scf, basis):
    
    '''
    This function provides the data that is selected, without asking user input
    2 is subtracted from the basis.
    
    INPUTS:
        
        method [int]: DFT or GW
        scf [int]: BLYP, B3LYP or BHLYP
        basis [int]: pVDZ, pVTZ, pVQZ or pV5Z
        
    OUTPUTS:
        
        data [2D numpy array]: The unsorted dataset with all the excitations
        
    '''   
    
    # Different simulation settings
    method_list = ['DFT', 'GW']
    scf_list = ['BLYP','B3LYP','BHLYP']
    basis_list = ['cc-pVDZ','cc-pVTZ','cc-pVQZ','cc-pV5Z']

    data_string = f'Energies_{scf_list[scf]}_{basis_list[basis-2]}_{method_list[method]}.npy'
    data = np.load(data_string)

    return data

def user_input_plot():
    
    '''
    Function that uses the console to ask what dataset should be loaded.
    
    INPUTS:
        
        None
        
    OUTPUTS:
        
        data [2D numpy array]: The unsorted dataset with all the excitations
        data_string: The corresponding string of what datafile was loaded
        
    '''   
    # Different simulation settings
    method_list = ['DFT', 'GW']
    scf_list = ['BLYP','B3LYP','BHLYP']
    basis_list = ['cc-pVDZ','cc-pVTZ','cc-pVQZ','cc-pV5Z']
    
    print('We are going to make one figure')
    method = int(input('What method do you want to visualize? [0,1] '))
    print(f'{method_list[method]} was chosen.')
    scf = int(input('What scf do you want to use? [0,1,2] '))
    print(f'{scf_list[scf]} was chosen.')
    basis = int(input('What basis do you want to use? [2,3,4,5] '))-2
    print(f'{basis_list[basis]} was chosen.')
    #excited_num = int(input('What excited state do you want to use? '))
    
    data_string = f'Energies_{scf_list[scf]}_{basis_list[basis]}_{method_list[method]}.npy'
    print(data_string)
    data = np.load(data_string)
    
    return data, data_string
    

def ground_state():
    
    '''
    This function is used to make a plot of the groundstate of the CO molecule.
    It plots the bondlength on the horizontal axis and the groundstate energy on the
    vertical axis.
    
    INPUTS
        none
    OUPUTS
        none
    
    
    '''
    # First we import the data
    data_string = f'Energies_B3LYP_cc-pV5Z_GW.npy'
    data = np.load(data_string)
    
    # Then we define the bondlength
    bondlength = np.linspace(1,1.4,400)
    
    # Than we plot the data
    plt.figure()
    plt.title('CO groundstate simulation using B3LYP and cc-pV5Z')
    plt.plot(bondlength, data[:,0], '.')
    plt.xlabel(r'Bond length [$\AA$]')
    plt.ylabel('Energy [eV]')














