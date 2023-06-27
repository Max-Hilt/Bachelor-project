# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 10:23:42 2023

@author: Max Hilt
"""

# Importing libraries
import matplotlib.pyplot as plt
import numpy as np


def thio_groundstate():
    '''
    INPUTS
        none
    OUPUTS
        none

    This function plots all the barplots for the ground state of thioformaldehyde
    '''
    # The simulation data
    
    simulation_CS = [1.619, 1.60797,1.6125]
    method_CS = ['CC3', '$G_0W_0$, B3LYP','Expt']
    error_CS = [0,0,0.0015] 
    simulation_CS_plot = simulation_CS[1:] - np.ones(len(simulation_CS[1:]))*simulation_CS[0]
    
    
    simulation_CH = [1.083, 1.08541, 1.0945]
    method_CH = ['CC3','$G_0W_0$, B3LYP','Expt']
    error_CH = [0,0,0.0015]
    simulation_CH_plot = simulation_CH[1:] - np.ones(len(simulation_CH[1:]))*simulation_CH[0]
     
    simulation_HCH = [116.1, 120.00745826655526, 116.55] 
    method_HCH = ['CC3', '$G_0W_0$, B3LYP', 'Expt']
    error_HCH = [0, 0, 0.35]
    simulation_HCH_plot = simulation_HCH[1:] - np.ones(len(simulation_HCH[1:]))*simulation_HCH[0]
    
    print(simulation_HCH_plot)
    
    simulation_oop_cs = [0.000, 0.0, 0.000]
    method_oop_cs = ['CC3', '$G_0W_0$, B3LYP', 'Expt']
    error_oop_cs = [0, 0, 0]
    visualizing_aid = [0.03,0.03]
    simulation_oop_cs_plot = simulation_oop_cs[1:] - np.ones(len(simulation_oop_cs[1:]))*simulation_oop_cs[0] + visualizing_aid
    
    # Plotting the data  
    # Set the width of each bar
    bar_width = 0.4
    
    # Create the combined bar graph
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    
    # Plotting the C-S bond length
    ax1.bar(np.arange(len(method_CS[1:])), simulation_CS_plot, width=bar_width, color='C0', yerr=error_CS[1:])
    ax1.set_ylabel('Difference with CC3 [$\AA$]')
    ax1.set_xlabel('Method')
    ax1.set_title('Thioformaldehyde ground-state')
    
    # Plotting the C-H1 bond length
    ax1.bar(np.arange(len(method_CH[1:])) + bar_width, simulation_CH_plot, width=bar_width, color='C1', yerr=error_CH[1:])
    
    # Plotting the H-C-H angle
    ax2.bar(np.arange(len(method_HCH[1:])) + 2*bar_width, simulation_HCH_plot, width=bar_width, color='C2', yerr=error_HCH[1:])
    ax2.set_ylabel('Difference with CC3 [Degrees]')
    
    # Plotting the C-S out of plane angle
    ax2.bar(np.arange(len(method_oop_cs[1:])) + 3*bar_width, simulation_oop_cs_plot, width=bar_width, color='C3', yerr=error_oop_cs[1:])
    
    # Set the x-axis tick positions and labels
    tick_positions = np.arange(len(method_CS[1:])) + bar_width*1.5
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(method_CS[1:])
    
    # Create a legend
    ax1.legend(['$C-S$ bond length', '$C-H_1$ bond length'], loc='upper left')
    ax2.legend([r'$\angle H_1$C$H_2$', '$C-S$ out of plane angle'], loc='upper right')
    
    # Find the maximum absolute value among the data points
    max_abs_value = max(max(abs(simulation_CS_plot)), max(abs(simulation_CH_plot)), max(abs(simulation_HCH_plot)), max(abs(simulation_oop_cs_plot)))
    
    # Set the y-limits of both axes
    y_limit_1 = 0.006 * max_abs_value  # Adjust the factor (1.2) as needed to control the y-axis range
    y_limit_2 = 1.3 * max_abs_value
    ax1.set_ylim(-y_limit_1*0.5, y_limit_1)
    ax2.set_ylim(-y_limit_2*0.5, y_limit_2)
    
    # Display the plot
    plt.tight_layout()
    plt.show()
    
    
    
    # Two seperate plots
    
    # Plot for C-S and C-H #####
    
    plt.figure()
    plt.bar(np.arange(len(method_CS[1:])), simulation_CS_plot, width=bar_width, color='C0', yerr=error_CS[1:])
    
    # Plotting the C-H1 bond length
    plt.bar(np.arange(len(method_CH[1:])) + bar_width, simulation_CH_plot, width=bar_width, color='C1', yerr=error_CH[1:])
    
    # Correcting the names of both sets of bars
    tick_positions = np.arange(len(method_CS[1:])) + bar_width*1
    plt.xticks(tick_positions, method_CS[1:])
    
    plt.ylabel('Difference with CC3 [$\AA$]')
    plt.xlabel('Method')
    plt.title('Thioformaldehyde ground-state')
    plt.legend(['$C-S$ bond length', '$C-H_1$ bond length'], loc='upper right')
    plt.tight_layout()
    
    
    
    # Plot for HCH and CS out of plane angle
    # Plotting the HCH angle
    plt.figure()
    plt.bar(np.arange(len(method_HCH[1:])), simulation_HCH_plot, width=bar_width, color='C2', yerr=error_HCH[1:])
    
    # Plotting the C-H1 bond length
    plt.bar(np.arange(len(method_oop_cs[1:])) + bar_width, simulation_oop_cs_plot, width=bar_width, color='C3', yerr=error_oop_cs[1:])
    
    # Correcting the names of both sets of bars
    tick_positions = np.arange(len(method_CS[1:])) + bar_width*1
    plt.xticks(tick_positions, method_CS[1:])
    
    plt.ylabel('Difference with CC3 [Degrees]')
    plt.xlabel('Method')
    plt.title('Thioformaldehyde ground-state')
    plt.legend([r'$\angle H_1$C$H_2$', '$C-S$ out of plane angle'], loc='upper right')
    plt.tight_layout()
        

def thio_excitedstate():
    '''
    INPUTS
        none
    OUPUTS
        none

    This function plots all the barplots for excited state of thioformaldehyde
    '''
    simulation_CS = [1.709,1.69790, 1.697]
    method_CS = ['CC3', '$G_0W_0$, B3LYP','Expt']
    error_CS = [0,0,0.015] 
    simulation_CS_plot = simulation_CS[1:] - np.ones(len(simulation_CS[1:]))*simulation_CS[0]
    
    simulation_CH = [1.078,1.08012,1.085]
    method_CH = ['CC3','G_0W_0$, B3LYP','Expt']
    error_CH = [0,0,0.008]
    simulation_CH_plot = simulation_CH[1:] - np.ones(len(simulation_CH[1:]))*simulation_CH[0]
     
    simulation_HCH = [120.2,120.00710298234027,119.2] 
    method_HCH = ['CC3', 'G_0W_0$, B3LYP', 'Expt']
    error_HCH = [0,0,2.4]
    simulation_HCH_plot = simulation_HCH[1:] - np.ones(len(simulation_HCH[1:]))*simulation_HCH[0]
    
    simulation_oop_cs = [0.0, 0.0, 4.45]
    method_oop_cs = ['CC3', '$G_0W_0$, B3LYP', 'Expt']
    error_oop_cs = [0,0,4.45]
    visualization_aid = [0.1,0]
    simulation_oop_cs_plot = simulation_oop_cs[1:] - np.ones(len(simulation_oop_cs[1:]))*simulation_oop_cs[0] + visualization_aid
    # Note that we add the visualization aid such that the bar gets drawn in the first place


    # Plotting the data  
    # Set the width of each bar
    bar_width = 0.2
    
    # Create the combined bar graph
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    
    # Plotting the C-S bond length
    ax1.bar(np.arange(len(method_CS[1:])), simulation_CS_plot, width=bar_width, color='C0', yerr=error_CS[1:])
    ax1.set_ylabel('Difference with CC3 [$\AA$]')
    ax1.set_xlabel('Method')
    ax1.set_title('Thioformaldehyde first excited-state')
    
    # Plotting the C-H1 bond length
    ax1.bar(np.arange(len(method_CH[1:])) + bar_width, simulation_CH_plot, width=bar_width, color='C1', yerr=error_CH[1:])
    
    # Plotting the H-C-H angle
    ax2.bar(np.arange(len(method_HCH[1:])) + 2*bar_width, simulation_HCH_plot, width=bar_width, color='C2', yerr=error_HCH[1:])
    ax2.set_ylabel('Difference with CC3 [Degrees]')
    
    # Plotting the C-S out of plane angle
    ax2.bar(np.arange(len(method_oop_cs[1:])) + 3*bar_width, simulation_oop_cs_plot, width=bar_width, color='C3', yerr=error_oop_cs[1:])
    
    # Set the x-axis tick positions and labels
    tick_positions = np.arange(len(method_CS[1:])) + bar_width*1.5
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(method_CS[1:])
    
    # Create a legend
    ax1.legend(['$C-S$ bond length', '$C-H_1$ bond length'], loc='upper left')
    ax2.legend([r'$\angle H_1$C$H_2$', '$C-S$ out of plane angle'], loc='upper right')
    
    # Find the maximum absolute value among the data points
    max_abs_value = max(max(abs(simulation_CS_plot)), max(abs(simulation_CH_plot)), max(abs(simulation_HCH_plot)), max(abs(simulation_oop_cs_plot)))
    
    # Set the y-limits of both axes
    y_limit_1 = 0.014 * max_abs_value  # Adjust the factor (1.2) as needed to control the y-axis range
    y_limit_2 = 2.1 * max_abs_value
    ax1.set_ylim(-y_limit_1*0.5, y_limit_1)
    ax2.set_ylim(-y_limit_2*0.5, y_limit_2)
    
    # Display the plot
    plt.tight_layout()
    plt.show()
    
    
    # Two seperate plots
    bar_width = 0.4
    # Plot for C-S and C-H #####
    
    plt.figure()
    plt.bar(np.arange(len(method_CS[1:])), simulation_CS_plot, width=bar_width, color='C0', yerr=error_CS[1:])
    
    # Plotting the C-H1 bond length
    plt.bar(np.arange(len(method_CH[1:])) + bar_width, simulation_CH_plot, width=bar_width, color='C1', yerr=error_CH[1:])
    
    # Correcting the names of both sets of bars
    tick_positions = np.arange(len(method_CS[1:])) + bar_width*1
    plt.xticks(tick_positions, method_CS[1:])
    
    plt.ylabel('Difference with CC3 [$\AA$]')
    plt.xlabel('Method')
    plt.title('Thioformaldehyde first excited-state')
    plt.legend(['$C-S$ bond length', '$C-H_1$ bond length'], loc='upper right')
    plt.tight_layout()
    
    
    
    # Plot for HCH and CS out of plane angle
    # Plotting the HCH angle
    plt.figure()
    plt.bar(np.arange(len(method_HCH[1:])), simulation_HCH_plot, width=bar_width, color='C2', yerr=error_HCH[1:])
    
    # Plotting the C-H1 bond length
    plt.bar(np.arange(len(method_oop_cs[1:])) + bar_width, simulation_oop_cs_plot, width=bar_width, color='C3', yerr=error_oop_cs[1:])
    
    # Correcting the names of both sets of bars
    tick_positions = np.arange(len(method_CS[1:])) + bar_width*1
    plt.xticks(tick_positions, method_CS[1:])
    
    plt.ylabel('Difference with CC3 [Degrees]')
    plt.xlabel('Method')
    plt.title('Thioformaldehyde first excited-state')
    plt.legend([r'$\angle H_1$C$H_2$', '$C-S$ out of plane angle'], loc='upper right')
    plt.tight_layout()
    
    
def CO_groundstate():
    '''
    INPUTS
        none
    OUPUTS
        none

    This function plots all the barplots for ground state of carbon monoxide
    '''
    
    simulation_bondlenght = [1.1303225, 1.135338346,1.123308271,1.110275689,1.128]
    method = ['CC3','BLYP', 'B3LYP', 'BHLYP','Expt']
    differences = simulation_bondlenght - np.ones(len(simulation_bondlenght))*simulation_bondlenght[0]

    
    method_plot = method[1:]
    difference_plot = differences[1:]
    
    plt.figure()
    plt.title('Bond length ground-state carbon monoxide')
    plt.bar(method_plot, difference_plot, color=['C0','C0','C0','C1'])
    plt.xlabel('Computational method')
    plt.ylabel('Difference with CC3 [$\AA$]')
    plt.tight_layout()
    plt.show()
    

def CO_excitedstate():
    '''
    INPUTS
        none
    OUPUTS
        none

    This function plots all the barplots for excited state of carbon monoxide
    '''
    
    simulation_bondlength = [1.224, 1.271679198, 1.238596491, 1.212531328, 1.21754386, 1.235]
    method = ['EOM-CCSD','$G_0W_0$, BLYP', '$G_0W_0$, B3LYP', '$G_0W_0$, BHLYP','DFT, BHLYP','Expt']
    differences = simulation_bondlength - np.ones(len(simulation_bondlength))*simulation_bondlength[0]
    
    method_plot = method[1:]
    difference_plot = differences[1:]
    
    plt.figure()
    plt.title('Bond length first excited-state carbon monoxide')
    plt.bar(method_plot, difference_plot, color=['C0','C0','C0','C0','C1'])
    plt.xlabel('Computational method')
    plt.ylabel('Difference with EOM-CCSD [$\AA$]')
    plt.tight_layout()
    
    
    plt.show()
    
def CO_self_consistant():
    '''
    INPUTS
        none
    OUPUTS
        none

    This function plots all the barplots for excited state of carbon monoxide of GW self consistant
    '''   
 
    # Setting the excited state values
    simulation_bondlength_excited = [1.224,1.267669173,1.234586466,1.215538847, 1.271679198, 1.238596491, 1.212531328, 1.235]
    method_excited = ['EOM-CCSD','$BLYP$','$B3LYP$','$BHLYP$','$BLYP $', '$B3LYP $', '$BHLYP $','Expt']
    differences_excited = simulation_bondlength_excited - np.ones(len(simulation_bondlength_excited))*simulation_bondlength_excited[0]
    
    method_plot = method_excited[1:]
    difference_plot = differences_excited[1:]

    
    plt.figure()
    plt.title('Bond length first excited-state carbon monoxide')
    bars = plt.bar(method_plot, difference_plot, color=['C2', 'C2', 'C2', 'C0', 'C0', 'C0', 'C1'])
    plt.xlabel('Computational method')
    plt.ylabel('Difference with EOM-CCSD [$\AA$]')
    plt.tight_layout()
    
    # Create a custom legend
    custom_legend = [
        plt.Rectangle((0, 0), 1, 1, color='C2', label='$G_9W_9$'),
        plt.Rectangle((0, 0), 1, 1, color='C0', label='$G_0W_0$'),
        plt.Rectangle((0, 0), 1, 1, color='C1', label='Expt')
    ]
    
    # Add the custom legend to the plot
    plt.legend(handles=custom_legend)
    
    # Adjust the layout to avoid overlapping labels
    plt.tight_layout()
    
    # Display the plot
    plt.show()

    
    return

def CO_forces():
    '''
    INPUTS
        none
    OUPUTS
        none

    This function plots all the barplots of excited state forces the first excited state carbon monoxide
    
    '''   
    
    forces = [-9.875095388,-10.11359903,-7.54003187705798,-10.6259327894717]
    methods = ['$G_0W_0$, BLYP','$G_0W_0$, B3LYP','$G_0W_0$, BHLYP','DFT BHLYP']

    
    plt.figure()
    plt.title('Excited-state forces of first excited-state carbon monoxide')
    bars = plt.bar(methods, forces)
    plt.xlabel('Computational method')
    plt.ylabel(r'Force $\left[\frac{eV}{\AA} \right]$')
    plt.tight_layout()

def thio_forces():
    '''
    INPUTS
        none
    OUPUTS
        none

    This function plots all the barplots for excited state forces of thioformaldehyde
    
    '''   
    
    print('Units of this are still a bit unclear to me')
    
    
    forces = [-2.88429096, 0.11073035, 0.11073435, 0.01372497, 0.01372597, 0.]
    what_derivative = ['C - S Bond', 'C - $H_1$ Bond', 'C - $H_2$ Bond', r'$\angle$SC$H_1$', r'$\angle$SC$H_2$', r'$\theta$($H_2CS$,$CSH_1$)']
    
    # Set the width of each bar
    bar_width = 0.35
    gap = 0.1  # Gap between the bars
    
    # Create the figure and the first axis
    fig, ax1 = plt.subplots()
    ax1.set_ylabel(r'Forces [$\mathrm{eV/\AA}$]')
    
    # Calculate the positions for ax1 bars
    num_bars1 = len(forces[:3])
    bar_positions1 = np.arange(num_bars1) * (bar_width + gap)
    
    # Plotting the first three data points on the first axis
    bars1 = ax1.bar(bar_positions1, forces[:3], width=bar_width, color='C0')
    
    # Create the second axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('Torque [$\mathrm{eV/deg}$]')
    
    # Calculate the positions for ax2 bars
    num_bars2 = len(forces[3:])
    bar_positions2 = np.arange(num_bars2) * (bar_width + gap) + (num_bars1 * (bar_width + gap))
    
    # Plotting the remaining three data points on the second axis
    bars2 = ax2.bar(bar_positions2, forces[3:], width=bar_width, color='C0')
    
    # Set the zero points of both axes at the same height
    min_y = 1.2*min(min(forces[:3]), min(forces[3:]), 0)
    max_y = 1.2*max(max(forces[:3]), max(forces[3:]), 0) + 2
    
    ax1.set_ylim(bottom=min_y, top=max_y)
    ax2.set_ylim(bottom=0.01*min_y, top=0.01*max_y)
    
    # Set the x-axis tick positions and labels
    tick_positions = np.concatenate([bar_positions1, bar_positions2])
    tick_labels = what_derivative[:3] + what_derivative[3:]
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(tick_labels, rotation=45, ha='right')
    
    ax1.set_title('Thioformaldehyde excited-state forces')
    
    # Adjust the spacing between subplots
    fig.tight_layout()
    
    # Display the plot
    plt.show()


    # Two seperate plots
    # Create the figure and the first axis
    fig1, ax1 = plt.subplots()
    ax1.set_ylabel(r'Force $\left[\frac{eV}{\AA}\right]$')
    
    # Calculate the positions for ax1 bars
    num_bars1 = len(forces[:3])
    bar_positions1 = np.arange(num_bars1) * (bar_width + gap)
    
    # Plotting the first three data points on the first axis
    bars1 = ax1.bar(bar_positions1, forces[:3], width=bar_width, color='C0')
    
    # Set the zero points of the axis at the same height
    min_y = 1.2 * min(min(forces[:3]), 0) 
    max_y = 1.2 * max(max(forces[:3]), 0) + 2
    
    ax1.set_ylim(bottom=min_y, top=max_y)
    
    # Set the x-axis tick positions and labels
    ax1.set_xticks(bar_positions1)
    ax1.set_xticklabels(what_derivative[:3], rotation=45, ha='right')
    
    ax1.set_title('Thioformaldehyde excited-state forces')
    
    # Adjust the spacing between subplots
    fig1.tight_layout()
    
    # Display the plot
    plt.show()
    

    # The second plot corresponding to the 'torques'
    fig2, ax2 = plt.subplots()
    ax2.set_ylabel(r'Torque $\left[\frac{eV}{Degrees}\right]$') 
    # Calculate the positions for ax2 bars
    num_bars2 = len(forces[3:])
    bar_positions2 = np.arange(num_bars2) * (bar_width + gap)
    
    # Plotting the remaining three data points on the second axis
    bars2 = ax2.bar(bar_positions2, forces[3:], width=bar_width, color='C0')
    
    # Set the zero points of the axis at the same height
    min_y = 1.2 * min(min(forces[3:]), 0)
    max_y = 1.2 * max(max(forces[3:]), 0)
    
    ax2.set_ylim(bottom=min_y, top=max_y)
    
    # Set the x-axis tick positions and labels
    ax2.set_xticks(bar_positions2)
    ax2.set_xticklabels(what_derivative[3:], rotation=45, ha='right')
    
    ax2.set_title('Thioformaldehyde excited-state forces')
    
    # Adjust the spacing between subplots
    fig2.tight_layout()
    
    # Display the plot
    plt.show()



plt.close('all')

plt.rc('axes', labelsize=14)
plt.rc('xtick', labelsize = 12)

#CO_groundstate()
#CO_excitedstate()
#CO_self_consistant()
#thio_groundstate()
thio_excitedstate()
#CO_forces()
#thio_forces()



























