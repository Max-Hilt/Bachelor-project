In order to visualize the energy surface of thioformaldehyde these scripts were created. Using these scripts 2 of the 6 parameters will be varied while the rest remain fixed at a chosen value. This allows to create a 3D plot which 'visualizes' 2 of the dimensions of the energy surface. The following files were needed for these simulations:

'gw.in'
'bse.in'
'bash.sub'
'Thio.zmat'
'generate_scripts.py'
'Get_data_2D.py'
'utils.py'


'gw.in' is the submission file for the molgw simulations. In this file the simulation settings can be changed. For example the used functional, basis and amount of GnWn iterations can be set.

'bse.in' is the second molgw simulation script. This script is for the BSE step. Again this file contains the simulation settings that are used.

'bash.sub' is the bash file which automatically runs the simulation after the files have been generated with 'generate_scripts.py'. The bash file goes over all the different configurations and saves the simulations results in the folder 'Simulations/location{number}'.

'Thio.zmat' All coordinates that are not varied get the value from this configuration file. So it can be said we are 'looking in the neighborhood' of this geometry.

'generate_scripts.py' This is the script which lets you choose the variables that are gonna get varied, in between what values they are gonna get varied and in how many steps they are gonna be varied. All of these questions are asked in the console once the script is run. There will also be a file called 'Simulation_settings.txt' which saves the selected parameters.
! Known limitation: you can twice select the same variable to be changed but this has not been accounted for in the code and can yield problems. If one is interested in still doing so it should be a relatively simple addition.

'Get_data_2D.py' After the simulation has been completed this python script can be run to save the ground state + some excited state energies (the exact amount can be set by changing the parameter num_excitations) at all grid positions into a numpy file called 'Surface_energy.npy'. This file can later be used to visualize the data.

'utils.py' Contains the various functions that are used by several of the other scripts. It mainly contains all funcitons related to converting internal coordiantes to atomic coordinates and back.