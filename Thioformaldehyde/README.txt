For the gradient descent algorithm with Thioformaldehyde the following 8 files are used:

'Thio.zmat'
'gw.in'
'bse.in'
'bash.sub'
'utils.py'
'Make_derivates_1.py'
'Gradient_descent.py'
'Get_data.py'


These files have the following functionallity:

'Thio.zmat' has the internal coordinates of the current iteration. When starting the simulation this file has the starting geometry. Old Thio.zmat files will be stored in a 'gradient_path' folder that is created.

'gw.in' is the submission file for the molgw simulations. In this file the simulation settings can be changed. For example the used functional, basis and amount of GnWn iterations can be set.

'bse.in' is the second molgw simulation script. This script is for the BSE step. Again this file contains the simulation settings that are used.

'bash.sub' This is the sbatch file which automatically runs the python scripts and submits the molgw scripts to perform the gradient descent. If all the files are present in the folder submitting this file should be everything that needs to be done.

During the simulation two folders will be generated with subfolders for all iterations.
The first folder is 'gradient_path' the gradient path folder contains the configuration of the molecule on each step of the simulation together with all the molgw outputs for that configuration of the molecule. Currently this data is only used to find the total energy as a function of iteration but more interesting properties could be gathered from this data.
The second folder generated is 'Simulations' folder contains for each iteration the simulations that needed to be performed to determine the gradient. The gradient is needed to perform the gradient descent script.


After the threshold (set in the gradient descent script) has been reached a file will automatically be generated called simulation.stop Once bash.sub sees this file the simulation is stopped and Get_data.py is automatically run.

'utils.py' Contains the various functions that are used by several of the other scripts. It mainly contains all funcitons related to converting internal coordiantes to atomic coordinates and back.

'Make_verivatives_1.py' Creates the geometries that need to be simulated with MOLGW to find the gradient. How much the various parameters are changed can be set, one value for distance, angle and dihedral angle. 

'Gradient_descent.py' This script performes the gradient descent script. In this script several parameters can and must be set. The first is the dimensions which is the degrees of the molecule under consideration. For thioformaldehyde the parameter is set to 6. Furtheremore the threshold must be set in this script which determines when the script is 'converged'. Finally the 'delta' distance and angle parameters must be set which are how far the parameters were varied to calculate the gradient.

'Get_data.py' This file is run last to collect the energies of the molecule at each step in the simulation. This can be helpful to see if the gradient descent algorithm is working properly. Note that this function can be easily adapted to collect other relevant information. The file created is 'path_energies.py'