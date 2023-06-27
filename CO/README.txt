The simulations for the CO molecule where the first that were performed.
For these simulations two methods of getting to the minima where used.




###### Steps ######
The first method was by doing a lot of simulations, the corresponding scripts are located in the 'Steps' folder. The steps folder has two important script needed to run the simulations. 

'Make_molgw11.py'
'Get_data6.py'


'Make_molgw11.py' is the script which generates all the scripts required to perform the simulation. This script generates all the submission files required by molgw and when running the program several parameters can be changed in the console such as the amount of steps, the min bond length and the maximum bondlength can be set. This script can with some adaptation easily be used for other diatomic atoms. However for more complex molecules this script is of little use. Note when running 'Make_molgw11.py' 'Get_data6.py' must be present in the same folder! The 'Get_data6.py' file gets copied to all relevant folder and locations.

Once the Make_molgw11.py script has been run the simulation can be started. This is not difficult since the bash scripts are also automatically generated. Simply sbatch either the DFT_BSE.sub or the GW_BSE.sub.

After the simulation has been completed the 'Get_data6.py' is located in all the subfolders for each of the seperate simulation setting combitnations. Running the 'Getdata6.py' script in these folder creates a file 'Energies_*functional*_*basis*_*DFT or GW*.npy for each of the different folders which have a completed simulation in them. These filesh contains all the energy of the ground state and some of the excited states. By changing one value in the 'Get_data6.py' script the exact amount of excitation energies that are included can easily be changed to the desired amount.


After the simulation the created energies file is structured with the first column all groundstate energies, the second the first excitation energies etc.
And the rows denote which of the steps the energy belongs to. The python scripts that were used to analyze this data can be found in the folder 'visualization', however writing new python scripts to visualize this data should be straight forward.


Finally a third script 'Make_molgw_evGW.py' can be used instead of 'Make_molgw11.py' if several GW iterations need to be performed. So for GnWn simulations.
The Steps_GW parameter can be changed for the amount of iterations that must be performed.

###### Gradient descent for CO ######

!!! Important !!!
This code was only written to verify the performace of the CO molecule. The gradient descent scripts of Thioformaldehyde where still improved on later and those scripts are 'nicer'.

The gradient descent script has also been run on the carbon monoxide molecule to verify its performance. The corresponding scripts are located in the 'Gradient descent' folder.

The gradient descent algorithm requires more scripts than the Steps method. The scripts required are:

'gw.in'
'bse.in'
'CO.zmat'
'bash.sub'
'gcutil.py'
'Make_derivatives_1.py'
'Gradient_descent.py'
'Get_data.py'
'Get_co_bondlength.py'


First of all there are 'gw.in' and 'bse.in'. Since the .xyz script of the gradient descent script is generated seperately the submission scripts for molgw do not change during the simulation.
Because they remain the same they are now both not generated using a python script anymore. One can enter these two files to switch the functional, basis, frozencore, number of GnWn steps or any other settings within the molgw program. 
The next file required is 'CO.zmat' which is the initial geometry of the as of yet not optimized CO molecule in internal coordinates. An example file is also added in the folder.
Then we have 'bash.sub' which must be run to start the molgw simulation.
'gcutil.py' contains a variation of functions used in the other scripts. A large part of the function is related with converting the .xyz and .zmat files into one another. This function also contains work written by Robert A Shaw.
'Make_derivatives_1.py' looks at the current .zmat file and creates for each degree of freedom two new zmat files. Where each of these new zmat files has one of the parameters slightly larger of slightly lower than the current configuration. These files are then automatically evalueated in the next steps of the 'bash.sub' script.
'Gradient_descent.py' Takes the all these evaluated data points and creates a new .zmat file of the new position. The learning rate can be set in this document.
'Get_data.py' Is a script that find the total energy at each step in the simulation and saves this in a 'path_energies.npy' file that can be evaluated later. This function is nice to see if the gradient descent algorithm is working as intended and to verify it is not overshooting.
'Get_co_bondlength.py' has been written to retriee the bondlength of the co molecule at each iteration of the simulation and save it to a 'CO_bondlengths.npy' file. This data was used to compare the gradient descent algorithm with the Steps method discussed above.



