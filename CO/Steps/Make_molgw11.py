# Author: Max Hilt

# This script is specially written for the CO molecule.
# Importing libraries
import numpy as np
import os
import yaml
import shutil

def make_gw_bse_bash():
    # Then we also make the submit file for the DFT, BSE job

    f = open(f'GW_BSE.sub', 'w+')

    f.write('#!/bin/bash\n')
    f.write('#SBATCH --time=2-23:00:00\n')
    f.write('#SBATCH -p ccp22\n')
    f.write('#SBATCH -N 1 --exclusive --ntasks-per-node 2\n')
    f.write('#SBATCH -J CO\n')
    f.write('#SBATCH --output=log%j\n')
    f.write('#SBATCH --ntasks-per-core=1\n')
    f.write('#SBATCH --error=error%j\n\n')
    f.write('module purge\n')
    f.write('module load python3\n')
    f.write('module load molGW/latest\n\n')
    f.write('ulimit -s unlimited\n\n')
    f.write('for folder in */;')
    f.write('do \n')
    f.write('   # Going over all the bond lengths of interest\n')
    f.write(f'  for i in {{0..{steps-1}..1}} \n')
    f.write('   do\n')
    f.write('        echo "Running configuration $i"\n\n')
    f.write('       # First we do the GW calculation and save that\n')
    f.write('       mpirun molgw ${folder}runfiles/GW/gw$i.in > gw.out\n')
    f.write('       mkdir -p ${folder}simulations/sim$i/GW\n')
    f.write('       mkdir -p ${folder}simulations/sim$i/BSE\n')
    f.write('       mv molgw.yaml ${folder}simulations/sim$i/GW/molgw.yaml\n')
    f.write('       mv gw.out ${folder}simulations/sim$i/GW/gw.out\n')
    f.write('       # Performing the BSE calculations\n')
    f.write('       mpirun molgw ${folder}runfiles/BSE/bse$i.in > bse.out\n')
    f.write('       mv molgw.yaml ${folder}simulations/sim$i/BSE/molgw.yaml\n')
    f.write('       mv bse.out ${folder}simulations/sim$i/BSE/bse.out\n')
    f.write('       mv RESTART ${folder}simulations/sim$i/BSE/RESTART\n')
    f.write('       mv ENERGY_QP ${folder}simulations/sim$i/BSE/ENERGY_QP\n')
    f.write('       mv dynamical_dipole_polarizability.dat ${folder}simulations/sim$i/BSE/dynamical_dipole_polarizability.dat \n')
    f.write('       mv photoabsorption_cross_section.dat ${folder}simulations/sim$i/BSE/photoabsorption_cross_section.dat\n')
    f.write('   done\n')
    f.write('   # Finally we collect the data to a good readable file\n')
    f.write('   mv ${folder}Post_processing.yaml Post_processing.yaml\n')
    f.write('   python3 Get_data6.py\n')
    f.write('   mv Post_processing.yaml ${folder}Post_processing.yaml\n')
    f.write('done')
    f.close()
    
    return


# Defining functions that write the files
def make_gw_bse(make_folder, path, scf_index, basis_index):

    # make_folder bool
    # path string must include the / at the end.

    if make_folder:
        os.mkdir(f'./{path}')
        path = path + '/'
    
    # Create the path to store the run files
    os.mkdir(f'./{path}runfiles')
    os.mkdir(f'./{path}runfiles/GW')
    os.mkdir(f'./{path}runfiles/BSE')

    # Choosing the correct scf and basis
    scf = scf_list[scf_index]
    basis = basis_list[basis_index]
    auxil_basis = auxil_basis_list[basis_index]

    # Creating all the files required to start simulations
    for i in range(steps):
        f = open(f'{path}runfiles/GW/gw{i}.in', 'w+')

        # Writing the input to the GW file
        f.write('&molgw\n')
        f.write('  comment=\'CO GW\'    ! an optional plain text here\n\n')
        f.write(f'  scf={scf}\n\n')
        f.write(f'  basis={basis}\n')
        f.write(f'  auxil_basis={auxil_basis}\n\n')
        f.write(f'  postscf={postscf}\n\n')
        f.write(f'  selfenergy_state_range={selfenergy_state_range}\n')
        f.write(f'  frozencore={frozencore}              ! accurate approximation: O1s will not be included\n')
        f.write('                                ! in the RPA/GW calculation \n\n')
        f.write(f'  natom={natom}\n')
        f.write(f'/\n')
        f.write('C      0.000000  0.000000  0.000000\n')
        f.write(f'O      0.000000  0.000000  {distances[i] :.6f}\n')

        f.close()

        f = open(f'{path}runfiles/BSE/bse{i}.in', 'w+')
        
        # Writing the input to the BSE file
        f.write('&molgw\n')
        f.write('  comment=\'CO BSE\'\n\n')
        f.write(f'  scf={scf}\n')
        f.write(f'  read_restart={read_restart}    ! read the RESTART file so to skip the DFT part\n')
        f.write(f'                        ! just to save time\n\n')
        f.write(f'  basis={basis}\n')
        f.write(f'  auxil_basis={auxil_basis}\n')
        f.write(f'  postscf=\'bse\'\n')
        f.write(f'  frozencore={frozencore}\n\n')
        f.write(f'  eta={eta}\n\n')
        f.write(f'  natom={natom}\n')
        f.write(f'/\n')
        f.write('C      0.000000  0.000000  0.000000\n')
        f.write(f'O      0.000000  0.000000  {distances[i] :.6f}\n')

        f.close()

    # Finally we create a file that gives all the simulation settings that were used during the simulation
    f = open(f'{path}Simulation_settings.txt', 'w+')
    f.write('##### This simulation had the following settings #####\n\n')
    f.write(f'Was GW used:  Yes\n\n')
    f.write(f'Minlength: {minlength}\n')
    f.write(f'Maxlength: {maxlength}\n')
    f.write(f'Steps: {steps}\n\n')
    f.write(f'scf: {scf}\n')
    f.write(f'Basis: {basis}\n')
    f.write(f'Auxil_basis: {auxil_basis}\n')
    f.write(f'postscf: {postscf}\n')
    f.write(f'selfenergy_state_range: {selfenergy_state_range}\n')
    f.write(f'frozencore: {frozencore}\n')
    f.write(f'natom: {natom}\n')

    f.close()

    # Also we write a small YAML file for assisting the data anaylsis file

    data = dict(
        Steps = steps,
        Scf = scf_text_list[scf_index],
        Basis = basis_text_list[basis_index],
        Path = path,
        Method = 'GW'
    )

    with open(f'{path}Post_processing.yaml', 'w') as outfile:
        yaml.dump(data, outfile)


    if make_folder:
        # And we copy the get data file to each subfolder
        src_path = 'Get_data6.py'
        dst_path = f'{path}Get_data6.py'
        shutil.copy(src_path, dst_path)
    return

def make_dft_bse_bash():
    # Then we also make the submit file for the DFT, BSE job

    f = open(f'DFT_BSE.sub', 'w+')

    f.write('#!/bin/bash\n')
    f.write('#SBATCH --time=2-23:00:00\n')
    f.write('#SBATCH -p ccp22\n')
    f.write('#SBATCH -N 1 --exclusive --ntasks-per-node 2\n')
    f.write('#SBATCH -J CO\n')
    f.write('#SBATCH --output=log%j\n')
    f.write('#SBATCH --ntasks-per-core=1\n')
    f.write('#SBATCH --error=error%j\n\n')
    f.write('module purge\n')
    f.write('module load python3\n')
    f.write('module load molGW/latest\n\n')
    f.write('ulimit -s unlimited\n\n')
    f.write('for folder in */;')
    f.write('do \n')
    f.write('   # Going over all the bond lengths of interest\n')
    f.write(f'  for i in {{0..{steps-1}..1}} \n')
    f.write('   do\n')
    f.write('        echo "Running configuration $i"\n\n')
    f.write('       # Performing DFT + BSE calculations\n')
    f.write('       mpirun molgw ${folder}runfiles/DFT_BSE/dft_bse$i.in > dft_bse.out\n')
    f.write('       mkdir -p ${folder}simulations/sim$i/DFT_BSE\n')
    f.write('       mv molgw.yaml ${folder}simulations/sim$i/DFT_BSE/molgw.yaml\n')
    f.write('       mv dft_bse.out ${folder}simulations/sim$i/DFT_BSE/gw.out\n')
    f.write('       mv RESTART ${folder}simulations/sim$i/DFT_BSE/RESTART\n')
    f.write('       mv dynamical_dipole_polarizability.dat ${folder}simulations/sim$i/DFT_BSE/dynamical_dipole_polarizability.dat \n')
    f.write('       mv photoabsorption_cross_section.dat  ${folder}simulations/sim$i/DFT_BSE/photoabsorption_cross_section.dat\n')
    f.write('   done\n')
    f.write('   # Finally we collect the data to a good readable file\n')
    f.write('   mv ${folder}Post_processing.yaml Post_processing.yaml\n')
    f.write('   python3 Get_data6.py\n')
    f.write('   mv Post_processing.yaml ${folder}Post_processing.yaml\n')
    f.write('done')
    f.close()

    return

def make_dft_bse(make_folder, path, scf_index, basis_index):

        # make_folder bool
    # path string must include the / at the end.

    if make_folder:
        os.mkdir(f'./{path}')
        path = path + '/'

    
    # Checking if there is another subpath to be added

    # Create the path to store the run files
    os.mkdir(f'./{path}runfiles')
    os.mkdir(f'./{path}runfiles/DFT_BSE')

    # Choosing the correct scf and basis
    scf = scf_list[scf_index]
    basis = basis_list[basis_index]
    auxil_basis = auxil_basis_list[basis_index]

    for i in range(steps):
    # Writing the input to the BSE file
        f = open(f'{path}runfiles/DFT_BSE/dft_bse{i}.in', 'w+')

        f.write('&molgw\n')
        f.write('  comment=\'CO BSE\'\n\n')
        f.write(f'  scf={scf}\n')
        f.write(f'  read_restart={read_restart_DFT}    ! read the RESTART file so to skip the DFT part\n')
        f.write(f'                        ! just to save time\n\n')
        f.write(f'  basis={basis}\n')
        f.write(f'  auxil_basis={auxil_basis}\n')
        f.write(f'  postscf=\'bse\'\n')
        f.write(f'  frozencore={frozencore}\n\n')
        f.write(f'  eta={eta}\n\n')
        f.write(f'  natom={natom}\n')
        f.write(f'/\n')
        f.write('C      0.000000  0.000000  0.000000\n')
        f.write(f'O      0.000000  0.000000  {distances[i] :.6f}\n')

        f.close()

    # Finally we create a file that gives all the simulation settings that were used during the simulation
    f = open(f'{path}Simulation_settings.txt', 'w+')
    f.write('##### This simulation had the following settings #####\n\n')
    f.write(f'Was GW used:  No\n\n')
    f.write(f'Minlength: {minlength}\n')
    f.write(f'Maxlength: {maxlength}\n')
    f.write(f'Steps: {steps}\n\n')
    f.write(f'scf: {scf}\n')
    f.write(f'Basis: {basis}\n')
    f.write(f'Auxil_basis: {auxil_basis}\n')
    f.write(f'postscf: \'bse\'\n')
    f.write(f'selfenergy_state_range: {selfenergy_state_range}\n')
    f.write(f'frozencore: {frozencore}\n')
    f.write(f'natom: {natom}\n')

    f.close()

    # Also we write a small YAML file for assisting the data anaylsis file

    data = dict(
        Steps = steps,
        Scf = scf_text_list[scf_index],
        Basis = basis_text_list[basis_index],
        Path = path,
        Method = 'DFT'
    )

    with open(f'{path}Post_processing.yaml', 'w') as outfile:
        yaml.dump(data, outfile)

    
    if make_folder:
        # And we copy the get data file to each subfolder
        src_path = 'Get_data6.py'
        dst_path = f'{path}Get_data6.py'
        shutil.copy(src_path, dst_path)
    return


    


# Setting parameters
scf_list = ['\'BLYP\'','\'B3LYP\'','\'BHLYP\'']
basis_list = ['\'cc-pVDZ\'','\'cc-pVTZ\'','\'cc-pVQZ\'','\'cc-pV5Z\''] # In theory PV6Z should also be there but I ran into some problems.
auxil_basis_list=['\'cc-pVDZ-RI\'','\'cc-pVTZ-RI\'','\'cc-pVQZ-RI\'','\'cc-pV5Z-RI\'']
methods = ['GW', 'DFT']
postscf= '\'G0W0\''
selfenergy_state_range=20           # How many layers above and below HOMO are GW corrected
frozencore='\'no\''
read_restart = '\'yes\''
read_restart_DFT = '\'no\''
eta = 0.01
natom = 2

# Creating text versions of the files to get proper path names
scf_text_list = ['BLYP','B3LYP','BHLYP']
basis_text_list = ['cc-pVDZ','cc-pVTZ','cc-pVQZ','cc-pV5Z']


# Letting the user pick parameters once the script is run
minlength = float(input('What should the minimum bondlength be in Angstrom?  '))
maxlength = float(input('What should the maximum bondlength be in Angstrom?  '))
steps = int(input('How many steps do you want to use?  '))
gw_dft = input('Do you want to use GW? [y,n]  ')
all_variations = input('Do you want to do all simulations? [y,n]  ')


# Creating an array with all the distances that will be simulated
distances = np.linspace(minlength,maxlength,steps)

if all_variations == 'y':
    # We want to create all possible configurations for either GW or DFT calculations

    # To count how many files have been created
    file_count = 0

    if gw_dft =='y':
    # Here we want to use GW calculation
        for i, scf_i in enumerate(scf_text_list):
            for j, basis_i in enumerate(basis_text_list):
                # We are now creating a GW file
                path = 'GW' + '_' + scf_i + '_' + basis_i
                make_gw_bse(True, path, i, j)

                # Keeping track of how many files have been created
                file_count += 2*steps + 1
        # Finally we also write our submit file
        file_count += 1
        make_gw_bse_bash()
    else:
    # Here we want to use DFT calculation
        for i, scf_i in enumerate(scf_text_list):
            for j, basis_i in enumerate(basis_text_list):
                # We are now creating a DFT file
                path = 'DFT' + '_' + scf_i + '_' + basis_i
                make_dft_bse(True, path, i, j)

                # Keeping track of how many files have been created
                file_count += steps + 1

        # Finally we also write our submit file
        file_count += 1
        make_dft_bse_bash()

    # Printing how many files have been created
    print(f'{len(scf_text_list)*len(basis_text_list)*len(methods)} Folders and {file_count} files have been created')
   
elif gw_dft == 'y':
    # We are only creating one file, and we want to use GW

    # Now we must specify what scf and basis is desired
    selected_scf = int(input('What scf do you want to use?  '))
    print(f'scf {scf_list[selected_scf]} was selected')
    selected_basis = int(input('What basis do you want to use?  ')) - 2
    print(f'basis {basis_list[selected_basis]} was selected')

    # Finally we create the files
    make_gw_bse(True, 'GW', selected_scf, selected_basis)

    # Final message to user indicating how many files were created
    print(f'{2*steps + 2} files have been created.')

else:
    # We are only creating one file, and want to use DFT

    # Now we must specify what scf and basis is desired
    selected_scf = int(input('What scf do you want to use?  '))
    print(f'scf {scf_list[selected_scf]} was selected')
    selected_basis = int(input('What basis do you want to use?  ')) - 2
    print(f'basis {basis_list[selected_basis]} was selected')

    # Finally we create the files
    make_dft_bse(True, 'DFT', selected_scf, selected_basis)

    # Final message to user indicating how many files were created
    print(f'{steps + 2} files have been created.')
    make_dft_bse_bash()
    



