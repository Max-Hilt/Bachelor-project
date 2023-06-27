# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:36:03 2023

@author: Max Hilt
"""

'''
This entire script is only written to make the plot comparing the different methods

Since automating the analyzing of the results was difficult and prone to error,
Hence the results are copied here.
'''

# Importing libraries
import numpy as np
import matplotlib.pyplot as plt

# Filling in the data points

# The order is DFT_BYLP, DFT_B3LYP, DFT_BHLYP, GW_BYLP, GW_B3LYP, GW_BHLYP
# Where each of these has 4 points corresponding to pVDZ, pVTZ, pVQZ, pV5Z. In this order
Ground_state_scf = np.array([1.147368421,1.137343358,1.135338346,1.135338346,1.13433584,1.126315789,1.123308271,1.123308271,1.120300752,1.113283208,1.111278195,1.110275689,1.147368421,1.137343358,1.135338346,1.135338346,1.13433584,1.126315789,1.123308271,1.123308271,1.120300752,1.113283208,1.111278195,1.110275689])
Excited_state_scf = np.array([None,None,None,None,None,None,None,None,None,1.224561404,1.21954887218045,1.21754385964912,1.29573934837092,1.28270676691729,1.274686717,1.271679198,1.26867167919799,1.24962406015037,1.24160401002506,1.23859649122807,1.24060150375939,1.22556390977443,1.21553884711779,1.2125313283208])

stacked = np.stack((Ground_state_scf,Excited_state_scf))

# We define the corresponding plot titles and lables
Names = np.array(['Comparison various methods for calculating groundstate bond length    ', 'Excited state','Ground state with mbpt energies', 'First excited state with mbpt energies'])
Basis = np.array(['cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-pV5Z'])
Functionals = np.array(['BLYP', 'B3LYP', 'BHLYP'])
Method = np.array(['DFT', 'GW'])

# Now we plot this data
plt.close('all')

# We iterate over 2 Excited and 2 ground state datasets
for i in range(2): 
    # This loops over DFT and GW
    for k in range(2):
        # Setting the plot title
        plt.figure()
        plt.title(Names[i])# + f' Using {Method[k]}'
        
        # We plot the theoretical value inside the plot
        if i == 0 or i == 2:
            # We are plotting an ground state so the theory value is 1.128
            #plt.axhline(y = 1.128, color = 'r', linestyle = '-', label='Experiment')
            print('no exp')
        else:
            # We are plotting an excited state so the theory value is 1.235
            #plt.axhline(y = 1.235, color = 'r', linestyle = '-', label='Experiment')
            print('no exp')
            
        # This loops over BLYP B3LYP and BHLYP
        for j in range(3):
            plt.plot(Basis, stacked[i,(j*4 + k*12):((j+1)*4 + k*12)],marker = '.',linestyle = '--', label= f'{Method[k]} {Functionals[j]}')
            print(f'{(j*4 + k*12)}:{((j+1)*4 + k*12)}')
        
        # Adding the legend and axis
        plt.legend()
        plt.ylabel('Distance [$\AA$]')
        plt.xlabel('Basis set')








