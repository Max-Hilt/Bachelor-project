# -*- coding: utf-8 -*-
"""
Created on Mon May 22 11:44:15 2023

@author: Max Hilt
"""

# Importing other script
import Visualizing_functions as vs

# Importing libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,mark_inset)

# First we close all the old files
plt.close('all')




# Distance parameterization
min_ang = 1.0
max_ang = 1.4
steps_ang = 400
angstrom = np.linspace(min_ang,max_ang,steps_ang)
dr = angstrom[1] - angstrom[0]


# Getting plots and data
data, string = vs.user_input_plot()
correct_data = vs.correct_data(data)


# This piece of code is only needed if correct data failed to nicely distinct
# The different excited states.

# maxnumber = 282
# minnumber = maxnumber - 26
# correct_data = correct_data[minnumber:maxnumber]
# angstrom = angstrom[minnumber:maxnumber]


min_ground, min_excited,force = vs.find_minima(correct_data, angstrom)
string = 'CO simulation using B3LYP and cc-pV5Z'
vs.one_excitation(correct_data, angstrom, string, min_ground, min_excited)


print(f'Minimum ground: {min_ground}, min_excited: {min_excited}, force: {force}')



print(f'Used angstrom grid {min_ang} - {max_ang} with {steps_ang} steps')






