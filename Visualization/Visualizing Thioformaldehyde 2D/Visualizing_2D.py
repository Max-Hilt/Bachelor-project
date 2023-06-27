# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 12:18:48 2023

@author: Max Hilt


This script is specifically written for the verification of convergence section to
'visualize' the region around the found minima.

It creates one plot with four corresponding subplots

"""

# Importing libraries
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm

plt.close('all')

# Loading the ground state data

#data1 = np.load('ground_0_1.npy')
#data2 = np.load('ground_1_2.npy')
#data3 = np.load('ground_0_3.npy')
#data4 = np.load('ground_4_5.npy')

# Loading the excited state data
data1 = np.load('excited_0_1.npy')
data2 = np.load('excited_1_2.npy')
data3 = np.load('excited_0_3.npy')
data4 = np.load('excited_4_5.npy')


# Adding the groundstate energy to the excitation energy
data1[:,:,1] = data1[:,:,1] + data1[:,:,0]
data2[:,:,1] = data2[:,:,1] + data2[:,:,0]
data3[:,:,1] = data3[:,:,1] + data3[:,:,0]
data4[:,:,1] = data4[:,:,1] + data4[:,:,0]

data1[:,:,2] = data1[:,:,2] + data1[:,:,0]
data2[:,:,2] = data2[:,:,2] + data2[:,:,0]
data3[:,:,2] = data3[:,:,2] + data3[:,:,0]
data4[:,:,2] = data4[:,:,2] + data4[:,:,0]

# Setting the x and y axis for all four subplots

x1 = np.linspace(1.2,2.1,50)
y1 = np.linspace(0.8,1.4,50)

x2 = np.linspace(0.8,1.4,50)
y2 = np.linspace(0.8,1.4,50)

x3 = np.linspace(1.2,2.1,50)
y3 = np.linspace(40,200,50)

x4 = np.linspace(40,200,50)
y4 = np.linspace(120,240,50)

# Now we use meshgrid to get the grids which we plot
X1, Y1 = np.meshgrid(x1,y1)
X2, Y2 = np.meshgrid(x2,y2)
X3, Y3 = np.meshgrid(x3,y3)
X4, Y4 = np.meshgrid(x4,y4)

# Create a figure and subplots
fig = plt.figure(figsize=(16, 6))
fig2 = plt.figure(figsize=(16, 6))

# First subplot
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
ax1.plot_surface(X1, Y1, data1[:,:,1], cmap='plasma')
ax1.plot_surface(X1, Y1, data1[:,:,2], cmap='Wistia', alpha = 0.5)
ax1.set_xlabel('C - S bondlength [$\AA$]')
ax1.set_ylabel('C - $H_1$ bondlength [$\AA$]')
ax1.set_zlabel('Energy [eV]')

# Second subplot
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
ax2.plot_surface(X2, Y2, data2[:,:,1], cmap='plasma')
ax2.plot_surface(X2, Y2, data2[:,:,2], cmap='Wistia', alpha = 0.5)
ax2.set_xlabel('C - $H_1$ bondlength [$\AA$]')
ax2.set_ylabel('C - $H_2$ bondlength [$\AA$]')
ax2.set_zlabel('Energy [eV]')

# Third subplot figure 2
ax3 = fig2.add_subplot(1, 2, 1, projection='3d')
ax3.plot_surface(X3, Y3, data3[:,:,1], cmap='plasma')
ax3.plot_surface(X3, Y3, data3[:,:,2], cmap='Wistia', alpha = 0.5)
ax3.set_xlabel('C - S bondlength [$\AA$]')
ax3.set_ylabel(r'$\angle$SC$H_1$ [degrees]')
ax3.set_zlabel('Energy [eV]')

# Fourth subplot figure 2
ax4 = fig2.add_subplot(1, 2, 2, projection='3d')
ax4.plot_surface(X4, Y4, data4[:,:,1], cmap='plasma')
ax4.plot_surface(X4, Y4, data4[:,:,2], cmap='Wistia', alpha = 0.5)
ax4.set_xlabel(r'$\angle$SC$H_2$ [degrees]')
ax4.set_ylabel(r'$\theta$($H_2CS$,$CSH_1$) [degrees]')
ax4.set_zlabel('Energy [eV]')

# Adjust the spacing between subplots
fig.tight_layout()
fig2.tight_layout()

# Display the plot
plt.show()



































