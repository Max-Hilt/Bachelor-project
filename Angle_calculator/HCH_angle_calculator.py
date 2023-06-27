# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 16:30:28 2023

@author: Max Hilt


This function is specially written to find the H-C-H angle.

Since the Z-matrix does not have this angle we calculate it using this small script.
"""

import numpy as np

# Defining position vectors

# Ground state
C_B3LYP_gr = np.array([0,0,0])
S_B3LYP_gr = np.array([1.60797,0.00000,0.00000])
H1_B3LYP_gr = np.array([-0.54264,0.94003,0.00000])
H2_B3LYP_gr = np.array([-0.54265,-0.94003,-0.00000])


# Excited state
C_B3LYP_ex = np.array([0,0,0])
S_B3LYP_ex = np.array([1.69790,0.00000,0.00000])
H1_B3LYP_ex = np.array([-0.54000,0.93545,0.00000])
H2_B3LYP_ex = np.array([-0.54001,-0.93545,-0.00000])

########### We start with calculating the H-C-H angle #############
# Now we define our vectors
a1 = H1_B3LYP_gr - C_B3LYP_gr
a2 = H2_B3LYP_gr - C_B3LYP_gr

b1 = H1_B3LYP_ex - C_B3LYP_ex
b2 = H2_B3LYP_ex - C_B3LYP_ex

# Then we calculate the dotproduct

dot1 = np.dot(a1,a2)
dot2 = np.dot(b1,b2)

# Finally we find the angle
angle1 = np.arccos(dot1 / (np.linalg.norm(a1)*np.linalg.norm(a2)))/(2*np.pi)*360
angle2 = np.arccos(dot2 / (np.linalg.norm(b1)*np.linalg.norm(b2)))/(2*np.pi)*360

print(f'The groundstate angle is {angle1} the excited state angle is {angle2}')

########## Now for the C-S out of plane angle #########

# We can re use the vectors already defined and then calculate the cross product
# To find the normal vector.


norm1 = np.cross(a1, a2)
norm2 = np.cross(b1, b2)

# After the normal vector is found we can do the exact same calculation as
# above once more.

dot_di_1 = np.dot(norm1, S_B3LYP_gr)
dot_di_2 = np.dot(norm2, S_B3LYP_ex)

# Note the 90 degrees and the minus sign since we calculate the angle
# With the normal vector whilst we want the angle with the plane
out_of_plane_angle_1 = 90 - np.arccos(dot_di_1 /(np.linalg.norm(norm1)*np.linalg.norm(S_B3LYP_gr)))/(2*np.pi)*360
out_of_plane_angle_2 = 90 - np.arccos(dot_di_2 /(np.linalg.norm(norm2)*np.linalg.norm(S_B3LYP_ex)))/(2*np.pi)*360

print(f'The ground state dihedral angle is {out_of_plane_angle_1} and the excited state dihedral angle is {out_of_plane_angle_2}')


