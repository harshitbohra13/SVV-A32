import numpy as np
import matplotlib.pyplot as plt
from Functions_Interpolation import find_interpolant,compute_interpolant

#Opening and reading the file text of the Fokker 100
file_f100 = open("aerodynamicloadf100.dat","r")
lines = file_f100.readlines()

#Variables
Nz = 81
Nx = 41
C_a = 0.505
l_a = 1.611

#Inserting the values of the text file to a matrix 81 times 41
matrix_data = np.zeros((Nz,Nx))
for line in lines:
    row = lines.index(line)
    line = line.replace("\n","")
    magnitude_list = line.split(",")
    matrix_data[row,:] = magnitude_list

#Z-Coordinate
theta_z = np.zeros(Nz+1)
for i in range(1,Nz+2):
    theta_z[i-1] = (i-1)*np.pi/Nz

coor_z = np.zeros(Nz)
for i in range(1,Nz+1):
    coor_z[i-1] = -1/2*(C_a/2*(1-np.cos(theta_z[i-1]))+C_a/2*(1-np.cos(theta_z[i])))

#X-Coordinate
theta_x = np.zeros(Nx+1)
for i in range(1,Nx+2):
    theta_x[i-1] = (i-1)*np.pi/Nx

coor_x = np.zeros(Nx)
for i in range(1,Nx+1):
    coor_x = 1/2*(l_a/2*(1-np.cos(theta_x[i-1]))+l_a/2*(1-np.cos(theta_x[i])))

#Find Coefficients of the first line chordwise
coeff_matrix=find_interpolant(coor_z,matrix_data[:,0])

#Plotting the splines
new_nodes,new_loading = compute_interpolant(coor_z,coeff_matrix,0.01)
plt.plot(new_nodes,new_loading,'r-')

#Plotting the points given from the .dat file
plt.plot(coor_z,matrix_data[:,0],'o')

plt.grid()
plt.show()