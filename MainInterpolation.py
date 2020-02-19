from Functions_Interpolation import find_interpolant,compute_interpolant
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()

#Opening and reading the file text of the Fokker 100
file_f100 = open("aerodynamicloadf100.dat","r")
lines = file_f100.readlines()

#Variables
Nz = 81
Nx = 41
C_a = 0.505
l_a = 1.611
resolution_span = 1 #mm
resolution_chord = 1 #mm

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
    coor_x[i-1] = 1/2*(l_a/2*(1-np.cos(theta_x[i-1]))+l_a/2*(1-np.cos(theta_x[i])))

#Finds the maximum spacing between nodes
max_diff_coor_z = np.amax(np.abs(np.diff(coor_z)))

#Number of steps needed between two MAIN NODES
internal_z_points= int(np.ceil(max_diff_coor_z/(resolution_chord*10**(-3))))
number_z = internal_z_points*(Nz)

#Finds the maximum spacing between nodes
max_diff_coor_x = np.amax(np.abs(np.diff(coor_x)))

#Number of steps needed between two MAIN NODES
internal_x_points=int(np.ceil(max_diff_coor_x/(resolution_span*10**(-3))))
number_x = internal_x_points*(Nx)

#Create a new aero_data
new_matrix_data= np.zeros((number_z,Nx))

#Increasing the grid-mesh chordwise
for chord in range(0,Nx):
    coeff_matrix_chord=find_interpolant(coor_z,matrix_data[:,chord])
    new_nodes_z,y_line_z = compute_interpolant(coor_z,coeff_matrix_chord,resolution_chord)
    new_matrix_data[:,chord]=y_line_z

#Create a new aero_data
aero_data= np.zeros((number_z,number_x))

#Increasing the grid-mesh chordwise
for span in range(0,number_z):
    coeff_matrix_span=find_interpolant(coor_x,new_matrix_data[span,:])
    new_nodes_x,y_line_x = compute_interpolant(coor_x,coeff_matrix_span,resolution_span)
    aero_data[span,:]=y_line_x

#Plotting a surface of the new aerodynamic loading
X,Z = np.meshgrid(new_nodes_x,new_nodes_z)
Y=aero_data

#Plotting the surfaces
plt.figure(1)
cp = plt.contour(X,Z,Y)
plt.colorbar(cp)
ax = plt.axes(projection='3d')
ax.plot_surface(X,Z,Y,cmap='magma')

ax.set_xlabel('X-Axis [m] ~ Spanwise')
ax.set_ylabel('Z-Axis [m] ~ Chordwise ')
ax.set_zlabel('Aerodynamic Loading [kPa]')

plt.show()

#Plotting the wireframe
plt.figure(2)
cp = plt.contour(X,Z,Y)
plt.colorbar(cp)
ax = plt.axes(projection='3d')
ax.plot_wireframe(X,Z,Y,cmap='magma')

ax.set_xlabel('X-Axis [m] ~ Spanwise')
ax.set_ylabel('Z-Axis [m] ~ Chordwise ')
ax.set_zlabel('Aerodynamic Loading [kPa]')

plt.show()

print("Runtime: %f seconds" % (time.time()-start_time))



