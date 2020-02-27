import numpy as np
from Functions_Interpolation import compute_value_interpolation, find_interpolant, integrate_spline
import matplotlib.pyplot as plt

#Opening and reading the file text of the Fokker 100
file_f100 = open("aerodynamicloadf100.dat", "r")
lines = file_f100.readlines()


#Variables
Ca = 0.505     # chord length aileron [m]
theta = 30*np.pi/180 #rad
Nz = 81
Nx = 41

#Inserting the values of the text file to a matrix 81 times 41
matrix_data = np.zeros((Nz,Nx))
for line in lines:
    row = lines.index(line)
    line = line.replace("\n","")
    magnitude_list = line.split(",")
    matrix_data[row,:] = magnitude_list

#Aerodynamic data
test_data = matrix_data[:,20] #chordline

#z-coordinates
theta_z = np.zeros(Nz+1)
for i in range(1,Nz+2):
    theta_z[i-1] = (i-1)*np.pi/Nz

coor_z = np.zeros(Nz)
for i in range(1,Nz+1):
    coor_z[i-1] = -1/2*(Ca/2*(1-np.cos(theta_z[i-1]))+Ca/2*(1-np.cos(theta_z[i])))

#computing the interpoling coefficients
Coeff_matrix = find_interpolant(coor_z,test_data)

# computing the spline integral
area = integrate_spline(coor_z,Coeff_matrix,coor_z[-1])
print(area)


plt.plot(coor_z, test_data)
plt.show()


