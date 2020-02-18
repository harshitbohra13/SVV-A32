import numpy as np
import matplotlib.pyplot as plt

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

coor_x = np.zeros(Nz)
for i in range(1,Nx+1):
    coor_x = 1/2*(l_a/2*(1-np.cos(theta_x[i-1]))+l_a/2*(1-np.cos(theta_x[i])))

#Interpolation for the first chord (there are 41)
magnitude_array = np.array(matrix_data[:,0]) # 0 is the actual column, we will have to do it for each column

#Computes the function difference (f_{i+1}-f_{i})
diff_magnitude = np.diff(magnitude_array)

#Computes the step-size h_{i}=z_{i+1}-z_{i}
diff_z=np.diff(coor_z)

#Computes the h_{i-1}/6*M_{i-1}+(h_{i-1}+h_{i})/3*M_{i}+h_{i}/6*M_{i+1}=(f_{i+1}-f_{i})\h_{i}(f_{i}-f_{i-1})/h_{i-1}
a_matrix = np.zeros((Nz-2,Nz-2)) #left hand side
b_matrix = np.zeros(Nz-2) #right hand side

for row in range(0,Nz-2):
    a_matrix[row,row] = (diff_z[row]+diff_z[row+1])/3
    if row!=Nz-3:
        a_matrix[row,row+1] = diff_z[row+1]/6
    if row!=0:
        a_matrix[row,row-1] = diff_z[row]/6

    b_matrix[row] = (magnitude_array[row+2]-magnitude_array[row+1])/diff_z[row+1]-(magnitude_array[row+1]-magnitude_array[row])/diff_z[row]

#Calculate the Constraint Matrix (M_0,M_1,M_2) --> A*x=B
m_matrix = np.linalg.tensorsolve(a_matrix,b_matrix)

#Adds M_{0} and M_{n} equal to 0 (just adding two elements, one at the beginning and one at the end)
m_matrix = np.insert(np.append(m_matrix,[0]),0,0)

#Computes coefficients a_{i}, b_{i}, c_{i} and d_{i}
coeff_a_matrix = np.diff(m_matrix)/(6*diff_z)

coeff_b_matrix = m_matrix/2
coeff_b_matrix = np.delete(coeff_b_matrix,-1)

coeff_c_matrix = np.zeros(Nz-1)
for element in range(0,Nz-1):
    coeff_c_matrix[element]=diff_magnitude[element]/diff_z[element]-diff_z[element]/3*m_matrix[element]-diff_z[element]/6*m_matrix[element+1]

coeff_d_matrix = np.delete(np.array(magnitude_array),-1)

#Coeff_matrix is a matrix containing all the coefficients
coeff_matrix = [coeff_a_matrix,coeff_b_matrix,coeff_c_matrix,coeff_d_matrix]
coeff_matrix= np.transpose(np.array(coeff_matrix))

#Plotting the points given from the .dat file
plt.plot(coor_z,magnitude_array,'o')
x_axis = []
y_axis = []

#Plotting the splines
for j in range(0,coeff_matrix.shape[0]):
    z_i_spline=coor_z[j]
    diff_z_spline=diff_z[j]
    step_spline=np.linspace(0,diff_z[j],2) #Taking 100 steps for each spline
    for step in range(len(step_spline)-1):
        spline_function = coeff_matrix[j,0]*(step_spline[step])**3+coeff_matrix[j,1]*(step_spline[step])**2+coeff_matrix[j,2]*(step_spline[step])+coeff_matrix[j,3]
        y_axis.append(spline_function)
        x_axis.append(z_i_spline+step_spline[step])
plt.plot(x_axis,y_axis,'r')
plt.grid()
plt.show()