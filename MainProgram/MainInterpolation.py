from Functions_Interpolation import find_interpolant, compute_value_interpolation,integrate_spline,doubleintegrate_spline, \
    quadrupleintegrate_spline,positive
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import time
from shearforces import shearstress
import SVV_structural_properties as prop
from shearcentercalc import get_sc
import Data as data

################DATA FOKKERF100##################
Ca = data.Ca #m
la = data.la #m
x1 = data.x1 #m
x2 = data.x2 #m
x3 = data.x3 #m
xa = data.xa #m
ha = data.ha #m
d1 = data.d1 # m
d3 = data.d3  # m
P = data.P  # N
theta = data.theta #rad
t_sk = data.t_sk   # skin thickness [m]
t_sp = data.t_sp  # spar thickness [m]
n_st = data.n_st     # number of stiffeners [-]
t_st = data.t_st   # thickness of stiffener [m]
h_st = data.h_st   # height of stiffener [m]
w_st = data.w_st  # width of stiffener [m]
t_sk = data.t_sk
E = data.E
G = data.G

z_hat = get_sc()[0] #m
z_cent = prop.z_cent
Izz = prop.I_zz
Iyy = prop.I_yy
J = prop.J
xI = x2 - xa / 2
xII = x2 + xa / 2
################################################

start_time = time.time()

#Opening and reading the file text of the Fokker 100
file_f100 = open("aerodynamicloadf100.dat", "r")
lines = file_f100.readlines()

#Variables
Nz = 81
Nx = 41

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
    coor_z[i-1] = -1/2*(Ca/2*(1-np.cos(theta_z[i-1]))+Ca/2*(1-np.cos(theta_z[i])))

#X-Coordinate
theta_x = np.zeros(Nx+1)
for i in range(1,Nx+2):
    theta_x[i-1] = (i-1)*np.pi/Nx

coor_x = np.zeros(Nx)
for i in range(1,Nx+1):
    coor_x[i-1] = 1/2*(la/2*(1-np.cos(theta_x[i-1]))+la/2*(1-np.cos(theta_x[i])))

#Diff_nodes
diff_x = np.diff(coor_x)
diff_z = np.diff(coor_z)

#Plotting a surface of the new aerodynamic loading
X,Z = np.meshgrid(coor_x,coor_z)
Y = matrix_data

#Plotting the surfaces
plt.figure(1)
cp = plt.contour(X,Z,Y)
plt.colorbar(cp)
ax = plt.axes(projection='3d')
ax.plot_surface(X,Z,Y,cmap='magma')

ax.set_xlabel('X-Axis [m] ~ Spanwise')
ax.set_ylabel('Z-Axis [m] ~ Chordwise ')
ax.set_zlabel('Aerodynamic Loading [kPa]')

#Plotting the wireframe
plt.figure(2)
cp = plt.contour(X,Z,Y)
plt.colorbar(cp)
ax = plt.axes(projection='3d')
ax.plot_wireframe(X,Z,Y,cmap='magma')

ax.set_xlabel('X-Axis [m] ~ Spanwise')
ax.set_ylabel('Z-Axis [m] ~ Chordwise ')
ax.set_zlabel('Aerodynamic Loading [kPa]')

#########INTERPOLATION################
#Plot the interpolation
#Plotting the scatter points with a resolution of maximum of 1 mm
fig = plt.figure(3)
ax = plt.axes(projection = '3d')

for chord in range(0,Nx):
    coeff_matrix=find_interpolant(coor_z,matrix_data[:,chord])
    y_line = []
    z_line = np.linspace(coor_z[0],coor_z[-1],500)
    for step in range(len(z_line)):
        value = compute_value_interpolation(coor_z,coeff_matrix,z_line[step])[0]
        y_line.append(value)
    nodes_length=len(z_line)
    x_line=coor_x[chord]*np.ones(nodes_length)
    ax.scatter3D(x_line,z_line,y_line,marker='o')

for span in range(0,Nz):
    coeff_matrix=find_interpolant(coor_x,matrix_data[span,:])
    y_line = []
    x_line = np.linspace(coor_x[0],coor_x[-1],500)
    for step in range(len(z_line)):
        value = compute_value_interpolation(coor_x,coeff_matrix,x_line[step])[0]
        y_line.append(value)
    nodes_length=len(x_line)
    z_line=coor_z[span]*np.ones(nodes_length)
    ax.scatter3D(x_line,z_line,y_line,marker='^')

ax.set_xlabel('X-Axis [m]')
ax.set_ylabel('Z-Axis [m]')
ax.set_zlabel('Aerodynamic Loading [kPa]')
#plt.show()

####Integration chordwise
Area_chord= []
for chord in range(0,Nx):
    coeff_matrix = find_interpolant(coor_z,matrix_data[:,chord])
    Area_singlechord = -1000*integrate_spline(coor_z,coeff_matrix,coor_z[-1])
    Area_chord.append(Area_singlechord)
Area_chord = np.array(Area_chord)

print(doubleintegrate_spline(coor_x,find_interpolant(coor_x,Area_chord),coor_x[-1]))

center_pressure = []
for chord in range(0,Nx):
    sum_z_load = []
    sum_load = []
    for span in range(0,Nz):
        sum_z_load.append(matrix_data[span,chord]*coor_z[span])
        sum_load.append(matrix_data[span, chord])
    center_pressure.append(np.sum(sum_z_load)/np.sum(sum_load))
center_pressure = np.array(center_pressure)

arm = np.array(center_pressure-z_hat)

torque_chord = []
for chord in range(0,Nx):
    torque_chord.append(Area_chord[chord]*arm[chord])
torque_chord = np.array(torque_chord)

######MainMatrix########

## Row 1: Shear Sy(la)=0
row1 = np.array([1,0,-np.sin(theta),1,0,1,0,0,0,0,0,0])

## Row 2: Shear Sz(la)=0
row2 = np.array([0,1,-np.cos(theta),0,1,0,1,0,0,0,0,0])

## Row 3: My(la)=0
row3 = np.array([0,(la-x1),-np.cos(theta)*(la-(xI)),0,(la-x2),0,(la-x3),0,0,0,0,0])

## Row 4: Mz(la)=0
row4 = np.array([la-x1,0,-np.sin(theta)*(la-(xI)),la-x2,0,la-x3,0,0,0,0,0,0])

## Row 5: T(la)=0
row5 = np.array([-(np.abs(z_hat)-ha/2),0,np.sin(theta)*(np.abs(z_hat)) - np.cos(theta)*ha/2,-(np.abs(z_hat)-ha/2),0,-(np.abs(z_hat)-ha/2),0,0,0,0,0,0])

## Row 6: v(x1) + theta(x1)*(z_hat+ha/2)=d1*cos(theta0)
row6 = np.array([-(1/(6*E*Izz))*(x1-x1)**3 - (1/(G*J))*(z_hat + ha/2)*(x1-x1)*(np.abs(z_hat) - ha/2),0,0,0,0,0,0,0,0,x1,1,(z_hat + ha/2)])

## Row 7: v(x2) + theta(x2)*(z_hat+ha/2)=0
row7 = np.array([-(1/(6*E*Izz))*(x2-x1)**3 - (1/(G*J))*(z_hat + ha/2)*(x2-x1)*(np.abs(z_hat) - ha/2),0,-(1/(6*E*Izz))*-np.sin(theta)*(x2-(xI))**3 + (1/(G*J))*(np.sin(theta)*np.abs(z_hat)*(x2-(xI)) - np.cos(theta)*ha/2*(x2-(xI)))*(z_hat + ha/2),-(1/(6*E*Izz))*(x2-x2)**3 - (1/(G*J))*(np.abs(z_hat) - ha/2)*(x2-x2)*(z_hat + ha/2),0,0,0,0,0,x2,1,(z_hat + ha/2)])

## Row 8: v(x3) + theta(x3)*(z_hat+ha/2)=d3*cos(theta0)
row8 = np.array([-(1/(6*E*Izz))*(x3-x1)**3 - (1/(G*J))*(z_hat + ha/2)*(x3-x1)*(np.abs(z_hat) - ha/2),0,-(1/(6*E*Izz))*-np.sin(theta)*(x3-(xI))**3 + (1/(G*J))*(np.sin(theta)*np.abs(z_hat)*(x3-(xI)) - np.cos(theta)*ha/2*(x3-(xI)))*(z_hat + ha/2),-(1/(6*E*Izz))*(x3-x2)**3 - (1/(G*J))*(np.abs(z_hat) - ha/2)*(x3-x2)*(z_hat + ha/2),0,-(1/(6*E*Izz))*(x3-x3)**3 - (1/(G*J))*(z_hat + ha/2)*(x3-x3)*(np.abs(z_hat) - ha/2),0,0,0,x3,1,(z_hat + ha/2)])

## Row 9: w(x1) = -d1*sin(theta0)
row9 = np.array([0,-(1/(6*E*Iyy))*(x1-x1)**3,0,0,0,0,0,x1,1,0,0,0])

## Row 10: w(x2) = 0
row10 = np.array([0,-(1/(6*E*Iyy))*(x2-x1)**3,(1/(6*E*Iyy))*np.cos(theta)*(x2-(xI))**3,0,-(1/(6*E*Iyy))*(x2-x2)**3,0,0,x2,1,0,0,0])

## Row 11: w(x3) = -d3*sin(theta0)
row11 = np.array([0,-(1/(6*E*Iyy))*(x3-x1)**3,(1/(6*E*Iyy))*np.cos(theta)*(x3-(xI))**3,0,-(1/(6*E*Iyy))*(x3-x2)**3,0,-(1/(6*E*Iyy))*(x3-x3)**3,x3,1,0,0,0])

## Row 12: w(xI)*cos(theta) + v(xI)*sin(theta) - twist(xI)*abs(zsc)*sin(theta) = 0
row12 = np.array([-(1/(6*E*Izz))*np.sin(theta)*(xI-x1)**3 + (1/(G*J))*np.sin(theta)*np.abs(z_hat)*(np.abs(z_hat)-ha/2),-(1/(6*E*Iyy))*np.cos(theta)*(xI-x1)**3,-(1/(6*E*Iyy))*np.cos(theta)*np.cos(theta)*(xI-xI)**3 - (1/(6*E*Izz))*np.sin(theta)*np.sin(theta)*(xI-xI)**3 + (1/(G*J))*np.abs(z_hat)*np.sin(theta)*np.sin(theta)*np.abs(z_hat)*(xI-xI),0,0,0,0,xI*np.cos(theta),np.cos(theta),xI*np.sin(theta),np.sin(theta),-np.abs(z_hat)*np.sin(theta)])

A_matrix = np.array([row1, row2, row3, row4, row5, row6, row7, row8, row9, row10, row11, row12])

B_matrix = np.array([[P*np.sin(theta) + integrate_spline(coor_x,find_interpolant(coor_x,Area_chord),coor_x[-1])], #+q
			   [P*np.cos(theta)],
			   [P*np.cos(theta)*(la-(xII))],
			   [P*np.sin(theta)*(la-(xII))+doubleintegrate_spline(coor_x,find_interpolant(coor_x,Area_chord),coor_x[-1])], #+q
			   [P*(np.cos(theta)*ha/2 - np.abs(z_hat)*np.sin(theta))+integrate_spline(coor_x,find_interpolant(coor_x,torque_chord),coor_x[-1])],#-q
			   [d1*np.cos(theta)-1/(E*Izz)*quadrupleintegrate_spline(coor_x,find_interpolant(coor_x,Area_chord),x1)+1/(G*J)*doubleintegrate_spline(coor_x,find_interpolant(coor_x,torque_chord),x1)*(z_hat+ha/2)],#+(1/(E*Izz))*q -(1/(G*J))*q Add AeroTorque
			   [-1/(E*Izz)*quadrupleintegrate_spline(coor_x,find_interpolant(coor_x,Area_chord),x2)+1/(G*J)*doubleintegrate_spline(coor_x,find_interpolant(coor_x,torque_chord),x2)*(z_hat+ha/2)], #[(1/(E*Izz))*q - (1/(G*J))*q] Add Aero
			   [d3*np.cos(theta) - (1/(6*E*Izz))*P*np.sin(theta)*(x3-(xII))**3 - (1/(G*J))*P*(z_hat + ha/2)*(x3-(xII))*(np.abs(z_hat)*np.sin(theta)-ha/2*np.cos(theta))-1/(E*Izz)*quadrupleintegrate_spline(coor_x,find_interpolant(coor_x,Area_chord),x3)+1/(G*J)*doubleintegrate_spline(coor_x,find_interpolant(coor_x,torque_chord),x3)*(z_hat+ha/2)], # + (1/(E*Izz))*q  - (1/(G*J))*q Add Aero
			   [-d1*np.sin(theta)],
			   [0],
			   [-d3*np.sin(theta) - (1/(6*E*Iyy))*P*np.cos(theta)*(x3-(xII))**3],
			   [-1/(E*Izz)*quadrupleintegrate_spline(coor_x,find_interpolant(coor_x,Area_chord),xI)*np.sin(theta)-1/(G*J)*np.abs(z_hat)*np.sin(theta)*doubleintegrate_spline(coor_x,find_interpolant(coor_x,torque_chord),xI)]]) #-(1/(G*J))*q Add Aero

x_matrix = np.dot(np.linalg.inv(A_matrix), B_matrix)

def Sy(x):
    return float(x_matrix[0] * positive(x - x1, 0) - x_matrix[2] * np.sin(theta) * positive(x - (xI), 0) + x_matrix[3] * positive(x - x2,0) - P * np.sin(theta) * positive(x - (xII), 0) + x_matrix[5] * positive(x - x3, 0) - integrate_spline(coor_x,find_interpolant(coor_x,Area_chord),x))

def Sz(x):
    return float(x_matrix[1] * positive(x - x1, 0) - x_matrix[2] * np.cos(theta) * positive(x - (xI), 0) + x_matrix[4] * positive(x - x2,0) - P * np.cos(theta) * positive(x - (xII), 0) + x_matrix[6] * positive(x - x3, 0))

def Moment_y(x):
    return float(x_matrix[1] * positive(x - x1, 1) - x_matrix[2] * np.cos(theta) * positive(x - (xI), 1) + x_matrix[4] * positive(x - x2,1) - P * np.cos(theta) * positive(x - (xII), 1) + x_matrix[6] * positive(x - x3, 1))

def Moment_z(x):
    return float(x_matrix[0] * positive(x - x1, 1) - x_matrix[2] * np.sin(theta) * positive(x - (xI), 1) + x_matrix[3] * positive(x - x2,1) - P * np.sin(theta) * positive(x - (xII), 1) + x_matrix[5] * positive(x - x3, 1) - doubleintegrate_spline(coor_x,find_interpolant(coor_x,Area_chord),x))

def T(x):
    return float(-x_matrix[0] * (np.abs(z_hat) - ha / 2) * positive(x - x1, 0) + x_matrix[2] * (np.sin(theta) * np.abs(z_hat) - np.cos(theta) * ha / 2) * positive(x - (xI), 0) - x_matrix[3] * (np.abs(z_hat) - ha / 2) * positive(x - x2, 0) + P * (np.sin(theta) * np.abs(z_hat) - np.cos(theta) * ha / 2) * positive(x - (xII), 0) - x_matrix[5] * (np.abs(z_hat) - ha / 2) * positive(x - x3, 0) - integrate_spline(coor_x,find_interpolant(coor_x,torque_chord),x))

def v(x):
    return float(-(1 / (6 * E * Izz)) * (x_matrix[0] * positive(x - x1, 3) - x_matrix[2] * np.sin(theta) * positive(x - (xI), 3) + x_matrix[3] * positive(x - x2,3) - P * np.sin(theta) * positive(x - (xII), 3) + x_matrix[5] * positive(x - x3, 3) - 6*quadrupleintegrate_spline(coor_x,find_interpolant(coor_x,Area_chord),x)) + x_matrix[9] * x + x_matrix[10])

def w(x):
    return float(-(1 / (6 * E * Iyy)) * (x_matrix[1] * positive(x - x1, 3) - x_matrix[2] * np.cos(theta) * positive(x - (xI), 3) + x_matrix[4] * positive(x - x2,3) - P * np.cos(theta) * positive(x - (xII), 3) + x_matrix[6] * positive(x - x3, 3)) + x_matrix[7] * x + x_matrix[8])

def Twist(x):
    return float((1 / (G * J)) * (-x_matrix[0] * (np.abs(z_hat) - ha / 2) * positive(x - x1, 1) + x_matrix[2] * (np.sin(theta) * np.abs(z_hat) - np.cos(theta) * ha / 2) * positive(x - (xI), 1) - x_matrix[3] * (np.abs(z_hat) - ha / 2) * positive(x - x2, 1) + P * (np.sin(theta) * np.abs(z_hat) - np.cos(theta) * ha / 2) * positive(x - (xII), 1) -x_matrix[5] * (np.abs(z_hat) - ha / 2) * positive(x - x3, 1) - doubleintegrate_spline(coor_x,find_interpolant(coor_x,torque_chord),x)) - x_matrix[11])
