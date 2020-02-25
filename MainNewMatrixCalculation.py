import numpy as np
import math as m
import matplotlib.pyplot as plt
import polarmoment as po
import SVV_structural_properties as prop

Ca = 0.505 #m
la = 1.611 #m
x1 = 0.125 #m
x2 = 0.498 #m
x3 = 1.494 #m
xa = 0.245 #m
ha = 0.161 #m
d1 = 0.00389 # m
d3 = 0.01245  # m
P = 49.2*1000  # N
theta = 30*np.pi/180 #rad
z_hat = -0.08553893540215983 #m
Izz = 4.753851442684436e-06
#prop.I_zz #
Iyy = 4.5943507864451845e-05
#prop.I_yy #
J = 7.748548555816593e-06
#po.J #
E = 72.9*10**9
G = 27.1*10**9
xI = x2 - xa / 2
xII = x2 + xa / 2

# Order of the results matrix
# R=[R1y,R1z,RI,R2y,R2z,R3y,R3z,C1,C2,C3,C4,C5]

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

## Row 6: Boundary condition 1  v(x1) + theta(x1)*(z_hat+ha/2)=d1*cos(theta0)
row6 = np.array([-(1/(6*E*Izz))*(x1-x1)**3 - (1/(G*J))*(z_hat + ha/2)*(x1-x1)*(np.abs(z_hat) - ha/2),0,0,0,0,0,0,0,0,x1,1,(z_hat + ha/2)])

## Row 7: Boundary condition 2  vy(x2) + theta(x2)*(z_hat+ha/2)=0
row7 = np.array([-(1/(6*E*Izz))*(x2-x1)**3 - (1/(G*J))*(z_hat + ha/2)*(x2-x1)*(np.abs(z_hat) - ha/2),0,-(1/(6*E*Izz))*-np.sin(theta)*(x2-(xI))**3 + (1/(G*J))*(np.sin(theta)*np.abs(z_hat)*(x2-(xI)) - np.cos(theta)*ha/2*(x2-(xI)))*(z_hat + ha/2),-(1/(6*E*Izz))*(x2-x2)**3 - (1/(G*J))*(np.abs(z_hat) - ha/2)*(x2-x2)*(z_hat + ha/2),0,0,0,0,0,x2,1,(z_hat + ha/2)])

## Row 8: Boundary condition 3  vy(x3) + theta(x3)*(z_hat+ha/2)=d3*cos(theta0)
row8 = np.array([-(1/(6*E*Izz))*(x3-x1)**3 - (1/(G*J))*(z_hat + ha/2)*(x3-x1)*(np.abs(z_hat) - ha/2),0,-(1/(6*E*Izz))*-np.sin(theta)*(x3-(xI))**3 + (1/(G*J))*(np.sin(theta)*np.abs(z_hat)*(x3-(xI)) - np.cos(theta)*ha/2*(x3-(xI)))*(z_hat + ha/2),-(1/(6*E*Izz))*(x3-x2)**3 - (1/(G*J))*(np.abs(z_hat) - ha/2)*(x3-x2)*(z_hat + ha/2),0,-(1/(6*E*Izz))*(x3-x3)**3 - (1/(G*J))*(z_hat + ha/2)*(x3-x3)*(np.abs(z_hat) - ha/2),0,0,0,x3,1,(z_hat + ha/2)])

## Row 9: Boundary condition 4  vz(x1) = -d1*sin(theta0)
row9 = np.array([0,-(1/(6*E*Iyy))*(x1-x1)**3,0,0,0,0,0,x1,1,0,0,0])

## Row 10: Boundary condition 5 vz(x2) = 0
row10 = np.array([0,-(1/(6*E*Iyy))*(x2-x1)**3,(1/(6*E*Iyy))*np.cos(theta)*(x2-(xI))**3,0,-(1/(6*E*Iyy))*(x2-x2)**3,0,0,x2,1,0,0,0])

## Row 11: Boundary condition 6 vz(x3) = -d3*sin(theta0)
row11 = np.array([0,-(1/(6*E*Iyy))*(x3-x1)**3,(1/(6*E*Iyy))*np.cos(theta)*(x3-(xI))**3,0,-(1/(6*E*Iyy))*(x3-x2)**3,0,-(1/(6*E*Iyy))*(x3-x3)**3,x3,1,0,0,0])

## Row 12: Boundary condition 7 w(xI)*cos(theta) + v(xI)*sin(theta) - twist(xI)*abs(zsc)*sin(theta) = 0
row12 = np.array([-(1/(6*E*Izz))*np.sin(theta)*(xI-x1)**3 + (1/(G*J))*np.sin(theta)*np.abs(z_hat)*(np.abs(z_hat)-ha/2),-(1/(6*E*Iyy))*np.cos(theta)*(xI-x1)**3,-(1/(6*E*Iyy))*np.cos(theta)*np.cos(theta)*(xI-xI)**3 - (1/(6*E*Izz))*np.sin(theta)*np.sin(theta)*(xI-xI)**3 + (1/(G*J))*np.abs(z_hat)*np.sin(theta)*np.sin(theta)*np.abs(z_hat)*(xI-xI),0,0,0,0,xI*np.cos(theta),np.cos(theta),xI*np.sin(theta),np.sin(theta),-np.abs(z_hat)*np.sin(theta)])

A_matrix = np.array([row1, row2, row3, row4, row5, row6, row7, row8, row9, row10, row11, row12])

B_matrix = np.array([[P*np.sin(theta)], #+q
			   [P*np.cos(theta)],
			   [P*np.cos(theta)*(la-(xII))],
			   [P*np.sin(theta)*(la-(xII))], #+q
			   [P*(np.cos(theta)*ha/2 - np.abs(z_hat)*np.sin(theta))],#-q
			   [d1*np.cos(theta)],#+(1/(E*Izz))*q -(1/(G*J))*q Add AeroTorque
			   [0], #[(1/(E*Izz))*q - (1/(G*J))*q] Add Aero
			   [d3*np.cos(theta) - (1/(6*E*Izz))*P*np.sin(theta)*(x3-(xII))**3 - (1/(G*J))*P*(z_hat + ha/2)*(x3-(xII))*(np.abs(z_hat)*np.sin(theta)-ha/2*np.cos(theta)) ], # + (1/(E*Izz))*q  - (1/(G*J))*q Add Aero
			   [-d1*np.sin(theta)],
			   [0],
			   [-d3*np.sin(theta) - (1/(6*E*Iyy))*P*np.cos(theta)*(x3-(xII))**3],
			   [0]]) #-(1/(G*J))*q Add Aero

x_matrix = np.dot(np.linalg.inv(A_matrix), B_matrix)


def positive(x, power):
    if power > 0 and x > 0:
        return x ** power
    elif power == 0 and x > 0:
        return 1
    else:
        return 0


def Sy(x):
    return x_matrix[0] * positive(x - x1, 0) - x_matrix[2] * np.sin(theta) * positive(x - (xI), 0) + x_matrix[3] * positive(x - x2,0) - P * np.sin(theta) * positive(x - (xII), 0) + x_matrix[5] * positive(x - x3, 0)

def Sz(x):
    return x_matrix[1] * positive(x - x1, 0) - x_matrix[2] * np.cos(theta) * positive(x - (xI), 0) + x_matrix[4] * positive(x - x2,0) - P * np.cos(theta) * positive(x - (xII), 0) + x_matrix[6] * positive(x - x3, 0)

def Moment_y(x):
    return x_matrix[1] * positive(x - x1, 1) - x_matrix[2] * np.cos(theta) * positive(x - (xI), 1) + x_matrix[4] * positive(x - x2,1) - P * np.cos(theta) * positive(x - (xII), 1) + x_matrix[6] * positive(x - x3, 1)

def Moment_z(x):
    return x_matrix[0] * positive(x - x1, 1) - x_matrix[2] * np.sin(theta) * positive(x - (xI), 1) + x_matrix[3] * positive(x - x2,1) - P * np.sin(theta) * positive(x - (xII), 1) + x_matrix[5] * positive(x - x3, 1)

def T(x):
    return -x_matrix[0] * (np.abs(z_hat) - ha / 2) * positive(x - x1, 0) + x_matrix[2] * (np.sin(theta) * np.abs(z_hat) - np.cos(theta) * ha / 2) * positive(x - (xI), 0) - x_matrix[3] * (np.abs(z_hat) - ha / 2) * positive(x - x2, 0) + P * (np.sin(theta) * np.abs(z_hat) - np.cos(theta) * ha / 2) * positive(x - (xII), 0) - x_matrix[5] * (np.abs(z_hat) - ha / 2) * positive(x - x3, 0)

def v(x):
    return -(1 / (6 * E * Izz)) * (x_matrix[0] * positive(x - x1, 3) - x_matrix[2] * np.sin(theta) * positive(x - (xI), 3) + x_matrix[3] * positive(x - x2,3) - P * np.sin(theta) * positive(x - (xII), 3) + x_matrix[5] * positive(x - x3, 3)) + x_matrix[9] * x + x_matrix[10]

def w(x):
    return -(1 / (6 * E * Iyy)) * (x_matrix[1] * positive(x - x1, 3) - x_matrix[2] * np.cos(theta) * positive(x - (xI), 3) + x_matrix[4] * positive(x - x2,3) - P * np.cos(theta) * positive(x - (xII), 3) + x_matrix[6] * positive(x - x3, 3)) + x_matrix[7] * x + x_matrix[8]

def Twist(x):
    return (1 / (G * J)) * (-x_matrix[0] * (np.abs(z_hat) - ha / 2) * positive(x - x1, 1) + x_matrix[2] * (np.sin(theta) * np.abs(z_hat) - np.cos(theta) * ha / 2) * positive(x - (xI), 1) - x_matrix[3] * (np.abs(z_hat) - ha / 2) * positive(x - x2, 1) + P * (np.sin(theta) * np.abs(z_hat) - np.cos(theta) * ha / 2) * positive(x - (xII), 1) -x_matrix[5] * (np.abs(z_hat) - ha / 2) * positive(x - x3, 1)) - x_matrix[11]

plt.figure(2)
xcheck = np.linspace(0, la, 100)

Shearplot = np.array([])
for i in range(len(xcheck)):
    Syplot = Moment_y(xcheck[i])
    Shearplot = np.append(Shearplot, Syplot)

plt.plot(xcheck, Shearplot)
plt.grid()
plt.show()
