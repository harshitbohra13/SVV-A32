import numpy as np

Ca = 0.505 #m
la = 1.611 #m
x1 = 0.125 #m
x2 = 0.498 #m
x3 = 1.494 #m
xa = 24.5*10**(-2) #m
ha = 0.161 #m
d1 = 0.389*10**(-2) #m
d3 = 1.245*10**(-2) #m
P = 49.2*10**(3) #N
theta = 30*np.pi()/180 #rad
z_centroid = -0.1 #m
Izz = 1 #
Iyy = 1 #
J = 1 #
E = 72.9*10**9
G = 27.1*10**9

#X matrix
#x_matrix = np.transpose(np.array([R1y,R1z,Fa,R2y,R2z,R3y,R3z,C1,C2,C3,C4,C5]))

#B_matrix
b_matrix = np.array([[-P*np.cos(theta)*(la-x2-xa/2)],
                    [-P*np.sin(theta)*(la-x2-xa/2)],#- moment q
                    [-P*np.sin(theta)*(0-z_centroid)*1+P*np.cos(0-z_centroid)*1],# - torque q
                    [-P*np.cos(theta)],
                    [-P*np.sin(theta)],# - shear q
                    [d1*np.sin(theta)],
                    [0],
                    [0],
                    [d3*np.sin(theta)+P/(6*E*Iyy)*np.cos(theta)*(x3-x2-xa/2)**3],
                    [d1*np.cos(theta)], #+ v and theta in q
                     [0], #+ v and theta in q
                     [-P/(G*J)*z_centroid*np.sin(theta)*(x3-x2-xa/2)*ha/2+P*np.sin(theta)/(6*E*Izz)*(x-x2-xa/2)**3] #+ v and theta in q
                     ])

#Moment My(la)=0
row_1 = np.array([0,-(la-x1),-np.cos(theta)*(la-x2+xa/2),0,-(la-x2),0,-(la-x3),0,0,0,0,0])

#Moment Mz(la)=0
row_2 = np.array([-(la-x1),0,-np.sin(theta)*(la-x2+xa/2),-(la-x2),0,-(la-x3),0,0,0,0,0,0])

#Torque T(la) = 0
row_3 = np.array([-(-ha/2-z_centroid)*1,0,(-np.sin(theta)*(0-z_centroid)+np.cos(theta)*ha/2)*1,-(-ha/2-z_centroid)*1,0,-(-ha/2-z_centroid)*1,0,0,0,0,0,0])

#Shear Sz(la)=0
row_4 = np.array([0,-1,-np.cos(theta),0,-1,0,-1,0,0,0,0,0])

#Shear Sy(la)=0
row_5 = np.array([-1,0,-np.sin(theta),-1,0,-1,0,0,0,0,0,0])

#Curvature w(x1) = D1*sin(theta)
row_6 = np.array([0,0,0,0,0,0,0,0,0,x1,1,0])

#Curvature w(x2-xa/2) = 0
row_7 = np.array([0,1/(6*E*Iyy)*(x2-xa/2-x1)**3,0,0,0,0,0,0,0,x2-xa/2,1,0])

#Curvature w(x2) = 0
row_8 = np.array([0,1/(6*E*Iyy)*(x2-x1)**3,1/(6*E*Iyy)*np.cos(theta)*(x2-x2+xa/2)**3,0,0,0,0,0,0,x2,1,0])

#Curvature w(x3) = D3*sin(theta)
row_9 = np.array([0,1/(6*E*Iyy)*(x3-x1)**3,1/(6*E*Iyy)*np.cos(theta)*(x3-x2+xa/2)**3,0,1/(6*E*Iyy)*(x3-x2)**3,0,0,0,0,x3,1,0])

#Curvature v(x1)+theta(x1)*(z_centroid)=d1*cos(theta)
row_10 = np.array([0,0,0,0,0,0,0,x1,1,0,0,z_centroid])

#Curvature v(x2)+theta(x2)*(z_centroid) = 0
row_11 = np.array([1/(6*E*Izz)*(x2-x1)**3-z_centroid/(G*J)*(-ha/2-z_centroid)*(x2-x1),0,
                   1/(6*E*Izz)*(x2-x2+xa/2)**3+z_centroid/(G*J)*(-np.sin(theta)(0-z_centroid)+np.cos(theta)*(ha/2))*(x2-x2+xa/2),
                  0,0,0,0,x2,1,0,0,z_centroid])

#Curvature v(x3)+theta(x3)*(z_centroid) = 0
row_12 = np.array([1/(6*E*Izz)*(x3-x1)**3-z_centroid/(G*J)*(-ha/2-z_centroid)*(x3-x1),0,
                   1/(6*E*Izz)*(x3-x2+xa/2)**3+z_centroid/(G*J)*(-np.sin(theta)(0-z_centroid)+np.cos(theta)*(ha/2))*(x3-x2+xa/2),
                  1/(6*E*Izz)*(x3-x2)**3-z_centroid/(G*J)*(-ha/2-z_centroid)*(x3-x2),0,
                   0,0,x3,1,0,0,z_centroid])

Big_matrix = np.array([row_1,row_2,row_3,row_4,row_5,row_6,row_7,row_8,row_9,row_10,row_11,row_12])

x_matrix = np.linalg.tensorsolve(Big_matrix,b_matrix)
