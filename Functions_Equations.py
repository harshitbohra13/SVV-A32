import numpy as np

Ca = 0.505 #m
la = 1.611 #m
x1 = 0.125 #m
x2 = 0.498 #m
x3 = 1.494 #m
xa = 24.5*10**2 #m
ha = 0.161 #m
P = 49.2*10**3 #N
theta = 30*np.pi/180 #rad
z_centroid = -0.1 #m
Izz = 1 #
Iyy = 1 #
J = 1 #
E = 72.9*10**9
G = 27.1*10**9

#Unknowns
R1y = 0
R1z = 0
Ractuator = 0
R2y = 0
R2z = 0
R3y = 0
R3z = 0
C1 = 0
C2 = 0
C3 = 0
C4 = 0
C5 = 0

def positive(n,power):
    if power>0 and n>0:
        return n**power
    elif power==0 and n>0:
        return 1
    else:
        return 0


def Moment_y(x):
    return -R1z*positive(x-x1,1)+P*np.cos(theta)*positive(x-x2-xa/2,1)-R2z*positive(x-x2,1)-R3z*positive(x-x3,1)-Ractuator*np.cos(theta)*positive(x-x2+xa/2,1)

def Moment_z(x): #+q
    return -R1y*positive(x-x1,1)-R2y*positive(x-x2,1)-R3y*(x-x3)+P*np.sin(theta)*positive(x-x2-xa/2,1)-Ractuator*np.sin(theta)*positive(x-x2+xa/2,1)

def v(x): #+q
    return -1/(6*E*Izz)*(-R1y*positive(x-x1,3)-R2y*positive(x-x2,3)-R3y*positive(x-x3,3)+P*np.sin(theta)*positive(x-x2-xa/2,3)
                         -Ractuator*np.sin(theta)*positive(x-x2+xa/2,3))+C1*x+C2

def w(x):
    return -1/(6*E*Iyy)*(-R1z*positive(x-x1,3)+P*np.cos(theta)*positive(x-x2-xa/2,3)-R2z*positive(x-x2,3)-R3z*positive(x-x3,3)-Ractuator*np.cos(theta)*positive(x-x2+xa/2,3))+C3*x+C4

def Sz(x):
    return -R1z*positive(x-x1,0)-R2z*positive(x-x2,0)-R3z*positive(x-x3,0)-Ractuator*np.cos(theta)*positive(x-x2+xa/2,0)+P*np.cos(theta)*positive(x-x2-xa/2,0)

def Sy(x): #+q
    return -R1y*positive(x-x1,0)-R2y*positive(x-x2,0)-R3y*positive(x-x3,0)-Ractuator*np.sin(theta)*positive(x-x2+xa/2,0)+P*np.sin(theta)*positive(x-x2-xa/2,0)

def T(x): #-q
    return -R1y*positive(x-x1,0)*(-ha/2-z_centroid)-R2y*positive(x-x2,0)*(-ha/2-z_centroid)-R3y*positive(x-x3,0)*(-ha/2-z_centroid)+P*np.sin(theta)*positive(x-x2-xa/2,0)*(-z_centroid)-P*np.cos(theta)*positive(x-x2-xa/2,0)*ha/2-Ractuator*np.sin(theta)*positive(x-x2+xa/2,0)*(-z_centroid)

def twist(x):
    return 1/(G*J)*(-R1y*positive(x-x1,1)*(-ha/2-z_centroid)-R2y*positive(x-x2,1)*(-ha/2-z_centroid)-R3y*positive(x-x3,1)*(-ha/2-z_centroid)+P*np.sin(theta)*positive(x-x2-xa/2,1)*(-z_centroid)-P*np.cos(theta)*positive(x-x2-xa/2,1)*ha/2-Ractuator*np.sin(theta)*positive(x-x2+xa/2,1)*(-z_centroid))+C5