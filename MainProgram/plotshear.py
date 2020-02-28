from Functions_Interpolation import find_interpolant, compute_value_interpolation,integrate_spline,doubleintegrate_spline, \
    quadrupleintegrate_spline,positive
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import time
from shearforces import shearstress
import SVV_structural_properties as prop
from scipy.integrate import quad
from shearcentercalc import get_sc
import Data as data
from MainInterpolation import Sy, Sz, T

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


####################AIRFOIL - STRESSES##########################


location = 0.4 #[m] #Location of the cross section
h = ha/2 #cm

#Step for each point
step=1000

#Third Section
z3=np.linspace(h,Ca,step)
y3=Ca*h/(Ca-h)-z3*h/(Ca-h)

#Fourth Section
z4=z3
y4=-y3

#Second Section
y2=np.linspace(0,h,step)
z2=np.ones(len(y2))*h

#Fifth Section
y5=-y2
z5=z2

#First Section
z1=np.linspace(0,h,step)
y1=np.sqrt(-(z1-h)**2+h**2)

#Sixth Section
z6=z1
y6=-y1

#############SHEAR-FLOWS######################
q1,q2,q3,q4,q5,q6 = shearstress(location,Sy(location),Sz(location),T(location))
q4 = q4[::-1]
q6 = q6[::-1]


plt.figure(7)
#Plot all the sections
marker_size = 15
plt.scatter(-z3,y3,marker_size,q3)
plt.scatter(-z4,y4,marker_size,q4)
plt.scatter(-z2,y2,marker_size,q2)
plt.scatter(-z5,y5,marker_size,q5)
plt.scatter(-z1,y1,marker_size,q1)
plt.scatter(-z6,y6,marker_size,q6)
cbar = plt.colorbar()
plt.title('Shear Flow Distribution')
plt.ylabel('y[m]')
plt.xlabel('-z[m]')
cbar.set_label('q[N/m]')
plt.grid()

plt.show()