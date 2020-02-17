import numpy as np
from math import *


#___________ Data Fokker 100 _________________________

c_a = 0.505     # chord length aileron [m]
l_a = 1.611     # span of the aileron [m]
h_a = 16.1e-2   # aileron height [m] 

t_sk = 1.1e-3   # skin thickness [m]
t_sp = 2.4e-3   # spar thickness [m]

n_st = 11.      # number of stiffeners [-]
t_st = 1.2e-3   # thickness of stiffener [m]
h_st = 1.3e-2   # height of stiffener [m]
w_st = 1.7e-2   # width of stiffener [m]


#________ Calculate stiffener locations____________

circ = pi*0.5*h_a + 2*sqrt((c_a - 0.5*h_a)**2. + (0.5*h_a)**2.) # circumference of aileron
s_st = circ/n_st # stiffener spacing [m]

# set up arrays for x and y locations of boom areas
B_z = np.zeros(11)
B_y = np.zeros(11)

theta = s_st / (0.5*h_a)
B_y[1] = (0.5*h_a)*sin(theta)
B_z[1] = 0.5*h_a - (0.5*h_a)*cos(theta)
B_y[10] = -(0.5*h_a)*sin(theta)
B_z[10] = 0.5*h_a - (0.5*h_a)*cos(theta)

space3b = 2*s_st - (pi*0.5*h_a)/2  
slope = - ((h_a/2 - 0)/(c_a-h_a/2 - 0))
theta_2 = atan((h_a/2)/(c_a - h_a/2))

B_y[2]=(h_a/2)-(space3b*sin(theta_2))
B_z[2]=(h_a/2)+(space3b*cos(theta_2))
B_y[9]=-((h_a/2)-(space3b*sin(theta_2)))
B_z[9]=(h_a/2)+(space3b*cos(theta_2))
B_y[3]=B_y[2]-(s_st*sin(theta_2))
B_z[3]=B_z[2]+(s_st*cos(theta_2))
B_y[8]=-(B_y[2]-(s_st*sin(theta_2)))
B_z[8]=B_z[2]+(s_st*cos(theta_2))
B_y[4]=B_y[3]-(s_st*sin(theta_2))
B_z[4]=B_z[3]+(s_st*cos(theta_2))
B_y[7]=-(B_y[3]-(s_st*sin(theta_2)))
B_z[7]=B_z[3]+(s_st*cos(theta_2))
B_y[5]=B_y[4]-(s_st*sin(theta_2))
B_z[5]=B_z[4]+(s_st*cos(theta_2))
B_y[6]=-(B_y[4]-(s_st*sin(theta_2)))
B_z[6]=B_z[4]+(s_st*cos(theta_2))

#______________ Calculate centroid location ____________

z_spar = 0.5*h_a
z_circ = 2*(0.5*h_a)/pi
z_tri = 0.5*h_a + 0.5* (c_a-0.5*h_a)
z_stiff = sum(B_z)/len(B_z)


