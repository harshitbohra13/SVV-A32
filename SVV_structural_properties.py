import numpy as np
from math import *
from matplotlib import pyplot as plt
from matplotlib.patches import Wedge

#================= Data Fokker 100 =======================

c_a = 0.505     # chord length aileron [m]
l_a = 1.611     # span of the aileron [m]
h_a = 16.1e-2   # aileron height [m] 

t_sk = 1.1e-3   # skin thickness [m]
t_sp = 2.4e-3   # spar thickness [m]

n_st = 11       # number of stiffeners [-]
t_st = 1.2e-3   # thickness of stiffener [m]
h_st = 1.3e-2   # height of stiffener [m]
w_st = 1.7e-2   # width of stiffener [m]
x = 1000



#============= Calculate stiffener locations ===========================

circ = pi*0.5*h_a + 2*sqrt((c_a - 0.5*h_a)**2. + (0.5*h_a)**2.) # circumference of aileron
s_st = circ/n_st # stiffener spacing [m]

# set up arrays for x and y locations of stiffener areas
B_z = np.zeros(11)
B_y = np.zeros(11)

angle = s_st / (0.5*h_a) # polar coordinate of the stiffener in the circular part
B_y[1] = (0.5*h_a)*sin(angle) 
B_z[1] = 0.5*h_a - (0.5*h_a)*cos(angle)
B_y[10] = -(0.5*h_a)*sin(angle)
B_z[10] = 0.5*h_a - (0.5*h_a)*cos(angle)
# calculate 
space = 2*s_st - (pi*0.5*h_a)/2  
slope = - ((h_a/2 - 0)/(c_a-h_a/2 - 0))
theta = atan((h_a/2)/(c_a - h_a/2))

B_y[2]=(h_a/2)-(space*sin(theta))
B_z[2]=(h_a/2)+(space*cos(theta))
B_y[9]=-((h_a/2)-(space*sin(theta)))
B_z[9]=(h_a/2)+(space*cos(theta))
B_y[3]=B_y[2]-(s_st*sin(theta))
B_z[3]=B_z[2]+(s_st*cos(theta))
B_y[8]=-(B_y[2]-(s_st*sin(theta)))
B_z[8]=B_z[2]+(s_st*cos(theta))
B_y[4]=B_y[3]-(s_st*sin(theta))
B_z[4]=B_z[3]+(s_st*cos(theta))
B_y[7]=-(B_y[3]-(s_st*sin(theta)))
B_z[7]=B_z[3]+(s_st*cos(theta))
B_y[5]=B_y[4]-(s_st*sin(theta))
B_z[5]=B_z[4]+(s_st*cos(theta))
B_y[6]=-(B_y[4]-(s_st*sin(theta)))
B_z[6]=B_z[4]+(s_st*cos(theta))
B_z = -1*B_z

#============= Calculate centroid location ======================

# z component centroid location individual components [m]
z_spar = -1*(0.5*h_a)
z_circ = -1*(0.5*h_a - h_a/pi)
z_tri = -1*(0.5*h_a + 0.5* (c_a-0.5*h_a))
z_stiff = (sum(B_z)/len(B_z))
# areas individual components [m^2]
A_spar = t_sp*h_a
A_circ = pi*0.5*h_a*t_sk
A_tri = 2*t_sk*sqrt((c_a-0.5*h_a)**2 + (0.5*h_a)**2)
A_stiff = n_st*t_st*(h_st+w_st-t_st) # total area of stiffeners

# Centroid location aileron 
z_cent = (z_spar*A_spar + z_circ*A_circ + z_tri*A_tri + z_stiff*A_stiff ) / (A_spar + A_circ + A_tri + A_stiff)
y_cent = 0.
print()
print("Centroid location (z,y) = (",z_cent,",",y_cent,")")

#=================== Calculate second moments of area =======================

# contribution to I due to circular part
I_z_circ = (1/16)*pi*(h_a**3)*t_sk 
I_y_circ = (1/16)*pi*(h_a**3)*t_sk -A_circ*(h_a/pi)**2+A_circ*(z_circ-z_cent)**2 
#I_y_circ = (1/16)*pi*(h_a**3)*t_sk + A_circ*(z_circ-z_cent)**2 

# Contribution to I due to spar
I_z_spar = (1/12)*t_sp*(h_a**3)

# contribution to I due to circular part
I_z_circ = (1/16)*pi*(h_a**3)*t_sk
I_y_circ = (1/16)*pi*(h_a**3)*t_sk -A_circ*(h_a/pi)**2+ A_circ*(z_circ-z_cent)**2 
#I_y_circ = (1/16)*pi*(h_a**3)*t_sk + A_circ*(z_circ-z_cent)**2 

# Contribution to I due to spar
I_z_spar = (1/12)*t_sp*(h_a**3)

I_y_spar = 0 + A_spar*(z_spar-z_cent)**2 # only steiner term because thin-walled assumption

# Contribution to I due to triangular part
len_sk = sqrt((c_a-0.5*h_a)**2 + (0.5*h_a)**2)
I_z_tri = 2 * ((1/12)*t_sk*((len_sk)**3)*(sin(theta))**2 + ((0.25*h_a)**2)*t_sk*len_sk)
I_y_tri = 2 * ((1/12)*t_sk*((len_sk)**3)*(cos(theta))**2 + ((z_tri-z_cent)**2)*t_sk*len_sk)

# Contribution to I due to stiffeners
I_z_stiff = 0
I_y_stiff = 0
for i in range(len(B_y)): 
    I_z_stiff = I_z_stiff + t_st*(h_st+w_st-t_st)*(
    B_y[i])**2          # add up all individual area moments
    I_y_stiff = I_y_stiff + t_st*(h_st+w_st-t_st)*(B_z[i]-z_cent)**2   # from the stiffeners

# Computation of total area moments about local z and y axes
I_zz = I_z_circ + I_z_spar + I_z_tri + I_z_stiff # [m^4]
I_yy = I_y_circ + I_y_spar + I_y_tri + I_y_stiff # [m^4]
print()
print('I_yy, I_zz = ',I_yy,', ',I_zz)

#============================ Polar moment of Inertia =====================

T = 1
r = h_a/2
A1 = 0.5*pi*r**2
A2 = (c_a-h_a/2)*h_a/2

A = np.array([[(h_a/(2*A1*t_sp))+(pi*r)/(2*A1*t_sk),-h_a/(2*A1*t_sp),-1],
              [-h_a/(2*A2*t_sp),(h_a/(2*A2*t_sp))+((2*len_sk)/(2*A2*t_sk)),-1],
              [2*A1,2*A2,0]])

B = np.array([0,0,T])

x = np.linalg.solve(A,B)
print("Matrix solution =",x)

Gdsigmadz = x[2]
J = T/(Gdsigmadz)
print("J =",J)

#============================ PLOTTING ====================================

fig, ax = plt.subplots()

ycirc = np.linspace(0,0.5*h_a,1000)
zcirc = np.sqrt(-(ycirc-0.5*h_a)**2+(0.5*h_a)**2)

zspar = np.full((1000, 1), -0.5*h_a)
yspar = np.linspace(-0.5*h_a,0.5*h_a,1000)

zskin = np.linspace(-0.5*h_a,-c_a,1000)
ytskin = -slope*(zskin+0.5*h_a)+0.5*h_a
ybskin = slope*(zskin+0.5*h_a)-0.5*h_a

ax.plot(-ycirc,zcirc,'k')  # plot upper part cemicircle
ax.plot(-ycirc,-zcirc,'k') # plot bottom part semicircle
ax.plot(zspar,yspar,'k')   # plot spar
ax.plot(zskin,ytskin,'k')  # plot 
ax.plot(zskin,ybskin,'k')
ax.scatter(B_z,B_y,label='Stiffeners')
ax.plot(z_cent,y_cent,'x', label='Centroid')
ax.arrow(0,0,0.035,0)
ax.arrow(0,0,0,0.035)
ax.text(0.035,0.01,"z'")
ax.text(0.01,0.045,"y'")
#matplotlib.pyplot.arrow(x, y, dx, dy, **kwargs)[source]¶

plt.title('Cross-section of aileron')
plt.ylabel('y [m]')
plt.xlabel('z [m]')
ax.legend()
ax.set_xlim(0.1*c_a, -1.05*c_a)
ax.set_ylim(-h_a, h_a)
# plt.show()

#========================== Verification ==================================

