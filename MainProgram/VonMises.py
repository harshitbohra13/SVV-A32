from MainInterpolation import *
from plotshear import *

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


location = 1 #[m] #Location of the cross section
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
y5= -y2
z5=z2

#First Section
z1=np.linspace(0,h,step)
y1=np.sqrt(-(z1-h)**2+h**2)

#Sixth Section
z6=z1
y6=-y1


##############DIRECT-STRESSES###############
def direct_stress(My,Mz,z,y):
    sigma_xx_z = My * (z - z_hat) / Iyy
    sigma_xx_y = Mz * y / Izz
    sigma_xx_total = sigma_xx_y + sigma_xx_z
    return sigma_xx_total

sigma_1 = []
sigma_2 = []
sigma_3 = []
sigma_4 = []
sigma_5 = []
sigma_6 = []

for n in range(0,step):
    sigma_3.append(direct_stress(Moment_y(location),Moment_z(location),-(z3[n]),y3[n]))
    sigma_4.append(direct_stress(Moment_y(location), Moment_z(location), -(z4[n]), y4[n]))
    sigma_2.append(direct_stress(Moment_y(location), Moment_z(location), -(z2[n]), y2[n]))
    sigma_5.append(direct_stress(Moment_y(location), Moment_z(location), -(z5[n]), y5[n]))
    sigma_1.append(direct_stress(Moment_y(location), Moment_z(location), -(z1[n]), y1[n]))
    sigma_6.append(direct_stress(Moment_y(location), Moment_z(location), -(z6[n]), y6[n]))

plt.figure(8)
#Plot all the sections
marker_size = 15
plt.scatter(-z3,y3,marker_size,sigma_3)
plt.scatter(-z4,y4,marker_size,sigma_4)
plt.scatter(-z2,y2,marker_size,sigma_2)
plt.scatter(-z5,y5,marker_size,sigma_5)
plt.scatter(-z1,y1,marker_size,sigma_1)
plt.scatter(-z6,y6,marker_size,sigma_6)
cbar = plt.colorbar()
plt.title('Direct Stress Distribution')
plt.ylabel('y[m]')
plt.xlabel('-z[m]')
cbar.set_label('$sigma_{xx}$[N/m^2]')
plt.grid()




def VonMises(My, Mz,tau_yz,z, y):

    sigma_xx_z = My * (z - z_hat) / Iyy
    sigma_xx_y = Mz * y / Izz

    sigma_xx_total = sigma_xx_y + sigma_xx_z
    sigma_vm =  np.sqrt(0.5 * ((sigma_xx_total**2) + ((-sigma_xx_total)**2)) + 3*(tau_yz ** 2))

    return sigma_vm

vonMises_1 = []
vonMises_2 = []
vonMises_3 = []
vonMises_4 = []
vonMises_5 = []
vonMises_6 = []

for n in range(0,step):
    vonMises_3.append(VonMises(Moment_y(location),Moment_z(location),q3[n]/t_sk,-(z3[n]),y3[n]))
    vonMises_4.append(VonMises(Moment_y(location), Moment_z(location),q4[n]/t_sk, -(z4[n]), y4[n]))
    vonMises_2.append(VonMises(Moment_y(location), Moment_z(location),q2[n]/t_sp, -(z2[n]), y2[n]))
    vonMises_5.append(VonMises(Moment_y(location), Moment_z(location),q5[n]/t_sp, -(z5[n]), y5[n]))
    vonMises_1.append(VonMises(Moment_y(location), Moment_z(location),q1[n]/t_sk, -(z1[n]), y1[n]))
    vonMises_6.append(VonMises(Moment_y(location), Moment_z(location),q6[n]/t_sk, -(z6[n]), y6[n]))

plt.figure(9)
#Plot all the sections
marker_size = 15
plt.scatter(-z3,y3,marker_size,vonMises_3)
plt.scatter(-z4,y4,marker_size,vonMises_4)
plt.scatter(-z2,y2,marker_size,vonMises_2)
plt.scatter(-z5,y5,marker_size,vonMises_5)
plt.scatter(-z1,y1,marker_size,vonMises_1)
plt.scatter(-z6,y6,marker_size,vonMises_6)
cbar = plt.colorbar()
plt.title('Von Mises Stress Distribution')
plt.ylabel('y[m]')
plt.xlabel('-z[m]')
cbar.set_label('$sigma_{vm}$[N/m^2]')
plt.grid()

print("Runtime: %f seconds" % (time.time()-start_time))

plt.show()
