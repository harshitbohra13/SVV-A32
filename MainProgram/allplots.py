from MainInterpolation import *



# plt.figure(3)
fig, axs = plt.subplots(2,2)
x_axis = np.linspace(coor_x[0], coor_x[-1], 100)
v_axis = np.array([])
mz_axis = np.array([])
sy_axis = np.array([])
for i in range(len(x_axis)):
    v_axis = np.append(v_axis, v(x_axis[i]))
    mz_axis = np.append(mz_axis,Moment_z(x_axis[i]))
    sy_axis = np.append(sy_axis,Sy(x_axis[i]))
slopev_axis =(np.diff(v_axis)/np.diff(x_axis))
axs[0,0].plot(x_axis,v_axis)
axs[0,0].set(title='Deflection in y', xlabel='x', ylabel='v(x)')
axs[0,0].grid()
axs[0,1].plot(x_axis[1:],slopev_axis)
axs[0,1].set(title='Slope in y', xlabel='x', ylabel='dv/dx(x)')
axs[0,1].grid()
axs[1,0].plot(x_axis,mz_axis)
axs[1,0].set(title='Bending Moment about z', xlabel='x', ylabel='M_z(x)')
axs[1,0].grid()
axs[1,1].plot(x_axis,sy_axis)
axs[1,1].set(title='Shear Force in y', xlabel='x', ylabel='S_y(x)')
axs[1,1].grid()

# plt.figure(4)
fig, axs = plt.subplots(2,2)
x_axis = np.linspace(coor_x[0], coor_x[-1], 100)
w_axis = np.array([])
my_axis = np.array([])
sz_axis = np.array([])
for i in range(len(x_axis)):
    w_axis = np.append(w_axis, w(x_axis[i]))
    my_axis = np.append(my_axis,Moment_y(x_axis[i]))
    sz_axis = np.append(sz_axis,Sz(x_axis[i]))
slopew_axis =(np.diff(w_axis)/np.diff(x_axis))
axs[0,0].plot(x_axis,w_axis)
axs[0,0].set(title='Deflection in z', xlabel='x', ylabel='w(x)')
axs[0,0].grid()
axs[0,1].plot(x_axis[1:],slopew_axis)
axs[0,1].set(title='Slope in z', xlabel='x', ylabel='dw/dx(x)')
axs[0,1].grid()
axs[1,0].plot(x_axis,my_axis)
axs[1,0].set(title='Bending Moment about y', xlabel='x', ylabel='M_y(x)')
axs[1,0].grid()
axs[1,1].plot(x_axis,sz_axis)
axs[1,1].set(title='Shear Force in z', xlabel='x', ylabel='S_z(x)')
axs[1,1].grid()

# plt.figure(5)
fig, axs = plt.subplots(1,2)
x_axis = np.linspace(coor_x[0], coor_x[-1], 100)
torque_axis = np.array([])
twist_axis = np.array([])
for i in range(len(x_axis)):
    torque_axis = np.append(torque_axis, T(x_axis[i]))
    twist_axis = np.append(twist_axis,Twist(x_axis[i]))
axs[0].plot(x_axis,torque_axis)
axs[0].set(title='Torque', xlabel='x', ylabel='T(x)')
axs[0].grid()
axs[1].plot(x_axis,twist_axis)
axs[1].set(title='Twist', xlabel='x', ylabel='$theta$(x)')
axs[1].grid()


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
    sigma_vm =  np.square(0.5 * ((sigma_xx_total**2) + ((-sigma_xx_total)**2)) + 3*(tau_yz ** 2))

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
