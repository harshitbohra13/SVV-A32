from MainInterpolation import Sy,coor_x,Moment_z,v,Sz,Moment_y,w,Twist,T
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import Verification.main as verification


##########Comparing Plots##############
x_axis = np.linspace(coor_x[0], coor_x[-1], 100)
v_axis = np.array([])
mz_axis = np.array([])
sy_axis = np.array([])
w_axis = np.array([])
my_axis = np.array([])
sz_axis = np.array([])
t_axis = np.array([])
twist_axis = np.array([])
for i in range(len(x_axis)):
    v_axis = np.append(v_axis, v(x_axis[i]))
    mz_axis = np.append(mz_axis,Moment_z(x_axis[i]))
    sy_axis = np.append(sy_axis,Sy(x_axis[i]))
    w_axis = np.append(w_axis, w(x_axis[i]))
    my_axis = np.append(my_axis,Moment_y(x_axis[i]))
    sz_axis = np.append(sz_axis,Sz(x_axis[i]))
    twist_axis = np.append(twist_axis,Twist(x_axis[i]))
    t_axis = np.append(t_axis,T(x_axis[i]))
slopev_axis =(np.diff(v_axis)/np.diff(x_axis))
slopew_axis =(np.diff(w_axis)/np.diff(x_axis))

#plt.subplot(2,2,1)
#plt.plot(x_axis,v_axis, 'r')
#
#plt.subplot(2,2,2)
#plt.plot(x_axis[1:],slopev_axis, 'r')
#
#plt.subplot(2,2,2)
#plt.plot(x_axis,mz_axis, 'r')
#
#plt.subplot(2,2,4)
#plt.plot(x_axis,sy_axis, 'r')

plt.figure(8)
x_axis = np.linspace(coor_x[0], coor_x[-1], 100)
Mzverification = []
Myverification = []
Vverification = []
Wverification = []
Tverification = []
Twistverification = []
for i in range(len(x_axis)):
    Mzverification.append(mz_axis[i]-verification.aileron.Mz(x_axis[i]))
    Myverification.append(my_axis[i]-verification.aileron.My(x_axis[i]))
    Vverification.append(v_axis[i]-verification.aileron.eval(x_axis[i])[0])
    Wverification.append(w_axis[i]-verification.aileron.eval(x_axis[i])[1])
    Tverification.append(t_axis[i]-verification.aileron.T(x_axis[i]))
    Twistverification.append(twist_axis[i]-verification.aileron.eval(x_axis[i])[2])
    
plt.subplot(2,2,1)
plt.plot(x_axis,Mzverification,'k',label='Mz error')
plt.plot(x_axis,Myverification,'r',label='My error')
plt.xlabel('x')
plt.ylabel('Bending Moment')
plt.legend()
plt.grid()

plt.subplot(2,2,2)
plt.plot(x_axis,Vverification,'k',label='v error')
plt.plot(x_axis,Wverification,'r',label='w error')
plt.xlabel('x')
plt.ylabel('Displacement')
plt.legend()
plt.grid()

plt.subplot(2,2,3)
plt.plot(x_axis,Tverification,'r',label='T error')
plt.xlabel('x')
plt.ylabel('Torque')
plt.legend()
plt.grid()

plt.subplot(2,2,4)
plt.plot(x_axis,Twistverification,'r',label='Twist error')
plt.xlabel('x')
plt.ylabel('Theta')
plt.legend()
plt.grid()

plt.show()
    