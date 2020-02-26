import numpy as np
import functionshearforce as sf
import SVV_structural_properties as prop 
from shearcentercalc import get_sc, h, lsk


Sy = 1
Sz = 1

z_c = prop.z_cent - h #z-location from the spar

# z_sc, y_sc = get_sc()

fv1 = lambda y: prop.t_sk*h**2*np.sin(y)
fv2 = lambda y: prop.t_sp*y
fv3 = lambda s: prop.t_sk*h - (prop.t_sk*h*s)/lsk
fv4 = lambda s: prop.t_sk*-h*s/lsk
fv5 = lambda y: prop.t_sp*y
fv6 = lambda y: prop.t_sk*h**2*np.sin(y)

fh1 = lambda z: -h**2*prop.t_sk +h**2*prop.t_sk* np.cos(z)+z_c*h*prop.t_sk
fh2 = lambda z: -h*prop.t_sp-prop.t_sp*z_c
fh3 = lambda s: -h*prop.t_sk -z_c*prop.t_sk-((prop.c_a-h)*prop.t_sk*s)/lsk
fh4 = lambda s: -prop.c_a*prop.t_sk-z_c*prop.t_sk + ((prop.c_a-h)*prop.t_sk*s)/lsk
fh5 = lambda z: -prop.t_sp*h - prop.t_sp*z_c
fh6 = lambda z: -h**2*prop.t_sk+ h**2*np.cos(z)*prop.t_sk - z_c*prop.t_sk*h

qy1 = sf.get_qy(1, Sy, prop.I_zz, fv1, 1000, 0, np.pi/2) #vertical shear flow for section 1
qz1 = sf.get_qz(1, Sz, prop.I_yy, fh1, 1000, 0, np.pi/2) #horizontal shear flow for section 1

qy2 = sf.get_qy(2, Sy, prop.I_zz, fv2, 1000, 0, h) #vertical shear flow for section 2
qz2 = sf.get_qz(2, Sz, prop.I_yy, fh2, 1000, 0, h) #horizontal shear flow for section 2

qy3 = sf.get_qy(3, Sy, prop.I_zz, fv3, 1000, 0, lsk) #vertical shear flow for section 3
qz3 = sf.get_qz(3, Sz, prop.I_yy, fh3, 1000, 0, lsk, (qy1+qz1+qy2+qz2)) #horizontal shear flow for section 3

qy4 = sf.get_qy(4, Sy, prop.I_zz, fv4, 1000, 0, lsk) #vertical shear flow for section 4
qz4 = sf.get_qz(4, Sz, prop.I_yy, fh4, 1000, 0, lsk, (qy3 + qz3)) #horizontal shear flow for section 4

qy5 = sf.get_qy(5, Sy, prop.I_zz, fv5, 1000, 0, -h) #vertical shear flow for section 5
qz5 = sf.get_qz(5, Sz, prop.I_yy, fh5, 1000, 0, -h) #horizontal shear flow for section 5

qy6 = sf.get_qy(6, Sy, prop.I_zz, fv6, 1000, np.pi/-2, 0) #vertical shear flow for section 6
qz6 = sf.get_qz(6, Sz, prop.I_yy, fh6, 1000, np.pi/-2 ,0, (qy4+qz4+qz5+qy5)) #horizontal shear flow for section 6


theta = np.linspace(0,np.pi/2, num = 1000)
y_points = np.linspace(0,prop.h_a/2, num = 1000) #y-points for section 2
s_points = np.linspace(0, lsk, num = 1000) #s points for section 3 and 4
y_points2 = np.linspace(0,prop.h_a/-2, num=1000) #y-points for section 5
theta2 = np.linspace(np.pi/-2,0,num = 1000) # theta for section 6

q1_ar = [] #total shear, horizontal + vertical for section 1
q2_ar = [] #total shear, horizontal + vertical for section 2
q3_ar = [] #total shear, horizontal + vertical for section 3
q4_ar = [] #total shear, horizontal + vertical for section 4
q5_ar = [] #total shear, horizontal + vertical for section 5
q6_ar = [] #total shear, horizontal + vertical for section 6

for i in range(1000):
    q1_ar.append(sf.get_qy(1, Sy, prop.I_zz, fv1, 100, 0, theta[i]) + sf.get_qz(1, Sz, prop.I_yy, fh1, 100, 0, theta[i]))
    q2_ar.append(sf.get_qy(2, Sy, prop.I_zz, fv2, 100, 0, y_points[i]) + sf.get_qz(2, Sz, prop.I_yy, fh2, 100, 0, y_points[i]))
    q3_ar.append(sf.get_qy(3, Sy, prop.I_zz, fv3, 100, 0, s_points[i]) + sf.get_qz(3, Sz, prop.I_yy, fh3, 100, 0, s_points[i], (qy1+qz1+qy2+qz2)))
    q4_ar.append(sf.get_qy(4, Sy, prop.I_zz, fv4, 100, 0, s_points[i]) + sf.get_qz(4, Sz, prop.I_yy, fh4, 100, 0, s_points[i], (qy3 + qz3)))
    q5_ar.append(sf.get_qy(5, Sy, prop.I_zz, fv5, 100, 0, y_points2[i]) + sf.get_qz(5, Sz, prop.I_yy, fh5, 100, 0, y_points2[i]))
    q6_ar.append((sf.get_qy(6, Sy, prop.I_zz, fv6, 100, np.pi/-2, theta2[i]))+sf.get_qz(6, Sz, prop.I_yy, fh6, 100, np.pi/-2 ,theta2[i], (qy4+qz4+qz5+qy5)))

print(q5_ar[999], qy5+qz5)

