import numpy as np
import functionshearforce as sf
import SVV_structural_properties as prop 
from shearcentercalc import get_sc, h, lsk

    
Sy = 1
Sz = 1

z_c = prop.z_cent - h #z-location from the spar

z_sc, y_sc = get_sc()

fh1 = lambda x: -h**2*t +h**2*t* np.cos(x)+z_c*h*t
fh2 = lambda x: -h*t-t*z_c
fh3 = lambda x: -h*t -z_c*t-((prop.c_a-h)*t*s/lsk)
fh4 = lambda x: -prop.c_a*t-z_c*t + (prop.c_a-h)*t*s/lsk
fh5 = lambda x: t*h + t*z_c
fh6 = lambda x: -h**2*t + h**2*np.cos(x)*t + z_c*t*h

qy1 = sf.get_qy(1, Sy, prop.I_zz, sf.sintegrate, prop.t_sk, h**2, 1000, 0, np.pi/2) #vertical shear flow for section 1
qz1 = sf.get_qz(1, Sz, prop.I_yy, fh1, 1000, 0 ,np.pi/2) #horizontal shear flow for section 1

qy2 = sf.get_qy(2, Sy, prop.I_zz, sf.get_integral, prop.t_sp, h, 1000, 0, h) #vertical shear flow for section 2
qz2 = sf.get_qz(2, Sz, prop.I_yy, fh2, 1000, 0, h) #horizontal shear flow for section 2

qy3 = sf.get_qy(3, Sy, prop.I_zz, sf.sintegrate, prop.t_sk, h**2, 1000, 0, lsk) #vertical shear flow for section 3
qz3 = sf.get_qz(3,Sz, prop.I_yy, fh3, 1000, 0, lsk, (qy1+qz1+qy2+qz2)) #horizontal shear flow for section 3

qy4 = sf.get_qy(4, Sy, prop.I_zz, sf.sintegrate, prop.t_sk, h**2, 1000, 0, lsk) #vertical shear flow for section 4
qz4 = sf.get_qz(4,Sz, prop.I_yy, fh4, 1000, 0, lsk, (qy3 + qz3)) #horizontal shear flow for section 4

qy5 = sf.get_qy(5, Sy, prop.I_zz, sf.sintegrate, prop.t_sk, h**2, 1000, 0, h) #vertical shear flow for section 5
qz5 = sf.get_qz(5, Sz, prop.I_yy, fh5, 1000, 0, h) #horizontal shear flow for section 5

qy6 = sf.get_qy(6, Sy, prop.I_zz, sf.sintegrate, prop.t_sk, h**2, 1000, 0, np.pi/-2) #vertical shear flow for section 6
qz6 = sf.get.qz(6,Sz, prop.I_yy, fh6, 1000, np.pi/-2 ,0, (qy4+qz4+qz5,qy5)) #horizontal shear flow for section 6

print("saga")
print(qz6)

