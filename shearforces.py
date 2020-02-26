import numpy as np
import functionshearforce as sf
import SVV_structural_properties as prop 
from shearcentercalc import get_sc, h, lsk

    
Sy = 1
Sz = 1
 
z_sc, y_sc = get_sc()

qy1 = sf.get_qy(1, Sy, prop.I_zz, sf.sintegrate, prop.t_sk, h**2, 1000, 0, np.pi/2)
qy2 = sf.get_qy(2, Sy, prop.I_zz, sf.get_integral, prop.t_sp, h, 1000, 0, np.pi/2)
qy3 = sf.get_qy(3, Sy, prop.I_zz, sf.sintegrate, prop.t_sk, h**2, 1000, 0, np.pi/2)
qy4 = sf.get_qy(4, Sy, prop.I_zz, sf.sintegrate, prop.t_sk, h**2, 1000, 0, np.pi/2)
qy5 = sf.get_qy(5, Sy, prop.I_zz, sf.sintegrate, prop.t_sk, h**2, 1000, 0, np.pi/2)
qy5 = sf.get_qy(6, Sy, prop.I_zz, sf.sintegrate, prop.t_sk, h**2, 1000, 0, np.pi/2)    
    