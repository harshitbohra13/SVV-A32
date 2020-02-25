import numpy as np
import SVV_structural_properties as prop

def forcetoshearcenter(F, z, y, theta=0):
    M = 0
    Fz = F*np.sin(theta)
    Fy = F*np.cos(theta)

    M += F*z
    M += F*y

    return(Fz, Fy, M)

def getq(Fx, Fy, T):
    q01 = T/(2*prop.pi*0.5*(prop.h_a/2)**2)
    q02 = T/(2*0.5*(prop.c_a - prop.h_a/2)*prop.h_a)
    q = [q01, q02]
    return (q)