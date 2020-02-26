import numpy as np
import SVV_structural_properties as prop

def forcetoshearcenter(F, z, y, theta=0):
    M = 0
    Fz = F*np.sin(theta)
    Fy = F*np.cos(theta)

    M += F*z
    M += F*y

    return(Fz, Fy, M)

def sintegrate(N, a, b):
    f = lambda x: np.sin(x)
    
    num0 = 0
    num1 = 0

    for i in range(1, N + 1):
        num0 += f(a + (i-(1/2))*((b-a)/N))
    num1 = ((b-a)/N)*num0

    return num1

def costegrate(N, a, b):
    f = lambda x: np.cos(x)
    
    num0 = 0
    num1 = 0

    for i in range(1, N + 1):
        num0 += f(a + (i-(1/2))*((b-a)/N))
    num1 = ((b-a)/N)*num0

    return num1
    
def get_integral(func, N, a, b):
    num0 = 0
    num1 = 0
    for i in range(1, N + 1):
        num0 += func(a + (i-(1/2))*((a - b)/N))
    num1 = ((a - b)/N)*num0
    
    return num1

def get_q0(Fx, Fy, T):
    q01 = T/(2*prop.pi*0.5*(prop.h_a/2)**2)
    q02 = T/(2*0.5*(prop.c_a - prop.h_a/2)*prop.h_a)
    q = [q01, q02]
    return (q)

def get_qboom(section, dir):
    if(section == 1):
        if(dir == "y"):
             return(prop.B_y[1] * prop.A_stiff)
        else:
             return(0.5* prop.B_z[0] * prop.A_stiff + prop.B_z[1] * prop.A_stiff)

    if (section == 6):
        if(dir == "y"):
             return(prop.B_y[10] * prop.A_stiff)
        else:
             return(0.5* prop.B_z[0] * prop.A_stiff + prop.B_z[10] * prop.A_stiff)

    if (section == 3):
        qboom = 0
        if(dir == "y"):
            for i in range(2, 6, 1):
                qboom += prop.B_y[i]*prop.A_stiff
            return(qboom)
        else:
            for i in range(2, 6, 1):
                qboom += prop.B_z[i]*prop.A_stiff
            return(qboom)
  
    if (section == 4):
        qboom = 0
        if(dir == "y"):
            for i in range(6, 10, 1):
                qboom += prop.B_y[i]*prop.A_stiff
            return(qboom)
        else:
            for i in range(6, 10, 1):
                qboom += prop.B_z[i]*prop.A_stiff
            return(qboom)
    else:
        return(0)


def get_qy(section, Sy, Izz, f, N, a, b, qs0 = 0):
    q = (get_integral(f,N,a,b) + get_qboom(section, "y"))
    return((-Sy/Izz  *  q) + qs0)

def get_qz(section, Sz, Iyy, f, N, a, b, qs0 = 0):
    q = (get_integral(f,N,a,b) + get_qboom(section, "z"))
    return((-Sz/Iyy *  q) + qs0)
