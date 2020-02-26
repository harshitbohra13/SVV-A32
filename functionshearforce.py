import numpy as np
import SVV_structural_properties as prop
import MainNewMatrixCalculation as Mat 
from shearcentercalc import h, get_sc


z_sc, y_sc = get_sc()

class forces:
    Sy = np.zeros(41)
    Sz = np.zeros(41)
    Ta = np.zeros(41)

    def get_allvals(self): 
         j = 0 
         vals = np.linspace(0, Mat.la, num = 41)
         for i in vals:
            self.Sy[j] =  Mat.Sy(i)
            self.Sz[j] =  Mat.Sz(i)
            self.Ta[j] =  self.Sy[j]*(h + z_sc)
            self.Ta[j] += Mat.T(i)
            j+=1
      
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

def get_q0(T):
    x = [[(2*0.5*(prop.c_a - prop.h_a/2)*prop.h_a), (2*prop.pi*0.5*(prop.h_a/2)**2)], 
        [-1, 0.783]]
    b = [T, 0]    
    return (np.linalg.solve(x,b))

def get_qboom(section, dir, s):
    omega = np.tanh(h/(prop.c_a - h))
    if(section == 1):
        if((h-h*np.cos(s)) < prop.B_z[1]):
            if(dir == "y"):
                return(prop.B_y[1] * prop.A_stiff)
            else:
                return(0.5* prop.B_z[0] * prop.A_stiff + prop.B_z[1] * prop.A_stiff)
        elif(dir == "z"):
            return(0.5* prop.B_z[0] * prop.A_stiff)

    if (section == 6):
        if((h-h*np.cos(s)) < prop.B_z[1]):
            if(dir == "y"):
                return(prop.B_y[10] * prop.A_stiff)
            else:
                return(0.5* prop.B_z[0] * prop.A_stiff + prop.B_z[10] * prop.A_stiff)
        elif(dir == "z"):
            return(0.5* prop.B_z[0] * prop.A_stiff)
    
    if (section == 3):
        qboom = 0
        if(dir == "y"):
            for i in range(2, 6, 1):
                if(prop.B_z[i] > (prop.c_a - h - s*np.cos(omega))):
                    qboom += prop.B_y[i]*prop.A_stiff
            return(qboom)
        else:
            for i in range(2, 6, 1):
                if(prop.B_z[i] > (prop.c_a - h - s*np.cos(omega))):
                    qboom += prop.B_z[i]*prop.A_stiff
            return(qboom)
  
    if (section == 4):
        qboom = 0
        if(dir == "y"):
            for i in range(6, 10, 1):
                if(prop.B_z[i] > (prop.c_a - h - s*np.cos(omega))):
                    qboom += prop.B_y[i]*prop.A_stiff
            return(qboom)
        else:
            for i in range(6, 10, 1):
                if(prop.B_z[i] > (prop.c_a - h - s*np.cos(omega))):
                    qboom += prop.B_z[i]*prop.A_stiff
            return(qboom)
    else:
        return(0)


def get_qy(section, Sy, Izz, f, N, a, b, qs0 = 0):
    q = (get_integral(f,N,a,b) + get_qboom(section, "y", b))
    return((-Sy/Izz  *  q) + qs0)

def get_qz(section, Sz, Iyy, f, N, a, b, qs0 = 0):
    q = (get_integral(f,N,a,b) + get_qboom(section, "z", b))
    return((-Sz/Iyy *  q) + qs0)

