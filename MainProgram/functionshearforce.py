import numpy as np
import MainNewMatrixCalculation as Mat, SVV_structural_properties as prop
from shearcentercalc import h, get_sc

z_c = prop.z_cent 

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
    
def get_integral(f, n, a, b):
    h = float(b-a)/n
    result = 0
    for i in range(n):
        result += f((a + h/2.0) + i*h)
    result *= h
    return result

def get_q0(T):
    x = np.array([[(2*0.5*(prop.Ca - prop.ha/2)*prop.ha), (2*prop.pi*0.5*(prop.ha/2)**2)],
        [-1, 1.0125]])
    b = np.array([T, 0])
    return (np.linalg.solve(x,b))

def get_qboom(section, dir, s):
    omega = np.tanh(h/(prop.Ca - h))

    if(section == 1):
        if(-(h-h*np.cos(s)) < prop.B_z[1]):
            if(dir == "y"):
                return(prop.B_y[1] * prop.A_stiff/prop.n_st)
            else:
                return(0.5*(prop.B_z[0] - z_c )* prop.A_stiff/prop.n_st + (prop.B_z[1]- z_c) * prop.A_stiff/prop.n_st)
        elif(dir == "z"):
            return(0.5*(prop.B_z[0]-z_sc) * prop.A_stiff)

    if (section == 6):
        if((h-h*np.cos(s)) < prop.B_z[1]):
            if(dir == "y"):
                return(prop.B_y[prop.n_st -1] * prop.A_stiff/prop.n_st)
            else:
                return(0.5* (prop.B_z[0] - z_c) * prop.A_stiff -(prop.B_z[prop.n_st-1] - z_c )  * prop.A_stiff/prop.n_st)
        elif(dir == "z"):
            return(0.5* -(prop.B_z[prop.n_st-1] - z_c ) * prop.A_stiff)
    
    if (section == 3):
        qboom = 0
        if(dir == "y"):
            for i in range(2, int((prop.n_st + 1)/2), 1):
                if(prop.B_z[i] > (prop.Ca - h - s*np.cos(omega))):
                    qboom += prop.B_y[i]*prop.A_stiff/prop.n_st
            return(qboom)
        else:
            for i in range(2, int((prop.n_st + 1)/2), 1):
                if(prop.B_z[i] > (prop.Ca - h - s*np.cos(omega))):
                    qboom += (prop.B_z[i] - z_c)*prop.A_stiff/prop.n_st
            return(qboom)
  
    if (section == 4):
        qboom = 0
        if(dir == "y"):
            for i in range(int((prop.n_st + 1)/2), prop.n_st-1, 1):
                if(prop.B_z[i] > (prop.Ca - h - s*np.cos(omega))):
                    qboom += prop.B_y[i]*prop.A_stiff/prop.n_st
            return(qboom)
        else:
            for i in range(int((prop.n_st + 1)/2), prop.n_st-1, 1):
                if(prop.B_z[i] > (prop.Ca - h - s*np.cos(omega))):
                    qboom += (prop.B_z[i] - z_c)*prop.A_stiff/prop.n_st
            return(qboom)
    else:
        return(0)


def get_qy(section, Sy, Izz, f, N, a, b, qs0 = 0):
    if(abs(b)>abs(a)):
         q = (get_integral(f,N,a,b) + get_qboom(section, "y", abs(b)))
    else:
         q = (get_integral(f,N,a,b) + get_qboom(section, "y", abs(a)))
    return((-Sy/Izz *  q) + qs0)

def get_qz(section, Sz, Iyy, f, N, a, b, qs0 = 0):  
    if(abs(b)>abs(a)):
         q = (get_integral(f,N,a,b)+ get_qboom(section, "z", abs(b)))
        #  q = (get_integral(f,N,a,b))
    else:
         q = (get_integral(f,N,a,b) + get_qboom(section, "z", abs(a)))
        #  q = (get_integral(f,N,a,b))
    return((-Sz/Iyy *  q) + qs0)

