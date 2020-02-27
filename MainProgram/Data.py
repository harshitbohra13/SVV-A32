import numpy as np
#================= Data Fokker 100 =======================

Ca = 0.505     # chord length aileron [m]
la = 1.611     # span of the aileron [m]
ha = 16.1e-2   # aileron height [m]
x1 = 0.125 #m
x2 = 0.498 #m
x3 = 1.494 #m
xa = 0.245 #m
d1 = 0.00389 # m
d3 = 0.01245  # m
P = -49.2*1000  # N
theta = 30*np.pi/180 #rad

t_sk = 1.1e-3   # skin thickness [m]
t_sp = 2.4e-3   # spar thickness [m]

n_st = 11       # number of stiffeners [-]
t_st = 1.2e-3   # thickness of stiffener [m]
h_st = 1.3e-2   # height of stiffener [m]
w_st = 1.7e-2   # width of stiffener [m]

E = 72.9*10**9
G = 27.1*10**9