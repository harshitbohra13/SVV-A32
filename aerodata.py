from math import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default='browser'
import pandas as pd




# Parameters
Ca = 0.505 #[m]
La = 1.611 #[m]

# Reading aerodata -- write as 81x41 matrix

ref_file = open("aerodynamicloadf100.dat", "r")
lines = ref_file.readlines()
ref_file.close()

aerodata = np.mat(np.zeros((81, 41)))


for line in lines:
    idx = lines.index(line)
    line = line.replace("\n","")
    values = line.split(',')
    values = np.matrix(values)
    
    aerodata[idx,:] = values
  


# Calculating x and z coordinates of data points 
theta_x = []  # List of theta_z_i angles
theta_z = []  # List of theta_x_i angles
coord_x = []  # List of x coordinates
coord_z = []  # List of z coordinates


N_x = 41
N_z = 81

for i in range(1,N_x+2):
    theta_i = (i-1)*np.pi/N_x
    theta_x.append(theta_i)
    

for i in range(1,N_x+1):
    x_i = -0.5*(   0.5*La*(1-np.cos(theta_x[i-1])) + 0.5*La*(1-np.cos(theta_x[i]))  )
    coord_x.append(x_i)


for i in range(1,N_z+2):
    theta_i = (i-1)*np.pi/N_z
    theta_z.append(theta_i) 
    
    
for i in range(1,N_z+1):
    z_i = -0.5*(   0.5*Ca*(1-np.cos(theta_z[i-1])) + 0.5*Ca*(1-np.cos(theta_z[i]))  )
    coord_z.append(z_i)
  
    

# Make all x-coordinates positive (positive x-axis starting from root towards tip )
coord_x = np.abs(coord_x)

# Plotting

header = np.arange(81)+1

df = pd.DataFrame(aerodata,header)
df.columns=np.arange(41)+1


layout = go.Layout(
    
title='Aerodynamic Loading on the Aileron Surface',
scene=dict(
xaxis=dict(range=[-0.5,2],title='X Axis'),
yaxis=dict(range=[0.1,-0.7],title='Z Axis'),
zaxis=dict(range=[0,30],title='Pressure Distribution [kPa]')
),
autosize=False,
width=1000,
height=850,
margin=dict(
l=65,
r=50,
b=65,
t=90
)
)

fig = go.Figure(data=[go.Surface(z=df.values,x=coord_x,y=coord_z)],layout=layout)
fig.show()



  

 

    
    
    
    
    
    
    
    
    
    