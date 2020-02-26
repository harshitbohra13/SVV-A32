import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default='browser'
import pandas as pd

### ==========  IMPORT VALIDATION DATA FROM FILES =============================

from Validation_datagathering import node_loc, element_nodes
from Validation_datagathering import dfbending, dfjambent, dfjamstraight

""" 

"""


HL_dy_bending = dfbending[dfbending.z==0]
HL_dy_bending = HL_dy_bending[HL_dy_bending.y==0]
HL_dy_bending = HL_dy_bending.sort_values(by=['x']) # sort the nodes by x value (necessary for plotting)


HL_dy_jambent = dfjambent[dfjambent.z==0]
HL_dy_jambent = HL_dy_jambent[HL_dy_jambent.y==0]
HL_dy_jambent = HL_dy_jambent.sort_values(by=['x']) # sort the nodes by x value (necessary for plotting)

HL_dy_jamstraight = dfjamstraight[dfjamstraight.z==0]
HL_dy_jamstraight = HL_dy_jamstraight[HL_dy_jamstraight.y==0]
HL_dy_jamstraight = HL_dy_jamstraight.sort_values(by=['x']) # sort the nodes by x value (necessary for plotting)
      
## plot the spanwise deflection in y direction
plt.plot(HL_dy_bending['x'],HL_dy_bending['y+Uy'],label='bending')
plt.plot(HL_dy_jambent['x'],HL_dy_jambent['y+Uy'],label='jam_bent')
plt.plot(HL_dy_jamstraight['x'],HL_dy_jamstraight['y+Uy'],label='jam_straight')
plt.legend()
plt.grid()
plt.ylabel('Deflection in Y [mm]')
plt.xlabel('Spanwise location [mm]')
#plt.show()



