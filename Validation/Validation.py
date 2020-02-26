import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default='browser'
import pandas as pd

### ==========  IMPORT VALIDATION DATA FROM FILES =============================

from Validation_datagathering import node_loc, element_nodes
from Validation_datagathering import bending, jambent, jamstraight

""" 
To use the data for three seperate load cases (bending, jambent and jamstraight)

POSSIBLE OUTPUT DATA:
    
stress    : von mises and shear stresses for region 1 (skin) and region 2 (spar)
U          : nodal displacement
Uassembly  : nodal displacement of assmebly
RF         : nodal reaction force
RFassembly : nodal reaction force of assembly       

To use output data for a given loadcase, use e.g. bending["stress1"]  
"""

### IMPORT TABLE CONTAINING INFO OF WHICH ELEMENTS ARE CONNECTED TO WHICH NODE

f = open('nodes.txt','r')
data = f.readlines()
f.close()
node_elements = []
for line in data: 
    line = line.rstrip('\n')
    line = line.split(',')
    line = [float(j) for j in line if j]
    if line:
        node_elements.append(line)
#node_elements = np.array(node_elements)


### ======== SET UP DATAFRAMES FOR EACH LOAD CASE =============================



Smises_avg = []
Ss12_avg   = []
for node in node_elements:
    Smisestab = []
    Ss12tab   = []
    for element in node[1:]:
        index = np.where(bending["stress"][:,0]==element)
        index = index[0][0]
        Smises = bending["stress"][index,1]
        Ss12   = bending["stress"][index,2]
        Smisestab.append(Smises)
        Ss12tab.append(Ss12)
    average_Smises = np.average(np.array(Smisestab))
    average_Ss12 = np.average(np.array(Ss12tab))
    Smises_avg.append(average_Smises)
    Ss12_avg.append(average_Ss12)
    
        
        
        
#        print(bending["stress1"][idx])


### ==== PLOT SPANWISE DEFLECTION IN Y dir DUE TO BENDING =====================

# the following lines create an array of all nodes with corresponding 
# x, y, z locations, and delete all nodes with non-zero z coordinates
# to create an array of all the leading edge (LE) nodes
LEcoordinates = np.array(node_loc)
zloc  = node_loc[:,3]
idx   = 0
for z in zloc:
    if z != 0:        
        LEcoordinates = np.delete(LEcoordinates,idx,0)   
        idx += -1
    idx += 1

# this piece of code generates an array of all LE nodes 
# with y displacement due to the bending load case
LE_dy = LEcoordinates
U_y = bending["U"][:,3]
idx = 0
for node in bending["U"]:
    if node[0] in LEcoordinates[:,0]: # if the node from nodal displacement data is on the LE:
        LE_dy[idx,2] += node[3] # add y deflection to y coordinate of LE node
        idx += 1
        
# put data in dataframe for the sake of easy processing
# the nodes are not sorted yet
LE_dy = pd.DataFrame({
        'node': LE_dy[:,0], 
        'x'   : LE_dy[:,1], 
        'y'   : LE_dy[:,2],
        'z'   : LE_dy[:,3]})
LE_dy = LE_dy.sort_values(by=['x']) # sort the nodes by x value (necessary for plotting)

# plot the spanwise deflection in y direction
plt.plot(LE_dy['x'],LE_dy['y'])
#plt.show()



