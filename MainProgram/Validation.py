import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default='browser'
import pandas as pd
import plotly.express as px

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

### ==========  IMPORT VALIDATION DATA FROM FILES =============================
from Validation_datagathering import node_loc, element_nodes
from Validation_datagathering import dfbending, dfjambent, dfjamstraight
import pandas as pd
import plotly.express as px

#
#import sys
#
#sys.path.append('../')

from MainInterpolation import *
#from MainProgram.MainInterpolation import *

""" 

"""
xrange = dfbending[dfbending.z==0]
xrange = xrange[xrange.y==0]
xrange = xrange['x']/1000
xrange = xrange.sort_values()


####################AIRFOIL - STRESSES##########################
#location = 0.4 #[m] #Location of the cross section
h = ha/2 #cm

#Step for each point
step=1000
step=20

#Third Section
z3=np.linspace(h,Ca,step)
y3=Ca*h/(Ca-h)-z3*h/(Ca-h)

#Fourth Section
z4=z3
y4=-y3

#Second Section
y2=np.linspace(0,h,step)
z2=np.ones(len(y2))*h

#Fifth Section
y5=-y2
z5=z2

#First Section
z1=np.linspace(0,h,step)
y1=np.sqrt(-(z1-h)**2+h**2)

#Sixth Section
z6=z1
y6=-y1


output = np.array([0,0,0,0,0])


for x in xrange[1:-1]:
#for x in [0.4,0.8,1.2]:

    #############SHEAR-FLOWS############################################
    q1,q2,q3,q4,q5,q6 = shearstress(x,Sy(x),Sz(x),T(x))
    q4 = q4[::-1]
    q6 = q6[::-1]
    
    
    ##############DIRECT-STRESSES AND VON MISES STRESSES ###############
    def direct_stress(My,Mz,z,y):
        sigma_xx_z = My * (z - z_hat) / Iyy
        sigma_xx_y = Mz * y / Izz
        sigma_xx_total = sigma_xx_y + sigma_xx_z
        return sigma_xx_total
    
    def VonMises(My, Mz,tau_yz,z, y):

        sigma_xx_z = My * (z - z_hat) / Iyy
        sigma_xx_y = Mz * y / Izz
    
        sigma_xx_total = sigma_xx_y + sigma_xx_z
        sigma_vm =  np.sqrt(0.5 * ((sigma_xx_total**2) + ((-sigma_xx_total)**2)) + 3*(tau_yz ** 2))
    
        return sigma_vm
    
    sec_1 = []
    sec_2 = []
    sec_3 = []
    sec_4 = []
    sec_5 = []
    sec_6 = []
    
    for n in range(0,step):
        sec_3.append([x,y3[n],z3[n],direct_stress(Moment_y(x), Moment_z(x), -(z3[n]), y3[n]),VonMises(Moment_y(x), Moment_z(x),q3[n]/t_sk, -(z3[n]), y3[n])])
        sec_4.append([x,y4[n],z4[n],direct_stress(Moment_y(x), Moment_z(x), -(z4[n]), y4[n]),VonMises(Moment_y(x), Moment_z(x),q4[n]/t_sk, -(z4[n]), y4[n])])
        sec_2.append([x,y2[n],z2[n],direct_stress(Moment_y(x), Moment_z(x), -(z2[n]), y2[n]),VonMises(Moment_y(x), Moment_z(x),q2[n]/t_sp, -(z2[n]), y2[n])])
        sec_5.append([x,y5[n],z5[n],direct_stress(Moment_y(x), Moment_z(x), -(z5[n]), y5[n]),VonMises(Moment_y(x), Moment_z(x),q5[n]/t_sp, -(z5[n]), y5[n])])
        sec_1.append([x,y1[n],z1[n],direct_stress(Moment_y(x), Moment_z(x), -(z1[n]), y1[n]),VonMises(Moment_y(x), Moment_z(x),q1[n]/t_sk, -(z1[n]), y1[n])])
        sec_6.append([x,y6[n],z6[n],direct_stress(Moment_y(x), Moment_z(x), -(z6[n]), y6[n]),VonMises(Moment_y(x), Moment_z(x),q6[n]/t_sk, -(z6[n]), y6[n])])

    output = np.vstack((np.array(output),np.array(sec_1)))
    output = np.vstack((output,np.array(sec_2)))
    output = np.vstack((output,np.array(sec_3)))
    output = np.vstack((output,np.array(sec_4)))
    output = np.vstack((output,np.array(sec_5)))
    output = np.vstack((output,np.array(sec_6)))

output = np.delete(output,(0),axis=0)

df = pd.DataFrame({
'x' : output[:,0],
'y' : output[:,1],
'z' : output[:,2],
'Sdirect': output[:,3],
'Smises' : output[:,4]
})

    
print("Runtime: %f seconds" % (time.time()-start_time))


fig = px.scatter_3d(df, x='x', y='y', z='z', color='Smises')
fig.update_layout(title='von mises stress')
fig.show()





#HL_dy_bending = dfbending[dfbending.z==0]
#HL_dy_bending = HL_dy_bending[HL_dy_bending.y==0]
#HL_dy_bending = HL_dy_bending.sort_values(by=['x']) # sort the nodes by x value (necessary for plotting)
#
#
#HL_dy_jambent = dfjambent[dfjambent.z==0]
#HL_dy_jambent = HL_dy_jambent[HL_dy_jambent.y==0]
#HL_dy_jambent = HL_dy_jambent.sort_values(by=['x']) # sort the nodes by x value (necessary for plotting)
#
#HL_dy_jamstraight = dfjamstraight[dfjamstraight.z==0]
#HL_dy_jamstraight = HL_dy_jamstraight[HL_dy_jamstraight.y==0]
#HL_dy_jamstraight = HL_dy_jamstraight.sort_values(by=['x']) # sort the nodes by x value (necessary for plotting)
#      
### plot the spanwise deflection in y direction
#plt.clf()
#plt.plot(HL_dy_bending['x'],HL_dy_bending['y+Uy'],label='bending')
#plt.plot(HL_dy_jambent['x'],HL_dy_jambent['y+Uy'],label='jam_bent')
#plt.plot(HL_dy_jamstraight['x'],HL_dy_jamstraight['y+Uy'],label='jam_straight')
#plt.legend()
#plt.grid()
#plt.ylabel('Deflection in Y [mm]')
#plt.xlabel('Spanwise location [mm]')
##plt.show()
#
#fig = px.scatter_3d(dfbending, x='x+Ux', y='y+Uy', z='z+Uz', color='Smises')
#fig.update_layout(title='Shear stresses due to pure bending')
#fig.show()

