### VALIDATION ####

""" THIS SCRIPT IS NOT FINISHED YET """

""" this script imports data from the B737.rpt file 
    and creates arrays of the tabular data from the file """ 

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import csv
import pandas as pd


B737rpt_data = open('B737_data/B737.rpt','r')
lines = B737rpt_data.readlines()
B737rpt_data.close()


B737inp_data = open("B737_data/B737.inp","r")
lines_inp = B737inp_data.readlines()

data_points = []
for i in range(9,6597): #get the coordinates of the data points
    split_lines = lines_inp[i].split(",")
    nodelabel = int(split_lines[0])
    x_loc = float(split_lines[1])
    y_loc = float(split_lines[2])
    z_loc = float(split_lines[3])-102.5 # place origin at leading edge to match reference frame of validation model with that of numerical model
    z_loc = float(split_lines[3]) # place origin at leading edge to match reference frame of validation model with that of numerical model
    data_points.append([nodelabel, x_loc, y_loc, z_loc])
node_loc = np.array(data_points)

node_sets = []
for j in range(6598,13232): #gives the points belonging to each node
    split_lines = lines_inp[j].split(",")
    elementlabel = int(split_lines[0])
    n1 = int(split_lines[1])
    n2 = int(split_lines[2])
    n3 = int(split_lines[3])
    n4 = int(split_lines[4].rstrip('\n'))
    node_sets.append([elementlabel, n1, n2, n3, n4])
element_nodes = np.array(node_sets)

##=============  STRESSES BENDING REGION 1  ================

elementlabel    = []
S_mises         = []
S_s12           = []
for i in range(20,5798):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    line.remove(line[1])    
    Smises = np.average([line[1],line[2]]) # computes average of S_mises @loc1 and S_mises @Loc2
    Ss12   = np.average([line[3],line[4]]) # computes average of S_S12   @Loc1 and S_S12   @Loc2    
    elementlabel.append(line[0])
    S_mises.append(Smises)
    S_s12.append(Ss12)
bending_region1 = np.transpose(np.array([elementlabel,S_mises,S_s12]))


##=============  STRESSES BENDING REGION 2  ================

elementlabel    = []
S_mises         = []
S_s12           = []
for i in range(5816,6672):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    line.remove(line[1])    
    Smises = np.average([line[1],line[2]]) # computes average of S_mises @loc1 and S_mises @Loc2
    Ss12   = np.average([line[3],line[4]]) # computes average of S_S12   @Loc1 and S_S12   @Loc2    
    elementlabel.append(line[0])
    S_mises.append(Smises)
    S_s12.append(Ss12)
bending_region2 = np.transpose(np.array([elementlabel,S_mises,S_s12]))


##=============  STRESSES JAM BENT REGION 1  ================

elementlabel    = []
S_mises         = []
S_s12           = []
for i in range(6705,12483):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    line.remove(line[1])    
    Smises = np.average([line[1],line[2]]) # computes average of S_mises @loc1 and S_mises @Loc2
    Ss12   = np.average([line[3],line[4]]) # computes average of S_S12   @Loc1 and S_S12   @Loc2    
    elementlabel.append(line[0])
    S_mises.append(Smises)
    S_s12.append(Ss12)
jambent_region1 = np.transpose(np.array([elementlabel,S_mises,S_s12])) 


##=============  STRESSES JAM BENT REGION 2  ================

elementlabel    = []
S_mises         = []
S_s12           = []
for i in range(12501,13357):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    line.remove(line[1])    
    Smises = np.average([line[1],line[2]]) # computes average of S_mises @loc1 and S_mises @Loc2
    Ss12   = np.average([line[3],line[4]]) # computes average of S_S12   @Loc1 and S_S12   @Loc2    
    elementlabel.append(line[0])
    S_mises.append(Smises)
    S_s12.append(Ss12)
jambent_region2 = np.transpose(np.array([elementlabel,S_mises,S_s12])) 


##=============  STRESSES JAM STRAIGHT REGION 1  ================

elementlabel    = []
S_mises         = []
S_s12           = []
for i in range(13390,19168):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    line.remove(line[1])    
    Smises = np.average([line[1],line[2]]) # computes average of S_mises @loc1 and S_mises @Loc2
    Ss12   = np.average([line[3],line[4]]) # computes average of S_S12   @Loc1 and S_S12   @Loc2    
    elementlabel.append(line[0])
    S_mises.append(Smises)
    S_s12.append(Ss12)
jamstraight_region1 = np.transpose(np.array([elementlabel,S_mises,S_s12])) 


##=============  STRESSES JAM STRAIGHT REGION 2  ================

elementlabel    = []
S_mises         = []
S_s12           = []
for i in range(19186,20042):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    line.remove(line[1])    
    Smises = np.average([line[1],line[2]]) # computes average of S_mises @loc1 and S_mises @Loc2
    Ss12   = np.average([line[3],line[4]]) # computes average of S_S12   @Loc1 and S_S12   @Loc2    
    elementlabel.append(line[0])
    S_mises.append(Smises)
    S_s12.append(Ss12)
jamstraight_region2 = np.transpose(np.array([elementlabel,S_mises,S_s12])) 




##=============  NODAL DISPLACEMENT BENDING  ================

nodelabel=[]; U_magnitude=[]; U_u1=[]; U_u2=[]; U_u3=[]
for i in range(20074,26662):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    u_magnitude = line[1]
    u_u1 = line[2]; u_u2 = line[3]; u_u3 = line[4]
    nodelabel.append(node)
    U_magnitude.append(u_magnitude); 
    U_u1.append(u_u1); U_u2.append(u_u2); U_u3.append(u_u3)
bending_nodal_U = np.transpose(np.array([nodelabel,U_magnitude,U_u1,U_u2,U_u3]))


##=============  NODAL DISPLACEMENT BENDING (ASSEMBLY) ================

nodelabel=[]; U_magnitude=[]; U_u1=[]; U_u2=[]; U_u3=[]
for i in range(26678,26694):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    u_magnitude = line[1]
    u_u1 = line[2]; u_u2 = line[3]; u_u3 = line[4]
    nodelabel.append(node)
    U_magnitude.append(u_magnitude); 
    U_u1.append(u_u1); U_u2.append(u_u2); U_u3.append(u_u3)
bending_nodal_U_assembly = np.transpose(np.array([nodelabel,U_magnitude,U_u1,U_u2,U_u3]))


##=============  NODAL DISPLACEMENT JAM BENT  ================

nodelabel=[]; U_magnitude=[]; U_u1=[]; U_u2=[]; U_u3=[]
for i in range(26724,33312):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    u_magnitude = line[1]
    u_u1 = line[2]; u_u2 = line[3]; u_u3 = line[4]
    nodelabel.append(node)
    U_magnitude.append(u_magnitude); 
    U_u1.append(u_u1); U_u2.append(u_u2); U_u3.append(u_u3)
jambent_nodal_U = np.transpose(np.array([nodelabel,U_magnitude,U_u1,U_u2,U_u3]))


##=============  NODAL DISPLACEMENT JAM BENT (ASSEMBLY)  ================

nodelabel=[]; U_magnitude=[]; U_u1=[]; U_u2=[]; U_u3=[]
for i in range(33328,33344):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    u_magnitude = line[1]
    u_u1 = line[2]; u_u2 = line[3]; u_u3 = line[4]
    nodelabel.append(node)
    U_magnitude.append(u_magnitude); 
    U_u1.append(u_u1); U_u2.append(u_u2); U_u3.append(u_u3)
jambent_nodal_U_assembly = np.transpose(np.array([nodelabel,U_magnitude,U_u1,U_u2,U_u3]))


##=============  NODAL DISPLACEMENT JAM STRAIGHT  ================

nodelabel=[]; U_magnitude=[]; U_u1=[]; U_u2=[]; U_u3=[]
for i in range(33374,39962):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    u_magnitude = line[1]
    u_u1 = line[2]; u_u2 = line[3]; u_u3 = line[4]
    nodelabel.append(node)
    U_magnitude.append(u_magnitude); 
    U_u1.append(u_u1); U_u2.append(u_u2); U_u3.append(u_u3)
jamstraight_nodal_U = np.transpose(np.array([nodelabel,U_magnitude,U_u1,U_u2,U_u3]))


##=============  NODAL DISPLACEMENT JAM STRAIGHT (ASSEMBLY)  ================

nodelabel=[]; U_magnitude=[]; U_u1=[]; U_u2=[]; U_u3=[]
for i in range(39978,39994):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    u_magnitude = line[1]
    u_u1 = line[2]; u_u2 = line[3]; u_u3 = line[4]
    nodelabel.append(node)
    U_magnitude.append(u_magnitude); 
    U_u1.append(u_u1); U_u2.append(u_u2); U_u3.append(u_u3)
jamstraight_nodal_U_assembly = np.transpose(np.array([nodelabel,U_magnitude,U_u1,U_u2,U_u3]))


##=============  NODAL REACTION FORCE BENDING  ================

nodelabel=[]; RF_magnitude=[]; RF1=[]; RF2=[]; RF3=[]
for i in range(40024,46612):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    rf_magnitude = line[1]
    rf1 = line[2]; rf2 = line[3]; rf3 = line[4]
    nodelabel.append(node)
    RF_magnitude.append(rf_magnitude); 
    RF1.append(rf1); RF2.append(rf2); RF3.append(rf3)
bending_nodal_RF = np.transpose(np.array([nodelabel,RF_magnitude,RF1,RF2,RF3]))


##=============  NODAL REACTION FORCE BENDING (ASSEMBLY)  ================

nodelabel=[]; RF_magnitude=[]; RF1=[]; RF2=[]; RF3=[]
for i in range(46628,46644):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    rf_magnitude = line[1]
    rf1 = line[2]; rf2 = line[3]; rf3 = line[4]
    nodelabel.append(node)
    RF_magnitude.append(rf_magnitude); 
    RF1.append(rf1); RF2.append(rf2); RF3.append(rf3)
bending_nodal_RF_assembly = np.transpose(np.array([nodelabel,RF_magnitude,RF1,RF2,RF3]))


##=============  NODAL REACTION FORCE JAM BENT  ================

nodelabel=[]; RF_magnitude=[]; RF1=[]; RF2=[]; RF3=[]
for i in range(46674,53262):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    rf_magnitude = line[1]
    rf1 = line[2]; rf2 = line[3]; rf3 = line[4]
    nodelabel.append(node)
    RF_magnitude.append(rf_magnitude); 
    RF1.append(rf1); RF2.append(rf2); RF3.append(rf3)
jambent_nodal_RF = np.transpose(np.array([nodelabel,RF_magnitude,RF1,RF2,RF3]))


##=============  NODAL REACTION FORCE JAM BENT (ASSEMBLY)  ================

nodelabel=[]; RF_magnitude=[]; RF1=[]; RF2=[]; RF3=[]
for i in range(53278,53294):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    rf_magnitude = line[1]
    rf1 = line[2]; rf2 = line[3]; rf3 = line[4]
    nodelabel.append(node)
    RF_magnitude.append(rf_magnitude); 
    RF1.append(rf1); RF2.append(rf2); RF3.append(rf3)
jambent_nodal_RF_assembly = np.transpose(np.array([nodelabel,RF_magnitude,RF1,RF2,RF3]))


##=============  NODAL REACTION FORCE JAM STRAIGHT  ================

nodelabel=[]; RF_magnitude=[]; RF1=[]; RF2=[]; RF3=[]
for i in range(53324,59912):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    rf_magnitude = line[1]
    rf1 = line[2]; rf2 = line[3]; rf3 = line[4]
    nodelabel.append(node)
    RF_magnitude.append(rf_magnitude); 
    RF1.append(rf1); RF2.append(rf2); RF3.append(rf3)
jamstraight_nodal_RF = np.transpose(np.array([nodelabel,RF_magnitude,RF1,RF2,RF3]))


##=============  NODAL REACTION JAM STRAIGHT (ASSEMBLY)  ================

nodelabel=[]; RF_magnitude=[]; RF1=[]; RF2=[]; RF3=[]
for i in range(59928,59944):
    line = lines[i].split(' ')
    line = [float(j) for j in line if j]
    node = line[0]
    rf_magnitude = line[1]
    rf1 = line[2]; rf2 = line[3]; rf3 = line[4]
    nodelabel.append(node)
    RF_magnitude.append(rf_magnitude); 
    RF1.append(rf1); RF2.append(rf2); RF3.append(rf3)
jamstraight_nodal_RF_assembly = np.transpose(np.array([nodelabel,RF_magnitude,RF1,RF2,RF3]))


bending_stress = np.vstack((bending_region1,bending_region2))
jambent_stress = np.vstack((jambent_region1,jambent_region2))
jamstraight_stress = np.vstack((jamstraight_region1,jamstraight_region2))

bending     = {"stress":bending_stress,  
               "U":bending_nodal_U, "Uassembly":bending_nodal_U_assembly,
               "RF":bending_nodal_RF,"RFassembly":bending_nodal_RF_assembly}
          
jambent     = {"stress":jambent_stress, 
               "U":jambent_nodal_U, "Uassembly":jambent_nodal_U_assembly,
               "RF":jambent_nodal_RF,"RFassembly":jambent_nodal_RF_assembly}

jamstraight = {"stress":jamstraight_stress, 
               "U":jamstraight_nodal_U, "Uassembly":jamstraight_nodal_U_assembly,
               "RF":jamstraight_nodal_RF,"RFassembly":jamstraight_nodal_RF_assembly}

#tab = []
#for node in node_loc[:,0]:
#    eltab = []
#    eltab = [int(node)]
#    for element in element_nodes:
#        if node in element[1:5]:
#            eltab.append(element[0])
#    tab.append(eltab)
#
#with open('nodes.txt', mode='w') as nodes_file:
#    nodes_writer = csv.writer(nodes_file, delimiter=',')
#    for i in range(len(tab)):
#        nodes_writer.writerow(tab[i])

   
        
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



### ------------------------------------------------------------------------###
###          load case: bending     
### ------------------------------------------------------------------------###

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

dfbending = pd.DataFrame({
        'node'      : node_loc[:,0], 
        'x'         : node_loc[:,1], 
        'y'         : node_loc[:,2],
        'z'         : node_loc[:,3],
        'Smises'    : Smises_avg,
        'Ss12'      : Ss12_avg,
        'U'         : bending["U"][:,1],
        'Ux'        : bending["U"][:,2],
        'Uy'        : bending["U"][:,3],
        'Uz'        : bending["U"][:,4],
        'RF'        : bending["RF"][:,1],
        'RFx'       : bending["RF"][:,2],
        'RFy'       : bending["RF"][:,3],
        'RFz'       : bending["RF"][:,4],
        'x+Ux'      : node_loc[:,1]+bending["U"][:,2],
        'y+Uy'      : node_loc[:,2]+bending["U"][:,3],
        'z+Uz'      : node_loc[:,3]+bending["U"][:,4]
        })
    
### ------------------------------------------------------------------------###
###          load case: jam bent     
### ------------------------------------------------------------------------###
    
Smises_avg = []
Ss12_avg   = []
for node in node_elements:
    Smisestab = []
    Ss12tab   = []
    for element in node[1:]:
        index = np.where(jambent["stress"][:,0]==element)
        index = index[0][0]
        Smises = jambent["stress"][index,1]
        Ss12   = jambent["stress"][index,2]
        Smisestab.append(Smises)
        Ss12tab.append(Ss12)
    average_Smises = np.average(np.array(Smisestab))
    average_Ss12 = np.average(np.array(Ss12tab))
    Smises_avg.append(average_Smises)
    Ss12_avg.append(average_Ss12)

dfjambent = pd.DataFrame({
        'node'      : node_loc[:,0], 
        'x'         : node_loc[:,1], 
        'y'         : node_loc[:,2],
        'z'         : node_loc[:,3],
        'Smises'    : Smises_avg,
        'Ss12'      : Ss12_avg,
        'U'         : jambent["U"][:,1],
        'Ux'        : jambent["U"][:,2],
        'Uy'        : jambent["U"][:,3],
        'Uz'        : jambent["U"][:,4],
        'RF'        : jambent["RF"][:,1],
        'RFx'       : jambent["RF"][:,2],
        'RFy'       : jambent["RF"][:,3],
        'RFz'       : jambent["RF"][:,4],
        'x+Ux'      : node_loc[:,1]+jambent["U"][:,2],
        'y+Uy'      : node_loc[:,2]+jambent["U"][:,3],
        'z+Uz'      : node_loc[:,3]+jambent["U"][:,4]
        })

    
### ------------------------------------------------------------------------###
###          load case: jam straight     
### ------------------------------------------------------------------------###
    
Smises_avg = []
Ss12_avg   = []
for node in node_elements:
    Smisestab = []
    Ss12tab   = []
    for element in node[1:]:
        index = np.where(jamstraight["stress"][:,0]==element)
        index = index[0][0]
        Smises = jamstraight["stress"][index,1]
        Ss12   = jamstraight["stress"][index,2]
        Smisestab.append(Smises)
        Ss12tab.append(Ss12)
    average_Smises = np.average(np.array(Smisestab))
    average_Ss12 = np.average(np.array(Ss12tab))
    Smises_avg.append(average_Smises)
    Ss12_avg.append(average_Ss12)

dfjamstraight = pd.DataFrame({
        'node'      : node_loc[:,0], 
        'x'         : node_loc[:,1], 
        'y'         : node_loc[:,2],
        'z'         : node_loc[:,3],
        'Smises'    : Smises_avg,
        'Ss12'      : Ss12_avg,
        'U'         : jamstraight["U"][:,1],
        'Ux'        : jamstraight["U"][:,2],
        'Uy'        : jamstraight["U"][:,3],
        'Uz'        : jamstraight["U"][:,4],
        'RF'        : jamstraight["RF"][:,1],
        'RFx'       : jamstraight["RF"][:,2],
        'RFy'       : jamstraight["RF"][:,3],
        'RFz'       : jamstraight["RF"][:,4],        
        'x+Ux'      : node_loc[:,1]+jamstraight["U"][:,2],
        'y+Uy'      : node_loc[:,2]+jamstraight["U"][:,3],
        'z+Uz'      : node_loc[:,3]+jamstraight["U"][:,4]
        })

### ===========================================================================

