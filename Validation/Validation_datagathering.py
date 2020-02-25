### VALIDATION ####

""" THIS SCRIPT IS NOT FINISHED YET """

""" this script imports data from the B737.rpt file 
    and creates arrays of the tabular data from the file """ 

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


B737rpt_data = open('B737_data/B737.rpt','r')
lines = B737rpt_data.readlines()
B737rpt_data.close()


B737inp_data = open("B737_data/B737.inp","r")
lines_inp = B737inp_data.readlines()

data_points = []
for i in range(9,6597): #get the coordinates of the data points
    split_lines = lines_inp[i].split(",")
    x_loc = float(split_lines[1])
    y_loc = float(split_lines[2])
    z_loc = float(split_lines[3])
    data_points.append([x_loc, y_loc, z_loc])
print("done with data points")

node_sets = []
for j in range(6598,13232): #gives the points belonging to each node
    split_lines = lines_inp[j].split(",")
    n1 = split_lines[1]
    n2 = split_lines[2]
    n3 = split_lines[3]
    n4 = split_lines[4].rstrip('\n')
    node_sets.append([n1, n2, n3, n4])
data_array = np.array(node_sets)
print("done with nodes")

##=============  STRESSES BENDING REGION 1  ================

elementlabel    = []
S_mises         = []
S_s12           = []
for i in range(20,5797):
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


