import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go

B737_data = open("B737.inp","r")
lines = B737_data.readlines()

data_points = []
for i in range(9,6597): #get the coordinates of the data points
    split_lines = lines[i].split(",")
    x_loc = float(split_lines[1])
    y_loc = float(split_lines[2])
    z_loc = float(split_lines[3])
    data_points.append([x_loc, y_loc, z_loc])
print("done with data points")

node_sets = []
for j in range(6598,13232): #gives the points belonging to each node
    split_lines = lines[j].split(",")
    n1 = split_lines[1]
    n2 = split_lines[2]
    n3 = split_lines[3]
    n4 = split_lines[4].rstrip('\n')
    node_sets.append([n1, n2, n3, n4])
data_array = np.array(node_sets)
print("done with nodes")




#
# fig = go.Figure(data=[go.Scatter3d(x=data_array[:,0], y=data_array[:,1], z=data_array[:,2],
#                                    mode='markers',marker=dict(
#         size=2,
#         opacity=1))])
# fig.show()

