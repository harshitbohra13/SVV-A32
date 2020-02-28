import numpy as np
from Functions_Interpolation import find_interpolant, compute_interpolant
from MainInterpolation import coor_x,coor_z,matrix_data,diff_x
import matplotlib.pyplot as plt


#approach:
# first try to find where in the grid you are
# then, find which coefficient you need to interpolate between points
# then linear

#setting up an example function to interpolate
x_tab = np.arange(0,3.001,0.001)
y_tab = []

# for i in range(len(x_tab)):
#     y = np.sin(x_tab[i])
#     y_tab.append(y)
#
# plt.plot(x_tab,y_tab)
# plt.show()

def interpol_func(x):
    y = np.sin(x)
    return y


#integrating function, linear approximation between two point
h = 0.001

#amount of sections:
total_sections = len(x_tab) - 1
area_total = 0
for i in range(total_sections):
    y1 = interpol_func(x_tab[i])
    y2 = interpol_func(x_tab[i+1])

    #area calculation(trapezoid approximation)
    area = 0.5*(y1 + y2)*h
    area_total += area

l_chord = 0.505
l_span = 1.611
#makes a mesh with n amount of nodes
n = 1000
mesh = np.arange(0, l_chord, 1/n)
check = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]


#function that identifies between which nodes an interval lies
def interval_check(v1, v2, reference):
    for i in range(len(reference)-1):
        if v1 >= reference[i] and v1<= reference[i+1]: #checks if it fits in an interval
            I1_index = i
            break

    for i in range(len(reference)-1):
        if v2 >= reference[i] and v2<= reference[i+1]:
            I2_index = i
            break

    return I1_index, I2_index

def integrate(nodes, coeff_matrix): #nodes is a list of the x or y position. data is aerodynamic data
    area = 0
    for i in range(coeff_matrix.size[0] - 1):
        Ac = coeff_matrix[i][0]
        Bc = coeff_matrix[i][1]
        Cc = coeff_matrix[i][2]
        Dc = coeff_matrix[i][3]

        x1 = nodes[i]
        x2 = nodes[i+1]
        A1 = Ac*(1/4)*(x1**4) + Bc*(1/3)*(x1**3) + Cc*(1/2)*(x1**2) + Dc*x1
        A2 = Ac * (1 / 4) * (x2 ** 4) + Bc * (1 / 3) * (x2 ** 3) + Cc * (1 / 2) * (x2 ** 2) + Dc * x2

        dA = A1 - A2
        area += dA

    return area

print(interval_check(0.3, 0.6, check))


#######Check Interpolation##########

plt.figure(11)
plt.plot(coor_x,matrix_data[0,:],'o')
new_nodes, new_loading = compute_interpolant(coor_x,find_interpolant(coor_x,matrix_data[0,:]),1)
plt.plot(new_nodes,new_loading,'r')
plt.grid()
plt.xlabel('x')
plt.ylabel('Aerodynamic Loading [kPa]')
plt.show()

##############Check Second Derivative##############
coeff = find_interpolant(coor_x,matrix_data[0,:])
matrix_value = np.zeros((len(coor_x)-2,2))
for i in range(0,len(coor_x)-2):
    matrix_value[i,0]=coeff[i,0]*6*(diff_x[i])+2*coeff[i,1]
    matrix_value[i,1]=coeff[i+1,1]*2
print(matrix_value)



