import numpy as np
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

print(interval_check(0.3, 0.6, check))




