import numpy as np
import matplotlib.pyplot as plt


#approach:
# first try to find where in the grid you are
# then, find which coefficient you need to interpolate between points
# then linear

#setting up an example function to interpolate
x_tab = np.arange(0,7.001,0.001)
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
h = 0.0005
#amount of sections:
total_sections = len(x_tab) - 1
for i in range(total_sections):
    y1 = interpol_func(x_tab[i])
    y2 = interpol_func(x_tab[i+1])

    #area calculation(trapezoid approximation)
    area = 0.5*(y1 + y2)*h


