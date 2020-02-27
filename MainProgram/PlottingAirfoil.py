import numpy as np
import matplotlib.pyplot as plt
from MainProgram import shearforces as shear

#Variables
Ca = 0.515 #m
h = 0.248/2 #cm

#Step for each point
step=0.00001

#Third Section
z3=np.arange(h,Ca+step,step)
y3=Ca*h/(Ca-h)-z3*h/(Ca-h)
q3 = shear.shearstress(loca)

#Fourth Section
z4=z3
y4=-y3

#Second Section
y2=np.arange(0,h+step,step)
z2=np.ones(len(y2))*h

#Fifth Section
y5=-y2
z5=z2

#First Section
z1=np.arange(0,h+step,step)
y1=np.sqrt(-(z1-h)**2+h**2)

#Sixth Section
z6=z1
y6=-y1

#

plt.figure(6)
#Plot all the sections
#plt.plot(z3,y3,'bo')
#plt.plot(z4,y4,'ro')
#plt.plot(z2,y2,'go')
#plt.plot(z5,y5,'ko')
#plt.plot(z1,y1,'ko')
#plt.plot(z6,y6,'go')
plt.grid()
plt.show()