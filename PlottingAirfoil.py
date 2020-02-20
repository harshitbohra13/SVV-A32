import numpy as np
import matplotlib.pyplot as plt

#Variables
Ca = 0.515 #m
h = 0.248/2 #cm

#Step for each point
step=0.00001

#First Section
z0=np.arange(h,Ca+step,step)
y0=Ca*h/(Ca-h)-z0*h/(Ca-h)

#Second Section
z1=z0
y1=-y0

#Third Section
y2=np.arange(0,h+step,step)
z2=np.ones(len(y2))*h

#Fourth Section
y3=-y2
z3=z2

#FIfth Section
z4=np.arange(0,h+step,step)
y4=np.sqrt(-(z4-h)**2+h**2)

#Sixth Section
z5=z4
y5=-y4

#Plot all the sections
plt.plot(z0,y0,'bo')
plt.plot(z1,y1,'ro')
plt.plot(z2,y2,'go')
plt.plot(z3,y3,'ko')
plt.plot(z4,y4,'ko')
plt.plot(z5,y5,'go')
plt.grid()
plt.show()