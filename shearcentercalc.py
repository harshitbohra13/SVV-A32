#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:41:42 2020

@author: daanwitte
"""
import numpy as np
import math as m
r = 0.0805
t = 0.0011
Vy = 1
x = 10

#The following will discretize the arc, by dividing it into x amount of equally long straight lines
#It will also find the location of these lines which can be used for shear centre calculations
#Dividing it into 10 lines for example, means 5 lines per quarter circle
#The first coordinates in the upcoming forloop, equals the coordinates of the line at the top of the arc
#The coordinates are measured from the center of the circle/arc
RadiansPerLine = m.pi/x

#The following matrix contains the x and y coordinates of each line at their starting and ending points
LocationMatrix = np.zeros((x,4))

for i in range(1,x+1):
    LocationMatrix[i-1,0] = r*m.cos(RadiansPerLine*i+-RadiansPerLine+m.pi/2)
    LocationMatrix[i-1,1] = r*m.cos(RadiansPerLine*i+m.pi/2)
    LocationMatrix[i-1,2] = r*m.sin(RadiansPerLine*i-RadiansPerLine+m.pi/2)
    LocationMatrix[i-1,3] = r*m.sin(RadiansPerLine*i+m.pi/2)

