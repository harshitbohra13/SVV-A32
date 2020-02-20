#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:41:42 2020

@author: daanwitte
"""
import numpy as np
import math
import SVV_structural_properties as prop

#to compute shear center 
lsk = np.sqrt((prop.c_a-0.5*prop.h_a)**2 + (0.5*prop.h_a)**2)
Sy = 1 

def get_qarc1():
    theta = np.arange(0, np.pi/2, 0.01)
    qarc1 = 0 
    for dtheta in theta: 
        qarc1 = qarc1 + prop.t_sk * prop.h_a/2 * np.sin(dtheta) * prop.h_a/2 * 0.01 
        
    print(qarc1)
        
    qarc1 = qarc1 * Sy/prop.I_zz
    return (qarc1)

def get_qsec2():
    qsec2 = 0
    y = np.arange(0, prop.h_a/2, 0.01)
    for dy in y:
        qsec2 = qsec2 + prop.t_sp * dy * 0.01
    
    qsec2 = qsec2*Sy/prop.I_zz
    return (qsec2)

#the slope after semi circle "\"
def get_qsec3():
    s = np.arange(0, lsk, 0.01)
    qsec3 = 0
    for ds in s:
        qsec3 = qsec3 + prop.t_sk*(prop.h_a/2 - (prop.h_a/2)*ds/lsk)*0.01
                                   
    qsec3 = qsec3*Sy/prop.I_zz + get_qarc1() + get_qsec2()    
    return qsec3()

def get_qsec4():
    s = np.arange(0, lsk, 0.01)
    qsec4 = 0
    for ds in s:
        qsec4 = qsec4 + prop.t_sk*((-prop.h_a/2*ds)/lsk)*0.01
                                   
    qsec4 = qsec4*Sy/prop.I_zz + get_qsec3()    
    return (qsec4)    

def get_qsec5():
    qsec5 = 0
    y = np.arange(0, prop.h_a/2, 0.01)
    for dy in y:
        qsec5 = qsec5 + prop.t_sp * dy * 0.01
    
    qsec5 = -qsec5
    
    qsec5 = qsec5*Sy/prop.I_zz
    return (qsec5)

def get_qsec6():
    theta = np.arange(0, -(np.pi/2), 0.01)
    qarc6 = 0 
    for dtheta in theta: 
        qarc6 = qarc6 + prop.t_sk * prop.h_a/2 * np.sin(dtheta) * prop.h_a/2 * 0.01 
    qarc6 = -qarc6
    qarc6 = qarc6 *Sy/prop.I_zz + get_qsec4() - get_qsec5()
    return (qarc6)




















# r = 0.0805
# t = 0.0011
# Vy = 1
# x = 10

# #The following will discretize the arc, by dividing it into x amount of equally long straight lines
# #It will also find the location of these lines which can be used for shear centre calculations
# #Dividing it into 10 lines for example, means 5 lines per quarter circle
# #The first coordinates in the upcoming forloop, equals the coordinates of the line at the top of the arc
# #The coordinates are measured from the center of the circle/arc
# RadiansPerLine = m.pi/x

# #The following matrix contains the x and y coordinates of each line at their starting and ending points
# LocationMatrix = np.zeros((x,4))

# for i in range(1,x+1):
#     LocationMatrix[i-1,0] = r*m.cos(RadiansPerLine*i+-RadiansPerLine+m.pi/2)
#     LocationMatrix[i-1,1] = r*m.cos(RadiansPerLine*i+m.pi/2)
#     LocationMatrix[i-1,2] = r*m.sin(RadiansPerLine*i-RadiansPerLine+m.pi/2)
#     LocationMatrix[i-1,3] = r*m.sin(RadiansPerLine*i+m.pi/2)

