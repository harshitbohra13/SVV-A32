#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:41:42 2020

@author: daanwitte
@editor: harshitbohra

"""
import numpy as np
import math
import SVV_structural_properties as prop

#to compute shear center 
lsk = np.sqrt((prop.c_a-0.5*prop.h_a)**2 + (0.5*prop.h_a)**2)
Sy = 1 
area1 =  np.pi/2*(prop.h_a/2)**2
area2 =  1/2*prop.h_a*(prop.c_a - prop.h_a/2)

# def get_redq():
    # q1 = Symbol('q1')
    # q2 = Symbol('q2')
    # solve([2*q1*prop.h_a/prop.t_sp = 2 * prop.h_a*q2/prop.t_sp + q2*2*lsk/prop.t_sk, q1*(prop.h_a*np.pi/prop.t_sk + 2*h_a/tsk)], q1, )

def get_qbooms():
    q_booms = np.zeros(prop.n_st)
    for i in range(prop.n_st):
        q_booms[i] = 1/prop.I_zz * prop.B_y[i] * prop.A_stiff/prop.n_st
    return (q_booms)

#sec - 1 (^
def get_qsec1():
    theta = np.arange(0, np.pi/2, 0.01)
    qsec1 = 0 
    for dtheta in theta: 
        qsec1 = qsec1 + prop.t_sk * prop.h_a/2 * np.sin(dtheta) * prop.h_a/2 * 0.01 
                
    qsec1 = qsec1 * Sy/prop.I_zz
    return (qsec1)
#sec - 2 |^
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
                                   
    qsec3 = qsec3*Sy/prop.I_zz + get_qsec1() + get_qsec2()    
    return (qsec3)

#the slope after semi circle "/"
def get_qsec4():
    s = np.arange(0, lsk, 0.01)
    qsec4 = 0
    for ds in s:
        qsec4 = qsec4 + prop.t_sk*((-prop.h_a/2*ds)/lsk)*0.01
                                   
    qsec4 = qsec4*Sy/prop.I_zz + get_qsec3()    
    return (qsec4)    

#spar bot "|"
def get_qsec5():
    qsec5 = 0
    y = np.arange(0, prop.h_a/2, 0.01)
    for dy in y:
        qsec5 = qsec5 + prop.t_sp * dy * 0.01
    
    qsec5 = -qsec5
    
    qsec5 = qsec5*Sy/prop.I_zz
    return (qsec5)

#bot semi circle
def get_qsec6():
    theta = np.arange(0, -(np.pi/2), 0.01)
    qarc6 = 0 
    for dtheta in theta: 
        qarc6 = qarc6 + prop.t_sk * prop.h_a/2 * np.sin(dtheta) * prop.h_a/2 * 0.01 
    qarc6 = -qarc6
    qarc6 = qarc6 *Sy/prop.I_zz + get_qsec4() - get_qsec5()
    return (qarc6)


def get_qbcell1():
    qbcell1 = []
    qbcell1.append(get_qsec1())
    qbcell1.append(get_qsec2())
    qbcell1.append(get_qsec3())
    qbcell1.append(get_qsec4())
    qbcell1.append(get_qsec5())
    qbcell1.append(get_qsec6())
    return(qbcell1)

def get_qs0():
    qb, ds_t = get_intqb()
    qs0 = sum(qb)/sum(ds_t)
    return (qs0)

def get_intqb():
    qb = get_qbcell1()
    ds_t = np.zeros(6)
    for i in range(1,7):
        if(i == 1 or i == 6):
            ds_t[i-1] = (np.pi/2)*(prop.h_a/2)/prop.t_sk
            qb[i-1] = qb[i-1]*ds_t[i-1]
        if(i == 2 or i == 5):
            ds_t[i-1] = (prop.h_a/2)/prop.t_sp
            qb[i-1] = qb[i-1]*ds_t[i-1]
        if(i == 3 or i == 4):
            ds_t[i-1] = (lsk)/prop.t_sk
            qb[i-1] = qb[i-1]*ds_t[i-1]
    return (qb, ds_t)
    
def get_sc():
    rht = 2 * area1 * get_qs0() + 2 * area2 * get_qs0()
    theta  = np.arange(0, np.pi/2, 0.01)
    theta1 = np.arange(-np.pi/2, 0, 0.01)
    qbooms = get_qbooms()  
    qbo = 0
    for i in range(len(qbooms)):
        qbo = qbo + (prop.B_z[i] * qbooms[i])    
    lht =[sum(prop.h_a/2 * get_qsec1() * prop.h_a/2 * dtheta for dtheta in theta),
          sum(prop.h_a/2 * get_qsec6() *  prop.h_a/2 * dtheta for dtheta in theta1),
          qbo]
    return(sum(lht)+ rht)


print(get_qs0())
print(get_sc())
    















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

