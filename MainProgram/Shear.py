#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:41:42 2020
@author: daanwitte
@editor: harshitbohra
This file consists shear flow calculations.
Stiffeners are reffered as booms since they are treated as point masses with areas
"""
# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:41:42 2020
@author: daanwitte
@editor: harshitbohra
This file consists shear flow calculations.
Stiffeners are reffered as booms since they are treated as point masses with areas
"""
import numpy as np
from MainProgram import SVV_structural_properties as prop

# to compute shear center
lsk = np.sqrt((prop.Ca - 0.5 * prop.ha) ** 2 + (0.5 * prop.ha) ** 2)
Sy = 1
area1 = np.pi / 2 * (prop.ha / 2) ** 2
area2 = 1 / 2 * prop.ha * (prop.Ca - prop.ha / 2)


# def get_redq():
# q1 = Symbol('q1')
# q2 = Symbol('q2')
# solve([2*q1*prop.ha/prop.t_sp = 2 * prop.ha*q2/prop.t_sp + q2*2*lsk/prop.t_sk, q1*(prop.ha*np.pi/prop.t_sk + 2*ha/tsk)], q1, )

def get_qbooms():
    '''
    This function  is used to calculate the q in all stiffeners
    '''
    q_booms = np.zeros(prop.n_st)
    for i in range(prop.n_st):
        q_booms[i] = 1 / prop.I_zz * prop.B_y[i] * prop.A_stiff / prop.n_st
    dq_booms = np.zeros(prop.n_st)

    dq_booms[10] = q_booms[1] - q_booms[10]

    for i in range(3, prop.n_st - 1, 1):
        if (i < 6):
            dq_booms[i] = q_booms[i] + q_booms[i - 1]
        if (i >= 6 and i < 10):
            dq_booms[i] = q_booms[i] - q_booms[i - 1]
    q_booms = dq_booms + q_booms

    return (q_booms)


# sec - 1 (^
def get_qsec1():
    delta = 0.01
    theta = np.arange(0, np.pi / 2, delta)
    qsec1 = 0
    for dtheta in theta:
        qsec1 = qsec1 + prop.t_sk * prop.ha / 2 * np.sin(dtheta) * prop.ha / 2 * delta

    qsec1 = qsec1 * Sy / prop.I_zz
    return (qsec1)


# sec - 2 |^
def get_qsec2():
    qsec2 = 0
    delta = 0.01
    y = np.arange(0, prop.ha / 2, delta)
    for dy in y:
        qsec2 = qsec2 + prop.t_sp * dy * delta

    qsec2 = qsec2 * Sy / prop.I_zz
    return (qsec2)


# the slope after semi circle "\"
def get_qsec3():
    delta = 0.01
    s = np.arange(0, lsk, delta)
    qsec3 = 0
    for ds in s:
        qsec3 = qsec3 + prop.t_sk * (prop.ha / 2 - (prop.ha / 2) * ds / lsk) * delta

    qsec3 = qsec3 * Sy / prop.I_zz + get_qsec1() + get_qsec2()
    return (qsec3)


# the slope after semi circle "/"
def get_qsec4():
    delta = 0.01
    s = np.arange(0, lsk, delta)
    qsec4 = 0
    for ds in s:
        qsec4 = qsec4 + prop.t_sk * ((-prop.ha / 2 * ds) / lsk) * delta

    qsec4 = qsec4 * Sy / prop.I_zz + get_qsec3()
    return (qsec4)


# spar bot "|"
def get_qsec5():
    delta = 0.01
    qsec5 = 0
    y = np.arange(0, prop.ha / 2, delta)
    for dy in y:
        qsec5 = qsec5 + prop.t_sp * dy * delta
    qsec5 = -qsec5

    qsec5 = qsec5 * Sy / prop.I_zz
    return (qsec5)


# bot semi circle
def get_qsec6():
    delta = 0.01
    theta = np.arange(0, -(np.pi / 2), delta)
    qarc6 = 0
    for dtheta in theta:
        qarc6 = qarc6 + prop.t_sk * prop.ha / 2 * np.sin(dtheta) * prop.ha / 2 * delta
    qarc6 = -qarc6
    qarc6 = qarc6 * Sy / prop.I_zz + get_qsec4() - get_qsec5()
    return (qarc6)


def get_qbcell():
    '''
    Finds qb of both the cells
    '''
    qbcell = []
    qbcell.append(get_qsec1())
    qbcell.append(get_qsec2())
    qbcell.append(get_qsec3())
    qbcell.append(get_qsec4())
    qbcell.append(get_qsec5())
    qbcell.append(get_qsec6())
    return (qbcell)


def get_qs0():
    '''
    Basically finds the redundant shear flow created by stiffeners
    and qb of sections.
    '''
    qb, ds_t = get_intqb()
    qbooms = get_qbooms()
    qs0 = (sum(qb) + sum(qbooms)) / (sum(ds_t) + prop.A_stiff)
    return (qs0)


def get_intqb():
    '''
    This function calculates the integral of qb with respect to ds
    Also ds as t is not constant.
    This is done for the compatibility equation, around shear center
    '''
    qb = get_qbcell()
    ds_t = np.zeros(6)
    for i in range(1, 7):
        if (i == 1 or i == 6):
            ds_t[i - 1] = (np.pi / 2) * (prop.ha / 2) / prop.t_sk
            qb[i - 1] = qb[i - 1] * ds_t[i - 1]
        if (i == 2 or i == 5):
            ds_t[i - 1] = (prop.ha / 2) / prop.t_sp
            qb[i - 1] = qb[i - 1] * ds_t[i - 1]
        if (i == 3 or i == 4):
            ds_t[i - 1] = (lsk) / prop.t_sk
            qb[i - 1] = qb[i - 1] * ds_t[i - 1]
    return (qb, ds_t)


def get_sc():
    '''
    This function calculates the sc from of z from center of the spar
    It basically uses the shear equation, lht is integral of (p*q_b*ds)
    rht is summation of (2*a*qs0) with respect to areas
    '''
    rht = 2 * area1 * get_qs0() + 2 * area2 * get_qs0()
    qbooms = get_qbooms()
    theta = np.zeros(round(np.pi / 2 * 1000))
    theta.fill(0.001)
    theta1 = np.zeros(round(np.pi / 2 * 1000))
    theta1.fill(0.001)
    qbo = 0
    for i in range(len(qbooms)):
        qbo = qbo + (prop.B_z[i] * qbooms[i])
    lht = [sum(prop.ha / 2 * get_qsec1() * prop.ha / 2 * dtheta for dtheta in theta),
           sum(prop.ha / 2 * get_qsec6() * prop.ha / 2 * -dtheta for dtheta in theta1),
           get_qsec2() * prop.ha / 2 * prop.t_sp * 0.5 * prop.ha,
           get_qsec5() * prop.ha / 2 * prop.t_sp * 0.5 * prop.ha,
           get_qsec3() * lsk * prop.t_sk * (prop.Ca - prop.ha / 2) * prop.ha / 2 / lsk,
           get_qsec4() * lsk * prop.t_sk * (prop.Ca - prop.ha / 2) * prop.ha / 2 / lsk,
           qbo]
    return (sum(lht) + rht, 0)


print(get_qs0())
print("sc(y,z) = ", get_sc())

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

print(get_qs0())
print("sc(y,z) = ", get_sc())