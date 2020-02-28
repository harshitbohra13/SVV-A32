#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 16:41:42 2020

@author: daanwitte
@editor: harshitbohra

This file consists shear flow calculations.
Stiffeners are reffered as booms since they are treated as point masses with areas

"""
import numpy as np
import SVV_structural_properties as prop
import Data as data

#to compute shear center 
h = prop.ha/2

lsk = np.sqrt((prop.Ca-h)**2 + (h)**2)
Sy = 1 
area1 =  np.pi/2*(h)**2
area2 =  1/2*prop.ha*(prop.Ca - h)
boareas = np.zeros(11)
boareas.fill(prop.A_stiff/prop.n_st) 
G = data.G #pascals

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
        q_booms[i] = 1/prop.I_zz * prop.B_y[i] * prop.A_stiff/prop.n_st
    dq_booms = np.zeros(prop.n_st)
    
    dq_booms[10] = q_booms[1] - q_booms[10] 

    for i in range(3, prop.n_st-1, 1):
        if(i<6):
            dq_booms[i] = q_booms[i] + q_booms[i-1]
        if(i>=6 and i<10):
            dq_booms[i] = q_booms[i] - q_booms[i-1]
    q_booms = dq_booms + q_booms

    return(q_booms)
'''
Area  1
'''
#--------------------------------------------------------------------------------------------

def sintegrate(N, a, b):
    def f(x):
        return np.sin(x)
    num0 = 0
    num1 = 0

    for i in range(1, N + 1):
        num0 += f(a + (i-(1/2))*((b-a)/N))
    num1 = ((b-a)/N)*num0

    return num1
    
def get_integral(N, a, b):
    def f(x):
        return x
    num0 = 0
    num1 = 0
    for i in range(1, N + 1):
        num0 += f(a + (i-(1/2))*((a - b)/N))
    num1 = ((a - b)/N)*num0
    
    return num1

#sec - 1 (^
def get_qsec1():
    # delta = 0.01
    # theta = np.arange(0,np.pi/2,delta)
    qsec1 = 0 
    qsec1 = prop.t_sk* h * h * sintegrate(1000, 0, np.pi/2) 
    qsec1 =  ((-1) * Sy/prop.I_zz) *(qsec1 + prop.B_y[1]*prop.A_stiff)
    return (qsec1)

#sec - 2 |^
def get_qsec2():
    qsec2 = 0
    qsec2 = (-1) * Sy/prop.I_zz * prop.t_sp *get_integral(1000, 0, h)  
    return (qsec2)

#spar bot "|"
def get_qsec5():
    qsec5 = 0

    qsec5 = (-1)*Sy/prop.I_zz * prop.t_sp* get_integral(1000, 0, -h)
    return (qsec5)

#bot semi circle
def get_qsec6():
    
    qarc6 = 0 
    qarc6 = prop.t_sk* h * h * sintegrate(1000,  -np.pi/2 ,0) 
    qarc6 =( (-1) *Sy/prop.I_zz *(qarc6 + prop.B_y[prop.n_st - 1] * prop.A_stiff/prop.n_st)) + (get_qsec4() - get_qsec5())
    
    return (qarc6)

#\
def get_qsec3():
    
    qsec3 = prop.t_sk* (h - h/lsk)* get_integral(1000, 0, lsk)
    for i in range(2, int(((prop.n_st)+1)/2), 1):
        qsec3 += prop.B_y[i]*prop.A_stiff/prop.n_st

    qsec3 = (-1)*Sy/prop.I_zz* qsec3 + get_qsec1() +  get_qsec2()    
    return (qsec3)

#the slope after semi circle "/"
def get_qsec4():
    qsec4 = 0 
    qsec4 = prop.t_sk *(h/lsk)* get_integral(1000, 0, lsk)
    for i in range(int((prop.n_st+1)/2), prop.n_st - 1, 1):
        qsec4 += prop.B_y[i]*prop.A_stiff/prop.n_st

    qsec4 = (-1)*Sy/prop.I_zz * qsec4 + get_qsec3()    
    return (qsec4)   



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
    return(qbcell)

def get_qs0():
    '''
    Basically finds the redundant shear flow created by stiffeners
    and qb of sections.
    '''
    qb1, ds1, qb2, ds2 = get_intqb()
    # qbooms = get_qbooms()
    # q01 = (qb1/ds1)
    # q02 = (qb2/ds2)

    x1 = (h)*((np.pi/2)*2) + prop.ha
    x2 = -1*(prop.ha)
    x3 = -1*(prop.ha)
    x4 = prop.ha + 2*lsk

    b1 = (h)*(get_qsec1()*(np.pi/2) + get_qsec6()*(np.pi/2))  + -1*get_qsec2()*(h) + -1*get_qsec5()*h
    b2 = get_qsec2()*h + get_qsec5()*h + get_qsec3()*lsk + get_qsec4()*lsk
    
    b = [-b1,-b2]
        
    matrix = np.array([[x1,x2],[x3,x4]])
    
    X = np.linalg.solve(matrix, b)


    return (X[0], X[1])

def get_intqb():
    '''
    This function calculates the integral of qb with respect to ds
    Also ds as t is not constant. 
    This is done for the compatibility equation, around shear center
    '''
    qb = get_qbcell()
    qb1 = 0
    qb2 = 0
    ds_t = np.zeros(6)
    ds1 = 0
    ds2 = 0

    for i in range(1,7):
        if(i == 1 or i == 6):
            ds_t[i-1] = (np.pi/2)*(h)/G*prop.t_sk
            ds1 += (np.pi/2)*(h)/G*prop.t_sk
            qb[i-1] = qb[i-1]*ds_t[i-1]
            qb1 += qb[i-1]*ds_t[i-1]

        if(i == 2 or i == 5):
            ds_t[i-1] = (h)/G*prop.t_sp
            qb[i-1] = qb[i-1]*ds_t[i-1]
            ds1 += (h)/G*prop.t_sp
            ds2 += (h)/G*prop.t_sp
            qb1 += qb[i-1]*ds_t[i-1]
            qb2 += qb[i-1]*ds_t[i-1]
            
        if(i == 3 or i == 4):
            ds_t[i-1] = (lsk)/G*prop.t_sk
            qb[i-1] = qb[i-1]*ds_t[i-1]
            qb2 +=  qb[i-1]*ds_t[i-1]
            ds2 += (lsk)/G*prop.t_sk 

    return (qb1, ds1, qb2, ds2)
    
def get_sc():
    '''
    This function calculates the sc from of z from center of the spar
    It basically uses the shear equation, lht is integral of (p*q_b*ds)
    rht is summation of (2*a*qs0) with respect to areas
    # '''
    q01, q02 = get_qs0()
    # rht = 2 * area1 * q01 + 2 * area2 * q02
    # qbooms = get_qbooms()
    # delta = 0.01
    # theta = np.arange(0, np.pi/2, delta)
    # theta1 = np.arange(0, np.pi/2, delta)
    # qbo = 0
    # for i in range(len(qbooms)):
    #     qbo = qbo + (prop.B_z[i] * qbooms[i])    
    # lht =[sum((h*get_qsec1() * delta * h * dtheta for dtheta in theta)),
    #       sum((h*get_qsec6() * delta * h * dtheta for dtheta in theta1)),
    #       get_qsec2()*h*prop.t_sp*0.5*prop.ha,
    #       get_qsec5()*h*prop.t_sp*0.5*prop.ha,
    #       get_qsec3()*lsk*prop.t_sk*(prop.Ca - h)*h/lsk,
    #       get_qsec4()*lsk*prop.t_sk*(prop.Ca - h)*h/lsk,
    #       ]
    # return((sum(lht)+ rht),0)

    q1_shear_tot = (get_qsec1() + q01)*(np.pi*h*0.5)*h
    q3_shear_tot = (get_qsec3() + q02 )*((prop.Ca-h)/2 + h)*lsk
    q4_shear_tot = (get_qsec4()+ q02 )*((prop.Ca-h)/2 + h)*lsk
    q6_shear_tot = (get_qsec6()+ q01)*(np.pi*h*0.5)*h

    Moment = q1_shear_tot + q3_shear_tot + q4_shear_tot + q6_shear_tot
    
    return(Moment-h, 0)
    
print("Shear center(y,z)", get_sc())


