#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 18:48:12 2020

@author: daanwitte
"""
import numpy as np
from math import *

T = 1
t_sk = 0.0011
t_sp = 2.4e-3
h_a = 16.1e-2
r = h_a/2
c_a = 0.505
len_sk = len_sk = sqrt((c_a-0.5*h_a)**2 + (0.5*h_a)**2)
A1 = 0.5*pi*r**2
A2 = (c_a-h_a/2)*h_a/2

A = np.array([[(h_a/(2*A1*t_sp))+(pi*r)/(2*A1*t_sk),-h_a/(2*A1*t_sp),-1],
              [-h_a/(2*A2*t_sp),(h_a/(2*A2*t_sp))+((2*len_sk)/(2*A2*t_sk)),-1],
              [2*A1,2*A2,0]])

B = np.array([0,0,1])

x = np.linalg.solve(A,B)
print(x)

Gdsigmadz = x[2]
J = T/(dsigmadz)
print(J)