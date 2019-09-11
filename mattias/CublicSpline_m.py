# -*- coding: utf-8 -*-
"""
Created on %(date)s
@author: %(username)s
"""
from scipy import *
from pylab import *

#class CubicSpline:
    #__init__
    #__call__

    #def plot:

# @index, 
    
u_vector = array([1., 2. ,4. ,5.5 ,8 ])


def createBasisFunction(j, u): 
    return 0
    

def getBasisFirstDegree(j, u):
        if u_vector[j-1] == u_vector[j]:
            return 0
        elif (u == u_vector[j-1] or u == u_vector[j]):
            return 1
        else:
            return 0

    
u = array([1., 2. ,4. ,5.5 ,8 ])
getBasisFirstDegree(1, 4)
    
    