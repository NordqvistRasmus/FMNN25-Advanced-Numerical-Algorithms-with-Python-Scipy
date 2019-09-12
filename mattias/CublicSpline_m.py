# -*- coding: utf-8 -*-
"""
Created on %(date)s
@author: %(username)s
"""
from scipy import *
from pylab import *
from matplotlib.pyplot import*

class CubicSpline:
    def __init__(self, u_vector = None):
        if u_vector == None:
            self.u_vector = linspace(0, 1.0, 9)
        print("Setup input vector", self.u_vector)
    
    def evaluateAt(self, u, I):
        return 99
    
    
    def __call__(self, u_vector):
            u_vector.sort()
            
            result = empty(len(u_vector))
        
            for i,u in enumerate(u_vector):
                I = (u < u_vector).argmax()
                print(I)
                result[i] = self.evaluateAt(u, I)
            

    def plot():
        return 0

input = array([0.1, 0.3, 0.6, 0.9])
Spline = CubicSpline()
res = Spline(input)




# @index, 
""" 

u_vector = array([1., 2. ,4. ,5.5 ,8 ])
N = array([])

def createBasisFunction(j, u): 
    return 0    

def getBasisFirstDegree(j, u):
        if u_vector[j-1] == u_vector[j]:
            return 0
        elif (u == u_vector[j-1] or u == u_vector[j]):
            return 1
        else:
            return 0

"""  


