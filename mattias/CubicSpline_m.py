# -*- coding: utf-8 -*-
"""
Created on %(date)s
@author: %(username)s
"""
from scipy import *
from pylab import *
from matplotlib import pyplot as plt

class CubicSpline:
    """
    Input: Control points and u grid (optional, otherwise default grid)
    """
    def __init__(self, cp, u_vector = None):
        self.cp = cp
        if u_vector is None:
            self.u_vector = linspace(0,10,len(cp) - 2) #default grid, two more points than the control points
            print("No input grid, set default:", self.u_vector)
        print("Control point vector:", cp)
    
    """
    Input: Point u to evaluate, I is hot index
    """
    

    
    def addClamp(self, u_vector):
        return r_[u_vector[0], u_vector[0], u_vector, u_vector[-1], u_vector[-1]]

    
    
    """
    Sort, add clamp (padding), and evaluate every u in u_vector
    """
    def __call__(self): #grid vector as parameter?
        self.u_vector.sort() 
        print("Grid:", self.u_vector)
        self.u_vector = self.addClamp(self.u_vector)
        print("Clamped:", self.u_vector)
        result = empty(len(self.u_vector)) #prep for input
    
        for i,u in enumerate(self.u_vector): #input count, element 
            I = (u < self.u_vector).argmax()
            print("Count:" , i, " Hot Index:" , I, "Value: ", u)
            result[i] = self.evaluateAt(u, I) #store result in resultvector
        return result
    
    """
    Here we run the algoritm for evaluating the spline
    """
    def evaluateAt(self, u, I): 
        return None
        
        

        
    def plot(self):
        plt.plot(self.cp[:,0], self.cp[:,1], '-.*') 
        plt.plot(self.u_vector, zeros(len(self.u_vector)), '.')
        plt.show()




