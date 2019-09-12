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
    Input: Control points and kont grid (optional, otherwise default grid)
    """
    def __init__(self, cp, knots = None):
        self.cp = cp
        if knots is None:
            self.knots = linspace(0,10,len(cp) - 2) #default grid, two more points than the control points
            print("No input grid, set default:", self.knots)
        self.knots.sort()
        print("Knots:", self.knots)
        self.u_vector = self.addClamp(self.knots)
        print("Clamped knots:", self.knots)
        print("Control point vector:", cp)
    
    """
    Input: Point u to evaluate, I is hot index
    """
    

    
    def addClamp(self, knots):
        return r_[knots[0], knots[0], knots, knots[-1], knots[-1]]

    
    
    """
    Sort, add clamp (padding), and evaluate every u in u_vector
    num - number of u to evaluate, increases precision. 
    """
    def __call__(self, num): #grid vector as parameter?
        u_vector = linspace(self.knots[0], self.knots[-1], num)
        result = empty(len(u_vector)) #prep for input
    
        for i,u in enumerate(u_vector): #input count, element 
            I = (u < self.knots).argmax()
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




