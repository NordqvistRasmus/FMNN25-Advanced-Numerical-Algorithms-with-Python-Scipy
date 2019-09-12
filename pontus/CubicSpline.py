import scipy.linalg as sl
import numpy as np

# -*- coding: utf-8 -*-

class CubicSpline:
    def __init__(self, u, d):
        self.x=x
        self.y=y
        
        uhot=hotIntervall(u)
    
    #def __call__
    
    #Maybe dosent work
    def hotIntervall(u):
        uHot = np.where(u != 0)[0]
        return uHot
    
    def spanIntervall(u):
        #Finds inceses of where array is nonzero
        #uNonZeroIndeces = np.where(u != 0)[0]
        uSpan=[l for l in u if l!=0]
        return uSpan
    
    def alpha(u):
        uSpan=spanIntervall(u)
        alpha=(u[-1]-u)/(u[-1]-u[0])
        return alpha
    #def __call__

    #def plot:
    
    print("test")
