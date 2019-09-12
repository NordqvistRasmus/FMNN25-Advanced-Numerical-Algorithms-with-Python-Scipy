import scipy.linalg as sl
import numpy as np

# -*- coding: utf-8 -*-

class CubicSpline:
    def __init__(self, x, y):
        self.x=x
        self.y=y
    
    def hotIntervall(u):
        uHot = np.where(u != 0)[0]
        return uHot
    
    #def __call__

    #def plot:
    
    print("test")
