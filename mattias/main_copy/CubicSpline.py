#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:42:55 2019

@author: Arvid
"""
from scipy import *
from matplotlib.pyplot import *

class CubicSpline:
    
    """
    A class used to represent a CubicSpline with the purpose of curve design in 2D. 

    ...

    Attributes
    ----------
    KNOTS : float
        Sequence of floats
    CONTROL : List of tuples
        Control points, used for shaping the curve.

    Methods
    -------
    __init__(knots, control):
        Sorts the knot grid and adds padding. 
    
    plot(plot_poly = True, precision = 150):
        Calculates and plot points using De Boors algorithm. 
        If 'plot_poly' is set True (default), then also plots the Control Points.  
        To increase plot resolution, increase 'precision' parameter.
        
    _point_eval(u):
        Helph method to 'plot', used to evaluate point 'u' with the De Boors algorithm.
        
    basis_fuction(knots, i):
        Takes a knot sequence 'knots' and index 'i' as input and returns a 
        python function to evaluate the i:th B-spline basis function. 
    """
    
    
    def __init__(self, knots, control):
        self.knots = knots
        self.knots.sort()
        self.knots = r_[knots[0], knots[0], knots, knots[-1], knots[-1]]
        self.control = array(control)
        
    def __call__(self):
        pass
    
    def plot(self, plot_poly = True, precision = 150):
        
        points = array([self.point_eval(point) for point in linspace(0, 1, precision)])
        #print(points)
        plot(points[:,0],points[:,1])
        if plot_poly:
            plot(self.control[:,0],self.control[:,1], 'yx--')
    
    def point_eval(self, u):
        indx = int(self.knots.searchsorted([u]))-1
        blossoms = array([self.control[i] for i in range(indx-2, indx+2)])
        knots = array([self.knots[i] for i in range(indx-2, indx+4)])
        
        alphas = array([((knots[i+3] - u)/(knots[i+3] - knots[i])) for i in range(3)])
        blossoms = array([alphas[i]*blossoms[i] + (1 - alphas[i])*blossoms[i+1]
                    for i in range(3)])
        
        alphas = array([(knots[i+3] - u)/(knots[i+3] - knots[i+1]) for i in range(2)])
        blossoms = array([alphas[i]*blossoms[i] + (1 - alphas[i])*blossoms[i+1]
                    for i in range(2)])
        
        alphas = (knots[3]-u)/(knots[3]-knots[2])
        blossoms = alphas*blossoms[0] + (1-alphas)*blossoms[1]
        
        return blossoms

    
    def basis_function(self, knots, i):
        
        def basis(u, i, k):
            if k == 0:
                if knots[i-1] == knots[i]:
                    return 0
                elif knots[i-1] <= u < knots[i]:
                    return 1
                else:
                    return 0
            else:
                if knots[i+k-1] == knots[i-1]: # and u == knots[i-1]:
                    coeff1 = 0
                else:
                    coeff1 = (u - knots[i-1])/(knots[i+k-1]-knots[i-1])
                if knots[i+k] == knots[i]: # and u == knots[i+k]:
                    coeff2 = 0
                else:
                    coeff2 = (knots[i+k] - u)/(knots[i+k] - knots[i])
                
                return basis(u, i, k-1)*coeff1 + basis(u, i+1, k-1)*coeff2
        def j_basis(u):
            return basis(u, i, 3)
            
        return j_basis


if __name__ == '__main__':
    
    CONTROL = [(-12.73564, 9.03455),
                   (-26.77725, 15.89208),
                   (-42.12487, 20.57261),
                   (-15.34799, 4.57169),
                   (-31.72987, 6.85753),
                   (-49.14568, 6.85754),
                   (-38.09753, -1e-05),
                   (-67.92234, -11.10268),
                   (-89.47453, -33.30804),
                   (-21.44344, -22.31416),
                   (-32.16513, -53.33632),
                   (-32.16511, -93.06657),
                   (-2e-05, -39.83887),
                   (10.72167, -70.86103),
                   (32.16511, -93.06658),
                   (21.55219, -22.31397),
                   (51.377, -33.47106),
                   (89.47453, -33.47131),
                   (15.89191, 0.00025),
                   (30.9676, 1.95954),
                   (45.22709, 5.87789),
                   (14.36797, 3.91883),
                   (27.59321, 9.68786),
                   (39.67575, 17.30712),
                   (39.67575, 17.30712)] 
        
    KNOTS = linspace(0, 1, 21)
    print(KNOTS)
    c = CubicSpline(KNOTS, CONTROL)
    #print(c.basis_function(KNOTS,0))
    print(c.point_eval(0.2))  
    c.plot()
    print("Debug")
        
    
        
        
        
    
        
        
        


