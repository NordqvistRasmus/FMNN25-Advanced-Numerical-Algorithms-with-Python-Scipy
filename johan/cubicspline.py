#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 11:38:09 2019

@author: johanliljegren
"""
from  scipy import *
from  pylab import *

class CubicSpline:
    def __init__(self, cp, p = None, knots = None):
        self.d = array(cp)
        if p is None:
            self.p = 3
        if knots is None:
            knots = linspace(0, 10, len(cp)-2)
        else:
            knots.sort()
        self.knots = r_[knots[0], knots[0], knots, knots[-1], knots[-1]]    
        
        
        
    def __call__(self, t):
        index = self.hot_interval(t)
        control_points = self.d[(index-3):(index+1)]
        interval = self.knots[index-3:index+3]
        s = self.blossomv2(t, control_points, interval)
        return s
    """
    Method for the blossom algorithm
    param:  t, given point around which we want to creata a spline
            d, the control points.
            dh, the hot interval, knot vector.
            
    returns:    d[0], i.e s(u) if the size of the hot interval are 1
                which occurs when the algorithm is complete
                
                else returns the same function, blossomv2, in a 
                recursive manner to again do the calculations.
    """
    def blossomv2(self, t, d, dh):
        if dh.size == 0:
            return d[0]
        else:
            mid = shape(dh)[0]//2
            left = dh[:mid]
            right = dh[mid:]
            alpha = zeros((shape(d)[0]-1,2))
            for i in range (shape(alpha)[0]):
                alpha[i] = (right[i] - t)/(right[i] - left[i])
            new_d = zeros((shape(d)[0]-1,2))
            for j in range (shape(new_d)[0]):
                new_d[j] = alpha[j]*d[j] + (1-alpha[j])*d[j-1]
            return self.blossomv2(t, new_d, dh[1:-1])
    
    """
    Finds the "hot interval" given the specific point.
    Then extract this part of the grid and returns it
    param: t, given point around which we want to creata a spline
    returns: interval, the hot interval, a point on the grid u
    """
    def hot_interval(self, t): 
        idx = self.knots.searchsorted([t])-1
        index = idx[0]
        return index
    
    def plot(self, control_poly = True):
        points = array([self.__call__(point) for point in linspace(0.1,0.8,150)])
        plot(points[:,0],points[:,1])
        if control_poly:
            plot(self.d[:,0],self.d[:,1], 'yx--')
        pass
    
    def interpolation(self):
        xi=(u[:-2]+u[1:-1]+u[2:])/3.
        
   #cp = array([[1, 3], [3.5, 8], [5, 4],[5.4, 5],[6.6, 2], [7, 4], [8, 3]])
