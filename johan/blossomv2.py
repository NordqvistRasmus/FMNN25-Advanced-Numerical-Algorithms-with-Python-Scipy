# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 12:40:55 2019
@author: pontusnordqvist
"""
from  scipy import *
from  pylab import *

"""
    Method for the blossom algorithm
    param:  t, given point around which we want to creata a spline
            d, the control points.
            u, the hot interval, knot vector.
            
    returns:    d[0], i.e s(u) if the size of the hot interval are 1
                which occurs when the algorithm is complete
                
                else returns the same function, blossomv2, in a 
                recursive manner to again do the calculations.
    """
    def blossomv2(self, t, d, u):
        if u.size == 1:
            return d[0]
        else:
            mid = u.size//2
            left = u[:mid]
            right = u[mid:]
            new_d = zeros(d.size-1)
            alpha = zeros(d.size-1)
            for i in range (alpha.size):
                alpha[i] = (right[i] - t)/(right[i] - left[i])
            for j in range (new_d.size):
                new_d[j] = alpha[j]*d[j] + (1-alpha[j])*d[j-1]
            return self.blossomv2(t, new_d, u[1:-1])