#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:42:55 2019

@author: Arvid
"""
from scipy import *
from matplotlib.pyplot import *


class CubicSpline:
    
    
    def __init__(self, knots, control = None, interpol_points = None,
                 interpolate = False):
        #self.knots = r_[knots[0], knots[0], knots, knots[-1], knots[-1]]
        self.knots = knots
        self.knots.sort()
        control = array(control)
        self.control = control
        
    def __call__(self, knots, control = None, interpol_points = None):
        pass
    
    def plot(self, plot_poly = True):
        
        points = array([self.point_eval(point) for point in linspace(0,1,150)])
        print(points)
        plot(points[:,0],points[:,1])
        if plot_poly:
            plot(self.control[:,0],self.control[:,1], 'yx--')
    
    def point_eval(self, u):
#        if u < self.grid_u[0] or u > self.grid_u[-1]:
#            raise ValueError(f"""u = {u} is not contained in
#                             the grid [{self.grid_u[0]}, {self.grid_u[-1]}]
#                             """)
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
         # Because of padding we want to remove the first and last 2 elements in knots
        
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
#        def j_basis(u):
#            return basis(u, i, 3)   # Write as lambda instead
#                                    # lambda u: basis(u,i,3)
        return lambda u: basis(u, i, 3)
    
    def interpolate(self, interpol_points):
        x_points = interpol_points[:,0]
        y_points = interpol_points[:,1]
        bases = [self.basis_function(self.knots, i)
                 for i in range(len(self.knots - 2))]
        


#if __name__ == '__main__':
#    c = CubicSpline(KNOTS, CONTROL)
#    print(c.point_eval(0.2))  
#    c.plot()
#        
        
if __name__ == '__main__':
    
    KNOTS = array([0.  , 0.  , 0.  , 0.12, 0.16, 0.2 , 0.24, 0.28, 0.32, 0.36, 0.4 ,
       0.44, 0.48, 0.52, 0.56, 0.6 , 0.64, 0.68, 0.72, 0.76, 0.8 , 0.84,
       0.88, 1.  , 1.  , 1.  ])
            
    c = CubicSpline(KNOTS)
    
    basis = [c.basis_function(KNOTS, i) for i in range(len(KNOTS)-2)]
    X = linspace(0, 1, 200)
    Y = array([[N(x) for x in X] for N in basis])
        
    
        
        
        
    
        
        
        


