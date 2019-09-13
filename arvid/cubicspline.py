#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:42:55 2019

@author: Arvid
"""
from scipy import *
from matplotlib.pyplot import *

class CubicSpline:
    
    
    def __init__(self, knots, control):
        self.knots = knots
        self.knots.sort()
        control = array(control)
        self.control = control
        
    def __call__(self):
        pass
    
    def plot(self, plot_poly = True):
        grid = linspace(self.grid_u[0], self.grid_u[-1],10)
        points = array([self.point_eval(point) for point in grid])
        plot(points[:,0],points[:,1])
        if plot_poly:
            plot(self.control_d[:,0],self.control_d[:,1], 'yx--')
    
    def point_eval(self, u):
#        if u < self.grid_u[0] or u > self.grid_u[-1]:
#            raise ValueError(f"""u = {u} is not contained in
#                             the grid [{self.grid_u[0]}, {self.grid_u[-1]}]
#                             """)
        indx = int(self.knots.searchsorted([u]))
        print(indx)
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
        
KNOTS = array([0.  , 0.  , 0.  , 0.12, 0.16, 0.2 , 0.24, 0.28, 0.32, 0.36, 0.4 ,
       0.44, 0.48, 0.52, 0.56, 0.6 , 0.64, 0.68, 0.72, 0.76, 0.8 , 0.84,
       0.88, 1.  , 1.  , 1.  ])
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
 (39.67575, 17.30712)]

if __name__ == '__main__':
    c = CubicSpline(KNOTS, CONTROL)
    c.point_eval(0.2)  
        
        
    
        
        
        
    
        
        
        


