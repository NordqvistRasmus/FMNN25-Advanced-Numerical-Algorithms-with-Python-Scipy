#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:42:55 2019

@author: Arvid
"""
from scipy import *
from matplotlib.pyplot import *

class CubicSpline:
    
    
    def __init__(self, grid_u, control_d):
        self.grid_u = grid_u
        self.control_d = control_d
        
    def __call__(self):
        pass
    
    def plot(self, plot_poly = False):
        grid = linspace(self.grid_u[0], self.grid_u[-1],10)
        points = array([self.point_eval(point) for point in grid])
        plot(points[:,0],points[:,1])
        if plot_poly:
            plot(self.control_d, 'bx')
    
    def point_eval(self, x):
#        if u < self.grid_u[0] or u > self.grid_u[-1]:
#            raise ValueError(f"""u = {u} is not contained in
#                             the grid [{self.grid_u[0]}, {self.grid_u[-1]}]
#                             """)
        indx = (x > self.grid_u).argmax()
        
        blossoms = [self.control_d[i] for i in range(indx-2, indx+2)]
        u = [self.grid_u[i] for i in range(indx-2, indx+4)]
        
        alphas = [((u[i+3] - x)/(u[i+3] - u[i])) for i in range(3)]
        blossoms = [alphas[i]*blossoms[i] + (1 - alphas[i])*blossoms[i+1]
                    for i in range(3)]
        
        alphas = [(u[i+3] - x)/(u[i+3] - u[i+1]) for i in range(2)]
        blossoms = [alphas[i]*blossoms[i] + (1 - alphas[i])*blossoms[i+1]
                    for i in range(2)]
        
        alphas = (u[3]-x)/(u[3]-u[2])
        blossoms = alphas*blossoms[0] + (1-alphas)*blossoms[1]
        
        return blossoms
        
        
        
        
    
        
        
        
    
        
        
        


