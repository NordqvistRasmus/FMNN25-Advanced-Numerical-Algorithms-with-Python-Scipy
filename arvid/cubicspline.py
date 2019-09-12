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
        
    def __call__(self,):
        pass
    
    def plot(self, plot_poly = False):
        pass
    
    def point_eval(self, u):
        if u < self.grid_u[0] or u > self.grid_u[-1]:
            raise ValueError(f"""u = {u} is not contained in
                             the grid [{self.grid_u[0]}, {self.grid_u[-1]}]
                             """)
        indx = (u > self.grid_u).argmax()
        ctrl_points = array([self.control_d[:,i] for i in range(indx-1, indx+2)])
        


