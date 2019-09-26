#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 21:31:59 2019

@author: johanliljegren
"""
from  scipy import *
from  pylab import *

class  OptimizationProblem():
    
    def __init__ (self,objective_function, x_0, gradient = None, hessian = None):
        self.objective_function = objective_function
        self.x_0 = x_0
        self.gradient = gradient
        self.hessian = hessian
        
    
    def __call__(self):
        """
        Might be redundant
        """
        pass