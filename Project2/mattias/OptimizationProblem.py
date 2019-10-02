#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 21:31:59 2019

@author: Mattias Lundström, Pontus Nordqvist, Johan Liljegren, Arvid Rolander, Antonio Alas
"""
from  scipy import *
from  pylab import *

class  OptimizationProblem():
    """
    Problem class for defining a problem.
    
    ...
    
    Attributes
    ----------
        objective_function: The problem function.
        x_0: The initial guess for the problem.
        gradient: The gradient of the problem. Defaults to None.
        hessian: The hessian of the problem. Defaults to None. 
    """
    
    def __init__ (self,objective_function, x_0, gradient = None, hessian = None):
        """
        Initializes a optimazation problem. Optional to give a predetermined 
        gradient or hessian, defults to none.
        """
        self.objective_function = objective_function
        self.x_0 = x_0
        self.gradient = gradient
        self.hessian = hessian
        
    
    def plot_problem(self):
        """
        Plot method for the the problem giver to show the problem. However,
        not implemented due to time constraints.
        """
        pass
    