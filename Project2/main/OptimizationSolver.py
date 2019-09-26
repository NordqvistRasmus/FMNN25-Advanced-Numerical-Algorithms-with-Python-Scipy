
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 21:12:56 2019
@author: Mattias Lundström1
"""
from  scipy import *
from  pylab import *
import OptimizationProblem.py

"""
Abstract solver class
"""
class OptimizationSolver():
    
    #def __init__(self, obj_function, x_0, gradient = None, hessian = None):
        
        #pass
    
   
class Solver(OptimizationSolver):
    
    def __init__(self):
        pass
    
    def __call__(self):
        pass

    def newton(self):
        pass
    
    def _newton_step(self):
        pass
    
    def exact_line_search(self):
        pass
    
    def inexact_line_search(self):
        pass
    
    def _hessian(self):

        
    
class GoodBroydenSolver(Solver):
    
    def _hessian(self, hessian):
        pass
        
class BadBroydenSolver(Solver):
    
    def _hessian(self, hessian):
        pass
        
class DFP2Solver(Solver):
    
    def _hessian(self, hessian):
        pass
    
class BFGS2Solver(Solver):
    
    def _hessian(self, hesssian):
        pass
    
    
        
    
    
