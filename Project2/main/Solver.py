
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 21:12:56 2019
@author: Mattias Lundström1
"""
from  scipy import *
from  pylab import *
import OptimizationProblem.py
    
class Solver:
    
    def __init__(self, problem):
        self.problem = problem
        self.function = problem.objective_function
        self.n = problem.x_0.shape[0]
        self.delta_values = diag([1e-8 for i in range(n)])
       
        
    
    def __call__(self):
        pass

    def newton(self, mode='default', tol = 1e-08, maxIteration = 15):
        iterations = 0
        x_k = self.problem.x_0
        x_next = self._newton_step(x_k)
        
        while(norm(x_k-x_next)>tol and iterations < maxIteration):
            x_k = x_next
            x_next = self._newton_step(x_k, mode)
            iterations +=1
        if norm(x_k-x_next) > tol:
            pass
            #raise exception didn't converge 
        
        return x_next
    
    def _newton_step(self, x_k, mode):
        if mode == 'default':
            alpha = 1
        
        elif mode == 'exact':
            alpha = self.exact_line_search(x_k)
            
        elif mode == 'inexact':
            alpha = self.inexact_line_search(x_k)
            
        if problem.gradient is not None:
            gradient = self.problem.gradient(x_k)
        else:
            gradient = self._gradient(x_k)
        hessian = self._hessian(x_k)
        
        x_next = solve(hessian, hessian@x_k - gradient)
        
        return x_next
    
    def _gradient(self, x_k):
        gradient = array([(self.function(x_k)
                            -self.function(x_k-self.delta_values[:,i]))*1e8
                            for i in range(n)])
        return gradient
    
    def _search_dir(self, x_k):
        pass    
    
    def exact_line_search(self):
        pass
    
    def inexact_line_search(self):
        pass
    
    def _hessian(self, x_k):

        
    
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
    
    
        
    
    
