#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 17:55:09 2019

@author: johanliljegren
"""
from  pylab import array
import numpy as np
from scipy.optimize import minimize, rosen, fmin_bfgs
from Solver import Solver, QuasiNewtonSolver, BFGS2Solver, DFP2Solver, BadBroydenSolver, GoodBroydenSolver
from OptimizationProblem import OptimizationProblem
from chebyquad_problem import chebyquad, gradchebyquad


import unittest


class Test_Optimization(unittest.TestCase):    
    
    #Default problem is rosen and guess [230, 30]    
    def setUp(self, func = rosen, x_0 = array([230,30]), gradient = None):
        self.problem = OptimizationProblem(func, x_0, gradient)
        
    def tearDown(self):
        self.problem = None
 
    def setGuess(self, guess):
        self.problem.x_0 = guess
        
    def setFunction(self, function):
        self.problem.objective_function = function
        
    def test_newton(self):
        s = Solver(self.problem)
        result = s.newton(mode = 'default')
        expected = array([1., 1.])
        self.assertAlmostEqual(0, norm(result-expected))
        
        """
        Set new zeros. New zeros is [a, a^2]
        """
    def test_newton_rosen_newzeros(self):
        self.setUp()
        a = 7
        f = lambda x: 100*((x[1]-x[0]**2)**2)+(a-x[0])**2
        self.setFunction(f)
        self.setGuess(array([10, 31]))
        s = Solver(self.problem)
        result = s.newton(mode = 'default', maxIteration = 100)
        expected = array([a, a**2])
        self.assertAlmostEqual(0, norm(result-expected))
 
        """
        Tests convergence fail with a bad guess and low maxIterations
        """    
    def test_newton_converge_fail(self):
        self.setGuess(array([699, 30300]))
        s = Solver(self.problem)
        result = s.newton(mode = 'default', maxIteration = 50)
        expected = array([1., 1.])
        close = np.isclose(result, expected) #Check that the result is 'Not almost equal' the expected value
        self.assertFalse(close[0])
        self.assertFalse(close[1])

    def test_exact_line_search(self):
        s = Solver(self.problem)
        result = s.newton(mode = 'exact', maxIteration = 400)
        expected = array([1., 1.])
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])     
      
    def test_chebyquad_4(self):
        guess = array([0.25, 0.5, 0.2, 0.9])
        self.setUp(chebyquad, guess, gradchebyquad)
        
        s = BadBroydenSolver(self.problem)
        
        result = s.newton(mode = 'inexact', maxIteration = 1000)  #HERE WE SHPULD HAVE OTHER MODE OR TEST BFGS, REGULAR DOES NOT WORK
        x = linspace(0,1, 4)
       
        xmin = fmin_bfgs(chebyquad,x,gradchebyquad)  
        print('result is:', result)
        print('xmin is:', xmin)
        self.assertAlmostEqual(0, norm(result-xmin), delta = 10e-3)
        
    """
    def test_exact_line_search(self):
        s = Solver(self.problem)
        result = s.newton(mode = 'inexact')
        expected = array([1., 1.])
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])  
     """
     
    def test_sanity(self):
        self.tearDown()
        with self.assertRaises(AttributeError):
            s = Solver(self.problem)
           
        
if __name__ == '__main__':
    unittest.main()
    