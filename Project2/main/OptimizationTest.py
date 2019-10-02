#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 17:55:09 2019

@author: Mattias Lundström, Pontus Nordqvist, Johan Liljegren, Arvid Rolander, Antonio Alas
"""
from  pylab import array
import numpy as np
from scipy.optimize import minimize, rosen, fmin_bfgs
from Solver import Solver, QuasiNewtonSolver, BFGS2Solver, DFP2Solver, BadBroydenSolver, GoodBroydenSolver
from OptimizationProblem import OptimizationProblem
from chebyquad_problem import chebyquad, gradchebyquad


import unittest


class Test_Optimization(unittest.TestCase):    
    
    """
    Test class which test Solver with OptimizationProblem
    """
   
    #Default problem is rosen and guess [230, 30]    
    def setUp(self, func = rosen, x_0 = array([230,30]), gradient = None):
        """
        Sets up the test class. 
        
        Parameters
        ----------
        
        func : Function to be tested, default is rosenbrock function
        x_0 : Guess, default guess is [230, 30] for the rosenbrock function
        gradient :  Optional parameter for gradient. Else it is solved. 
        """
        self.problem = OptimizationProblem(func, x_0, gradient)
        
      
    def tearDown(self):
        """
        Tears down the problem. 
        """
        self.problem = None
        
        
    def setGuess(self, guess):
        """
        Set a new guess for the test function
        """
        self.problem.x_0 = guess
        
    def setFunction(self, function):
        """
        Set a new function for the test function
        """
        self.problem.objective_function = function
        
    def test_newton(self):
        """
        Test rosenbrock function
        """
        s = Solver(self.problem)
        result = s.newton(mode = 'default')
        expected = array([1., 1.])
        self.assertAlmostEqual(0, norm(result-expected))
        
        
    def test_newton_rosen_newzeros(self):
        """
        Test rosenbrock with new zeros. New zeros is [a, a^2]
        """
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
        """
        Test Exact line search 
        """
        s = Solver(self.problem)
        result = s.newton(mode = 'exact', maxIteration = 400)
        expected = array([1., 1.])
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])   
        
        
    def test_good_broyden(self):
        """
        Test GoodBroyden with simple function
        """
        f = lambda x: (x[0]-10)**2 + (x[1]-10)**2
        self.setUp(f, array([20, 5]))
        
        s = GoodBroydenSolver(self.problem) 
        result = s.newton(mode = 'inexact', maxIteration = 400)
        
        expected = array([10, 10])
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])   
        
    def test_bad_broyden(self):
        """
        Test BadBroyden with simple function
        """
        f = lambda x: (x[0]-10)**2 + (x[1]-10)**2
        self.setUp(f, array([20, 5]))
        
        s = BadBroydenSolver(self.problem) 
        result = s.newton(mode = 'inexact', maxIteration = 400)
        
        expected = array([10, 10])
        self.assertAlmostEqual(result[0], expected[0], delta = 10e-4)
        self.assertAlmostEqual(result[1], expected[1], delta = 10e-4)
        
    def test_DFP2Solver(self):
        """
        Test DFP2 with simple function
        """
        f = lambda x: (x[0]-10)**2 + (x[1]-10)**2
        self.setUp(f, array([20, 5]))
        
        s = DFP2Solver(self.problem) 
        result = s.newton(mode = 'inexact', maxIteration = 400)
        
        expected = array([10, 10])
        self.assertAlmostEqual(result[0], expected[0], delta = 10e-4)
        self.assertAlmostEqual(result[1], expected[1], delta = 10e-4)
        
    def test_BFGS2Solver(self):
        """
        Test BFGS2 with simple function
        """
        f = lambda x: (x[0]-10)**2 + (x[1]-10)**2
        self.setUp(f, array([20, 5]))
        
        s = BFGS2Solver(self.problem) 
        result = s.newton(mode = 'inexact', maxIteration = 400)
        
        expected = array([10, 10])
        self.assertAlmostEqual(result[0], expected[0], delta = 10e-4)
        self.assertAlmostEqual(result[1], expected[1], delta = 10e-4)
      
    def test_chebyquad_4(self): #Does not currently work, issues with too small steps.
        """
        Test chebyquad function with inexact line serach (Bad broyden)
        """
        
        guess = array([0.25, 0.5, 0.2, 0.9])
        self.setUp(chebyquad, guess, gradchebyquad)
        
        s = BadBroydenSolver(self.problem)
        
        result = s.newton(mode = 'inexact', maxIteration = 1000)  
        x = linspace(0,1, 4)
       
        xmin = fmin_bfgs(chebyquad,x,gradchebyquad)  
        print('result is:', result)
        print('xmin is:', xmin)
        self.assertAlmostEqual(0, norm(result-xmin), delta = 10e-3)
        
   
    def test_sanity(self):
         """
         Simple sanity test. Test the test class
         """
         self.tearDown()
         with self.assertRaises(AttributeError):
             s = Solver(self.problem)
             
           
        
if __name__ == '__main__':
    unittest.main()
    