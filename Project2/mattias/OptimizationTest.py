#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 17:55:09 2019

@author: johanliljegren
"""
from  scipy import *
from  pylab import *
import numpy as np
from scipy.optimize import minimize, rosen
from Solver import Solver
from OptimizationProblem import OptimizationProblem

import unittest


#class exact_line(TestCase):
    
    #def setUp(self):
    #   Någonstans här ska man väl inputa en gissning?
        #prob = OptimizationProblem(rosen)
        #s = Solver(prob)
        #res = s.newton(mode='exact', tol = 1e-6)
        #expected = [1,1]
        #self.assertAlmostEqual(res, expected, delta = 1e-6)


class Test_Optimization(unittest.TestCase):    
    
    #Default problem is rosen and guess [230, 30]    
    def setUp(self):
        self.problem = OptimizationProblem(rosen, array([230,30]))
        
    def setUpRosen(self):
        pass
            
            
    def tearDown(self):
        self.problem = None
 
    def SetGuess(self, x, y):
        self.problem.x_0 = array([x, y])
    
    def test_newton(self):
        s = Solver(self.problem)
        result = s.newton(mode = 'default')
        expected = array([1., 1.])
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])
        
    def test_newton_converge_fail(self):
        self.SetGuess(699, 699**2)
        s = Solver(self.problem)
        result = s.newton(mode = 'default')
        expected = array([1., 1.])
        close = np.isclose(result, expected) #Check that the result is 'Not almost equal' the expected value
        self.assertFalse(close[0])
        self.assertFalse(close[1])

    def test_exact_line_search(self):
        s = Solver(self.problem)
        result = s.newton(mode = 'exact')
        expected = array([1., 1.])
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])     
      
    """
    def test_exact_line_search(self):
        s = Solver(self.problem)
        result = s.newton(mode = 'inexact')
        expected = array([1., 1.])
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])  
      """
        
    def test_sanity(self):
        expected = 1
        result = 1
        self.assertAlmostEqual(result, expected)
        
if __name__ == '__main__':
    unittest.main()
    

        
    
        
    