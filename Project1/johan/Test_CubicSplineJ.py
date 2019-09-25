# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 12:56:35 2019
@author: pontusnordqvist
"""
from  scipy import *
from  pylab import *


import unittest
from cubicspline import CubicSpline

import random

#__unittest = True

class TestCubicSpline(unittest.TestCase):
    def setUp(self):
        self.cp = [(-12.73564, 9.03455),
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
        self.knots = linspace(0, 1, 26)
    
    def test_almostEqual(self):
        s = CubicSpline(cp, knots = knots)
        result = s(0.2)
        expected = array([-31.90219167, 6.47655833])
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])
    
    def test_order(self):
        random.shuffle(knots)
        s = CubicSpline(cp, knots = knots)
        result = s(0.2)
        expected = array([-31.90219167, 6.47655833])
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])
    
    def test_dimensions(self):
        with self.assertRaises(TypeError):
            cp_1D = [i[0] for i in cp]
            s = CubicSpline(cp_1D, knots = knots)
            
    
    def test_null(self):
        self.tearDown()
        s = CubicSpline(cp, knots = knots)
        with self.assertRaises(TypeError):
            s(0.2)
            
    def test_blossom(self):
        s = CubicSpline(cp, knots = knots)
        s.basis_function(knots, 0.2)
        bases = [s.basis_function(self.knots, i) for i in range(len(self.knots))]
        bases_in_point = [N(0.2) for N in bases]
        result = sum(bases_in_point)
        expected = s(0.2)
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])
        
    def test_basis_sum(self):
        s = CubicSpline(cp, knots = knots)
        bases = [s.basis_function(self.knots, i) for i in range(len(self.knots))]
        x = linspace(0,1, 200)
        y = [[N(val) for val in x] for N in bases]
        Y = zeros(1, len(x))
        for element in y:
            Y += element
        self.assertEqual([y for y in Y], 1)
        pass
    
    def test_padded(self):
        s = CubicSpline(self.cp, self.knots)
        checker_left = checker_right = False
        if (s.knots[0] == s.knots[1] and s.knots[1] == s.knots[2]):
            checker_left = True
        if (s.knots[-1] == s.knots[-2] and s.knots[-2] == s.knots[-3]):
            checker_right = True        
        self.assertTrue(checker_left)
        self.assertTrue(checker_right)
        
        
    def tearDown(self):
        cp = None
        knots = None
        s = None
    

            
        
if __name__=='__main__':
    unittest.main()