# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 11:44:57 2019
@author: Pontus Nordqvist
"""
from  scipy import *
from  pylab import *
from CubicSpline import CubicSpline


import  unittest
import random

class Test_CubicSpline(unittest.TestCase):
    
    def setUp(self):
        #Might be cleaned up using txt document
        self.CONTROL = [(-12.73564, 9.03455),
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
        
        self.KNOTS = linspace(0, 1, 26)
        self.KNOTS[ 1] = self.KNOTS[ 2] = self.KNOTS[ 0]
        self.KNOTS[-3] = self.KNOTS[-2] = self.KNOTS[-1]
        
        cs = CubicSpline(self.KNOTS, self.CONTROL)
    
    def tearDown(self):
        CONTROL = None
        KNOTS = None
    
    def test_exampleTest(self):

        result = cs.point_eval(0.2)
        expected = array([-31.90219167, 6.47655833])

        self.assertAlmostEqual(result[0],expected[0])
        self.assertAlmostEqual(result[1],expected[1])

    def test_order(self):
        random.shuffle(self.KNOTS)
        cs = CubicSpline(self.KNOTS, self.CONTROL)
        result = cs.point_eval(0.2)
        expected = array([-31.90219167, 6.47655833])
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])
    
    """
    def test_dimensions(self):
        with self.assertRaises(TypeError):
            cp_1D = [i[0] for i in self.CONTROL]
            cs = CubicSpline(self.KNOTS, self.CONTROL)
    """
    
    def test_null(self):
        self.tearDown()
        cs = CubicSpline(self.KNOTS, self.CONTROL)
        with self.assertRaises(TypeError):
            cs.point_eval(0.2)
    
    """
    Uppg 4
    
    def test_blossom(self):
        #s = CubicSpline(cp, knots = knots)
        cs.basis_function(KNOTS, 0.2)
        bases = [cs.basis_function(self.KNOTS, i) for i in range(len(self.KNOTS)-4)]
        bases_in_point = [N(0.2) for N in bases]
        result = sum(bases_in_point)
        expected = cs.point_eval(0.2)
        self.assertAlmostEqual(result[0], expected[0])
        self.assertAlmostEqual(result[1], expected[1])
    
    def test_padded(self):
        #s = CubicSpline(self.cp, self.knots)
        checker_left = checker_right = False
        if (cs.knots[0] == cs.knots[1] and cs.knots[1] == cs.knots[2]):
            checker_left = True
        if (cs.knots[-1] == cs.knots[-2] and cs.knots[-2] == cs.knots[-3]):
            checker_right = True        
        self.assertTrue(checker_left)
        self.assertTrue(checker_right)
    
    def test_basis_sum(self):
        #s = CubicSpline(cp, knots = knots)
        bases = [cs.basis_function(self.KNOTS, i) for i in range(len(self.KNOTS)-4)]
        x = linspace(0,1, 200)
        y = [[N(val) for val in x] for N in bases]
        Y = zeros((1, len(x)))
        for element in y:
            Y += element
        #self.assertEqual([y for y in Y], 1)
        one = ones(len(Y))
        for i in range(0,len(Y)):
            self.assertAlmostEqual(Y[i],one[i])
    """
        
        
    
    def test_basisFunctionPositive(self):
        positiveValue = True
        
        bases = [c.basis_function(KNOTS,i) for i in range(len(KNOTS))]
        
        for i in range(0,len(KNOTS)-4):
            for j in linspace(0,1,200):

                positiveValue = bases[i](j) > 0
                
                if bases[i](j) == -0.0:
                    positiveValue = True
                    
                if positiveValue == False:
                    #print(bases[i](j))
                    break
                
        self.assertTrue(positiveValue)
        
if __name__ == '__main__':
    unittest.main()

