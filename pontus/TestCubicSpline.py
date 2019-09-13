# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 11:44:57 2019
@author: Pontus Nordqvist
"""
from  scipy import *
from  pylab import *
from CubicSpline import CubicSpline


import  unittest

class TestCubicSpline(unittest.TestCase):
    
    #def setUP(self):
    
    """    
    def exampleTest(self):
        #Can be cleaned up using txt file
        
        CONTROL = [(-12.73564, 9.03455),(-26.77725, 15.89208),(-42.12487, 20.57261),
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
        
        CONTROL = [(0.,0.), (0.,0.)]
        KNOTS = linspace(0, 1, 26)
        KNOTS[ 1] = KNOTS[ 2] = KNOTS[ 0]
        KNOTS[-3] = KNOTS[-2] = KNOTS[-1]
        
        c=CubicSpline(KNOTS, CONTROL)
        result = c.point_eval(0.2)
        expected = array([-31.90219167, 6.47655833]) 
        #expected = array([-0.8, 6])
        
        self.assertEqual(result , expected)
        """
    #Test method to test unittest
    def simpleMath(self):
        self.assertEqual(5+9, 14)
        
if __name__ == '__main__':
    unittest.main()
#if  __name__ == '__main__ ':
#    unittest.main()
