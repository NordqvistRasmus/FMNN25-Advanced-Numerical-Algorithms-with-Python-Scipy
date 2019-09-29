#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 17:55:09 2019

@author: johanliljegren
"""
from  scipy import *
from  pylab import *
from scipy.optimize import minimize, rosen
from OptimizationProblem import OptimizationProblem
from Solver import Solver

import unittest


class exact_line(TestCase):
    def setUp(self):
        #Någonstans här ska man väl inputa en gissning?
        prob = OptimizationProblem(rosen)
        s = Solver(prob)
        res = s.newton(mode='exact', tol = 1e-6)
        expected = [1,1]
        self.assertAlmostEqual(res, expected, delta = 1e-6)
      
        
    