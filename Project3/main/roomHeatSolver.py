#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 15:09:42 2019

@author: johanliljegren
"""
from  scipy import *
from  pylab import *

class roomHeatSolver:
    
    def __init__(self, problem):
        self.problem = problem
        self.dx = dx
        self.BC_normal = problem.wall
        self.BC_heater = problem.heater
        self.BC_window = problem.window
        self.n = int(1/dx)
    
    def solveRoom(self):
        pass
    
    def getBound(self):
        pass
    
    def updateBound(self):
        pass
    
    def relax(self):
        pass
    
    def getMatrix(self):
        pass
    