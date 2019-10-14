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
        self.dx = problem.dx
        self.wall = problem.wall
        self.heater = problem.heater
        self.window = problem.window
        self.n = int(1/self.dx)
    
    def solveRoom(self):
        pass
    
    def getBound(self):
        pass
    
    def getDerives(self):
        pass
    
    def updateBound(self):
        pass
    
    def relax(self):
        pass
    
    def getMatrix(self):
        pass
    