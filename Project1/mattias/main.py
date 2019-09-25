
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:19:55 2019
@author: Mattias Lundström1
"""
from  scipy import *
from  pylab import *
from CubicSpline_m import CubicSpline as spline

grid = array([1, 3, 5, 7, 8, 9]) # K = 6 
cp = array([[1, 3], [3.5, 8], [5, 4],[5.4, 5],[6.6, 2], [7, 4], [8, 3]]) # L = K - 2 = 4

s = spline(cp)
res = s(50)

s.plot()