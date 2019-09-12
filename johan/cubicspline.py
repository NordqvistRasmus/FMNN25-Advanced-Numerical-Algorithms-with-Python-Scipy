#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 11:38:09 2019

@author: johanliljegren
"""
from  scipy import *
from  pylab import *

class CubicSpline:
    def __init__(self, cp, p, u = None):
        self.d = array[cp]
        self.p = p
        If u == None:
            u = linspace(0, 10, len(cp)-1)
        else:
            u.sort()
        self.u = r_[u[0], u[0], u, u[-1], u[-1]]    
        
        
        
    def __call__(self, u):
        print('hej')
        pass
     '''
     input:     control points: d
                knot vector: u
                the point: t
                
                
     '''
    def blossom(t):
        index = hot_intervall(t, u)
        hot_int = zeros(p+1)
        for k in range(-2,2):
            hot_int[k+2] = u(index + k)
        for i in range(0,p):
            for j in range(i,p):
                a = (u[j-i]-t) / (u[j-i]-u[j-p])
                d[j] = a*d[j]*(1-a)*d[j+1]
        s = d[2]
        return s
    
    def hot_intervall(t): """ U är instoppade värdet... """
        u.sort()
        idx = u.searchsorted([t]])
        index = idx[0]
        return index
    
    def plot():
        pass
    
sasdas dasd asd asd asd
    