#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:42:55 2019

@author: Arvid
"""
from scipy import *
from matplotlib.pyplot import *

class CubicSpline:
    
    
    def __init__(self, knots, control):
        #self.knots = r_[knots[0], knots[0], knots, knots[-1], knots[-1]]
        self.knots = knots
        self.knots.sort()
        control = array(control)
        self.control = control
        
    def __call__(self):
        pass
    
    def plot(self, u_i = None, plot_poly = True, precision = 200):
        if u_i is not None:
            markIndx = searchsorted(self.knots, u_i)
            self.markControl = []
            self.markBloss = []
            self.points = array([self.point_eval(point, markIndx) for point in linspace(0,1,precision)])
            plot([y[0] for y in self.markControl], [y[1] for y in self.markControl], 'bs') #Control points
            plot([y[0] for y in self.markBloss], [y[1] for y in self.markBloss], 'ro')
            print(self.markBloss)
        else:
            self.points = array([self.point_eval(point) for point in linspace(0,1,precision)])
            
        if plot_poly:
            plot(self.control[:,0], self.control[:,1], 'yx--')
        plot(self.points[:,0], self.points[:,1])
        show()
    
    def point_eval(self, u, markIndx = None):
            
        indx = int(self.knots.searchsorted([u]))-1   

        blossoms0 = array([self.control[i] for i in range(indx-2, indx+2)])
        knots = array([self.knots[i] for i in range(indx-2, indx+4)])
        
        alphas = array([((knots[i+3] - u)/(knots[i+3] - knots[i])) for i in range(3)])
        blossoms1 = array([alphas[i]*blossoms0[i] + (1 - alphas[i])*blossoms0[i+1]
                    for i in range(3)])

        alphas = array([(knots[i+3] - u)/(knots[i+3] - knots[i+1]) for i in range(2)])
        blossoms2 = array([alphas[i]*blossoms1[i] + (1 - alphas[i])*blossoms1[i+1]
                    for i in range(2)])
        
        alphas = (knots[3]-u)/(knots[3]-knots[2])
        blossoms3 = alphas*blossoms2[0] + (1-alphas)*blossoms2[1]
       
        #print(blossoms)
        if markIndx == indx:
            self.markControl.append(blossoms0[0])
            self.markControl.append(blossoms0[1])
            self.markControl.append(blossoms0[2])
            self.markControl.append(blossoms0[3])
            self.markBloss.append(blossoms1[0])
            self.markBloss.append(blossoms1[1])
            self.markBloss.append(blossoms1[2])
            self.markBloss.append(blossoms2[0])
            self.markBloss.append(blossoms2[1])

        return blossoms3
        
        
    def basis_function(self, knots, i):
        def basis(u,i, k):
            if k == 0:
                if knots[i-1] == knots[i]:
                    return 0
                elif knots[i-1] <= u < knots[i]:
                    return 1
                else:
                    return 0
            else:
                if knots[i+k-1] == knots[i-1] and u == knots[i-1]:
                    coeff1 = 0
                else:
                    coeff1 = (u - knots[i-1])/(knots[i+k-1]-knots[i-1])
                if knots[i+k] == knots[i] and u == knots[i+k]:
                    coeff2 = 0
                else:
                    coeff2 = (knots[i+k] - u)/(knots[i+k] - knots[i])
                
                return basis(u, i, k-1)*coeff1 + basis(u, i+1, k-1)*coeff2
        def j_basis(u):
            return basis(u, i, 3)
            
        return j_basis
                
        
KNOTS = array([0.  , 0.  , 0.  , 0.12, 0.16, 0.2 , 0.24, 0.28, 0.32, 0.36, 0.4 ,
       0.44, 0.48, 0.52, 0.56, 0.6 , 0.64, 0.68, 0.72, 0.76, 0.8 , 0.84,
       0.88, 1.  , 1.  , 1.  ])
    
 
CONTROL = [(-12.73564, 9.03455),
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

if __name__ == '__main__':
    c = CubicSpline(KNOTS, CONTROL)
   # pts = c(0.25)
    c.plot(0.5, True, 200) 


        
        
    
        
        
        
    
        
        
        


