#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 14:42:55 2019

@author: Arvid
"""
from scipy import *
from matplotlib.pyplot import *
from timeit import default_timer as timer


class CubicSpline:
    
    """
    CubicSpline
    
    A class used to represent a CubicSpline with the purpose of curve design in 2D. 

    ...

    Attributes
    ----------
    KNOTS : float
        Sequence of floats
    CONTROL : tuples
        List of tuples. Control points, used for shaping the curve.

    Methods 
    -------
    __init__(control, knots = None):
        Sorts the knot grid and adds padding. 
        If no knot sequence is passed, the default will be a equidistant vector from 0 to 1.
    
    __call__(u):
        See method point_eval(u). 
    
    plot(u_i = None, plot_poly = True, precision = 150):
        Calculates and plot points using De Boors algorithm.
        If 'plot_poly' is set True (default), then also plots the Control Points.  
        To increase plot resolution, increase 'precision' parameter.
        
        u_i is used for Add-on task. Adds parts of the blossoms points to the plot. 
        
    point_eval(u):
        Help method to method 'plot', used to evaluate point 'u' with the De Boors algorithm and calculate splines. 
        Input: A point 'u' to be evaluated, needs to be inside the knot sequence.
        Return: Blossom point(in R^2) 
            
    basis_fuction(knots, i):
        Takes a knot sequence 'knots' and index 'i' as input and returns a 
        python function to evaluate the i:th B-spline basis function. 
        
    ------
    """
    
    """
    Sorts the knot grid and adds padding. 
    If no knot sequence is passed, the default will be a equidistant vector from 0 to 1.
    """
        
    def __init__(self, control, knots = None, Inter=False):
        
        if control is None:
            raise TypeError('Control points are not defiened properly')
       
        if shape(control[0]) != (2,): #Exception for input shape
            raise ValueError('Control points are not two dimensional')
        
        self.control = array(control)
        
        if knots is None:
            knots = linspace(0,1, len(control) - 2) #Padding adds 5 points, and we need knots length to be 2 more than control
            self.knots = r_[knots[0], knots[0], knots, knots[-1], knots[-1], knots[-1]]
        else: #Knots are provided, need to sort. 
            self.knots = knots
            self.knots.sort()
        
        if len(self.knots) != len(self.control) + 3: #Exception for input length
            raise IndexError('The length of control points, K, and knot grid, L, \
            need to match with K+3=L. "\n"Current lengths are K =', len(self.control), '\
            and L = ', len(self.knots))
            
        if Inter:
            self.control=control
            knots=np.linspace(0,1,len(self.control)-2)
            self.knots = np.concatenate([[0,0],knots,[1,1]])
            self.interpol(self.control)           
    
   
    def __call__(self, u):
        if not min(self.knots) <= u <= max(self.knots): 
            raise IndexError('Input in is not a value inside the knot sequence.')
        else:
            return self.point_eval(u)
    
    """
    Calculates and plot points using De Boors algorithm. 
    If 'plot_poly' is set True (default), then also plots the Control Points.  
    To increase plot resolution, increase 'precision' parameter
    """
    def plot(self, u_i = None, plot_poly = True, precision = 150):
        figure(1)
        
        #   Add-on 1:
        if u_i is not None:
            if u_i not in self.knots: 
                raise IndexError('Input in is not a value of the knot sequence.')
            markIndx = searchsorted(self.knots, u_i)
            self.markControl = []
            self.markBloss = []
            self.points = array([self.point_eval(point, markIndx) for point in linspace(0,1,precision)])
            plot([x[0] for x in self.markControl], [y[1] for y in self.markControl], 'bs') #Control points
            plot([x[0] for x in self.markBloss], [y[1] for y in self.markBloss], 'r1') #Blossoms
        #   Add-on 1:
        
        points = array([self.point_eval(point) for point in linspace(0,1, precision)])
        #print(points)
        plot(points[:,0],points[:,1])
        if plot_poly:
            plot(self.control[:,0], self.control[:,1], 'yx--')
        
    """
    Help method to method 'plot', used to evaluate point 'u' with the De Boors algorithm and calculate splines. 
    """
    def point_eval(self, u, markIndx = None):
        indx = int(self.knots.searchsorted([u])) - 1
        
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
        
        # Add on attemp 
        if markIndx == indx:
            [self.markControl.append(x) for x in blossoms0]
            [self.markBloss.append(x) for x in blossoms1]
            [self.markBloss.append(x) for x in blossoms2]
        
        return blossoms3

    """
    Takes a knot sequence 'knots' and index 'i' as input and returns a 
    python function to evaluate the i:th B-spline basis function. 
    """
    def basis_function(self, knots, i):  
        def basis(u, i, k):
            if k == 0:
                if knots[i-1] == knots[i]:
                    return 0
                elif knots[i-1] <= u < knots[i]:
                    return 1
                else:
                    return 0
            else:
                if knots[i+k-1] == knots[i-1]: # and u == knots[i-1]:
                    coeff1 = 0
                else:
                    coeff1 = (u - knots[i-1])/(knots[i+k-1]-knots[i-1])
                if knots[i+k] == knots[i]: # and u == knots[i+k]:
                    coeff2 = 0
                else:
                    coeff2 = (knots[i+k] - u)/(knots[i+k] - knots[i])
                
                return basis(u, i, k-1)*coeff1 + basis(u, i+1, k-1)*coeff2
        return lambda u: basis(u, i, 3)
    
    def interpol(self,points):

        d = points #Our interpolationpoints
        grevabs=np.zeros([len(self.knots)-2])

        for i in range(len(self.knots)-2):
            grevabs[i]=(self.knots[i]+self.knots[i+1]+self.knots[i+2])/3 #Our Greville abscissae vector

        VanderMatrix = np.zeros((len(grevabs),len(grevabs))) #Empty matrix for our A in Ax=b

        for i in range(len(grevabs)):

            for j in range(len(grevabs)):
                controlBase = zeros([len(self.knots)-2,2]) #For creating our b-spline basis
                controlBase[j]=[1,1]
                bspline=CubicSpline(controlBase,self.knots) #Creates the b-spline base    
                u=(np.array([grevabs[i]])) 
                VanderMatrix[i,j] = bspline(u)[0] #Evalues the point u using our bspline

        x = linalg.solve(VanderMatrix,d[:,0]) #Solves for our vector of x coords
        y = linalg.solve(VanderMatrix,d[:,1]) ##Solves for our vector of y coords

    

        controlpoints = np.zeros((len(x),2))

        for i in range(len(x)): #Gatheres the coordinates of the controlpoints
            controlpoints[i,0]=x[i]
            controlpoints[i,1]=y[i]
        self.control=controlpoints
        alfa=CubicSpline(controlpoints) #Cubic Spline the controlpoints

    """
    Example program if class executed as main.
    """
if __name__ == '__main__':
    
    #print(c.basis_function(KNOTS,0))
    #print(c.point_eval(0.2))  
    #c.plot()
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
        
    KNOTS = linspace(0, 1, 26)
    KNOTS[ 1] = KNOTS[ 2] = KNOTS[ 0]
    KNOTS[-3] = KNOTS[-2] = KNOTS[-1]
    KNOTS = r_[KNOTS, KNOTS[-1]] #Fix index error with basis function
    
    
    c = CubicSpline(CONTROL, KNOTS)
    #c = CubicSpline(CONTROL)
    #c.plot()
    start = timer()
    c.plot(0.52, True, 200)
    end = timer()
    print(end-start, 'seconds.')
    
    #Cant get the last basis function to show (-2 instead of -3), index error
    
    basis = [c.basis_function(KNOTS, i) for i in range(len(KNOTS)-3)]
    X = linspace(0, 1, 200)
    Y = array([[N(x) for x in X] for N in basis])
    figure(2); [plot(X, y) for y in Y]
    
    #Interpol test
    #points = np.array([[0,3],[3,6],[5.5,6.5],[4,0],[2,-3],[0,-6],
    #                 [-2,-3],[-4,0],[-5.5,6.5],[-3,5.5],[0,3],[0,3]]) 
    # b = CubicSpline(points, Inter = True)
    #b.plot(plot_poly=False) #plot it
