# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 20:16:48 2019
@author: Mattias Lundström1
"""
#from mpi4py import MPI
from numpy import diag, ones, zeros, array
import numpy as np
from scipy.linalg import block_diag, solve
from Problem import Problem

class roomHeatSolverSmall:
    
    def __init__(self, problem):
        self.problem = problem
        self.dx = problem.dx
        self.BC_normal = problem.wall
        self.BC_heater = problem.heater
        self.BC_window = problem.window
        self.n = int(1/self.dx)
        ####
        self.solvesSmallRoom(self.n)            #Denna behövs ej i init?
        self.u = 0 #Old room 
        self.interface_1 = 20
        self.interface_2 = 20


    def getBoundaries(self, room):
        if self.u == 0:
            print("Room not setup")
        if room == "room1": # We want boundaries for room 1 -> solve room 2 and return the derivates
            return self.interface_1
        elif room == "room3":
            return self.interface_2
            #Calculate derivates at boundary 
        else:
            #raise exception
            print("Room not defined")
            return None

#    def updateBoundaries(self, interface = interface1):
#        n = self.n
#        if interface:
#            pass



    def solveSmallRoom(self, interface_E, interface_W):
        n = self.n

        A1 = diag(-4*ones(n-1)) + diag(ones(n-2), 1) + diag(ones(n-2), -1) 
        A = block_diag(A1, A1) #First block
        for i in range(2,n-1): # All other building blocks (to 2n-1 for 2x1 block)
            A = block_diag(A, A1)

        A = A + diag(ones(A.shape[0]-(n-1)), n-1) +  diag(ones(A.shape[0]-(n-1)), -(n-1))
        N = A.shape[0] # Size of A

        # Create arrays of boundry conditions. North, West, East, South
        BC_N = array([[*self.BC_heater*ones(n-1), *zeros(N - (n-1))]]).T #OK
        BC_S = array([[*zeros(N - (n-1)), *self.BC_window*ones(n-1)]]).T #OK

        BC_W = 25*array([zeros(N)]).T  #OK
        for i in range(0, int(np.ceil(N/2)), n-1): #Sets interface also to BC_normal temperature
            BC_W[i] = self.BC_normal

        BC_E = array([zeros(N)]).T  #OK
        for i in range(N-1, int(np.ceil(N/2)), -(n-1)): #Sets interface also to BC_normal temperature
            BC_E[i] = self.BC_normal

        # Collect dirichlet boundry conditions in b, then solve Au = b
        b = - BC_N - BC_S - BC_W - BC_E
        A = A/(self.dx**2) 
        b = b/(self.dx**2)
        self.u = solve(A, b)
        print(self.u.reshape(2*n-1, n-1))
        return self.u


    


#p = Problem(0.2)
#solver = roomHeatSolver(p)




