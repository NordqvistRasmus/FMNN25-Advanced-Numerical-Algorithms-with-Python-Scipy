# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 20:16:48 2019
@author: Mattias Lundström1
"""
from mpi4py import MPI
from numpy import diag, ones, zeros, array
import numpy as np
from scipy.linalg import block_diag, solve
from Problem import Problem

class roomHeatSolver:
    
    def __init__(self, problem):
        self.problem = problem
        self.dx = problem.dx
        self.BC_normal = problem.wall
        self.BC_heater = problem.heater
        self.BC_window = problem.window
        self.n = int(1/self.dx)
        self.interface_1 = 20*ones(self.n-1) 
        self.interface_2 = 20*ones(self.n-1)
        self.u = 0 #Old room 
        #### for testing
        self.solveLargeRoom()
        ## Set interface value guess to 20, size of interface is n-1.
        


    def getBoundaries(self, room):
        if self.u == 0:
            print("Room not setup")

        if room == "room1": # We want boundaries for room 1 -> solve room 2 and return the derivates
            return self.interface_1
        elif room == "room2":
            return (self.interface_1, self.interface_2)
        elif room == "room3":
            return self.interface_2
            
        else:
            #raise exception
            print("Room not defined")
            return None

    def updateBoundaries(self, interface):
        pass

    def getNeumannBC(self, interface):
        n = self.n
        u = self.u
        if (interface == "interface1"):
            #Calculate derivates at boundary 
            u_matrix = u.reshape(2*n-1, n-1)
            u_i = u_matrix[:,0] #get first row, since interface1
            u_i = u_i[n:] #slice array, get only elements along interface
            derivates = (self.interface_1 - u_i) / self.dx
            return(derivates)
        ###
        #mooooore

            
            

    def solveLargeRoom(self):
        n = self.n

        A1 = diag(-4*ones(n-1)) + diag(ones(n-2), 1) + diag(ones(n-2), -1) 
        A = block_diag(A1, A1) #First block
        for i in range(2,2*n-1): # All other building blocks (to 2n-1 for 2x1 block) to build matrix A,, one block for each line
            A = block_diag(A, A1)

        A = A + diag(ones(A.shape[0]-(n-1)), n-1) +  diag(ones(A.shape[0]-(n-1)), -(n-1))
        N = A.shape[0] # Size of A

        # Create arrays of boundry conditions. North, West, East, South
        BC_N = array([[*self.BC_heater*ones(n-1), *zeros(N - (n-1))]]).T #OK
        BC_S = array([[*zeros(N - (n-1)), *self.BC_window*ones(n-1)]]).T #OK

        BC_W = 25*array([zeros(N)]).T  #OK
        for i in range(0, int(np.ceil(N/2)), n-1): 
            BC_W[i] = self.BC_normal
        BC_W_M = BC_W.reshape(2*n-1, n-1) #Reshape matrix and insert interface values
        BC_W_M[n:,0] = self.interface_1
        BC_W = array([BC_W.reshape(N)]).T #OK

        BC_E = array([zeros(N)]).T  #OK
        for i in range(N-1, int(np.ceil(N/2)), -(n-1)): 
            BC_E[i] = self.BC_normal
        BC_E_M = BC_E.reshape(2*n-1, n-1)
        BC_E_M[:n-1,-1] = self.interface_2 #Reshape matrix and insert interface values
        BC_E = array([BC_E.reshape(N)]).T #OK

        # Collect dirichlet boundry conditions in b, then solve Au = b
        b = - BC_N - BC_S - BC_W - BC_E
        print(b.reshape(2*n-1, n-1))
        A = A/(self.dx**2) 
        b = b/(self.dx**2)
        self.u = solve(A, b)
        return self.u


    


p = Problem(0.2)
solver = roomHeatSolver(p)
solver.getNeumannBC("interface1")

#u.reshape(2*n-1, n-1)




