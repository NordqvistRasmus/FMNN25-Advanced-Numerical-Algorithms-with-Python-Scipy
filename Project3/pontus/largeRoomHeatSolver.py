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

from roomHeatSolver import roomHeatSolver

class largeRoomHeatSolver(roomHeatSolver):
    
    def __init__(self, problem):
        super().__init__(problem)
        self.solveLargeRoom(self.n)
        self.u = 0
        
    def solveLargeRoom(self, n):
        A1 = diag(-4*ones(n-1)) + diag(ones(n-2), 1) + diag(ones(n-2), -1) 
        A = block_diag(A1, A1) #First block
        for i in range(2,2*n-1): # All other building blocks (to 2n-1 for 2x1 block)
            A = block_diag(A, A1)

        A = A + diag(ones(A.shape[0]-(n-1)), n-1) +  diag(ones(A.shape[0]-(n-1)), -(n-1))
        N = A.shape[0] # Size of A

        # Create arrays of boundry conditions. North, West, East, South
        BC_N = array([[*self.heater*ones(n-1), *zeros(N - (n-1))]]).T #OK
        BC_S = array([[*zeros(N - (n-1)), *self.window*ones(n-1)]]).T #OK

        BC_W = 25*array([zeros(N)]).T  #OK
        for i in range(0, int(np.ceil(N/2)), n-1): #Sets interface also to BC_normal temperature
            BC_W[i] = self.wall

        BC_E = array([zeros(N)]).T  #OK
        for i in range(N-1, int(np.ceil(N/2)), -(n-1)): #Sets interface also to BC_normal temperature
            BC_E[i] = self.wall

        # Collect dirichlet boundry conditions in b, then solve Au = b
        b = - BC_N - BC_S - BC_W - BC_E
        A = A/(self.dx**2) 
        b = b/(self.dx**2)
        self.u = solve(A, b)
        #print(u)
        
    def getMatrix(self):
        return self.u.reshape(2*self.n-1,self.n-1)
        

        

#if __name__ == '__main__':
p = Problem(0.2)
rhs = roomHeatSolver(p)
lrhs = largeRoomHeatSolver(rhs)
lrhs.solveLargeRoom(lrhs.n)
print(lrhs.getMatrix())


