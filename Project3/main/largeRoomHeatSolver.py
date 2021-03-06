# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 20:16:48 2019
@author: Mattias Lundström, Arvid Rolander, Pontus Nordqvist, Johan Liljegren, Antonio Alas
"""
from mpi4py import MPI
from numpy import diag, ones, zeros, array, block
import numpy as np
from scipy.linalg import block_diag, solve
from Problem import Problem

class LargeRoomHeatSolver:
    
    def __init__(self, problem):
        self.problem = problem
        self.dx = problem.dx
        self.n = int(1/self.dx)
        self.heater = problem.heater
        self.wall = problem.wall
        self.window = problem.window
        self.interfaceArray1 = 20*ones(self.n-1) # Interface 1 and 2 start with 20 degress. 
        self.interfaceArray2 = 20*ones(self.n-1)
        self.u = None
        self.u_prev = None # Previous room, used for relaxation
        #### for testing
        #self.solveLargeRoom()
        ## Set interface value guess to 20, size of interface is n-1.
        


    def getBound(self, interface = None):
        if self.u is None:
            print("Room needs setup")
        if (interface == "interface1"): # We want boundaries for room 1 -> solve room 2 and return the derivates
            return self.interfaceArray1
        if (interface == "interface2"):
            return self.interfaceArray2  
        else:
            raise ValueError("Can only get boundaries at interface1 or interface2.")


    def updateBound(self, interface, boundArray):
        if (boundArray.shape !=  self.interfaceArray1.shape):
            print("BoundArray not same length as interfaceArray. Not possible to update.")
            raise IndexError()
        if (interface == "interface1"):
            self.interfaceArray1 = boundArray
            print("Boundaries at interface1 updated.\n")
            return
        if (interface == "interface2"):
            self.interfaceArray2 = boundArray
            print("Boundaries at interface2 updated.\n")
            return
        else:
            raise ValueError("Can only update boundaries at interface1 or interface2.")
        
    def getMatrix(self):
        if self.u is None:
            raise TypeError('U is not defined')
        
        room = np.zeros((2*self.n+1, self.n+1))
        size = room.shape
        u_matrix = self.u.reshape(2*self.n -1, self.n - 1)
        R = array([*self.interfaceArray2, *self.wall*ones(self.n)])
        L = array([*self.wall*ones(self.n), *self.interfaceArray1])
               
        room[0, :] = self.problem.heater*ones(self.n+1)
        room[1:-1, 0] = L
        room[1:-1, 1:-1]=u_matrix
        room[1:-1, -1] = R
        room[-1, :] = self.problem.window*ones(self.n+1)

        return room

    def getDerives(self, interface):
        """
        Calculates the derivates along interface1 or interface2 and returns them as a vector. 
        """
        n = self.n
        u = self.u

        if (interface == "interface1"):
            #Calculate derivates at boundary 
            u_matrix = u.reshape(2*n-1, n-1)
            u_i = u_matrix[:,0] #get first row, since interface1
            u_i = u_i[n:] #slice array, get only elements along interface1
            derivates = (self.interfaceArray1 - u_i) / self.dx
            return(derivates)
        if (interface == "interface2"):
            u_matrix = u.reshape(2*n-1, n-1)
            u_i = u_matrix[:,-1] #get last row, since interface1
            u_i = u_i[:n-1] #slice array, get only elements along interface2
            derivates = (u_i - self.interfaceArray2) / self.dx  #WANT POSITIVE OR NEGATIVE DERIVATES?
            return(derivates)
        else:
            raise ValueError('Interface not defined properly. Choose interface1 or interface2') 


    def solveLargeRoom(self):
        """
        Solves the large room with a linear system using interfaceArray1 and interfaceArray2 values along the room interfaces.
        Wall, heater and window temperatures are defined in the problem class. 

        Returns u, the solved inner room temperaturesm as a vector which can be reshaped back into a room matrix.
        """
        n = self.n

        A1 = diag(-4*ones(n-1)) + diag(ones(n-2), 1) + diag(ones(n-2), -1) 
        A = block_diag(A1, A1) #First block
        for i in range(2,2*n-1): # All other building blocks (to 2n-1 for 2x1 block) to build matrix A,, one block for each line
            A = block_diag(A, A1)

        A = A + diag(ones(A.shape[0]-(n-1)), n-1) +  diag(ones(A.shape[0]-(n-1)), -(n-1))
        N = A.shape[0] # Size of A

        # Create arrays of boundry conditions. North, West, East, South
        BC_N = array([[*self.heater*ones(n-1), *zeros(N - (n-1))]]).T #OK
        BC_S = array([[*zeros(N - (n-1)), *self.window*ones(n-1)]]).T #OK

        BC_W = 25*array([zeros(N)]).T  #OK
        for i in range(0, int(np.ceil(N/2)), n-1): 
            BC_W[i] = self.wall
        BC_W_M = BC_W.reshape(2*n-1, n-1) #Reshape matrix and insert interface values
        BC_W_M[n:,0] = self.interfaceArray1
        BC_W = array([BC_W.reshape(N)]).T #OK

        BC_E = array([zeros(N)]).T  #OK
        for i in range(N-1, int(np.ceil(N/2)), -(n-1)): 
            BC_E[i] = self.wall
        BC_E_M = BC_E.reshape(2*n-1, n-1)
        BC_E_M[:n-1,-1] = self.interfaceArray2 #Reshape matrix and insert interface values
        BC_E = array([BC_E.reshape(N)]).T #OK

        # Collect dirichlet boundry conditions in b, then solve Au = b
        b = - BC_N - BC_S - BC_W - BC_E

        A = A*(1/(self.dx**2)) 
        b = b*(1/(self.dx**2))
        self.u = solve(A, b)
        #print(self.u)
        return self.u

if __name__ == '__main__':
    p = Problem(1/5)
    solver = LargeRoomHeatSolver(p)
    solver.solveLargeRoom()
    #print(solver.getDerives("interface2"))
    #print(solver.getMatrix())





