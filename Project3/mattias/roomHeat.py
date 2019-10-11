# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 20:16:48 2019
@author: Mattias Lundström1
"""
from mpi4py import MPI
from numpy import diag, ones, zeros, array
import numpy as np
from scipy.linalg import block_diag, solve




dx = 1/6
n = int(1/dx)
BC_normal = 15
BC_heater = 40
BC_window = 5

#Create matrix A for Au=b with 2nd finite diferences.
# Big room, 2x1. Starting in top left corner
# We get an (2n-1)X(n-1) inner matrix, and an (2n+1)X(n+1) outer matrix with BP

A1 = diag(-4*ones(n-1)) + diag(ones(n-2), 1) + diag(ones(n-2), -1) 
A = block_diag(A1, A1) #First block
for i in range(2,2*n-1): # All other building blocks (to 2n-1 for 2x1 block)
    A = block_diag(A, A1)

A = A + diag(ones(A.shape[0]-(n-1)), n-1) +  diag(ones(A.shape[0]-(n-1)), -(n-1))
N = A.shape[0] # Size of A

# Create arrays of boundry conditions. North, West, East, South
BC_N = array([[*BC_heater*ones(n-1), *zeros(N - (n-1))]]).T #OK
BC_S = array([[*zeros(N - (n-1)), *BC_window*ones(n-1)]]).T #OK

BC_W = 25*array([zeros(N)]).T  #OK
for i in range(0, int(np.ceil(N/2)), n-1): #Sets interface also to BC_normal temperature
    BC_W[i] = BC_normal

BC_E = array([zeros(N)]).T  #OK
for i in range(N-1, int(np.ceil(N/2)), -(n-1)): #Sets interface also to BC_normal temperature
    BC_E[i] = BC_normal

# Collect dirichlet boundry conditions in b, then solve Au = b
b = - BC_N - BC_S - BC_W - BC_E
A = A/(dx**2) 
b = b/(dx**2)
u = solve(A, b)
#print()
print(u.reshape(2*n-1, n-1)) #For 2x1 room inner points
