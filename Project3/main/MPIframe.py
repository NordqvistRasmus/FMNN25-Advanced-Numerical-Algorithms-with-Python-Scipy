#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 13:02:49 2019

@author: johanliljegren
"""
#from  scipy import *
#from  pylab import 
from numpy import zeros, array, diag, ones

from mpi4py import MPI
from roomHeatSolver import roomHeatSolver
from smallRoomHeatSolver import smallRoomHeatSolver
from largeRoomHeatSolver import largeRoomHeatSolver
from Problem import Problem 
import seaborn as sns; sns.set()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

dx = (1/20)
n = 20

prob = Problem(dx)
room1 = smallRoomHeatSolver( 'east', zeros(n - 1), prob,"room1")
room2 = largeRoomHeatSolver(prob)
room3 = smallRoomHeatSolver('west', zeros(n - 1), prob,"room3")

nbrits = 2 #hur många gånger vi vill iterera

for i in range(nbrits):
    if rank == 0: #room2
        if(i==0): #first iteration
            room2.solveLargeRoom()
        else:
            comm.recv(bound1, source = 1)
            comm.recv(bound3, source = 2)
            room2.updateBound('interface1', bound1)
            room2.updateBound('interface2', bound3)
            room2.solveLargeRoom()
        bounds_r1 = room2.getDerives('interface1')
        bounds_r3 = room2.getDerives('interface2')
        comm.send(bounds_r1, source = 1)
        comm.send(bounds_r3, source = 2)
    
    if rank == 1  and i != 0: #room1
        comm.recv(bounds_r1, source = 0)
        u, bound1 = room1.solve_system(bounds_r1)
        comm.send(bound1, source = 0)
    
    if rank == 2  and i != 0: #room3
        comm.recv(bounds_r3, source = 0)
        u, bound3 = room3.solve_system(bounds_r3)
        comm.send(bound3, source = 0)
        
#A = room1.getMatrix()
B = room2.getMatrix()
""" #C = room3.getMatrix()

fig, ax = plt.subplots()

upper_left = zeros((A.shape[0],A.shape[1]))
lower_right = zeros((C.shape[0],C.shape[1]))
upper_left.fill(None)
lower_right.fill(None)

splitted_room2 = vsplit(array(B), 2)
total = block([[upper_left, splitted_room2[0], C, A, splitted_room2[1], lower_right]])
ax = sns.heatmap(total, cmap = "YlOrRd")
plt.show()

 """