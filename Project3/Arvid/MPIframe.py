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
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

dx = (1/20)
n = 20

prob = Problem(dx)
room1 = smallRoomHeatSolver( 'east', zeros(n - 1), prob,"room1")
room2 = largeRoomHeatSolver(prob)
room3 = smallRoomHeatSolver('west', zeros(n - 1), prob,"room3")

omega  = 0.8
bound1_old = 20*ones(n-1)
bound3_old = 20*ones(n-1)

nbrits = 2 #hur många gånger vi vill iterera

for i in range(nbrits):
    if rank == 0: #room2
        if(i==0): #first iteration
            room2.solveLargeRoom()
        else:
            bound1 = comm.recv(source = 1)
            bound3 = comm.recv(source = 2)
            room2.updateBound('interface1', bound1)
            room2.updateBound('interface2', bound3)
            room2.solveLargeRoom()
            #print('Room2: {}'.format(room2.getMatrix()))
            
        bounds_r1 = room2.getDerives('interface1')
        bounds_r3 = room2.getDerives('interface2')
        comm.send(bounds_r1, dest = 1)
        comm.send(bounds_r3, dest = 2)
    
    if rank == 1: #room1
        bounds_r1 = comm.recv(source = 0)
        u, bound1 = room1.solve_system(bounds_r1)
        bound1 = omega*bound1 +(1-omega)*bound1_old
        comm.send(bound1, dest = 0)
        bound1_old = bound1
        #print('Room1: {}'.format(room1.getMatrix()))
    
    if rank == 2: #room3
        bounds_r3 = comm.recv(source = 0)
        print('bounds_r3: {}'.format(bounds_r3))
        u, bound3 = room3.solve_system(bounds_r3)
        bound3 = omega*bound3 + (1-omega)*bound3_old
        bound3_old = bound3
        comm.send(bound3, dest = 0)
        #print('Room3: '.format(room3.getMatrix()))
    A = room1.getMatrix()
        
A = room1.getMatrix()
B = room2.getMatrix()
C = room3.getMatrix()

fig, ax = plt.subplots()

upper_left = zeros((A.shape[0],A.shape[1]))
lower_right = zeros((C.shape[0],C.shape[1]))
upper_left.fill(None)
lower_right.fill(None)

splitted_room2 = vsplit(array(B), 2)
total = block([[upper_left, splitted_room2[0], C, A, splitted_room2[1], lower_right]])
ax = sns.heatmap(total, cmap = "YlOrRd")
plt.show()

