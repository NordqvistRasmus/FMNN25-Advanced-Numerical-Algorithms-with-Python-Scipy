#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 13:02:49 2019

@author: johanliljegren
"""
from  scipy import *
from  pylab import *

from mpi4py import MPI
from roomHeatSolver import roomHeatSolver
from Problem import Problem 

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

dx = int(1/10)

prob = Problem(dx)
room1 = roomHeatSolver(prob)
room2 = roomHeatSolver(prob)
room3 = roomHeatSolver(prob)

nbrits = 10 #hur många gånger vi vill iterera

for i in range(nbrits):
    if rank == 0: #room2
        if(i==0): #first iteration
            room2.solveLargeRoom()
        else:
            comm.recv(bound1, source = 1)
            comm.recv(bound3, source = 2)
            room2.updateBound(bound1, bound3, room = 'room2')
            room2.solveLargeRoom(()
        comm.send(room2.getBound(), soruce = 1)
        comm.send(room2.getBound(), soruce = 2)
    
    if rank == 1: #room1
        comm.recv(bound2, source = 0)
        room1.updateBound(bound2, room = 'room1')
        room1.solveSmallRoom()
        comm.send(room1.getBound(), source = 0)
    
    if rank == 2: #room3
        comm.recv(bounds, source = 0)
        room3.updateBound(bound2, room = 'room3')
        room3.solveSmallRoom()
        comm.send(room3.getBound(), source = 0)
        
A = room1.getMatrix()
B = room2.getMatrix()
C = room3.getMatrix()