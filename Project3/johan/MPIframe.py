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

dx = int(1/20)

prob = Problem(dx)
room1 = smallRoomHeatSolver(prob)
room2 = largeRoomHeatSolver(prob)
room3 = smallRoomHeatSolver(prob)

nbrits = 10 #hur många gånger vi vill iterera

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
        comm.send(room2.getDerives('interface1'), source = 1)
        comm.send(room2.getDerives('interface2'), source = 2)
    
    if rank == 1: #room1
        comm.recv(bounds, source = 0)
        room1.solveSmallRoom(bounds)
        comm.send(room1.getBound(), source = 0)
    
    if rank == 2: #room3
        comm.recv(bounds, source = 0)
        [] = room3.solveSmallRoom(bounds)
        comm.send(room3.getBound(), source = 0)
        
A = room1.getMatrix()
B = room2.getMatrix()
C = room3.getMatrix()