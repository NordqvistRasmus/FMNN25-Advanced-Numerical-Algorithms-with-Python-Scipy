#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 13:02:49 2019

@author: Mattias Lundström, Arvid Rolander, Pontus Nordqvist, Johan Liljegren, Antonio Alas
"""
from numpy import zeros, array, diag, ones, vsplit, block

from mpi4py import MPI
#from roomHeatSolver import roomHeatSolver
from smallRoomHeatSolver import SmallRoomHeatSolver
from largeRoomHeatSolver import LargeRoomHeatSolver
from Problem import Problem 
import matplotlib.pyplot as plt
import seaborn as sns; sns.set()


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

dx = (1/20)
n = 20

prob = Problem(dx)
room1 = SmallRoomHeatSolver('east', zeros(n - 1), prob, 'room1')
room2 = LargeRoomHeatSolver(prob)
room3 = SmallRoomHeatSolver('west', zeros(n - 1), prob, 'room3')

omega  = 0.8
bound1_old = 20*ones(n-1)
bound3_old = 20*ones(n-1)

nbrits = 10 #hur många gånger vi vill iterera

for i in range(nbrits):
    if rank == 0: #room2
        print("Rank is 0")
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
        print("Rank is 1")
        bounds_r1 = comm.recv(source = 0)
       
        u, bound1 = room1.solve_system(bounds_r1)
        bound1 = omega*bound1 +(1-omega)*bound1_old
        comm.send(bound1, dest = 0)
        bound1_old = bound1
        #print('Room1: {}'.format(room1.getMatrix()))
    
    if rank == 2: #room3
        print("Rank is 2")
        bounds_r3 = comm.recv(source = 0)
       
        #print('bounds_r3: {}'.format(bounds_r3))
        u, bound3 = room3.solve_system(bounds_r3)
        bound3 = omega*bound3 + (1-omega)*bound3_old
        bound3_old = bound3
        comm.send(bound3, dest = 0)
        #print('Room3: '.format(room3.getMatrix()))
    if(i == nbrits-1):
        if rank == 0:
            B = room2.getMatrix()
            comm.send(B, dest=3, tag=2)
        if rank == 1:
            A = room1.getMatrix()
            comm.send(A, dest=3, tag=1)
        if rank == 2:
            C = room3.getMatrix()
            comm.send(C, dest=3, tag=3)
if rank == 3:
    A = comm.recv(source = 1, tag=1)
    C = comm.recv(source = 2, tag=3)
    B = comm.recv(source=0, tag=2)

    
    fig, ax = plt.subplots()

    upper_left = zeros([A.shape[0] - 1, A.shape[1]])
    lower_right = zeros([C.shape[0] - 1,C.shape[1]])
    upper_left.fill(None)
    lower_right.fill(None)
    

    First = block([[upper_left],[A]])
    Second = B
    Third = block([[C], [lower_right]])
    print(First)
    print(B)
    print(Third)
    total = block([First, Second, Third])
    print('-------------------------------------------------')
    print(total)
    ax = sns.heatmap(total, cmap = "YlOrRd")
    #ax.invert_yaxis()
    plt.show()
    
    
    

