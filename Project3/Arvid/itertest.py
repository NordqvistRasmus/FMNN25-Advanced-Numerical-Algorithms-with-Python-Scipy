#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 20:56:07 2019

@author: Arvid
"""
from scipy import *
from matplotlib.pyplot import *
from numpy import zeros, array, diag, ones

from mpi4py import MPI
from roomHeatSolver import roomHeatSolver
from smallRoomHeatSolver import smallRoomHeatSolver
from largeRoomHeatSolver import largeRoomHeatSolver
from Problem import Problem 
import seaborn as sns; sns.set()

dx = (1/20)
n = 20

prob = Problem(dx)
room1 = smallRoomHeatSolver( 'east', zeros(n - 1), prob,"room1")
room2 = largeRoomHeatSolver(prob)
room3 = smallRoomHeatSolver('west', zeros(n - 1), prob,"room3")

omega  = 0.8
bound1_old = 20*ones(n-1)
bound3_old = 20*ones(n-1)

nbrits = 10 #hur många gånger vi vill iterera

for i in range(nbrits):
    
    # Room 2
    if i == 0:
        room2.solveLargeRoom()
    else:
        room2.updateBound('interface1', bound1)
        room2.updateBound('interface2', bound3)
        room2.solveLargeRoom()
    bounds_r1 = room2.getDerives('interface1')
    bounds_r3 = room2.getDerives('interface2')
    
    # Room1
    u1, bound1 = room1.solve_system(bounds_r1)
    bound1 = omega*bound1 +(1-omega)*bound1_old
    bound1_old = bound1
    # Room3
    u2, bound3 = room3.solve_system(bounds_r3)
    bound3 = omega*bound3 + (1-omega)*bound3_old
    bound3_old = bound3
    
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
    
    
    
    