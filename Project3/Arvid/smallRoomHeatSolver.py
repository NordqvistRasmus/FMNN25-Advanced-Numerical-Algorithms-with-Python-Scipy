#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 11 15:49:31 2019

@author: Arvid
"""
from scipy import *
from numpy import *
from matplotlib.pyplot import *
from scipy.linalg import lu_factor, lu_solve

from Problem import Problem


class smallRoomHeatSolver():
    
    def __init__(self, interface_dir, interface_vals, problem, room):
     #            geom=(1,1), heater=40, normal_wall=15):
        #super().__init__(problem)
        
        self.interface_dir = interface_dir
        self.interface_vals = interface_vals
        #Change with global geometry
        self.x_len = problem.geometry[room][0]
        self.y_len = problem.geometry[room][1]
        self.dx = problem.dx
        self.heater = problem.heater
        self.normal_wall = problem.wall
        self.n_rows = round(self.x_len/self.dx) -1 # Defines the number of rows in the coordinate mesh.
        self.n_cols = round(self.y_len/self.dx)    # Defines the number of columns in the coordinate mesh.  
        self.size = (self.n_rows, self.n_cols)
        self.N_elements = self.n_rows*self.n_cols # Number of points in which to calculate u 
        BC, neu_ind = self._make_boundaries()
        self.BC = BC
        self.neu_ind = neu_ind
        self.A = self._make_matrix()
        lu, piv = lu_factor(self.A)
        self.lu = lu
        self.piv = piv
        self.solution = None
        
    def _make_boundaries(self):
        
        BC_W = zeros(self.size)
        BC_E = zeros(self.size)
        BC_N = zeros(self.size)
        BC_S = zeros(self.size)
        BC_N[0,0:] = self.normal_wall
        BC_S[-1,0:] = self.normal_wall
        
        if self.interface_dir == 'west':
            BC_E[:,-1] = self.heater/(self.dx**2)
            BC_W[:,0] = self.interface_vals/self.dx
            neumann_ind = nonzero(BC_W.reshape(self.N_elements))
        elif self.interface_dir == 'east':
            BC_E[:,-1] = self.interface_vals/self.dx
            BC_W[:,0] = self.heater/self.dx**2
            neumann_ind = nonzero(BC_E.reshape(self.N_elements))
        BC_tot = BC_W + BC_E + BC_N/self.dx**2 + BC_S/self.dx**2
        
        BC_tot = BC_tot.reshape(self.size[0]*self.size[1])
        return BC_tot, neumann_ind
    
    def _update_boundaries(self, interface_vals):
        self.interface_vals = interface_vals
        BC, neu = self._make_boundaries()
        self.BC = -BC
    
    def _make_matrix(self):
        A = (diag(-4*ones(self.N_elements))
            + diag(ones(self.N_elements-1), -1)
            + diag(ones(self.N_elements-1), 1)
            + diag(ones(self.N_elements-self.n_cols), self.n_cols)
            + diag(ones(self.N_elements-self.n_cols), -self.n_cols))
        for ind in self.neu_ind:
            A[ind,ind] = -3
        
        for i in range(self.n_cols-1, self.N_elements-1, self.n_cols):
            A[i, i+1] = 0
        for i in range(self.n_cols, self.N_elements-1, self.n_cols):
            A[i,i-1] = 0
        return A*(1/self.dx**2)
    
    def solve_system(self, interface_vals):
        self._update_boundaries(interface_vals)
        u = lu_solve((self.lu, self.piv), self.BC)
        mesh_vals = u.reshape(self.n_rows,self.n_cols)
        if self.interface_dir == 'east':
            interface_vals = mesh_vals[:,-1]
        elif self.interface_dir == 'west':
            interface_vals = mesh_vals[:,0]
        self.solution = u    
        return u, interface_vals
        
    def getMatrix(self):
        room = zeros((self.n_rows+2, self.n_cols+1))
        if self.interface_dir == 'east':
            room[1:-1,0:-1] = flip(self.solution.reshape(self.size)) #Might have flipped to much heh (mirror flip?)
            room[0, :] = self.normal_wall*ones(self.n_cols+1)
            room[:, -1] = self.heater*ones(self.n_rows+2)
            room[-1, :] = self.normal_wall*ones(self.n_cols+1) 
        elif self.interface_dir == 'west':
            room[1:-1, 1:] = flip(self.solution.reshape(self.size)) #Might have flipped to much heh (mirror flip?)
            room[0,:] = room[-1,:] = self.normal_wall*ones(self.n_cols+1)
            room[:, 0] = self.heater*ones(self.n_rows+2)
        print('Complete room is: {}'.format(room))
        return room
    
    
if __name__ == '__main__':
    p = Problem(1/4)
    print(p.geometry)
    interface_vals = array([20,20,20])
    s = smallRoomHeatSolver('east', interface_vals, p, 'room1')
    #BC, neumann_ind = s._make_boundaries()
    A=s._make_matrix()
    s.solve_system(interface_vals)
    print(s.getMatrix())
        