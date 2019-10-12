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


class smallRoomHeatSolver:
    
    def __init__(self, interface_dir, dx, interface_vals,
                 geom=(1,1), heater=40, normal_wall=15):
        
        self.interface_dir = interface_dir
        self.interface_vals = interface_vals
        self.x_len = geom[0]
        self.y_len = geom[1]
        self.heater = heater
        self.normal_wall = normal_wall
        self.n_rows = round(self.x_len/dx) -1 # Defines the number of rows in the coordinate mesh.
        self.n_cols = round(self.y_len/dx)    # Defines the number of columns in the coordinate mesh.  
        self.size = (self.n_rows, self.n_cols)
        self.N_elements = self.n_rows*self.n_cols # Number of points in which to calculate u 
        BC, neu_ind = self._make_boundaries()
        self.BC = BC
        self.neu_ind = neu_ind
        self.A = self._make_matrix()
        lu, piv = lu_factor(self.A)
        self.lu = lu
        self.piv = piv
        
    def _make_boundaries(self):
        
        BC_W = zeros(self.size)
        BC_E = zeros(self.size)
        BC_N = zeros(self.size)
        BC_S = zeros(self.size)
        BC_N[0,0:] = self.normal_wall
        BC_S[-1,0:] = self.normal_wall
        
        if self.interface_dir == 'west':
            BC_E[:,-1] = self.heater
            BC_W[:,0] = self.interface_vals
            neumann_ind = nonzero(BC_W.reshape(self.N_elements))
        elif self.interface_dir == 'east':
            BC_E[:,-1] = self.interface_vals
            BC_W[:,0] = self.heater
            neumann_ind = nonzero(BC_E.reshape(self.N_elements))
        BC_tot = BC_W + BC_E + BC_N + BC_S
        
        BC_tot = BC_tot.reshape(self.size[0]*self.size[1])
        return BC_tot, neumann_ind
    
    def _update_boundaries(self, interface_vals):
        self.interface_vals = interface_vals
        BC, neu = self._make_boundaries()
        self.BC = BC
    
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
        return A
    
    def solve_system(self, interface_vals):
        self._update_boundaries(interface_vals)
        u = lu_solve((self.lu, self.piv, self.BC))
        mesh_vals = u.reshape(self.n_rows,self.n_cols)
        if self.interface_dir == 'east':
            interface_vals = mesh_vals[:,-1]
        elif self.interface_dir == 'west':
            interface_vals = mesh_vals[:,0]
            
        return u, interface_vals
        
        
    
    
if __name__ == '__main__':
    s = smallRoomHeatSolver('east', 1/4, array([20,20,20]))
    #BC, neumann_ind = s._make_boundaries()
    A=s._make_matrix()
        