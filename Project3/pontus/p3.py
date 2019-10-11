
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 20:16:48 2019
@author: Mattias Lundström1
"""
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocessors = comm.Get_size()
print ("helloworld from process", rank)
