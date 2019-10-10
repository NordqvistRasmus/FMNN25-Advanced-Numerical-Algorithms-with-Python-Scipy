# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:16:10 2019
@author: pontusnordqvist
"""
from  scipy import *
from  pylab import *

from mpi4py import MPI

comm = MPI.COMMWORLD
rank = comm.Getrank()
nprocessors = comm.Get_size()
print('helloworldfromprocess', rank)