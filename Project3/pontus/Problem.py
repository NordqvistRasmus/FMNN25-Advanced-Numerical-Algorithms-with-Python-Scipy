# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 10:35:23 2019
@author: pontusnordqvist
"""
from  scipy import *
from  pylab import *


class Problem:
    
   def  __init__(self, dx, room1 = [1,1], room2 = [1,2], room3 = [1,1]):
        wall = 15
        heat = 40
        wind = 5
        
        self.geometry = {'room1': room1, 'room2': room2, 'room3': room3}
        self.boundary = {'room1': [wall, interface, wall, heater],
                         'room2': [heater, ]}
        
        for r in self.geometry.keys():
            self.geometry[r]= self.init_room(dx, self.geometry[r])
        
   def init_room(self, dx, dimensions, north = 15, east = 15, south = 15, west = 15):
        x = arange(0, dimensions[0],dx) 
        y = arange(0, dimensions[1],dx)
        room = zeros((len(x), len(y)))
        for i in range(room.shape[1]):
            room[0,i] = north
        for i in range(room.shape[0]):
            room[i,-1] = east
        for i in range(room.shape[1]):
            room[-1,i] = south
        for i in range(room.shape[0]):
            room[i,0] = west
             
        return room
        
if __name__ == '__main__':
    p = Problem(0.1)
    print('room1',p.geometry['room1'])
    print('room2',p.geometry['room2'])
    print('room3',p.geometry['room3'])
    #for i in p.geometry:
    #    print(p.ge)
        