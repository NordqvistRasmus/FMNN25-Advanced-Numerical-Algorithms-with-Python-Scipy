# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 10:35:23 2019
@author: pontusnordqvist
"""
from  scipy import *
from  pylab import *

import matplotlib.pyplot as plt
import seaborn as sns; sns.set()

class Problem:
    
    def  __init__(self, dx, wall = 15, heater = 40, window = 5):
        self.dx = dx
        self.wall = wall
        self.heater = heater
        self.window = window
        room1 = [1,1]
        room2 = [1,2]
        room3 = [1,1]
        self.geometry = {'room1': room1, 'room2': room2, 'room3': room3}
        self.boundary = {'room1': [wall, wall, wall, heater],
                         'room2': [heater, wall, window, wall],
                         'room3': [wall, heater, wall, wall]}
        
        for r in self.geometry.keys():
            self.geometry[r]= self.init_room(dx, self.geometry[r], self.boundary[r])
    
    def init_room(self, dx, dimensions, bound = [15,15,15,15]):
        x = arange(0, dimensions[1],dx) 
        y = arange(0, dimensions[0],dx)
        room = zeros((len(x), len(y)))
        for i in range(room.shape[1]):
            room[0,i] = bound[0]
        for i in range(room.shape[0]):
            room[i,-1] = bound[1]
        for i in range(room.shape[1]):
            room[-1,i] = bound[2]
        for i in range(room.shape[0]):
            room[i,0] = bound[3]
             
        return room 
    
    def plot(self):
        fig, ax = plt.subplots()

        upper_left = zeros((self.geometry['room1'].shape[0],self.geometry['room1'].shape[1]))
        lower_right = zeros((self.geometry['room3'].shape[0],self.geometry['room3'].shape[1]))
        upper_left.fill(None)
        lower_right.fill(None)

        splitted_room2 = vsplit(array(self.geometry['room2']), 2)
        total = block([[upper_left, splitted_room2[0], self.geometry['room3']],
                       [self.geometry['room1'], splitted_room2[1], lower_right]])
        
        ax = sns.heatmap(total, cmap = "YlOrRd")
        ax.invert_yaxis()
        plt.show()
        
if __name__ == '__main__':
    p = Problem(0.1)
    p.plot()
    #print('room1',p.geometry['room1'])
    #print('room2',p.geometry['room2'])
    #print('room3',p.geometry['room3'])
    #for i in p.geometry:
    #    print(p.ge)

        