# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 10:33:00 2015

@author: fsc
"""

import numpy as np

class Parameters():
    def __init__(self):
        self.meshfile='sphere.obj'
        self.filename='sphere-line'
        self.init_node=np.array([-1.0 ,0., 0.])
        self.second_node=np.array([-0.964,  0.00,  0.266      ])
        self.init_length=0.5
#Number of iterations (generations of branches)
        self.N_it=10
#Median length of the branches
        self.length=.3
        self.branch_angle=0.15
        self.w=0.1
#Length of the segments (approximately, because the lenght of the branch is random)
        self.l_segment=.01

        self.Fascicles=True
###########################################
# Fascicles data
###########################################
        self.fascicles_angles=[-1.5,.2] #rad
        self.fascicles_length=[.5,.5]
# Save data?
        self.save=True