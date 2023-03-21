# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 11:02:34 2015

@author: fsc
"""
import logging
from fractal_tree.FractalTree import Fractal_Tree_3D
from fractal_tree.parameters import Parameters

logging.basicConfig(level=logging.INFO)
param = Parameters()
branches, nodes = Fractal_Tree_3D(param)
