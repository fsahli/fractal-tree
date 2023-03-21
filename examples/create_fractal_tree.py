# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 11:02:34 2015

@author: fsc
"""
import logging
from fractal_tree.tree import FractalTree3D
from fractal_tree.mesh import Mesh
from fractal_tree.parameters import Parameters

logging.basicConfig(level=logging.INFO)
param = Parameters()
mesh = Mesh(param.meshfile)
branches, nodes = FractalTree3D(mesh, param)
