# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 11:02:34 2015

@author: fsc
"""
import logging
from fractal_tree import generate_fractal_tree, FractalTreeParameters, Mesh

logging.basicConfig(level=logging.INFO)
param = FractalTreeParameters(filename="sphere-line")
mesh = Mesh("sphere.obj")
branches, nodes = generate_fractal_tree(mesh, param)
