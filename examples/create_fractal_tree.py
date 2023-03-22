# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 11:02:34 2015

@author: fsc
"""
import logging
from fractal_tree import generate_fractal_tree, FractalTreeParameters, Mesh

logging.basicConfig(level=logging.INFO)
mesh = Mesh.from_file("sphere.obj")
param = FractalTreeParameters(
    filename="sphere-line3",
    N_it=10,
    second_node=mesh.verts[10, :],
)

branches, nodes = generate_fractal_tree(mesh, param)
