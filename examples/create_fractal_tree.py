# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 11:02:34 2015

@author: fsc
"""
import logging
import meshio
from fractal_tree import generate_fractal_tree, FractalTreeParameters, Mesh

logging.basicConfig(level=logging.INFO)
fname = "sphere.obj"
msh = meshio.read(fname)
mesh = Mesh(verts=msh.points, connectivity=msh.cells[0].data)

param = FractalTreeParameters(
    filename="sphere-line3",
    N_it=10,
    second_node=mesh.verts[10, :],
)

branches, nodes = generate_fractal_tree(mesh, param)
