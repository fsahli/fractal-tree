.. . documentation master file, created by
   sphinx-quickstart on Mon Nov 30 12:48:33 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Fractal-Tree's documentation!
========================================

This code is to create a fractal tree over a surface discretized by triangles. It was developed to create a representation of the Purkinje network in the ventricles of the human heart. 

The details of the algorithm are presented in this `article <https://www.sciencedirect.com/science/article/pii/S0021929015007332>`_. If you are going to use this code, please cite:

|	Generating Purkinje networks in the human heart.
|	F. Sahli Costabal, D. Hurtado and E. Kuhl.
|	Journal of Biomechanics, accepted for publication.
|

**Pre-requisites:**

* Numpy
* Scipy
* Mayavi, if you want to export Paraview files for visualization.

You will need .obj mesh file to create the tree. A very nice software to manipulate the mesh and export it to .obj is `MeshLab <http://meshlab.sourceforge.net>`_. Please check if the mesh has duplicated vertex or faces before running the code. Also the orientation of the normals can change your results, because the angles will be fliped. To visualize the output, the best alternative is `Paraview <http://www.paraview.org>`_.

To define the mesh file and the parameters of the tree to use, edit the parameters.py file and then run:

.. code-block:: python

	from FractalTree import *
	from parameters import Parameters

	param=Parameters()

	branches, nodes = Fractal_Tree_3D(param)

If you have questions you can contact me at francisco.sahli  at  gmail.com

Contents:

.. toctree::
   :maxdepth: 4

   Branch3D
   FractalTree
   Mesh
   parameters


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

