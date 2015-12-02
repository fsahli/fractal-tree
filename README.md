# fractal-tree

This code is to create a fractal tree over a surface discretized by triangles. It was developed to create a representation of the Purkinje network in the ventricles of the human heart. 

The details of the algorithm are presented in this [article](http://biomechanics.stanford.edu/paper/JBIOM16.pdf). If you are going to use this code, please cite:

	Generating Purkinje networks in the human heart.
	F. Sahli Costabal, D. Hurtado and E. Kuhl.
	Journal of Biomechanics, accepted for publication.


**Pre-requisites:**

* Numpy
* Scipy
* Mayavi, if you want to export Paraview files for visualization.

You will need .obj mesh file to create the tree. A very nice software to manipulate the mesh and export it to .obj is [MeshLab](http://meshlab.sourceforge.net). Please check if the mesh has duplicated vertex or faces before running the code. Also the orientation of the normals can change your results, because the angles will be fliped. To visualize the output, the best alternative is [Paraview](http://www.paraview.org).

To define the mesh file and the parameters of the tree to use, edit the parameters.py file and then run:

```
	from FractalTree import *
	from parameters import Parameters

	param=Parameters()

	branches, nodes = Fractal_Tree_3D(param)
```

If you have questions you can contact me at francisco.sahli  at  gmail.com
