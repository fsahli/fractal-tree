# fractal-tree

This code is to create a fractal tree over a surface discretized by triangles. It was developed to create a representation of the Purkinje network in the ventricles of the human heart.

Read the documentation in [fractal-tree.readthedocs.org](http://fractal-tree.readthedocs.org/en/latest/).

The details of the algorithm are presented in this [article](http://www.sciencedirect.com/science/article/pii/S0021929015007332). If you are going to use this code, please cite:

> Generating Purkinje networks in the human heart. F. Sahli Costabal, D. Hurtado and E. Kuhl. Journal of Biomechanics, doi:10.1016/j.jbiomech.2015.12.025


**Pre-requisites:**

* Numpy
* Scipy

You will need .obj mesh file to create the tree. A very nice software to manipulate the mesh and export it to .obj is [MeshLab](http://meshlab.sourceforge.net). Please check if the mesh has duplicated vertex or faces before running the code. Also the orientation of the normals can change your results, because the angles will be flipped. To visualize the output, the best alternative is [Paraview](http://www.paraview.org).

To define the mesh file and the parameters of the tree to use, edit the parameters.py file and then run:

```python
from fractal_tree.tree import FractalTree3D
from fractal_tree.mesh import Mesh
from fractal_tree.parameters import Parameters

param = Parameters()
mesh = Mesh(param.meshfile)
branches, nodes = FractalTree3D(mesh, param)
```

If you have questions you can contact me at francisco.sahli  at  gmail.com
