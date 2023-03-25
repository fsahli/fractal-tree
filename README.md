# fractal-tree

This code is to create a fractal tree over a surface discretized by triangles. It was developed to create a representation of the Purkinje network in the ventricles of the human heart.

Note that this is a rewrite of the original code found at https://github.com/fsahli/fractal-tree

The details of the algorithm are presented in this [article](http://www.sciencedirect.com/science/article/pii/S0021929015007332). If you are going to use this code, please cite:

> Generating Purkinje networks in the human heart. F. Sahli Costabal, D. Hurtado and E. Kuhl. Journal of Biomechanics, doi:10.1016/j.jbiomech.2015.12.025

- Source code: https://github.com/finsberg/fractal-tree
- Documentation: https://github.com/finsberg/fractal-tree

## Install
You can install the library with pip
```
python3 -m pip install fractal-tree
```
Note that you also need a way to load the mesh from e.g gmsh or another meshing tool. For this we recommend to use [`meshio`](https://github.com/nschloe/meshio) as it support the most common formats.

## Getting started

The following illustrates a minimal example, assuming that you have surface mesh called `sphere.obj` in your current directory.

```python
import meshio
import numpy as np
from fractal_tree import generate_fractal_tree, FractalTreeParameters, Mesh

msh = meshio.read("sphere.obj")
mesh = Mesh(verts=msh.points, connectivity=msh.cells[0].data)
param = FractalTreeParameters(
    filename="sphere-line",
    N_it=10,
)
branches, nodes = generate_fractal_tree(mesh, param)
```

For a more elaborate example you can checkout the [gmsh example](https://github.com/finsberg/fractal-tree/examples/demo_gmsh.html).


## License
MIT

## Need help or having issues
Please submit an [issue](https://github.com/finsberg/fractal-tree/issues)
