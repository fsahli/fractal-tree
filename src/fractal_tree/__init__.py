from importlib.metadata import metadata

from . import tree
from . import branch
from . import mesh
from . import viz

from .tree import generate_fractal_tree, FractalTreeParameters
from .mesh import Mesh

meta = metadata("fractal-tree")
__version__ = meta["Version"]
__author__ = meta["Author"]
__license__ = meta["License"]
__email__ = meta["Author-email"]
__program_name__ = meta["Name"]

__all__ = [
    "tree",
    "branch",
    "mesh",
    "viz",
    "generate_fractal_tree",
    "Mesh",
    "FractalTreeParameters",
]
