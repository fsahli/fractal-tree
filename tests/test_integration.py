import numpy as np
from pathlib import Path
import pytest


from fractal_tree.tree import FractalTree3D
from fractal_tree.mesh import Mesh
from fractal_tree.parameters import Parameters


here = Path(__file__).absolute().parent


@pytest.fixture
def filename():
    path = here / "tmp-sphere"
    yield path.as_posix()
    for p in here.iterdir():
        if p.name.startswith(path.name):
            p.unlink()


def test_integration(filename):

    param = Parameters()
    param.meshfile = here / ".." / "examples" / "sphere.obj"
    param.filename = filename

    # Read Mesh
    mesh = Mesh(param.meshfile)

    np.random.seed(1234)
    branches, nodes = FractalTree3D(mesh, param)

    assert len(branches) == 893
    assert len(nodes.nodes) == 5822
    assert np.allclose(nodes.nodes[0], [-1, 0, 0])
    assert np.allclose(nodes.nodes[2], [-0.9990874, 0.0, 0.01975606])
    assert branches[0].nodes == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    assert branches[0].triangles == [
        3122,
        3122,
        3122,
        3122,
        3122,
        3122,
        3122,
        614,
        614,
        614,
    ]
    assert branches[0].tri == 614
