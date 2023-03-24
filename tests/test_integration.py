import numpy as np
from pathlib import Path
import pytest
import meshio

from fractal_tree.tree import generate_fractal_tree, FractalTreeParameters
from fractal_tree.mesh import Mesh


here = Path(__file__).absolute().parent


@pytest.fixture
def filename():
    path = here / "tmp-sphere"
    yield path.as_posix()
    for p in here.iterdir():
        if p.name.startswith(path.name):
            p.unlink()


def test_integration(filename):
    param = FractalTreeParameters(
        filename=filename, second_node=np.array([-0.964, 0.00, 0.266])
    )

    # Read Mesh
    fname = here / ".." / "examples" / "sphere.obj"
    msh = meshio.read(fname)
    mesh = Mesh(verts=msh.points, connectivity=msh.cells[0].data)

    np.random.seed(1234)
    branches, nodes = generate_fractal_tree(mesh, param)

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
