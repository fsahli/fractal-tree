"""
This module contains the function that creates the fractal tree.
"""


from dataclasses import dataclass
import sys
from typing import Optional, NamedTuple
import numpy as np
import logging

from .branch import Nodes, Branch
from .mesh import Mesh
from .viz import write_line_VTU

logger = logging.getLogger(__name__)


def grow_fascicles(branches, parameters, mesh, nodes, lines, last_branch):
    brother_nodes = []
    brother_nodes += branches[0].nodes
    for i_branch in range(len(parameters.fascicles_angles)):
        last_branch += 1
        angle = parameters.fascicles_angles[i_branch]
        branches[last_branch] = Branch(
            mesh,
            branches[0].nodes[-1],
            branches[0].dir,
            branches[0].tri,
            parameters.fascicles_length[i_branch],
            angle,
            0.0,
            nodes,
            brother_nodes,
            int(parameters.fascicles_length[i_branch] / parameters.l_segment),
        )
        brother_nodes += branches[last_branch].nodes

        for i_n in range(len(branches[last_branch].nodes) - 1):
            lines.append(
                [
                    branches[last_branch].nodes[i_n],
                    branches[last_branch].nodes[i_n + 1],
                ]
            )
    branches_to_grow = list(range(1, len(parameters.fascicles_angles) + 1))
    return branches_to_grow, last_branch


def run_generation(
    branches_to_grow, parameters, branches, last_branch, mesh, nodes, lines
):
    choices = 2 * np.random.randint(2, size=len(branches_to_grow)) - 1
    lengths = np.random.normal(0, parameters.std_length, size=2 * len(branches_to_grow))
    k = 0
    k1 = 0
    np.random.shuffle(branches_to_grow)
    new_branches_to_grow = []
    for g in branches_to_grow:
        angle = -parameters.branch_angle * choices[k]
        k += 1
        for j in range(2):
            brother_nodes = []
            brother_nodes += branches[g].nodes
            if j > 0:
                brother_nodes += branches[last_branch].nodes

            # Add new branch
            last_branch += 1
            logger.debug(last_branch)
            l = parameters.length + lengths[k1]
            k1 += 1
            if l < parameters.min_length:
                l = parameters.min_length
            branches[last_branch] = Branch(
                mesh,
                branches[g].nodes[-1],
                branches[g].dir,
                branches[g].tri,
                l,
                angle,
                parameters.w,
                nodes,
                brother_nodes,
                int(parameters.length / parameters.l_segment),
            )
            # Add nodes to lines
            for i_n in range(len(branches[last_branch].nodes) - 1):
                lines.append(
                    [
                        branches[last_branch].nodes[i_n],
                        branches[last_branch].nodes[i_n + 1],
                    ]
                )

            # Add to the new array
            if branches[last_branch].growing:
                new_branches_to_grow.append(last_branch)

            branches[g].child[j] = last_branch
            angle = -angle
    branches_to_grow = new_branches_to_grow
    return branches, nodes, lines, branches_to_grow, lines, last_branch


def save_tree(parameters, nodes, lines):
    if parameters.save_paraview:

        logger.info("Finished growing, writing paraview file")
        xyz = np.zeros((len(nodes.nodes), 3))
        for i in range(len(nodes.nodes)):
            xyz[i, :] = nodes.nodes[i]
        write_line_VTU(xyz, lines, parameters.filename + ".vtu")

    np.savetxt(parameters.filename + "_lines.txt", lines, fmt="%d")
    np.savetxt(parameters.filename + "_xyz.txt", xyz)
    np.savetxt(parameters.filename + "_endnodes.txt", nodes.end_nodes, fmt="%d")


@dataclass
class FractalTreeParameters:
    """Class to specify the parameters of the fractal tree.

    Attributes:
        filename (str):
            name of the output files.
        init_node (numpy array):
            the first node of the tree.
        second_node (numpy array):
            this point is only used to calculate the
            initial direction of the tree and is not
            included in the tree. Please avoid selecting
            nodes that are connected to the init_node by a
            single edge in the mesh, because it causes numerical issues.
        init_length (float):
            length of the first branch.
        N_it (int):
            number of generations of branches.
        length (float):
            average length of the branches in the tree.
        branch_angle (float):
            angle with respect to the direction of
            the previous branch and the new branch.
        w (float):
            repulsivity parameter.
        l_segment (float):
            length of the segments that compose one branch
            (approximately, because the length of the branch is random).
            It can be interpreted as the element length in a finite element mesh.
        Fascicles (bool):
            include one or more straight branches with different lengths and
            angles from the initial branch. It is motivated by the fascicles
            of the left ventricle.
        fascicles_angles (list):
            angles with respect to the initial branches of the fascicles.
            Include one per fascicle to include.
        fascicles_length (list):
            length  of the fascicles. Include one per fascicle to include.
            The size must match the size of fascicles_angles.
        save (bool):
            save text files containing the nodes, the connectivity and end
            nodes of the tree.
        save_paraview (bool):
            save a .vtu paraview file. The tvtk module must be installed.

    """

    filename: str = "results"
    second_node: np.ndarray = np.array([-0.964, 0.00, 0.266])
    init_length: float = 0.1
    N_it: int = 10  # Number of iterations (generations of branches)
    length: float = 0.1  # Median length of the branches
    branch_angle: float = 0.15
    w: float = 0.1
    l_segment: float = 0.01  # Length of the segments (approximately, because
    # the length of the branch is random)
    generate_fascicles: bool = True
    fascicles_angles: tuple[float, float] = (-1.5, 0.2)  # rad
    fascicles_length: tuple[float, float] = (0.5, 0.5)
    save: bool = True
    save_paraview: bool = True

    @property
    def std_length(self) -> float:
        """Standard deviation of the length.
        Set to zero to avoid random lengths."""
        return np.sqrt(0.2) * self.length

    @property
    def min_length(self) -> float:
        """Minimum length of the branches.
        To avoid randomly generated negative lengths."""
        return self.length / 10.0

    def as_dict(self):
        return {k: v for k, v in self.__dict__.items()}


class FractalTreeResult(NamedTuple):
    branches: dict[int, Branch]
    nodes: Nodes


def generate_fractal_tree(
    mesh: Mesh, parameters: Optional[FractalTreeParameters] = None
) -> FractalTreeResult:
    """This function creates the fractal tree.

    Args:
        parameters (Optional[FractalTreeParameters]):
            This object contains all the parameters that
            define the tree. See the parameters module documentation for details:

    Returns:
        FractalTreeResult: branches containes a dictionary that contains
        all the branches objects, and nodes is the object that
        contains all the nodes of the tree.
    """
    if parameters is None:
        parameters = FractalTreeParameters()

    # Define the initial direction
    initial_direction = (parameters.second_node - mesh.init_node) / np.linalg.norm(
        parameters.second_node - mesh.init_node
    )

    # Initialize the nodes object, contains the nodes and all the distance functions
    nodes = Nodes(mesh.init_node)
    # Project the first node to the mesh.
    point = mesh.project_new_point(nodes.nodes[0])
    if point.triangle_index >= 0:
        initial_triangle = point.triangle_index
    else:
        logger.error("initial point not in mesh")
        sys.exit(0)
    # Initialize the dictionary that stores the branches objects
    branches = {}
    last_branch = 0
    # Compute the first branch
    branches[last_branch] = Branch(
        mesh=mesh,
        initial_node=0,
        initial_direction=initial_direction,
        initial_triangle=initial_triangle,
        length=parameters.init_length,
        angle=0.0,
        repulsitivity=0.0,
        nodes=nodes,
        brother_nodes=[0],
        num_segments=int(parameters.init_length / parameters.l_segment),
    )
    branches_to_grow = []
    branches_to_grow.append(last_branch)

    lines = []
    for i_n in range(len(branches[last_branch].nodes) - 1):
        lines.append(
            [branches[last_branch].nodes[i_n], branches[last_branch].nodes[i_n + 1]]
        )
    # To grow fascicles
    if parameters.generate_fascicles:
        branches_to_grow, last_branch = grow_fascicles(
            branches, parameters, mesh, nodes, lines, last_branch
        )

    for _ in range(parameters.N_it):
        branches, nodes, lines, branches_to_grow, lines, last_branch = run_generation(
            branches_to_grow, parameters, branches, last_branch, mesh, nodes, lines
        )

    if parameters.save:
        save_tree(parameters, nodes, lines)

    return FractalTreeResult(branches, nodes)
