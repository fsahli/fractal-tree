"""
This module contains the function that creates the fractal tree.
"""
from __future__ import annotations
from dataclasses import dataclass
import sys
from typing import Optional, NamedTuple, Union
import numpy as np
import logging
import tqdm
from pathlib import Path

from .branch import Nodes, Branch
from .mesh import Mesh
from .viz import write_line_VTU

logger = logging.getLogger(__name__)


def grow_fascicles(
    branches: dict[int, Branch],
    parameters: FractalTreeParameters,
    mesh: Mesh,
    nodes: Nodes,
    lines: list[tuple[int, int]],
    last_branch: int,
):
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
                (
                    branches[last_branch].nodes[i_n],
                    branches[last_branch].nodes[i_n + 1],
                )
            )
    branches_to_grow = list(range(1, len(parameters.fascicles_angles) + 1))
    return branches_to_grow, last_branch


def run_generation(
    branches_to_grow: list[int],
    parameters: FractalTreeParameters,
    branches: dict[int, Branch],
    last_branch: int,
    mesh: Mesh,
    nodes: Nodes,
    lines: list[tuple[int, int]],
):
    choices = 2 * np.random.randint(2, size=len(branches_to_grow)) - 1
    lengths = np.random.normal(0, parameters.std_length, size=2 * len(branches_to_grow))
    k = 0
    k1 = 0
    np.random.shuffle(branches_to_grow)
    new_branches_to_grow = []
    for branch_index in branches_to_grow:
        branch = branches[branch_index]
        angle = -parameters.branch_angle * choices[k]
        k += 1
        for j in range(2):
            brother_nodes = branch.nodes.copy()
            if j > 0:
                brother_nodes += branches[last_branch].nodes

            # Add new branch
            last_branch += 1
            logger.debug(last_branch)
            total_length = parameters.length + lengths[k1]
            k1 += 1
            if total_length < parameters.min_length:
                total_length = parameters.min_length
            branches[last_branch] = Branch(
                mesh=mesh,
                initial_node=branch.nodes[-1],
                initial_direction=branch.dir,
                initial_triangle=branch.tri,
                length=total_length,
                angle=angle,
                repulsitivity=parameters.repulsitivity,
                nodes=nodes,
                brother_nodes=brother_nodes,
                num_segments=int(parameters.length / parameters.l_segment),
            )
            # Add nodes to lines
            for n1, n2 in zip(
                branches[last_branch].nodes[:-1], branches[last_branch].nodes[1:]
            ):
                lines.append((n1, n2))

            # Add to the new array
            if branches[last_branch].growing:
                new_branches_to_grow.append(last_branch)

            branch.child[j] = last_branch
            angle = -angle

    branches_to_grow = new_branches_to_grow
    return branches, nodes, lines, branches_to_grow, lines, last_branch


def save_tree(
    filename: Union[Path, str],
    nodes: Nodes,
    lines: list[tuple[int, int]],
    save_paraview: bool = True,
):
    if save_paraview:
        logger.info("Finished growing, writing paraview file")
        xyz = np.zeros((len(nodes.nodes), 3))
        for i in range(len(nodes.nodes)):
            xyz[i, :] = nodes.nodes[i]
        write_line_VTU(xyz, lines, Path(filename).with_suffix(".vtu"))
    name = Path(filename).name
    np.savetxt(
        Path(filename).with_name(name + "_lines").with_suffix(".txt"), lines, fmt="%d"
    )
    np.savetxt(Path(filename).with_name(name + "_xyz").with_suffix(".txt"), xyz)
    np.savetxt(
        Path(filename).with_name(name + "_endnodes").with_suffix(".txt"),
        nodes.end_nodes,
        fmt="%d",
    )


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
            If no node is provided, a random node will be selected
        init_length (float):
            length of the first branch.
        N_it (int):
            number of generations of branches.
        length (float):
            average length of the branches in the tree.
        branch_angle (float):
            angle with respect to the direction of
            the previous branch and the new branch.
        repulsitivity (float):
            repulsitivity parameter.
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
    second_node: Optional[np.ndarray] = None
    initial_direction: Optional[np.ndarray] = None
    init_length: float = 0.1
    N_it: int = 10  # Number of iterations (generations of branches)
    length: float = 0.1  # Median length of the branches
    branch_angle: float = 0.15
    repulsitivity: float = 0.1
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


def node_direction(src: np.ndarray, target: Optional[np.ndarray] = None) -> np.ndarray:
    """Return the direction from src to target.

    Parameters
    ----------
    src : np.ndarray
        Source node
    target : Optional[np.ndarray]
        Target node

    Returns
    -------
    np.ndarray
        The unit vector from src to node 2
    """
    return (target - src) / np.linalg.norm(target - src)


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

    if parameters.initial_direction is None:
        # Get the second node to define the initial direction
        second_node = parameters.second_node
        if second_node is None:
            # If no node is specified, lets just pick a random node
            second_node = mesh.verts[np.random.choice(mesh.valid_nodes), :]

        # Define the initial direction
        initial_direction = node_direction(src=mesh.init_node, target=second_node)
    else:
        initial_direction = parameters.initial_direction

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

    lines = [
        (n1, n2)
        for n1, n2 in zip(
            branches[last_branch].nodes[:-1], branches[last_branch].nodes[1:]
        )
    ]

    # To grow fascicles
    if parameters.generate_fascicles:
        branches_to_grow, last_branch = grow_fascicles(
            branches, parameters, mesh, nodes, lines, last_branch
        )

    for _ in tqdm.tqdm(range(parameters.N_it)):
        branches, nodes, lines, branches_to_grow, lines, last_branch = run_generation(
            branches_to_grow, parameters, branches, last_branch, mesh, nodes, lines
        )

    if parameters.save:
        save_tree(
            filename=parameters.filename,
            nodes=nodes,
            lines=lines,
            save_paraview=parameters.save_paraview,
        )

    return FractalTreeResult(branches, nodes)
