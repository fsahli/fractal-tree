"""
This module contains the mesh class. This class is the
triangular surface where the fractal tree is grown.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, NamedTuple
import numpy as np
from scipy.spatial import cKDTree
import collections


def closest_point_projection(
    triangles, pre_projected_point, verts, connectivity, normals
):
    return (
        (pre_projected_point - verts[connectivity[triangles, 0], :])
        * normals[triangles, :]
    ).sum(-1)


def get_node_to_triangle(connectivity: np.ndarray) -> dict[int, list[int]]:
    node_to_tri = collections.defaultdict(list)
    for i in range(len(connectivity)):
        for j in range(3):
            node_to_tri[connectivity[i, j]].append(i)
    return node_to_tri


def compute_normals(connectivity: np.ndarray, verts: np.ndarray) -> np.ndarray:
    U = verts[connectivity[:, 1], :] - verts[connectivity[:, 0], :]
    V = verts[connectivity[:, 2], :] - verts[connectivity[:, 0], :]
    N = np.cross(U, V)
    normals = (N.T / np.linalg.norm(N, axis=1)).T
    return normals


class ProjectedPoint(NamedTuple):
    point: np.ndarray
    triangle_index: int


class InvalidNodeError(Exception):
    pass


@dataclass
class Mesh:
    """Class that contains the mesh where fractal tree is grown.
    It must be Wavefront .obj file. Be careful on how the normals
    are defined. It can change where an specified angle will go.


    Args:
        filename (str):
            the path and filename of the .obj file with the mesh.

    Attributes:
        verts (array):
            a numpy array that contains all the nodes of the mesh.
            verts[i,j], where i is the node index and j=[0,1,2] is
            the coordinate (x,y,z).
        connectivity (array):
            a numpy array that contains all the connectivity of the
            triangles of the mesh. connectivity[i,j], where i is
            the triangle index and j=[0,1,2] is node index.
        normals (array):
            a numpy array that contains all the normals of the triangles
            of the mesh. normals[i,j], where i is the triangle index
            and j=[0,1,2] is normal coordinate (x,y,z).
        node_to_tri (dict):
            a dictionary that relates a node to the triangles
            that it is connected. It is the inverse relation of connectivity.
            The triangles are stored as a list for each node.
        tree (scipy.spatial.cKDTree):
            a k-d tree to compute the distance from any point
            to the closest node in the mesh.
    """

    verts: np.ndarray
    connectivity: np.ndarray
    init_node: Optional[np.ndarray] = None
    normals: np.ndarray = field(init=False)
    node_to_tri: dict[int, list[int]] = field(init=False)
    tree: cKDTree = field(init=False)

    def __post_init__(self):
        self.verts = np.array(self.verts)
        self.connectivity = np.array(self.connectivity)

        self.normals = compute_normals(self.connectivity, verts=self.verts)
        self.node_to_tri = get_node_to_triangle(connectivity=self.connectivity)
        self.valid_nodes = np.array(tuple(self.node_to_tri.keys()))

        self.tree = cKDTree(self.verts)

        if self.init_node is None:
            min_node = self.verts[self.valid_nodes, 0].argmin()
            self.init_node = self.verts[self.valid_nodes[min_node], :]
        self.init_node = np.array(self.init_node)

    def project_new_point(self, point: np.ndarray) -> ProjectedPoint:
        """This function projects any point to
        the surface defined by the mesh.

        Args:
            point (array):
                coordinates of the point to project.

        Returns:
            ProjectedPoint: Contains projected_point which is the coordinates
            of the projected point that lies in the surface. triangle_index is
            the index of the triangle where the projected point lies.
            If the point is outside surface, triangle_index=-1.

        """
        # Get the closest point
        node = self.tree.query(point)[1]

        # Get triangles connected to that node
        triangles = self.node_to_tri[node]
        if len(triangles) == 0:
            raise InvalidNodeError(
                f"node {node} with point {point} not connected to triangles, check your mesh"
            )

        # Compute the vertex normal as the avergage of the triangle normals.
        vertex_normal = np.sum(self.normals[triangles], axis=0)
        # Normalize
        vertex_normal = vertex_normal / np.linalg.norm(vertex_normal)
        # Project to the point to the closest vertex plane
        pre_projected_point = point - vertex_normal * np.dot(
            point - self.verts[node], vertex_normal
        )
        # Calculate the distance from point to plane (Closest point projection)
        CPP = closest_point_projection(
            triangles, pre_projected_point, self.verts, self.connectivity, self.normals
        )
        triangles = np.array(triangles)
        # Sort from closest to furthest
        order = np.abs(CPP).argsort()

        return self.check_in_triangle(order, triangles, pre_projected_point, CPP)

    def check_in_triangle(
        self, order, triangles, pre_projected_point, CPP
    ) -> ProjectedPoint:
        for o in order:
            i = triangles[o]

            projected_point = pre_projected_point - CPP[o] * self.normals[i, :]

            u = (
                self.verts[self.connectivity[i, 1], :]
                - self.verts[self.connectivity[i, 0], :]
            )
            v = (
                self.verts[self.connectivity[i, 2], :]
                - self.verts[self.connectivity[i, 0], :]
            )
            w = projected_point - self.verts[self.connectivity[i, 0], :]

            vxw = np.cross(v, w)
            vxu = np.cross(v, u)
            uxw = np.cross(u, w)
            sign_r = np.dot(vxw, vxu)
            sign_t = np.dot(uxw, -vxu)

            if sign_r >= 0 and sign_t >= 0:
                r = np.linalg.norm(vxw) / np.linalg.norm(vxu)
                t = np.linalg.norm(uxw) / np.linalg.norm(vxu)

                if r <= 1 and t <= 1 and (r + t) <= 1.001:
                    return ProjectedPoint(point=projected_point, triangle_index=i)

        return ProjectedPoint(point=projected_point, triangle_index=-1)
