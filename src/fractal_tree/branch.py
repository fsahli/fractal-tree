"""
This module contains the Branch class (one branch of the tree)
and the Nodes class name
"""

import numpy as np
import logging
from scipy.spatial import cKDTree

logger = logging.getLogger(__name__)


class Branch:
    """Class that contains a branch of the fractal tree

    Args:
        mesh:
            an object of the mesh class, where the fractal tree will grow
        initial_node (int):
            initial node to grow the branch. This is an index that refers
            to a node in the nodes.nodes array.
        initial_direction (array):
            initial direction to grow the branch. In general, it refers to
            the direction of the last segment of the mother brach.
        initial_triangle (int):
            the index of triangle of the mesh where the initial_node sits.
        l (float):
            total length of the branch
        angle (float):
            angle (rad) with respect to the initial_direction
            in the plane of the initial_triangle triangle
        repulsitivity (float):
            repulsitivity parameter. Controls how much the branches repel each other.
        nodes:
            the object of the class nodes that contains all the
            nodes of the existing branches.
        brother_nodes (list):
            the nodes of the brother and mother branches, to be excluded
            from the collision detection between branches.
        num_segments (int):
            number of segments to divide the branch.


    Attributes:
        child (list):
            contains the indexes of the child branches.
            It is not assigned when created.
        dir (array):
            vector direction of the last segment of the branch.
        nodes (list):
            contains the node indices of the branch. The node coordinates can
            be retrieved using nodes.nodes[i]
        triangles (list):
            contains the indices of the triangles from the mesh where every
            node of the branch lies.
        tri (int):
            triangle index where last node sits.
        growing (bool):
            False if the branch collide or is out of the surface. True otherwise.

    """

    def __init__(
        self,
        mesh,
        initial_node,
        initial_direction,
        initial_triangle,
        length,
        angle,
        repulsitivity,
        nodes,
        brother_nodes,
        num_segments,
    ):

        self.child = [0, 0]
        self.dir = np.array([0.0, 0.0, 0.0])
        self.nodes = []
        self.triangles = []

        self.queue = []
        self.growing = True
        init_normal = mesh.normals[initial_triangle]
        nodes.update_collision_tree(brother_nodes)

        inplane = -np.cross(initial_direction, init_normal)
        dir = np.cos(angle) * initial_direction + np.sin(angle) * inplane
        dir = dir / np.linalg.norm(dir)
        self.nodes.append(initial_node)
        self.queue.append(nodes.nodes[initial_node])
        self.triangles.append(initial_triangle)
        grad = nodes.gradient(self.queue[0])
        dir = (dir + repulsitivity * grad) / np.linalg.norm(dir + repulsitivity * grad)

        for i in range(1, num_segments):
            intriangle = self.add_node_to_queue(
                mesh, self.queue[i - 1], dir * length / num_segments
            )

            if not intriangle:
                logger.info(f"Point {i} not in triangle")
                self.growing = False
                break
            collision = nodes.collision(self.queue[i])
            if collision[1] < length / 5.0:
                logger.debug(f"Collision {i}: {collision}")
                self.growing = False
                self.queue.pop()
                self.triangles.pop()
                collision[0]
                break
            grad = nodes.gradient(self.queue[i])
            normal = mesh.normals[self.triangles[i], :]
            # Project the gradient to the surface
            grad = grad - (np.dot(grad, normal)) * normal
            dir = (dir + repulsitivity * grad) / np.linalg.norm(
                dir + repulsitivity * grad
            )
        nodes_id = nodes.add_nodes(self.queue[1:])
        [self.nodes.append(x) for x in nodes_id]
        if not self.growing:
            nodes.end_nodes.append(self.nodes[-1])
        self.dir = dir

        self.tri = self.triangles[-1]

    # Uncomment the following lines for a closed network
    #   if shared_node is not -1:
    #      self.nodes.append(shared_node)

    def add_node_to_queue(self, mesh, initial_node, dir):
        """Functions that projects a node in the mesh surface
        and it to the queue is it lies in the surface.

        Args:
            mesh:
                an object of the mesh class, where the fractal tree will grow
            initial_node (array):
                vector that contains the coordinates of the
                last node added in the branch.
            dir (array):
                vector that contains the direction from the initial_node
                to the node to project.

        Return:
            success (bool):
                true if the new node is in the triangle.

        """
        point, triangle = mesh.project_new_point(initial_node + dir)

        if triangle >= 0:
            self.queue.append(point)
            self.triangles.append(triangle)
            success = True
        else:
            success = False

        return success


class Nodes:
    """A class containing the nodes of the branches plus some
    functions to compute distance related quantities.

    Args:
        initial_node (array):
            an array with the coordinates of the initial node of the first branch.

    Attributes:
        nodes (list):
            list of arrays containing the coordinates of the nodes
        last_node (int):
            last added node.
        end_nodes (list):
            a list containing the indices of all end nodes
            (nodes that are not connected) of the tree.
        tree (scipy.spatial.cKDTree):
            a k-d tree to compute the distance from any point
            to the closest node in the tree. It is updated once a branch is finished.
        collision_tree (scipy.spatial.cKDTree):
            a k-d tree to compute the distance from any point to the closest node
            in the tree, except from the brother and mother branches.
            It is used to check collision between branches.

    """

    def __init__(self, initial_node):
        self.nodes = []
        self.nodes.append(initial_node)
        self.last_node = 0
        self.end_nodes = []
        self.tree = cKDTree(self.nodes)

    def add_nodes(self, queue):
        """This function stores a list of nodes of a branch and
        returns the node indices. It also updates the tree to compute distances.

        Args:
            queue (list):
            a list of arrays containing the coordinates of the nodes of one branch.

        Returns:
            nodes_id (list):
            the indices of the added nodes.
        """
        nodes_id = []
        for point in queue:
            self.nodes.append(point)
            self.last_node += 1
            nodes_id.append(self.last_node)
        self.tree = cKDTree(self.nodes)
        return nodes_id

    def distance_from_point(self, point):
        """This function returns the distance from any
        point to the closest node in the tree.

        Args:
            point (array):
                the coordinates of the point to calculate the distance from.

        Returns:
            d (float):
            the distance between point and the closest node in the tree.
        """
        d, node = self.tree.query(point)
        #      distance=pool.map(lambda a: np.linalg.norm(a-point),self.nodes.values())
        return d

    def distance_from_node(self, node):
        """This function returns the distance from any
        node to the closest node in the tree.

        Args:
            node (int):
                the index of the node to calculate the distance from.

        Returns:
            d (float):
                the distance between specified node and the closest node in the tree.
        """
        d, node = self.tree.query(self.nodes[node])
        #     distance=pool.map(lambda a: np.linalg.norm(a-self.nodes[node]),self.nodes.values())
        return d

    def update_collision_tree(self, nodes_to_exclude):
        """This function updates the collision_tree excluding a
        list of nodes from all the nodes in the tree. If all the
        existing nodes are excluded, one distant node is added.

        Args:
            nodes_to_exclude (list):
                contains the nodes to exclude from the tree.
                Usually it should be the mother and the brother branch nodes.

        Returns:
            none
        """
        nodes = set(range(len(self.nodes)))
        nodes = nodes.difference(nodes_to_exclude)
        nodes_to_consider = [self.nodes[x] for x in nodes]
        self.nodes_to_consider_keys = [x for x in nodes]
        if len(nodes_to_consider) == 0:
            nodes_to_consider = [
                np.array([-100000000000.0, -100000000000.0, -100000000000.0])
            ]
            self.nodes_to_consider_keys = [100000000]
            logger.info("no nodes to consider")
        self.collision_tree = cKDTree(nodes_to_consider)

    def collision(self, point):
        """This function returns the distance between one point and
        the closest node in the tree and the index of the closest
        node using the collision_tree.

        Args:
            point (array):
                the coordinates of the point to calculate the distance from.

        Returns:
            collision (tuple):
                (distance to the closest node, index of the closest node)
        """
        d, node = self.collision_tree.query(point)
        collision = (self.nodes_to_consider_keys[node], d)
        return collision

    def gradient(self, point):
        """This function returns the gradient of the distance from
        the existing points of the tree from any point. It uses a
        central finite difference approximation.

        Args:
            point (array):
                the coordinates of the point to calculate the gradient of the distance from.

        Returns:
            grad (array):
                (x,y,z) components of gradient of the distance.
        """
        delta = 0.01
        dx = np.array([delta, 0, 0])
        dy = np.array([0.0, delta, 0.0])
        dz = np.array([0.0, 0.0, delta])
        distx_m = self.distance_from_point(point - dx)
        distx_p = self.distance_from_point(point + dx)
        disty_m = self.distance_from_point(point - dy)
        disty_p = self.distance_from_point(point + dy)
        distz_m = self.distance_from_point(point - dz)
        distz_p = self.distance_from_point(point + dz)
        grad = np.array(
            [
                (distx_p - distx_m) / (2 * delta),
                (disty_p - disty_m) / (2 * delta),
                (distz_p - distz_m) / (2 * delta),
            ]
        )
        return grad
