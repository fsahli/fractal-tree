# -*- coding: utf-8 -*-
"""
This module contains the mesh class. This class is the triangular surface where the fractal tree is grown. 
"""
import numpy as np
from scipy.spatial import cKDTree
import collections

class Mesh:
    """Class that contains the mesh where fractal tree is grown. It must be Wavefront .obj file. Be careful on how the normals are defined. It can change where an specified angle will go.
    
    
    Args:    
        filename (str): the path and filename of the .obj file with the mesh.

    Attributes:
        verts (array): a numpy array that contains all the nodes of the mesh. verts[i,j], where i is the node index and j=[0,1,2] is the coordinate (x,y,z).
        connectivity (array): a numpy array that contains all the connectivity of the triangles of the mesh. connectivity[i,j], where i is the triangle index and j=[0,1,2] is node index.
        normals (array): a numpy array that contains all the normals of the triangles of the mesh. normals[i,j], where i is the triangle index and j=[0,1,2] is normal coordinate (x,y,z).
        node_to_tri (dict): a dictionary that relates a node to the triangles that it is connected. It is the inverse relation of connectivity. The triangles are stored as a list for each node.
        tree (scipy.spatial.cKDTree): a k-d tree to compute the distance from any point to the closest node in the mesh.
        
    """
    def __init__(self,filename):
        verts, connectivity = self.loadOBJ(filename)
        self.verts=np.array(verts)
        self.connectivity=np.array(connectivity)
        self.normals=np.zeros(self.connectivity.shape)
        self.node_to_tri=collections.defaultdict(list)
        for i in range(len(self.connectivity)):
            for j in range(3):
                self.node_to_tri[self.connectivity[i,j]].append(i)
            u=self.verts[self.connectivity[i,1],:]-self.verts[self.connectivity[i,0],:]
            v=self.verts[self.connectivity[i,2],:]-self.verts[self.connectivity[i,0],:]
            n=np.cross(u,v)
            self.normals[i,:]=n/np.linalg.norm(n)
        self.tree=cKDTree(verts)
        
    def loadOBJ(self,filename):  
        """This function reads a .obj mesh file
        
        Args:
            filename (str): the path and filename of the .obj file.
            
        Returns:
             verts (array): a numpy array that contains all the nodes of the mesh. verts[i,j], where i is the node index and j=[0,1,2] is the coordinate (x,y,z).
             connectivity (array): a numpy array that contains all the connectivity of the triangles of the mesh. connectivity[i,j], where i is the triangle index and j=[0,1,2] is node index.
        """
        numVerts = 0  
        verts = []  
        norms = []   
        connectivity=[]
        for line in open(filename, "r"):  
            vals = line.split()
            if len(vals)>0:
                if vals[0] == "v":  
                    v = map(float, vals[1:4])  
                    verts.append(v)  
                if vals[0] == "vn":  
                    n = map(float, vals[1:4])  
                    norms.append(n)  
                if vals[0] == "f": 
                    con=[]
                    for f in vals[1:]:  
                        w = f.split("/")  
  #                      print w
                        # OBJ Files are 1-indexed so we must subtract 1 below  
                        con.append(int(w[0])-1)
                        numVerts += 1  
                    connectivity.append(con)
        return verts, connectivity
    

        
    def project_new_point(self,point):
        """This function projects any point to the surface defined by the mesh.
        
        Args:
            point (array): coordinates of the point to project.
            
        Returns:
             projected_point (array): the coordinates of the projected point that lies in the surface.
             intriangle (int): the index of the triangle where the projected point lies. If the point is outside surface, intriangle=-1.
        """
        #Get the closest point
        d, node=self.tree.query(point)
        #print d, node
        #Get triangles connected to that node
        triangles=self.node_to_tri[node]
        #print triangles
        #Compute the vertex normal as the avergage of the triangle normals.
        vertex_normal=np.sum(self.normals[triangles],axis=0)
        #Normalize
        vertex_normal=vertex_normal/np.linalg.norm(vertex_normal)
        #Project to the point to the closest vertex plane
        pre_projected_point=point-vertex_normal*np.dot(point-self.verts[node],vertex_normal)
        #Calculate the distance from point to plane (Closest point projection)
        CPP=[]
        for tri in triangles:
            CPP.append(np.dot(pre_projected_point-self.verts[self.connectivity[tri,0],:],self.normals[tri,:]))
        CPP=np.array(CPP)
     #   print 'CPP=',CPP
        triangles=np.array(triangles)
        #Sort from closest to furthest
        order=np.abs(CPP).argsort()
       # print CPP[order]
        #Check if point is in triangle
        intriangle=-1
        for o in order:
            i=triangles[o]
      #      print i
            projected_point=(pre_projected_point-CPP[o]*self.normals[i,:])
      #      print projected_point
            u=self.verts[self.connectivity[i,1],:]-self.verts[self.connectivity[i,0],:]
            v=self.verts[self.connectivity[i,2],:]-self.verts[self.connectivity[i,0],:]
            w=projected_point-self.verts[self.connectivity[i,0],:]
       #     print 'check ortogonality',np.dot(w,self.normals[i,:])
            vxw=np.cross(v,w)
            vxu=np.cross(v,u)
            uxw=np.cross(u,w)
            sign_r=np.dot(vxw,vxu)
            sign_t=np.dot(uxw,-vxu)
        #    print sign_r,sign_t            
            if sign_r>=0 and sign_t>=0:
                r=np.linalg.norm(vxw)/np.linalg.norm(vxu)
                t=np.linalg.norm(uxw)/np.linalg.norm(vxu)
             #   print 'sign ok', r , t
                if r<=1 and t<=1 and (r+t)<=1.001:
              #      print 'in triangle',i
                    intriangle = i
                    break
        return projected_point, intriangle
                
                
            
        
        
        
        
    
