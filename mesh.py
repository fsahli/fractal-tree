# -*- coding: utf-8 -*-
"""
Created on Thu May  7 20:16:20 2015

@author: fsc
"""
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree
import collections

class Mesh:
    def __init__(self,filename):
        verts, norms,connectivity = self.loadOBJ(filename)
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
        return verts, norms,connectivity
    
    def plot_surface(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')        
        ax.plot_trisurf(self.verts[:,0], self.verts[:,1],self.connectivity, self.verts[:,2], cmap=cm.jet, linewidth=0.2)        
        plt.show()
        
    def project_new_point(self,point):
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
                
                
            
        
        
        
        
    
