# -*- coding: utf-8 -*-
"""
Created on Fri May  8 18:33:11 2015

@author: fsc
"""

import numpy as np
from multiprocessing.dummy import Pool as ThreadPool
from operator import itemgetter
from scipy.spatial import cKDTree

pool = ThreadPool(16) 

class Branch:
    def __init__(self,mesh,init_node,init_dir,init_tri,l,angle,w,nodes,brother_nodes,Nsegments):
        self.nnodes=0
        self.child = [0,0]
        self.dir = np.array([0.0,0.0,0.0])
        self.nodes=[]
        self.triangles=[]
        self.normal=np.array([0.0,0.0,0.0])
        self.queue=[]
        self.growing=True
        shared_node=-1
        init_normal=mesh.normals[init_tri]
        nodes.update_collision_tree(brother_nodes)
        global_nnodes=len(nodes.nodes)
        
      #  R=np.array([[np.cos(angle),-np.sin(angle)],[ np.sin(angle), np.cos(angle)]])
        inplane=-np.cross(init_dir,init_normal)
        dir=np.cos(angle)*init_dir+np.sin(angle)*inplane
        dir=dir/np.linalg.norm(dir)
        self.nodes.append(init_node)
        self.queue.append(nodes.nodes[init_node])
        self.triangles.append(init_tri)
        grad=nodes.gradient(self.queue[0])
        dir=(dir+w*grad)/np.linalg.norm(dir+w*grad)
    #    print nodes.nodes[init_node]+dir*l/Nsegments
        for i in range(1,Nsegments):
            intriangle=self.add_node_to_queue(mesh,self.queue[i-1],dir*l/Nsegments)
            #print 'intriangle',intriangle
            if not intriangle:
                print 'Point not in triangle',i
#                print self.queue[i-1]+dir*l/50.
                self.growing=False
                break
            collision=nodes.collision(self.queue[i])
            if collision[1]<l/5.:
                print "Collision",i, collision
                self.growing=False
                self.queue.pop()
                self.triangles.pop()
                shared_node=collision[0]
                break
            grad=nodes.gradient(self.queue[i])
            normal=mesh.normals[self.triangles[i],:]
            #Project the gradient to the surface
            grad=grad-(np.dot(grad,normal))*normal
            dir=(dir+w*grad)/np.linalg.norm(dir+w*grad)
        nodes_id=nodes.add_nodes(self.queue[1:])
        [self.nodes.append(x) for x in nodes_id]
        if not self.growing:            
            nodes.end_nodes.append(self.nodes[-1])
        self.dir=dir
       # #print self.triangles
        self.tri=self.triangles[-1]
    #Uncomment the following lines for a closed network
     #   if shared_node is not -1:
      #      self.nodes.append(shared_node)
        
    def add_node_to_queue(self,mesh,init_node,dir):
       # print 'node trying to project', init_node+dir
        point, triangle=mesh.project_new_point(init_node+dir)
       # print 'Projected point', point, 'dist', np.linalg.norm(point-init_node)
        if triangle>=0:
            self.queue.append(point)
            self.triangles.append(triangle)
            success=True
        else:
#            print point, triangle
            success=False
        #print 'Success? ',success
        return success

class Nodes:
    def __init__(self,init_node):
        self.nodes=[]
        self.nodes.append(init_node)
        self.last_node=0
        self.end_nodes=[]
        self.tree=cKDTree(self.nodes)
    def add_nodes(self,queue):
        nodes_id=[]
        for point in queue:
            self.nodes.append(point)
            self.last_node+=1
            nodes_id.append(self.last_node)
        self.tree=cKDTree(self.nodes)
        return nodes_id
    def distance_from_point(self,point):
        d,node=self.tree.query(point)
  #      distance=pool.map(lambda a: np.linalg.norm(a-point),self.nodes.values())
        return d
    def distance_from_node(self,node):
        d, node = self.tree.query(self.nodes[node])
   #     distance=pool.map(lambda a: np.linalg.norm(a-self.nodes[node]),self.nodes.values())
        return d
    def update_collision_tree(self,nodes_to_exclude):
        nodes=set(range(len(self.nodes)))
        nodes=nodes.difference(nodes_to_exclude)
        nodes_to_consider=[self.nodes[x] for x in nodes]
        self.nodes_to_consider_keys=[x for x in nodes]
        if len(nodes_to_consider)==0:
            nodes_to_consider=[np.array([-100000000000.0,-100000000000.0,-100000000000.0])]
            self.nodes_to_consider_keys=[100000000]
            print "no nodes to consider"
        self.collision_tree=cKDTree(nodes_to_consider)
    def collision(self,point):
        d,node=self.collision_tree.query(point)
        collision=(self.nodes_to_consider_keys[node],d)
        return collision
    def gradient(self,point):
        delta=0.01
        dx=np.array([delta,0,0])
        dy=np.array([0.0,delta,0.0])
        dz=np.array([0.0,0.0,delta])
        distx_m=self.distance_from_point(point-dx)
        distx_p=self.distance_from_point(point+dx)
        disty_m=self.distance_from_point(point-dy)
        disty_p=self.distance_from_point(point+dy)
        distz_m=self.distance_from_point(point-dz)
        distz_p=self.distance_from_point(point+dz)
        grad=np.array([(distx_p-distx_m)/(2*delta),(disty_p-disty_m)/(2*delta),(distz_p-distz_m)/(2*delta)])
        return grad
