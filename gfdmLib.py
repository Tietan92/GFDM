# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 17:05:07 2022

@author: Stefan
"""

import numpy as np


###############################################################################

def readMesh(filename):
    
    physicalNames = []
    Nodes = []
    innerElements = []
    boundaryElements = []
      
    mesh = open(filename,'r')  
    loop = True 
 
    while loop == True:    
        line = mesh.readline()
    
        if line[:-1] == '$PhysicalNames':

            line = mesh.readline()
            numPhysicalNames = int(line[:-1])
            for i in range(0,numPhysicalNames):
                line = mesh.readline()
                values = line[:-1].split(' ')
                physicalNames.append([int(values[0]),int(values[1]), values[2][1:-1]])
            
        elif line[:-1] == '$Nodes':
            line = mesh.readline()
            numNodes = int(line[:-1])
            Nodes = np.zeros([numNodes,3])
            for i in range(0,numNodes):
                line = mesh.readline()
                values = line[:-1].split(' ')
                Nodes[i,:] = [float(values[1]), float(values[2]), float(values[3])]

        
        elif line[:-1] == '$Elements':
                line = mesh.readline()
                numElements = int(line[:-1])
                for i in range(0,numElements): 
                    line = mesh.readline()
                    values = line[:-1].split(' ')
                    nodeTags = []
                    for j in range(0,len(values)-5):
                        nodeTags.append(int(values[5+j])-1)
                    element = [int(values[3]), nodeTags]
                    if int(values[1]) == 1:
                        boundaryElements.append(element)
                    else:
                        innerElements.append(element)
                    
        elif line[:-1] == '$EndElements':
                loop = False
                
    return physicalNames, Nodes, innerElements, boundaryElements


###############################################################################

def computeFaceNormals(Nodes, innerElements, boundaryElements):
    
    faceNormals = []
    N =len(innerElements)
    
    for i in range(0,len(boundaryElements)):

        
        boundaryNodes = boundaryElements[i][1]

        for j in range(0,N):
            
            Nodesj = innerElements[j][1]
       

            difference = set(Nodesj) - set(boundaryNodes)


            if len(difference) == 1:
                
                
                n3 = list(difference)[0]
                
                v1 = Nodes[boundaryNodes[0],:2]
                v2 = Nodes[boundaryNodes[1],:2]
                v3 = Nodes[n3,:2]
                
                v12 = v1 - v2
                v13 = v3 - v1
                v23 = v3 - v2
                coord = v2 + 0.5 *v12
                
                normal = np.cross([v12[0], v12[1], 0],[0, 0, 1])[:-1]
                
                normal = normal/np.linalg.norm(normal) 
                
                if np.dot(normal,v13) > 0:
                    if np.dot(normal,v23) > 0:
                        normal = -normal
                        
                faceNormals.append([[coord[0],coord[1]],[normal[0],normal[1]]])
    
                break
        
    return faceNormals

#################################################################################
            

def computeSupport(Nodes):
    
    dGlobal = np.zeros(len(Nodes))

    for i in range(0,len(Nodes)):
        vector = Nodes[i,:2]
        diff = np.linalg.norm(Nodes[:,:2]-vector, axis = 1)
        diff[i] = np.inf
        dGlobal[i] = np.min(diff)
        
    h_rmax = np.max(dGlobal)
    h = h_rmax/0.45

    print('h = '+str(h))
    h_rmin = np.min(dGlobal)
    if h*0.2 > h_rmin:
        print('schlechte Verteilung der Punkte')
        
    support = []    

    for i in range(0,len(Nodes)):

        vector = Nodes[i,:2]
        diff = np.linalg.norm(Nodes[:,:2]-vector, axis = 1)
        diff[i] = np.inf
        comp = diff <= h
        I = np.where(comp)[0]
        supportI = I.tolist()
        supportI.insert(0,i)
        support.append(supportI)
        
    return h, support

def computeNodeTags(numOfNodes, boundaryElements):
    NodeTags = []
    for i in range(0,numOfNodes): NodeTags.append([])
    
    for i in range(0,len(boundaryElements)):
        tag = boundaryElements[i][0]
        
        I1 = boundaryElements[i][1][0]
        if tag not in NodeTags[I1]:
            NodeTags[I1].append(tag)
            
        I2 = boundaryElements[i][1][1]
        if tag not in NodeTags[I2]:
            NodeTags[I2].append(tag)
            
    for i in range(0,numOfNodes):
        if NodeTags[i] == []:
            NodeTags[i].append(0)
    
    return NodeTags

def ComputeNodeNormals(boundaryElements,faceNormals):
    NodeNormals = np.zeros([len(faceNormals),2])
    
    for i in range(0,len(faceNormals)):
        faceNormal = faceNormals[i][1]
        
        I1 = boundaryElements[i][1][0]
        NodeNormals[I1] += faceNormal
        
        I2 = boundaryElements[i][1][1]
        NodeNormals[I2] += faceNormal
        
    NodeNormalsNorm = np.array([np.linalg.norm(NodeNormals, axis = 1)])
    
    return NodeNormals/NodeNormalsNorm.transpose()
        
            
        
        
            
        
                
                



        
    
     


            
    
