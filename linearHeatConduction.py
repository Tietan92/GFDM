# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 18:29:08 2022

@author: Stefan
"""

import gfdmLib as gl
import numpy as np 
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import bicgstab

class model:
    
    def __init__(self,f, tend, T0, mshFile):
        
        self.f = f          # Simulationsfrequenz in Hz
        self.tend = tend     # Simulationszeit in s
        self.t = np.arange(0,self.tend,1/self.f) 
        self.condy = 2000.0             # Wärmeleitfähigkeit [W/(m K)]
        self.rho    = 6450           # Dichte [kg/m^3]
        self.heatcap = 800           # spez. Wärmekapazität [J/kg K]
        self.physicalNames, self.Nodes, self.innerElements, self.boundaryElements = gl.readMesh(mshFile)
        self.faceNormals = gl.computeFaceNormals(self.Nodes, self.innerElements, self.boundaryElements)
        self.h, self.support = gl.computeSupport(self.Nodes)
        self.normals = gl.ComputeNodeNormals(self.boundaryElements, self.faceNormals)
        self.N = self.Nodes.shape[0]
        self.NodeTags = gl.computeNodeTags(self.N, self.boundaryElements)
        
        self.dirichlets  = np.zeros([self.N,self.t.size])
        self.dirichlets[:,0] = T0
        self.neumanns    = np.zeros([self.N,self.t.size])
        self.boundaryType = np.zeros(self.N)
        
    def setBoundary(self, boundary_old, boundary_new):
        for i in range(0,self.N):
            if set(self.NodeTags[i]) == set(boundary_old):
                self.NodeTags[i] = boundary_new
        
    def defineBoundaryCondition(self, physicalTag, boundaryType, Value):
        for i in range(0,self.N):
            if self.NodeTags[i] == physicalTag:
                if boundaryType == 1:
                    self.dirichlets[i] = Value 
                    self.boundaryType[i] = 1
                elif boundaryType == 2:
                    self.neumanns[i] = Value
                    self.boundaryType[i] = 2
                    
    def solve(self):
        
        data = []
        row_ind = []
        col_ind = []

        beta  = np.zeros(self.N)
        gamma = np.zeros(self.N)
        zeta  = np.zeros(self.N)
        
        a = self.condy/(self.rho*self.heatcap)
        
        G_omega = np.zeros([1,6])
        G_omega[0,0]  = self.f 
        G_omega[0,3]  = - a 
        G_omega[0,5]  = - a

        
        for i in range(0,self.N):
            
            if self.boundaryType[i] == 1:
            
                    data.append(1) 
                    row_ind.append(i)
                    col_ind.append(i)
                    zeta[i] = 1
                    
            else:
                
                sup = self.support[i]
                n = len(sup)
            
                D = np.zeros([n, 6])
                Wd = np.zeros([n,n])
                
                for j in range(0,n):
                    
                    diff =  self.Nodes[sup[j],:] - self.Nodes[i,:]
                
                    D[j,0] = 1
                    D[j,1] = diff[0] #xj - xi 
                    D[j,2] = diff[1] #yj - yi
                    D[j,3] = 0.5*diff[0]**2 #1/2*(xj - xi)**2 
                    D[j,4] = diff[0]*diff[1]#(xj - xi )*(yj - yi )
                    D[j,5] = 0.5*diff[1]**2 #*(yj - yi)**2
                
                    Wd[j,j] = (np.exp(-5*np.linalg.norm(diff)**2/(2*self.h**2)))**2
                
                W_omega = np.zeros([n+1,n+1])
                W_omega[:n,:n] = Wd
                W_omega[-1,-1] = 1 # Wichtungsfaktor Innere Gebietspunkte       
                
                M_omega = np.append(D,G_omega,axis=0)
                    
                if self.boundaryType[i] == 0:
                        
                    M = M_omega
                    W = W_omega
                
                    c = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(M),W),M)),np.transpose(M)),W)
                    beta[i]  = c[0,n]
                        
                elif self.boundaryType[i] == 2:
                        
                    W_omega_n = np.zeros([n+2,n+2])
                    W_omega_n[:n+1,:n+1] = W_omega
                    W_omega_n[-1,-1]  = 1 # Wichtungsfaktor Neumann Randbedingungen
                        
                    nx = self.normals[i][0]
                    ny = self.normals[i][1]
                        
                    G_omega_n = np.zeros([1,6])
                    G_omega_n[0,1] = -self.condy*nx
                    G_omega_n[0,2] = -self.condy*ny
                    
                    M_omega_n = np.append(M_omega,G_omega_n,axis=0)

                    M = M_omega_n
                    W = W_omega_n
                
                    c = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(M),W),M)),np.transpose(M)),W)
                    beta[i]  = c[0,n]
                    gamma[i] = c[0,n+1]
                        
                data.append(1- c[0,0])
                data += list(-c[0,1:n])
                row_ind += [i]*n
                col_ind += sup
              
        A = csc_matrix((data, (row_ind, col_ind)),shape=(self.N,self.N))
            
        Y=np.zeros([self.N,self.t.size]) # Temperatur
        Y[:,0] = self.dirichlets[:,0]
        for k in range(2,self.t.size):
            
            qa = self.neumanns[:,k]
            Td = self.dirichlets[:,k]

            B = beta * self.f * Y[:,k-1] + gamma * qa  + zeta * Td 
            
            Y[:,k] = bicgstab(A, B, x0=np.copy(Y[:,k-1]))[0]   
            print(k)
        self.Y = Y