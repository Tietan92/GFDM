#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 14:39:51 2022

@author: stefan
"""
import gfdmLib as gl
import numpy as np 
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import bicgstab
from scipy.sparse.linalg import spsolve

class model:
    
    def __init__(self,mshFile, f, tend):
        
        self.f = f          # Simulationsfrequenz in Hz
        self.tend = tend     # Simulationszeit in s
        self.t = np.arange(0,self.tend,1/self.f) 
        self.v = 0.3            # Querkontraktionszahl
        self.E = 210 * 10**9    # E-Modul [N/m]
        self.rho = 7850         # Dichte [kg/m³]

        self.physicalNames, self.Nodes, self.innerElements, self.boundaryElements = gl.readMesh(mshFile)
        self.faceNormals = gl.computeFaceNormals(self.Nodes, self.innerElements, self.boundaryElements)
        self.h, self.support = gl.computeSupport(self.Nodes)
        self.normals = gl.ComputeNodeNormals(self.boundaryElements, self.faceNormals)
        self.N = self.Nodes.shape[0]
        self.NodeTags = gl.computeNodeTags(self.N, self.boundaryElements)
        
        self.Px = np.zeros([self.N,self.t.size])
        self.Py = np.zeros([self.N,self.t.size])
        self.barU = np.zeros([self.N,self.t.size])
        self.barV = np.zeros([self.N,self.t.size])
        
        self.boundaryType = np.zeros(self.N)
        
        self.Y=np.zeros([2*self.N,self.t.size])
        
        
    def setBoundary(self, boundary_old, boundary_new):
        for i in range(0,self.N):
            if set(self.NodeTags[i]) == set(boundary_old):
                self.NodeTags[i] = boundary_new
        
    def defineBoundaryCondition(self, physicalTag, boundaryType, Value):
        for i in range(0,self.N):
            if self.NodeTags[i] == physicalTag:
                if boundaryType == 1:
                    self.barU[i] = Value[0]
                    self.barV[i] = Value[1]
                    self.boundaryType[i] = 1
                elif boundaryType == 2:
                    self.Px[i] = Value[0]
                    self.Py[i] = Value[1]
                    self.boundaryType[i] = 2
                    
    def solve(self):
        
        data = []
        #dataUx = []
        #dataUy = []
        #dataVx = []
        #dataVy = []
        row_ind = []
        col_ind = []
        #row_ind_spat = []
        #col_ind_spat = []

        k01  = np.zeros(2*self.N)
        k02 =  np.zeros(2*self.N)
        k11  = np.zeros(2*self.N)
        k12 =  np.zeros(2*self.N)
        k21  = np.zeros(2*self.N)
        k22 =  np.zeros(2*self.N)

        
        lame1 = self.v/(1-2*self.v)*1/(1+self.v)*self.E
        lame2 = 1/2*1/(1+self.v)*self.E
        

        
        a0 = (self.rho*self.f**2) / (lame1 + 2 * lame2)
        a1 = (lame2) / (lame1 + 2 * lame2)
        a2 = (lame1 + lame2) / (lame1 + 2 * lame2)
        a3 = (lame1) / (lame1 + 2 * lame2)
        a4 = (1) / (lame1 + 2 * lame2)
        
        

        
        G_omega = np.zeros([2,12])
        G_omega[0,0]  = a0
        G_omega[0,3]  = -1 # u_xx
        G_omega[0,5]  = -a1         # u_yy
        G_omega[0,10] = -a2   # v_xy
        G_omega[1,4]  = -a2   # u_xy
        G_omega[1,6]  = a0
        G_omega[1,9]  = -a1           # v_xx
        G_omega[1,11] = -1  # v_yy
        
        G_omega_d = np.zeros([2,12])
        G_omega_d[0,0] = 1
        G_omega_d[1,6] = 1
        
        G_omega_n = np.zeros([2,12])
        
        
        for i in range(0,self.N):
            
            if self.boundaryType[i] == 1:
                
                data.append(1) 
                row_ind.append(i)
                col_ind.append(i)
                
                data.append(1) 
                row_ind.append(i+self.N)
                col_ind.append(i+self.N)
                  
                k11[i]        = 1
                k12[i+self.N] = 1
                
                
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
                        
                    Wd[j,j] = (np.exp(-6*np.linalg.norm(diff)**2/(2*self.h**2)))**2
                
                D_z = np.zeros([2*n,12])
                D_z[:n,:6]         = D
                D_z[n:2*n,6:12]    = D
                 
                W_z = np.zeros([2*n,2*n])
                W_z[:n,:n] = Wd
                W_z[n:2*n,n:2*n] = Wd
                
                M_omega = np.append(D_z,G_omega, axis = 0)
                
                W_omega = np.zeros([2*n+2,2*n+2])
                W_omega[:2*n,:2*n] = W_z

                W_omega[-2,-2] = 1 # Wichtung u
                W_omega[-1,-1] = 1 # Wichtung v
                
                if self.boundaryType[i] == 0:
                        
                    M = M_omega
                    W = W_omega
                
                    c = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(M),W),M)),np.transpose(M)),W)
                
                    k01[i]        = c[0,2*n] * a0
                    k01[i+self.N] = c[6,2*n] * a0
                    k02[i]        = c[0,2*n+1] * a0
                    k02[i+self.N] = c[6,2*n+1] * a0
                
                
                elif self.boundaryType[i] == 1: # Dirichlet Randbedingugen
            
                    M_omega_d = np.append(M_omega,G_omega_d, axis = 0)
                    W_omega_d = np.zeros([2*n+4,2*n+4])
                    W_omega_d[:2*n+2,:2*n+2] = W_omega
            
                    W_omega_d[-2,-2] = 1 # Wichtung bar_U
                    W_omega_d[-1,-1] = 1 # Wichtung bar_v
                
                    M = M_omega_d
                    W = W_omega_d
                    
                    c = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(M),W),M)),np.transpose(M)),W)
                
                    k01[i]        = c[0,2*n] * a0
                    k01[i+self.N] = c[6,2*n] * a0
                    k02[i]        = c[0,2*n+1] * a0
                    k02[i+self.N] = c[6,2*n+1] * a0
                
                    k11[i]        = c[0,2*n+2] 
                    k11[i+self.N] = c[6,2*n+2] 
                    k12[i]        = c[0,2*n+3] 
                    k12[i+self.N] = c[6,2*n+3] 
                
                        
                elif self.boundaryType[i] == 2: # Neumann Randbedingugen
            
                    nx = self.normals[i][0]
                    ny = self.normals[i][1]
            
                    G_omega_n[0,1] = nx   # u_x
                    G_omega_n[0,2] = a1*ny              # u_y
                    G_omega_n[0,7] = a1*ny                # v_x
                    G_omega_n[0,8] = a3*nx                # v_y
                
                    G_omega_n[1,1] = a3*ny                # u_x
                    G_omega_n[1,2] = a1*nx                 # u_y
                    G_omega_n[1,7] = a1*nx                # v_x
                    G_omega_n[1,8] = ny   # v_y
                
                    M_omega_n = np.append(M_omega,G_omega_n, axis = 0)
                    W_omega_n = np.zeros([2*n+4,2*n+4])
                    W_omega_n[:2*n+2,:2*n+2] = W_omega
            
                    W_omega_n[-2,-2] = 1 # Wichtung p_x
                    W_omega_n[-1,-1] = 1 # Wichtung p_y
                    
                    M = M_omega_n
                    W = W_omega_n
                    
                    c = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(M),W),M)),np.transpose(M)),W)
                
                    k01[i]        = c[0,2*n] * a0
                    k01[i+self.N] = c[6,2*n] * a0
                    k02[i]        = c[0,2*n+1] * a0
                    k02[i+self.N] = c[6,2*n+1] * a0
                
                    k21[i]        = c[0,2*n] * a4
                    k21[i+self.N] = c[6,2*n] * a4
                    k22[i]        = c[0,2*n+1] * a4
                    k22[i+self.N] = c[6,2*n+1] * a4
                
                data.append(1- c[0,0])    # 1
                data += list(-c[0,1:n])   # n-1
                row_ind += [i]*n          # n 
                col_ind += sup            # n
                
                data += list(-c[0,n:2*n])
                row_ind += [i]*n
                col_ind += list(np.array(sup)+self.N)
                
                data.append(1- c[6,n])
                data += list(-c[6,n+1:2*n])
                row_ind += [i+self.N]*n
                col_ind += list(np.array(sup)+self.N)
                
                data += list(-c[6,:n])
                row_ind += [i+self.N]*n
                col_ind += sup
            
                #dataUx += list(c[1,:n])
                #dataUy += list(c[2,:n])
                #dataVx += list(c[7,:n])
                #dataVy += list(c[8,:n])
                #row_ind_spat += [i]*n
                #col_ind_spat += sup   
            
                #dataUx += list(c[1,n:2*n])
                #dataUy += list(c[2,n:2*n])
                #dataVx += list(c[7,n:2*n])
                #dataVy += list(c[8,n:2*n])
                #row_ind_spat += [i]*n
                #col_ind_spat += list(np.array(sup)+self.N)
                
        A = csc_matrix((data, (row_ind, col_ind)),shape=(2*self.N,2*self.N))

        self.coordsNewX = np.zeros([self.N,self.t.size-2])
        self.coordsNewY = np.zeros([self.N,self.t.size-2])


        
        for k in range(2,self.t.size):
            
            uk1 = np.append(self.Y[:self.N,k-1],self.Y[:self.N,k-1])
            uk2 = np.append(self.Y[:self.N,k-2],self.Y[:self.N,k-2])
            vk1 = np.append(self.Y[self.N:,k-1],self.Y[self.N:,k-1])
            vk2 = np.append(self.Y[self.N:,k-2],self.Y[self.N:,k-2])
            barUk = np.append(self.barU[:,k],self.barU[:,k])
            barVk = np.append(self.barV[:,k],self.barV[:,k])
            Pxk   =  np.append(self.Px[:,k],self.Px[:,k])
            Pyk   =  np.append(self.Py[:,k],self.Py[:,k])
            
            B = k01 * (2 * uk1 - uk2) + k02 * (2* vk1 - vk2) + k11 * barUk + k12 * barVk
            
            #self.Y[:,k] = bicgstab(A, B, x0=np.copy(self.Y[:,k-1]),tol=1e-50)[0] 
            self.Y[:,k] = spsolve(A,B)
            self.coordsNewX[:,k-2] = self.Y[:self.N,k]+self.Nodes[:,0]
            self.coordsNewY[:,k-2] = self.Y[self.N:,k]+self.Nodes[:,1]
            
            print(k)
            
            
                    
