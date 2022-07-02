#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 16:32:25 2022

@author: stefan
"""

import gfdmLib as gl
import numpy as np 
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import bicgstab

class model:
    
    def __init__(self,mshFile):
        
        self.v = 0.3            # Querkontraktionszahl
        self.E = 210 * 10**9    # E-Modul [N/m]

        self.physicalNames, self.Nodes, self.innerElements, self.boundaryElements = gl.readMesh(mshFile)
        self.faceNormals = gl.computeFaceNormals(self.Nodes, self.innerElements, self.boundaryElements)
        self.h, self.support = gl.computeSupport(self.Nodes)
        self.normals = gl.ComputeNodeNormals(self.boundaryElements, self.faceNormals)
        self.N = self.Nodes.shape[0]
        self.NodeTags = gl.computeNodeTags(self.N, self.boundaryElements)
        
        self.Px = np.zeros(self.N)
        self.Py = np.zeros(self.N)
        self.barU = np.zeros(self.N)
        self.barV = np.zeros(self.N)
        
        self.boundaryType = np.zeros(self.N)
        
        
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
        dataUx = []
        dataUy = []
        dataVx = []
        dataVy = []
        row_ind = []
        col_ind = []
        row_ind_spat = []
        col_ind_spat = []

        beta  = np.zeros(2*self.N)
        gamma = np.zeros(2*self.N)
        zeta  = np.zeros(2*self.N)
        
        lame1 = self.v/(1-2*self.v)*1/(1+self.v)*self.E
        lame2 = 1/2*1/(1+self.v)*self.E
        
        a0 = lame2/(lame1+2*lame2)
        a1 = (lame1+lame2)/(lame1+2*lame2)
        a2 = lame1/(lame1+2*lame2)
        a3 = 1/(lame1+2*lame2)
        
        G_omega = np.zeros([2,12])
        G_omega[0,3]  = 1  # u_xx
        G_omega[0,5]  = a0              # u_yy
        G_omega[0,10] = a1      # v_xy
        G_omega[1,4]  = a1      # u_xy
        G_omega[1,9]  = a0           # v_xx
        G_omega[1,11] = 1  # v_yy
        
        G_omega_d = np.zeros([2,12])
        G_omega_d[0,0] = 1
        G_omega_d[1,6] = 1
        
        G_omega_n = np.zeros([2,12])
        
        B = np.zeros(self.N*2)
        Bux = np.zeros(self.N)
        Buy = np.zeros(self.N)
        Bvx = np.zeros(self.N)
        Bvy = np.zeros(self.N)
         
        
        for i in range(0,self.N):
                
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
                
            elif self.boundaryType[i] == 1: # Dirichlet Randbedingugen
            
                M_omega_d = np.append(M_omega,G_omega_d, axis = 0)
                W_omega_d = np.zeros([2*n+4,2*n+4])
                W_omega_d[:2*n+2,:2*n+2] = W_omega
            
                W_omega_d[-2,-2] = 1 # Wichtung bar_U
                W_omega_d[-1,-1] = 1 # Wichtung bar_v
                
                M = M_omega_d
                W = W_omega_d
                    
                c = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(M),W),M)),np.transpose(M)),W)
                
                B[i]        = c[0,-2] * self.barU[i] +  c[0,-1] * self.barV[i]
                B[i+self.N] = c[6,-2] * self.barU[i] +  c[6,-1] * self.barV[i]
                
                Bux[i]        = c[1,-2] * self.barU[i] +  c[1,-1] * self.barV[i]
                Buy[i]        = c[2,-2] * self.barU[i] +  c[2,-1] * self.barV[i]
                Bvx[i]        = c[7,-2] * self.barU[i] +  c[7,-1] * self.barV[i]
                Bvy[i]        = c[8,-2] * self.barU[i] +  c[8,-1] * self.barV[i]

            

                        
            elif self.boundaryType[i] == 2: # Neumann Randbedingugen
            
                nx = self.normals[i][0]
                ny = self.normals[i][1]
            
                G_omega_n[0,1] = nx   # u_x
                G_omega_n[0,2] = a0*ny              # u_y
                G_omega_n[0,7] = a0*ny                # v_x
                G_omega_n[0,8] = a2*nx                # v_y
                
                G_omega_n[1,1] = a2*ny                # u_x
                G_omega_n[1,2] = a0*nx                 # u_y
                G_omega_n[1,7] = a0*nx                # v_x
                G_omega_n[1,8] = ny   # v_y
                
                M_omega_n = np.append(M_omega,G_omega_n, axis = 0)
                W_omega_n = np.zeros([2*n+4,2*n+4])
                W_omega_n[:2*n+2,:2*n+2] = W_omega
            
                W_omega_n[-2,-2] = 1 # Wichtung p_x
                W_omega_n[-1,-1] = 1 # Wichtung p_y
                    
                M = M_omega_n
                W = W_omega_n
                    
                c = np.matmul(np.matmul(np.linalg.inv(np.matmul(np.matmul(np.transpose(M),W),M)),np.transpose(M)),W)
                    
                B[i]        = (c[0,-2] * self.Px[i] +  c[0,-1] * self.Py[i])*a3
                B[i+self.N] = (c[6,-2] * self.Px[i] +  c[6,-1] * self.Py[i])*a3
                
                Bux[i] = (c[1,-2] * self.Px[i] +  c[1,-1] * self.Py[i])*a3
                Buy[i] = (c[2,-2] * self.Px[i] +  c[2,-1] * self.Py[i])*a3
                Bvx[i] = (c[7,-2] * self.Px[i] +  c[7,-1] * self.Py[i])*a3
                Bvy[i] = (c[8,-2] * self.Px[i] +  c[8,-1] * self.Py[i])*a3
                    
                    
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
            
            dataUx += list(c[1,:n])
            dataUy += list(c[2,:n])
            dataVx += list(c[7,:n])
            dataVy += list(c[8,:n])
            row_ind_spat += [i]*n
            col_ind_spat += sup   
            
            dataUx += list(c[1,n:2*n])
            dataUy += list(c[2,n:2*n])
            dataVx += list(c[7,n:2*n])
            dataVy += list(c[8,n:2*n])
            row_ind_spat += [i]*n
            col_ind_spat += list(np.array(sup)+self.N)
                
        A = csc_matrix((data, (row_ind, col_ind)),shape=(2*self.N,2*self.N))
        
        Aux = csc_matrix((dataUx, (row_ind_spat, col_ind_spat)),shape=(self.N,2*self.N))
        Auy = csc_matrix((dataUy, (row_ind_spat, col_ind_spat)),shape=(self.N,2*self.N))
        
        Avx = csc_matrix((dataVx, (row_ind_spat, col_ind_spat)),shape=(self.N,2*self.N))
        Avy = csc_matrix((dataVy, (row_ind_spat, col_ind_spat)),shape=(self.N,2*self.N))
        
        
        sol = bicgstab(A,B)[0]
        
        
        self.coordsNew = np.zeros([self.N,2])
        self.coordsNew[:,0] = sol[:self.N]+self.Nodes[:,0]
        self.coordsNew[:,1] = sol[self.N:]+self.Nodes[:,1]
        
        
        uX = Aux.dot(sol) + Bux
        uY = Auy.dot(sol) + Buy
        vX = Avx.dot(sol) + Bvx
        vY = Avy.dot(sol) + Bvy
        
        self.epsilonXX = uX
        self.epsilonXY = 1/2*(uY+vX)
        self.epsilonYY = vY
        
        self.sigmaXX = self.epsilonXX * (lame1+2*lame2) + self.epsilonYY * lame1
        self.sigmaXY = self.epsilonXY * 2 * lame2
        self.sigmaYY = self.epsilonYY * (lame1+2*lame2) + self.epsilonXX * lame1
        
        self.sol = sol
        
            
                
                

                
                
                
                
                
        
                        

        
        

        
        