# -*- coding: utf-8 -*-
"""
Created on Mon May  3 15:06:04 2021

@author: Stefan
"""
import numpy as np
from numpy import linalg as LA

class nonlinearHeatConduct:
    
    
    def __init__(self,pointcloud, f, tend, T0):
        self.pointcloud = pointcloud # Punktwolke
        self.f = f          # Abtastfrequenz in Hz
        self.tend = tend        # Simulationszeit in s
        self.t = np.arange(0,self.tend,1/self.f)
        
        
        # Materialparameter ###################################################
        
        self.Lambda = []         # Wärmeleitfähigkeit  
        self.c = []              # spezifische Wärmekapazität
        self.rho = []            # Dichte
        self.alpha = [0,0]          # Wärmeübergangskoeffizient
        self.K0 = -273.15           # Absoluter Nullpunkt
        self.T0 = T0 #-self.K0               # Anfangstemperatur

        
        self.dirichlets  = np.zeros([self.pointcloud.size,self.t.size])
        self.dirichlets[:,0] = self.T0
        self.neumanns    = np.zeros([self.pointcloud.size,self.t.size])

        
    def defineBoundaryCondition(self, physicalTag, boundaryType, Value):
        for i in np.arange(0,self.pointcloud.size):
            if self.pointcloud[i].physical_id == physicalTag:
                self.pointcloud[i].boundary_type = boundaryType
                if boundaryType == 1:
                    self.dirichlets[i] = Value #- self.K0
                elif boundaryType == 2:
                    self.neumanns[i] = Value
                elif boundaryType == 3:
                    self.neumanns[i] = Value #- self.K0
                    
                    
    def mergeBoundary(self, boundary_old, boundary_new):
        for i in np.arange(0,self.pointcloud.size):
            if np.all(np.isin(boundary_old,self.pointcloud[i].physical_id)) == True:
                self.pointcloud[i].physical_id = boundary_new
                
    def addBoundary(self, node_id, physical_id):
        for i in np.arange(0,self.pointcloud.size):
            if self.pointcloud[i].node_id == node_id:
                self.pointcloud[i].physical_id = physical_id
    
        
        
    
    def solve(self, solver ='pjg'):
        
        
        from scipy.sparse import csc_matrix
        from scipy.sparse.linalg import bicgstab

        N = self.pointcloud.size
        
        # Vor initialisierung #########################################################################################
        ###############################################################################################################
        
        for i in np.arange(0,N):
            
            support = self.pointcloud[i].support
            
            
            n = support.size
            
            xi = self.pointcloud[i].coord[0]
            yi = self.pointcloud[i].coord[1]
            
            D = np.zeros([n, 6])
            
            Wd = np.zeros([n,n])
            
            for j in np.arange(0,support.size):
                
              
                xj = self.pointcloud[support[j]].coord[0]
                yj = self.pointcloud[support[j]].coord[1]
                
                D[j,0] = 1
                D[j,1] = xj - xi 
                D[j,2] = yj - yi
                D[j,3] = 1/2*(xj - xi)**2 
                D[j,4] = (xj - xi )*(yj - yi )
                D[j,5] = 1/2*(yj - yi)**2
                
                h_i = self.pointcloud[i].h
                h_j = self.pointcloud[j].h
                
                Wd[j,j] = (np.exp(-5*LA.norm(np.array([xj, yj])-np.array([xi, yi]))**2/(h_i**2+h_j**2)))**2
                
                W_omega = np.zeros([n+1,n+1])
                W_omega[:n,:n] = Wd
                W_omega[-1,-1] = 1 # Wichtungsfaktor Innere Gebietspunkte
                
                W_omega_n = np.zeros([n+2,n+2])
                W_omega_n[:n+1,:n+1] = W_omega
                W_omega_n[-1,-1]  = 1 # Wichtungsfaktor Neumann Randbedingungen
                
            self.pointcloud[i].D = D
            self.pointcloud[i].W_omega = W_omega
            self.pointcloud[i].W_omega_n = W_omega_n
            
            

        Y=np.zeros([N,self.t.size]) # Temperatur
        Y[:,0] = self.dirichlets[:,0]
        I = np.arange(1,self.t.size,dtype='int')
        
        for k in I:
            
            #A = np.zeros([N,N])
            data = []
            row_ind = []
            col_ind = []

            beta  = np.zeros(N)
            gamma = np.zeros(N)
            zeta  = np.zeros(N)
            
            Lambda = np.interp( Y[:,k-1],#+self.K0 , 
                                   self.Lambda[0], 
                                   self.Lambda[1] , 
                                   left = self.Lambda[1,0], 
                                   right = self.Lambda[1,-1])
            
            heatcap = np.interp( Y[:,k-1],#+self.K0 , 
                                        self.c[0] , 
                                        self.c[1] , 
                                        left = self.c[1,0], 
                                        right = self.c[1,-1]) 
            
            LambdaT = np.interp( Y[:,k-1],#+self.K0 ,
                                 self.LambdaT[0] ,
                                 self.LambdaT[1] , 
                                 left = self.LambdaT[1,0], 
                                 right = self.LambdaT[1,-1])
            
            
            
            alpha = np.interp( Y[:,k-1],#+self.K0,
                               self.alpha[0],
                               self.alpha[1],
                               left = self.alpha[1,0],
                               right = self.alpha[1,-1])
        
            for i in np.arange(0,N):
                
                
                boundary_type = self.pointcloud[i].boundary_type
            
                support = self.pointcloud[i].support
            
                n = support.size
                
                # T,x und T,y ##########################################################################################
                ########################################################################################################
                
                Tx = 0
                Ty = 0
                
                if k > 1 and self.pointcloud[i].boundary_type != 1:
                    
                    c = self.pointcloud[i].c
                    
                    for j in np.arange(0,n):
                        
                        index_j = int(support[j])
                        
                        Tx += c[1,j] * Y[index_j,k-1]
                        Ty += c[2,j] * Y[index_j,k-1] 

                     
                           
                    Tx += c[1,n] * self.f  * Y[index_j,k-2]
                    Ty += c[2,n] * self.f  * Y[index_j,k-2]

            
                # D_i und Wd_i ############################################################################################
                ###########################################################################################################
           
                D =  self.pointcloud[i].D 

            
            # G_omega #################################################################################################
            ###########################################################################################################
                '''
                G_omega = np.array([[self.f,
                                 - LambdaT[i] * Tx  / (self.rho * heatcap[i]), 
                                 - LambdaT[i] * Ty  / (self.rho * heatcap[i]),
                                 - Lambda[i] / (self.rho * heatcap[i]),
                                 0,
                                 - Lambda[i] / (self.rho * heatcap[i])]]) 
                '''
                G_omega = np.zeros([1,6])
                G_omega[0,0]  = self.f 
                G_omega[0,1]  = - LambdaT[i] * Tx  / (self.rho * heatcap[i]) 
                G_omega[0,2]  = - LambdaT[i] * Ty  / (self.rho * heatcap[i]) 
                G_omega[0,3]  = - Lambda[i] / (self.rho * heatcap[i]) 
                G_omega[0,5]  = - Lambda[i] / (self.rho * heatcap[i]) 
                

            
            # innere Gebietspunkte ####################################################################################
            ###########################################################################################################
            
                M_omega = np.append(D,G_omega,axis=0)
                W_omega = self.pointcloud[i].W_omega
            
            
            # Zuweisung der Gleichungen zum Punkt #####################################################################
            ###########################################################################################################
        
                if boundary_type == 0: # Innere Gebietspunkte
            
                    M = M_omega
                    W = W_omega
                
                    c = np.matmul(np.matmul(LA.inv(np.matmul(np.matmul(np.transpose(M),W),M)),np.transpose(M)),W)
                    self.pointcloud[i].c = c
                    self.pointcloud[i].M = M
                    self.pointcloud[i].W = W
                
                    beta[i]  = c[0,n]
            
                elif boundary_type == 2: # Neumannn Randbedingungen
            
                    nx = self.pointcloud[i].normal[0]
                    ny = self.pointcloud[i].normal[1]
                    
                    #G_omega_n = np.array([[0,Lambda[i]*nx,Lambda[i]*ny,0,0,0]])

                    
                    G_omega_n = np.zeros([1,6])
                    G_omega_n[0,1] = Lambda[i]*nx
                    G_omega_n[0,2] = Lambda[i]*ny
                    
                
                    M_omega_n = np.append(M_omega,G_omega_n,axis=0)
                    W_omega_n = self.pointcloud[i].W_omega_n
                
                    M = M_omega_n
                    W = W_omega_n
                
                    c = np.matmul(np.matmul(LA.inv(np.matmul(np.matmul(np.transpose(M),W),M)),np.transpose(M)),W)
                    self.pointcloud[i].c = c
                    self.pointcloud[i].M = M
                    self.pointcloud[i].W = W
                
                    beta[i]  = c[0,n]
                    gamma[i] = c[0,n+1]
                    
                elif boundary_type == 3:   # Übergangsbedingungen
                    
                    nx = self.pointcloud[i].normal[0]
                    ny = self.pointcloud[i].normal[1]
                    
                    #G_omega_n = np.array([[alpha[i],Lambda[i]*nx,Lambda[i]*ny,0,0,0]])
                    
                    G_omega_n = np.zeros([1,6])
                    G_omega_n[0,0] = alpha[i]
                    G_omega_n[0,1] = Lambda[i]*nx
                    G_omega_n[0,2] = Lambda[i]*ny
                    
                    M_omega_n = np.append(M_omega,G_omega_n,axis=0)                  
                    W_omega_n = self.pointcloud[i].W_omega_n
                
                    M = M_omega_n
                    W = W_omega_n
                
                    c = np.matmul(np.matmul(LA.inv(np.matmul(np.matmul(np.transpose(M),W),M)),np.transpose(M)),W)
                    self.pointcloud[i].c = c
                    self.pointcloud[i].M = M
                    self.pointcloud[i].W = W
                
                    beta[i]  = c[0,-2]
                    gamma[i] = c[0,-1] * alpha[i] 
                    
                    
            # Einarbeitung in Matrix A
            
            
                if boundary_type == 1:
                
                    data.append(1) 
                    row_ind.append(i)
                    col_ind.append(i) 
                    zeta[i] = 1
                    
                else:
                    
                    data.append(1- c[0,0]) 
                    row_ind.append(i) 
                    col_ind.append(i)
                
                    for j in np.arange(1,n):
                    
                        index_j = int(support[j])
                        data.append(-c[0,j]) 
                        row_ind.append(i) 
                        col_ind.append(index_j)
                        
            
            A = csc_matrix((data, (row_ind, col_ind)),shape=(N,N))

            
            qa = self.neumanns[:,k]
            Td = self.dirichlets[:,k]
            #print(gamma[:10])
            #print(qa[:10])
            B = beta * self.f * Y[:,k-1] + gamma * qa  + zeta * Td 
            
            Y[:,k] = bicgstab(A, B, x0=np.copy(Y[:,k-1]))[0]   

            print(k)

        for i in np.arange(0,N):
            self.pointcloud[i].output = Y[i]#+self.K0
        return Y #+ self.K0