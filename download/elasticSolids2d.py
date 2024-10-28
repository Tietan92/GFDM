#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 18:40:24 2023

@author: stefan
"""


import os 
import numpy as np
from gfdm2d import PreprocessorParent, PostprocessorParent

                
class Preprocessor(PreprocessorParent):
    
    # private attributes ------------------------------------------------------
    
    __massDensity    = [] # Specific Mass of the material
    __elasticModulus = [] # elastic Modulus of the material
    __poissonRatio   = [] # Poisson Ratio of the material
    
    # privatemethods ----------------------------------------------------------  
      
    def __init__(self, projectName):
        # Constructor of the preprocessor class. Calls PreprocessorParent(projectName)
        
        super().__init__(projectName)
        self._PreprocessorParent__importSolversList(["elasticSolids2d"])    
        
        
        
    def __checkIfComplete(self):
        # Controls if all the necessary simulation data was already set up.
        
        if self._PreprocessorParent__frequency == []:
            
            raise Exception('Missing Parameter: frequency')
        
        if self._PreprocessorParent__timeEnd == []:
            
            raise Exception('Missing Parameter: timeEnd')
        
        if self._Preprocessor__massDensity == []:
            
            raise Exception('Missing Parameter: massDensity')
        
        if self._Preprocessor__elasticModulus == []:
            
            raise Exception('Missing Parameter: elastic Modulus')
            
        if self._Preprocessor__poissonRatio == []:
                
            raise Exception('Missing Parameter: Poisson Ratio')
        
        for i in range(0,len(self._PreprocessorParent__BoundaryConditions.values)):
            
            if self._PreprocessorParent__BoundaryConditions.values[i] == []:
                
                raise Exception('Boundary Condition is missing:' + str(self._PreprocessorParent__BoundaryConditions.physicalTag[i]))
    
    
    # public -----------------------------------------------------------------
            
    def setMassDensity(self, massDensity):
        # Sets the mass density used in the simulation
        
        self.__massDensity = massDensity
        
    def setElasticModulus(self, elasticModulus):
         # Sets the  Young Modulus used in the simulation
            
        self.__elasticModulus = elasticModulus
        
    def setPoissonRatio(self, poissonRatio):
        # Sets the  Young Modulus used in the simulation
                
        self.__poissonRatio = poissonRatio
        
        
    def setBoundaryCondition(self, physicalTag, boundaryType, values):
        # This function is used to implement boundary conditions for elastic solids.
        # The function is making a list out of the input values and call the function 
        # `PreprocessorParent.setBoundConditionBase` with a the list as a input parameter
        
            
        if boundaryType == 'velocity' or boundaryType == 'displacement':
                
            boundaryTypeID = 1
                
        elif boundaryType == 'stress':
                
            boundaryTypeID = 2
                
        else:
                
            raise Exception('The following boundary type does not exist:_' + boundaryType)
        

        self._PreprocessorParent__setBoundConditionBase(physicalTag, boundaryType=boundaryTypeID, listValues=values, listValuesDim=[1,1])
            
    
    def setInitialCondition(self,vel,stress):
        # Sets the initial values for the velocity field of the simulation. 
        # Calls the function `PreprocessorParent.setInitialConditionSkalar`.
        
        self._PreprocessorParent__setInitialConditionBase(values=[vel[0],vel[1],stress[0],stress[1],stress[2]])
        
                
    def runSolver(self):
        # Checks if all needed requirements are fulfilled, send the simulation data via pipe and runs the solver. 
        # Gives an live output of the simulation progress.
        
        import subprocess
        import time
        
        print("Solver Running...")
        self.__checkIfComplete()
        
        solverPathSplitted = self._PreprocessorParent__solver.split("/")
        solverName = solverPathSplitted[-1]
        solverFolder = ""
        
        for i in range(0,len(solverPathSplitted)-1):
            solverFolder += solverPathSplitted[i] + "/"
        
        args = "cd " + solverFolder + " ; " + " stdbuf -oL -eL " + "./" + solverName + " " +  self._PreprocessorParent__projectName + " " + "\u0022" + os.getcwd() + "\u0022"
        
        start_time = time.time()
        process = subprocess.Popen(args, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
        
        # Send simulation parameter
        process.stdin.write((self._PreprocessorParent__mesh._Mesh__mshFile + '\n').encode('utf-8'))                               # Name of msh-File
        process.stdin.write((str(self._PreprocessorParent__frequency) + "\n").encode('utf-8'))                                    # Frequency
        process.stdin.write((str(self._PreprocessorParent__timeEnd) + "\n").encode('utf-8'))                                      # End time of the simulation
        process.stdin.write((str(self._PreprocessorParent__mesh._Mesh__Nodes.numNodes) + "\n").encode('utf-8'))                   # Number of Nodes
        process.stdin.write((str(self._PreprocessorParent__numOfTimeSteps) + "\n").encode('utf-8'))                               # Number if time-steps
        process.stdin.write((str(self._PreprocessorParent__mesh._Mesh__PhysicalNames.numPhysicalNames-1) + "\n").encode('utf-8')) # Number of boundaries
        process.stdin.write((str(self._PreprocessorParent__numThreads) + "\n").encode('utf-8'))                                   # Number of threads
        process.stdin.write((str(self._PreprocessorParent__solverTolerance) + "\n").encode('utf-8'))                              # Solver tolerance

        
        # Send material parameter
        process.stdin.write((str(self.__massDensity) + "\n").encode('utf-8'))     # Mass Densitiy
        process.stdin.write((str(self.__elasticModulus) + "\n").encode('utf-8'))  # elastic Modulus
        process.stdin.write((str(self.__poissonRatio) + "\n").encode('utf-8'))    # Poisson Ratio
    
        # Initial Conditions
        for i in range(0,self._PreprocessorParent__mesh._Mesh__Nodes.numNodes):
            
            process.stdin.write((str(self._PreprocessorParent__initialValues[0][i]) + "\n").encode('utf-8')) # initial velocity values in x-direction
            process.stdin.write((str(self._PreprocessorParent__initialValues[1][i]) + "\n").encode('utf-8')) # initial velocity values in y-direction
            process.stdin.write((str(self._PreprocessorParent__initialValues[2][i]) + "\n").encode('utf-8')) # stress sigma_xx
            process.stdin.write((str(self._PreprocessorParent__initialValues[3][i]) + "\n").encode('utf-8')) # stress sigma_xy
            process.stdin.write((str(self._PreprocessorParent__initialValues[4][i]) + "\n").encode('utf-8')) # stress sigma_yy
            
        # Boundary Conditions
        for i in range(0,len(self._PreprocessorParent__BoundaryConditions.physicalTag)):
                       
            process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.physicalTag[i]) + "\n").encode('utf-8'))  # Physical Tag
            process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.boundaryType[i]) + "\n").encode('utf-8')) # Boundary Type 
            
            for j in range(0,self._PreprocessorParent__numOfTimeSteps):
                
                process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.values[i][0][j]) + "\n").encode('utf-8')) # boundary values x-direction
                process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.values[i][1][j]) + "\n").encode('utf-8')) # boundary values y-direction
                    
        process.stdin.flush() # move the buffered data to console


        while process.poll() is None:
            
            nextline = process.stdout.readline()
            
            if nextline == '':
                
                continue
            
            print(nextline.strip().decode('ascii'))
            
            if nextline.strip().decode('ascii')=="100%":
                
                break

        print("Time needed: " + str(time.time() - start_time) + " seconds")
        
        
class Postprocessor(PostprocessorParent):
    
    # private methods ---------------------------------------------------------
    
    def __init__(self, projectName):
        # Constructor of the Postprocessor class. Calls `PostprocessorParent(projectName)`
        
        super().__init__(projectName)
        
    # public methods ----------------------------------------------------------
     
    def triangulationPlotAnimated(self,fig,ax, var ,name, colormap,valMin = [], valMax = []):
        # Calls PostprocessorParent.triangulationPlotAnimatedBase(fig,ax,PostprocessorParent.outputData,name)
        
        if var == 'sigma_xx':
            
            varInt = 0
            
        elif var == 'sigma_xy':
            
            varInt = 1
            
        elif var == 'sigma_yy':
            
            varInt = 2
        
        elif var == 'mises':
                
            varInt = 3
            
        else:
            
            raise Exception('The variable name is not defined')
        
        self._PostprocessorParent__triangulationPlotAnimatedUpdatedLagrangeBase(fig,ax,varInt,name, colormap, valMin,valMax)
        
        
    def plotResults(self,fig,ax,var,tVal,colormap, valMin = [], valMax = []):
        # Calls PostprocessorParent.plotResultsBase(fig,ax,PostprocessorParent.outputData,tVal)
        
        if var == 'sigma_xx':
            
            varInt = 0
            
        elif var == 'sigma_xy':
            
            varInt = 1
            
        elif var == 'sigma_yy':
            
            varInt = 2
            
        elif var == 'mises':
                
            varInt = 3
            
        else:
            
            raise Exception('The variable name is not defined')
            
        
        self._PostprocessorParent__plotResultsUpdatedLagrangeBase(fig,ax,varInt,tVal,colormap,valMin,valMax)
        
        
