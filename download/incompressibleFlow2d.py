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
    
    __massDensity  = [] # Specific Mass of the material
    __dynViscosity = [] # Dynamic Viscosity oft the material
    
    # privatemethods ----------------------------------------------------------  
      
    def __init__(self, projectName):
        # Constructor of the preprocessor class. Calls PreprocessorParent(projectName)
        
        super().__init__(projectName)
        self._PreprocessorParent__importSolversList(["incompressibleFlow2d_BiCGSTAB"])    
        
        
        
    def __checkIfComplete(self):
        # Controls if all the necessary simulation data was already set up.
        
        if self._PreprocessorParent__frequency == []:
            
            raise Exception('Missing Parameter: frequency')
        
        if self._PreprocessorParent__timeEnd == []:
            
            raise Exception('Missing Parameter: timeEnd')
        
        if self._Preprocessor__massDensity == []:
            
            raise Exception('Missing Parameter: massDensity')
        
        if self._Preprocessor__dynViscosity == []:
            
            raise Exception('Missing Parameter: dynViscosity')
        
        
        for i in range(0,len(self._PreprocessorParent__BoundaryConditions.values)):
            
            if self._PreprocessorParent__BoundaryConditions.values[i] == []:
                
                raise Exception('Boundary Condition is missing:' + str(self._PreprocessorParent__BoundaryConditions.physicalTag[i]))
    
    
    # public -----------------------------------------------------------------
            
    def setMassDensity(self, massDensity):
        # Sets the mass density used in the simulation
        
        self.__massDensity = massDensity
        
    def setDynViscosity(self, dynViscosity):
         # Sets the dynamic viscosity used in the simulation
            
        self.__dynViscosity = dynViscosity
        
        
    def setBoundaryCondition(self, physicalTag, boundaryType, values = [0,0]):
        # This function is used to implement boundary conditions for the imcompressible flow problem.
        # The function is making a list out of the input values and call the function 
        # `PreprocessorParent.setBoundConditionBase` with a the list as a input parameter
        
            
        if boundaryType == 'velocity':
                
            boundaryTypeID = 1
                
        elif boundaryType == 'outflow':
                
            boundaryTypeID = 2
                
        elif boundaryType == 'slip':
                    
            boundaryTypeID = 3
                
        else:
                
            raise Exception('The following boundary type does not exist:_' + boundaryType)
        

        self._PreprocessorParent__setBoundConditionBase(physicalTag, boundaryType=boundaryTypeID, listValues=values, listValuesDim=[1,1])
            
    
    def setInitialCondition(self,v0):
        # Sets the initial values for the velocity field of the simulation. 
        # Calls the function `PreprocessorParent.setInitialConditionSkalar`.
        
        self._PreprocessorParent__setInitialConditionBase(values=v0)
        
                
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
        process.stdin.write((str(self.__massDensity) + "\n").encode('utf-8'))  # Mass Densitiy
        process.stdin.write((str(self.__dynViscosity) + "\n").encode('utf-8')) # Dynamic Viscosity
        
        # Initial Conditions
        for i in range(0,self._PreprocessorParent__mesh._Mesh__Nodes.numNodes):
            
            process.stdin.write((str(self._PreprocessorParent__initialValues[0][i]) + "\n").encode('utf-8')) # initial velocity values in x-direction
            process.stdin.write((str(self._PreprocessorParent__initialValues[1][i]) + "\n").encode('utf-8')) # initial velocity values in y-direction
            
        # Boundary Conditions
        for i in range(0,len(self._PreprocessorParent__BoundaryConditions.physicalTag)):
                       
            process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.physicalTag[i]) + "\n").encode('utf-8'))     # Physical Tag
            process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.boundaryType[i]) + "\n").encode('utf-8')) # Boundary Type 
            
            for j in range(0,self._PreprocessorParent__numOfTimeSteps):
                
                process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.values[i][0][j]) + "\n").encode('utf-8')) # boundary values x-direction
                process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.values[i][1][j]) + "\n").encode('utf-8')) # boundary values in y-direction
                    
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
        
        if var == 'vel_x':
            
            varInt = 0
            
        elif var == 'vel_y':
            
            varInt = 1
            
        elif var == '|vel|':
            
            varInt = 2
            
        else:
            
            raise Exception('The variable name is not defined')
        
        self._PostprocessorParent__triangulationPlotAnimatedBase(fig,ax,varInt,name, colormap, valMin,valMax)
        
        
    def plotResults(self,fig,ax,var,tVal,colormap, valMin = [], valMax = []):
        # Calls PostprocessorParent.plotResultsBase(fig,ax,PostprocessorParent.outputData,tVal)
        
        if var == 'vel_x':
            
            varInt = 0
            
        elif var == 'vel_y':
            
            varInt = 1
            
        elif var == '|vel|':
            
            varInt = 2
            
        else:
            
            raise Exception('The variable name is not defined')
            
        
        self._PostprocessorParent__plotResultsBase(fig,ax,varInt,tVal,colormap,valMin,valMax)
        
        
