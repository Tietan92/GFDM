#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 18:40:24 2023

@author: stefan
"""


import os 
from gfdm2d import PreprocessorParent, PostprocessorParent
                
class Preprocessor(PreprocessorParent):
    
    # private attributes ------------------------------------------------------
    
    __specHeatCapacity    = [] # specific heat capacity of the material
    __thermalConductivity = [] # Thermal conductivity of the material
    __massDensity         = [] # Specific Mass of the material
    
    # privatemethods ----------------------------------------------------------  
      
    def __init__(self, projectName):
        # Constructor of the preprocessor class. Calls PreprocessorParent(projectName)
        
        super().__init__(projectName)
        self._PreprocessorParent__importSolversList(["linearHeatConduction2d_BiCGSTAB"])    
        
        
        
    def __checkIfComplete(self):
        # Controls if all the necessary simulation data was already set up.
        
        if self._PreprocessorParent__frequency == []:
            
            raise Exception('Missing Parameter: frequency')
        
        if self._PreprocessorParent__timeEnd == []:
            
            raise Exception('Missing Parameter: timeEnd')
        
        if self._Preprocessor__massDensity == []:
            
            raise Exception('Missing Parameter: massDensity')
        
        if self._Preprocessor__specHeatCapacity == []:
            
            raise Exception('Missing Parameter: specHeatCapacity')
        
        if self._Preprocessor__thermalConductivity == []:
            
            raise Exception('Missing Parameter: thermalConductivity')
        
        for i in range(0,len(self._PreprocessorParent__BoundaryConditions.values)):
            
            if self._PreprocessorParent__BoundaryConditions.values[i] == []:
                
                raise Exception('Boundary Condition is missing:' + str(self._PreprocessorParent__BoundaryConditions.physicalTag[i]))
    
    
    # public -----------------------------------------------------------------
    
    def setSpecHeatCapacity(self, specHeatCapacity):
        # Sets the specific heat capacity used for the simulation
        
        self.__specHeatCapacity = specHeatCapacity
        
        
    def setThermalConductivity(self, thermalConductivity):
        # Sets the thermal conductivity used in the simulation
        
        self.__thermalConductivity = thermalConductivity
        
        
    def setMassDensity(self, massDensity):
        # Sets the mass density used in the simulation
        
        self.__massDensity = massDensity
        
        
    def setBoundaryCondition(self, physicalTag, boundaryType, temp=0, filmCoeff=0, heatflux=0 ):
        # This function is used to implement boundary conditions for the linear heat conduction problem.
        # The function is making a list out of the input values and call the function
        #PreprocessorParent.setBoundConditionBase` with a the list as a input parameter
        
        self._PreprocessorParent__setBoundConditionBase(physicalTag, boundaryType, [temp, filmCoeff, heatflux], [1, 0, 1])
        
    
    def setInitialCondition(self,T0):
        # Sets the initial values for the temperature field of the simulation.
        # Calls the function `PreprocessorParent.setInitialConditionSkalar`.
        
        self._PreprocessorParent__setInitialConditionSkalar(T0)
        
                
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
        process.stdin.write((str(self.__massDensity) + "\n").encode('utf-8'))         # Mass Densitiy
        process.stdin.write((str(self.__specHeatCapacity) + "\n").encode('utf-8'))    # Specific heat capacity
        process.stdin.write((str(self.__thermalConductivity) + "\n").encode('utf-8')) # Thermal conductivity

        
        # Initial Conditions
        for i in range(0,self._PreprocessorParent__mesh._Mesh__Nodes.numNodes):
            
            process.stdin.write((str(self._PreprocessorParent__initialValues[i]) + "\n").encode('utf-8'))
            
        # Boundary Conditions
        for i in range(0,len(self._PreprocessorParent__BoundaryConditions.physicalTag)):
                       
            process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.physicalTag[i]) + "\n").encode('utf-8'))  # Physical Tag
            process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.boundaryType[i]) + "\n").encode('utf-8')) # Boundary Type
            
            for j in range(0,self._PreprocessorParent__numOfTimeSteps):
                
                process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.values[i][0][j]) + "\n").encode('utf-8')) # Temperature Values

            if self._PreprocessorParent__BoundaryConditions.boundaryType[i] == 2:
                
                process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.values[i][1]) + "\n").encode('utf-8')) # Film coefficient
            
                for j in range(0,self._PreprocessorParent__numOfTimeSteps):
                    
                    process.stdin.write((str(self._PreprocessorParent__BoundaryConditions.values[i][2][j]) + "\n").encode('utf-8')) # Heat Flux values
                    
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
     
    def triangulationPlotAnimated(self,fig,ax,name, colormap):
        # Calls PostprocessorParent.triangulationPlotAnimatedBase(fig,ax,PostprocessorParent.outputData,name)
        
        self._PostprocessorParent__triangulationPlotAnimatedBase(fig,ax,self._PostprocessorParent__outputData,name, colormap)
        
        
    def plotResults(self,fig,ax,tVal):
        # Calls PostprocessorParent.plotResultsBase(fig,ax,PostprocessorParent.outputData,tVal)
        
        self._PostprocessorParent__plotResultsBase(fig,ax,self._PostprocessorParent__outputData,tVal)
        
        
