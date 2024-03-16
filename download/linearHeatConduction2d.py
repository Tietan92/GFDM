#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 18:40:24 2023

@author: stefan
"""


import os 
from gfdm2d import PreprocessorParent, PostprocessorParent
from variableNames import VariableNames

                
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
        
          
    def __writeInputFile(self):
        # Writes all the simulation data into an input file for the C++ solver
        
        varNames = VariableNames()
        fileName = self._PreprocessorParent__projectName + "_input.txt"
        
        if os.path.exists(fileName):
            
            os.remove(fileName)
        
        inputFile = open(fileName, "w")
        
        # Section Mesh
        inputFile.write(varNames.ImportMesh.strImportMesh + "\n")
        
        inputFile.write(varNames.ImportMesh.strLoadMesh + "\n")
        inputFile.write(self._PreprocessorParent__mesh._Mesh__mshFile + '\n')
        
        inputFile.write(varNames.ImportMesh.strImportMeshEnd + "\n")
        
        inputFile.write("\n")
        
        # Section simulation parameter
        inputFile.write(varNames.SimulationParameters.strSimulationParameters + "\n")
        
        inputFile.write(varNames.SimulationParameters.strFrequency + "\n")
        inputFile.write(str(self._PreprocessorParent__frequency) + "\n")
        
        inputFile.write(varNames.SimulationParameters.strTimeEnd + "\n")
        inputFile.write(str(self._PreprocessorParent__timeEnd) + "\n")
        
        inputFile.write(varNames.SimulationParameters.strNumNodes + "\n")
        inputFile.write(str(self._PreprocessorParent__mesh._Mesh__Nodes.numNodes) + "\n")
        
        inputFile.write(varNames.SimulationParameters.strNumTimesteps + "\n")
        inputFile.write(str(self._PreprocessorParent__numOfTimeSteps) + "\n")
        
        inputFile.write(varNames.SimulationParameters.strNumBoundaries + "\n")
        inputFile.write(str(self._PreprocessorParent__mesh._Mesh__PhysicalNames.numPhysicalNames-1) + "\n")
        
        if self._PreprocessorParent__numThreads != []:
            
            inputFile.write(varNames.SimulationParameters.strNumThreads + "\n")
            inputFile.write(str(self._PreprocessorParent__numThreads) + "\n")
            
        if self._PreprocessorParent__solverTolerance != []:
            
            inputFile.write(varNames.SimulationParameters.strSolverTolerance + "\n")
            inputFile.write(str(self._PreprocessorParent__solverTolerance) + "\n")  
            
        inputFile.write(varNames.SimulationParameters.strSimulationParametersEnd + "\n")
        
        inputFile.write("\n")
        
        # Section material Parameters
        inputFile.write(varNames.MaterialParameters.strMaterialParameters + "\n")
        
        inputFile.write(varNames.MaterialParameters.strMassDensity + "\n")
        inputFile.write(str(self.__massDensity) + "\n")
        
        inputFile.write(varNames.MaterialParameters.strSpecHeatCapacity + "\n")
        inputFile.write(str(self.__specHeatCapacity) + "\n")
        
        inputFile.write(varNames.MaterialParameters.strThermalConductivity + "\n")
        inputFile.write(str(self.__thermalConductivity) + "\n")
        
        inputFile.write(varNames.MaterialParameters.strMaterialParametersEnd + "\n")
        
        inputFile.write("\n")
        
        # Section initial conditions 
        inputFile.write(varNames.InitialConditions.strInitialConditions + "\n")
        symbolCount = 0
        
        for i in range(0,self._PreprocessorParent__mesh._Mesh__Nodes.numNodes):
            
            currentString = str(self._PreprocessorParent__initialValues[i])+ ', '
            inputFile.write(currentString)
            symbolCount += len(currentString)
            
            if symbolCount >= 200:
                inputFile.write("\n")
                symbolCount=0
                
        inputFile.write("\n")
        inputFile.write(varNames.InitialConditions.strInitialConditionsEnd + "\n")
        
        inputFile.write("\n")
        
        # Section Boundary Conditions
        inputFile.write(varNames.BoundaryConditions.strBoundaryConditions + "\n")
        
        for i in range(0,len(self._PreprocessorParent__BoundaryConditions.physicalTag)):
                       
            inputFile.write(varNames.BoundaryConditions.strDefineBoundaryCondition + "\n")
            
            inputFile.write(varNames.BoundaryConditions.strPhysicalTag + "\n")
            inputFile.write(str(self._PreprocessorParent__BoundaryConditions.physicalTag[i]) + "\n")
            
            inputFile.write(varNames.BoundaryConditions.strBoundaryType + "\n")
            inputFile.write(str(self._PreprocessorParent__BoundaryConditions.boundaryType[i]) + "\n")
            
            inputFile.write(varNames.BoundaryConditions.strTemperature + "\n")
            symbolCount = 0
            
            for j in range(0,self._PreprocessorParent__numOfTimeSteps):
                
                currentString = str(self._PreprocessorParent__BoundaryConditions.values[i][0][j]) + ', '
                inputFile.write(currentString)
                symbolCount += len(currentString)
                
                if symbolCount >= 200:
                    symbolCount = 0
                    inputFile.write("\n")
            
            if self._PreprocessorParent__BoundaryConditions.boundaryType[i] == 2:
                
                inputFile.write("\n")
                
                inputFile.write(varNames.BoundaryConditions.strFilmCoefficient + "\n")
                inputFile.write(str(self._PreprocessorParent__BoundaryConditions.values[i][1]) + "\n")
                
                inputFile.write(varNames.BoundaryConditions.strHeatflux + "\n")
                symbolCount = 0
                
                for j in range(0,self._PreprocessorParent__numOfTimeSteps):
                    
                    currentString = str(self._PreprocessorParent__BoundaryConditions.values[i][2][j]) + ', '
                    inputFile.write(currentString)
                    symbolCount += len(currentString)
                    
                    if symbolCount >= 200:
                        symbolCount = 0
                        inputFile.write("\n") 
                
            inputFile.write("\n")
            inputFile.write(varNames.BoundaryConditions.strDefineBoundaryConditionEnd + "\n")
                       
        inputFile.write(varNames.BoundaryConditions.strBoundaryConditionsEnd + "\n")
        inputFile.close()
        
        
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
        # Checks if all needed requirements are fulfilled, writes the simulation data into the input file and runs the solver. 
        # Gives an live output of the simulation progress.
        
        import subprocess
        import time
        
        print("Solver Running...")
        self.__checkIfComplete()
        self.__writeInputFile()
        
        solverPathSplitted = self._PreprocessorParent__solver.split("/")
        solverName = solverPathSplitted[-1]
        solverFolder = ""
        
        for i in range(0,len(solverPathSplitted)-1):
            solverFolder += solverPathSplitted[i] + "/"
        
        args = "cd " + solverFolder + " ; " + " stdbuf -oL -eL " + "./" + solverName + " " +  self._PreprocessorParent__projectName + " " + "\u0022" + os.getcwd() + "\u0022"
        start_time = time.time()
        
        process = subprocess.Popen(args, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

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
        
        
