#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 07:24:40 2023

@author: stefan
"""
import numpy as np

class Mesh():
    
    # private attributes-------------------------------------------------------
    
    class __PhysicalNames():
        
        numPhysicalNames = [] # Number of physical groups
        dimension        = [] # Dimension of the topology of the physical group:
        physicalTag      = [] # Physical Tag is the identification number of the physical group
        name             = [] # Identification string of the physical group
        
    class __Nodes():
        
        numNodes = [] # Number of nodes
        coord    = [] # Array with the coordinates of the nodes in 3 dimensions
        
    class __InnerElements():
        
        numInnerElements = [] # Number of inner elements (triangles)
        elementTag       = [] # Number of inner elements (triangles)
        nodeTag          = [] # The index of the nodes for each triangle
        
    class __BoundaryElements():
        
        numBoundaryElements = [] # Number of boundary elements (lines)
        elementTag          = [] # The physical tag of each line element
        nodeTag             = [] # The index of the nodes included in the line element
        
    __mshFile = [] # Name of the msh file
    
    # private methods ---------------------------------------------------------
        
    def __init__(self, mshFile):
        # Konstruktor of the mesh class. Reads an msh-Files and stores all the nesessary informations inside the class attributes.
        
        self.__mshFile = mshFile
        mesh = open(mshFile,'r')  
        loop = True 
     
        while loop == True:    
            line = mesh.readline()
        
            if line[:-1] == '$PhysicalNames':

                line = mesh.readline()
                
                self.__PhysicalNames.numPhysicalNames = int(line[:-1])
                self.__PhysicalNames.dimension        = np.zeros(self.__PhysicalNames.numPhysicalNames, dtype='int')
                self.__PhysicalNames.physicalTag      = np.zeros(self.__PhysicalNames.numPhysicalNames, dtype = 'int')
                self.__PhysicalNames.name             = np.zeros(self.__PhysicalNames.numPhysicalNames, dtype='<U127')
                
                for i in range(0,self.__PhysicalNames.numPhysicalNames):
                    
                    line = mesh.readline()
                    values = line[:-1].split(' ')
                    
                    self.__PhysicalNames.dimension[i]   = int(values[0])
                    self.__PhysicalNames.physicalTag[i] = int(values[1])      
                    self.__PhysicalNames.name[i]        = values[2][1:-1]
                
            elif line[:-1] == '$Nodes':
                
                line = mesh.readline()
                
                self.__Nodes.numNodes = int(line[:-1])
                self.__Nodes.coord    = np.zeros([self.__Nodes.numNodes,3])
                
                for i in range(0,self.__Nodes.numNodes):
                    
                    line = mesh.readline()
                    values = line[:-1].split(' ')
                    
                    self.__Nodes.coord[i,:] = [float(values[1]), float(values[2]), float(values[3])]
            
            elif line[:-1] == '$Elements':
                
                    line = mesh.readline()
                    numElements = int(line[:-1])
                    
                    boundaryElementsNum        = 0
                    boundaryElementsElementTag = []
                    boundaryElementsNodeTag    = []
                        
                    innerElementsNum        = 0
                    innerElementsElementTag = []
                    innerElementsNodeTag    = []
                    
                    for i in range(0,numElements): 
                        
                        line = mesh.readline()
                        values = line[:-1].split(' ')
                        
                        nodeTags = []
                        
                        for j in range(0,len(values)-5):
                            
                            nodeTags.append(int(values[5+j])-1)
                            
                        if int(values[1]) == 1:
                        
                            boundaryElementsNum += 1
                            boundaryElementsElementTag.append(int(values[3]))
                            boundaryElementsNodeTag.append(nodeTags)
                            
                        else:
                            
                            innerElementsNum += 1
                            innerElementsElementTag.append(int(values[3]))
                            innerElementsNodeTag.append(nodeTags)                   
                        
            elif line[:-1] == '$EndElements':
                
                    loop = False
        
        self.__BoundaryElements.numBoundaryElements = np.array(boundaryElementsNum)
        self.__BoundaryElements.elementTag          = np.array(boundaryElementsElementTag)
        self.__BoundaryElements.nodeTag             = np.array(boundaryElementsNodeTag)
        
        self.__InnerElements.numInnerElements = np.array(innerElementsNum)
        self.__InnerElements.elementTag       = np.array(innerElementsElementTag)
        self.__InnerElements.nodeTag          = np.array(innerElementsNodeTag)
        
        
    # public methods ----------------------------------------------------------
        
    def plotBoundary(self, ax, lw=2):
        # Function which can be used to plot the boundary together with the associated labels
        
        import matplotlib as mpl

        coordX = self.__Nodes.coord[:,0]
        coordY = self.__Nodes.coord[:,1]
        physicalIds = []
        physicalNames = []
        boundCoordX = []
        boundCoordY = []
        
        for i in range(0,self.__PhysicalNames.numPhysicalNames):
            
            if self.__PhysicalNames.dimension[i] == 1:
                
                physicalIds.append(self.__PhysicalNames.physicalTag[i])
                physicalNames.append(self.__PhysicalNames.name[i])
                boundCoordX.append([])
                boundCoordX.append([])
                boundCoordY.append([])
                boundCoordY.append([])
        
        colorRange = np.linspace(0,1,len(physicalIds))

        colors = mpl.colormaps['tab20b']
        
        for i in range(0,self.__BoundaryElements.numBoundaryElements):
            
            index = physicalIds.index(self.__BoundaryElements.elementTag[i])
            
            boundCoordX[2*index].append(coordX[self.__BoundaryElements.nodeTag[i,0]])
            boundCoordX[2*index+1].append(coordX[self.__BoundaryElements.nodeTag[i,1]])
            boundCoordY[2*index].append(coordY[self.__BoundaryElements.nodeTag[i,0]])
            boundCoordY[2*index+1].append(coordY[self.__BoundaryElements.nodeTag[i,1]])
            
        for i in range(0,len(physicalIds)):
            
            label = physicalNames[i] + ' (' + str(physicalIds[i]) + ')'
            color = colors(colorRange[i])
            
            segs = np.zeros((len(boundCoordX[2*i]),2, 2))
            segs[:, :, 1] = np.array([boundCoordY[2*i],boundCoordY[2*i+1]]).transpose()
            segs[:, :, 0] = np.array([boundCoordX[2*i],boundCoordX[2*i+1]]).transpose()
            
            line_segments = mpl.collections.LineCollection(segs,linewidths=lw, colors = color, label=label)
            ax.add_collection(line_segments)
            
        ax.axis('equal')
        
        
    def plotMesh(self,ax,lw=0.5):
        # Function which can be used to plot a triangulation mesh.
            
        import matplotlib as mpl
            
        coordX = self.__Nodes.coord[:,0]
        coordY = self.__Nodes.coord[:,1]
        
        triangles = []
        
        for i in range(0,self.__InnerElements.numInnerElements):
                
            triangles.append(self.__InnerElements.nodeTag[i])
                
        triang = mpl.tri.Triangulation(coordX, coordY, triangles)
        ax.triplot(triang, 'k-', lw=lw)
        ax.axis('equal')
        
        
class PreprocessorParent():
    
    # private attributes ------------------------------------------------------
    
    __projectName     = [] # Name which will be used for the simulation project
    __mesh            = [] # Mesh class object which stores the information for the mesh
    __solvers         = [] # A list with the names of the available solvers
    __numOfTimeSteps  = [] # The number of time-steps needed for the specified frequency and final time.
    __frequency       = [] # Frequency of the simulation. Describes how many time steps are processed during one time unit
    __timeEnd         = [] # The final time when the simulation stops
    __numThreads      = [] # The number of threads used for multiphreading implementation
    __solverTolerance = [] # Tolerance of the solver
    __initialValues   = [] # List with values which satisfy the initial conditions.
            
    class __BoundaryConditions():
        
        physicalTag  = [] # List with the defined physical tags of the boundaries which are saved in PreprocessorParent.mesh
        name         = [] # List of the defined names of the boundaries which are saved in PreprocessorParent.mesh
        boundaryType = [] # List which stores the boundary type as an integer for each boundary. 0 = "Dirichlet", 1 = "Neumann"
        values       = [] # List which stores the boundary condition data for each defined boundary
        possibleConditions = ['Dirichlet','Neumann'] # List with the possible boundary conditions.
        
    # private Methods -----------------------------------------------------------------    
    
    def __init__(self, projectName):
        # Constructor of the preprocessor parent class. The only attribute which is needed for the initialisation is the 
        # name of the project.
        
        self.__projectName = projectName

    def __importSolversList(self,solverList):
        # Imports a list with the names of all possible solvers for the specific problem
        
        self.__solvers = solverList
        
        
    def __setBoundConditionBase(self, physicalTag, boundaryType ,listValues, listValuesDim):
        # The base function which is used to set boundary conditions.
        
        listValuesApproved = []
            
        if self.__numOfTimeSteps == []:
            
            raise Exception('boundary condition error: please define "timeEnd" and "frequency" first')
        
        if self.__BoundaryConditions.physicalTag == []:
            
            raise Exception('boundary condition error: please import a mesh first')
            
        for i in range(0,len(listValuesDim)):
            
            # listValuesDim[i] = 0 --> listValues[i] must be skalar
            # listValuesDim[i] = 1 --> listValues[i] must be an 1-dimensional array
            # listValuesDim[i] = 2 --> listValues[i] must be an 2-dimensional array 
            
            if listValuesDim[i] == 0:
                
                if type(listValues[i]) != int and type(listValues[i]) == float:
                    
                    raise Exception('Dimension of input is wrong: skalar value expected')
                    
                listValuesApproved.append(listValues[i])
                    
            elif listValuesDim[i] == 1:
                
                if type(listValues[i]) == int or type(listValues[i]) == float:
                    
                    listValuesApproved.append([listValues[i]]*self.__numOfTimeSteps)
                    
                else:
                    
                    if self.__numOfTimeSteps == len(listValues[i]):
                        
                        listValuesApproved.append(listValues[i])
                        
                    else:
                        
                        raise Exception('boundary condition error: the length of the boundary values must be equal to "numOfTimeSteps"')
                
            elif listValuesDim[i] == 2:
                
                pass
                    
        if type(physicalTag) == str:
                        
            if physicalTag not in self.__BoundaryConditions.name:
                
                raise Exception('boundary condition error: the physical Tag does not exist')
                            
            index = self.__BoundaryConditions.name.index(physicalTag)            
                        
        elif type(physicalTag) == int:
            
            if physicalTag not in self.__BoundaryConditions.physicalTag:
                
                raise Exception('boundary condition error: the physical Tag does not exist')
                        
            index = self.__BoundaryConditions.physicalTag.index(physicalTag)
                    
        self.__BoundaryConditions.values[index] = listValuesApproved
                
        if boundaryType not in self.__BoundaryConditions.possibleConditions:
            
            raise Exception('boundary condition error: this boundary type is not defined')
                            
        self.__BoundaryConditions.boundaryType[index] = self.__BoundaryConditions.possibleConditions.index(boundaryType)+1
        
                    
    def __setInitialConditionSkalar(self,values):
        # A function to implement initial conditions for skalar fields (e.g. Temperature)
        
        if self.__mesh  == []:
            
            raise Exception('initial condition error: please import a mesh first')
        
        if type(values) == int or type(values) == float:
                
            values = [values] * self.__mesh._Mesh__Nodes.numNodes
                
        if self.__mesh._Mesh__Nodes.numNodes != len(values):
            
            raise Exception('initial condition error: the number of values must be equal to the number of Nodes')
                
        self.__initialValues = values    
        
    
    # public methods ----------------------------------------------------------
        
    def setFrequency(self, frequency):
        # Set the simulation frequency.
        
        self.__frequency = frequency
        
        if self.__timeEnd != []:
            
            self.__numOfTimeSteps = int(self.__frequency*self.__timeEnd)+1
            
        
    def setTimeEnd(self, timeEnd):
        # Set the end time of the simulation.
        
        self.__timeEnd = timeEnd
        
        if self.__frequency != []:
            
            self.__numOfTimeSteps = int(self.__frequency*self.__timeEnd)+1
            
            
    def setNumThreads(self, numThreads):
        # Set the number of threads used for multiphreading.
        
        self.__numThreads = numThreads
        
        
    def setSolverTolerance(self, solverTolerance):
        # Optional function. Sets the tolerance of the solver. If not called the tolerance is set to 1e-09.
        
        self.__solverTolerance = solverTolerance
        
            
    def linkMesh(self,mesh):
        # Links a mesh object to the Preprocessor
        
        self.__mesh = mesh
        
        for i in range(0,self.__mesh._Mesh__PhysicalNames.numPhysicalNames):
            
            if self.__mesh._Mesh__PhysicalNames.dimension[i] == 1:
                
                self.__BoundaryConditions.physicalTag.append(self.__mesh._Mesh__PhysicalNames.physicalTag[i])
                self.__BoundaryConditions.name.append(self.__mesh._Mesh__PhysicalNames.name[i])
                self.__BoundaryConditions.values.append([])
                self.__BoundaryConditions.boundaryType.append([])
                
                
    def setSolver(self,solverPath):
        # This function sets the Path for the solver that should be used.
        
        solverName = solverPath.split("/")[-1]
        
        if solverName in self.__solvers:
            
            self.__solver = solverPath
            
        else:
            
            raise Exception('Solver is not in Database')

        
class PostprocessorParent():
    
    # private attributes ------------------------------------------------------
    
    __projectName = [] # The name of the existing project the Postprocessor should analyse.
    __outputData  = [] # The output data imported as an array from the output file
    __tValues     = [] # An array which saves the time values for each time-step performed in the simulation
    __mesh        = [] # The Mesh class object which corresponds to the output data
    
    # private methods ---------------------------------------------------------
             
    def __init__(self, projectName):
        # Konstructor of the PostprocessorParent class. Reads the simulation output file specified 
        # in the parameter and saves the data in the corresponding arrays.
        
        self.__projectName = projectName
        
        outputData = []
        tValues = []
        fileName = self.__projectName + "_output.txt"
        outputFile = open(fileName,'r')
        loop = True 
             
        while loop == True:    
            line = outputFile.readline()
                    
            if line[:5] == "<t = ":
                tValues.append(float(line[5:-2]))
                loop2 = True
                valuesTemp = []
                line = outputFile.readline()
                while loop2 == True:
                                
                    values = line[:-1].split(', ')
                    for i in range(0,len(values)-1):
                        valuesTemp.append(float(values[i]))
                                
                    line = outputFile.readline()
                                     
                    if line[:-1] == "<timestepEnd>":
                        loop2 = False
                        outputData.append(valuesTemp)
                                
            if line[:-1] == "<$SimulationOutputEnd>":
                loop = False
                        
        self.__outputData = np.array(outputData).transpose()
        self.__tValues = tValues
        
                        
    def __triangulationPlotAnimatedBase(self,fig,ax,data,name, colormap):
        # Generates an animated triangulation plot of skalar time dependent data and saves it as an .mp4 file.
            
        import matplotlib.tri as mtri
        import matplotlib.animation as animation
        from matplotlib import pyplot as plt
        from IPython.display import Video, display
            
        print("Generate Animation...")
        triangles = []
        for i in range(0,self.__mesh._Mesh__InnerElements.numInnerElements):
                
            triangles.append(self.__mesh._Mesh__InnerElements.nodeTag[i])
                
        triang = mtri.Triangulation(self.__mesh._Mesh__Nodes.coord[:,0], self.__mesh._Mesh__Nodes.coord[:,1], triangles)
            
        valMax = np.max(np.max(data))
        valMin = np.min(np.min(data))
        
        frequency = 1/(self.__tValues[1]-self.__tValues[0])
            
        ax.axis('equal')
        
        def update_plot(frame_number):
                
            ax.clear()
            ax.tripcolor(triang, data[:,frame_number], vmin=valMin, vmax=valMax, cmap=colormap, shading = 'gouraud')
            
        plot = ax.tripcolor(triang, data[:,0], vmin=valMin, vmax=valMax, cmap=colormap, shading = 'gouraud')
        fig.colorbar(plot)
        plt.close()
        
        ani = animation.FuncAnimation(fig, update_plot, interval=1/frequency*1000,  save_count=int(frequency*self.__tValues[-1]))
            
        writer = animation.FFMpegWriter(fps=frequency, metadata=dict(artist='Me'), bitrate=-1)
        ani.save(name + ".mp4", writer=writer)
        
        display(Video(name + ".mp4"))
        

    def __plotResultsBase(self,fig,ax,data,tVal):
        # Generates different sub-plots for fixed time points
            
        ax = np.array(ax)
            
        import matplotlib.tri as mtri
        triangles = []
        for i in range(0,self.__mesh._Mesh__InnerElements.numInnerElements):
                
            triangles.append(self.__mesh._Mesh__InnerElements.nodeTag[i])
                
        triang = mtri.Triangulation(self.__mesh._Mesh__Nodes.coord[:,0], self.__mesh._Mesh__Nodes.coord[:,1], triangles)
            
        valMax = np.max(np.max(data))
        valMin = np.min(np.min(data))
        ax = ax.flatten()
        frequency = 1/(self.__tValues[1]-self.__tValues[0])
            
        for i in range(0,len(tVal)):
            k = round(tVal[i]*frequency)
            ax[i].axis('equal')
            plot = ax[i].tripcolor(triang, data[:,k], vmin=valMin, vmax=valMax, cmap='rainbow', shading = 'gouraud')
            fig.colorbar(plot)
            ax[i].set_title("t=" + str(tVal[i]))
            
            
    # public methods ----------------------------------------------------------
            
    def linkMesh(self,mesh):
        # Links a mesh object to the Postprocessor. 
        # The number of nodes from the simulation output file must be equal to the number of nodes of the imported mesh.
                
        self.__mesh = mesh
        
