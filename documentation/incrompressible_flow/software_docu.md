## Mesh

| Module Informations |     |
| --- | --- |
| Associated Python File | `/source_code/python/gfdm2d.py` |
| Libraries needed | `numpy`,`matplotlib` |

### 1\. Description

The Mesh class is used to import a mesh in the msh format and save it in fitting Python containers. It should also offer basic functionalities to visualise the mesh and to see which boundaries exists and what are the labels for the boundaries.

### 2\. Public Methods

**Mesh(mshFile)**
Constructor of the mesh class. Reads an msh-File and stores all the necessary information inside the class attributes.

| Parameter |     |
| --- | --- |
| `mshFile` : str | Name of the .msh input file. Should have the structure "filename.msh" |

&nbsp;

**plotBoundary(ax, lw=2)**
Function which can be used to plot the boundary together with the associated labels

| Parameter |     |
| --- | --- |
| `ax` : matplotlib axes object | The new plot is added to this axes |
| `lw`: int or float | Optional: linewidths of the boundaries |

&nbsp;

**plotMesh(ax,lw=0.5)**
Function which can be used to plot a triangulation mesh.

| Parameter |     |
| --- | --- |
| `ax` : matplotlib axes object | The new plot is added to this axes |
| `lw`: int or float | Optional: linewidths of the connecting lines |

## Preprocessor

### 1\. Description
The Preprocessor class is used to build the model based on the imported mesh. The material parameters, initial conditions and boundary conditions as well as basic simulation parameters should be defined inside an object of this class. All the data is send to the solver via pipe. To improve the general usage, two classes are defined: PreprocessorParent is the parent class and provides general functions that can be used by several different solvers. The Preprocessor class already specifies functions that are required for a specific problem. The inheritance scheme is as follows:

### 2.\. Public Methods
**Preprocessor.setFrequency(frequency)** 
Set the simulation frequency.

| Parameter |     |
| --- | --- |
| `frequency` : int or float | Value for the simulation frequency |

&nbsp;

**Preprocessor.setTimeEnd(timeEnd)** 
Set the end time of the simulation.

| Parameter |     |
| --- | --- |
| `timeEnd` : int or float | Value for the end time of the simulation |

&nbsp;

**Preprocessor.setNumThreads(numThreads)**
Set the number of threads used for multiphreading.

| Parameter |     |
| --- | --- |
| `numThreads` : int | Value for the number of threads |

&nbsp;

**Preprocessor.setSolverTolerance(solverTolerance)**
Optional function. Sets the tolerance of the solver. If not called the tolerance is set to 1e-09.

| Parameter |     |
| --- | --- |
| `solverTolerance` : float | Value of the solver Tolerance |

&nbsp;

**Preprocessor.linkMesh(mesh)**
Links a mesh object to the Preprocessor

| Parameter |     |
| --- | --- |
| `mesh` : Mesh Object | The mesh that should be linked to the Preprocessor |

&nbsp;

**Preprocessor.setSolver(solverPath)**
This function sets the Path for the solver that should be used.

| Parameter |     |
| --- | --- |
| `solverPath` : str | Path of the solver. The solver name must be included in the path and must be a valid one which can be found in `PreprocessorParent.solvers` |

**Preprocessor.getNumNodes()**  
Returns an integer with the number of nodes

**Preprocessor.getNumTimeSteps()** 
Returns an integer with the number of time-steps
