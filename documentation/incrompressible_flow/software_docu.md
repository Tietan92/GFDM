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
The Preprocessor class is used to build the model based on the imported mesh. The material parameters, initial conditions and boundary conditions as well as basic simulation parameters should be defined inside an object of this class. All the data is send to the solver via pipe. To improve the general usage, two classes are defined: PreprocessorParent is the parent class and provides general functions that can be used by several different solvers. 

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

**Preprocessor.setMassDensity(massDensity)**
Sets the mass density used in the simulation

| Parameter |     |
| --- | --- |
| `massDensity` : int or float | The value for mass density of the material |

&nbsp;

**Preprocessor.setDynViscosity(dynViscosity**)
Sets the dynamic viscosity used in the simulation

| Parameter |     |
| --- | --- |
| `dynViscosity` : int or float | The value for the dynamic viscosity of the material |

&nbsp;

**Preprocessor.setBoundaryCondition(physicalTag,  boundaryType, values)**
This function is used to implement boundary conditions for the incompressible flow problem. The function is making a list out of the input values and call the function `PreprocessorParent.setBoundConditionBase` with a the list as a input parameter

<table class="jop-noMdConv"><thead class="jop-noMdConv"><tr class="jop-noMdConv"><th class="jop-noMdConv">Parameter</th><th class="jop-noMdConv"></th></tr></thead><tbody class="jop-noMdConv"><tr class="jop-noMdConv"><td class="jop-noMdConv"><code class="inline-code jop-noMdConv" spellcheck="false">physicalTag</code> : int or str</td><td class="jop-noMdConv">The corresponding physical tag of the boundary the conditions should be specified. If the input is from type <code spellcheck="false" class="jop-noMdConv">int</code> it must be included in&nbsp; <code spellcheck="false" class="jop-noMdConv">PreprocessorParent.BoundaryConditions.physicalTag</code>. If the input is from type <code spellcheck="false" class="jop-noMdConv">str</code> it must be included in <code spellcheck="false" class="jop-noMdConv">PreprocessorParent.BoundaryConditions.name</code></td></tr><tr class="jop-noMdConv"><td class="jop-noMdConv"><code spellcheck="false" class="jop-noMdConv">boundaryType</code> : List of str</td><td class="jop-noMdConv"><p>The input must be a list with two elements, each of which defines a boundary condition for a direction. The following boundary conditions are possible:</p><ul class="jop-noMdConv"><li class="jop-noMdConv"><code spellcheck="false" class="jop-noMdConv">'velocity'</code> : the velocity values will be given on the boundary for the defined direction</li><li class="jop-noMdConv"><code spellcheck="false" class="jop-noMdConv">'outflow'</code> : homogeneous Neumann boundary conditions used to map the outflow of the domain . If this boundary condition is applied, no parameter for values has to be specified.&nbsp;</li><li class="jop-noMdConv"><code spellcheck="false">'slip'</code> : Slip boundary condition: The velocity component in the normal direction of the boundary is 0. The velocity component in tangential direction is not constraints. If this boundary condition is applied, no parameter for <code spellcheck="false">values</code> has to be specified.</li></ul></td></tr><tr class="jop-noMdConv"><td class="jop-noMdConv"><code spellcheck="false" class="jop-noMdConv">values</code>: List of float or array</td><td class="jop-noMdConv"><p>Only needed in case velocity boundary conditions are applied. The input must be a list with two elements, containing the velocity values in x- and y-direction. If a single value is specified as an element of the list, it is assumed that the boundary condition is constant over time. If an array is specified, it is assumed that the array corresponds to the boundary values for each time-step. Therefore the length of the array must be equal to <code spellcheck="false" class="jop-noMdConv">PreprocessorParent.numOfTimeSteps</code></p></td></tr></tbody></table>

&nbsp;

**Preprocessor.setInitialCondition(v0)**  
Sets the initial values for the velocity field of the simulation. Calls the function `PreprocessorParent.setInitialConditionSkalar`.

| Parameter |     |
| --- | --- |
| `v0` : List of float or array | The values for the initial velocity field. In case it is constant in can be set by a single value for the field in x- and  a single value for the field in y- direction. In case it is different for each node it can be set via an array with the length equal to the number of nodes. An element from the array with the index i sets the velocity for the node with the index i |

&nbsp;

**Preprocessor.runSolver()**
Checks if all needed requirements are fulfilled, send the simulation data via pipe and runs the solver. Gives an live output of the simulation progress.

**Preprocessor.getNumTimeSteps()** 
Returns an integer with the number of time-steps


## Postprocessor


### 1\. Description

The Postprocessor Class reads the generated output file and offers different options to visualise the results and also make them accessible for further analysis. To improve the general usage, two classes are defined: PreprocessorParent is the parent class and provides general functions that can be used by several different solvers. The Preprocessor class already specifies functions that are required for a specific problem.

### 2\. Public Methods

**Postprocessor(projectName)**
Constructor of the Postprocessor class. Calls `PostprocessorParent(projectName)`

| Parameter |     |
| --- | --- |
| `projectName` : str | The name of the project of the Postprocessor should analyse. Must be equal to the project name defined in the Preprocessor, otherwise it will not be able to find the output file |

&nbsp;

**Postprocessor.triangulationPlotAnimated(fig,ax,name,colormap)**
Calls `PostprocessorParent.triangulationPlotAnimatedBase(fig,ax,var,name,colormap,valMin=[],valMax=[])`

<table><thead><tr><th>Parameter</th><th></th></tr></thead><tbody><tr><td><code class="inline-code" spellcheck="false">fig</code> : matplotlib Figure object</td><td>The matplotlib figure used for the plot</td></tr><tr><td><code class="inline-code" spellcheck="false">ax</code> : matplotlib Axes object</td><td>The matplotlib figure used for the plot</td></tr><tr><td><code spellcheck="false">var</code> : str</td><td><p>The variable that should be plotted:</p><ul><li><code spellcheck="false">'vel_x'</code> : velocity in x-direction</li><li><code spellcheck="false">'vel_y'</code> : velocity in y-direction</li><li><code spellcheck="false">'|vel|'</code> : norm of the velocity vector&nbsp;</li></ul></td></tr><tr><td><code class="inline-code" spellcheck="false">name</code> : str</td><td>The name of the .mp4 file to generate</td></tr><tr><td><code spellcheck="false">colormap</code> : str</td><td>The colormap used for the plot&nbsp;</td></tr><tr><td><code spellcheck="false">valMin</code> : float</td><td>Minimum value used for the colormap. If not defined, the global minimum value for the chosen variable is taken</td></tr><tr><td><code spellcheck="false">valMax</code>: float</td><td>Maximum value used for the colormap. If not defined, the global maximum value for the chosen variable is taken</td></tr></tbody></table>

&nbsp;

**Postprocessor.plotResults(fig,ax,tVal)**
Calls `PostprocessorParent.plotResultsBase(fig,ax,var,tVal,colormap,valMin=[],valMax=[])`

<table><thead><tr><th>Parameter</th><th></th></tr></thead><tbody><tr><td><code class="inline-code" spellcheck="false">fig</code> : matplotlib Figure object</td><td>The matplotlib figure used for the plot</td></tr><tr><td><code class="inline-code" spellcheck="false">ax</code> : matplotlib Axes object</td><td>The matplotlib figure used for the plot</td></tr><tr><td><code class="inline-code" spellcheck="false">var</code> : str</td><td><p>The variable that should be plotted:</p><ul><li><code class="inline-code" spellcheck="false">'vel_x'</code> : velocity in x-direction</li><li><code class="inline-code" spellcheck="false">'vel_y'</code> : velocity in y-direction</li><li><code class="inline-code" spellcheck="false">'|vel|'</code> : norm of the velocity vector</li></ul></td></tr><tr><td><code class="inline-code" spellcheck="false">tVal</code> : Array of flot</td><td>The time values at which the data should be plotted. The values must be in an range between 0 and the simulation end time</td></tr><tr><td><code class="inline-code" spellcheck="false">colormap</code> : str</td><td>The colormap used for the plot</td></tr><tr><td><code class="inline-code" spellcheck="false">valMin</code> : float</td><td>Minimum value used for the colormap. If not defined, the global minimum value for the chosen variable is taken</td></tr><tr><td><code class="inline-code" spellcheck="false">valMax</code>: float</td><td>Maximum value used for the colormap. If not defined, the global maximum value for the chosen variable is taken</td></tr></tbody></table>

**PostprocessorParent.linkMesh(mesh)**
Links a mesh object to the Postprocessor. The number of nodes from the simulation output file must be equal to the number of nodes of the imported mesh.

| Parameter |     |
| --- | --- |
| `mesh` : Mesh Object | The mesh that should be linked to the Preprocessor |
