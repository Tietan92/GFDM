# gfdm-Toolkit : A physical solver library based on the generalized finite difference method
The generalized finite difference method is a meshfree approach to solve partial differential equations numericaly. Instead of a mesh only a pointcloud of the computitional domain is needed. This offers many advantages specially for non-linear field problems in continuum mechanics. The gfdm-Toolkit consists solvers for the following problems:
- **Dynamics of incompressible fluids (2D)**
- **Linear Heat Conduction (2D)**

## Features:
- **User-friendly pointcloud generation**: The solver is able to read finite element meshes from [Gmsh](https://gmsh.info/) and generate the point cloud out of it. 
- **Python User Interface**: The user interface was implemented entirely in Python and therefore allows intuitive and simple handling. A simulation model can be built from just a few lines of code and does not require in-depth knowledge of numerical methods.
- **Fast and efficient Solver implementation**:The entire computationally intensive part of the code was programmed in C++ with the [ViennaCL](https://viennacl.sourceforge.net/) and [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page) libraries for linear algebra. In addition, most of the code was parallelised using [openMP](https://www.openmp.org/) and can be executed on multiple processor cores. The result is a very fast solver that can also solve high-resolution simulations with more than 100,000 particles on normal computers.
- **Easy post processing using the whole possibilities with Python**: The simulation results can be imported into Python and therefore allow simple further processing of the data or graphical visualisation. In addition, a post-processor in Python based on [Matplotlib](https://matplotlib.org/) was implemented, which enables quick and easy visualisation of the simulation results.

## Documentation:
- **gfdm-Toolkit Concept**: [Link](https://github.com/Tietan92/GFDM/blob/main/documentation/solver_concept.md)
- **Dynamics of incompressible fluids**
    - Mathematical Foundation: [Link](https://raw.githack.com/Tietan92/GFDM/main/documentation/incrompressible_flow/Mathematical%20Foundation.html)
    - Software Documentation: [Link](https://github.com/Tietan92/GFDM/blob/main/documentation/incrompressible_flow/software_docu.md)
- **Linear Heat Conduction**
    - Mathematical Foundation: [Link](https://rawcdn.githack.com/Tietan92/GFDM/3c39fe69ffd24cbd256ec1d5a6d3c92fc7151adf/documentation/gfdm/1.%20Mathematical%20Foundation.html)
    - Software Documentation: [Link](https://rawcdn.githack.com/Tietan92/GFDM/3c39fe69ffd24cbd256ec1d5a6d3c92fc7151adf/documentation/gfdm/3.Python%20Software%20Architecture.html) 

## Examples:

**Turbulent flow along airfoil**

https://github.com/user-attachments/assets/9cf57511-0f52-4552-a163-a2a85eb8a975

- Go to the code page: [Link](https://github.com/Tietan92/GFDM/blob/main/examples/flow%20around%20airfoil/flow_around_airfoil.ipynb)

**Turbulent flow along square**

https://github.com/user-attachments/assets/52af2088-f746-4df5-8642-62523d965e47

- Go to the code page: [Link](https://github.com/Tietan92/GFDM/blob/main/examples/flow%20along%20square/flow%20around%20square%20big%20area.ipynb)


**CPU cooler heatflow simulation**

https://github.com/Tietan92/GFDM/assets/108583448/5428ffb2-e368-4d86-b4f7-984f10ee1599

- Go to the code page: [Link](https://github.com/Tietan92/GFDM/blob/main/examples/cpu%20cooler/cpu_cooler.ipynb)
- Step by step guide: [Link](https://rawcdn.githack.com/Tietan92/GFDM/3c39fe69ffd24cbd256ec1d5a6d3c92fc7151adf/documentation/gfdm/Example_%20Heat%20conduction%20in%20a%20cpu%20cooler.html)

**Thermal design of a micro heating element**

https://github.com/Tietan92/GFDM/assets/108583448/d494b75b-d693-488f-a1ac-787c3a4b13a2

- Go to the code page: [Link](https://github.com/Tietan92/GFDM/blob/main/examples/micro%20heating%20element/micro_heating_element.ipynb)
- Step by step guide: [Link](https://rawcdn.githack.com/Tietan92/GFDM/3c39fe69ffd24cbd256ec1d5a6d3c92fc7151adf/documentation/gfdm/Example_%20Heat%20conduction%20in%20a%20micro%20heater%20element.html)

## Download

All the required Python libs and the compiled solvers can be found in the download folder. 

## Operation system
Currently the solver is only compiled for the Linux platform. Windows will follow soon. 



