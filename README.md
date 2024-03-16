# 2D linear heat conduction solver based on the generalised finite difference method
The general finite difference method is a meshfree numerical approach to solve physical field problems. Instead of a mesh only a pointcloud of the computitional domain is needed. This offers many advantes specially for non-linear field problems in continuum mechanics. Currently, only a solver for linear heat conduction problems is implemented, but more are planned. 


## Features:
- **User-friendly pointcloud generation**: The solver is able to read finite element meshes from [Gmsh](https://gmsh.info/) and generate the point cloud out of it. 
- **Python User Interface**: The user interface was implemented entirely in Python and therefore allows intuitive and simple handling. A simulation model can be built from just a few lines of code and does not require in-depth knowledge of numerical methods.
- **Fast and efficient Solver implementation**:The entire computationally intensive part of the code was programmed in C++ with the [Eigen3 library](https://eigen.tuxfamily.org/index.php?title=Main_Page) for linear algebra. In addition, most of the code was parallelised using [openMP](https://www.openmp.org/) and can be executed on multiple processor cores. The result is a very fast solver that can also solve high-resolution simulations with more than 100,000 particles in just a few minutes.
- **Easy post processing using the whole possibilities with Python**: The simulation results can be imported into Python and therefore allow simple further processing of the data or graphical visualisation. In addition, a post-processor in Python based on [Matplotlib](https://matplotlib.org/) was implemented, which enables quick and easy visualisation of the simulation results.

## Documentation:
- **Mathematical foundation** of the generalised finite difference method applied to heat conduction problems: [Link](https://rawcdn.githack.com/Tietan92/GFDM/3c39fe69ffd24cbd256ec1d5a6d3c92fc7151adf/documentation/gfdm/1.%20Mathematical%20Foundation.html)
- **System Architecture**: [Link](https://rawcdn.githack.com/Tietan92/GFDM/3c39fe69ffd24cbd256ec1d5a6d3c92fc7151adf/documentation/gfdm/2.%20System%20Archtitecture.html)
- **Python Software Architecture**: [Link](https://rawcdn.githack.com/Tietan92/GFDM/3c39fe69ffd24cbd256ec1d5a6d3c92fc7151adf/documentation/gfdm/3.Python%20Software%20Architecture.html)
- **C++ Software Architecture**: [Link](https://rawcdn.githack.com/Tietan92/GFDM/3c39fe69ffd24cbd256ec1d5a6d3c92fc7151adf/documentation/gfdm/4.%20C++%20Software%20Architecture.html)

## Examples:

**CPU cooler heatflow simulation**: 

https://github.com/Tietan92/GFDM/assets/108583448/5428ffb2-e368-4d86-b4f7-984f10ee1599

- Go to the code page [Link](https://github.com/Tietan92/GFDM/blob/main/examples/cpu%20cooler/cpu_cooler.ipynb)
- Step by step guide: [Link](https://rawcdn.githack.com/Tietan92/GFDM/3c39fe69ffd24cbd256ec1d5a6d3c92fc7151adf/documentation/gfdm/Example_%20Heat%20conduction%20in%20a%20cpu%20cooler.html)

**Thermal design of a micro heating element**

https://github.com/Tietan92/GFDM/assets/108583448/d494b75b-d693-488f-a1ac-787c3a4b13a2

- Go to the code page [Link](https://github.com/Tietan92/GFDM/blob/main/examples/micro%20heating%20element/micro_heating_element.ipynb)
- Step by step guide: [Link](https://rawcdn.githack.com/Tietan92/GFDM/3c39fe69ffd24cbd256ec1d5a6d3c92fc7151adf/documentation/gfdm/Example_%20Heat%20conduction%20in%20a%20micro%20heater%20element.html)

## Download

All the required Python libs and the compiled solvers can be found in the download folder. 



