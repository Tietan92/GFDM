# 2D linear heat conduction solver based on the generalised finite difference method
The general finite difference method is a meshfree numerical approach to solve physical field problems. Instead of a mesh only a pointcloud of the computitional domain is needed. This offers many advantes specially for non-linear field problems in continuum mechanics. Currently, only a solver for linear heat conduction problems is implemented, but more are planned. 


## Features:
- **User-friendly pointcloud generation**: The solver is able to read finite element meshes from [Gmsh](https://gmsh.info/) and generate the point cloud out of it. 
- **Python User Interface**: The user interface was implemented entirely in Python and therefore allows intuitive and simple handling. A simulation model can be built from just a few lines of code and does not require in-depth knowledge of numerical methods.
- **Fast and efficient Solver implementation**:The entire computationally intensive part of the code was programmed in C++ with the [Eigen3 library](https://eigen.tuxfamily.org/index.php?title=Main_Page) for linear algebra. In addition, most of the code was parallelised using [openMP](https://www.openmp.org/) and can be executed on multiple processor cores. The result is a very fast solver that can also solve high-resolution simulations with more than 100,000 particles in just a few minutes.
- **Easy post processing using the whole possibilities with Python**: The simulation results can be imported into Python and therefore allow simple further processing of the data or graphical visualisation. In addition, a post-processor in Python based on [Matplotlib](https://matplotlib.org/) was implemented, which enables quick and easy visualisation of the simulation results.

## Download:



## Examples:


## Documentation:

