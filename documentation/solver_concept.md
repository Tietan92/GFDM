# gfdm-Toolik: Concept
![system_architecture_concept](https://github.com/user-attachments/assets/9bcc3024-38c6-4d07-9169-4f579c8ea6fa)

## Mesh Class (Python)
The Mesh class is used to import a mesh in the msh format and save it in fitting Python containers. It offers basic functionalities to visualise the mesh and to see which boundaries exists and what are the labels for the boundaries. This makes it later easier to implement the desired boundary coundditions.

## Preprocessor Class (Python)
The simulation model is constructed entirely in Python. For this a Preprocessor Class is provided. The steps here are as follows:
- Definition of simulation parameters
- Definition of material parameters
- Mesh Import
- Definition of initial conditions
- Definition of boundary conditions
Once all the preconditions have been met, the simulation can be started. To do this, the solver executable is called and all necessary data is transferred via pipe.

## Solver (C++)
The solver is programmed in C++ and uses implicit time integration. The resulting system of equation is solved using iterative solvers (GMRES and BiCGSTAB depending on the problem).
The output variables are written to the harddrive in binary form.

## Postprocessor
The Postprocessor Class reads the generated output file and offers different options to visualise the results and also make them accessible for further analysis.
