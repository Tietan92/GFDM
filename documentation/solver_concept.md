![system_architecture_concept](https://github.com/user-attachments/assets/9bcc3024-38c6-4d07-9169-4f579c8ea6fa)

The simulation model is constructed entirely in Python. For this a Preprocessor Class is provided. The steps here are as follows:
- Definition of simulation parameters
- Definition of material parameters
- Mesh Import
- Definition of initial conditions
- Definition of boundary conditions

Once all the preconditions have been met, the simulation can be started. To do this, the solver executable is called and all necessary data is transferred via pipe.
The solver is programmed in C++ and uses implicit time integration. The resulting system of equation is solved using iterative solvers (GMRES and BiCGSTAB depending on the problem).
The output variables are written to the harddrive in binary form.

For the visualization of the results a Postprocessor Class for Python is provided. It can plot the results at a specific time point and can also generate animations of the whole simulation.  
