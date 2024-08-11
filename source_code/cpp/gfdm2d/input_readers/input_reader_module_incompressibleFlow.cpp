#include "input_reader_module_incompressibleFlow.h"

structModelInput readParameter()
{
    structModelInput modelInput;

    // Read simulation parameters
    char tmp[101];
    scanf("%s", tmp);
    modelInput.meshFile = tmp;                                      // Name of msh-File
    scanf("%lf", &modelInput.simulationParameters.frequency);       // Frequency
    scanf("%lf", &modelInput.simulationParameters.timeEnd);         // End time of the simulation
    scanf("%d",  &modelInput.simulationParameters.numNodes);        // Number of Nodes
    scanf("%d",  &modelInput.simulationParameters.numTimesteps);    // Number if time-steps
    scanf("%d",  &modelInput.simulationParameters.numBoundaries);   // Number of boundaries
    scanf("%d",  &modelInput.simulationParameters.numThreads);      // Number of threads
    scanf("%lf", &modelInput.simulationParameters.solverTolerance); // Solver tolerance

    // Read material parameters
    scanf("%lf", &modelInput.materialParameters.massDensity);  // Mass Densitiy
    scanf("%lf", &modelInput.materialParameters.dynViscosity); // Dynamic Viscosity

    // Read initial condition
    modelInput.initialConditions.startVelocity.resize(modelInput.simulationParameters.numNodes,2);
    for (int i=0; i<modelInput.simulationParameters.numNodes; i++)
    {
        scanf("%lf", &modelInput.initialConditions.startVelocity(i,0)); // start velocity v_x at node i
        scanf("%lf", &modelInput.initialConditions.startVelocity(i,1)); // start velocity v_x at node i
    }

    // Read boundary conditions
    modelInput.boundaryConditions.boundaryConditionPhysicalTag.resize(modelInput.simulationParameters.numBoundaries,1);
    modelInput.boundaryConditions.boundaryConditionType.resize(modelInput.simulationParameters.numBoundaries,1);
    modelInput.boundaryConditions.values1.resize(modelInput.simulationParameters.numBoundaries,modelInput.simulationParameters.numTimesteps);
    modelInput.boundaryConditions.values2.resize(modelInput.simulationParameters.numBoundaries,modelInput.simulationParameters.numTimesteps);

    for (int i=0; i<modelInput.simulationParameters.numBoundaries; i++)
    {
        scanf("%d", &modelInput.boundaryConditions.boundaryConditionPhysicalTag(i)); // Physical Tag
        scanf("%d", &modelInput.boundaryConditions.boundaryConditionType(i));        // boundary type

        for (int j=0; j<modelInput.simulationParameters.numTimesteps; j++)
        {
            scanf("%lf", &modelInput.boundaryConditions.values1(i,j)); // boundary values in direction 1
            scanf("%lf", &modelInput.boundaryConditions.values2(i,j)); // velocity values in direction 2
        }
    }
    return modelInput;
}
