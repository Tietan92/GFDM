#include "input_reader_module_linearHeatConduction.h"

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
    scanf("%lf", &modelInput.materialParameters.massDensity);         // Mass Densitiy
    scanf("%lf", &modelInput.materialParameters.specHeatCapacity);    // Specific heat capacity
    scanf("%lf", &modelInput.materialParameters.thermalConductivity); // Thermal conductivity

    // Read initial condition
    modelInput.initialConditions.startTemperature.resize(modelInput.simulationParameters.numNodes,1);
    for (int i=0; i<modelInput.simulationParameters.numNodes; i++)
    {
        scanf("%lf", &modelInput.initialConditions.startTemperature(i)); // Start temperature at node i
    }

    // Read boundary conditions
    modelInput.boundaryConditions.boundaryConditionPhysicalTag.resize(modelInput.simulationParameters.numBoundaries,1);
    modelInput.boundaryConditions.boundaryConditionType.resize(modelInput.simulationParameters.numBoundaries,1);
    modelInput.boundaryConditions.filmCoefficient = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(modelInput.simulationParameters.numBoundaries,1);
    modelInput.boundaryConditions.heatflux.resize(modelInput.simulationParameters.numBoundaries,modelInput.simulationParameters.numTimesteps);
    modelInput.boundaryConditions.temperature.resize(modelInput.simulationParameters.numBoundaries,modelInput.simulationParameters.numTimesteps);

    for (int i=0; i<modelInput.simulationParameters.numBoundaries; i++)
    {
        scanf("%d", &modelInput.boundaryConditions.boundaryConditionPhysicalTag(i)); // Physical Tag
        scanf("%d", &modelInput.boundaryConditions.boundaryConditionType(i));        //  Boundary Type

        for (int j=0; j<modelInput.simulationParameters.numTimesteps; j++)
        {
            scanf("%lf", &modelInput.boundaryConditions.temperature(i,j)); // # Temperature Values
        }

        if (modelInput.boundaryConditions.boundaryConditionType(i) == 2) // only for Neumann boundary conditions: film coefficients and heat flux as to be imported
        {
            scanf("%lf", &modelInput.boundaryConditions.filmCoefficient(i)); // Film coefficient

            for (int j=0; j<modelInput.simulationParameters.numTimesteps; j++)
                {
                    scanf("%lf", &modelInput.boundaryConditions.heatflux(i,j)); // Heat Flux values
                }
        }
    }
    return modelInput;
}
