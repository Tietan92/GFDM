#ifndef INPUT_READER_MODULE_H_INCLUDED
#define INPUT_READER_MODULE_H_INCLUDED

#include <string>     // for string objects
#include "Eigen/Core" // Eigen Core Module for Eigen::Matrix and Eigen::Array Objects
#include <unistd.h>   // For the scanf function

struct structSimulationParameters{
    double frequency;
    double timeEnd;
    int numThreads;
    double solverTolerance;
    int numNodes;
    int numTimesteps;
    int numBoundaries;
};

struct structMaterialParameters{
    double massDensity;
    double dynViscosity;
};

struct structInitialConditions{
    Eigen::Array<double, Eigen::Dynamic, 2> startVelocity;
};

struct structBoundaryConditions{
    Eigen::Array<int, Eigen::Dynamic, 1> boundaryConditionPhysicalTag;
    Eigen::Array<int, Eigen::Dynamic, 1> boundaryConditionType;
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> values1;
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> values2;
};

struct structModelInput{
    std::string meshFile;
    structSimulationParameters simulationParameters;
    structMaterialParameters materialParameters;
    structInitialConditions initialConditions;
    structBoundaryConditions boundaryConditions;
};


structModelInput readParameter();

#endif // INPUT_READER_MODULE_H_INCLUDED
