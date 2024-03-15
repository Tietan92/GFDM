#ifndef INPUT_READER_MODULE_H_INCLUDED
#define INPUT_READER_MODULE_H_INCLUDED

#include <string>     // for string objects
#include <fstream>    // input output stream to operate on files
#include <sstream>    // string streams
#include <vector>     // for std::vector object
#include "Eigen/Core" // Eigen Core Module for Eigen::Matrix and Eigen::Array Objects

struct structSimulationParameters{
    double frequency;
    double timeEnd;
    int numThreads{1};
    double solverTolerance{1e-09};
    int numNodes;
    int numTimesteps;
    int numBoundaries;
};

struct structMaterialParameters{
    double massDensity;
    double specHeatCapacity;
    double thermalConductivity;
};

struct structInitialConditions{
    Eigen::Array<double, Eigen::Dynamic, 1> startTemperature;
};

struct structBoundaryConditions{
    Eigen::Array<int, Eigen::Dynamic, 1> boundaryConditionPhysicalTag;
    Eigen::Array<int, Eigen::Dynamic, 1> boundaryConditionType;
    Eigen::Array<double, Eigen::Dynamic, 1> filmCoefficient;
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> temperature;
    Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> heatflux;
};

struct structModelInput{
    std::string meshFile;
    structSimulationParameters simulationParameters;
    structMaterialParameters materialParameters;
    structInitialConditions initialConditions;
    structBoundaryConditions boundaryConditions;
};


structModelInput importParameterSet(std::string filePath);

#endif // INPUT_READER_MODULE_H_INCLUDED
