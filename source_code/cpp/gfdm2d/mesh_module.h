#ifndef MESH_H_INCLUDED
#define MESH_H_INCLUDED

#include <string>     // for string objects
#include <vector>     // for std::vector object
#include <fstream>    // input output stream to operate on files
#include "Eigen/Core" // Eigen Core Module for Eigen::Matrix and Eigen::Array Objects

struct structPhysicalNames{
    int numPhysicalNames;
    Eigen::Array<int,Eigen::Dynamic,1> dimension;
    Eigen::Array<int,Eigen::Dynamic,1> physicalTag;
    Eigen::Array<std::string,Eigen::Dynamic,1> name;
};

struct structNodes{
    int numNodes;
    Eigen::Array<double,Eigen::Dynamic,3> coord;
};

struct structElements{
    int numElements = 0;
    Eigen::Array<int,Eigen::Dynamic,1> elementTag;
    Eigen::Array<int,Eigen::Dynamic,Eigen::Dynamic> nodeTag;
};


class Mesh
{
public:
    // Attributes
    structPhysicalNames physicalNames;
    structNodes nodes;
    structElements innerElements;
    structElements boundaryElements;

    // Methods
    Mesh(std::string filePath);
};
#endif // MESH_H_INCLUDED
