#ifndef POINTCLOUD_MODULE_H
#define POINTCLOUD_MODULE_H

// external headers:
#include "Eigen/Core" // Eigen Core Module for Eigen::Matrix and Eigen::Array Objects
#include <iostream>   // Input Output Stream
#include <algorithm>  // for std::find
#include <iterator>   // for begin(), end(),...
#include <omp.h>      // open-mp parallel computing
#include <math.h>     // for INFINITY

// gfdm headers
#include "mesh_module.h"


class Pointcloud
{
public:

    // Attributes:
    int numInnerElements;
    Eigen::Array<int,Eigen::Dynamic,3> boundaryTriangles;
    Eigen::Matrix<double, Eigen::Dynamic, 2> coord;
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> innerElements;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> faceNormals;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeNormals;
    Eigen::Array<std::vector<int>,Eigen::Dynamic,1> support;
    double h;

    // Methods
    Pointcloud(Mesh mesh);
    void computeFaceNormals();
    void computeNodeNormals();
    void computeSupport();
};

#endif
