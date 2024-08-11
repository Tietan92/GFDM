#ifndef POINTCLOUD_MODULE_H
#define POINTCLOUD_MODULE_H

// external headers:
#include "Eigen/Core" // Eigen Core Module for Eigen::Matrix and Eigen::Array Objects
#include <iostream>   // Input Output Stream
#include <algorithm>  // for std::find and std::cort
#include <iterator>   // for begin(), end(),...
#include <omp.h>      // open-mp parallel computing
#include <math.h>     // for INFINITY
#include <cmath>      // for squareroot

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
    Eigen::Array<int,Eigen::Dynamic,1> supportSize;
    double h;

    // Methods
    Pointcloud(Mesh mesh);
    void computeFaceNormals();
    void computeNodeNormals();
    void computeSupport();
};

#endif
