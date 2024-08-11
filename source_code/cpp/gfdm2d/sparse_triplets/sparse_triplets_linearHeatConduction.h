#ifndef SPARSE_TRIPLETS_LINEARHEATCONDUCTION_H_INCLUDED
#define PARSE_TRIPLETS_LINEARHEATCONDUCTION_H_INCLUDED

// external headers:
#include <vector>       // for std::vector object
#include "Eigen/Core"   // Eigen Core Module for Eigen::Matrix and Eigen::Array Objects
#include "Eigen/Dense"  // For algebra operations on dense matrices
#include "Eigen/Sparse" // For algebra operations on sparse matrices
#include <omp.h>        // open-mp parallel computing
#include <math.h>       // for exp function

// gfdm headers
#include "input_reader_module_linearHeatConduction.h"
#include "pointcloud_module.h"

extern int N;                                                               // Number of Rows off the Output Matrix (equal to the number of Nodes)
extern int K;                                                               // Number if Columns if the Output Matrix (equal to the number of timesteps)
extern int numBoundNodes;                                                   // Number if Boundary Nodes
extern Eigen::Array<double, Eigen::Dynamic, 1> B1;                          // Array with coefficients for inner nodes
extern Eigen::Array<double, Eigen::Dynamic, 1> B2;
extern Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> boundAllocMat;
extern double frequency;

typedef Eigen::Triplet<double> Tri;


std::vector<Tri> calcTripletList(structModelInput &modelInput,Pointcloud &pointcloud, Eigen::Array<unsigned int, Eigen::Dynamic, 1> &nodeTypes);

#endif
