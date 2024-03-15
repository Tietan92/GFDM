#ifndef BOUNDARY_CONDITION_MODULE_H_INCLUDED
#define BOUNDARY_CONDITION_MODULE_H_INCLUDED

// external headers:
#include "Eigen/Core" // Eigen Core Module for Eigen::Matrix and Eigen::Array Objects
#include <omp.h>      // open-mp parallel computing

// gfdm headers
#include "mesh_module.h"
#include "input_reader_module_linearHeatConduction.h"

extern Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> boundAllocMat;
extern int N;
extern int numBoundNodes;


void allocateBoundaryConditions(Mesh &mesh, structModelInput &modelInput, Eigen::Array<unsigned int, Eigen::Dynamic, 1> &nodeTypes);

#endif // BOUNDARY_CONDITION_MODULE_H_INCLUDED
