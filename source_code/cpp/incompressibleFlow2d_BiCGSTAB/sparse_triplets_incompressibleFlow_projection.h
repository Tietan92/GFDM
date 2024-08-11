#ifndef SPARSE_TRIPLETS_LINEARHEATCONDUCTION_H_INCLUDED
#define PARSE_TRIPLETS_LINEARHEATCONDUCTION_H_INCLUDED

// external headers:
#include <vector>       // for std::vector object
#include "Eigen/Core"   // Eigen Core Module for Eigen::Matrix and Eigen::Array Objects
#include "Eigen/Dense"  // For algebra operations on dense matrices
#include "Eigen/Sparse" // For algebra operations on sparse matrices
#include <omp.h>        // open-mp parallel computing
#include <math.h>       // for exp function

// ViennaCL headers
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/io/matrix_market.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"

// gfdm headers
#include "input_reader_module_incompressibleFlow.h"
#include "pointcloud_module.h"

// Global Variables
// ============================================================================

// skalar parameters
extern int N;               // Number of Nodes
extern int K;               // Number of time-steps
extern int numBoundNodes;   // Numbber of noundary Nodes
extern double frequency;    // Frequency of the simulation
extern double massDensity;  // Mass density
extern double dynViscosity; // dynamic viscosity
extern double hSquared;     //
extern std::size_t sizeDouble;

// Eigen3 Matrix Objects
extern Eigen::Matrix<double, Eigen::Dynamic, 1> solverOutput;
extern Eigen::Matrix<double, Eigen::Dynamic, 1> B_predVel;
extern Eigen::Matrix<double, Eigen::Dynamic, 1> B_predVel_x;
extern Eigen::Matrix<double, Eigen::Dynamic, 1> B_predVel_y;
extern Eigen::Matrix<double, Eigen::Dynamic, 1> corrP_x;
extern Eigen::Matrix<double, Eigen::Dynamic, 1> corrP_y;
extern Eigen::Matrix<double, Eigen::Dynamic, 2> nodeCoord;                      // Node coordinates
extern Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeNormals; // Node Normals

// Eigen3 Array Objects
extern Eigen::Array<std::vector<int>,Eigen::Dynamic,1> support;               // Support for each node
extern Eigen::Array<int,Eigen::Dynamic,1> supportSize;
extern Eigen::Array<unsigned int, Eigen::Dynamic, 1> nodeTypes;
extern Eigen::Array<double, Eigen::Dynamic,6> corrP_rhsCoeff;
extern Eigen::Array<double, Eigen::Dynamic,3> velPred_boundNodes_rhsCoeffX; // right hand side coefficient for boundary nodes (x-direction)
extern Eigen::Array<double, Eigen::Dynamic,3> velPred_boundNodes_rhsCoeffY; // right hand side coefficient for boundary nodes (y-direction)
extern Eigen::Array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,Eigen::Dynamic,1> innerNodesTailorCoeff;
extern Eigen::Array<Eigen::Matrix<double, Eigen::Dynamic, 1>,Eigen::Dynamic,1> innerNodesWeights;
extern Eigen::Array<int,Eigen::Dynamic,1> bufferOffsets;
extern Eigen::Array<int,Eigen::Dynamic,1> posDiagElements;

// viennaCL Objects
extern viennacl::compressed_matrix<double> vcl_A_predVel;
extern viennacl::compressed_matrix<double> vcl_A_predVel_x;
extern viennacl::compressed_matrix<double> vcl_A_predVel_y;

//std::vector Objects
extern std::vector< std::map< unsigned int, double>> corrP_triplets;
extern std::vector< std::map< unsigned int, double>> corrP_x_triplets;
extern std::vector< std::map< unsigned int, double>> corrP_y_triplets;
extern std::vector< std::map< unsigned int, double>> corrP_xx_triplets;
extern std::vector< std::map< unsigned int, double>> corrP_xy_triplets;
extern std::vector< std::map< unsigned int, double>> corrP_yy_triplets;

void computeTailorCoeffInnerNodes();
void fillMat_predVel_boundNodes();
void fillMat_predVel_innerNodes(int k);
void fillMat_corrP();
#endif
