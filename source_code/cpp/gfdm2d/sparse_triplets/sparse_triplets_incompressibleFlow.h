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


// gfdm headers
#include "input_reader_module_incompressibleFlow.h"
#include "pointcloud_module.h"

// Global Variables
// ============================================================================

// skalar parameters
extern int N;                                   // Number of Nodes
extern int K;                                   // Number of time-steps
extern int numBoundNodes;                       // Numbber of noundary Nodes
extern double frequency;                        // Frequency of the simulation
extern double massDensity;                      // Mass density
extern double dynViscosity;                     // dynamic viscosity
extern double hSquared;                         // Square of the support radius
extern std::size_t sizeDouble; // size of double in bytes

// Eigen3 Matrix Objects
extern Eigen::Matrix<double, Eigen::Dynamic, 1> solverOutput;             // Saves the final values of the simulation output for the current timestep
extern Eigen::Matrix<double, Eigen::Dynamic, 1> B_predVel;                // RHS of the LES to calculate the predicted velocity
extern Eigen::Matrix<double, Eigen::Dynamic, 1> B_predVel_x;              // RHS of the LES to calculate the partial derivative in x-direction for the predicted velocity
extern Eigen::Matrix<double, Eigen::Dynamic, 1> B_predVel_y;              // RHS of the LES to calculate the partial derivative in y-direction for the predicted velocity
extern Eigen::Matrix<double, Eigen::Dynamic, 1> corrP_x;                  // pressure correction x-derivative
extern Eigen::Matrix<double, Eigen::Dynamic, 1> corrP_y;                  // pressure correction y-derivative
extern Eigen::Matrix<double, Eigen::Dynamic, 2> nodeCoord;                // Node coordinates
extern Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeNormals; // Node Normals

// Eigen3 Array Objects
extern Eigen::Array<std::vector<int>,Eigen::Dynamic,1> support;                                                    // Support for each node
extern Eigen::Array<int,Eigen::Dynamic,1> supportSize;                                                             // The size of the support for each node
extern Eigen::Array<unsigned int, Eigen::Dynamic, 1> nodeTypes;                                                    // the boundary condition type for each node
extern Eigen::Array<double, Eigen::Dynamic,6> corrP_rhsCoeff;                                                      // the  pressure correction RHS coefficients
extern Eigen::Array<double, Eigen::Dynamic,3> velPred_boundNodes_rhsCoeffX;                                        // predicted velocity RHS coefficients for boundary nodes (x-direction)
extern Eigen::Array<double, Eigen::Dynamic,3> velPred_boundNodes_rhsCoeffY;                                        // predicted velocity RHS coefficients for boundary nodes (y-direction)
extern Eigen::Array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,Eigen::Dynamic,1> innerNodesTailorCoeff; // Taylor series coefficients for inner nodes
extern Eigen::Array<Eigen::Matrix<double, Eigen::Dynamic, 1>,Eigen::Dynamic,1> innerNodesWeights;                  // weight function coefficients for inner nodes
extern Eigen::Array<int,Eigen::Dynamic,1> bufferOffsets;                                                           // offsets bit which indicates the start of a new sparse matrix row in the buffer
extern Eigen::Array<int,Eigen::Dynamic,1> posDiagElements;                                                         // for each node, the index of the own position in the support is saved

// viennaCL Objects
extern viennacl::compressed_matrix<double> vcl_A_predVel;   // predicted velocity sparse coefficient matrix
extern viennacl::compressed_matrix<double> vcl_A_predVel_x; // predicted velocity x-derivative sparse coefficient matrix
extern viennacl::compressed_matrix<double> vcl_A_predVel_y; // predicted velocity xyderivative sparse coefficient matrix

//std::vector Objects
extern std::vector< std::map< unsigned int, double>> corrP_triplets;    // corrected pressure sparse matrix triplets
extern std::vector< std::map< unsigned int, double>> corrP_x_triplets;  // corrected pressure x-derivate sparse matrix triplets
extern std::vector< std::map< unsigned int, double>> corrP_y_triplets;  // corrected pressure y-derivate sparse matrix triplets
extern std::vector< std::map< unsigned int, double>> corrP_xx_triplets; // corrected pressure xx-derivate sparse matrix triplets
extern std::vector< std::map< unsigned int, double>> corrP_xy_triplets; // corrected pressure xy-derivate sparse matrix triplets
extern std::vector< std::map< unsigned int, double>> corrP_yy_triplets; // corrected pressure xx-derivate sparse matrix triplets

void computeTailorCoeffInnerNodes();
void fillMat_predVel_boundNodes();
void fillMat_predVel_innerNodes(int k);
void fillMat_corrP();
#endif
