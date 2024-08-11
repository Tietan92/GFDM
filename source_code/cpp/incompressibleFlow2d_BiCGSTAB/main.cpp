// stl headers
// ============================================================================
#include <string>          // for string objects
#include <fstream>         // input output stream to operate on files
#include  <iostream>       // std::cout
#include <map>             // std::map

// omp and openCL headers
// ============================================================================
#include <omp.h>           // open-mp parallel computing

// Preprocessor Variables
// ============================================================================
#define VIENNACL_WITH_EIGEN 1  // Enable Interface between Eigen3 and ViennaCL
#define VIENNACL_WITH_OPENMP 1 // Using ViennaCL with OpenMP
#define EIGEN_NO_DEBUG 1       // Disables Eigen3 assertions

// Eigen headers
// ============================================================================
#include "../Eigen/Core"   // Eigen Core Module for Eigen::Matrix and Eigen::Array Objects
#include "../Eigen/Sparse" // For algebra operations on sparse matrices

// viennaCL headers
// ============================================================================
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/io/matrix_market.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"

// gfdm headers
// ============================================================================
#include "input_reader_module_incompressibleFlow.h"
#include "mesh_module.h"
#include "pointcloud_module.h"
#include "boundary_condition_module.h"
#include "output_writer_module.h"
#include "sparse_triplets_incompressibleFlow.h"

// Global Variables
// ============================================================================

// skalar parameters
int numVars = 3;                         // Number of Variables: Velocity in x and y direction and pressure
int N;                                   // Number of Nodes
int K;                                   // Number of time-steps
int numBoundNodes;                       // Numbber of noundary Nodes
double frequency;                        // Frequency of the simulation
double massDensity;                      // Mass density
double dynViscosity;                     // dynamic viscosity
double hSquared;                         // Square of the support radius
std::size_t sizeDouble = sizeof(double); // size of double in bytes

// Eigen3 Matrix Objects
Eigen::Matrix<double, Eigen::Dynamic, 1> solverOutput;             // Saves the final values of the simulation output for the current timestep
Eigen::Matrix<double, Eigen::Dynamic, 1> B_predVel;                // RHS of the LES to calculate the predicted velocity
Eigen::Matrix<double, Eigen::Dynamic, 1> B_predVel_x;              // RHS of the LES to calculate the partial derivative in x-direction for the predicted velocity
Eigen::Matrix<double, Eigen::Dynamic, 1> B_predVel_y;              // RHS of the LES to calculate the partial derivative in y-direction for the predicted velocity
Eigen::Matrix<double, Eigen::Dynamic, 1> corrP_x;                  // pressure correction x-derivative
Eigen::Matrix<double, Eigen::Dynamic, 1> corrP_y;                  // pressure correction y-derivative
Eigen::Matrix<double, Eigen::Dynamic, 2> nodeCoord;                // Node coordinates
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeNormals; // Node Normals

// Eigen3 Array Objects
Eigen::Array<std::vector<int>,Eigen::Dynamic,1> support;                                                    // Support for each node
Eigen::Array<int,Eigen::Dynamic,1> supportSize;                                                             // The size of the support for each node
Eigen::Array<unsigned int, Eigen::Dynamic, 1> nodeTypes;                                                    // the boundary condition type for each node
Eigen::Array<double, Eigen::Dynamic,6> corrP_rhsCoeff;                                                      // the  pressure correction RHS coefficients
Eigen::Array<double, Eigen::Dynamic,3> velPred_boundNodes_rhsCoeffX;                                        // predicted velocity RHS coefficients for boundary nodes (x-direction)
Eigen::Array<double, Eigen::Dynamic,3> velPred_boundNodes_rhsCoeffY;                                        // predicted velocity RHS coefficients for boundary nodes (y-direction)
Eigen::Array<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,Eigen::Dynamic,1> innerNodesTailorCoeff; // Taylor series coefficients for inner nodes
Eigen::Array<Eigen::Matrix<double, Eigen::Dynamic, 1>,Eigen::Dynamic,1> innerNodesWeights;                  // weight function coefficients for inner nodes
Eigen::Array<int,Eigen::Dynamic,1> bufferOffsets;                                                           // offsets bit which indicates the start of a new sparse matrix row in the buffer
Eigen::Array<int,Eigen::Dynamic,1> posDiagElements;                                                         // for each node, the index of the own position in the support is saved

// viennaCL Objects
viennacl::compressed_matrix<double> vcl_A_predVel;   // predicted velocity sparse coefficient matrix
viennacl::compressed_matrix<double> vcl_A_predVel_x; // predicted velocity x-derivative sparse coefficient matrix
viennacl::compressed_matrix<double> vcl_A_predVel_y; // predicted velocity xyderivative sparse coefficient matrix

//std::vector Objects
std::vector< std::map< unsigned int, double>> corrP_triplets;    // corrected pressure sparse matrix triplets
std::vector< std::map< unsigned int, double>> corrP_x_triplets;  // corrected pressure x-derivate sparse matrix triplets
std::vector< std::map< unsigned int, double>> corrP_y_triplets;  // corrected pressure y-derivate sparse matrix triplets
std::vector< std::map< unsigned int, double>> corrP_xx_triplets; // corrected pressure xx-derivate sparse matrix triplets
std::vector< std::map< unsigned int, double>> corrP_xy_triplets; // corrected pressure xy-derivate sparse matrix triplets
std::vector< std::map< unsigned int, double>> corrP_yy_triplets; // corrected pressure xx-derivate sparse matrix triplets

// other object types
std::ofstream outputFile;

// main function
// ============================================================================
int main(int argc, char *argv[])
{
    std::string name             = argv[1];               // The first console argument tells the project name
    std::string workingDirectory = argv[2];               // The second console argument tells the working directory
    std::string outputName       = name + "_output.dat";  // naming convention for the output file

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> boundAllocMat; // boundary allocation matrix
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> boundValMatX;  // boundary value matrix x-direction
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> boundValMatY;  // boundary values matrix y-direction

    viennacl::linalg::bicgstab_tag solverTagBicgstab; // solver tag for bicgstab solver --> used for LES of predicted velocity
    viennacl::linalg::gmres_tag solverTagGmres;       // solver tag for gmres solver --> used for LES of pressure correction

    {   // Begin Local Scope
        structModelInput modelInput = readParameter(); // read input

        Mesh mesh(workingDirectory + "/" + modelInput.meshFile); // load mesh

        omp_set_num_threads(modelInput.simulationParameters.numThreads); // set the number of threads for multiphreading with open mp

        // Assign Parameters
        N             = mesh.nodes.numNodes;
        K             = modelInput.simulationParameters.numTimesteps;
        numBoundNodes = mesh.boundaryElements.numElements;
        frequency     = modelInput.simulationParameters.frequency;
        massDensity   = modelInput.materialParameters.massDensity;
        dynViscosity  = modelInput.materialParameters.dynViscosity;
        boundValMatX  = modelInput.boundaryConditions.values1;
        boundValMatY  = modelInput.boundaryConditions.values2;

        // Setting up solver tags
        solverTagBicgstab = viennacl::linalg::bicgstab_tag(modelInput.simulationParameters.solverTolerance);
        solverTagGmres = viennacl::linalg::gmres_tag(modelInput.simulationParameters.solverTolerance);


        // resizing of global variables
        solverOutput.resize(3*N);
        B_predVel.resize(2*N);
        B_predVel_x.resize(2*N);
        B_predVel_y.resize(2*N);
        corrP_x.resize(N);
        corrP_y.resize(N);
        velPred_boundNodes_rhsCoeffX.resize(2*numBoundNodes,3);
        velPred_boundNodes_rhsCoeffY.resize(2*numBoundNodes,3);
        corrP_rhsCoeff.resize(N,6);
        innerNodesTailorCoeff.resize(N-numBoundNodes);
        innerNodesWeights.resize(N-numBoundNodes);
        bufferOffsets.resize(2*N);
        posDiagElements.resize(N);
        corrP_triplets.resize(N);
        corrP_x_triplets.resize(N);
        corrP_y_triplets.resize(N);
        corrP_xx_triplets.resize(N);
        corrP_xy_triplets.resize(N);
        corrP_yy_triplets.resize(N);

        allocateBoundaryConditions(mesh, modelInput, nodeTypes, boundAllocMat); // Allocate the boundary conditions

        // load pointcloud
        Pointcloud pointcloud(mesh);
        support = pointcloud.support;
        nodeCoord = pointcloud.coord;
        nodeNormals = pointcloud.nodeNormals;
        hSquared = pointcloud.h*pointcloud.h;
        supportSize = pointcloud.supportSize;

        // Setting Initial Conditions
        solverOutput(Eigen::seqN(0,N)) = modelInput.initialConditions.startVelocity(Eigen::all,0);
        solverOutput(Eigen::seqN(N,N)) = modelInput.initialConditions.startVelocity(Eigen::all,1);
        corrP_x = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(N);
        corrP_y = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(N);

        // Analyze Matrix Patterns
        std::vector< std::map< unsigned int, double>> predVel_SparsePatterns(2*N);
        for (int i = 0; i<N; i++)
        {
            for (int j : support(i))
            {
                predVel_SparsePatterns[i][j]     = 1;
                predVel_SparsePatterns[i][j+N]   = 1;
                predVel_SparsePatterns[i+N][j]   = 1;
                predVel_SparsePatterns[i+N][j+N] = 1;
            }
        }
        viennacl::copy(predVel_SparsePatterns, vcl_A_predVel);
        viennacl::copy(predVel_SparsePatterns, vcl_A_predVel_x);
        viennacl::copy(predVel_SparsePatterns, vcl_A_predVel_y);

        //Calculate Buffer Offsets
        bufferOffsets(0) = 0;
        for (int i=0; i<N; i++)
        {
            bufferOffsets(i+1) = bufferOffsets(i) + supportSize(i) * 2;
        }
        for (int i=0; i<N-1; i++)
        {
            bufferOffsets(i+1+N) = bufferOffsets(i+N) + supportSize(i) * 2;
        }
        bufferOffsets = bufferOffsets * sizeDouble;

        // calculate the positions of the diagonal elements
        for (int i=0; i<N; i++)
        {
            for (int j=0; j<supportSize(i); j++)
            {

                if (support(i)[j] == i)
                {
                    posDiagElements(i) = j;
                    break;
                }
            }
        }
    } // End of local Scope

    // Initialise Eigen Arrays
    Eigen::Matrix<double, Eigen::Dynamic, 1> corrP(N);      // pressure correction
    Eigen::Matrix<double, Eigen::Dynamic, 1> corrP_xx(N);   // pressure correction xx-derivative
    Eigen::Matrix<double, Eigen::Dynamic, 1> corrP_xy(N);   // pressure correction xy-derivative
    Eigen::Matrix<double, Eigen::Dynamic, 1> corrP_yy(N);   // pressure correction yy-derivative

    Eigen::Matrix<double, Eigen::Dynamic, 1> B_corrP(N);    // pressure correction RHS of LES
    Eigen::Matrix<double, Eigen::Dynamic, 1> B_corrP_x(N);  // pressure correction x-derivative RHS of LES
    Eigen::Matrix<double, Eigen::Dynamic, 1> B_corrP_y(N);  // pressure correction x-derivative RHS of LES
    Eigen::Matrix<double, Eigen::Dynamic, 1> B_corrP_xx(N); // pressure correction xx-derivative RHS of LES
    Eigen::Matrix<double, Eigen::Dynamic, 1> B_corrP_xy(N); // pressure correction xy-derivative RHS of LES
    Eigen::Matrix<double, Eigen::Dynamic, 1> B_corrP_yy(N); // pressure correction yy-derivative RHS of LES

    Eigen::Matrix<double, Eigen::Dynamic, 1> predVel(2*N);   // predicted velocity
    Eigen::Matrix<double, Eigen::Dynamic, 1> predVel_x(2*N); // predicted velocity x-derivative
    Eigen::Matrix<double, Eigen::Dynamic, 1> predVel_y(2*N); // predicted veloctiy y-derivative

    Eigen::Matrix<double, Eigen::Dynamic, 1> corrVel(2*N);   // velocity correction
    Eigen::Matrix<double, Eigen::Dynamic, 1> corrVel_x(2*N); // velocity correction x-derivative
    Eigen::Matrix<double, Eigen::Dynamic, 1> corrVel_y(2*N); // velocity correction y-derivative

    Eigen::Array<double, Eigen::Dynamic, 1> boundValuesX = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(numBoundNodes); // saves the boundary values if the current timestep in x-direction
    Eigen::Array<double, Eigen::Dynamic, 1> boundValuesY = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(numBoundNodes); // saves the boundary values if the current timestep in y-direction

    Eigen::Array<double, Eigen::Dynamic, 1> valuesMinMax = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(6); // [vMax, wMax, pMax, vMin, wMin, pMin]
    Eigen::Array<double, 6, 1> valuesMinMaxTemp = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(6); // [vMax, wMax, pMax, vMin, wMin, pMin]

    fillMat_predVel_boundNodes(); // fill the predVel matrices for boundary nodes

    computeTailorCoeffInnerNodes(); // Pre-computing if the tailer coeffićients of the inner nodes

    // Initialise viennaCL Objects
    viennacl::vector<double> vcl_corrP(N);         // pressure correction
    viennacl::vector<double> vcl_corrP_x(N);       // pressure correction x-derivative
    viennacl::vector<double> vcl_corrP_y(N);       // pressure correction y-derivative
    viennacl::vector<double> vcl_corrP_xx(N);      // pressure correction xx-derivative
    viennacl::vector<double> vcl_corrP_xy(N);      // pressure correction xy-derivative
    viennacl::vector<double> vcl_corrP_yy(N);      // pressure correction yy-derivative

    viennacl::vector<double> vcl_B_corrP(N);       // pressure correction RHS
    viennacl::vector<double> vcl_B_corrP_x(N);     // pressure correction x-derivative RHS
    viennacl::vector<double> vcl_B_corrP_y(N);     // pressure correction y-derivative RHS
    viennacl::vector<double> vcl_B_corrP_xx(N);    // pressure correction xx-derivative RHS
    viennacl::vector<double> vcl_B_corrP_xy(N);    // pressure correction xy-derivative RHS
    viennacl::vector<double> vcl_B_corrP_yy(N);    // pressure correction yy-derivative RHS

    viennacl::vector<double> vcl_predVel(2*N);     // predicted velocity
    viennacl::vector<double> vcl_predVel_x(2*N);   // predicted velocity x-derivative
    viennacl::vector<double> vcl_predVel_y(2*N);   // predicted velocity y-derivative

    viennacl::vector<double> vcl_B_predVel(2*N);   // predicted velocity RHS
    viennacl::vector<double> vcl_B_predVel_x(2*N); // predicted velocity x-derivative RHS
    viennacl::vector<double> vcl_B_predVel_y(2*N); // predicted velocity y-derivative RHS

    viennacl::compressed_matrix<double> vcl_A_corrP(N,N);    // pressure correction sparse matrix
    viennacl::compressed_matrix<double> vcl_A_corrP_x(N,N);  // pressure correction x-derivative sparse matrix
    viennacl::compressed_matrix<double> vcl_A_corrP_y(N,N);  // pressure correction y-derivatie sparse matrix
    viennacl::compressed_matrix<double> vcl_A_corrP_xx(N,N); // pressure correction xx-derivative  sparse matrix
    viennacl::compressed_matrix<double> vcl_A_corrP_xy(N,N); // pressure correction xy-derivative sparse matrix
    viennacl::compressed_matrix<double> vcl_A_corrP_yy(N,N); // pressure correction yy-derivative sparse matrix

    // Fill Matrices for pressure correction
    fillMat_corrP();
    viennacl::copy(corrP_triplets,    vcl_A_corrP);
    viennacl::copy(corrP_x_triplets,  vcl_A_corrP_x);
    viennacl::copy(corrP_y_triplets,  vcl_A_corrP_y);
    viennacl::copy(corrP_xx_triplets, vcl_A_corrP_xx);
    viennacl::copy(corrP_xy_triplets, vcl_A_corrP_xy);
    viennacl::copy(corrP_yy_triplets, vcl_A_corrP_yy);

    // Set Up Preconditioners:
    viennacl::linalg::jacobi_precond< viennacl::compressed_matrix<double>> vcl_predVel_Precond(vcl_A_predVel, viennacl::linalg::jacobi_tag()); // Jacobi preconditioner for LES of predicted velocity
    viennacl::linalg::jacobi_precond< viennacl::compressed_matrix<double>> vcl_coorP_Precond(vcl_A_corrP, viennacl::linalg::jacobi_tag());     // Jacobi preconditioner for LES of pressure correction

    // setting up solvers
    viennacl::linalg::bicgstab_solver<viennacl::vector<double>> vlc_predVelSolver(solverTagBicgstab); // bicgstab solver for LES of predicted velocity
    viennacl::linalg::gmres_solver<viennacl::vector<double>> vlc_corrPSolver(solverTagGmres);         // gmres solver for LES of pressure correction

    // Setting Up Output File
    generateOutputFile(workingDirectory + "/" + outputName);
    int percentage = 10; // Percentage of the simulation progress

    for (int k=0; k<K; k++)
    {
        if (k>0)
        {
            // Predictor Step: Calculate Velocity Prediction ----------------------------------------------------------

            // Calculate right-hand-side  coefficients for boundary nodes:
            boundValuesX = (boundAllocMat*boundValMatX(Eigen::all,k)).array();
            boundValuesY = (boundAllocMat*boundValMatY(Eigen::all,k)).array();

            #pragma omp parallel for
            for (int i = 0; i<numBoundNodes; i++)
            {
                B_predVel(i)     = velPred_boundNodes_rhsCoeffX(i,0)              * boundValuesX(i) + velPred_boundNodes_rhsCoeffY(i,0)               * boundValuesY(i);
                B_predVel(i+N)   = velPred_boundNodes_rhsCoeffX(i+numBoundNodes,0)* boundValuesX(i) + velPred_boundNodes_rhsCoeffY(i+numBoundNodes,0) * boundValuesY(i);
                B_predVel_x(i)   = velPred_boundNodes_rhsCoeffX(i,1)              * boundValuesX(i) + velPred_boundNodes_rhsCoeffY(i,1)               * boundValuesY(i);
                B_predVel_x(i+N) = velPred_boundNodes_rhsCoeffX(i+numBoundNodes,1)* boundValuesX(i) + velPred_boundNodes_rhsCoeffY(i+numBoundNodes,1) * boundValuesY(i);
                B_predVel_y(i)   = velPred_boundNodes_rhsCoeffX(i,2)              * boundValuesX(i) + velPred_boundNodes_rhsCoeffY(i,2)               * boundValuesY(i);
                B_predVel_y(i+N) = velPred_boundNodes_rhsCoeffX(i+numBoundNodes,2)* boundValuesX(i) + velPred_boundNodes_rhsCoeffY(i+numBoundNodes,2) * boundValuesY(i);
            }
            // Fill matrices and calculate right-hand side coefficients for inner nodes
            fillMat_predVel_innerNodes(k);

            // Copy eigen to vcl:
            viennacl::copy(B_predVel, vcl_B_predVel);
            viennacl::copy(B_predVel_x, vcl_B_predVel_x);
            viennacl::copy(B_predVel_y, vcl_B_predVel_y);

            // Solve the linear Systems
            vcl_predVel_Precond.init(vcl_A_predVel);
            vlc_predVelSolver.set_initial_guess(vcl_predVel);
            vcl_predVel   = vlc_predVelSolver(vcl_A_predVel,vcl_B_predVel,vcl_predVel_Precond);
            vcl_predVel_x = viennacl::linalg::prod(vcl_A_predVel_x, vcl_predVel) + vcl_B_predVel_x;
            vcl_predVel_y = viennacl::linalg::prod(vcl_A_predVel_y, vcl_predVel) + vcl_B_predVel_y;

            // Copy back to eigen
            viennacl::copy(vcl_predVel, predVel);
            viennacl::copy(vcl_predVel_x, predVel_x);
            viennacl::copy(vcl_predVel_y, predVel_y);

            // First Corrector Step: Calculate first Pressure Correction ----------------------------------------------
            B_corrP(Eigen::seqN(0,numBoundNodes)) = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(numBoundNodes); // right-hand-side coefficients for boundary nodes

            #pragma omp parallel for
            for (int i = numBoundNodes; i<N;i++) // right-hand-side coefficients for inner nodes
            {
                double rhs = massDensity * frequency * (predVel_x(i)  + predVel_y(i+N));
                B_corrP(i)    = corrP_rhsCoeff(i,0) * rhs;
                B_corrP_x(i)  = corrP_rhsCoeff(i,1) * rhs;
                B_corrP_y(i)  = corrP_rhsCoeff(i,2) * rhs;
                B_corrP_xx(i) = corrP_rhsCoeff(i,3) * rhs;
                B_corrP_xy(i) = corrP_rhsCoeff(i,4) * rhs;
                B_corrP_yy(i) = corrP_rhsCoeff(i,5) * rhs;
            }
            B_corrP = B_corrP.array() - B_corrP.sum()/N; // Subtracting the mean of the right-hand-side vector so ensure discrete compatibility

            // Copy eigen to vcl:
            viennacl::copy(B_corrP, vcl_B_corrP);
            viennacl::copy(B_corrP_x, vcl_B_corrP_x);
            viennacl::copy(B_corrP_y, vcl_B_corrP_y);
            viennacl::copy(B_corrP_xx, vcl_B_corrP_xx);
            viennacl::copy(B_corrP_xy, vcl_B_corrP_xy);
            viennacl::copy(B_corrP_yy, vcl_B_corrP_yy);

            // Solve the linear Systems for current timestep to calculate the first pressure correction and its spartial derivates
            vcl_corrP = vlc_corrPSolver(vcl_A_corrP,vcl_B_corrP,vcl_coorP_Precond);
            vcl_corrP_x = viennacl::linalg::prod(vcl_A_corrP_x, vcl_corrP) + vcl_B_corrP_x;
            vcl_corrP_y = viennacl::linalg::prod(vcl_A_corrP_y, vcl_corrP) + vcl_B_corrP_y;
            vcl_corrP_xx = viennacl::linalg::prod(vcl_A_corrP_xx, vcl_corrP) + vcl_B_corrP_xx;
            vcl_corrP_xy = viennacl::linalg::prod(vcl_A_corrP_xy, vcl_corrP) + vcl_B_corrP_xy;
            vcl_corrP_yy = viennacl::linalg::prod(vcl_A_corrP_yy, vcl_corrP) + vcl_B_corrP_yy;

            // copy back to eigen
            viennacl::copy(vcl_corrP, corrP);
            viennacl::copy(vcl_corrP_x, corrP_x);
            viennacl::copy(vcl_corrP_y, corrP_y);
            viennacl::copy(vcl_corrP_xx, corrP_xx);
            viennacl::copy(vcl_corrP_xy, corrP_xy);
            viennacl::copy(vcl_corrP_yy, corrP_yy);

            // Perform first velocity correction for inner nodes:
            #pragma omp parallel for
            for (int i = numBoundNodes; i<N;i++)
            {
                corrVel(i)    = predVel(i)    - 1/(massDensity*frequency) * corrP_x(i);
                corrVel_x(i)  = predVel_x(i)  - 1/(massDensity*frequency) * corrP_xx(i);
                corrVel_y(i)  = predVel_y(i)  - 1/(massDensity*frequency) * corrP_xy(i);

                corrVel(i+N)    = predVel(i+N)    - 1/(massDensity*frequency) * corrP_y(i);
                corrVel_x(i+N)  = predVel_x(i+N)  - 1/(massDensity*frequency) * corrP_xy(i);
                corrVel_y(i+N)  = predVel_y(i+N)  - 1/(massDensity*frequency) * corrP_yy(i);
            }
            solverOutput(Eigen::seqN(2*N,N)) += corrP;
            // Second Corrector Step: Calculate second pressure correction and final velocity -------------------------

            B_corrP(Eigen::seqN(0,numBoundNodes)) = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(numBoundNodes); // right-hand-side coefficients for boundary nodes

            #pragma omp parallel for
            for (int i = numBoundNodes; i<N;i++) // right-hand-side coefficients for inner nodes
            {
                double rhs =  massDensity * (corrVel_x(i)*corrVel_x(i) + 2 * corrVel_x(i+N)*corrVel_y(i) + corrVel_y(i+N)*corrVel_y(i+N))
                            - massDensity * (predVel_x(i)*predVel_x(i)   + 2 * predVel_x(i+N)*predVel_y(i) + predVel_y(i+N)*predVel_y(i+N));

                B_corrP(i)    = corrP_rhsCoeff(i,0) * rhs;
                B_corrP_x(i)  = corrP_rhsCoeff(i,1) * rhs;
                B_corrP_y(i)  = corrP_rhsCoeff(i,2) * rhs;
            }

            B_corrP = B_corrP.array() - B_corrP.sum()/N; // Subtracting the mean of the right-hand-side vector so ensure discrete compatibility

            // Copy eigen to vcl:
            viennacl::copy(B_corrP, vcl_B_corrP);
            viennacl::copy(B_corrP_x, vcl_B_corrP_x);
            viennacl::copy(B_corrP_y, vcl_B_corrP_y);

            // Solve the linear Systems for current timestep to calculate the first pressure correction and its spartial derivates
            vcl_corrP = vlc_corrPSolver(vcl_A_corrP,vcl_B_corrP,vcl_coorP_Precond);
            vcl_corrP_x = viennacl::linalg::prod(vcl_A_corrP_x, vcl_corrP) + vcl_B_corrP_x;
            vcl_corrP_y = viennacl::linalg::prod(vcl_A_corrP_y, vcl_corrP) + vcl_B_corrP_y;

            // copy back to eigen
            viennacl::copy(vcl_corrP, corrP);
            viennacl::copy(vcl_corrP_x, corrP_x);
            viennacl::copy(vcl_corrP_y, corrP_y);

            // Update finale velocity and pressure
            solverOutput(Eigen::seqN(0,numBoundNodes)) = predVel(Eigen::seqN(0,numBoundNodes)); // boundary values dont need correction and are therefore taken from the prediction step
            solverOutput(Eigen::seqN(N,numBoundNodes)) = predVel(Eigen::seqN(N,numBoundNodes));
            #pragma omp parallel for
            for (int i = numBoundNodes; i<N;i++)
            {
                solverOutput(i)   = corrVel(i)   - 1/frequency * (corrVel(i) * corrVel_x(i)   + corrVel(i+N) * corrVel_y(i))
                                                 + 1/frequency * (predVel(i) * predVel_x(i)   + predVel(i+N) * predVel_y(i))
                                                 - 1/(massDensity*frequency) * corrP_x(i);

                solverOutput(i+N) = corrVel(i+N) - 1/frequency * (corrVel(i) * corrVel_x(i+N) + corrVel(i+N) * corrVel_y(i+N))
                                                 + 1/frequency * (predVel(i) * predVel_x(i+N) + predVel(i+N) * predVel_y(i+N))
                                                 - 1/(massDensity*frequency) * corrP_y(i);
            }
            //solverOutput(Eigen::seqN(2*N,N)) += corrP;
            for (int i = 0; i<N;i++)
            {
                Eigen::Matrix<double, 2, 1> diff;
                diff(0) = solverOutput(i);
                diff(1) = solverOutput(i+N);
                solverOutput(i+2*N) = diff.norm();
            }

        }
        // Calculate the min-max values for the current timestep
        #pragma omp parallel for
        for (int i=0; i<6; i++)
        {
            if (i<3)
            {
                valuesMinMaxTemp(i) = solverOutput(Eigen::seqN(i*N,N)).array().maxCoeff();

                if (valuesMinMaxTemp(i)> valuesMinMax(i))
                {
                    valuesMinMax(i) = valuesMinMaxTemp(i);
                }
            }
            else
            {
                valuesMinMaxTemp(i) = solverOutput(Eigen::seqN((i-3)*N,N)).array().minCoeff();

                if (valuesMinMaxTemp(i) < valuesMinMax(i))
                {
                    valuesMinMax(i) = valuesMinMaxTemp(i);
                }
            }
        }
        writeToOutputFile(k, percentage);
    }
    closeOutputFile(valuesMinMax);
}
