// stl headers
// ============================================================================
#include <string>          // for string objects
#include <fstream>         // input output stream to operate on files
#include  <iostream>       // std::cout
#include <map>             // std::map
#include <iomanip>

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
#include "Eigen/Dense"  // For algebra operations on dense matrices

// viennaCL headers
// ============================================================================
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/io/matrix_market.hpp"
#include "viennacl/linalg/jacobi_precond.hpp"

// gfdm headers
// ============================================================================
#include "input_reader_module_elasticSolids.h"
#include "mesh_module.h"
#include "pointcloud_module.h"
#include "boundary_condition_module.h"
#include "output_writer_module.h"
//#include "sparse_triplets_incompressibleFlow.h"

// Global Variables
// ============================================================================

// skalar parameters
int numVars = 6;                         // Number of Variables: Velocity in x and y direction and pressure
int N;                                   // Number of Nodes
int K;                                   // Number of time-steps
int numBoundNodes;                       // Numbber of noundary Nodes
double frequency;                        // Frequency of the simulation
double massDensity;                      // Mass density
double lame1;                            // first lame constant
double lame2;                            // second lame constant
double hSquared;                         // Square of the support radius
std::size_t sizeDouble = sizeof(double); // size of double in bytes

// Eigen3 Matrix Objects
Eigen::Matrix<double, Eigen::Dynamic, 1> solverOutput;             // Saves the final values of the simulation output for the current timestep
Eigen::Matrix<double, Eigen::Dynamic, 2> nodeCoord;                // Node coordinates
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeNormals; // Node Normals

// Eigen3 Array Objects
Eigen::Array<std::vector<int>,Eigen::Dynamic,1> support;                                                    // Support for each node
Eigen::Array<int,Eigen::Dynamic,1> supportSize;                                                             // The size of the support for each node
Eigen::Array<unsigned int, Eigen::Dynamic, 1> nodeTypes;                                                    // the boundary condition type for each node
Eigen::Array<int,Eigen::Dynamic,1> bufferOffsets;                                                           // offsets bit which indicates the start of a new sparse matrix row in the buffer
Eigen::Array<int,Eigen::Dynamic,1> posDiagElements;                                                         // for each node, the index of the own position in the support is saved

// viennaCL Objects
viennacl::compressed_matrix<double> vcl_disp_A;     // displacement field sparse coefficient matrix
viennacl::compressed_matrix<double> vcl_disp_x_A;   // displacement field derivative-x sparse coefficient matrix
viennacl::compressed_matrix<double> vcl_disp_y_A;   // displacement field derivative-y sparse coefficient matrix

double scale_mass = 1;
double scale_time = 1;
double scale_length = 1e3;
double scale_stress = scale_mass/(scale_length*scale_time*scale_time);

// other object types
std::ofstream outputFile;

void computeTailorCoeff( int i, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& M, Eigen::Matrix<double, Eigen::Dynamic, 1> &W)
{
    Eigen::Matrix<double, 2, 1> diff;

    for (int j=0; j<supportSize(i); j++)
    {
        diff = nodeCoord(support(i)[j],Eigen::all) - nodeCoord(i,Eigen::all);

        M(j,0)      = 1;
        M(j,1)      = diff(0);
        M(j,2)      = diff(1);
        M(j,3)      = 0.5*diff(0)*diff(0);
        M(j,4)      = diff(0)*diff(1);
        M(j,5)      = 0.5*diff(1)*diff(1);
        W(j)        = exp(-5.25*diff.squaredNorm()/(hSquared));
    }
}

void computeStencilCoeff(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& M,
                         Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& MdotW,
                         Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& C,
                         Eigen::Matrix<double, Eigen::Dynamic, 1> W)
{
    MdotW = M.transpose()*W.asDiagonal();
    C = (MdotW*M).llt().solve(MdotW);
    //Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Wdiag = W.asDiagonal();
    //std::cout << (MdotW*M).norm()*(MdotW*M).inverse().norm() << std::endl;
    //C = (Wdiag*M).completeOrthogonalDecomposition().solve(Wdiag);
}

void setMomentumEquation(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& M, int &n)
{
    // momentum equation x-direction
    M(2*n,0)  = (massDensity*frequency*frequency); // v
    M(2*n,3)  = -(lame1+2*lame2);                  // v_xx
    M(2*n,5)  = -lame2;                            // v_yy
    M(2*n,10) = -(lame1+lame2);                    // w_xy

    // momentum equation y-direction
    M(2*n+1,4)  = -(lame1+lame2);                    // v_xy
    M(2*n+1,6)  = (massDensity*frequency*frequency); // w
    M(2*n+1,9)  = -lame2;                            // w_xx
    M(2*n+1,11) = -(lame1+2*lame2);                  // w_yy
}

void fillMat(int& i, int& n, std::vector< std::map< unsigned int, double>> &triplets, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& C)
{
    for (int j=0; j<n; j++)
    {
        if (support(i)[j]==i)
        {
            triplets[i][i]     = 1-C(0,j);
            triplets[i][i+N]   =  -C(0,j+n);

            triplets[i+N][i]   =  -C(6,j);
            triplets[i+N][i+N] = 1-C(6,j+n);
        }
        else
        {
            triplets[i][support(i)[j]]     = -C(0,j);
            triplets[i][support(i)[j]+N]   = -C(0,j+n);

            triplets[i+N][support(i)[j]]   = -C(6,j);
            triplets[i+N][support(i)[j]+N] = -C(6,j+n);
        }
    }
}

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
    Eigen::Array<double, Eigen::Dynamic, 4> disp_coeff;                  // RHS coefficients for displacement field
    Eigen::Array<double, Eigen::Dynamic, 4> disp_x_coeff;                  // RHS coefficients for displacement x-derivative field
    Eigen::Array<double, Eigen::Dynamic, 4> disp_y_coeff;                  // RHS coefficients for displacement y-derivative field

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
        frequency     = modelInput.simulationParameters.frequency/scale_time;
        massDensity   = modelInput.materialParameters.massDensity*scale_mass/(scale_length*scale_length*scale_length);
        double elasticModulus = modelInput.materialParameters.elasticModulus*scale_stress;
        double poissonRatio   = modelInput.materialParameters.poissonRatio;
        lame1 = poissonRatio/(1-2*poissonRatio) * 1/(1+poissonRatio) * elasticModulus;
        lame2 = 0.5 * 1/(1+poissonRatio) * elasticModulus;
        boundValMatX  = modelInput.boundaryConditions.values1;
        boundValMatY  = modelInput.boundaryConditions.values2;

        // Setting up solver tags
        solverTagBicgstab = viennacl::linalg::bicgstab_tag(modelInput.simulationParameters.solverTolerance);
        solverTagGmres = viennacl::linalg::gmres_tag(modelInput.simulationParameters.solverTolerance);

        // resizing of global variables
        solverOutput.resize(6*N);

        allocateBoundaryConditions(mesh, modelInput, nodeTypes, boundAllocMat); // Allocate the boundary conditions

        // load pointcloud
        mesh.nodes.coord = mesh.nodes.coord*scale_length;
        Pointcloud pointcloud;
        pointcloud.importMesh(mesh);
        support = pointcloud.support;
        nodeCoord = pointcloud.coord;
        nodeNormals = pointcloud.nodeNormals;
        hSquared = pointcloud.h*pointcloud.h;
        supportSize = pointcloud.supportSize;

        // Setting Initial Conditions
        solverOutput(Eigen::seqN(0,N))   = pointcloud.coord(Eigen::all,0);                           // x-pos
        solverOutput(Eigen::seqN(N,N))   = pointcloud.coord(Eigen::all,1);                           // y-pos
        solverOutput(Eigen::seqN(2*N,N)) = modelInput.initialConditions.initialValues(Eigen::all,2); // stress sigma_xx
        solverOutput(Eigen::seqN(3*N,N)) = modelInput.initialConditions.initialValues(Eigen::all,3); // stress sigma_xy
        solverOutput(Eigen::seqN(4*N,N)) = modelInput.initialConditions.initialValues(Eigen::all,4); // stress sigma_yy
        solverOutput(Eigen::seqN(5*N,N)) = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(N,1);;      // stress sigma_yy

        // Fill Matrix
        disp_coeff   = Eigen::Array<double, Eigen::Dynamic, 4>::Zero(2*N,4);
        disp_x_coeff = Eigen::Array<double, Eigen::Dynamic, 4>::Zero(2*N,4);
        disp_y_coeff = Eigen::Array<double, Eigen::Dynamic, 4>::Zero(2*N,4);

        std::vector< std::map< unsigned int, double>> disp_triplets(2*N);   // Displacement Triplet List For Matrix filling
        std::vector< std::map< unsigned int, double>> disp_x_triplets(2*N); // Displacement x-derivative Triplet List For Matrix filling
        std::vector< std::map< unsigned int, double>> disp_y_triplets(2*N); // Displacement y-derivative Triplet List For Matrix filling
        #pragma omp parallel
            {
                int n;
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M;
                Eigen::Matrix<double, Eigen::Dynamic, 1> W;
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MdotW;
                Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C;

                #pragma omp for nowait
                for (int i=0; i<numBoundNodes; i++) // Boundary Points
                {
                    n = support(i).size();
                    M = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(2*n+4,12);
                    W.resize(2*n+4);
                    MdotW.resize(12,2*n+4);
                    C.resize(12,2*n+4);

                    // compute tailor-expansion coefficients
                    computeTailorCoeff(i,M,W);
                    M(Eigen::seqN(n,n),Eigen::seqN(6,6)) = M(Eigen::seqN(0,n),Eigen::seqN(0,6));
                    W(Eigen::seqN(n,n)) = W(Eigen::seqN(0,n));

                    switch(nodeTypes(i))
                    {
                        case 1: // Velocity boundary condition

                            setMomentumEquation(M,n);

                            M(2*n+2,0) = 1; // velocity boundary condition x direction
                            M(2*n+3,6) = 1; // velocity boundary condition y direction

                            // Weight function
                            W(Eigen::seqN(2*n,4)) << 2, 2, 2, 2;

                            computeStencilCoeff(M,MdotW,C,W); // Compute Stencil coefficients

                            fillMat(i,n,disp_triplets,C); // Fill Triplet Lists
                            for (int j=0; j<n; j++)
                            {
                                disp_x_triplets[i][support(i)[j]]     = C(1,j);
                                disp_x_triplets[i][support(i)[j]+N]   = C(1,j+n);
                                disp_x_triplets[i+N][support(i)[j]]   = C(7,j);
                                disp_x_triplets[i+N][support(i)[j]+N] = C(7,j+n);

                                disp_y_triplets[i][support(i)[j]]     = C(2,j);
                                disp_y_triplets[i][support(i)[j]+N]   = C(2,j+n);
                                disp_y_triplets[i+N][support(i)[j]]   = C(8,j);
                                disp_y_triplets[i+N][support(i)[j]+N] = C(8,j+n);
                            }

                            // Set RHS coefficients:
                            disp_coeff(i,Eigen::seqN(0,2))   = C(0,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_coeff(i,Eigen::seqN(2,2))   = C(0,Eigen::seqN(2*n+2,2));
                            disp_coeff(i+N,Eigen::seqN(0,2)) = C(6,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_coeff(i+N,Eigen::seqN(2,2)) = C(6,Eigen::seqN(2*n+2,2));

                            disp_x_coeff(i,Eigen::seqN(0,2))   = C(1,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_x_coeff(i,Eigen::seqN(2,2))   = C(1,Eigen::seqN(2*n+2,2));
                            disp_x_coeff(i+N,Eigen::seqN(0,2)) = C(7,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_x_coeff(i+N,Eigen::seqN(2,2)) = C(7,Eigen::seqN(2*n+2,2));

                            disp_y_coeff(i,Eigen::seqN(0,2))   = C(2,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_y_coeff(i,Eigen::seqN(2,2))   = C(2,Eigen::seqN(2*n+2,2));
                            disp_y_coeff(i+N,Eigen::seqN(0,2)) = C(8,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_y_coeff(i+N,Eigen::seqN(2,2)) = C(8,Eigen::seqN(2*n+2,2));

                            break;

                        case 2: // Stress boundary condition

                            setMomentumEquation(M,n);

                            // Neumann boundary condition x-direction
                            M(2*n+2,1) = nodeNormals(i,0)*(lame1+2*lame2); // v_x
                            M(2*n+2,2) = (lame2 * nodeNormals(i,1));       // v_y
                            M(2*n+2,7) = (lame2 * nodeNormals(i,1));       // w_x
                            M(2*n+2,8) = (lame1 * nodeNormals(i,0));       // w_y

                            // Neumann boundary condition y-direction
                            M(2*n+3,1) = (lame1 * nodeNormals(i,1));       // v_x
                            M(2*n+3,2) = (lame2 * nodeNormals(i,0));       // v_y
                            M(2*n+3,7) = (lame2 * nodeNormals(i,0));       // w_x
                            M(2*n+3,8) = nodeNormals(i,1)*(lame1+2*lame2); // w_y

                            // Weight function
                            W(Eigen::seqN(2*n,4)) << 2, 2, 2, 2;

                            computeStencilCoeff(M,MdotW,C,W); // Compute Stencil coefficients

                            fillMat(i,n,disp_triplets,C); // Fill Triplet Lists
                            for (int j=0; j<n; j++)
                            {
                                disp_x_triplets[i][support(i)[j]]     = C(1,j);
                                disp_x_triplets[i][support(i)[j]+N]   = C(1,j+n);
                                disp_x_triplets[i+N][support(i)[j]]   = C(7,j);
                                disp_x_triplets[i+N][support(i)[j]+N] = C(7,j+n);

                                disp_y_triplets[i][support(i)[j]]     = C(2,j);
                                disp_y_triplets[i][support(i)[j]+N]   = C(2,j+n);
                                disp_y_triplets[i+N][support(i)[j]]   = C(8,j);
                                disp_y_triplets[i+N][support(i)[j]+N] = C(8,j+n);
                            }

                            disp_coeff(i,Eigen::seqN(0,2))   = C(0,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_coeff(i,Eigen::seqN(2,2))   = C(0,Eigen::seqN(2*n+2,2))*scale_stress;
                            disp_coeff(i+N,Eigen::seqN(0,2)) = C(6,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_coeff(i+N,Eigen::seqN(2,2)) = C(6,Eigen::seqN(2*n+2,2))*scale_stress;

                            disp_x_coeff(i,Eigen::seqN(0,2))   = C(1,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_x_coeff(i,Eigen::seqN(2,2))   = C(1,Eigen::seqN(2*n+2,2))*scale_stress;
                            disp_x_coeff(i+N,Eigen::seqN(0,2)) = C(7,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_x_coeff(i+N,Eigen::seqN(2,2)) = C(7,Eigen::seqN(2*n+2,2))*scale_stress;

                            disp_y_coeff(i,Eigen::seqN(0,2))   = C(2,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_y_coeff(i,Eigen::seqN(2,2))   = C(2,Eigen::seqN(2*n+2,2))*scale_stress;
                            disp_y_coeff(i+N,Eigen::seqN(0,2)) = C(8,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                            disp_y_coeff(i+N,Eigen::seqN(2,2)) = C(8,Eigen::seqN(2*n+2,2))*scale_stress;

                            break;
                    }
                }

                #pragma omp for nowait
                for (int i=numBoundNodes; i<N; i++) // inner nodes
                {
                    n = support(i).size();
                    M = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(2*n+2,12);
                    W.resize(2*n+2);
                    MdotW.resize(12,2*n+2);
                    C.resize(12,2*n+2);

                    // compute tailor-expansion coefficients
                    computeTailorCoeff(i,M,W);
                    M(Eigen::seqN(n,n),Eigen::seqN(6,6))    = M(Eigen::seqN(0,n),Eigen::seqN(0,6));
                    W(Eigen::seqN(n,n))   = W(Eigen::seqN(0,n));

                    setMomentumEquation(M,n);

                    // Weight function
                    W(Eigen::seqN(2*n,2)) << 2, 2;

                    computeStencilCoeff(M,MdotW,C,W); // Compute Stencil coefficients

                    fillMat(i,n,disp_triplets,C); // Fill Triplet Lists
                    for (int j=0; j<n; j++)
                    {
                        disp_x_triplets[i][support(i)[j]]     = C(1,j);
                        disp_x_triplets[i][support(i)[j]+N]   = C(1,j+n);
                        disp_x_triplets[i+N][support(i)[j]]   = C(7,j);
                        disp_x_triplets[i+N][support(i)[j]+N] = C(7,j+n);

                        disp_y_triplets[i][support(i)[j]]     = C(2,j);
                        disp_y_triplets[i][support(i)[j]+N]   = C(2,j+n);
                        disp_y_triplets[i+N][support(i)[j]]   = C(8,j);
                        disp_y_triplets[i+N][support(i)[j]+N] = C(8,j+n);
                    }

                    // Set RHS coefficients:
                    disp_coeff(i,Eigen::seqN(0,2))   = C(0,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                    disp_coeff(i+N,Eigen::seqN(0,2)) = C(6,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);

                    disp_x_coeff(i,Eigen::seqN(0,2))   = C(1,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                    disp_x_coeff(i+N,Eigen::seqN(0,2)) = C(7,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);

                    disp_y_coeff(i,Eigen::seqN(0,2))   = C(2,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                    disp_y_coeff(i+N,Eigen::seqN(0,2)) = C(8,Eigen::seqN(2*n,2))*(massDensity*frequency*frequency);
                }
            } // End of parallel Region
            viennacl::copy(disp_triplets, vcl_disp_A);
            viennacl::copy(disp_x_triplets, vcl_disp_x_A);
            viennacl::copy(disp_y_triplets, vcl_disp_y_A);
    } // End of local Scope

    // Initialise Eigen Arrays
    Eigen::Array<double, Eigen::Dynamic, 1> disp         = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(2*N);       // displacement
    Eigen::Array<double, Eigen::Dynamic, 1> disp_x       = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(2*N);      // displacement x-derivative
    Eigen::Array<double, Eigen::Dynamic, 1> disp_y       = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(2*N);      // displacement y-derivative
    Eigen::Array<double, Eigen::Dynamic, 1> disp_old     = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(2*N);       // solution of the global-equatio-system
    Eigen::Array<double, Eigen::Dynamic, 1> disp_B       = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(2*N);      // displacementright hand side of the global-equation system
    Eigen::Array<double, Eigen::Dynamic, 1> disp_x_B     = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(2*N);      // displacement x-derivative right hand side of the global-equation system
    Eigen::Array<double, Eigen::Dynamic, 1> disp_y_B     = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(2*N);      // displacement y-derivative right hand side of the global-equation system
    Eigen::Array<double, Eigen::Dynamic, 1> valuesMinMax = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(numVars*2); // [vMax, wMax, pMax, vMin, wMin, pMin]
    Eigen::Array<double, 12, 1> valuesMinMaxTemp         = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(numVars*2); // [vMax, wMax, pMax, vMin, wMin, pMin]

    // Initialise viennaCL Objects
    viennacl::vector<double> vcl_disp(2*N);
    viennacl::vector<double> vcl_disp_B(2*N);
    viennacl::vector<double> vcl_disp_x(2*N);
    viennacl::vector<double> vcl_disp_x_B(2*N);
    viennacl::vector<double> vcl_disp_y(2*N);
    viennacl::vector<double> vcl_disp_y_B(2*N);

    // Set Up Preconditioners:
    viennacl::linalg::jacobi_precond< viennacl::compressed_matrix<double>> vcl_Precond(vcl_disp_A, viennacl::linalg::jacobi_tag()); // Jacobi preconditioner for LES of predicted velocity

    // setting up solvers
    viennacl::linalg::bicgstab_solver<viennacl::vector<double>> vcl_solver(solverTagBicgstab);         // gmres solver for LES of pressure correction

    // Setting Up Output File
    frequency = frequency*scale_time;
    generateOutputFile(workingDirectory + "/" + outputName);
    int percentage = 10; // Percentage of the simulation progress

    Eigen::Array<double, Eigen::Dynamic, 1> boundValuesX = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(numBoundNodes); // saves the boundary values if the current timestep in x-direction
    Eigen::Array<double, Eigen::Dynamic, 1> boundValuesY = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(numBoundNodes); // saves the boundary values if the current timestep in y-direction

    for (int k=0; k<K; k++)
    {
        if (k>0)
        {
            // Set RHS coefficients
            boundValuesX = (boundAllocMat*boundValMatX(Eigen::all,k)).array();
            boundValuesY = (boundAllocMat*boundValMatY(Eigen::all,k)).array();

            #pragma omp parallel for
            for (int i = 0; i<numBoundNodes; i++)
            {
                disp_B(i)     = disp_coeff(i,0)     * (2*disp(i)-disp_old(i))
                              + disp_coeff(i,1)     * (2*disp(i+N)-disp_old(i+N))
                              + disp_coeff(i,2)     * boundValuesX(i)
                              + disp_coeff(i,3)     * boundValuesY(i);

                disp_B(i+N)   = disp_coeff(i+N,0)   * (2*disp(i)-disp_old(i))
                              + disp_coeff(i+N,1)   * (2*disp(i+N)-disp_old(i+N))
                              + disp_coeff(i+N,2)   * boundValuesX(i)
                              + disp_coeff(i+N,3)   * boundValuesY(i);

                disp_x_B(i)   = disp_x_coeff(i,0)   * (2*disp(i)-disp_old(i))
                              + disp_x_coeff(i,1)   * (2*disp(i+N)-disp_old(i+N))
                              + disp_x_coeff(i,2)   * boundValuesX(i)
                              + disp_x_coeff(i,3)   * boundValuesY(i);

                disp_x_B(i+N) = disp_x_coeff(i+N,0) * (2*disp(i)-disp_old(i))
                              + disp_x_coeff(i+N,1) * (2*disp(i+N)-disp_old(i+N))
                              + disp_x_coeff(i+N,2) * boundValuesX(i)
                              + disp_x_coeff(i+N,3) * boundValuesY(i);

                disp_y_B(i)   = disp_y_coeff(i,0)   * (2*disp(i)-disp_old(i))
                              + disp_y_coeff(i,1)   * (2*disp(i+N)-disp_old(i+N))
                              + disp_y_coeff(i,2)   * boundValuesX(i)
                              + disp_y_coeff(i,3)   * boundValuesY(i);

                disp_y_B(i+N) = disp_y_coeff(i+N,0) * (2*disp(i)-disp_old(i))
                              + disp_y_coeff(i+N,1) * (2*disp(i+N)-disp_old(i+N))
                              + disp_y_coeff(i+N,2) * boundValuesX(i)
                              + disp_y_coeff(i+N,3) * boundValuesY(i);
            }

            #pragma omp parallel for
            for (int i = numBoundNodes; i<N; i++)
            {
                disp_B(i)     = disp_coeff(i,0)     * (2*disp(i)-disp_old(i))
                              + disp_coeff(i,1)     * (2*disp(i+N)-disp_old(i+N));

                disp_B(i+N)   = disp_coeff(i+N,0)   * (2*disp(i)-disp_old(i))
                              + disp_coeff(i+N,1)   * (2*disp(i+N)-disp_old(i+N));

                disp_x_B(i)   = disp_x_coeff(i,0)   * (2*disp(i)-disp_old(i))
                              + disp_x_coeff(i,1)   * (2*disp(i+N)-disp_old(i+N));

                disp_x_B(i+N) = disp_x_coeff(i+N,0) * (2*disp(i)-disp_old(i))
                              + disp_x_coeff(i+N,1) * (2*disp(i+N)-disp_old(i+N));

                disp_y_B(i)   = disp_y_coeff(i,0)   * (2*disp(i)-disp_old(i))
                              + disp_y_coeff(i,1)   * (2*disp(i+N)-disp_old(i+N));

                disp_y_B(i+N) = disp_y_coeff(i+N,0) * (2*disp(i)-disp_old(i))
                              + disp_y_coeff(i+N,1) * (2*disp(i+N)-disp_old(i+N));
            }

            disp_old = disp;

            // Copy Eigen-Arrays to vcl arrays
            viennacl::copy(disp_B, vcl_disp_B);

            // Solve System
            vcl_solver.set_initial_guess(vcl_disp);
            vcl_disp = vcl_solver(vcl_disp_A,vcl_disp_B);
            vcl_disp_x = viennacl::linalg::prod(vcl_disp_x_A, vcl_disp) + vcl_disp_x_B;
            vcl_disp_y = viennacl::linalg::prod(vcl_disp_y_A, vcl_disp) + vcl_disp_y_B;

            // copy back to Eigen
            viennacl::copy(vcl_disp,disp);
            viennacl::copy(vcl_disp_x,disp_x);
            viennacl::copy(vcl_disp_y,disp_y);

            // Update node coordinates and calculate stress:
            #pragma omp parallel for
            for (int i=0; i<N; i++)
            {
                solverOutput(i)   = (nodeCoord(i,0) + disp(i))/scale_length;    // particle position x-direction
                solverOutput(i+N) = (nodeCoord(i,1) + disp(i+N))/scale_length;  // particle position y-direction

                solverOutput(i+2*N) = (disp_x(i) * (lame1+2*lame2) + 0.5 * lame1 * disp_y(i+N))/scale_stress; // sigma_xx
                solverOutput(i+3*N) = (lame2 * (disp_x(i+N) + disp_y(i)))/scale_stress;                       // sigma_xy
                solverOutput(i+4*N) = (disp_y(i+N) * (lame1+2*lame2) + 0.5 * lame1 * disp_x(i))/scale_stress; // sigma_yy
                solverOutput(i+5*N) = sqrt(solverOutput(i+2*N)*solverOutput(i+2*N) + 3 * solverOutput(i+3*N)*solverOutput(i+3*N) + solverOutput(i+4*N)*solverOutput(i+4*N)); // von-Mises stress
            }
        }
        // Calculate the min-max values for the current timestep
        #pragma omp parallel for
        for (int i=0; i<12; i++)
        {
            if (i<6)
            {
                valuesMinMaxTemp(i) = solverOutput(Eigen::seqN(i*N,N)).array().maxCoeff();

                if (valuesMinMaxTemp(i)> valuesMinMax(i))
                {
                    valuesMinMax(i) = valuesMinMaxTemp(i);
                }
            }
            else
            {
                valuesMinMaxTemp(i) = solverOutput(Eigen::seqN((i-6)*N,N)).array().minCoeff();

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
