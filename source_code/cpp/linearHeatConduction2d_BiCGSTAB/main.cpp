// external headers
#include <string>          // for string objects
#include <fstream>         // input output stream to operate on files
#include "../Eigen/Core"   // Eigen Core Module for Eigen::Matrix and Eigen::Array Objects
#include "../Eigen/Sparse" // For algebra operations on sparse matrices
#include <omp.h>           // open-mp parallel computing

// gfdm headers
#include "input_reader_module_linearHeatConduction.h"
#include "mesh_module.h"
#include "pointcloud_module.h"
#include "boundary_condition_module.h"
#include "output_writer_module.h"
#include "sparse_triplets_linearHeatConduction.h"


typedef Eigen::Triplet<double> Tri;// Matrix Triplet

// Arrays for linear equation system
int N;                                                               // Number of Rows off the Output Matrix (equal to the number of Nodes)
int K;                                                               // Number if Columns if the Output Matrix (equal to the number of timesteps)
int numBoundNodes;                                                   // Number if Boundary Nodes
Eigen::Array<double, Eigen::Dynamic, 1> B1;                          // Array with coefficients for inner nodes
Eigen::Array<double, Eigen::Dynamic, 1> B2;                          // Array with Boundary Condition coefficients
Eigen::Matrix<double, Eigen::Dynamic, 1> solverOutput;               // Solution of the equation System
Eigen::SparseMatrix<double,Eigen::RowMajor> A;                       // coefficient Matrix if the linear equation system
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> boundAllocMat; // Boundary Condition Allocation Matrix
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> boundValMat;   // Boundary Condition Value Matrix

double frequency; // frequency must be defined globaly, because it is used by the output writer module

// variables for output file
std::ofstream outputFile;
Eigen::Matrix<std::string, Eigen::Dynamic, 1> solverOutputString;

int main(int argc, char *argv[])
{
    std::string name             = argv[1];               // The first console argument tells the project name
    std::string workingDirectory = argv[2];               // The second console argument tells the working directory
    std::string inputName        = name + "_input.txt";   // naming convention for the input file
    std::string outputName       = name + "_output.txt";  // naming convention for the output file

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver; // initialising the slver

    // Scope for Coefficient Allocation
    {
        structModelInput modelInput = importParameterSet(workingDirectory + "/" + inputName); // load Input file

        // read simulation parameters
        frequency = modelInput.simulationParameters.frequency;                // The simulation frequency determines the size of a time-step during the simulation
        omp_set_num_threads(modelInput.simulationParameters.numThreads);      // set the number of threads for multiphreading with open mp
        solver.setTolerance(modelInput.simulationParameters.solverTolerance); // The tolerance used as a criteria so stop the iteration process

        Mesh mesh(workingDirectory + "/" + modelInput.meshFile); // load mesh

        // resizing arrays and fill them with 0 values, if needed
        N = mesh.nodes.numNodes;
        K = modelInput.simulationParameters.numTimesteps;
        numBoundNodes = mesh.boundaryElements.numElements;
        solverOutput.resize(N);
        solverOutputString.resize(N);
        A.resize(N,N);
        boundValMat.resize(modelInput.simulationParameters.numBoundaries, modelInput.simulationParameters.numTimesteps);

        // Allocate the boundary conditionsa
        Eigen::Array<unsigned int, Eigen::Dynamic, 1> nodeTypes;
        allocateBoundaryConditions(mesh, modelInput, nodeTypes);

        for (int i=0; i<modelInput.simulationParameters.numBoundaries; i++)
        {
            if (modelInput.boundaryConditions.boundaryConditionType(i) ==1)
            {
                boundValMat(i,Eigen::all) = modelInput.boundaryConditions.temperature(i,Eigen::all);
            }
            else
            {
                boundValMat(i,Eigen::all) = modelInput.boundaryConditions.filmCoefficient(i) * modelInput.boundaryConditions.temperature(i,Eigen::all) +
                                            modelInput.boundaryConditions.heatflux(i,Eigen::all);
            }
        }

        // load pointcloud
        Pointcloud pointcloud(mesh);

        // Setting Initial Conditions
        solverOutput = modelInput.initialConditions.startTemperature;

        std::vector<Tri> tripletList = calcTripletList(modelInput,pointcloud,nodeTypes);
        A.setFromTriplets(tripletList.begin(), tripletList.end());
    }
    solver.compute(A);

    Eigen::Matrix<double, Eigen::Dynamic, 1> B = Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(N,1);

    // Setting Up Output File
    generateOutputFile(workingDirectory + "/" + outputName);

    int outputFilePos;   // the actual position in the output file
    int percentage = 10; // Percentage of the simulation progress

    for (int k=0; k<K; k++)
    {
        if (k>0)
        {
            B = B1 * solverOutput.array();
            B(Eigen::seq(0,numBoundNodes-1)) += (B2 * (boundAllocMat*boundValMat(Eigen::all,k)).array()).matrix();
            solverOutput = solver.solveWithGuess(B,solverOutput);
        }
        writeToOutputFile(k,outputFilePos,percentage);
    }
    closeOutputFile();

    return 0;
}
