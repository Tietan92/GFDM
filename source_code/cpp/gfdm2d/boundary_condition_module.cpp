#include "boundary_condition_module.h"

void allocateBoundaryConditions(Mesh &mesh, structModelInput &modelInput, Eigen::Array<unsigned int, Eigen::Dynamic, 1> &nodeTypes)
{
    nodeTypes.resize(N,1);
    nodeTypes = nodeTypes * 0;

    boundAllocMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(numBoundNodes,modelInput.simulationParameters.numBoundaries);
    #pragma omp parallel for
    for (int i=0; i<numBoundNodes; i++)
    {
        for (int j=0; j<modelInput.simulationParameters.numBoundaries; j++)
        {
            if (modelInput.boundaryConditions.boundaryConditionPhysicalTag(j) == mesh.boundaryElements.elementTag(i))
            {
                boundAllocMat(mesh.boundaryElements.nodeTag(i,0),j) = 1;
                boundAllocMat(mesh.boundaryElements.nodeTag(i,1),j) = 1;
                break;
            }
        }
    }
    // If a node has two Boundaries: Boundary wheights have to be calculated
    #pragma omp parallel for
    for (int i = 0; i < numBoundNodes; i++)
    {
        std::vector<int> indices;
        for (int j=0; j<modelInput.simulationParameters.numBoundaries; j++)
        {
            if (boundAllocMat(i,j) == 1)
            {
                indices.push_back(j);
            }
        }
        if (indices.size() == 1)
        {
            nodeTypes(i) = modelInput.boundaryConditions.boundaryConditionType(indices[0]);
        }
        else if (indices.size() == 2)
        {
            int boundType1 = modelInput.boundaryConditions.boundaryConditionType(indices[0]);
            int boundType2 = modelInput.boundaryConditions.boundaryConditionType(indices[1]);

            if (boundType1 == boundType2)
            {
                boundAllocMat(i,indices[0]) = 0.5;
                boundAllocMat(i,indices[1]) = 0.5;
                nodeTypes(i) = boundType1;
            }
            else if (boundType1 == 1 && boundType2 == 2)
            {
                boundAllocMat(i,indices[1]) = 0;
                nodeTypes(i) = 1;
            }
            else if (boundType1 == 2 && boundType2 == 1)
            {
                boundAllocMat(i,indices[0]) = 0;
                nodeTypes(i) = 1;
            }
        }
    }
}
