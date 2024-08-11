#include "mesh_module.h"


template<typename T>
Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> vectorToEigenArray(int M, int N, std::vector<std::vector<T>> vectorInput)
{
    Eigen::Array<T,Eigen::Dynamic,Eigen::Dynamic> arrayOutput;
    arrayOutput.resize(M,N);
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
        {
            arrayOutput(i,j) = vectorInput[i][j];
        }
    }
    return arrayOutput;
}


Mesh::Mesh(std::string filePath)
{
    std::string currentLine;
    std::ifstream reader(filePath);

    while (std::getline(reader, currentLine))
    {
        if (currentLine.substr(0,currentLine.size()) == "$PhysicalNames")
        {
            getline(reader, currentLine);
            physicalNames.numPhysicalNames = stoi(currentLine);
            physicalNames.dimension.resize(physicalNames.numPhysicalNames);
            physicalNames.physicalTag.resize(physicalNames.numPhysicalNames);
            physicalNames.name.resize(physicalNames.numPhysicalNames);
            for (int i=0; i<physicalNames.numPhysicalNames; i++)
            {
                getline(reader, currentLine,' ');
                physicalNames.dimension(i) = stoi(currentLine);
                getline(reader, currentLine,' ');
                physicalNames.physicalTag(i) = stoi(currentLine);
                getline(reader, currentLine,'\n');
                physicalNames.name(i) = currentLine;
            }
        }
        else if (currentLine.substr(0,currentLine.size()) == "$Nodes")
        {
            getline(reader, currentLine);
            nodes.numNodes = stoi(currentLine);
            nodes.coord.resize(nodes.numNodes,3);
            for (int i=0; i<nodes.numNodes; i++)
            {
                getline(reader, currentLine,' ');
                getline(reader, currentLine,' ');
                nodes.coord(i,0) = stod(currentLine);
                getline(reader, currentLine,' ');
                nodes.coord(i,1) = stod(currentLine);
                getline(reader, currentLine);
                nodes.coord(i,2) = stod(currentLine);
            }
        }
        else if (currentLine.substr(0,currentLine.size()) == "$Elements")
        {
            std::vector<std::vector<int>> boundaryElementsElementTag;
            std::vector<std::vector<int>> innerElementsElementTag;
            std::vector<std::vector<int>> boundaryElementsNodeTag;
            std::vector<std::vector<int>> innerElementsNodeTag;

            getline(reader, currentLine);
            int m = stoi(currentLine);
            for (int i=0; i<m; i++)
            {
                getline(reader, currentLine,' ');
                getline(reader, currentLine,' ');
                int ElementType = stoi(currentLine);
                std::getline(reader, currentLine,' ');
                std::getline(reader, currentLine,' ');
                if (ElementType == 1)
                {
                    boundaryElements.numElements += 1;
                    std::vector<int> currentValue = {stoi(currentLine)};
                    boundaryElementsElementTag.push_back(currentValue);
                    getline(reader, currentLine,' ');
                    getline(reader, currentLine,' ');
                    int node1 = stoi(currentLine);
                    getline(reader, currentLine,'\n');
                    int node2 = stoi(currentLine);
                    boundaryElementsNodeTag.push_back({node1-1,node2-1});
                }
                else
                {
                    innerElements.numElements += 1;
                    innerElementsElementTag.push_back({stoi(currentLine)});
                    getline(reader, currentLine,' ');
                    getline(reader, currentLine,' ');
                    int node1 = stoi(currentLine);
                    getline(reader, currentLine,' ');
                    int node2 = stoi(currentLine);
                    getline(reader, currentLine,'\n');
                    int node3 = stoi(currentLine);
                    innerElementsNodeTag.push_back({node1-1,node2-1,node3-1});
                }
            }
            boundaryElements.elementTag = vectorToEigenArray<int>(boundaryElements.numElements,1,boundaryElementsElementTag);
            boundaryElements.nodeTag = vectorToEigenArray<int>(boundaryElements.numElements,2,boundaryElementsNodeTag);
            innerElements.elementTag = vectorToEigenArray<int>(innerElements.numElements,1,innerElementsElementTag);
            innerElements.nodeTag = vectorToEigenArray<int>(innerElements.numElements,3,innerElementsNodeTag);
        }
    }
}

