#include "pointcloud_module.h"


extern int N;
extern int numBoundNodes;

void printVector(std::vector<int> input)
{
    for (auto& str : input)
    {
        std::cout << str << ' ';
    }
    std::cout << '\n';
    std::cout << '\n';
}


void vectorUnion(std::vector<int> &v1, std::vector<int> &v2){
    //Make the result initially equal to v1
    std::vector<int> result = v1;

    //Remove the elements in v2 from the result
    for(const int &element: v2){
        const auto it = std::find(result.begin(), result.end(), element);
        if(it != result.end()){
            result.erase(it);
        }
    }
    v2.insert(v2.end(),result.begin(),result.end());

}


Pointcloud::Pointcloud(Mesh mesh)
{
    numInnerElements = mesh.innerElements.numElements;
    boundaryTriangles.resize(numBoundNodes,3);
    coord = mesh.nodes.coord(Eigen::all,Eigen::seq(0,1));
    innerElements = mesh.innerElements.nodeTag;
    faceNormals.resize(numBoundNodes,2);
    nodeNormals.resize(numBoundNodes,2);
    support.resize(N,1);
    supportSize.resize(N,1);

    #pragma omp parallel for
    for (int i=0; i<numBoundNodes; i++)
    {
        for (int j=0; j<mesh.innerElements.numElements; j++)
        {
            bool found1 = false;
            bool found2 = false;
            int n3 = 0;
            for (int k=0; k<3; k++)
            {
                if (mesh.innerElements.nodeTag(j,k) == mesh.boundaryElements.nodeTag(i,0))
                {
                    found1 = true;
                }
                else if (mesh.innerElements.nodeTag(j,k) == mesh.boundaryElements.nodeTag(i,1))
                {
                    found2 = true;
                }
                else
                {
                    n3 = mesh.innerElements.nodeTag(j,k);
                }
            }
            if (found1 == true && found2 == true)
            {
                boundaryTriangles(i,0) = mesh.boundaryElements.nodeTag(i,0);
                boundaryTriangles(i,1) = mesh.boundaryElements.nodeTag(i,1);
                boundaryTriangles(i,2) = n3;
                break;
            }
        }
    }
    computeFaceNormals();
    computeNodeNormals();
    computeSupport();
}

void Pointcloud::computeFaceNormals()
{
    Eigen::Matrix<double, 2, 1> v1;
    Eigen::Matrix<double, 2, 1> v2;
    Eigen::Matrix<double, 2, 1> v3;
    Eigen::Matrix<double, 2, 1> v12;
    Eigen::Matrix<double, 2, 1> v13;
    Eigen::Matrix<double, 2, 1> v23;
    Eigen::Matrix<double, 2, 1> faceNormal;

    for (int i=0; i<numBoundNodes; i++)
    {
        v1 = {coord(boundaryTriangles(i,0),0),coord(boundaryTriangles(i,0),1)};
        v2 = {coord(boundaryTriangles(i,1),0),coord(boundaryTriangles(i,1),1)};
        v3 = {coord(boundaryTriangles(i,2),0),coord(boundaryTriangles(i,2),1)};
        v12 = v1 - v2;
        v13 = v3 - v1;
        v23 = v3 - v2;
        faceNormal = {-v12[1], v12[0]};
        faceNormal = faceNormal/faceNormal.norm();
        if (faceNormal.dot(v13) > 0 && faceNormal.dot(v23) > 0)
        {
            faceNormal = - faceNormal;
        }
        faceNormals(i,Eigen::all) = faceNormal;
    }
}

void Pointcloud::computeNodeNormals()
{
    Eigen::Matrix<double, Eigen::Dynamic, 1> nodeNormalNorm(numBoundNodes);
    for (int i=0; i<numBoundNodes; i++)
    {
        int I1 = boundaryTriangles(i,0);
        int I2 = boundaryTriangles(i,1);

        nodeNormals(I1,0) += faceNormals(i,0);
        nodeNormals(I2,0) += faceNormals(i,0);

        nodeNormals(I1,1) += faceNormals(i,1);
        nodeNormals(I2,1) += faceNormals(i,1);
    }
    nodeNormalNorm = nodeNormals.rowwise().norm();
    Eigen::Array<double,Eigen::Dynamic,2> tempArray(numBoundNodes,2);
    tempArray(Eigen::all,0) = nodeNormals(Eigen::all,0).array()/nodeNormalNorm.array();
    tempArray(Eigen::all,1) = nodeNormals(Eigen::all,1).array()/nodeNormalNorm.array();
    nodeNormals = tempArray.matrix();
}

void Pointcloud::computeSupport()
{
    Eigen::Array<std::vector<int>,Eigen::Dynamic,1> direktFamily(N);
    Eigen::Array<int,Eigen::Dynamic,1> sizesOld(N);
    Eigen::Matrix<double, Eigen::Dynamic, 1> dGlobal(N);
    double vectorNorm;

    for (int i = 0; i< N; i++)
    {
        direktFamily(i).push_back(i);
        sizesOld(i) = 1;
    }

    // Calculate Node Family
    for (int i=0; i<numInnerElements; i++)
    {
        std::vector<int> tripel{innerElements(i,0),innerElements(i,1),innerElements(i,2)};
        for (int j=0; j<3; j++)
        {
            int index = innerElements(i,j);
            vectorUnion(tripel,direktFamily(index));
        }
    }
    Eigen::Array<std::vector<int>,Eigen::Dynamic,1> supportUnlimited = direktFamily;

    // Calculation of h
    #pragma omp parallel for
    for (int i=0; i<N; i++)
    {
        Eigen::Matrix<double, 2, 1> v1{coord(i,0),coord(i,1)};
        double diffMin = INFINITY;
        for (int& element : supportUnlimited(i))
        {
            if (element == i)
            {
                continue;
            }
            Eigen::Matrix<double, 2, 1> v2{coord(element,0),coord(element,1)};
            vectorNorm = (v2-v1).norm();
            if (vectorNorm < diffMin)
            {
                diffMin = vectorNorm;
            }
        }
        dGlobal(i) = diffMin;
    }
    h = dGlobal.maxCoeff()/0.45;

    // Adding Next Rows to supportUnlimited
    for (int iter = 1; iter<4; iter++)
    {
        #pragma omp parallel for
        for (int i = 0; i< N; i++)
        {
            int n = supportUnlimited(i).size();
            for (int j=sizesOld(i); j<n; j++)
            {
                int index = supportUnlimited(i)[j];
                vectorUnion(direktFamily(index),supportUnlimited(i));
            }
            sizesOld(i) = n;
        }
    }

    // Calulating support
    #pragma omp parallel for
    for (int i = 0; i< N; i++)
    {
        Eigen::Matrix<double, 2, 1> v1{coord(i,0),coord(i,1)};
        for (int& element : supportUnlimited(i))
        {
            Eigen::Matrix<double, 2, 1> v2{coord(element,0),coord(element,1)};
            vectorNorm = (v2-v1).norm();
            if (vectorNorm <= h)
            {
                support(i).push_back(element);
            }
        }
        std::sort(support(i).begin(), support(i).end());
        supportSize(i) = support(i).size();
    }
}

