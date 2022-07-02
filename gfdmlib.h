#ifndef GFDMLIB_H_INCLUDED
#define GFDMLIB_H_INCLUDED

#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <math.h>
#include <Eigen/Core>
#include "omp.h"


using namespace std;
using namespace Eigen;


struct structPhysicalNames{
    vector<int> dim;
    vector<int> num;
    vector<string> name;
};
struct structElements{
    vector<int> tag;
    vector<vector<int>> nodelist;
};


void readMesh(string filename,
              structPhysicalNames &physicalNames,
              Matrix<double, Dynamic, 2> &nodes,
              structElements &innerElements,
              structElements &boundaryElements)
{
    string currentLine;
    int m;
    //
    ifstream reader(filename);
    while (getline(reader, currentLine))
    {
        int n = currentLine.size();
        if (currentLine.substr(0,n) == "$PhysicalNames")
        {
            getline(reader, currentLine);
            m = stoi(currentLine);
            for (int i=0; i<m; i++)
            {
                getline(reader, currentLine,' ');
                physicalNames.dim.push_back(stoi(currentLine));
                getline(reader, currentLine,' ');
                physicalNames.num.push_back(stoi(currentLine));
                getline(reader, currentLine,'\n');
                physicalNames.name.push_back(currentLine);
            }
        }
        else if (currentLine.substr(0,n) == "$Nodes")
        {
            getline(reader, currentLine);
            m = stoi(currentLine);
            nodes.resize(m,2);
            for (int i=0; i<m; i++)
            {
                getline(reader, currentLine,' ');
                getline(reader, currentLine,' ');
                nodes(i,0) = stod(currentLine);
                getline(reader, currentLine,' ');
                nodes(i,1) = stod(currentLine);
                getline(reader, currentLine,'\n');
            }

        }
        else if (currentLine.substr(0,n) == "$Elements")
        {
            getline(reader, currentLine);
            m = stoi(currentLine);
            for (int i=0; i<m; i++)
            {
                getline(reader, currentLine,' ');
                getline(reader, currentLine,' ');
                int ElementType = stoi(currentLine);
                std::getline(reader, currentLine,' ');
                std::getline(reader, currentLine,' ');
                if (ElementType == 1)
                {

                    boundaryElements.tag.push_back(stoi(currentLine));
                    getline(reader, currentLine,' ');
                    getline(reader, currentLine,' ');
                    int node1 = stoi(currentLine);
                    getline(reader, currentLine,'\n');
                    int node2 = stoi(currentLine);
                    boundaryElements.nodelist.push_back({node1-1,node2-1});
                }
                else
                {
                    innerElements.tag.push_back(stoi(currentLine));
                    getline(reader, currentLine,' ');
                    getline(reader, currentLine,' ');
                    int node1 = stoi(currentLine);
                    getline(reader, currentLine,' ');
                    int node2 = stoi(currentLine);
                    getline(reader, currentLine,'\n');
                    int node3 = stoi(currentLine);
                    innerElements.nodelist.push_back({node1-1,node2-1,node3-1});
                }
            }
        }
    }
}

void computeFaceNormals(Matrix<double, Dynamic, 2> &nodes,
                        structElements &innerElements,
                        structElements &boundaryElements,
                        vector<vector<double>> &faceNormals )
{
    int N = innerElements.tag.size();
    int M = boundaryElements.tag.size();



    for (int i=0; i<M; i++)
    {
        vector<int> boundaryNodes = boundaryElements.nodelist[i];
        for (int j=0; j<N; j++)
        {

            vector<int> Nodesj = innerElements.nodelist[j];
            bool found1 = false;
            bool found2 = false;
            int n3;
            for (int k=0; k<3; k++)
            {
                if (Nodesj[k] == boundaryNodes[0])
                {
                    found1 = true;
                }
                else if (Nodesj[k] == boundaryNodes[1])
                {
                    found2 = true;
                }
                else
                {
                    n3 = Nodesj[k];
                }
            }
            if (found1 == true && found2 == true)
            {
                Matrix<double, 2, 1> v1{nodes(boundaryNodes[0],0),nodes(boundaryNodes[0],1)};
                Matrix<double, 2, 1> v2{nodes(boundaryNodes[1],0),nodes(boundaryNodes[1],1)};
                Matrix<double, 2, 1> v3{nodes(n3,0),nodes(n3,1)};
                Matrix<double, 2, 1> v12 = v1 - v2;
                Matrix<double, 2, 1> v13 = v3 - v1;
                Matrix<double, 2, 1> v23 = v3 - v2;
                Matrix<double, 2, 1> normal{-v12[1], v12[0]};
                normal = normal/normal.norm();
                if (normal.dot(v13) > 0 && normal.dot(v23) > 0)
                {
                    normal = - normal;
                }
                faceNormals.push_back({normal[0], normal[1]});
                break;
            }
        }
    }

}

void computeNodeNormals(structElements &boundaryElements,
                        vector<vector<double>> &faceNormals,
                        Matrix<double, Dynamic, 2> &nodeNormals)
{
    int n = faceNormals.size();
    Matrix<double, Dynamic, 1> nodeNormalNorm(n);
    nodeNormals.resize(n,2);
    for (int i=0; i<n; i++)
    {
        int I1 = boundaryElements.nodelist[i][0];
        int I2 = boundaryElements.nodelist[i][1];

        nodeNormals(I1,0) += faceNormals[i][0];
        nodeNormals(I2,0) += faceNormals[i][0];

        nodeNormals(I1,1) += faceNormals[i][1];
        nodeNormals(I2,1) += faceNormals[i][1];
    }
    nodeNormalNorm = nodeNormals.rowwise().norm();
    for (int i=0; i<n; i++)
    {
        nodeNormals(i,0) = nodeNormals(i,0)/nodeNormalNorm(i);
        nodeNormals(i,1) = nodeNormals(i,1)/nodeNormalNorm(i);
    }
}

void computeSupport(Matrix<double, Dynamic, 2> &nodes,
                    double &h, vector<vector<int>> &support)
{

    int n = nodes.rows();
    Matrix<double, Dynamic, 1> dGlobal(n);
    Array<double,Dynamic,1>  diff(n);
    Array<bool,Eigen::Dynamic,1> comp(n);


    for (int i=0; i<n; i++)
    {
        Matrix<double, 2, 1> v{nodes(i,0),nodes(i,1)};
        diff = (nodes.rowwise()- v.transpose()).rowwise().norm();
        diff(i) = INFINITY;
        dGlobal(i) = diff.minCoeff();
    }

    h = dGlobal.maxCoeff()/0.45;


    for (int i=0; i<n; i++)
    {
        Matrix<double, 2, 1> v{nodes(i,0),nodes(i,1)};
        diff = (nodes.rowwise()- v.transpose()).rowwise().norm();
        diff(i) = INFINITY;
        comp = diff <= h;
        vector<int> supportI;
        for (int j=0; j<n; j++)
        {
            if (comp(j) == true)
            {
                supportI.push_back(j);
            }
        }
        support.push_back(supportI);
    }
}

void computeNodeTags(int numOfNodes,
                     structElements &boundaryElements,
                     vector<vector<int>> &nodeTags)
{
    nodeTags.resize(numOfNodes);
    int n = boundaryElements.nodelist.size();
    for (int i=0; i<n; i++)
    {
        int tag = boundaryElements.tag[i];
        int n1 = boundaryElements.nodelist[i][0];
        int n2 = boundaryElements.nodelist[i][1];


        if (std::find(nodeTags[n1].begin(), nodeTags[n1].end(), tag) == nodeTags[n1].end())
        {
            nodeTags[n1].push_back(tag);
        }
        if (std::find(nodeTags[n2].begin(), nodeTags[n2].end(), tag) == nodeTags[n2].end())
        {
            nodeTags[n2].push_back(tag);
        }
    }
}

void changeNodeTags(vector<int> boundaryNew, vector<int> boundaryOld, vector<vector<int>> &nodeTags, int N)
{
    sort(boundaryOld.begin(), boundaryOld.end());
    for (int i=0; i<N; i++)
    {
        vector<int> tagI = nodeTags[i];
        sort(tagI.begin(), tagI.end());
        if (tagI == boundaryOld)
        {
            nodeTags[i] = boundaryNew;
        }
    }
}

void setBoundaryCondition(vector<int> physicalTag,
                          int boundaryType,
                          Array<double, Dynamic, 1> values,
                          vector<vector<int>> &nodeTags,
                          Array<double, Dynamic, Dynamic> &neumanns,
                          Array<double, Dynamic, Dynamic> &dirichlets,
                          Matrix<int, Dynamic, 1> &boundaryTypes,
                          int N
                           )
{
    for (int i=0; i<N; i++)
    {
        if (nodeTags[i] == physicalTag)
        {
            if (boundaryType == 1)
            {
                dirichlets(i,all) = values;
                boundaryTypes(i) = 1;
            }
            else if (boundaryType == 2)
            {
                neumanns(i,all) = values;
                boundaryTypes(i) = 2;
            }
        }
    }
}

double weightfunction(double h, Matrix<double, 2, 1> diff)
{
    return pow(exp(-5.25*pow(diff.norm(),2)/(2*h*h)),2);
}




#endif // GFDMLIB_H_INCLUDED
