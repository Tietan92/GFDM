#ifndef LINEARHEATCONDUCTION_H_INCLUDED
#define LINEARHEATCONDUCTION_H_INCLUDED
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "gfdmlib.h"
#include <cmath>
#include "omp.h"
#include <Eigen/Core>


class model
{
public:
    structPhysicalNames physicalNames;
    Matrix<double, Dynamic, 2> nodes;
    structElements innerElements;
    structElements boundaryElements;
    vector<vector<double>> faceNormals;
    Matrix<double, Dynamic, 2> nodeNormals;
    vector<vector<int>> support;
    vector<vector<int>> nodeTags;
    double h;
    Array<double, Dynamic, Dynamic> neumanns;
    Array<double, Dynamic, Dynamic> dirichlets;
    Matrix<int, Dynamic, 1> boundaryTypes;
    Matrix<double, Dynamic, Dynamic> Y;

    double f;
    double tend;
    double T0;
    double condy;
    double rho;
    double heatcap;
    int K;
    int N;


    model(string);
    void preprocessing();
    void setBoundary(vector<int>,vector<int>);
    void defineBoundaryCondition(vector<int>, int , Matrix<double, Dynamic, 1>);
    void solve();

};

model::model(string filename)
{
    readMesh(filename,physicalNames,nodes,innerElements,boundaryElements);
}

void model::preprocessing()
{
    K = round(f*tend);
    N = nodes.rows();
    neumanns   = Array<double, Dynamic, Dynamic>::Zero(N,K);
    dirichlets = Array<double, Dynamic, Dynamic>::Zero(N,K);
    Y.resize(N,K);
    dirichlets(all,0) = T0*Matrix<double, Dynamic, 1>::Ones(N);
    boundaryTypes.resize(N);
    computeFaceNormals(nodes,innerElements,boundaryElements,faceNormals);
    computeNodeNormals(boundaryElements, faceNormals,nodeNormals);
    computeSupport(nodes, h, support);
    computeNodeTags(N,boundaryElements,nodeTags);
}

void model::setBoundary(vector<int> boundaryNew, vector<int> boundaryOld)
{
    changeNodeTags(boundaryNew,boundaryOld,nodeTags,N);
}

void model::defineBoundaryCondition(vector<int> physicalTag, int boundaryType, Matrix<double, Dynamic, 1> values)
{
    setBoundaryCondition(physicalTag, boundaryType, values, nodeTags, neumanns, dirichlets, boundaryTypes, N);
}

void model::solve()
{
    omp_set_dynamic(0);   // fixed number of threads
    omp_set_num_threads(6);  // assign number of threads



    typedef Triplet<double> Tri;

    vector<Tri> tripletList;

    Array<double, Dynamic, 1> beta  = Array<double, Dynamic, 1>::Zero(N,1);
    Array<double, Dynamic, 1> gamma = Array<double, Dynamic, 1>::Zero(N,1);
    Array<double, Dynamic, 1> zeta  = Array<double, Dynamic, 1>::Zero(N,1);

    Matrix<double, Dynamic, 6> D;
    Matrix<double, Dynamic, Dynamic> Wd;
    Matrix<double, Dynamic, Dynamic> Womega;
    Matrix<double, Dynamic, Dynamic> WomegaN;
    Matrix<double, 6, Dynamic> c;
    Matrix<double, Dynamic, 6> Momega;
    Matrix<double, Dynamic, 6> MomegaN;


    double a = condy/(rho*heatcap);

    Matrix<double, 1, 6> Gomega = Matrix<double, 1, 6>::Zero(1,6);
    Gomega(0,0) = f;
    Gomega(0,3) = -a;
    Gomega(0,5) = -a;

    Matrix<double, 1, 6> GomegaN = Matrix<double, 1, 6>::Zero(1,6);

    SparseMatrix<double,RowMajor> A(N,N);
    Matrix<double, Dynamic, 1> B = Matrix<double, Dynamic, 1>::Zero(N,1);
    BiCGSTAB<SparseMatrix<double, RowMajor>> solver;
    solver.setTolerance(1e-09);


    for (int i=0; i<N; i++)
    {
        if (boundaryTypes(i) == 1)
        {
            tripletList.push_back(Tri(i,i,1));
            zeta(i) = 1;
        }
        else
        {
            vector<int> sup = support[i];
            int n = sup.size();
            D.resize(n,6);
            Wd = Matrix<double, Dynamic, Dynamic>::Zero(n,n);

            for (int j=0; j<n; j++)
            {
                int I = sup[j];
                Matrix<double, 2, 1> diff = nodes(I,all) - nodes(i,all);

                D(j,0) = 1;
                D(j,1) = diff(0);
                D(j,2) = diff(1);
                D(j,3) = 0.5*diff(0)*diff(0);
                D(j,4) = diff(0)*diff(1);
                D(j,5) = 0.5*diff(1)*diff(1);
                Wd(j,j) = weightfunction(h,diff);
            }

            Womega =  Matrix<double, Dynamic, Dynamic>::Zero(n+1,n+1);
            Womega(seq(0,n-1),seq(0,n-1)) = Wd;
            Womega(n,n) = 1;

            Momega.resize(n+1,6);
            Momega(seq(0,n-1),all) = D;
            Momega(n,all) = Gomega;

            if (boundaryTypes(i)== 0)

            {   c.resize(6,n+1);
                c = (Momega.transpose()*Womega*Momega).inverse()*Momega.transpose()*Womega;
                beta(i) = c(0,n);
            }
            else if (boundaryTypes(i) == 2)
            {
                WomegaN =  Matrix<double, Dynamic, Dynamic>::Zero(n+2,n+2);
                WomegaN(seq(0,n),seq(0,n)) = Womega;
                WomegaN(n+1,n+1) = 1;

                GomegaN(0,1) = -condy*nodeNormals(i,0);
                GomegaN(0,2) = -condy*nodeNormals(i,1);

                MomegaN.resize(n+2,6);
                MomegaN(seq(0,n),all) = Momega;
                MomegaN(n+1,all) = GomegaN;

                c.resize(6,n+2);
                c = (MomegaN.transpose()*WomegaN*MomegaN).inverse()*MomegaN.transpose()*WomegaN;
                beta(i) = c(0,n);
                gamma(i) = c(0,n+1);
            }

            tripletList.push_back(Tri(i,i,(1-c(0,0))));
            for (int j=1; j<n; j++)
            {
                int I = sup[j];
                tripletList.push_back(Tri(i,I,-c(0,j)));
            }
        }
    }

    A.setFromTriplets(tripletList.begin(), tripletList.end());
    solver.compute(A);
    Y(all,0) = dirichlets(all,0);

    for (int k=1; k<K; k++)
    {
        B = (beta * f * Y(all,k-1).array() + gamma*neumanns(all,k) + zeta*dirichlets(all,k)).matrix();
        Y(all,k) = solver.solveWithGuess(B,Y(all,k-1));
    }



}
#endif // LINEARHEATCONDUCTION_H_INCLUDED
