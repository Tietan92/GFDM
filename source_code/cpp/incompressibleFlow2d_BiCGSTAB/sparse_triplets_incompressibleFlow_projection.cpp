#include "sparse_triplets_incompressibleFlow_projection.h"
#include  <iostream>        // std::cout
#include <cmath>

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

void computeTailorCoeffInnerNodes()
{
    int n;
    for (int i = 0; i<(N-numBoundNodes); i++)
    {
        n = supportSize(i+numBoundNodes);
        innerNodesTailorCoeff(i).resize(n,6);
        innerNodesWeights(i).resize(n);
        computeTailorCoeff(i+numBoundNodes,innerNodesTailorCoeff(i),innerNodesWeights(i));
    }
}

void computeStencilCoeff(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& M,
                         Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& MdotW,
                         Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& C,
                         Eigen::Matrix<double, Eigen::Dynamic, 1> W)
{
    MdotW = M.transpose()*W.asDiagonal();
    C = (MdotW*M).llt().solve(MdotW);
}


void fillMatRow(viennacl::compressed_matrix<double>& A ,int row, Eigen::Matrix<double, Eigen::Dynamic, 1> stencilCoeff)
{
    int offset = bufferOffsets(row);
    for (double element : stencilCoeff)
    {
        viennacl::backend::memory_write(A.handle(),offset,sizeDouble,&element);
        offset += sizeDouble;
    }
}

void fillMat_predVel_boundNodes()
{
    #pragma omp parallel
    {
        // Setting Up Matrices for Least Square Optimization
        int n;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M;
        Eigen::Matrix<double, Eigen::Dynamic, 1> W;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MdotW;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C;

        #pragma omp for nowait
        for (int i=0; i<numBoundNodes; i++) // Boundary Points
        {
            n = support(i).size();
            M = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(2*n+3,12);
            W.resize(2*n+3);
            MdotW.resize(12,2*n+3);
            C.resize(12,2*n+3);

            computeTailorCoeff(i,M,W);
            M(Eigen::seqN(n,n),Eigen::seqN(6,6)) = M(Eigen::seqN(0,n),Eigen::seqN(0,6));
            W(Eigen::seqN(n,n)) = W(Eigen::seqN(0,n));

            switch(nodeTypes(i))
            {
                case 1: // Dirichlet Boundary Condition

                    // Dirichlet velocity boundary condition x direction
                    M(2*n,0)      = 1;

                    // Dirichlet velocity boundary condition x direction
                    M(2*n+1,6)    = 1;

                    break;

                case 2: // Outflow Boubndary Condition

                    // Outflow Boundary Condition x direction
                    M(2*n,1)    = nodeNormals(i,0);
                    M(2*n,2)    = nodeNormals(i,1);

                    // Outflow Boundary Condition y direction
                    M(2*n+1,7)  = nodeNormals(i,0);
                    M(2*n+1,8)  = nodeNormals(i,1);

                    break;

                case 3: // Slip Boundary condition

                    // Contraint in normal direction
                    M(2*n,0) = nodeNormals(i,0);
                    M(2*n,6) = nodeNormals(i,1);

                    // Constrain in tangential direction
                    M(2*n+1,1) = nodeNormals(i,0) * (- nodeNormals(i,1));
                    M(2*n+1,2) = nodeNormals(i,1) * (- nodeNormals(i,1));
                    M(2*n+1,7) = nodeNormals(i,0) * nodeNormals(i,0);
                    M(2*n+1,8) = nodeNormals(i,1) * nodeNormals(i,0);

                    break;
            }

            // Mass balance equation
            M(2*n+2,1)  = 1;
            M(2*n+2,8)  = 1;

            // Weight function
            W(2*n)   = 2;
            W(2*n+1) = 2;
            W(2*n+2) = 2;

            // Solve the system
            computeStencilCoeff(M,MdotW,C,W);

            // Set A
            C(0,posDiagElements(i))   -= 1;
            C(6,posDiagElements(i)+n) -= 1;
            fillMatRow(vcl_A_predVel, i,   -C(0,Eigen::seqN(0,2*n)));
            fillMatRow(vcl_A_predVel, i+N, -C(6,Eigen::seqN(0,2*n)));

            fillMatRow(vcl_A_predVel_x, i,   C(1,Eigen::seqN(0,2*n)));
            fillMatRow(vcl_A_predVel_x, i+N, C(7,Eigen::seqN(0,2*n)));
            fillMatRow(vcl_A_predVel_y, i,   C(2,Eigen::seqN(0,2*n)));
            fillMatRow(vcl_A_predVel_y, i+N, C(8,Eigen::seqN(0,2*n)));

            // calculate right hand side coefficients
            velPred_boundNodes_rhsCoeffX(i,0)                 = C(0,2*n);
            velPred_boundNodes_rhsCoeffX(i+numBoundNodes,0)   = C(6,2*n);
            velPred_boundNodes_rhsCoeffY(i,0)                 = C(0,2*n+1);
            velPred_boundNodes_rhsCoeffY(i+numBoundNodes,0)   = C(6,2*n+1);

            velPred_boundNodes_rhsCoeffX(i,1)                 = C(1,2*n);
            velPred_boundNodes_rhsCoeffX(i+numBoundNodes,1)   = C(7,2*n);
            velPred_boundNodes_rhsCoeffY(i,1)                 = C(1,2*n+1);
            velPred_boundNodes_rhsCoeffY(i+numBoundNodes,1)   = C(7,2*n+1);

            velPred_boundNodes_rhsCoeffX(i,2)                 = C(2,2*n);
            velPred_boundNodes_rhsCoeffX(i+numBoundNodes,2)   = C(8,2*n);
            velPred_boundNodes_rhsCoeffY(i,2)                 = C(2,2*n+1);
            velPred_boundNodes_rhsCoeffY(i+numBoundNodes,2)   = C(8,2*n+1);
        }
    }// End of parallel Region
}

void fillMat_predVel_innerNodes(int k)
{
    #pragma omp parallel
    {
        // Setting Up Matrices for Least Square Optimization
        int n;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M;
        Eigen::Matrix<double, Eigen::Dynamic, 1> W;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MdotW;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C;

        double rhs1 = 0;
        double rhs2 = 0;

        #pragma omp for nowait
        for (int i=numBoundNodes; i<N; i++) // inner nodes
        {
            n = support(i).size();
            M = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(2*n+2,12);
            W.resize(2*n+2);
            MdotW.resize(12,2*n+2);
            C.resize(12,2*n+2);

            // Import weight and values for tailer coefficients
            M(Eigen::seqN(0,n),Eigen::seqN(0,6))    = innerNodesTailorCoeff(i-numBoundNodes);
            M(Eigen::seqN(n,n),Eigen::seqN(6,6))    = M(Eigen::seqN(0,n),Eigen::seqN(0,6));
            W(Eigen::seqN(0,n))   = innerNodesWeights(i-numBoundNodes);
            W(Eigen::seqN(n,n))   = W(Eigen::seqN(0,n));

            // Navier Stokes equation x-direction
            M(2*n,0)    = 1;
            M(2*n,1)    = solverOutput(i)/frequency;
            M(2*n,2)    = solverOutput(i+N)/frequency;
            M(2*n,3)    = (- dynViscosity)/(massDensity * frequency);
            M(2*n,5)    = (- dynViscosity)/(massDensity * frequency);

            // Navier Stokes equation y-direction
            M(2*n+1,6)  = 1;
            M(2*n+1,7)  = solverOutput(i)/frequency;
            M(2*n+1,8)  = solverOutput(i+N)/frequency;
            M(2*n+1,9)  = (- dynViscosity)/(massDensity * frequency);
            M(2*n+1,11) = (- dynViscosity)/(massDensity * frequency);

            // Weight Function
            W(2*n)   = 2;
            W(2*n+1) = 2;

            // Solve the system
            computeStencilCoeff(M,MdotW,C,W);

            // Set A
            C(0,posDiagElements(i))   -= 1;
            C(6,posDiagElements(i)+n) -= 1;
            fillMatRow(vcl_A_predVel, i,   -C(0,Eigen::seqN(0,2*n)));
            fillMatRow(vcl_A_predVel, i+N, -C(6,Eigen::seqN(0,2*n)));

            fillMatRow(vcl_A_predVel_x, i,   C(1,Eigen::seqN(0,2*n)));
            fillMatRow(vcl_A_predVel_x, i+N, C(7,Eigen::seqN(0,2*n)));
            fillMatRow(vcl_A_predVel_y, i,   C(2,Eigen::seqN(0,2*n)));
            fillMatRow(vcl_A_predVel_y, i+N, C(8,Eigen::seqN(0,2*n)));

            rhs1 = solverOutput(i)   - 1/(frequency*massDensity) * corrP_x(i);
            rhs2 = solverOutput(i+N) - 1/(frequency*massDensity) * corrP_y(i);

            B_predVel(i)      = rhs1 * C( 0,2*n) + rhs2 * C(0,2*n+1);
            B_predVel(i+N)    = rhs1 * C( 6,2*n) + rhs2 * C(6,2*n+1);
            B_predVel_x(i)    = rhs1 * C(1, 2*n) + rhs2 * C(1, 2*n+1);
            B_predVel_x(i+N)  = rhs1 * C(7, 2*n) + rhs2 * C(7, 2*n+1);
            B_predVel_y(i)    = rhs1 * C(2, 2*n) + rhs2 * C(2, 2*n+1);
            B_predVel_y(i+N)  = rhs1 * C(8, 2*n) + rhs2 * C(8, 2*n+1);
        }
    }// End of parallel Region
}

void fillMat_corrP()
{
    #pragma omp parallel
    {
        // Setting Up Matrices for Least Square Optimization
        int n;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> M;
        Eigen::Matrix<double, Eigen::Dynamic, 1> W;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MdotW;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C;

        #pragma omp for nowait
        for (int i=0; i<numBoundNodes; i++) // Boundary Points
        {
            n = support(i).size();
            M = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(n+2,6);
            W.resize(n+2);
            MdotW.resize(6,n+2);
            C.resize(6,n+2);

            computeTailorCoeff(i,M,W);

            // Pressure Poisson Equation
            M(n,0) = 0.001;
            M(n,3) = 1;
            M(n,5) = 1;

            // Homogen Neumann boundaryy condititios
            M(n+1,1) = nodeNormals(i,0);
            M(n+1,2) = nodeNormals(i,1);

            // Weight function
            W(n)   = 2;
            W(n+1) = 2;

            // Solve the system
            computeStencilCoeff(M,MdotW,C,W);

            // Fill Triplet Lists
            C(0,posDiagElements(i))   -= 1;

            for (int j=0; j<n; j++)
            {
                corrP_triplets[i][support(i)[j]]    = -C(0,j);
                corrP_x_triplets[i][support(i)[j]]  = C(1,j);
                corrP_y_triplets[i][support(i)[j]]  = C(2,j);
                corrP_xx_triplets[i][support(i)[j]] = C(3,j);
                corrP_xy_triplets[i][support(i)[j]] = C(4,j);
                corrP_yy_triplets[i][support(i)[j]] = C(5,j);
            }
        }

        #pragma omp for nowait
        for (int i=numBoundNodes; i<N; i++) // inner nodes
        {
            n = support(i).size();
            M = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(n+1,6);
            W.resize(n+1);
            MdotW.resize(6,n+1);
            C.resize(6,n+1);

            computeTailorCoeff(i,M,W);

            // pressure poison equation
            M(n,3) = 1;
            M(n,5) = 1;

            // Weight function
            W(n)   = 2;

            // Solve the system
            computeStencilCoeff(M,MdotW,C,W);

            // Fill Triplet Lists
            C(0,posDiagElements(i))   -= 1;

            for (int j=0; j<n; j++)
            {
                corrP_triplets[i][support(i)[j]]    = -C(0,j);
                corrP_x_triplets[i][support(i)[j]]  = C(1,j);
                corrP_y_triplets[i][support(i)[j]]  = C(2,j);
                corrP_xx_triplets[i][support(i)[j]] = C(3,j);
                corrP_xy_triplets[i][support(i)[j]] = C(4,j);
                corrP_yy_triplets[i][support(i)[j]] = C(5,j);
            }
            // Save Right Hand Side Stencil Coefficients
            corrP_rhsCoeff(i,Eigen::all) = C(Eigen::all,n);
        }
    } // End of parallel Region
}

