#include "sparse_triplets_linearHeatConduction.h"


std::vector<Tri> calcTripletList(structModelInput &modelInput,Pointcloud &pointcloud, Eigen::Array<unsigned int, Eigen::Dynamic, 1> &nodeTypes)
{
    std::vector<Tri> tripletList;

    // resize B1 and B2
    B1 = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(N,1);
    B2 = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(numBoundNodes,1);

    // read physical parameters
    double massDensity         = modelInput.materialParameters.massDensity;             // mass density of the material
    double specHeatCapacity    = modelInput.materialParameters.specHeatCapacity;        // specific heat capacity of the material
    double thermalConductivity = modelInput.materialParameters.thermalConductivity;     // thermal conductivity of the material
    double thermalDiffusivity  = thermalConductivity/(massDensity * specHeatCapacity); // thermal diffussivity

    // read pointcloud arrays
    Eigen::Matrix<double, Eigen::Dynamic, 2> nodeCoord = pointcloud.coord;                      // Node coordinates
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeNormals = pointcloud.nodeNormals; // Node Normals
    Eigen::Array<std::vector<int>,Eigen::Dynamic,1> support = pointcloud.support;               // Support for each node
    double hSquared = pointcloud.h*pointcloud.h;

    // Calculate averaged film coefficients
    Eigen::Matrix<double,Eigen::Dynamic,1> filmCoefficientsAvg(numBoundNodes);
    filmCoefficientsAvg = boundAllocMat*modelInput.boundaryConditions.filmCoefficient.matrix();

    // Gomega is constant and therefore can be calculated before the parallel region
    Eigen::Matrix<double, 1, 6> Gomega = Eigen::Matrix<double, 1, 6>::Zero(1,6);
    Gomega(0,0) = frequency;
    Gomega(0,3) = - thermalDiffusivity;
    Gomega(0,5) = - thermalDiffusivity;

    #pragma omp parallel
    {
        // Setting Up Matrices for Least Square Optimization
        int n;
        Eigen::Matrix<double, 1, 6> GomegaN = Eigen::Matrix<double, 1, 6>::Zero(1,6);
        Eigen::Matrix<double, Eigen::Dynamic, 6> D;
        Eigen::Matrix<double, Eigen::Dynamic, 1> Wd;
        Eigen::Matrix<double, Eigen::Dynamic, 1> Womega;
        Eigen::Matrix<double, Eigen::Dynamic, 1> WomegaN;
        Eigen::Matrix<double, 6, Eigen::Dynamic> C;
        Eigen::Matrix<double, Eigen::Dynamic, 6> Momega;
        Eigen::Matrix<double, Eigen::Dynamic, 6> MomegaN;
        Eigen::Matrix<double, 6, Eigen::Dynamic> MdotW;

        std::vector<Tri> tripletListPrivate;

        #pragma omp for nowait
        for (int i=0; i<N; i++)
        {

            if (nodeTypes(i) == 1)
            {
                tripletListPrivate.push_back(Tri(i,i,1));
                B2(i) = 1;
            }
            else
            {
                n = support(i).size();
                D.resize(n,6);
                Wd.resize(n);
                for (int j=0; j<n; j++)
                {
                    Eigen::Matrix<double, 2, 1> diff = nodeCoord(support(i)[j],Eigen::all) - nodeCoord(i,Eigen::all);

                    D(j,0) = 1;
                    D(j,1) = diff(0);
                    D(j,2) = diff(1);
                    D(j,3) = 0.5*diff(0)*diff(0);
                    D(j,4) = diff(0)*diff(1);
                    D(j,5) = 0.5*diff(1)*diff(1);
                    Wd(j)  = exp(-5.25*diff.squaredNorm()/(hSquared));
                }
                Womega.resize(n+1);
                Womega(Eigen::seq(0,n-1)) = Wd;
                Womega(n) = 1;

                Momega.resize(n+1,6);
                Momega(Eigen::seq(0,n-1),Eigen::all) = D;
                Momega(n,Eigen::all) = Gomega;

                if (nodeTypes(i)== 0)
                {
                    MdotW.resize(6,n+1);
                    MdotW = Momega.transpose()*Womega.asDiagonal();

                    C.resize(6,n+1);
                    C = (MdotW*Momega).llt().solve(MdotW);

                    B1(i) = C(0,n) * frequency;
                }
                else if (nodeTypes(i) == 2)
                {
                    WomegaN.resize(n+2);
                    WomegaN(Eigen::seq(0,n)) = Womega;
                    WomegaN(n+1) = 1;

                    GomegaN(0,0) = filmCoefficientsAvg(i);
                    GomegaN(0,1) = thermalConductivity * nodeNormals(i,0);
                    GomegaN(0,2) = thermalConductivity * nodeNormals(i,1);

                    MomegaN.resize(n+2,6);
                    MomegaN(Eigen::seq(0,n),Eigen::all) = Momega;
                    MomegaN(n+1,Eigen::all) = GomegaN;

                    MdotW.resize(6,n+2);
                    MdotW = MomegaN.transpose()*WomegaN.asDiagonal();

                    C.resize(6,n+2);
                    C = (MdotW*MomegaN).llt().solve(MdotW);

                    B1(i) = C(0,n) * frequency;
                    B2(i) = C(0,n+1);
                }
                tripletListPrivate.push_back(Tri(i,i,(1-C(0,0))));
                for (int j=1; j<n; j++)
                {
                    tripletListPrivate.push_back(Tri(i, support(i)[j], -C(0,j)));
                }
            }
        }
        #pragma omp critical
        tripletList.insert(tripletList.end(), tripletListPrivate.begin(), tripletListPrivate.end());
    }
    return tripletList;
}
