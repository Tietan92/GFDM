#include <iostream>
//#include "mshreader.h"
#include <vector>
#include <string>
#include <iomanip>      // std::setprecision
#include "linearHeatConduction.h"
#include <cmath>
#include "omp.h"
#include <Eigen/Core>
#include "Eigen/IterativeLinearSolvers"



using namespace std;


double f = 100;
double tend = 10;
double T0 = 20;
int N;






int main()
{
    model model1("viereck.msh");
    model1.f = f;
    model1.tend = tend;
    model1.T0 = T0;
    model1.condy = 200.0;
    model1.rho = 6450.0;
    model1.heatcap = 800.0;
    model1.preprocessing();




    int K = model1.K;


    Array<double, Dynamic, 1> Taussen = 90*Array<double, Dynamic, 1>::Ones(K);


    model1.defineBoundaryCondition(vector<int> {5},1,Taussen);


    model1.solve();

    cout << model1.Y(33333,all) << endl;

    return 0;
}
