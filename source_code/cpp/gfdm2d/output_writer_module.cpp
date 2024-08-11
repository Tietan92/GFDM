#include "output_writer_module.h"


void generateOutputFile(std::string outputPath)
{
    outputFile.open(outputPath, std::ios::binary); // open the file in binary mode
    double doubleNumVars = static_cast<double>(numVars);
    double doubleN = static_cast<double>(N);
    double doubleK = static_cast<double>(K);
    outputFile.write(reinterpret_cast<char*>(&doubleNumVars),sizeDouble);
    outputFile.write(reinterpret_cast<char*>(&doubleN),sizeDouble);
    outputFile.write(reinterpret_cast<char*>(&doubleK),sizeDouble);
    outputFile.write(reinterpret_cast<char*>(&frequency),sizeDouble);
    printf("0%%\n");
}

void writeToOutputFile(int &k, int &percentage)
{
    for (int i=0; i<N*numVars; i++)
    {
        outputFile.write(reinterpret_cast<char*>(&solverOutput(i)),sizeDouble);
    }

    if ((k*100)/(K-1) == percentage)
    {
        printf("%d", percentage);
        printf("%%\n");
        percentage += 10;
    }
}


void closeOutputFile(Eigen::Array<double, Eigen::Dynamic, 1> &valuesMinMax)
{
    for (int i=0; i<2*numVars; i++)
    {
        outputFile.write(reinterpret_cast<char*>(&valuesMinMax(i)),sizeDouble);
    }
    outputFile.close();
}

