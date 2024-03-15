#include "output_writer_module.h"


void generateOutputFile(std::string outputPath)
{
    outputFile.open(outputPath);
    outputFile << "<$SimulationOutput>" << std::endl;
    printf("0%%\n");
}

void writeToOutputFile(int &k, int &outputFilePos, int &percentage)
{
    #pragma omp parallel
    {
        std::stringstream double2stringConverter;
        #pragma omp for nowait
        for (int i=0; i<N; i++)
        {
            double2stringConverter << std::setprecision(11) << solverOutput(i) << ", ";
            solverOutputString(i) = double2stringConverter.str();
            double2stringConverter.str(std::string());
        }
    }

    std::stringstream textInput_k;
    textInput_k << "<t = " << k/frequency  << ">"<< std::endl;
    outputFilePos = textInput_k.tellp();
    for (int i=0; i<N; i++)
    {
        textInput_k << solverOutputString(i);
        if (textInput_k.tellp() >= outputFilePos+200)
        {
            textInput_k << std::endl;
            outputFilePos += 200;
        }
    }
    textInput_k << std::endl << "<timestepEnd>" << std::endl;
    outputFile << textInput_k.str();

    if ((k*100)/(K-1)==percentage)
    {
        printf("%d", percentage);
        printf("%%\n");
        percentage += 10;
    }
}

void closeOutputFile()
{
    outputFile << "<$SimulationOutputEnd>" << std::endl;
    outputFile.close();
}

