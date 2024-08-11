#ifndef OUTPUT_WRITER_MODULE_H
#define OUTPUT_WRITER_MODULE_H

#include "Eigen/Core" // Eigen Core Module for Eigen::Matrix and Eigen::Array Objects
#include <iomanip>    // needed for setting decimal precision
#include <omp.h>      // open-mp parallel computing
#include <fstream>    // input output stream to operate on files
#include <cstdio>     // needed for c printf function

extern int numVars;
extern int N;
extern int K;
extern Eigen::Matrix<double, Eigen::Dynamic, 1> solverOutput;
extern double frequency;
extern std::ofstream outputFile;
extern std::size_t sizeDouble;

void generateOutputFile(std::string outputPath);

void writeToOutputFile(int &k, int &percentage);

void closeOutputFile(Eigen::Array<double, Eigen::Dynamic, 1> &valuesMinMax);

#endif
