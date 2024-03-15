#include "input_reader_module_linearHeatConduction.h"


template<typename T>
Eigen::Array<T, Eigen::Dynamic, 1> readArray(std::string& currentLine,std::ifstream& reader, std::string endLine, int arraySize)
{
    int index = 0;
    Eigen::Array<T, Eigen::Dynamic, 1> arrayOutput(arraySize);
    while (currentLine != endLine)
    {
        int stringLength = currentLine.length();
        std::stringstream sstr(currentLine.substr(0,stringLength-2));
        std::string currentString;
        while (getline(sstr, currentString,','))
        {
            arrayOutput(index) = stod(currentString);
            index += 1;
        }
        getline(reader, currentLine);
    }
    return arrayOutput;
}

template<typename T>
void readParameter(std::string& currentLine, std::ifstream& reader, std::string strParameter, T& outputParameter)
{
   if (currentLine == strParameter)
   {
       getline(reader, currentLine);
       std::stringstream ss;
       ss << currentLine;
       ss >> outputParameter;
   }
}

structModelInput importParameterSet(std::string filePath)
{
    structModelInput modelInput;

    std::string currentLine;
    std::ifstream reader(filePath);
    while (getline(reader, currentLine))
    {
        if (currentLine == "$ImportMesh") // Begin of Section Import Mesh
        {
            while(currentLine !="$ImportMeshEnd")
            {
                readParameter<std::string>(currentLine, reader,"<loadMesh>", modelInput.meshFile);
                getline(reader, currentLine);
            }
        }
        else if (currentLine == "$SimulationParameters") // Begin of Section Simulation Parameters
        {
            while(currentLine !="$SimulationParametersEnd")
            {
                readParameter<double>(currentLine, reader, "<frequency>",     modelInput.simulationParameters.frequency);
                readParameter<double>(currentLine, reader, "<timeEnd>"  ,     modelInput.simulationParameters.timeEnd);
                readParameter<int   >(currentLine, reader, "<numThreads>",    modelInput.simulationParameters.numThreads);
                readParameter<int   >(currentLine, reader, "<numNodes>",      modelInput.simulationParameters.numNodes);
                readParameter<int   >(currentLine, reader, "<numTimesteps>",  modelInput.simulationParameters.numTimesteps);
                readParameter<int   >(currentLine, reader, "<numBoundaries>", modelInput.simulationParameters.numBoundaries);

                getline(reader, currentLine);
            }
        }
        else if (currentLine == "$MaterialParameters") //  Begin of Section Physical Parameter
        {
            while(currentLine !="$MaterialParametersEnd")
            {
                readParameter<double>(currentLine, reader, "<massDensity>",         modelInput.materialParameters.massDensity);
                readParameter<double>(currentLine, reader, "<specHeatCapacity>",    modelInput.materialParameters.specHeatCapacity);
                readParameter<double>(currentLine, reader, "<thermalConductivity>", modelInput.materialParameters.thermalConductivity);

                getline(reader, currentLine);
            }
        }
        else if (currentLine == "$InitialConditions") // Begin of Section initial conditions
        {
            modelInput.initialConditions.startTemperature.resize(modelInput.simulationParameters.numNodes,1);
            getline(reader, currentLine); // start temperatue identifier
            modelInput.initialConditions.startTemperature = readArray<double>(currentLine, reader, "$InitialConditionsEnd", modelInput.simulationParameters.numNodes);
        }
        else if (currentLine == "$BoundaryConditions") // Begin of Section boundary conditions
        {
            modelInput.boundaryConditions.boundaryConditionPhysicalTag.resize(modelInput.simulationParameters.numBoundaries,1);
            modelInput.boundaryConditions.boundaryConditionType.resize(modelInput.simulationParameters.numBoundaries,1);
            modelInput.boundaryConditions.filmCoefficient = Eigen::Array<double, Eigen::Dynamic, 1>::Zero(modelInput.simulationParameters.numBoundaries,1);
            modelInput.boundaryConditions.heatflux.resize(modelInput.simulationParameters.numBoundaries,modelInput.simulationParameters.numTimesteps);
            modelInput.boundaryConditions.temperature.resize(modelInput.simulationParameters.numBoundaries,modelInput.simulationParameters.numTimesteps);

            int index = 0; // index is used to write the values at the correct array position

            while(currentLine !="$BoundaryConditionsEnd")
            {
                if (currentLine == "<defineBoundaryCondition>")
                {
                    getline(reader, currentLine); // physical Tag identifier
                    getline(reader, currentLine);  // physical Tag value
                    modelInput.boundaryConditions.boundaryConditionPhysicalTag(index) = stoi(currentLine);

                    getline(reader, currentLine); // boundary condition Type identifier
                    getline(reader, currentLine); // boundary condition Type value
                    modelInput.boundaryConditions.boundaryConditionType(index) = stoi(currentLine);

                    getline(reader, currentLine); // temperature identifier
                    getline(reader, currentLine); // temperature values first row
                    if (modelInput.boundaryConditions.boundaryConditionType(index) == 2) // Neumann Boundary Conditions: Temperature is not the last parameter
                    {
                        modelInput.boundaryConditions.temperature(index,Eigen::all) = readArray<double>(currentLine, reader, "<filmCoefficient>",modelInput.simulationParameters.numTimesteps);
                    }
                    else // Dirichlet Boundary Conditions: Temperature is the last parameter
                    {
                        modelInput.boundaryConditions.temperature(index,Eigen::all) = readArray<double>(currentLine, reader, "<defineBoundaryConditionEnd>",modelInput.simulationParameters.numTimesteps);

                    }

                    if (modelInput.boundaryConditions.boundaryConditionType(index) == 2) // only for Neumann boundary conditions: film coefficients and heat flux as to be imported
                    {
                        getline(reader, currentLine); // film Coefficient value
                        modelInput.boundaryConditions.filmCoefficient(index) = stoi(currentLine);

                        getline(reader, currentLine); // heat flux identifier
                        getline(reader, currentLine); // heat flux values first row
                        modelInput.boundaryConditions.heatflux(index,Eigen::all) = readArray<double>(currentLine, reader, "<defineBoundaryConditionEnd>",modelInput.simulationParameters.numTimesteps);
                    }
                    index +=1;
                }
                getline(reader, currentLine);
            }
        }
    }
    return modelInput;
}
