#ifndef ReactionDiffusionSimulation_h
#define ReactionDiffusionSimulation_h 1
#include <vector>
#include <map>
#include <utility>  
#include <string>
using namespace std;


class ReactionDiffusionSimulation
{
public: 
    void SetVerbose(int verboseLevel);// function for setting up the output information during the calculation process
    void MeshGeometry(double dX, double dY, double dZ, double gridSize ); // dimension unit in mm, function meshing the simulation geometry
    void CellStateInitialization(std::map<int, double> cX, std::map<int, double> cY,std::map<int, double> cZ, std::map<int, int>cState,std::map<int,double> secRMap);
    void DiffusionReactionCalculation(std::vector<std::vector<std::vector<double> > > &initialC,\
    double DiffC, double reactionRate, int Rt, double delta_T, double T); // concentration in molecules per grid,  in unit of mm^2/s, delta_T and T in unit of second

    std::vector<std::vector<std::vector<double> > > GetConcentration(); // function returning the Concentration at final time step
    std::vector<std::vector<std::vector<double> > > GetIntegralConcentration();// function returning the integral concentration at final time step
    
    void WriteConcentrationProfile(std::string path);// function write concentration profile in time domain

    int GetMeshDIMX();// function geting total grids number in x dimension
    int GetMeshDIMY();// function geting total grids number in y dimension
    int GetMeshDIMZ();// function geting total grids number in Z dimension
    
private:
    
    std::pair<std::map<int, std::vector<std::vector<std::vector<double> > > >,std::vector<std::vector<std::vector<double> > > > DiffusionReactionCalculation3D(std::vector<std::vector<std::vector<double> > > &initialC,\
    double DiffC, double reactionRate, int Rt, double delta_T, double T);
    std::pair<std::map<int, std::vector<std::vector<double> > >, std::vector<std::vector<double> > > DiffusionReactionCalculation2D(std::vector<std::vector<double> >&initialC,\
        double DiffC, double reactionRate, int Rt,double delta_T,double T
    );// concentration in molecules per grid, D in unit of mm^2/s, delta_T and T in unit of second
    
    int verbose;
    int gridNumX;
    int gridNumY;
    int gridNumZ;
    
    double xDim; // length in X dimension, unit in um
    double yDim;  // length in Y dimension, unit in um
    double zDim; // lenght in Z dimension, unit in um
    double d; // grid length, unit in um
    double deltaT; // time step size, unit in second

    struct IJKVec // define a struct to store the index of each cell of cellular automaton
    {
        int i;
        int j;
        int k;
    };
    
    struct IJVec
    {
        int i;
        int j;
    };
    
    
    std::vector<std::vector<IJVec> > gridIndex2D;
    std::map<int, IJVec> cellIndex2D;
    std::vector<std::vector<std::vector<IJKVec> > > gridIndex; // 3D vector containing the index of all the grids of geometry
    std::map<int, int> cellState; // A map containing the cell state, key is cell ID, value is cell state ID
    std::map<int, double>secreteRateMap; // a map containing the singal secretion rate of each cell, key is cell ID
    std::map<int, IJKVec> cellIndex; // A map containing the cell index in cellular automaton
    std::map<int, std::vector<std::vector<std::vector<double> > > > concentration; // A map containing the chemical Concentration at all grid cell
    std::vector<std::vector<std::vector<double> > > integralConcentration;
};

#endif