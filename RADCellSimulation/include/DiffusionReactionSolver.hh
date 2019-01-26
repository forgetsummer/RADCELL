#ifndef DiffusionReactionSolver_h 
#define DiffusionReactionSolver_h 1

#include <vector>
#include <map>
#include <string>
using namespace std;

class DiffusionReactionSolver
{
public: 
    void SetVerbose(int verboseLevel); // function to set up the output information, when verboseLevel=1, output information, verboseLevel=0 not
    void DiffusionReactionInitialization(double dX, double dY, double dZ, double gridSize,double DiffC, double reactionRate, \
        int receptorNum,double delta_T);//initilize the diffusion reaction system
    void CellStateUpdate(int cellID, double x, double y, double z, int state,double fmu);// set the cell stat
    void CellStateUpdateInitialize();// initilize cell state when cell state needed to be updated
    void DiffusionReactionCalculation(double deltaT_diffusion, int frequency);// diffusion reaction calculation
    double GetConcentration(double x, double y, double z);// function returning the chemical concentration at location (x,y,z)
    double GetIntegralConcentration(double x, double y, double z);
    void WriteConcentrationToFile(std::string path, int fileID);

    
private:
    int verbose;
    double dimX;
    double dimY;
    double dimZ;
    double d;
    double D;
    double r;
    int Rt;
    double mu;
    double dt;
    int N_X;
    int N_Y;
    int N_Z;
    std::map<int, double> cX;// map containing the cell position information, x
    std::map<int, double> cY;// map containing the cell position information, y
    std::map<int, double> cZ;// map containing the cell position information ,z
    std::map<int, int> cellState; // A map containing the cell state, key is cell ID, value is cell state ID
    std::map<int, double>secreteRate; // a map containing the cell signal emitting rate , key is cell id, value is mu
    std::vector<std::vector<std::vector<double> > >  concentration; // A vector containing the chemical Concentration at all grid cell
    std::vector<std::vector<std::vector<double> > >  integralConcentration;
    std::vector<std::vector<std::vector<double> > > initialC; 
};



#endif