#ifndef MonteCarloErrorAnalysis_h
#define MonteCarloErrorAnalysis_h 1
#include <vector>
#include <map>
using namespace std;

class  MonteCarloErrorAnalysis
{
public:
    void ImportSimulationResults(std::vector<std::vector<double> > result);
    std::vector<double> GetMeanOfSimulationResults();
    std::vector<double> GetStdOfSimulatioinResults();
private:
    std::vector<std::vector<double> > simResult; // a matrix storing all the simulation results for all each run history

};



#endif