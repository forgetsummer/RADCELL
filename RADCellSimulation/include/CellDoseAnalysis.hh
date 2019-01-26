#ifndef CellDoseAnalysis_h
#define CellDoseAnalysis_h 1


#include "projectDataTypeGlobals.hh"

#include <vector>
#include <map>
using namespace std;



class CellDoseAnalysis
{
public:
    void ReadEnergyDepositionFile(std::string filename);
    void WriteCellDoseAnalysisFile(std::string doseTallyFileName);
    void DoseTally();
    void ImportDoseTallyMapForWholeRun(EventDoseTallyMap doseTallyMap);
    void ImportTotalParticleNumber(int num);
    std::map<int, double> GetTotalCellDoseMeanMap(); // function returning mean cell dose map, key is cellID
    std::map<int, double> GetTotalCellDoseStdMap();
    std::map<int, double> GetTotalNucleusDoseMeanMap();
    std::map<int, double> GetTotalNucleusDoseStdMap();
    
    
private:
    int totalParticleNumber;
    EventDoseTallyMap doseTallyMapForWholeRun;
    
    std::map<int, double> cellDoseMeanMap;// map storing the mean cell total dose, the key is cellID, value is total dose of cell
    std::map<int, double> cellDoseStdMap; // map storing the statistical error of mean cell total dose, the key is cellID
    std::map<int, double>nucleusDoseMeanMap;
    std::map<int, double>nucleusDoseStdMap;
     
    
};
#endif 