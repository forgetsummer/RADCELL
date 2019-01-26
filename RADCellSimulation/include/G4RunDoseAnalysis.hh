#ifndef G4RunDoseAnalysis_h
#define G4RunDoseAnalysis_h 1

#include "DetectorConstruction.hh"
#include "projectDataTypeGlobals.hh"
#include <string>


using namespace std;

class G4RunDoseAnalysis
{
public:
    static void CollectEventDoseTallyMap(EventDoseTallyMap doseTallyMapPerEvent); // function collecting the dose tally infomation for each object of G4Event
    void CellDoseTally(const DetectorConstruction *detector);
    void ImportTotalParticleNumber(int num);
    void RadiationTransportDoseOutPut(std::string filename); // function for outputting the radiation dose tally information for whole Run
    void ClearRunData();
    int  GetSizeOfDoseTallyMap(){return doseTallyMap.size();};
    
    std::map<int, double> GetCellDoseMeanMap(){return cellDoseMeanMap;};
    std::map<int, double> GetCellDoseStdMap(){return cellDoseStdMap;};
    std::map<int, double> GetNucleusDoseMeanMap(){return nucleusDoseMeanMap;};
    std::map<int, double> GetNucleusDoseStdMap(){return nucleusDoseStdMap;};
private:
    static EventDoseTallyMap doseTallyMap;
    int totalParticleNumber;
    std::map<int, double> cellDoseMeanMap;// map storing the mean cell total dose, the key is cellID, value is total dose of cell
    std::map<int, double> cellDoseStdMap; // map storing the statistical error of mean cell total dose, the key is cellID
    std::map<int, double> nucleusDoseMeanMap; // map storing the mean dose of nuclues, key is cellID, value is nuclues mean dose
    std::map<int, double> nucleusDoseStdMap; //
};




#endif