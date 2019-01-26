#ifndef G4RunDNADamageAnalysis_h
#define G4RunDNADamageAnalysis_h 1


#include "projectDataTypeGlobals.hh"
#include <string>


using namespace std;

class G4RunDNADamageAnalysis
{
public:
    static void CollectEventDNADamageTallyMap(EventDNADamagesTallyInfo dnaDNADamagesTallyInfoPerEvent);
    
    void CellDNADamageTally(std::map<int, double> cellDoseMeanMap,std::map<int, double> cellDoseStdMap);
    void ImportTotalParticleNumber(int num);
    
    void RadiationTransportDNADamageOutPut(std::string filename);
    
    std::map<int, double> GetCellDSBMeanMap(){return cellDSBMeanMap;};
    std::map<int, double> GetCellDSBStdMap(){return cellDSBStdMap;};
    std::map<int, double> GetCellSSBMeanMap(){return cellSSBMeanMap;};
    std::map<int, double> GetCellSSBStdMap(){return cellSSBStdMap;};
    std::map<int, std::vector<double> > GetCellDSBDimensionMap(){return cellDSBDimensionMap;};
    
    std::map<int, double> GetDSBPerGyPerCellMeanMap() {return ratio_DSB_Gy_Gbp_MeanMap;};
    std::map<int, double> GetDSBPerGyPerCellStdMap() {return ratio_DSB_Gy_Gbp_StdMap;};
    std::map<int, double> GetSSBPerGyPerCellMeanMap() {return ratio_SSB_Gy_Gbp_MeanMap;};
    std::map<int, double> GetSSBPerGyPerCellStdMap() {return ratio_SSB_Gy_Gbp_StdMap;};
    
    void ClearRunData();
    
private:
    
    static EventDNASSBMap SSBTallyMap;
    static EventDNADSBMap DSBTallyMap;
    static int N;


    int totalParticleNumber;
    
    std::map<int, double> cellDSBMeanMap; // map storing the mean DSB number of cell, key is cellID, value is DSB number
    std::map<int, double> cellDSBStdMap; // map storing the DSB number std of cell, key is cell ID
    std::map<int, double> cellSSBMeanMap;
    std::map<int, double> cellSSBStdMap;
    std::map<int, std::vector<double> > cellDSBDimensionMap;
    
    std::map<int, double> ratio_DSB_Gy_Gbp_MeanMap;
    std::map<int, double> ratio_DSB_Gy_Gbp_StdMap;
   
    std::map<int, double> ratio_SSB_Gy_Gbp_MeanMap;
    std::map<int, double> ratio_SSB_Gy_Gbp_StdMap;
    
    std::map<int, double> cellDoseMeanMap;
    
};





#endif