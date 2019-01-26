#ifndef CellDNADamageAnalysis_h
#define CellDNADamageAnalysis_h 1
#include <iostream>
#include <vector>
#include <map>
#include "DBSCAN.hh"
#include "projectDataTypeGlobals.hh"

using namespace std;




class CellDNADamageAnalysis
{
public:
    void ReadEnergyDepositionFile(std::string filename);
    
    void ImportTotalParticleNumber(int num);
   
    void ImportDNADamagesTallyMapForWholeRun(EventDNASSBMap ssbMap, EventDNADSBMap dsbMap,int Num);
    void CellDNADamageTally();
    
    void WriteDNADamageTallyFile(std::string  outPutFileName);
    
    
///////////////////////////////////////////////////////////////////
    // The functions returning the DNA damages tally information after CellDNADamageAnalysis was implemented 
    
    std::map<int, double> GetCellSSBNumberMeamMap();
    std::map<int, double> GetCellSSBNumberStdMap();
    std::map<int, double> GetCellDSBNumberMeanMap();
    std::map<int, double> GetCellDSBNumberStdMap();
    std::map<int, std::vector<double> > GetCellDSBDimensionMap();
    
    
private:
    int totalParticleNumber;
    
    EventDNASSBMap SSBTallyMap;
    EventDNADSBMap DSBTallyMap;
    
    int N;
  
    std::map<int, double>cellSSBMeanMap;
    std::map<int, double>cellSSBStdMap;
    std::map<int, double>cellDSBMeanMap;
    std::map<int, double>cellDSBStdMap;
    
};



#endif 