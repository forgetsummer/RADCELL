#ifndef G4EventDNADamageAnalysis_h
#define G4EventDNADamageAnalysis_h 1

#include "projectDataTypeGlobals.hh"



class G4EventDNADamageAnalysis
{

public:
    
    EventDNADamagesTallyInfo EventDNADamageTally(EnergyDepositionEventMap edepEventMap,double SSBProb,double SSBEmin, double SSBEmax,double eps, int MinPts, double p); // function to tally the DNA damage for each particle event


private:
    EventDNASSBMap EventSSBTally(EventEdepMap edepMap, double SSBEmin, double SSBEmax);// function doing SSB tally, returning EvenDNASSBTallyMap
    void EventDSBTally(EventDNASSBMap &potentialSSBTallyMap,double probOfInSensitiveRegion, double eps, int MinPts, double p);
    DSBVector GetDSBCluster(std::map<int, std::vector<int> > & SSBClusterMap, SSBVector& SSBVec);
    double GetDSBDimension( std::vector<int>& SSBPts, SSBVector& SSBVec);
    
    double Rand2();// define a function to produce pseudo random number for monte carlo sampling 
    
    EventDNASSBMap SSBTallyMap;
    EventDNADSBMap DSBTallyMap;
};





#endif