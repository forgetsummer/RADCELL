#ifndef projectDataTypeGlobals_h
#define projectDataTypeGlobals_h 1

#include "EnergyDepositionEventInfo.hh"
#include "DNADamageTypes.hh"

#include <vector>
#include <map>
#include <utility> 
using namespace std;


// this class if for defining some data types used in the project

typedef std::vector<EnergyDepositionEventInfo> edepEventVector;
typedef std::map<int, edepEventVector> EnergyDepositionEventMap;// key is eventID

typedef std::map<std::string, double> OrganelleDoseMap; //  key is name of organelle, value is dose
typedef std::map<int,OrganelleDoseMap> CellDoseTallyMap; //  key is cellID
typedef std::map<int, CellDoseTallyMap> EventDoseTallyMap; //  key is eventID


typedef std::vector<EnergyDepositionPoint> EdepPoints; // define vector to store energy deposition point information
typedef std::map<int,EdepPoints> CellEdepMap; // define map of edep points of a cell, key is cell ID
typedef std::map<int, CellEdepMap> EventEdepMap; // define map of edep points of a event, key is event ID

typedef std::vector<SSB>SSBVector;//define a vector to store the SSBs
typedef std::vector<DSB>DSBVector;// define a vector to store the DSBs
typedef std::map<int,SSBVector> CellDNASSBMap; // define a map storing SSB for a cell, key is cell ID
typedef std::map<int,DSBVector> CellDNADSBMap; // define a map stroing DSB for a cell, key is cell ID
typedef std::map<int, CellDNASSBMap> EventDNASSBMap; //define a map storing the DNA SSB damages, the key is event ID
typedef std::map<int, CellDNADSBMap> EventDNADSBMap; //define a map storing the DNA DSB damages, the key is event ID


// typedef std::pair< std::map<int, std::map<int, double> >,std::map<int,std::map<int,double> > >EventSSBDSBPairMap; // define a pair to store SSB and DSB map


    
struct EventDNADamagesTallyInfo
{
    EventDNADSBMap eventDSBMap; // this a map storing ssb tally information
    EventDNASSBMap eventSSBMap; // this is a map storing dsb tally information
    int N; // this is edep points re-use factor, usually you can set this as 1, sometimes 10, or 100, etc
};


#endif
