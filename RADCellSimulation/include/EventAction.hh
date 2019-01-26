#ifndef EventAction_h
#define EventAction_h 1

#include "EnergyDepositionEventInfo.hh"
#include "G4UserEventAction.hh"
#include "projectDataTypeGlobals.hh"


#include <vector>
#include <map>
using namespace std;


class EventAction : public G4UserEventAction
{
public:
    EventAction(); // constructor
    virtual ~EventAction(); // destructor 
    
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    
    void  StoreEdepEvent(EnergyDepositionEventInfo edepEventInfo);//storing edep event into map 

private:

    EnergyDepositionEventMap edepEventMap;// edepEventMap , key is eventID
    
    void CellDoseTally(); // function for the dose tally for each particle
    void CellDNADamageTally(); //function for DNA damage tally for each particle
    
    
};





#endif
