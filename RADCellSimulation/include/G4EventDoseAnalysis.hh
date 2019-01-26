#ifndef G4EventDoseAnalysis_h
#define G4EventDoseAnalysis_h 1

#include "projectDataTypeGlobals.hh"



class G4EventDoseAnalysis
{
public:
    EventDoseTallyMap EventDoseTally(EnergyDepositionEventMap edepEventMap);// function to tally the dose for each particle event,
};







#endif 