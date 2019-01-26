#ifndef EnergyDepositionEventInfo_h
#define EnergyDepositionEventInfo_h 1
#include "globals.hh"

class EnergyDepositionEventInfo // define a class to store the energy deposition information
{
public:
  int eventID;// eventID indicating in which simulation event the edep happened
  int cellID; // cellID indicating in which cell energy deposition happened
  G4double pX; // position coordinate X
  G4double pY; // coordinate Y
  G4double pZ; // coordinate Z
  G4double edep;// energy deposited
  G4String affectedCellOrganelle;// cell organelle where energy deposition event happened
};



#endif

