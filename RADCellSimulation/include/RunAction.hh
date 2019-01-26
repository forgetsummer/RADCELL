// This program is largely amended from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// Reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// $ID$
/// \file RunAction.hh
/// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <iostream>
#include <map>
#include <vector>

#include "EnergyDepositionEventInfo.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class TrackingAction;



class RunAction : public G4UserRunAction
{
public:
  
  RunAction();
  virtual ~RunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);



private:

  /////////////////
  // Histogramming
  //
  void CreateHistogram();
  void WriteHistogram();

  /////////////////
  // Worker
  //
  void BeginMaster(const G4Run*);
  void EndMaster(const G4Run*);

  /////////////////
  // Worker
  //
  void InitializeWorker(const G4Run*);
  void BeginWorker(const G4Run*);
  void EndWorker(const G4Run*);

  /////////////////
  // Print Info
  //
  void PrintRunInfo(const G4Run* run);

  /////////////////
  // Attributes
  //
  TrackingAction* fpTrackingAction;
  bool fInitialized;
  bool fDebug;
  

  int totalNumberOfEvent;


};
#endif
