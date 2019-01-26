//
// This code is directly used from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// Reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// $ID$
/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "TrackingAction.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization(DetectorConstruction* detConstruction)
: G4VUserActionInitialization(),fDetConstruction(detConstruction) // using construct initializer here !!!
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
 // In MT mode, to be clearer, the RunAction class for the master thread might
 // be different than the one used for the workers.
 // This RunAction will be called before and after starting the
 // workers.
 // For more details, please refer to :
 //https://twiki.cern.ch/twiki/bin/view/Geant4/Geant4MTForApplicationDevelopers
 //

  SetUserAction(new RunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
  // G4cout << "Build for = "
  // << G4RunManager::GetRunManager()->GetRunManagerType()
  // << G4endl;

  SetUserAction(new PrimaryGeneratorAction);

  TrackingAction* trackingAction = new TrackingAction();
  SetUserAction(trackingAction);

  RunAction* runAction= new RunAction();
  SetUserAction(runAction);
  EventAction* eventAction=new EventAction();
  SetUserAction(eventAction);
  

  SetUserAction(new SteppingAction(fDetConstruction,eventAction));
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
