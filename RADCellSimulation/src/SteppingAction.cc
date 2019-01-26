
// This program is largely amended from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $ID$
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "Analysis.hh"
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include <vector>
#include <sstream>
#include "ArgumentInterpreter.hh"
#include "G4EventManager.hh" // using this class to get the current event ID

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* detectorConstruction,EventAction* eventAction)
  :G4UserSteppingAction(),fDetConstruction(detectorConstruction),fEventAction(eventAction)/// We should use construct initializer
{
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{ 
  G4double flagParticle=0.;
  G4double flagProcess=0.;
  G4double x,y,z,xp,yp,zp;

  // Process sub-types are listed in G4PhysicsListHelper.cc

  /*
 // The following method avoids the usage of string comparison 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ==
 G4Electron::ElectronDefinition())
    flagParticle = 1; 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ==
 G4Proton::ProtonDefinition())
    flagParticle = 2; 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ==
 G4Alpha::AlphaDefinition())
    flagParticle = 4; 

    G4DNAGenericIonsManager *instance;
    instance = G4DNAGenericIonsManager::Instance();
    G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
    G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
    G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
    G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
    G4ParticleDefinition* heliumDef = instance->GetIon("helium");

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ==
 instance->GetIon("hydrogen"))
    flagParticle = 3; 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ==
 instance->GetIon("alpha+"))
    flagParticle = 5; 

 if (step->GetTrack()->GetDynamicParticle()->GetDefinition() ==
 instance->GetIon("helium"))
    flagParticle = 6; 
   */

  const G4String& particleName = step->GetTrack()->GetDynamicParticle()->
      GetDefinition()->GetParticleName();
  const G4String& processName = step->GetPostStepPoint()->
      GetProcessDefinedStep()->GetProcessName();

  if (particleName == "e-")       flagParticle = 1;
  else if (particleName == "proton")   flagParticle = 2;
  else if (particleName == "hydrogen") flagParticle = 3;
  else if (particleName == "alpha")    flagParticle = 4;
  else if (particleName == "alpha+")   flagParticle = 5;
  else if (particleName == "helium")   flagParticle = 6;

  if (processName=="msc")      flagProcess =10;
  else if (processName=="e-_G4DNAElastic")    flagProcess =11;
  else if (processName=="e-_G4DNAExcitation")    flagProcess =12;
  else if (processName=="e-_G4DNAIonisation")    flagProcess =13;
  else if (processName=="e-_G4DNAAttachment")    flagProcess =14;
  else if (processName=="e-_G4DNAVibExcitation")  flagProcess =15;
  else if (processName=="eCapture")      flagProcess =16;

  else if (processName=="proton_G4DNAExcitation")  flagProcess =17;
  else if (processName=="proton_G4DNAIonisation")  flagProcess =18;
  else if (processName=="proton_G4DNAChargeDecrease")  flagProcess =19;

  else if (processName=="hydrogen_G4DNAExcitation")   flagProcess =20;
  else if (processName=="hydrogen_G4DNAIonisation")   flagProcess =21;
  else if (processName=="hydrogen_G4DNAChargeIncrease")flagProcess =22;

  else if (processName=="alpha_G4DNAExcitation")  flagProcess =23;
  else if (processName=="alpha_G4DNAIonisation")  flagProcess =24;
  else if (processName=="alpha_G4DNAChargeDecrease")  flagProcess =25;

  else if (processName=="alpha+_G4DNAExcitation")  flagProcess =26;
  else if (processName=="alpha+_G4DNAIonisation")  flagProcess =27;
  else if (processName=="alpha+_G4DNAChargeDecrease")  flagProcess =28;
  else if (processName=="alpha+_G4DNAChargeIncrease")  flagProcess =29;

  else if (processName=="helium_G4DNAExcitation")  flagProcess =30;
  else if (processName=="helium_G4DNAIonisation")  flagProcess =31;
  else if (processName=="helium_G4DNAChargeIncrease")  flagProcess =32;

  else if (processName=="hIoni")  flagProcess =33;
  else if (processName=="eIoni")  flagProcess =34;

  if (processName!="Transportation")
  {
    x=step->GetPreStepPoint()->GetPosition().x()/nanometer;
    y=step->GetPreStepPoint()->GetPosition().y()/nanometer;
    z=step->GetPreStepPoint()->GetPosition().z()/nanometer;
    xp=step->GetPostStepPoint()->GetPosition().x()/nanometer;
    yp=step->GetPostStepPoint()->GetPosition().y()/nanometer;
    zp=step->GetPostStepPoint()->GetPosition().z()/nanometer;

//     CommandLineParser* parser = CommandLineParser::GetParser();
//     Command* command(0);
//     if((command = parser->GetCommandIfActive("-out"))==0) return;
    G4String simulationMode = ArgumentInterpreter::GetSimulationMode();

    if (simulationMode=="gui") return; // if it is gui mode, then return, while if it is out , then continue
    
    
   //////////////////////////////////////////////
   
   // Collect energy and energy deposition position step by step

  // get volume of the current step
    G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
    G4double edep = step->GetTotalEnergyDeposit()/eV;// get deposited energy in eV
    
    int cellID=fDetConstruction->GetEdepCellInformation(volume).cellID;
    
    G4String affectedCellOrganelle=fDetConstruction->GetEdepCellInformation(volume).affectedCellOrganelle;
    int eventID=G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID(); // get evendID here
    
    EnergyDepositionEventInfo edepEventInfo;// make an instance of EnergyDepositionEventInfo class
    edepEventInfo.cellID=cellID;
    edepEventInfo.eventID=eventID;//store the eventID in to EnergyDepositionEventInfo class
    edepEventInfo.pX=x;
    edepEventInfo.pY=y;
    edepEventInfo.pZ=z;
    edepEventInfo.edep=edep;
    edepEventInfo.affectedCellOrganelle=affectedCellOrganelle;
    
    fEventAction->StoreEdepEvent(edepEventInfo);

    // get analysis manager
    //cout<<"the first argument is "<<argument_vector[0]<<endl;
    
     bool postProcessMarker=ArgumentInterpreter::GetPostProcessMarker();
     if (postProcessMarker==true) // if needed to write out root files, then do this
     {
        G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    //    if(!analysisManager->IsActive()) {return; }

        // fill ntuple
        analysisManager->FillNtupleDColumn(0, flagParticle);
        analysisManager->FillNtupleDColumn(1, flagProcess);
        analysisManager->FillNtupleDColumn(2, x);
        analysisManager->FillNtupleDColumn(3, y);
        analysisManager->FillNtupleDColumn(4, z);
        analysisManager->FillNtupleDColumn(5, step->GetTotalEnergyDeposit()/eV);
        analysisManager->FillNtupleDColumn(6,
                                        std::sqrt((x-xp)*(x-xp)+
                                            (y-yp)*(y-yp)+(z-zp)*(z-zp))/nm);
        analysisManager->FillNtupleDColumn(7, (step->GetPreStepPoint()->
            GetKineticEnergy() - step->GetPostStepPoint()->GetKineticEnergy())/eV );

        
        analysisManager->FillNtupleDColumn(8,cellID);
        analysisManager->AddNtupleRow();
     }

  }
    
  
  //cerr<<"DONE WITH STEPPING ACTION"<<endl;
  
}    
