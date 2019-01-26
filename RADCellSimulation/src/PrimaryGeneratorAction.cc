//
// This program is largely amended from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// Reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// $ID$
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction():
G4VUserPrimaryGeneratorAction()
{
  //int n_particle = 1;
  //~ fpParticleGun  = new G4ParticleGun(n_particle);
//~ 
  //~ // default gun parameters
  //~ fpParticleGun->SetParticleEnergy(10.*keV);
  //~ fpParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  //~ fpParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-0.1*mm));
  fpParticleGun = new G4GeneralParticleSource();
  
  fpParticleGun->SetVerbosity(0);// set Verbosity level for printing the particle information on screen
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fpParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  fpParticleGun->GeneratePrimaryVertex(anEvent);
}
