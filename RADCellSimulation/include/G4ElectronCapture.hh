// This program is direclty used from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// Reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// $ID$
/// \file G4ElectronCapture.hh
/// \brief Definition of the G4ElectronCapture class

#ifndef ElectronCapture_h
#define ElectronCapture_h 1

#include "G4VDiscreteProcess.hh"
#include "globals.hh"
#include "G4ParticleChangeForGamma.hh"

class G4Region;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ElectronCapture : public G4VDiscreteProcess
{
public:

  G4ElectronCapture(const G4String& regName, G4double ekinlimit);

  virtual ~G4ElectronCapture();

  void SetKinEnergyLimit(G4double);

  virtual void BuildPhysicsTable(const G4ParticleDefinition&);

  virtual G4bool IsApplicable(const G4ParticleDefinition&);

  virtual G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
      G4double previousStepSize,
      G4ForceCondition* condition);

  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&);

protected:

  virtual G4double GetMeanFreePath(const G4Track&, G4double,G4ForceCondition*);

private:

  // hide assignment operator as private
  G4ElectronCapture(const G4ElectronCapture&);
  G4ElectronCapture& operator = (const G4ElectronCapture &right);

  G4double fKinEnergyThreshold;
  G4String fRegionName;
  G4Region* fpRegion;
  G4ParticleChangeForGamma fParticleChange;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

