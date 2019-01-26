//
// This code is directly used from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// Reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// $ID$
/// \file RunInitObserver.cc
/// \brief Implementation of the RunInitObserver & RunInitManager classes

#include "RunInitObserver.hh"

G4ThreadLocal RunInitManager* RunInitManager::fgInstance(0);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunInitManager::RunInitManager()
{
  fgInstance = this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunInitManager::~RunInitManager()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunInitObserver::RunInitObserver()
{
  RunInitManager::Instance()->Insert(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunInitObserver::~RunInitObserver()
{
  // TODO Auto-generated destructor stub
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunInitManager* RunInitManager::Instance()
{
  if(fgInstance == 0) new RunInitManager();
  return fgInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
RunInitManager::Initialize()
{
  std::vector<RunInitObserver*>::iterator it = fObservers.begin();
  std::vector<RunInitObserver*>::iterator end = fObservers.end();
  for(; it != end ; it++) (*it)->Initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
RunInitManager::Insert(RunInitObserver* observer)
{
  fObservers.push_back(observer);
}
