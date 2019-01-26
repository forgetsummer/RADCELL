// This program is direclty from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// Reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// $ID$
/// \file RunInitObserver.hh
/// \brief Definition of the RunInitObserver & RunInitManager classes

#ifndef RUNINITOBSERVER_HH_
#define RUNINITOBSERVER_HH_

#include "globals.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunInitObserver
{
public:
  RunInitObserver();
  virtual
  ~RunInitObserver();

  virtual void Initialize() = 0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunInitManager
{
public:
  static RunInitManager* Instance();
  void Initialize();

protected:
  friend class RunInitObserver;
  RunInitManager();
  ~RunInitManager();
  void Insert(RunInitObserver*);

  std::vector<RunInitObserver*> fObservers;
  static G4ThreadLocal RunInitManager* fgInstance;
};

#endif /* RUNINITOBSERVER_HH_ */
