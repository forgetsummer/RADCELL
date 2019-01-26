// This program is largely amended from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// Reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// $ID$
/// \file TrackingAction.hh
/// \brief Definition of the TrackingAction class

#ifndef TrackingAction_h
#define TrackingAction_h

#include "G4UserTrackingAction.hh"
#include "RunInitObserver.hh"
#include <map>

class G4Region;
class G4ParticleDefinition;

class TrackingAction : public G4UserTrackingAction, public RunInitObserver
{
public:
    TrackingAction();
    ~TrackingAction();

    void Initialize();

    virtual void PreUserTrackingAction(const G4Track*);

    std::map<const G4ParticleDefinition*, int>& GetNParticlesCreatedInTarget()
    {
        return fNParticleInTarget;
    }

    std::map<const G4ParticleDefinition*, int>& GetNParticlesCreatedInWorld()
    {
        return fNParticleInWorld;
    }

private:
    G4Region* fpTargetRegion;
    std::map<const G4ParticleDefinition*, int> fNParticleInTarget;
    std::map<const G4ParticleDefinition*, int> fNParticleInWorld;
};


#endif
