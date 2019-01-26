#ifndef RadiationSourceSettings_h
#define RadiationSourceSettings_h 1
#include "G4String.hh"

using namespace std;
class RadiationSourceSettings
{
public:
    void RadiationPlaneSourceSetUp(std::string fileName, int particleNum, double particleEn, double cx, double cy, double cz, double xDim, double yDim);
    void RadiationPointSourceSetUp(std::string fileName, int particleNum, double particleEn,double cx, double cy, double cz);
    std::string GetSimulationSourceID();
    int GetSimulationParticleNumber();
    double GetSimulationParticleEnergy();
    
private:
    G4String simulationSourceID;
    int particleNumber;
    double particleEnergy;

};

#endif