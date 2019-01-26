#ifndef  DNADamageTypes_h
#define  DNADamageTypes_h 1

// This is a header file containing all the DNA damage types, such as SSB and DSB

class EnergyDepositionPoint
{
    public:
        double edep; // energy deposited in the energy deposition point
        double x;
        double y;
        double z;

};
class SSB
{
    public:
    EnergyDepositionPoint edepPoint; // edepPoint is edep information for this SSB
    int DNALabel; // DNALabel is for indicating which DNA strand is broke by this edep point, could be 0 or 1
};
class DSB
{
    public:
        double DSBTotalEnergy;
        double DSBDimension;
        double xC;// center of DSB, x
        double yC; // center of DSB,y
        double zC;// center of DSB , z
        int DSBType; // label for indicating DSB types, 1 is simple DSB, 2 is complex DSB
};



#endif