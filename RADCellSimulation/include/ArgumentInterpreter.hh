// 
// This is code is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015

#ifndef ArgumentInterpreter_h
#define ArgumentInterpreter_h 1
#include "G4String.hh"

struct stat;class ArgumentInterpreter
{
    public:
        static void SetArgument(G4String arg);
        static G4String GetSimulationMode(); // function returning the simulation mode, two modes: gui, out
        static void SetPostProcessMarker(bool postMarker);
        static bool GetPostProcessMarker();
        static G4String GetOutPutFileName();
    private:
        static G4String argument;
        static bool postProcessMarker;
        
};
#endif