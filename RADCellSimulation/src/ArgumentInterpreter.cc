#include "ArgumentInterpreter.hh"
#include <sstream>
#include <vector>
using namespace std;

G4String ArgumentInterpreter::argument;// Initialize the argument here
bool ArgumentInterpreter::postProcessMarker;


void ArgumentInterpreter::SetArgument(G4String arg)
{
    argument = arg;
}

G4String ArgumentInterpreter::GetSimulationMode()
{
    G4String argument_case;
    std::vector<G4String>argument_vector;
    stringstream ss(argument);
    while (ss>>argument_case)
    {
        argument_vector.push_back(argument_case);
    }
    
    return argument_vector[0];
}

G4String ArgumentInterpreter::GetOutPutFileName()
{
    G4String argument_case;
    std::vector<G4String>argument_vector;
    stringstream ss(argument);
    while (ss>>argument_case)
    {
        argument_vector.push_back(argument_case);
    }
    
    G4String filename;
    if (argument_vector[0]=="out" && argument_vector.size()==1)
    {
        filename = "microdosimetry";// default output filename is "microdosimetry"
    }
    if (argument_vector[0]=="out"&& argument_vector.size()==2)
    {
        filename = argument_vector[1]; // use the filename setup by user
    }
    
    return filename;

}

void ArgumentInterpreter::SetPostProcessMarker(bool postMarker)
{
    postProcessMarker = postMarker;
}
bool ArgumentInterpreter::GetPostProcessMarker()
{
    return postProcessMarker;
}



