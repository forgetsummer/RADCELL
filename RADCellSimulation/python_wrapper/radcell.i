                  
// Module Name
%module ("threads"=1) radcell
%import <G4Types.hh>
// ************************************************************
// Module Includes
// ************************************************************

// These are copied directly to the .cxx file and are not parsed
// by SWIG.  Include include files or definitions that are required
// for the module to build correctly. 

%apply char **STRING_ARRAY { char ** }
      
%{    
                
// #include <G4String.hh> 
    
//#include <RADCellSimulationPyWrapper.hh>   
#include <RADCellSimulation.hh>
#include <CellLayoutInitializer.hh>          
#include <CellDoseAnalysis.hh>
#include <CellDNADamageAnalysis.hh> 
#include <ReadRadiationTransportInfo.hh>
#include <DiffusionReactionSolver.hh>


     
     
 

// Namespaces
using namespace std;

#include <string>

%}

// #define XMLUTILS_EXPORT
  


// C++ std::string handling
%include "std_string.i"
%include <std_string.i>
  
 
// C++ std::map handling         
%include "std_map.i"
 
// C++ std::map handling
%include "std_vector.i"      

//C++ std::list handling
%include "std_list.i"

%template (vectorstring) std::vector<std::string>;



//%include <RADCellSimulationPyWrapper.hh>
%include <RADCellSimulation.hh>
%include <CellLayoutInitializer.hh>    
%include <CellDoseAnalysis.hh>
%include <CellDNADamageAnalysis.hh>
%include <ReadRadiationTransportInfo.hh>
%include <DiffusionReactionSolver.hh>

//%include <G4Types.hh>


// %include <G4String.hh>