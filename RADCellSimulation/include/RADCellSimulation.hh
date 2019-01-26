// 
// This code is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015
// This is a class for running the radiation transportation in cells

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "DetectorConstruction.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4Types.hh"
#include <string>

using namespace std;


class RADCellSimulation 
 {
      public: 
      RADCellSimulation();// the constructor
      ~RADCellSimulation(); // the destructor
      void CreateCell(string cellType, string cellOrganelle, string cellShape,string organelleColor); // the function for creating one type of cell
      void SetCellWorld(G4double X, G4double Y, G4double Z); // the function for seting the dimension of world containing all the cells
      void CellConstruction(int cellID,string cellType, G4double X, G4double Y, G4double Z, G4double cellRadius,G4double cellNucleusRadius); 
      void SetCellSimulationParameter(G4double defaultCut, string targetOrganelle);
      void RADCellSimulationInitialize(int argc, char** argv); // function initializing the radiation transport simulation
      void RADCellSimulationInitializePyWrapper(int argc,const std::vector<std::string> & _vec_str); // for python wrapping RADCellSimulationInitialize()
      void EnergyDistributionCalculation(string runArgument,string inputFile,bool postProMarker); // function for running the radiation transport, the arguments are the running mode.
      void UpdateGeometryInitialize(); // the function initializing the process of changing the cell geometry between different runs
      void UpdateGeometryFinalize(); // the function finalizing the process of changing the cell geometry between different runs

      
      
      private:
      
      DetectorConstruction* detector;
      G4VisManager* visManager;
        #ifdef G4MULTITHREADED
        G4MTRunManager* runManager;
        #else
        G4RunManager* runManager;
        #endif
      G4UImanager* UImanager;
      G4UIExecutive* ui; 

      
      

 };
	
