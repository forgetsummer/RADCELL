// 
// This code is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 

#include "RADCellSimulation.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
#include "ActionInitialization.hh"
#include "PhysicsList.hh"
#include "G4UIQt.hh"
#include "G4SystemOfUnits.hh"
#include "G4String.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RegionStore.hh"
#include "G4Material.hh"
#include "SteppingAction.hh"
#include "RunAction.hh"
#include "ArgumentInterpreter.hh"
#include <sstream>
#include <iostream>
using namespace std;

RADCellSimulation::RADCellSimulation()
{
    detector = new DetectorConstruction();
}


RADCellSimulation::~RADCellSimulation()
{
    delete visManager;
    delete runManager;
    delete ui;

}


void RADCellSimulation::SetCellWorld(G4double X, G4double Y, G4double Z)
{
    detector->DefineMaterials();//define materials of cell system
    G4double XWithUnit=X*mm;
    G4double YWithUnit=Y*mm;
    G4double ZWithUnit=Z*mm;
    detector->SetCellWorld(XWithUnit,YWithUnit,ZWithUnit);// set up a world which contains cells;
}

void RADCellSimulation::CellConstruction(int cellID, string cellType, G4double X, G4double Y, G4double Z, G4double cellRadius,G4double cellNucleusRadius)
{
  G4double XWithUnit=X*mm;//define the position of cell  in unit of mm
  G4double YWithUnit=Y*mm;
  G4double ZWithUnit=Z*mm;
  G4double cellRadiusWithUnit=cellRadius*um; // define the dimension of cell in unit of um
  G4double cellNucleusRadiusWithUnit=cellNucleusRadius*um; // define the dimension of nucleus in unit of micronmeter
 
  G4ThreeVector cellPosition=G4ThreeVector(XWithUnit,YWithUnit,ZWithUnit);
  detector->SetCellPlacement(cellID, cellType, cellPosition,cellRadiusWithUnit, cellNucleusRadiusWithUnit);//place cells inside cell world
     
}

void RADCellSimulation::CreateCell(string cellType, string cellOrganelle, string cellShape,string organelleColor)
{
  detector->CreateCell(cellType, cellOrganelle,cellShape,organelleColor);
}

void RADCellSimulation::SetCellSimulationParameter(G4double defaultCut, string targetOrganelle)
{
  detector->SetCellSimulationParameter(defaultCut,targetOrganelle);
}

void RADCellSimulation::RADCellSimulationInitializePyWrapper(int argc,const std::vector<std::string> & _vec_str)
{
    cerr<<"argc="<<argc<<endl;
    cerr<<"_vec_str.size()="<<_vec_str.size()<<endl;
    char** argv  = (char**)&_vec_str[0];
    cerr<<"RADCellSimulationInitializePyWrapper argv  [0]="<<argv[0]<<endl;
    
    for (int i = 0 ; i < _vec_str.size() ; ++i){
        cerr<<"argument["<<i<<"]="<<argv[i]<<endl;
    }
    
    RADCellSimulationInitialize(argc,argv);  
}
//----------------------------------------------------------------------
// Below is for function of EnergyDistributionCalculation()
//----------------------------------------------------------------------



void RADCellSimulation::RADCellSimulationInitialize(int argc, char** argv)
{
    // Choose the Random engine
  
  G4Random::setTheEngine(new CLHEP::RanecuEngine);


  #ifdef G4MULTITHREADED
   // G4MTRunManager* runManager = new G4MTRunManager;
    runManager= new G4MTRunManager;
    runManager->SetNumberOfThreads(4);// set up the number of number of threads for multithread simulation
  #else
   // G4RunManager* runManager = new G4RunManager;
    runManager= new G4RunManager;
  #endif
  
  runManager->SetUserInitialization(detector);//here we use the internal variable detector in class 


  runManager->SetUserInitialization(new PhysicsList);

  // User action initialization
  
  runManager->SetUserInitialization(new ActionInitialization(detector));
  


  // Initialize G4 kernel
  
    runManager->Initialize();
 
  // Initialize visualization 
  visManager = new G4VisExecutive;

  visManager->Initialize();  // comment this if you want to run radcell in CC3D
  ui = new G4UIExecutive(argc, argv); // 
  UImanager = G4UImanager::GetUIpointer();  // Get the pointer to the User Interface manager
  UImanager->ApplyCommand("/control/execute vis.mac");

}

void RADCellSimulation::UpdateGeometryInitialize()
{
   runManager->ReinitializeGeometry();
//    G4GeometryManager::GetInstance()->OpenGeometry();
//    G4RegionStore::GetInstance()->Clean();// delete the regions in store
   detector->ClearCellStore();
   
}
void RADCellSimulation::UpdateGeometryFinalize()
{
   
  runManager->DefineWorldVolume(detector->Construct());

  runManager->PhysicsHasBeenModified();

}


void RADCellSimulation::EnergyDistributionCalculation(string runArgument,string inputFile,bool postProMarker)
{

    stringstream argument(runArgument); // get run argument
    G4String argument_case;
    std::vector<G4String>argument_vector;
    stringstream commandArgument;
    commandArgument<<"/control/execute" <<" "<<inputFile;
    string runCommand;
    runCommand=commandArgument.str();
    
  
    
    while (argument>>argument_case)
    {
        argument_vector.push_back(argument_case);
       
    }
    ArgumentInterpreter::SetArgument(runArgument);
    ArgumentInterpreter::SetPostProcessMarker(postProMarker);
    
    
    if (argument_vector[0]=="gui")// running in gui interactive mode 
    {
          // Define UI terminal for interactive mode  
        //G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        
        UImanager->ApplyCommand(runCommand);
        ui->SessionStart(); 
       // delete ui;
    }
    if (argument_vector[0]=="out")// running in output file mode
    { 
        UImanager->ApplyCommand(runCommand);
        
//         UImanager->ApplyCommand("/control/execute microdosimetry.in");
        
    }

  
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
