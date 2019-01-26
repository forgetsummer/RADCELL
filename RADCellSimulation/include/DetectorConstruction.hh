// This program is largely amended from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// Reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// $ID$
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Ellipsoid.hh"
#include "G4UnionSolid.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4String.hh"
#include "globals.hh"
#include <map>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Region;
class Cell;

class EdepCellInformation
{
public:
  int cellID;
  G4String affectedCellOrganelle;
};
////////////////////////////////////////////////

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();

  virtual ~DetectorConstruction();

  virtual G4VPhysicalVolume* Construct();
  
  void DefineMaterials();

  G4Region* GetTargetRegion() {return fpRegion;}
    
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //// The first step: Register a general cell type for later palcement into tissue system
  
  G4String RegisterCellType(G4String cellType){return cellType;}// function for registering cell type
  G4String RegisterCellOrganelle(G4String cellOrganelle){return cellOrganelle;}//function for registeriing cell organs
  G4String RegisterCellShape(G4String cellShape){return cellShape;}// functions for registering cell shape
  G4String RegisterCellOrganelleColor(G4String organelleColor){return organelleColor;}// function registering color of organelle
  void   CreateCell(G4String cellType, G4String cellOraganelle, G4String cellShape,G4String organelleColor); // function creating a general cell and register it to generalCellMap
 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///// The second step: Place the cell into tissue system according to its position
  
  void SetCellWorld(G4double X, G4double Y, G4double Z);// create a cell world as tissue 
  G4AssemblyVolume* CreateCellAssembly (G4String cellType, G4double cellSize, G4double nucleusSize);// function creating cell assembly according to cell sizes
  void SetCellPlacement (int cellID,G4String cellType,G4ThreeVector position,G4double cellSize,G4double nucleusSize);// The size of cytoplasma and nuclues will change with cell cycle
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  /////////////////////////////////////////////////////////////////////////////////////////
  /////The Third step: Accessing energy deposition information in cell

  EdepCellInformation  GetEdepCellInformation(G4VPhysicalVolume* physicalVolume);// function getting enery deposition information in cell
 
  /////////////////////////////////////////////////////////////////////////////////////////////
  
  void SetCellSimulationParameter(G4double defaultCut, G4String targetOrganelle);// in this function set cell simulation parameters, such as energy cutting value, and fpRegion values
  
  void ClearCellStore();
  
  int GetTestNum()const {return generalCellMap.size();};
  
  std::map<int, double>  GetCellMassMap()const;
  std::map<int, double>  GetNucleusMassMap() const;

                         
protected:
  G4Material*        fpWaterMaterial;//assume the cell material is water
  G4Material*        fpAirMaterial; // assume the work material is air for single cell DNA damage simulation
  G4Region*          fpRegion;
  
  G4VPhysicalVolume* physicalWorld;//define a physical volume of world
  
  std::map<G4String, Cell> generalCellMap;// map storing a general cell, int is the cell type ID
  std::vector<G4String> targetOrganelleVector; // vector storing the target volume which will use Geant4DNA crosssection
  
  std::map<int, double> cellMassMap; // a map storing the mass of each cell created for simulation
  std::map<int, double> nucleusMassMap;
  

  
     
};


#endif
