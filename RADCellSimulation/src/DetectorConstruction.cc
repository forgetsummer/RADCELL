// This code is largely amended from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// Reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// $ID$
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "Randomize.hh"  // uisng random number
#include <math.h> //uisng math function
#include <stdlib.h>
#include "DetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "Cell.hh"
#include <iostream>
#include <sstream>

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),fpRegion(0),physicalWorld(0)
{
    DefineMaterials();// define the material in the constructor 
    fpRegion = new G4Region("Target");// declare a new region for target
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fpRegion;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()

{
//   DefineMaterials();
  return physicalWorld;//return the physical world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 

  // Water is defined from NIST material database
  G4NistManager * man = G4NistManager::Instance();
  G4Material * H2O = man->FindOrBuildMaterial("G4_WATER");
  G4Material *Air=man->FindOrBuildMaterial("G4_AIR");
  

  // Default materials in setup.
  fpWaterMaterial = H2O;
  fpAirMaterial = Air;

}

void DetectorConstruction::SetCellSimulationParameter(G4double defaultCut, G4String targetOrganelle)
{
     
     stringstream sOrganelle(targetOrganelle);
     G4String targetOrg;
     while (sOrganelle>>targetOrg)
     {
       targetOrganelleVector.push_back(targetOrg);// store the target organelle in to a vector
     }
     
    G4ProductionCuts* cuts = new G4ProductionCuts();
    G4double defCut = defaultCut*nanometer;// set up default cut value here
    cuts->SetProductionCut(defCut,"gamma");
    cuts->SetProductionCut(defCut,"e-");
    cuts->SetProductionCut(defCut,"e+");
    cuts->SetProductionCut(defCut,"proton");
    cuts->SetProductionCut(defCut,"alpha");
    
    fpRegion->SetProductionCuts(cuts);
    
}

void DetectorConstruction::SetCellWorld(G4double X, G4double Y, G4double Z)
{
    
    G4double worldSizeX=X;
    G4double worldSizeY=Y;
    G4double worldSizeZ=Z;
    G4VSolid* solidWorld = new G4Box("World",         //its name
    worldSizeX/2,
    worldSizeY/2,
    worldSizeZ/2);  //its size


   G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,  //its solid
                                    fpWaterMaterial,  //its material
                                    "World");    //its name
   physicalWorld = new G4PVPlacement(0,      //no rotation
                                  G4ThreeVector(),  //at (0,0,0)
                                  "World",    //its name
                                  logicWorld,    //its logical volume
                                  0,      //its mother  volume
                                  false,      //no boolean operation
                                  0);      //copy number
                                    // Visualization attributes
   G4VisAttributes* worldVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
   worldVisAtt->SetVisibility(true);
   logicWorld->SetVisAttributes(worldVisAtt);
   
}

void DetectorConstruction::CreateCell(G4String cellType, G4String cellOrganelle, G4String cellShape,G4String organelleColor)
{
  Cell myCell;
  myCell.CellConstruct(cellType,cellOrganelle,cellShape,organelleColor); // construct a cell
  if (generalCellMap.find(myCell.GetCellType())==generalCellMap.end())
  {
    generalCellMap[myCell.GetCellType()]=myCell;// store this type of cell into generalCellMap according to its cell type
  }
  else
  {
    G4cout<<"FOUND DUPLICATE CELL TYPE"<<myCell.GetCellShape()<<G4endl;
    exit(1);
  }
}


G4AssemblyVolume* DetectorConstruction::CreateCellAssembly(G4String cellType,G4double cellSize, G4double nucleusSize)
{
  G4AssemblyVolume* cellAssembly= new G4AssemblyVolume();//declare a cell assembly
  G4Transform3D Tr;
  std::map<G4String,Cell>::iterator mitr=generalCellMap.find(cellType);// find general cell according to cell type
  
  // creating a cell assembly according to cell shape and cell organell
  if (mitr->second.GetCellShape()=="Sphere")
  {
    std::vector<G4String> & cellOrganelleRef=mitr->second.GetCellOrganelle();
    if (cellOrganelleRef.size())//if it is empty, it means no cell organell
    {

      if (std::find(cellOrganelleRef.begin(),cellOrganelleRef.end(),"Cytoplasma")!=cellOrganelleRef.end())
      {
	G4Sphere* cytoplasma = new G4Sphere ("Cytoplasma",nucleusSize,cellSize,0,360,0,180);// define cytoplasma using sphere
	G4LogicalVolume* logicalCytoplasma = new G4LogicalVolume(cytoplasma,  //cell cytoplasma
                                                   fpWaterMaterial,  // material is water
                                                   "Cytoplasma");    // name of the volume
	cellAssembly->AddPlacedVolume(logicalCytoplasma,Tr); // add one cytoplasma volume to assembly
      }
      if (std::find(cellOrganelleRef.begin(),cellOrganelleRef.end(),"Nucleus")!=cellOrganelleRef.end())
      {
	G4Sphere* nucleus = new G4Sphere ("Nucleus",0,nucleusSize,0,360,0,180);// define nucleus using sphere
	G4LogicalVolume* logicalNucleus = new G4LogicalVolume (nucleus, // cell nucleus
                                                       fpWaterMaterial,
                                                       "Nucleus");
	cellAssembly->AddPlacedVolume(logicalNucleus,Tr);// add nucleus volume to assembly
      }
      // you can continue to check other cell organelle cell has other organelle, such as Mitochondial,...
    }
    else
    {
	G4cout<<"NO CELLORGANELLE WAS FOUND!"<<G4endl;
	exit(1);
    }
  } 

  else if (mitr->second.GetCellShape()=="Box")
  {
    std::vector<G4String> & cellOrganelleRef=mitr->second.GetCellOrganelle();
    if (cellOrganelleRef.size())//if it is empty, it means no cell organell
    {
       if (std::find(cellOrganelleRef.begin(),cellOrganelleRef.end(),"Cytoplasma")!=cellOrganelleRef.end())
       {
	 G4Box* cytoplasma = new G4Box("Cytoplasma",cellSize/2,cellSize/2,cellSize/2);// here just an example, the size could be set up according to specific cell
	 G4LogicalVolume* logicalCytoplasma = new G4LogicalVolume(cytoplasma,  //cell cytoplasma
                                                   fpWaterMaterial,  //material is water
                                                   "Cytoplasma");    // name of the volume
	 cellAssembly->AddPlacedVolume(logicalCytoplasma,Tr); // add one cytoplasma volume to assembly
       }
       if (std::find(cellOrganelleRef.begin(),cellOrganelleRef.end(),"Nucleus")!=cellOrganelleRef.end())
      {
	G4Sphere* nucleus = new G4Sphere ("Nucleus",0,nucleusSize,0,360,0,180);// define nucleus using sphere
	G4LogicalVolume* logicalNucleus = new G4LogicalVolume (nucleus, // cell nucleus
                                                       fpWaterMaterial,
                                                       "Nucleus");
	cellAssembly->AddPlacedVolume(logicalNucleus,Tr);// add nucleus volume to assembly
      }
      // could add other cell organelle in here
    }
    else
    {
      G4cout<<"NO CELLORGANELLE WAS FOUND!"<<G4endl;
      exit(1);
      
    }
    
  }
   
  else if (mitr->second.GetCellShape()=="Ellipsoid")
  {
    std::vector<G4String> & cellOrganelleRef=mitr->second.GetCellOrganelle();
    if (cellOrganelleRef.size())//if it is empty, it means no cell organell
    {
       if (std::find(cellOrganelleRef.begin(),cellOrganelleRef.end(),"Cytoplasma")!=cellOrganelleRef.end())
       {
	 G4Ellipsoid* cytoplasma = new G4Ellipsoid("Cytoplasma",cellSize/2,cellSize/3,cellSize/4,0,0);// here is also just an example
	 G4LogicalVolume* logicalCytoplasma = new G4LogicalVolume(cytoplasma,  //cell cytoplasma
                                                   fpWaterMaterial,  //material is water
                                                   "Cytoplasma");    // name of the volume
	 cellAssembly->AddPlacedVolume(logicalCytoplasma,Tr); // add one cytoplasma volume to assembly
       }
       if (std::find(cellOrganelleRef.begin(),cellOrganelleRef.end(),"Nucleus")!=cellOrganelleRef.end())
      {
	G4Sphere* nucleus = new G4Sphere ("Nucleus",0,nucleusSize,0,360,0,180);// define nucleus using sphere
	G4LogicalVolume* logicalNucleus = new G4LogicalVolume (nucleus, // cell nucleus
                                                       fpWaterMaterial,
                                                       "Nucleus");
	cellAssembly->AddPlacedVolume(logicalNucleus,Tr);// add nucleus volume to assembly
      }
      // could add other cell organelle in here
    }
    else
    {
      G4cout<<"NO CELLORGANELLE WAS FOUND!"<<G4endl;
      exit(1);
      
    }
  }
    
  else
  {
    G4cout<<"NO CELLSHAPE WAS FOUND!"<<G4endl;
    exit(1);
  }
  
    return cellAssembly;   // retrun the created cell assembly
     
}


void DetectorConstruction::ClearCellStore()// this function is used to clear all the registered cell when the detector gemoetry update
{
    generalCellMap.clear();
    targetOrganelleVector.clear();
    cellMassMap.clear();
}

void  DetectorConstruction::SetCellPlacement(int cellID, G4String cellType, G4ThreeVector position,G4double cellSize, G4double nucleusSize)
{
    stringstream sName;
    sName<<cellID;
    G4RotationMatrix Rm;
    G4Transform3D Tr;
    Tr= G4Transform3D(Rm,position); 
    G4AssemblyVolume* cellAssembly=CreateCellAssembly(cellType,cellSize,nucleusSize);
    cellAssembly->MakeImprint(physicalWorld->GetLogicalVolume(),Tr);// place the cell into tissue
   // cerr<< "The assembly id is :" <<cellAssembly->GetAssemblyID()<<endl;
    
    std::map<G4String,Cell>::iterator mitr=generalCellMap.find(cellType);// find general cell according to cell type
    std::vector<G4String> & cellOrganelleRef=mitr->second.GetCellOrganelle();
    std::vector<G4String> & organelleColorRef=mitr->second.GetCellOrganelleColor();
    
    std::vector<G4String>::iterator vitCyto=std::find(cellOrganelleRef.begin(),cellOrganelleRef.end(),"Cytoplasma");
    std::vector<G4String>::iterator vitNucl=std::find(cellOrganelleRef.begin(),cellOrganelleRef.end(),"Nucleus");
    // other cell organelle could also be added here
    
    if (vitCyto!=cellOrganelleRef.end())// When the cell has cytoplasma, then do the below stuff
    {
      G4Colour colorOfCytoplasma;
      G4String color;
      color=organelleColorRef.at(std::distance(cellOrganelleRef.begin(),vitCyto)); // get color of cytoplasma
      if (color=="Red")
      colorOfCytoplasma=G4Colour::Red();
      else if (color=="Gray")
      colorOfCytoplasma=G4Colour::Gray();
      else if (color=="Grey")
      colorOfCytoplasma=G4Colour::Grey();
      else if (color=="Black")	
      colorOfCytoplasma=G4Colour::Black();
      else if (color=="Green")
      colorOfCytoplasma=G4Colour::Green();
      else if (color=="Blue")
      colorOfCytoplasma=G4Colour::Blue();
      else if (color=="Cyan")
      colorOfCytoplasma=G4Colour::Cyan();
      else if (color=="Magenta")
      colorOfCytoplasma=G4Colour::Magenta();
      else if (color=="Yellow")
      colorOfCytoplasma=G4Colour::Yellow();
      else
      colorOfCytoplasma=G4Colour::Red();//default color of cell is red

      G4VisAttributes* cytoplasmaVisAtt = new G4VisAttributes(colorOfCytoplasma);// set the cell nuclues color 
      cytoplasmaVisAtt->SetVisibility (true);
      std::vector<G4VPhysicalVolume*>::iterator pVItr=cellAssembly->GetVolumesIterator();    
      G4VPhysicalVolume* pVCytoplasma= *(pVItr+std::distance(cellOrganelleRef.begin(),vitCyto));
      
    // cerr<<"The  name of PV is "<<pVCytoplasma->GetName()<<endl;
      
      pVCytoplasma->SetName(sName.str()+"_"+"Cytoplasma");// set the name of physicalVolume of Cytoplasma
      
    //  cerr<<"The  name of PV is: "<<pVCytoplasma->GetName()<<endl;
      G4LogicalVolume* logicalCytoplasma=pVCytoplasma->GetLogicalVolume();
      logicalCytoplasma->SetVisAttributes(cytoplasmaVisAtt);// display cell Nucleus
      if (std::find(targetOrganelleVector.begin(),targetOrganelleVector.end(),"Cytoplasma")!=targetOrganelleVector.end())
      {
	fpRegion->AddRootLogicalVolume(logicalCytoplasma); //add cell cytoplasma into region if it is set up in targetOrganelleVector
      }
      
      G4double massOfCytop=logicalCytoplasma->GetMass(); // get the mass of cytoplasm
   
      cellMassMap[cellID]=cellMassMap[cellID]+ massOfCytop/kilogram; //add cytoplasm mass to cellMassMap
      
      
    
    }
    
    
    if (vitNucl!=cellOrganelleRef.end())// when the cell has nuclues, then do the below stuff
    {
      G4Colour colorOfNucleus;
      G4String color;
      color=organelleColorRef.at(std::distance(cellOrganelleRef.begin(),vitNucl)); // get color of nucleus
      if (color=="Red")
      colorOfNucleus=G4Colour::Red();
      else if (color=="Gray")
	colorOfNucleus=G4Colour::Gray();
      else if (color=="Grey")
	colorOfNucleus=G4Colour::Grey();
      else if (color=="Black")
	colorOfNucleus=G4Colour::Black();
      else if (color=="Green")
	colorOfNucleus=G4Colour::Green();
      else if (color=="Blue")
	colorOfNucleus=G4Colour::Blue();
      else if (color=="Cyan")
	colorOfNucleus=G4Colour::Cyan();
      else if (color=="Magenta")
	colorOfNucleus=G4Colour::Magenta();
      else if (color=="Yellow")
	colorOfNucleus=G4Colour::Yellow();
      else
        colorOfNucleus=G4Colour::Red();//default color of cell is red
	
      G4VisAttributes* nucleusVisAtt = new G4VisAttributes(colorOfNucleus);// set the cell nuclues color 
      nucleusVisAtt->SetVisibility (true);
      
      std::vector<G4VPhysicalVolume*>::iterator pVItr=cellAssembly->GetVolumesIterator();// get iterator of vector storing physical volumes
     
      
      G4VPhysicalVolume* pVNucleus=  *(pVItr+std::distance(cellOrganelleRef.begin(),vitNucl));
    
      pVNucleus->SetName(sName.str()+"_"+"Nucleus");// set the name of physicalVolume of Nuclues
     // cerr<<"The name of PV is :"<<pVNucleus->GetName()<<endl;
      G4LogicalVolume* logicalNucleus=pVNucleus->GetLogicalVolume();
      logicalNucleus->SetVisAttributes(nucleusVisAtt);// display cell Nucleus
      if (std::find(targetOrganelleVector.begin(),targetOrganelleVector.end(),"Nucleus")!=targetOrganelleVector.end())
      {
	fpRegion->AddRootLogicalVolume(logicalNucleus);// add cell nucleus into region
      }
      
      G4double massOfNuc=logicalNucleus->GetMass();
      
      cellMassMap[cellID]=cellMassMap[cellID]+massOfNuc/kilogram; // add nucleus mass to cellMassMap
      nucleusMassMap[cellID]=massOfNuc/kilogram;
       
    }
    // can add the color of other organelle here
    

   
    
}

EdepCellInformation DetectorConstruction::GetEdepCellInformation(G4VPhysicalVolume* physicalVolume)
{
  G4String physicalVolumeName = physicalVolume->GetName();// get the name of this phyiscalVolume
  
  stringstream sName(physicalVolumeName); 
  G4String segment;
  std::vector<std::string> seglist;
  EdepCellInformation cellEdep;
  
  if (physicalVolumeName!="World")// if it is not in World, then this is energy deposition in cells
  {
    while(std::getline(sName, segment, '_'))
    {
      seglist.push_back(segment);
    }
    stringstream idBuffer(seglist[0]);
    int cellID;
    idBuffer>>cellID;
    cellEdep.cellID=cellID;
    cellEdep.affectedCellOrganelle=*(seglist.begin()+1);
  }
  else // if it is in World, then this is energy deposition in ECM
  {
    cellEdep.cellID=-1;
    cellEdep.affectedCellOrganelle="ECM";
  }

  return cellEdep;
  
}

map< int, double > DetectorConstruction::GetCellMassMap() const
{
    return cellMassMap;
}

map< int, double > DetectorConstruction::GetNucleusMassMap() const
{
    return nucleusMassMap;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
