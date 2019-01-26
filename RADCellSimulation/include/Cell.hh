// 
// This code is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015
#ifndef Cell_h
#define Cell_h 1

#include "G4String.hh"
#include <vector>
class Cell // creating a class named cell
{
public:
//   Cell(); // constructor
  void CellConstruct(G4String cType, G4String cOrganelle, G4String cShape,G4String cColor);
  G4String GetCellType(){return cellType;}
  G4String GetCellShape(){return cellShape;}
  
  std::vector<G4String> & GetCellOrganelle () {return cellOrganelle;}
  std::vector<G4String> & GetCellOrganelleColor() {return colorOfOrganelle;}

private:
  G4String cellType;// cell type
  std::vector<G4String> cellOrganelle; // a vector storing cell organelle
  std::vector<G4String> colorOfOrganelle;// a vector storing cell
  G4String cellShape; // cell shape
};
#endif
