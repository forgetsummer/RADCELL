// 
// This is code is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 

#include "Cell.hh"
#include <sstream>

using namespace std;


void Cell::CellConstruct(G4String cType, G4String cOrganelle, G4String cShape,G4String cColor)
{
  cellType=cType;
  cellShape=cShape;
  stringstream organList(cOrganelle);
  stringstream colorList(cColor);
  G4String organ;
  G4String color;
  while(organList>>organ)
  {
    cellOrganelle.push_back(organ);
  }
  while (colorList>>color)
  {
    colorOfOrganelle.push_back(color);
  }
  
}