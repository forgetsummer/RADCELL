#ifndef CellLayoutInitializer_h
#define CellLayoutInitializer_h 1
#include <vector>

class Index
{
  public:
    int i;
    int j;
    int k;
};
class CellLayoutInitializer
{
  public:
    CellLayoutInitializer();// constructor, maybe not needed somehow
    ~CellLayoutInitializer();
    void RectangularSlab(double xDim, double yDim, double zDim,int seededCellNumber);// function for seeding certain number of cells into slab
    void SetCellHomeParamter(double sizeX, double sizeY, double sizeZ);// function setting up the geometry parameter of cell home containing a cell
    double GetCellPositionX(int i){return cellPositionX[i];}
    double GetCellPositionY(int i){return cellPositionY[i];}
    double GetCellPositionZ(int i){return cellPositionZ[i];}
    int GetCellNumber(){return cellPositionX.size();} // function of obtaining the total seeded cells
    
    
    void Blob(double dDim, int seededCellNumber); // function for seeding certain number of cells in to a blob, dDim is diameter
    
    void Blob2D (double dDim, int seededCellNumber);// dDim is diameter
    
    void Heart(double dDim, int seededCellNumber);// function for seeding certain number of cells in to a heart shape object, purely for fun test, dDim is diameter
    
    void Heart2D (double dDim, int seededCellNumber);// dDim is diameter
  private:
     std::vector<double>cellPositionX; // define a vector to store the cell position
     std::vector<double>cellPositionY;
     std::vector<double>cellPositionZ;
     std::vector<Index> positionIndex;
    

    double cellHomeSizeX;
    double cellHomeSizeY;
    double cellHomeSizeZ; 
  
};
#endif
