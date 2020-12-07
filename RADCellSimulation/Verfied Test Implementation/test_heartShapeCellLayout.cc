
// This program is largely amended from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// This the test function for the project

#include "RADCellSimulation.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include "G4RegionStore.hh"
#include "CellLayoutInitializer.hh"
#include "CellDoseAnalysis.hh"
#include "CellDNADamageAnalysis.hh"
#include <fstream>
#include <time.h>
#include "ReadRadiationTransportInfo.hh"

int main(int argc,char** argv)
{
    int N=1500; // this line for seeding 300 cell
    
    double dim = 5;
    double xDim = dim;
    double yDim = dim;
    double zDim = dim;
    
        
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(0.01,0.01,0.01);// cell home size is in unit of mm
//         layout.RectangularSlab(0.05, 0.05, 0.05,N);// unit here is also mm, this line is for single cell simulation
   // layout.RectangularSlab(1,1,0,N); // this line is for multiple cell simulation
    
    layout.Heart2D(0.8*dim,N);
//     layout.Heart(0.8*dim,N);
//     exit(0);
    RADCellSimulation mySim;
//     mySim.SetCellWorld(0.1,0.1,0.1);//create cell world, dimension in unit of mm, this line is for testing single cell simulation

    mySim.SetCellWorld(1.1*xDim,1.1*yDim,1.1*zDim); // this line is for testing multiple cell simulation
    mySim.CreateCell("Epithelia","Cytoplasma Nucleus","Sphere","Red Green");
    mySim.SetCellSimulationParameter(1,"Nucleus");
    for (int i=0; i<layout.GetCellNumber(); i++)
    {
        cout<<"X="<<layout.GetCellPositionX(i)<<" Y="<<layout.GetCellPositionY(i)<<" Z="<<layout.GetCellPositionZ(i)<<endl;

        mySim.CellConstruction(i,"Epithelia",layout.GetCellPositionX(i),layout.GetCellPositionY(i),layout.GetCellPositionZ(i),20,5);
    }

    mySim.RADCellSimulationInitialize(argc, argv);

    mySim.EnergyDistributionCalculation("gui","Mono_Electron.in",false);
    mySim.EnergyDistributionCalculation("gui","Proton-1MeV.in",false);

    return 0;
}

