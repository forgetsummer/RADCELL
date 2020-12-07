
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
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProductionCuts.hh"
#include "Randomize.hh"  // uisng random number
#include <math.h> //uisng math function
#include <vector>
#include "G4RegionStore.hh"
#include "CellLayoutInitializer.hh"

int main(int argc,char** argv)
{
    G4int N=300;
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(0.05,0.05,0.05);// cell home size is in unit of mm
    layout.RectangularSlab(1, 1, 1,N);// unit here is also mm
    RADCellSimulation mySim;
    mySim.SetCellWorld(1.1,1.1,1.1);//create cell world, dimension in unit of mm
    mySim.CreateCell("Epithelia","Cytoplasma Nucleus","Sphere","Red Green");
    mySim.SetCellSimulationParameter(1," Nucleus");
    
    
    for (G4int i=0; i<layout.GetCellNumber(); i++)
    {

        mySim.CellConstruction(i,"Epithelia",layout.GetCellPositionX(i),layout.GetCellPositionY(i),layout.GetCellPositionZ(i),20,5);
    }
    cout<<"The registered region in first run is "<<G4RegionStore::GetInstance()->GetRegion("Target")->GetNumberOfRootVolumes()<<endl;   
    cout<<"The number of seeded cells is "<<layout.GetCellNumber()<<endl;
    mySim.RADCellSimulationInitialize(argc, argv);
    
//     mySim.EnergyDistributionCalculation( "out","microdosimetry.in",true);
    mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
//     mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
//     mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
    
    
    cout<<"run second time"<<endl;
    CellLayoutInitializer layout1;
    G4int N1=320;
    layout1.SetCellHomeParamter(0.05,0.05,0.05);// cell home size is in unit of mm
    layout1.Blob(1,N1);// unit here is also mm
    mySim.UpdateGeometryInitialize();
    cout<<"after UpdateGeometryInitialize"<<endl;
    mySim.SetCellWorld(1.1,1.1,1.1);//create cell world, dimension in unit of mm
    mySim.CreateCell("Epithelia","Cytoplasma Nucleus","Sphere","Blue Green");
    mySim.SetCellSimulationParameter(1," Nucleus");
    
    for (G4int i=0; i<layout1.GetCellNumber(); i++)
    {
        mySim.CellConstruction(i,"Epithelia",layout1.GetCellPositionX(i),layout1.GetCellPositionY(i),layout1.GetCellPositionZ(i),20,5);
    }
    cout<<"The registered region in second run is "<<G4RegionStore::GetInstance()->GetRegion("Target")->GetNumberOfRootVolumes()<<endl;
    cout<<"The number of seeded cells is "<<layout1.GetCellNumber()<<endl;
    cout<<"after rebuild all cells"<<endl;
    mySim.UpdateGeometryFinalize();
//     mySim.EnergyDistributionCalculation( "out","microdosimetry.in",true);
    mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
//     mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
//     mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
//     
    cout<<"run third time"<<endl;
    CellLayoutInitializer layout2;
    G4int N2=30;
    layout2.SetCellHomeParamter(0.05,0.05,0.05);// cell home size is in unit of mm
    layout2.RectangularSlab(1, 1, 0,N2);// unit here is also mm
    mySim.UpdateGeometryInitialize();
    cout<<"after UpdateGeometryInitialize"<<endl;

    mySim.SetCellWorld(1.1,1.1,1.1);//create cell world, dimension in unit of mm
    mySim.CreateCell("Epithelia","Cytoplasma Nucleus","Sphere","Blue Green");
    mySim.SetCellSimulationParameter(1," Nucleus");
    
    for (G4int i=0; i<layout2.GetCellNumber(); i++)
    {
        mySim.CellConstruction(i,"Epithelia",layout2.GetCellPositionX(i),layout2.GetCellPositionY(i),layout2.GetCellPositionZ(i),20,5);
    }
    
    cout<<"The number of seeded cells is "<<layout2.GetCellNumber()<<endl;
    cout<<"after rebuild all cells"<<endl;
    mySim.UpdateGeometryFinalize();
    cout<<"The registered region in second run is "<<G4RegionStore::GetInstance()->GetRegion("Target")->GetNumberOfRootVolumes()<<endl;
    mySim.EnergyDistributionCalculation( "out","microdosimetry.in",true);
//     mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
//     mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
//     mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
    
    return 0;
}

