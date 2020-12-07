
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
#include "ReactionDiffusionSimulation.hh"
#include "DiffusionReactionSolver.hh"
using namespace std;

int main(int argc,char** argv)
{
    double xDim=3.1415926;//unit in mm
    double yDim=3.1415926;//unit in mm
    double zDim=0;//unit in mm
    double d=0.01;//unit in mm
    double D=1E-5;//unit in mm^2/s 
//     double r=0; // unit in 1/#.s
//     int Rt=0;// unit in #
    double r=0.63E-17; // unit in 1/#.s
    int Rt=5000;// unit in #
    double mu=200; // unit in #
    int cellNum = 100000; //unit in number
    double deltaT_diffusion =0.01; // unit in second
    double deltaT_diffusionUpdate=1;// unit in second    
    
    DiffusionReactionSolver diffusionSolver;
    diffusionSolver.DiffusionReactionInitialization(xDim,yDim,zDim,d,D,r,Rt,deltaT_diffusion);
    diffusionSolver.SetVerbose(1);
    
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d,d,d);
    layout.RectangularSlab(xDim,yDim,zDim,cellNum);
    cout<<"The total cell number is "<<layout.GetCellNumber()<<endl;
    
//     double cX = 0+xDim/2.0; 
//     double cY = 0+yDim/2.0;
//     double cZ = 0+zDim/2.0;
//     int cellState = 2;
//     diffusionSolver.CellStateUpdate(1,cX,cY,cZ,cellState);
    
//     for(int i=0;i<cellNum;i++)
//     {
//         double cX = layout.GetCellPositionX(i)+xDim/2.0; 
//         double cY = layout.GetCellPositionY(i)+yDim/2.0;
//         double cZ = layout.GetCellPositionZ(i)+zDim/2.0;
//         if (i%2==0)
//         {
//             int cellState=1;// initial cell state is S1, which means all the cells are in healthy state 
//             diffusionSolver.CellStateUpdate(i,cX,cY,cZ,cellState);// get the initial state information for diffusion solver
//         }
//         else
//         {
//             int cellState=2;// initial cell state is S1, which means all the cells are in healthy state 
//             diffusionSolver.CellStateUpdate(i,cX,cY,cZ,cellState);// get the initial state information for diffusion solver
//         }
// 
//     }
    
    for (int n=0;n<500;n++)// loop time steps, t =[0,500] seconds
    {
        diffusionSolver.CellStateUpdateInitialize();
        for(int i =0;i<layout.GetCellNumber();i++)
        {
            double cX = layout.GetCellPositionX(i)+xDim/2.0; 
            double cY = layout.GetCellPositionY(i)+yDim/2.0;
            double cZ = layout.GetCellPositionZ(i)+zDim/2.0;
            int cellState = 2;
//             cout<<"cX="<<cX<<" "<<"cY="<<cY<<" "<<"cZ="<<cZ<<endl;
//             double fmu = cos(cX)*cos(cY)*(1+2*D*n); // this is for manufactured method which is used to verify the CAM method
            double fmu = cos(cX)*cos(cY)*(1+2*D*n)+r*Rt*cos(cX)*cos(cY)*n; // this is for manufactured method which is used to verify the CAM method
//             cout<<"cos(cX)="<<cos(cX)<<" "<<"cos(cY)="<<cos(cY)<<endl;
//             cout<<"fmu="<<fmu<<endl;
            diffusionSolver.CellStateUpdate(i,cX,cY,cZ,2,fmu);// get the initial state information for diffusion solver
        }
        
        diffusionSolver.DiffusionReactionCalculation(deltaT_diffusionUpdate,1);// diffusion solver progresses
        cout<<"The concentration at point "<<0.5<<" "<<0.6<<" "<<0<<" is "<<diffusionSolver.GetConcentration(0.5,0.6,0)<<endl;
        string path;
        path = "/home/ruirui/RAD-build-DiffusionVerification/concentration/";
        diffusionSolver.WriteConcentrationToFile(path,n); 
    }
    
    return 0;
}

