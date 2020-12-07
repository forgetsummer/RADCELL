
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
    double xDim=1;//unit in mm
    double yDim=1;//unit in mm
    double zDim=1;//unit in mm
    double d=0.01;//unit in mm
    double D=1E-5;//unit in mm^2/s 
    double r=0.63E-17; // unit in 1/#.s
    int Rt=5000;// unit in #
    double mu=200; // unit in #
    int cellNum = 100; //unit in number
    double deltaT_diffusion =1; // unit in second
    double deltaT_diffusionUpdate=1;// unit in second    
    
    DiffusionReactionSolver diffusionSolver;
    diffusionSolver.DiffusionReactionInitialization(xDim,yDim,zDim,d,D,r,Rt,deltaT_diffusion);
    diffusionSolver.SetVerbose(0);
    
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d,d,d);
    layout.RectangularSlab(xDim,yDim,zDim,cellNum);
    cout<<"The total cell number is "<<layout.GetCellNumber()<<endl;
    
    for(int i=0;i<cellNum;i++)
    {
        double cX = layout.GetCellPositionX(i)+xDim/2.0; 
        double cY = layout.GetCellPositionY(i)+yDim/2.0;
        double cZ = layout.GetCellPositionZ(i)+zDim/2.0;
        if (i%2==0)
        {
            int cellState=1;// initial cell state is S1, which means all the cells are in healthy state 
            diffusionSolver.CellStateUpdate(i,cX,cY,cZ,cellState,mu);// get the initial state information for diffusion solver
        }
        else
        {
            int cellState=2;// initial cell state is S1, which means all the cells are in healthy state 
            diffusionSolver.CellStateUpdate(i,cX,cY,cZ,cellState,mu);// get the initial state information for diffusion solver
        }

    }
    
    for (int i=0;i<1000;i++)
    {
        diffusionSolver.DiffusionReactionCalculation(deltaT_diffusionUpdate,1);// diffusion solver progresses
        cout<<"The concentration at point "<<0.5<<" "<<0.6<<" "<<0<<" is "<<diffusionSolver.GetConcentration(0.5,0.6,0)<<endl;
        string path;
        path = "/home/ruirui/RAD-build/concentration/";
        diffusionSolver.WriteConcentrationToFile(path,i); 
    }
   

    return 0;
}

