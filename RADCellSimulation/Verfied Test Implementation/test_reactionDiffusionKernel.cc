
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

#include <sstream>


int main(int argc,char** argv)
{
    ReactionDiffusionSimulation react;
    double xDim=M_PI; // unit in mm
    double yDim=M_PI;
    double zDim=0;
    double d = 0.01;//unit in mm
    double D = 1E-6; // D in unit of mm^2/s
    double r = 0.63E-17;
    int Rt = 5000;
    double mu = 200;
    double deltaT = 1;
    double T = 100;
    react.MeshGeometry(xDim, yDim, zDim, d );           
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d,d,d);
    int cellNum =100;
    layout.RectangularSlab(xDim,yDim,zDim,cellNum);
    
    std::map<int, double>cX;
    std::map<int, double>cY;
    std::map<int, double>cZ;
    std::map<int, int> cellState;
    std::map<int,double>  secretionRate;

    for (int i=0;i<cellNum;i++)
    {
        cX[i] = layout.GetCellPositionX(i)+xDim/2.0;
        cY[i] = layout.GetCellPositionY(i)+yDim/2.0;
        cZ[i] = layout.GetCellPositionZ(i)+zDim/2.0;
        secretionRate[i] = mu;
//         cellState[i]=2;
        if (i%2==0)
        {
            cellState[i]=1;
        }
        else{
            cellState[i]=2;
        }
    }
    
//     for (int i=0;i<50;i++)
//     {
//         cout<<"The cell location is "<<cX[i]<<" "<<cY[i]<<" "<<cZ[i]<<endl; 
//     }
//     
    react.CellStateInitialization(cX,cY,cZ,cellState,secretionRate);
    
    std::vector<std::vector<std::vector<double> > > initialC;
    
    int N_X=react.GetMeshDIMX();
    int N_Y=react.GetMeshDIMY();
    int N_Z=react.GetMeshDIMZ();
    cout<<"N_X N_Y N_Z="<< N_X <<" "<<N_Y<<" "<<N_Z<<endl;

    for (int i=0;i<N_X;i++ )
    {
        initialC.push_back(std::vector<std::vector<double> >());
        for (int j=0;j<N_Y;j++)
        {
            initialC[i].push_back(std::vector<double>());
            for (int k=0;k<N_Z;k++)
            {
                initialC[i][j].push_back(0); // initial concentration
            }
        }
    }

    react.DiffusionReactionCalculation(initialC,D,r,Rt,deltaT,T);
    react.WriteConcentrationProfile("/home/ruirui/RAD-build/concentration/");

    return 0;
}

