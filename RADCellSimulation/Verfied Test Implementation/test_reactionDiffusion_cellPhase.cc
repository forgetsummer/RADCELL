
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
#include "CellStateModel.hh"

#include <sstream>


int main(int argc,char** argv)
{
    double xDim=1;
    double yDim=1;
    double zDim=0;
    double d=0.01;
    int cellNum = 10;
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d,d,d);
    layout.RectangularSlab(xDim,yDim,zDim,cellNum);
    DiffusionReactionSolver diffusionSolver;
    double D=1E-6;
    double r=0.63E-17;
    int Rt=5000;
    double mu=200;
    double deltaT_diffusion=1;
    double T=1;
    int cellType=1;
    double tG1;
    double sigmaG1;
    double tS;
    double sigmaS;
    double tG2;
    double sigmaG2;
    double tM;
    double sigmaM;
    diffusionSolver.DiffusionReactionInitialization(xDim,yDim,zDim,d,D,r,Rt,mu,deltaT_diffusion,T);
    for(int i=0;i<cellNum;i++)
    {
        double cX = layout.GetCellPositionX(i)+xDim/2.0;
        double cY = layout.GetCellPositionY(i)+yDim/2.0;
        double cZ = layout.GetCellPositionZ(i)+zDim/2.0;
        int cellState=1;// initial cell state is S1, which means all the cells are in healthy state 
        diffusionSolver.CellStateUpdate(i,cX,cY,cZ,cellState);// get the initial state information for diffusion solver
    }
    
    CellStateModel myCellState;   
    myCellState.CellInformationSetup(cellType,tG1,sigmaG1,tS,sigmaS,tG2,sigmaG2,tM,sigmaM,Rt);
    for (int i =0;i<cellNum;i++) // seed cellNum cells with random cell phase
    {
        myCellState.CellPhaseInitializationRandom(i,cellType);
    }
    
//     for(int i=0;i<10000;i++)// Loop all the time steps, absolute time 
//     {
//         diffusionSolver.DiffusionReactionCalculation();
//         cout<<"The concentration at point "<<0.5<<" "<<0.6<<" "<<0<<" is "<<diffusionSolver.GetConcentration(0.5,0.6,0)<<endl;
//         string path;
//         path = "/home/ruirui/RAD-build/concentration/";
//         diffusionSolver.WriteConcentrationToFile(path,i);
//     }
//     

    cout<<"THE CELL PHASE INITIALLY IS "<<endl;
    std::map<int, std::string> phaseMap;
    std::map<int, double>durationMap;
    phaseMap = myCellState.GetCellPhase();
    durationMap = myCellState.GetCellPhaseDuration();
    int mitosisNum = 0;
    for (std::map<int, std::string>::iterator mitr_cell=phaseMap.begin();mitr_cell!=phaseMap.end();mitr_cell++)
    {
        cout<<"THE CELL PHASE AND DURATION OF CELL "<<mitr_cell->first<<" : "<<mitr_cell->second<<" "<<durationMap[mitr_cell->first] <<endl;

    }
    
    double deltaT = 0.01;
    int N = 1200;
    ofstream file;
    file.open("cellPhase.csv");
    
    for (int i=0;i<N;i++)
    {
        cout<<"PROCESSING STEP:"<<i<<endl;
        std::map<int, std::string> phaseMap;
        std::map<int, double>durationMap;
        std::map<int,double> ageMap;
        phaseMap = myCellState.GetCellPhase();
        durationMap = myCellState.GetCellPhaseDuration();
        ageMap = myCellState.GetCellAge();
//         cout<<"Cell Information for Step: "<<i<<endl;
// 
//         for (std::map<int,std::string>::iterator mitr_cell=phaseMap.begin();mitr_cell!=phaseMap.end();mitr_cell++)
//         {
//             cout<<"Cell: "<<mitr_cell->first<<" Phase: "<<mitr_cell->second<<" Age: "<<ageMap[mitr_cell->first]<<" Duration: "<<durationMap[mitr_cell->first]<<endl;
//         }
        int mitosisNum = 0;
        int G1Num=0;
        int G2Num=0;
        int SNum=0;
        for (std::map<int, std::string>::iterator mitr_cell=phaseMap.begin();mitr_cell!=phaseMap.end();mitr_cell++)
        {
            if (mitr_cell->second == "M")
            {
                mitosisNum = mitosisNum+1;
            }
            if (mitr_cell->second == "G1")
            {
                G1Num = G1Num+1;
            }
            if (mitr_cell->second == "G2")
            {
                G2Num = G2Num+1;
            }
            if (mitr_cell->second == "S")
            {
                SNum = SNum+1;
            }
        }
        file<<i<<","<<G1Num<<","<<SNum<<","<<G2Num<<","<<mitosisNum<<","<< phaseMap.size()<<endl;
        
        for (std::map<int, std::string>::iterator mitr_cell=phaseMap.begin();mitr_cell!=phaseMap.end();mitr_cell++)
        {            
            myCellState.CellPhaseUpdate(mitr_cell->first,cellType,true,deltaT,1);
        }

    }

    return 0;
}

