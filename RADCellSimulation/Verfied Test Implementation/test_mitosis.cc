
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
    CellStateModel myCellState;
    
//     double mu = 1.5;
//     double sigma = 0.25;
//     ofstream file1;
//     file1.open("testGaussian.csv");
//     
//     for (int i =0;i<10000;i++)
//     {
//         double t;
//         t = myCellState.GaussianSampling(mu,sigma);
// //         cout<<t<<","<<0<<endl;
//         file1<<t<<","<<0<<endl;
//     }
//     
    
    myCellState.CellInformationSetup(1,1.5,0.25,6.0,0.25,1.5,0.25,1.0,0.25,100);
    int cellNum = 10;
    for (int i =0;i<cellNum;i++) // seed 1000 cells with random cell phase
    {
        myCellState.CellPhaseInitializationRandom(i,1);
//         myCellState.CellPhaseInitialization(i,1,"G0");
    }
//     myCellState.CellPhaseInitialization(0,1,"M");
//     myCellState.CellPhaseInitialization(1,1,"G1"); 
//     myCellState.CellPhaseInitialization(2,1,"G0");
//     myCellState.CellPhaseInitialization(3,1,"G0");
//     myCellState.CellPhaseInitialization(4,1,"S");
//     myCellState.CellPhaseInitialization(5,1,"S");
//     myCellState.CellPhaseInitialization(6,1,"M");
//     myCellState.CellPhaseInitialization(7,1,"M");
//     myCellState.CellPhaseInitialization(8,1,"G1");
//     myCellState.CellPhaseInitialization(9,1,"M");
    
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
    
    double deltaT = 60;// time step is 1 min
    int N = 100000;
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
            myCellState.CellPhaseUpdate(mitr_cell->first,true,1,deltaT,1);
        }

    }

    return 0;
}

