
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
    double xDim=5;//unit in mm
    double yDim=5;//unit in mm
    double zDim=0;//unit in mm
    double d=0.05;//unit in mm
    double D=1E-5;//unit in mm^2/s 
    double r=0.63E-17; // unit in 1/#.s
    int Rt=5000;// unit in #
    double mu=200; // unit in #
    double deltaT_diffusion =1; // unit in second
    double deltaT_diffusionUpdate=60;// unit in second
    double deltaT_cellPhaseUpdate =60;//unit in second 
    double deltaT_cellStateUpdate =60;// unit in second
    int cellType=1;
//     double tG1 = 1.5;// unit in hour
//     double sigmaG1=0.25;//unit in hour
//     double tS=6;
//     double sigmaS=0.25;
//     double tG2=1.5;
//     double sigmaG2=0.25;
//     double tM=1;
//     double sigmaM=0.25;
    double tG1 = 9;// unit in hour
    double sigmaG1=1.8;//unit in hour
    double tS=11;
    double sigmaS=2.2;
    double tG2=1;
    double sigmaG2=0.2;
    double tM=1;
    double sigmaM=0.2;
    double alpha=0.2497; // unit in 1/DSB
    double beta=25.71; // unit in ml/pg
    double E1=0;
    double E2=18.15168;
    double E3=48.47;
    double sigma = 6.96;
    double fG1 = 1;
    double fG2 = 1.015306122;
    double fM = 1.015306122;
    double fS = 0.816326;
    double f1 = 0.62;
    double f2 = 0.38;
    double lambda1=3.31; // unit in 1/hour
    double lambda2=0.14; // unit in 1/hour
    double T=10; // time step for one MCS, unit in second
    int totalTimeStepNum=60480;
    int cellNum = 1000; //unit in number
    
//*************************************************************************************************
//**********************Above is the simulation parameters
//*************************************************************************************************


    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d,d,d);
    layout.RectangularSlab(xDim,yDim,zDim,cellNum);
    cout<<"The total cell number is "<<layout.GetCellNumber()<<endl;

    RADCellSimulation mySim;
    double xWorld = 1.1*xDim;
    double yWorld = 1.1*yDim;
    double zWorld = 1.1*xDim;
    mySim.SetCellWorld(xWorld,yWorld,zWorld); // this line is for testing multiple cell simulation
    mySim.CreateCell("Epithelia","Cytoplasma Nucleus","Sphere","Blue Green");// create biological cells
    mySim.SetCellSimulationParameter(1,"Nucleus");
    
    for (int i=0; i<layout.GetCellNumber(); i++)
    {
        cout<<layout.GetCellPositionX(i)<<","<<layout.GetCellPositionY(i)<<","<<layout.GetCellPositionZ(i)<<endl;
        mySim.CellConstruction(i,"Epithelia",layout.GetCellPositionX(i),layout.GetCellPositionY(i),layout.GetCellPositionZ(i),20,5);
    }

    mySim.RADCellSimulationInitialize(argc, argv);// initialize the radiation transport module
    mySim.EnergyDistributionCalculation("gui","microdosimetry.in",false); // use this to check cell layout condition
    
    double clock_cellPhaseUpdate=0;
    double clock_cellStateUpdate=0;
    double clock_diffusionUpdate=0;
    int cycle_cellPhaseUpdate=0;
    int cycle_cellStateUPdate=0;
    int cycle_diffusionUpdate=0;
    
    CellStateModel myCellState;
    myCellState.CellInformationSetup(cellType,tG1,sigmaG1,tS,sigmaS,tG2,sigmaG2,tM,sigmaM,100);
    myCellState.CellStateModelParameterSetup(alpha,beta,E1,E2,E3,sigma,fG1,fS,fG2,fM,f1,lambda1,f2, lambda2);
    myCellState.TissueGeometryInitialization(xDim,yDim,zDim,d);
   
    
    for (int i =0;i<cellNum;i++)//initialize the cell position 
    {
        double cX = layout.GetCellPositionX(i)+xDim/2.0; 
        double cY = layout.GetCellPositionY(i)+yDim/2.0;
        double cZ = layout.GetCellPositionZ(i)+zDim/2.0;
        myCellState.CellPositionInitialization(i,cX,cY,cZ);
    }

    for (int i =0;i<cellNum;i++) // seed cellNum cells with random cell phase
    {
        myCellState.CellPhaseInitializationRandom(i,1);
//         myCellState.CellPhaseInitialization(i,1,"G1");
    }
    
    for (int i =0;i<cellNum;i++)
    {
        myCellState.CellStateInitialization(i,"S1");
    }

//     myCellState.SetUpContactInhibition(false);
    ofstream file;
    ofstream file1;
    file.open("cellPhase.csv");
    file1.open("cellState.csv");

//*****************************************************************************************
//****************Above is cell system CellPhaseInitialization
//****************************************************************************************
    
    for(int i=0;i<totalTimeStepNum;i++)// Loop all the time steps, absolute time 
    {
        clock_cellPhaseUpdate = clock_cellPhaseUpdate + T;
        clock_cellStateUpdate = clock_cellStateUpdate + T;
        clock_diffusionUpdate = clock_diffusionUpdate + T; 

        std::map<int, std::string> phaseMap;// Here we use phase map to retract the dynamic cell index with respect to time
        phaseMap = myCellState.GetCellPhase();// update cellPhaseMap
        if (clock_cellPhaseUpdate>=deltaT_cellPhaseUpdate) // cell phase update
        {
            cycle_cellPhaseUpdate=cycle_cellPhaseUpdate+1;
            cout<<"PROCESSING STEP FOR UPDATING CELL PHASE: "<<cycle_cellPhaseUpdate<<endl;
            std::map<int, double> ageMap;
            std::map<int, double> durationMap;
            ageMap = myCellState.GetCellAge();
            durationMap = myCellState.GetCellPhaseDuration();
//             cout<<"Cell Information for Step: "<<i<<endl;
//             for (std::map<int,std::string>::iterator mitr_cell=phaseMap.begin();mitr_cell!=phaseMap.end();mitr_cell++)
//             {
//                 cout<<"Cell: "<<mitr_cell->first<<" Phase: "<<mitr_cell->second<<" Age: "<<ageMap[mitr_cell->first]<<" Duration: "<<durationMap[mitr_cell->first]<<endl;
//             }
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
            file<<cycle_cellPhaseUpdate<<","<<G1Num<<","<<SNum<<","<<G2Num<<","<<mitosisNum<<","<< phaseMap.size()<<endl;

            for (std::map<int, std::string>::iterator mitr_cell=phaseMap.begin();mitr_cell!=phaseMap.end();mitr_cell++)
            {            
                myCellState.CellPhaseUpdate(mitr_cell->first,1,true,deltaT_cellPhaseUpdate,1);
            }
            clock_cellPhaseUpdate=0;   
//             cout<<"The total cell number now is  "<<phaseMap.size()<<endl;
        }
        
  
    }
    cout<<"finish"<<endl;

    return 0;
}

