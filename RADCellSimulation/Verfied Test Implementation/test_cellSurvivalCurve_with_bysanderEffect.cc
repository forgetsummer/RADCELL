
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
#include "math.h"
#include <sstream>
#include <stdlib.h>


int main(int argc,char** argv)
{
    double xDim=1;//unit in mm
    double yDim=1;//unit in mm
    double zDim=0;//unit in mm
    double d=0.01;//unit in mm
    double D=1E-5;//unit in mm^2/s 
    double r=0.63E-17; // unit in 1/#.s
    int Rt=5000;// unit in #
    double mu=200; // unit in #
    double deltaT_diffusion =1; // unit in second
    double deltaT_diffusionUpdate=60;// unit in second
    double deltaT_cellPhaseUpdate =60;//unit in second 
    double deltaT_cellStateUpdate =60;// unit in second
    int cellType=1;
    double tG1 = 1.5;// unit in hour
    double sigmaG1=0.25;//unit in hour
    double tS=6;
    double sigmaS=0.25;
    double tG2=1.5;
    double sigmaG2=0.25;
    double tM=1;
    double sigmaM=0.25;
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
    int totalTimeStepNum=10000;
    int cellNum = 1000; //unit in number
        

        
//*************************************************************************************************
//**********************Above is the simulation parameters
//*************************************************************************************************


    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d,d,d);
    layout.RectangularSlab(xDim,yDim,zDim,cellNum);
    cout<<"The total cell number is "<<layout.GetCellNumber()<<endl;
    
    ofstream file;
    ofstream file1;
    file.open("cellPhase.csv");
    file1.open("cellState.csv");
    int doseCaseNum = 20;
    int maxDose = 8;
    double step_dose = (double)maxDose/doseCaseNum;
    cout<<step_dose<<endl;
//     exit(0);
    for (int n=0;n<=doseCaseNum;n++)
    {
        double clock_cellPhaseUpdate=0;
        double clock_cellStateUpdate=0;
        double clock_diffusionUpdate=0;
        int cycle_cellPhaseUpdate=0;
        int cycle_cellStateUPdate=0;
        int cycle_diffusionUpdate=0;
        cout<<"PROCESSING DOSE STEP: "<<n<<endl;
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
//             myCellState.CellPhaseInitialization(i,1,"G1");
        }
        
        for (int i =0;i<cellNum;i++)
        {
            myCellState.CellStateInitialization(i,"S1");
        }
        
        DiffusionReactionSolver diffusionSolver;
        diffusionSolver.DiffusionReactionInitialization(xDim,yDim,zDim,d,D,r,Rt,deltaT_diffusion);
        diffusionSolver.SetVerbose(0);
        
        for(int i=0;i<cellNum;i++)
        {
            double cX = layout.GetCellPositionX(i)+xDim/2.0; 
            double cY = layout.GetCellPositionY(i)+yDim/2.0;
            double cZ = layout.GetCellPositionZ(i)+zDim/2.0;
            int cellState=1;// initial cell state is S1, which means all the cells are in healthy state 
            diffusionSolver.CellStateUpdate(i,cX,cY,cZ,cellState,mu);// get the initial state information for diffusion solver
        }
        

//         cout<<"get here"<<endl;
    //*****************************************************************************************
    //****************Above is cell system CellPhaseInitialization
    //****************************************************************************************

        for(int i=0;i<totalTimeStepNum;i++)// Loop all the time steps, absolute time 
        {
            clock_cellStateUpdate = clock_cellStateUpdate + T;
            clock_cellPhaseUpdate = clock_cellPhaseUpdate + T;
            clock_diffusionUpdate = clock_diffusionUpdate + T;
            std::map<int, std::string> phaseMap;// Here we use phase map to retract the dynamic cell index with respect to time
            phaseMap = myCellState.GetCellPhase();// update cellPhaseMap
            if (clock_cellPhaseUpdate>=deltaT_cellPhaseUpdate) // cell phase update
            {
                cycle_cellPhaseUpdate=cycle_cellPhaseUpdate+1;
//                 cout<<"PROCESSING STEP FOR UPDATING CELL PHASE: "<<cycle_cellPhaseUpdate<<endl;
                std::map<int, double> ageMap;
                std::map<int, double> durationMap;
                ageMap = myCellState.GetCellAge();
                durationMap = myCellState.GetCellPhaseDuration();
//                 cout<<"Cell Information for Step: "<<i<<endl;
//                 for (std::map<int,std::string>::iterator mitr_cell=phaseMap.begin();mitr_cell!=phaseMap.end();mitr_cell++)
//                 {
//                     cout<<"Cell: "<<mitr_cell->first<<" Phase: "<<mitr_cell->second<<" Age: "<<ageMap[mitr_cell->first]<<" Duration: "<<durationMap[mitr_cell->first]<<endl;
//                 }

                for (std::map<int, std::string>::iterator mitr_cell=phaseMap.begin();mitr_cell!=phaseMap.end();mitr_cell++)
                {            
                    myCellState.CellPhaseUpdate(mitr_cell->first,true,1,deltaT_cellPhaseUpdate,1);
                }
                clock_cellPhaseUpdate=0;   
                
            }
            
            if(clock_diffusionUpdate>=deltaT_diffusionUpdate) // diffusion update
            {
                diffusionSolver.DiffusionReactionCalculation(deltaT_diffusionUpdate,1);// diffusion solver progresses
//                 cout<<"The concentration at point "<<0.5<<" "<<0.6<<" "<<0<<" is "<<diffusionSolver.GetConcentration(0.5,0.6,0)<<endl;
//                 string path;
//                 path = "/home/ruirui/RAD-build/concentration/";
//                 diffusionSolver.WriteConcentrationToFile(path,i);
                clock_diffusionUpdate = 0;
            }
                
            
            if (clock_cellStateUpdate>=deltaT_cellStateUpdate)// cell state update
            {
                cycle_cellStateUPdate = cycle_cellStateUPdate+1;
                cout<<"PROCESSING STEP FOR UPDATING CELL STATE: "<<cycle_cellStateUPdate<<endl;
                std::map<int, std::string> stateMap;
                std::map<int, double> durationMap;
                std::map<int, double> ageMap;
                stateMap = myCellState.GetCellState(); 
                durationMap = myCellState.GetCellStateDuration();
                ageMap = myCellState.GetCellStateAge();
//                 for (std::map<int, std::string>::iterator mitr_cell=stateMap.begin();mitr_cell!=stateMap.end();mitr_cell++)//phase map contains the real dynamic cell index
//                 {  
//                     cout<<"Cell: "<<mitr_cell->first<<" State: "<<stateMap[mitr_cell->first]<<" Age: "<<ageMap[mitr_cell->first]<<" Duration: "<<durationMap[mitr_cell->first]<<endl;
//                     
//                 }
                std::map<int, double> cellPositionX;
                std::map<int, double> cellPositionY;
                std::map<int, double> cellPositionZ;
                cellPositionX = myCellState.GetCellPositionX();// get cell position
                cellPositionY = myCellState.GetCellPositionY();
                cellPositionZ = myCellState.GetCellPositionZ();
                for (std::map<int, std::string>::iterator mitr_cell=stateMap.begin();mitr_cell!=stateMap.end();mitr_cell++)//phase map contains the real dynamic cell index
                {  
                    double cX = cellPositionX.at(mitr_cell->first);
                    double cY = cellPositionY.at(mitr_cell->first);
                    double cZ = cellPositionZ.at(mitr_cell->first);
                    double massOfBystanderSignal = 1.66E-12;// unit in pico gram
                    double volume = d*d*d; 
                    double concentration = diffusionSolver.GetConcentration(cX, cY, cZ);
//                     cout<<"the concentration is "<<concentration<<endl;
                    double concentrationInMass = concentration*massOfBystanderSignal/volume;
//                     double concentrationInMass = 1E+4*massOfBystanderSignal/volume;
                    if (cycle_cellStateUPdate==1)// at initial time, give radiation exposure
                    {
                        myCellState.CellStateUpdate(mitr_cell->first,cellType,n*step_dose*40,0,deltaT_cellStateUpdate,1);// update cell state
                    }
                    else
                    {
                        myCellState.CellStateUpdate(mitr_cell->first,cellType,0,concentrationInMass,deltaT_cellStateUpdate,1);// update cell state
                    }
                    
                }

                diffusionSolver.CellStateUpdateInitialize();
                for (std::map<int, std::string>::iterator mitr_cell = stateMap.begin();mitr_cell!=stateMap.end();mitr_cell++)
                {     
                    if (stateMap[mitr_cell->first]=="S1")
                    {
                        double cX = cellPositionX.at(mitr_cell->first);
                        double cY = cellPositionY.at(mitr_cell->first);
                        double cZ = cellPositionZ.at(mitr_cell->first);
                        diffusionSolver.CellStateUpdate(mitr_cell->first,cX, cY, cZ,1,mu);// update bystander signalling state
                        
                    }
                    if (stateMap[mitr_cell->first]=="S21")
                    {
                        double cX = cellPositionX.at(mitr_cell->first);
                        double cY = cellPositionY.at(mitr_cell->first);
                        double cZ = cellPositionZ.at(mitr_cell->first);
                        diffusionSolver.CellStateUpdate(mitr_cell->first,cX, cY, cZ,2,mu);
                    }
                    if (stateMap[mitr_cell->first]=="S22")
                    {
                        double cX = cellPositionX.at(mitr_cell->first);
                        double cY = cellPositionY.at(mitr_cell->first);
                        double cZ = cellPositionZ.at(mitr_cell->first);
                        if (ageMap[mitr_cell->first]<=3)
                        {
                            diffusionSolver.CellStateUpdate(mitr_cell->first,cX, cY, cZ,2,mu);// maximum bystander signal secreting time
                        }
                        else
                        {
                            diffusionSolver.CellStateUpdate(mitr_cell->first,cX, cY, cZ,1,mu);
                        }
   
                    }
                    if (stateMap[mitr_cell->first]=="S3")
                    {
                        double cX = cellPositionX.at(mitr_cell->first);
                        double cY = cellPositionY.at(mitr_cell->first);
                        double cZ = cellPositionZ.at(mitr_cell->first);
                        if (ageMap[mitr_cell->first]<=3)
                        {
                            diffusionSolver.CellStateUpdate(mitr_cell->first,cX, cY, cZ,3,mu);// maximum bystander signal secreting time
                        }
                        else
                        {
                            diffusionSolver.CellStateUpdate(mitr_cell->first,cX, cY, cZ,1,mu);
                        }
                    }
                }
                    
                clock_cellStateUpdate=0;  
            }

        }
        
        std::map<int, std::string> stateMap;
        std::map<int, double> durationMap;
        std::map<int, double> ageMap;
        stateMap = myCellState.GetCellState();
        ageMap = myCellState.GetCellStateAge();
        durationMap = myCellState.GetCellStateDuration();
        std::map<int, int> cellAncestryIDMap;
        cellAncestryIDMap = myCellState.GetCellAncestryID();
        std::map<int, int > cellColonySizeMap;
        for (int i =0;i<cellNum;i++)
        {
            cellColonySizeMap[i] = 0; 
        }

        for (std::map<int, std::string>::iterator mitr_cell = stateMap.begin();mitr_cell!=stateMap.end();mitr_cell++)
        {
            cout<<"cell id: "<<mitr_cell->first <<" GetCellAncestryID is: "<<cellAncestryIDMap[mitr_cell->first]<<endl;
            cellColonySizeMap[cellAncestryIDMap[mitr_cell->first]] = cellColonySizeMap[cellAncestryIDMap[mitr_cell->first]]+1;
            
        }
        
        int colonySizeThreshold = 2;
        int colonyNumber = 0;
        
        for (std::map<int, int >::iterator mitr_cell = cellColonySizeMap.begin();mitr_cell!=cellColonySizeMap.end();mitr_cell++)
        {
            cout<<"The colonyNumber for cell "<<mitr_cell->first<< " is "<<mitr_cell->second<<endl;
            if (mitr_cell->second>=colonySizeThreshold)
            {
                colonyNumber = colonyNumber + 1;
            }
        }

        file1<<n*step_dose<<","<<stateMap.size() <<","<<(double)colonyNumber/cellNum<<endl;
                    
    }

    return 0;
}

