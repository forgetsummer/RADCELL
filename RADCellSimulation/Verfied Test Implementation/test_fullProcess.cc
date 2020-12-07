
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
    clock_t t1,t2;
    t1=clock();
    
    double xDim=1;//unit in mm
    double yDim=1;//unit in mm
    double zDim=1;//unit in mm
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

    ofstream file;
    ofstream file1;
    file.open("cellPhase.csv");
    file1.open("cellState.csv");
    int doseCaseNum = 50;
    int maxDose = 50;
    double step_dose = (double)maxDose/doseCaseNum;
    cout<<step_dose<<endl;
//     exit(0);
    
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(d,d,d);// cell home size is in unit of mm
//     layout.RectangularSlab(xDim,yDim,zDim,cellNum);// unit here is also mm, this line is for single cell simulation
    layout.Blob(xDim,cellNum);
    cout<<"The total cell number is "<<layout.GetCellNumber()<<endl;
    RADCellSimulation mySim; // define an object of Geant4 radiation transport module
//         mySim.SetCellWorld(0.1,0.1,0.1);//create cell world, dimension in unit of mm, this line is for testing single cell simulation
    double xWorld = 1.1*xDim;
    double yWorld = 1.1*yDim;
    double zWorld = 1.1*xDim;
    mySim.SetCellWorld(xWorld,yWorld,zWorld); // this line is for testing multiple cell simulation
    mySim.CreateCell("Epithelia","Cytoplasma Nucleus","Sphere","Blue Green");// create biological cells
    mySim.SetCellSimulationParameter(1,"Nucleus");
    
    for (int i=0; i<layout.GetCellNumber(); i++)
    {
        mySim.CellConstruction(i,"Epithelia",layout.GetCellPositionX(i),layout.GetCellPositionY(i),layout.GetCellPositionZ(i),20,5);
    }

    mySim.RADCellSimulationInitialize(argc, argv);// initialize the radiation transport module
    
    string doseFileName;
    string DNADamageFileName;
    std::string line;
    int particleNumber = 1000000;
    double particleEnergy = 1; // set up the particle energy

    ifstream ifile; // read the microdosimetry input file
    G4String testRadiationName="microdosimetry";
    G4String radiationInputFileName;
    radiationInputFileName=testRadiationName+".in";// original input file, could be defined for different sources
    ifile.open(radiationInputFileName.c_str());
            
    stringstream ss_beamOn;
    stringstream ss_energy;
    stringstream ss_out;
    stringstream ss_OutFile;
    string beamOnArgument;
    string energyArgument;
    string inputFileName;
    string outFileName;
            
    ss_energy<<"/gps/energy" <<" "<<particleEnergy<<" "<<"MeV"; // set up the radiation energy
    ss_beamOn<<"/run/beamOn" <<" "<<particleNumber; // set up beam on partilce number  
    ss_out<<testRadiationName<<"_"<<"Eng"<<"_"<< particleEnergy<<"PNum"<<"_"<<particleNumber<<".in";
    ss_OutFile<<testRadiationName<<"_"<<"Eng"<<"_"<< particleEnergy<<"PNum"<<"_"<<particleNumber;
    inputFileName = ss_out.str();
    outFileName = ss_OutFile.str();
    ofstream ofile(inputFileName.c_str());// write the changed microdosimetry input file
    doseFileName = outFileName+"_dose.csv";
    DNADamageFileName = outFileName+"_DNADamage.csv";
    beamOnArgument=ss_beamOn.str();
    energyArgument=ss_energy.str();

    if (ifile.is_open())
    {
        while (!ifile.eof())
        {
            std::getline(ifile,line);
            std::istringstream ls(line);
            std::string firstString;
            std::string secondString;
            ls>>firstString>>secondString;
//                 cout<<"The first string is "<<firstString <<" The second string is "<<secondString<<endl;
            if (firstString=="/run/beamOn")
            {
                line=beamOnArgument;
            }
            if (firstString=="/gps/ene/mono")
            {
                line = energyArgument;
            }
            ofile<<line<<endl;
        }
        
    }

    ifile.close();
    stringstream ss2;
    string simulationArgument;
    
    ss2<<"out"<<" "<<testRadiationName<<"_"<<"Eng"<<"_"<< particleEnergy<<"PNum"<<"_"<<particleNumber; // run argument in out mode
    simulationArgument=ss2.str();
            
//     mySim.EnergyDistributionCalculation("gui","microdosimetry.in",false); // run in gui mode, for checking the initial geometry
//     mySim.EnergyDistributionCalculation(simulationArgument.c_str(),inputFileName,false);// implementing radiation transport
// 
    
//*******************************************************************************************************************
//*******************************************************************************************************************
    
    for (int n=0;n<=doseCaseNum;n++)// simulation for different prescribed dose
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

        for (int i =0;i<layout.GetCellNumber();i++)//initialize the cell position 
        { 
            double cX = layout.GetCellPositionX(i)+xDim/2.0; 
            double cY = layout.GetCellPositionY(i)+yDim/2.0;
            double cZ = layout.GetCellPositionZ(i)+zDim/2.0;
            myCellState.CellPositionInitialization(i,cX,cY,cZ);
        }

        for (int i =0;i<layout.GetCellNumber();i++) // seed cellNum cells with random cell phase
        {
            myCellState.CellPhaseInitializationRandom(i,1);
//             myCellState.CellPhaseInitialization(i,1,"G1");
        }
                       
        for (int i =0;i<layout.GetCellNumber();i++)
        {
            myCellState.CellStateInitialization(i,"S1");
        }
        
        DiffusionReactionSolver diffusionSolver; // define an object for diffusion solver
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


    //*****************************************************************************************
    //****************Above is cell system CellPhaseInitialization
    //****************************************************************************************

        for(int i=0;i<totalTimeStepNum;i++)// Loop all the time steps, absolute time 
        {         
            cout<<"Get Time Step "<<i<<endl;
//             exit(0);
            clock_cellStateUpdate = clock_cellStateUpdate + T;
            clock_cellPhaseUpdate = clock_cellPhaseUpdate + T;
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
                    if (cycle_cellStateUPdate==1)
                    {
                        ReadRadiationTransportInfo radiationTransResulut;// define object to read MC simulation results
                        radiationTransResulut.ReadDoseTallyOutPut(doseFileName);
                        radiationTransResulut.ReadDNADamageTallyOutPut(DNADamageFileName);
                        radiationTransResulut.GetAbsoluteResultsByCOIMethod(n*step_dose,"indirect","mean");
                        double dose_cell = radiationTransResulut.GetAbsDoseOfCell(mitr_cell->first);
                        double DSB_cell = radiationTransResulut.GetAbsDSBOfCell(mitr_cell->first);
//                         cout<<"doseFileName "<<doseFileName<<endl;
//                         cout<<"prescribed Dose is "<<n*step_dose<<"Simulated Dose "<<dose_cell<<" The DSB number "<<DSB_cell<<endl;
                        
                        myCellState.CellStateUpdate(mitr_cell->first,cellType,DSB_cell,0,deltaT_cellStateUpdate,1);
                    }
                    else
                    {
                        myCellState.CellStateUpdate(mitr_cell->first,cellType,0,concentrationInMass,deltaT_cellStateUpdate,1);
                    }
                    
                }
                diffusionSolver.CellStateUpdateInitialize();// update the bystander signal signalling state of cells
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
        
        std::map<int, double> cellPositionX;
        std::map<int, double> cellPositionY;
        std::map<int, double> cellPositionZ;
        cellPositionX = myCellState.GetCellPositionX();// get cell position
        cellPositionY = myCellState.GetCellPositionY();
        cellPositionZ = myCellState.GetCellPositionZ();
        
        cout<<"Check the cell colony formation: "<<endl;
        mySim.UpdateGeometryInitialize();
        cout<<"after UpdateGeometryInitialize"<<endl;
        mySim.SetCellWorld(xWorld,yWorld,zWorld);//create cell world, dimension in unit of mm
        mySim.CreateCell("Epithelia","Cytoplasma Nucleus","Sphere","Blue Green");
        mySim.CreateCell("Epithelia_Daughter","Cytoplasma Nucleus","Sphere","Red Green");
        mySim.SetCellSimulationParameter(1," Nucleus");
        
        for (std::map<int, double>::iterator mitr_cell=cellPositionX.begin();mitr_cell!=cellPositionX.end();mitr_cell++)
        {
            if (mitr_cell->first>=0 && mitr_cell->first<cellNum)// if it is the mother celll
            {
                mySim.CellConstruction(mitr_cell->first,"Epithelia",cellPositionX.at(mitr_cell->first)-xDim/2,\
                 cellPositionY.at(mitr_cell->first)-yDim/2,cellPositionZ.at(mitr_cell->first)-zDim/2,20,5);
            }
            else
            {
                mySim.CellConstruction(mitr_cell->first,"Epithelia_Daughter",cellPositionX.at(mitr_cell->first)-xDim/2,\
                 cellPositionY.at(mitr_cell->first)-yDim/2,cellPositionZ.at(mitr_cell->first)-zDim/2,20,5);
            }

        }

        cout<<"after rebuild all cells"<<endl;
        mySim.UpdateGeometryFinalize();
//         mySim.EnergyDistributionCalculation( "out","microdosimetry.in",true);
//         mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
    //     mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
//     mySim.EnergyDistributionCalculation("gui","microdosimetry.in",true);
                    
    }
   
    
    t2=clock();
    float diff ((float)t2-(float)t1);
    cout<<"The processing time for the simulation is "<<diff / CLOCKS_PER_SEC<<"seconds"<<endl;

    return 0;
}

