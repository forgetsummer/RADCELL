
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

    int N=1; // this line for seeding 300 cell
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(0.05,0.05,0.05);// cell home size is in unit of mm
    layout.RectangularSlab(0.05, 0.05, 0.05,N);// unit here is also mm, this line is for single cell simulation
    //         layout.RectangularSlab(1,1,1,N); // this line is for multiple cell simulation
    RADCellSimulation mySim;
    mySim.SetCellWorld(0.1,0.1,0.1);//create cell world, dimension in unit of mm, this line is for testing single cell simulation

    //         mySim.SetCellWorld(1.1,1.1,1.1); // this line is for testing multiple cell simulation
    mySim.CreateCell("Epithelia","Cytoplasma Nucleus","Sphere","Blue Green");
    //         mySim.SetCellSimulationParameter(1,"Nucleus");
    mySim.SetCellSimulationParameter(1,"Nucleus Cytoplasma");
    for (int i=0; i<layout.GetCellNumber(); i++)
    {

        mySim.CellConstruction(i,"Epithelia",layout.GetCellPositionX(i),layout.GetCellPositionY(i),layout.GetCellPositionZ(i),20,5);
    }

    mySim.RADCellSimulationInitialize(argc, argv);


    clock_t t1,t2;
    t1=clock();
    std::string line;
    
    int particleNumber = 100000;
    
    string particleType;
    particleType = "e";
    std::vector<double> particleEnergyVec; // the vector storing the particle energy for each run
    
    if (particleType=="p") // if the particle is proton
    {
        for (int i=2;i<=50;i++)
        {
            particleEnergyVec.push_back(i*2); // the second energy groups, delt_E=2MeV
        } 
    }
    if (particleType == "e") // if the particle is electron
    {
        for (int i=1;i<=4;i++)
        {
            particleEnergyVec.push_back(i*0.2); // delt_E=0.2MeV
        }
        
//         for (int i=6;i<=10;i++)
//         {
//             particleEnergyVec.push_back(i*1); //  delt_E=1MeV
//         }
    }


    for (int i=0;i<particleEnergyVec.size();i++) // loop energy groups
    {
        ifstream ifile; // read the microdosimetry input file
        G4String testRadiationName;
        if (particleType=="e")
        {
            testRadiationName="Mono_Electron";
        }
        if (particleType=="p")
        {
            testRadiationName="Mono_Proton";
        }
        G4String radiationInputFileName;
        radiationInputFileName=testRadiationName+".in";
        ifile.open(radiationInputFileName.c_str());
        
        stringstream ss_beamOn;
        stringstream ss_energy;
        stringstream ss1;
        string beamOnArgument;
        string energyArgument;
        string inputFileName;
        
        ss_energy<<"/gps/energy" <<" "<<particleEnergyVec[i]<<" "<<"MeV"; // set up the radiation energy
        ss_beamOn<<"/run/beamOn" <<" "<<particleNumber; // set up beam on partilce number
        
        ss1<<testRadiationName<<"_"<<"Eng"<<"_"<< particleEnergyVec[i]<<"PNum"<<"_"<<particleNumber<<".in";
        inputFileName=ss1.str();
        ofstream ofile(inputFileName.c_str());// write the changed microdosimetry input file

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
                if (firstString=="/gps/energy")
                {
                    line = energyArgument;
                }
            
                ofile<<line<<endl;
            }
            
        }
        cout<<"the particle number is "<<particleNumber<<endl;
        ifile.close();
        stringstream ss2;
        string simulationArgument;
        ss2<<"out"<<" "<<testRadiationName<<"_"<<"Eng"<<"_"<< particleEnergyVec[i]<<"PNum"<<"_"<<particleNumber; // run argument in out mode
        simulationArgument=ss2.str();
        mySim.EnergyDistributionCalculation(simulationArgument.c_str(),inputFileName,false);
//                 mySim.EnergyDistributionCalculation("gui",inputFileName,true); // run in gui mode
        
    }

    t2=clock();
    float diff ((float)t2-(float)t1);
    cout<<"The processing time for the simulation is "<<diff / CLOCKS_PER_SEC<<"seconds"<<endl;

    return 0;
}

