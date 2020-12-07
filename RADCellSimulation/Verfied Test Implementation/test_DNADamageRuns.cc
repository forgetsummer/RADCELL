
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
    int runmode=0;
    if (runmode==0)
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

        std::vector<int> particleNumberVec; //the vector storing the particle number for each run

        for (int i =1; i<2;i++) // running five different experiment, set up the particle number in here
        {
            particleNumberVec.push_back(i*10); // initialize the vector 
        }
        //////////////////////////////////////////////////////////////////////////////////
        string particleType = "e";
        std::vector<double> particleEnergyVec; // the vector storing the particle energy for each run
        if (particleType=="e")// this is for electron simulation
        {
            for (int i=1;i<=25;i++)
            {
                particleEnergyVec.push_back(i*0.2); // the first 25 energy groups, from 0.2MeV to 5MeV, delta_E=0.2MeV
            }
            
            for (int i=6;i<=100;i++)
            {
                particleEnergyVec.push_back(i); // from 6MeV to 100MeV
            }
        }
        if (particleType == "p") // this is for proton simulation
        {      

//             for (int i=1;i<=10;i++)
//             {
//                 particleEnergyVec.push_back(i*1); // the first 10 energy groups, from 0.1MeV to 1MeV, delta_E=0.1MeV
//             }
//             
            for (int i=2;i<=50;i++)
            {
                particleEnergyVec.push_back(i*2); // the second energy groups, from 1MeV to 20MeV, delta_E=1MeV
            } 
        }
        
        G4String testRadiationName;
        if (particleType=="e")
        {
            testRadiationName ="Mono_Electron";// choose the input file name
        }
        if (particleType == "p")
        {
            testRadiationName = "Mono_Proton";
        }
        string fileName;
        fileName = testRadiationName+".csv";
        ofstream ofileID(fileName.c_str());// output the input file names to a file, for later post simulation processing
        for (int i =0;i<particleEnergyVec.size();i++) // loop energy groups
        {
            for (int j=0;j<1;j++) // loop particle number tests
            {
                ifstream ifile; // read the microdosimetry input file
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
                ss_beamOn<<"/run/beamOn" <<" "<<particleNumberVec[j]; // set up beam on partilce number
              
                ss1<<testRadiationName<<"_"<<"Eng"<<"_"<< particleEnergyVec[i]<<"PNum"<<"_"<<particleNumberVec[j]<<".in";
                inputFileName=ss1.str();
                ofstream ofile(inputFileName.c_str());// write the changed microdosimetry input file
                
                ofileID<<inputFileName<<endl;
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
                cout<<"test"<<j<<endl; // output test case number
                cout<<"the particle number is "<<particleNumberVec[j]<<endl;
                ifile.close();
                stringstream ss2;
                string simulationArgument;
                
                ss2<<"out"<<" "<<testRadiationName<<"_"<<"Eng"<<"_"<< particleEnergyVec[i]<<"PNum"<<"_"<<particleNumberVec[j]; // run argument in out mode
                simulationArgument=ss2.str();
//                 mySim.EnergyDistributionCalculation(simulationArgument.c_str(),inputFileName,false);
                
                mySim.EnergyDistributionCalculation("gui",inputFileName,true); // run in gui mode
            }
            
        }

        t2=clock();
        float diff ((float)t2-(float)t1);
        cout<<"The processing time for the simulation is "<<diff / CLOCKS_PER_SEC<<"seconds"<<endl;
    
    
    }
    
    if (runmode==1)
    {
        clock_t t1, t2;
        t1=clock();
        cout<<"Begin Analysis:"<<endl;
        
        CellDoseAnalysis myAnalysis;
        G4String edepForDNA="microdosimetry_100_edep.csv";
        myAnalysis.ReadEnergyDepositionFile(edepForDNA);
        myAnalysis.DoseTally();
        
        std::map<int, double> cellDoseMeanMap;// map storing the mean cell total dose, the key is cellID, value is total dose of cell
        std::map<int, double> cellDoseStdMap; // map storing the statistical error of mean cell total dose, the key is cellID
        cellDoseMeanMap=myAnalysis.GetTotalCellDoseMeanMap();
        cellDoseStdMap=myAnalysis.GetTotalCellDoseStdMap();
        
        CellDNADamageAnalysis dnaAnalysis;
        dnaAnalysis.ReadEnergyDepositionFile(edepForDNA);
        dnaAnalysis.CellDNADamageTally();
        
        std::map<int, double> cellDSBMeanMap;
        std::map<int, double> cellDSBStdMap;
        
        std::map<int, double> cellSSBMeanMap;
        std::map<int, double> cellSSBStdMap;
        
        cellDSBMeanMap=dnaAnalysis.GetCellDSBNumberMeanMap();
        cellDSBStdMap=dnaAnalysis.GetCellDSBNumberStdMap();
    
        
        cellSSBMeanMap=dnaAnalysis.GetCellSSBNumberMeamMap();
        cellSSBStdMap=dnaAnalysis.GetCellSSBNumberStdMap();
        
        cout<<"the cell dose is "<<cellDoseMeanMap[268]<<endl;
        cout<<"the cell dose std is "<<cellDoseStdMap[0]<<endl;
        
        cout<<"the SSB mean of cell is "<<cellSSBMeanMap[0]<<endl;
        cout<<"the SSB std of cell is "<<cellSSBStdMap[0]<<endl;
        


        cout<<"the DSB mean of cell is "<<cellDSBMeanMap[0]<<endl;
        cout<<"the DSB std of cell is "<<cellDSBStdMap[0]<<endl;
        
        double ratio_DSB_Gy_Gbp;
        ratio_DSB_Gy_Gbp=cellDSBMeanMap[0]/(cellDoseMeanMap[0]*2.8656)*pow10(8);
        cout<<"the ratio_DSB_Gy_Gbp is "<<ratio_DSB_Gy_Gbp<<endl;
        
        double ratio_SSB_Gy_Gbp;
        ratio_SSB_Gy_Gbp=cellSSBMeanMap[0]/(cellDoseMeanMap[0]*2.8656)*pow10(8);
        cout<<"the ratio_SSB_Gy_Gbp is "<<ratio_SSB_Gy_Gbp<<endl;
        
        t2=clock();
        float diff ((float)t2-(float)t1);
        cout<<"The processing time for the simulation is "<<diff / CLOCKS_PER_SEC<<"seconds"<<endl;
        

        
    }
   
    
    
    
    
    ////////////////////////////////////////////
    

    
    

   
    return 0;
}

