
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

int main(int argc,char** argv)
{
    int N=4;
    CellLayoutInitializer layout;
    layout.SetCellHomeParamter(0.05,0.05,0.05);// cell home size is in unit of mm
    layout.RectangularSlab(0.05, 0.05, 0.05,N);// unit here is also mm
    RADCellSimulation mySim;
    mySim.SetCellWorld(0.1,0.1,0.1);//create cell world, dimension in unit of mm
    mySim.CreateCell("Epithelia","Cytoplasma Nucleus","Sphere","Blue Green");
    mySim.SetCellSimulationParameter(1,"Nucleus ");
    for (int i=0; i<layout.GetCellNumber(); i++)
    {

        mySim.CellConstruction(i,"Epithelia",layout.GetCellPositionX(i),layout.GetCellPositionY(i),layout.GetCellPositionZ(i),20,5);
    }
  
    mySim.RADCellSimulationInitialize(argc, argv);



    clock_t t1,t2;
    t1=clock();
    std::string line;

    std::vector<int> particleNumberVec; //the vector storing the particle number for each run

    for (int i =1; i<=1;i++) // running five different experiment, set up the particle number in here
    {

        particleNumberVec.push_back(i*40); // initialize the vector 
    }

    for (int j=0;j<particleNumberVec.size();j++)
    {
        ifstream ifile; // read the microdosimetry input file
        ifile.open(("microdosimetry.in"));
        
        stringstream ss;
        stringstream ss1;
        string beamOnArgument;
        string inputFileName;

        ss<<"/run/beamOn" <<" "<<particleNumberVec[j];
        ss1<<"microdosimetry"<<"_"<<"particleNumber"<<"_"<<particleNumberVec[j]<<".in";
        inputFileName=ss1.str();
        ofstream ofile(inputFileName.c_str());// write the changed microdosimetry input file

        beamOnArgument=ss.str();

        
        if (ifile.is_open())
        {
            while (!ifile.eof())
            {
                std::getline(ifile,line);
                std::istringstream s(line);
                std::string firstString;
                std::string secondString;
                s>>firstString>>secondString;
//                 cout<<"The first string is "<<firstString <<" The second string is "<<secondString<<endl;
                if (firstString=="/run/beamOn")
                {
                    line=beamOnArgument;
                }
              
                ofile<<line<<endl;
            }
             
        }
       cout<<"test"<<j<<endl; // output test case number
       cout<<"the particle number is "<<particleNumberVec[j]<<endl;
       ifile.close();
       stringstream ss2;
       string simulationArgument;
       
       ss2<<"out"<<" "<<"OutPutFile"<<"_"<<particleNumberVec[j];
       simulationArgument=ss2.str();
       

       mySim.EnergyDistributionCalculation(simulationArgument.c_str(),inputFileName,true);
       
 

    }
     
    t2=clock();
    float diff ((float)t2-(float)t1);
    cout<<"The processing time for the simulation is "<<diff / CLOCKS_PER_SEC<<"seconds"<<endl;


    return 0;
}

