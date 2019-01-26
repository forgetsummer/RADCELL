#include "RadiationSourceSettings.hh"
#include <sstream>
#include <fstream>
#include <stdlib.h>

void RadiationSourceSettings::RadiationPlaneSourceSetUp(string fileName, int particleNum, double particleEn,double cx, double cy, double cz,double xDim,double yDim)
{
    std::string line;
    particleNumber = particleNum;
    particleEnergy = particleEn; // set up the particle energy

    ifstream ifile; // read the microdosimetry input file
    G4String testRadiationName=fileName;
    G4String radiationInputFileName;
    radiationInputFileName=testRadiationName+".in";// original input file, could be defined for different sources
    ifile.open(radiationInputFileName.c_str());
            
    stringstream ss_beamOn;
    stringstream ss_energy;
    stringstream ss_halfX;
    stringstream ss_halfY;
    stringstream ss_center;
    stringstream ss_out;
    stringstream ss_OutFile;
    string beamOnArgument;
    string energyArgument;
    string planeXDimArgument;
    string planeYDimArgument;
    string planeCenterArgument;
    string inputFileName;
    string outFileName;
            
    ss_energy<<"/gps/energy" <<" "<<particleEnergy<<" "<<"MeV"; // set up the radiation energy
    ss_beamOn<<"/run/beamOn" <<" "<<particleNumber; // set up beam on partilce number 
    ss_halfX<<"/gps/pos/halfx"<<" "<<xDim/2<<" "<<"mm";
    ss_halfY<<"/gps/pos/halfy"<<" "<<yDim/2<<" "<<"mm";
    ss_center<<"/gps/pos/centre" <<" "<<cx<<" "<<cy<<" "<<cz<<" "<<"mm";
    ss_out<<testRadiationName<<"_"<<"Eng"<<"_"<< particleEnergy<<"PNum"<<"_"<<particleNumber<<".in";
    inputFileName = ss_out.str();
    ofstream ofile(inputFileName.c_str());// write the changed microdosimetry input file
    beamOnArgument=ss_beamOn.str();
    energyArgument=ss_energy.str();
    planeXDimArgument = ss_halfX.str();
    planeYDimArgument = ss_halfY.str();
    planeCenterArgument = ss_center.str();

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
            if (firstString=="/gps/ene/mono")// this is for mono energy source
            {
                line = energyArgument;
            }
            if(firstString=="/gps/energy")
            {
                line = energyArgument;
            }
            if (firstString == "/gps/pos/halfx")
            {
                line = planeXDimArgument;
            }
            if (firstString == "/gps/pos/halfy")
            {
                line = planeYDimArgument;
            }
            if (firstString == "/gps/pos/centre")
            {
                line = planeCenterArgument;
            }
            ofile<<line<<endl;
        }
        
    }

    ifile.close();
    stringstream ss2;
    
    ss2<<testRadiationName<<"_"<<"Eng"<<"_"<< particleEnergy<<"PNum"<<"_"<<particleNumber; // run argument in out mode
    simulationSourceID=ss2.str();

}

void RadiationSourceSettings::RadiationPointSourceSetUp(string fileName, int particleNum, double particleEn,double cx, double cy, double cz)
{
    std::string line;
    particleNumber = particleNum;
    particleEnergy = particleEn; // set up the particle energy
    

    ifstream ifile; // read the microdosimetry input file
    G4String testRadiationName=fileName;
    G4String radiationInputFileName;
    radiationInputFileName=testRadiationName+".in";// original input file, could be defined for different sources
    ifile.open(radiationInputFileName.c_str());
            
    stringstream ss_beamOn;
    stringstream ss_energy;
    stringstream ss_sourceCenter;
    stringstream ss_out;
    stringstream ss_OutFile;
    string beamOnArgument;
    string energyArgument;
    string sourceCenterArgument;
    string inputFileName;
    string outFileName;
            
    ss_energy<<"/gps/energy" <<" "<<particleEnergy<<" "<<"MeV"; // set up the radiation energy
    ss_beamOn<<"/run/beamOn" <<" "<<particleNumber; // set up beam on partilce number  
    ss_sourceCenter<<"/gps/pos/centre"<<" "<<cx<<" "<<cy<<" "<<cz<<" "<<"mm";
    ss_out<<testRadiationName<<"_"<<"Eng"<<"_"<< particleEnergy<<"PNum"<<"_"<<particleNumber<<".in";
    inputFileName = ss_out.str();
    ofstream ofile(inputFileName.c_str());// write the changed microdosimetry input file
    beamOnArgument=ss_beamOn.str();
    energyArgument=ss_energy.str();
    sourceCenterArgument = ss_sourceCenter.str();

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
            if (firstString=="/gps/ene/mono")// this is for mono energy source
            {
                line = energyArgument;
            }
            if(firstString=="/gps/energy")
            {
                line = energyArgument;
            }
            if(firstString == "/gps/pos/centre")
            {
                line = sourceCenterArgument;
            }
            ofile<<line<<endl;
        }
        
    }

    ifile.close();
    stringstream ss2;
    
    ss2<<testRadiationName<<"_"<<"Eng"<<"_"<< particleEnergy<<"PNum"<<"_"<<particleNumber; // run argument in out mode
    simulationSourceID=ss2.str();       

}

string RadiationSourceSettings::GetSimulationSourceID()
{
    return simulationSourceID;
}


double RadiationSourceSettings::GetSimulationParticleEnergy()
{
    return particleEnergy;
}
int RadiationSourceSettings::GetSimulationParticleNumber()
{
    return particleNumber;
}

