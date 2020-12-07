
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
        clock_t t1, t2;
        t1=clock();
        cout<<"Begin Analysis:"<<endl;
        
        std::vector<int> particleNumberVec;
        
        
        for (int i =1; i<=5;i++) // running five different experiment, set up the particle number in here
        {

            particleNumberVec.push_back(i*40); // initialize the vector 
        }
        
        for (int i=0;i<particleNumberVec.size();i++)
        {
            for (int j=0;j<=1;j++)
            {
                cout<<"test"<<i<<" and Thread "<<j<<endl;
                stringstream ss;
                ss<<"OutPutFile_"<<particleNumberVec[i]<<"_"<<"edep"<<j<<".csv";
                string edepForDNA=ss.str();
            
            

//         string edepForDNA="OutPutFile_4000_edep-2.csv";
        
                CellDNADamageAnalysis myDNAAnalysis;
                CellDoseAnalysis doseAnalysis;
                
                myDNAAnalysis.ReadEnergyDepostionFile(edepForDNA);
                doseAnalysis.ReadEnergyDepostionFile(edepForDNA);
                
                std::map<int, double> cellDoseMeanMap;// map storing the mean cell total dose, the key is cellID, value is total dose of cell
                std::map<int, double> cellDoseStdMap; // map storing the statistical error of mean cell total dose, the key is cellID
                doseAnalysis.DoseTally();
                cellDoseMeanMap=doseAnalysis.GetTotalCellDoseMeanMap();
                cellDoseStdMap=doseAnalysis.GetTotalCellDoseStdMap();
                
            
        //         myDNAAnalysis.testfunction();
                
                std::map<int, double> cellDSBMeanMap;
                std::map<int, double> cellDSBStdMap;
                
                std::map<int, double> cellSSBMeanMap;
                std::map<int, double> cellSSBStdMap;
        //         myDNAAnalysis.DNADamageTally(0.0143,5,57.63,3.2,2,2);// 57.63 is from old paper of firelind, 37.5 is new 
        //         myDNAAnalysis.DNADamageTally(0.0143,5,37.5,3.2,2,2);// 57.63 is from old paper of firelind, 37.5 is new 
        //         myDNAAnalysis.DNADamageTally(0.0143,5,57.63,2000,2,2);// 57.63 is from old paper of firelind, 37.5 is new 
                myDNAAnalysis.DNADamageTally(0.0143,10.79,57.63,3.2,2,2);// 57.63 is from old paper of firelind, 37.5 is new 

                cellDSBMeanMap=myDNAAnalysis.GetCellDSBNumberMeanMap();
                cellDSBStdMap=myDNAAnalysis.GetCellDSBNumberStdMap();
            
                
                cellSSBMeanMap=myDNAAnalysis.GetCellSSBNumberMeamMap();
                cellSSBStdMap=myDNAAnalysis.GetCellSSBNumberStdMap();
                
                cout<<"the cell dose is "<<cellDoseMeanMap[0]<<endl;
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
                
            }
           
        
        }

    
    t2=clock();
    float diff ((float)t2-(float)t1);
    cout<<"The processing time for the simulation is "<<diff / CLOCKS_PER_SEC<<"seconds"<<endl;

        
   
    return 0;
}

