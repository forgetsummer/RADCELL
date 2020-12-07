
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
#include "ReactionDiffusionSimulation.hh"cout<<"get here"<<endl;
#include "DiffusionReactionSolver.hh"
#include "CellStateModel.hh"
#include "math.h"
#include <sstream>
#include <stdlib.h>


int main(int argc,char** argv)
{
   ReadRadiationTransportInfo testRead;
   testRead.ReadDoseTallyOutPut("microdosimetry_Eng_1PNum_1000000_dose.csv");
   testRead.ReadDNADamageTallyOutPut("microdosimetry_Eng_1PNum_1000000_DNADamage.csv");
   int cellID = 99;
   cout<<"The MC dose of the cell is "<<testRead.GetMCDoseOfCell(cellID)<<" and Std dose is "<<testRead.GetMCDoseStdOfCell(cellID)<<endl;
   cout<<"The MC DSB of the cell is "<<testRead.GetMCDSBOfCell(cellID)<<" and Std DSB is "<<testRead.GetMCDSBStdOfCell(cellID)<<endl;
   testRead.GetAbsoluteResultsByCOIMethod(1,"direct","mean");
   
   cout<<"The abs dose of the cell is "<<testRead.GetAbsDoseOfCell(cellID)<<" and Std dose is "<<testRead.GetAbsDoseStdOfCell(cellID)<<endl;
   cout<<"The abs DSB of the cell is "<<testRead.GetAbsDSBOfCell(cellID)<<" and Std DSB is "<<testRead.GetAbsDSBStdOfCell(cellID)<<endl;
   
   for (int i =0;i<99;i++)
   {
        cout<<"The abs dose of the cell "<<i<< " is "<<testRead.GetAbsDoseOfCell(i)<<" and Std dose is "<<testRead.GetAbsDoseStdOfCell(i)<<endl;
       
   }
   
   for (int i =0;i<99;i++)
   {
//        cout<<"The abs DSB of the cell "<<i<< " is "<<testRead.GetAbsDSBOfCell(i)<<" and Std DSB is "<<testRead.GetAbsDSBStdOfCell(i)<<endl;
       
   }
   testRead.GetAbsoluteResultsByFluenceMethod(1E+8,"direct");
    
   cout<<"The abs dose of the cell is "<<testRead.GetAbsDoseOfCell(cellID)<<" and Std dose is "<<testRead.GetAbsDoseStdOfCell(cellID)<<endl;
   cout<<"The abs DSB of the cell is "<<testRead.GetAbsDSBOfCell(cellID)<<" and Std DSB is "<<testRead.GetAbsDSBStdOfCell(cellID)<<endl;
   //*********************************************************************************************************
   //**********get net results using single cell results as reference
   //*********************************************************************************************************
   cout<<endl;
   cout<<"get net results using single cell results as reference"<<endl;
   testRead.ReadSingleCellDoseAsReference("Mono_Electron_Eng_0.2PNum_10000_dose.csv");
   testRead.ReadSingleCellDNADamageAsReference("Mono_Electron_Eng_0.2PNum_10000_DNADamage.csv");
   testRead.GetAbsoluteResultsByCOIMethod(1,"indirect","mean");
//    testRead.GetAbsoluteResultsByFluenceMethod(1E+8,"indirect");
    for (int i =0;i<99;i++)
   {
       cout<<"The abs dose of the cell "<<i<< " is "<<testRead.GetAbsDoseOfCell(i)<<" and Std dose is "<<testRead.GetAbsDoseStdOfCell(i)<<endl;
       
   }
   
   for (int i =0;i<99;i++)
   {
       cout<<"The abs DSB of the cell "<<i<< " is "<<testRead.GetAbsDSBOfCell(i)<<" and Std DSB is "<<testRead.GetAbsDSBStdOfCell(i)<<endl;  
   }
   
   
   return 0;
}

