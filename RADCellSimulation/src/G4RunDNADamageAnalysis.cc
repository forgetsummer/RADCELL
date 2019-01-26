#include "G4RunDNADamageAnalysis.hh"

#include "CellDoseAnalysis.hh"
#include "CellDNADamageAnalysis.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>
#include <stdlib.h>
#include <cmath>

EventDNASSBMap G4RunDNADamageAnalysis::SSBTallyMap;
EventDNADSBMap G4RunDNADamageAnalysis::DSBTallyMap;
int G4RunDNADamageAnalysis::N;



void G4RunDNADamageAnalysis::CollectEventDNADamageTallyMap(EventDNADamagesTallyInfo dnaDNADamagesTallyInfoPerEvent)
{
    SSBTallyMap.insert(dnaDNADamagesTallyInfoPerEvent.eventSSBMap.begin(),dnaDNADamagesTallyInfoPerEvent.eventSSBMap.end());
    DSBTallyMap.insert(dnaDNADamagesTallyInfoPerEvent.eventDSBMap.begin(),dnaDNADamagesTallyInfoPerEvent.eventDSBMap.end());
    N=dnaDNADamagesTallyInfoPerEvent.N;


}


void G4RunDNADamageAnalysis::ImportTotalParticleNumber(int num)
{
    totalParticleNumber = num;

}



void G4RunDNADamageAnalysis::CellDNADamageTally(map< int, double > cellDoseMeanMap, map< int, double > cellDoseStdMap)
{
    CellDNADamageAnalysis dnaAnalysis;
    dnaAnalysis.ImportDNADamagesTallyMapForWholeRun(SSBTallyMap,DSBTallyMap,N);
    dnaAnalysis.ImportTotalParticleNumber(totalParticleNumber);
    dnaAnalysis.CellDNADamageTally();
    cellSSBMeanMap=dnaAnalysis.GetCellSSBNumberMeamMap();
    cellSSBStdMap=dnaAnalysis.GetCellSSBNumberStdMap();
    cellDSBMeanMap=dnaAnalysis.GetCellDSBNumberMeanMap();
    cellDSBStdMap=dnaAnalysis.GetCellDSBNumberStdMap();
    cellDSBDimensionMap=dnaAnalysis.GetCellDSBDimensionMap();
    
    for (std::map<int,double>::iterator mitr_cell=cellDSBMeanMap.begin();mitr_cell!=cellDSBMeanMap.end();mitr_cell++)
    {
        ratio_DSB_Gy_Gbp_MeanMap[mitr_cell->first]= cellDSBMeanMap[mitr_cell->first]/cellDoseMeanMap[mitr_cell->first]/double(6);
        ratio_DSB_Gy_Gbp_StdMap[mitr_cell->first] = 1/pow(cellDoseMeanMap[mitr_cell->first],2)*sqrt(pow(cellDoseMeanMap[mitr_cell->first],2)\
            *pow(cellDSBStdMap[mitr_cell->first],2)+pow(cellDSBMeanMap[mitr_cell->first],2)*pow(cellDoseStdMap[mitr_cell->first],2)
        )/double(6); // this calculating the error of DSB according the error propagation rule
    }
    
    for (std::map<int, double>::iterator mitr_cell=cellSSBMeanMap.begin();mitr_cell!=cellSSBMeanMap.end();mitr_cell++)
    {
        ratio_SSB_Gy_Gbp_MeanMap[mitr_cell->first] = cellSSBMeanMap[mitr_cell->first]/cellDoseMeanMap[mitr_cell->first]/double(6);
        ratio_SSB_Gy_Gbp_StdMap[mitr_cell->first] = 1/pow(cellDoseMeanMap[mitr_cell->first],2)*sqrt(pow(cellDoseMeanMap[mitr_cell->first],2)\
            *pow(cellSSBStdMap[mitr_cell->first],2)+pow(cellSSBMeanMap[mitr_cell->first],2)*pow(cellDoseStdMap[mitr_cell->first],2)
        )/double(6); // calculating the error of SSB
        
    }

}




void G4RunDNADamageAnalysis::RadiationTransportDNADamageOutPut(string filename)
{
    G4String outputDNADamageFileName=filename+"_DNADamage"+".csv"; // write DSB tally information out to a file
    ofstream file_DNA;
    file_DNA.open(outputDNADamageFileName);
    file_DNA<<"cellID"<<","<<"SSB"<<","<<"statisticalError"<<","<<"DSB"<<","<<"statisticalError"<<endl;
    
    for (std::map<int, double>::iterator mitr_cell=cellSSBMeanMap.begin();mitr_cell!=cellSSBMeanMap.end();mitr_cell++)
    {
        file_DNA<<mitr_cell->first<<","<<cellSSBMeanMap[mitr_cell->first]<<","<<cellSSBStdMap[mitr_cell->first]<<","\
        <<cellDSBMeanMap[mitr_cell->first]<<","<<cellDSBStdMap[mitr_cell->first]<<endl;
    }
    file_DNA.close();
    
    G4String outputDSBDimensionFileName= filename+"_DSBDimension"+".csv";
   
    ofstream file_DSBDimension;
    file_DSBDimension.open(outputDSBDimensionFileName);
   
    file_DSBDimension<<"cellID"<<","<<"DSBDimension"<<endl;
   
    for (std::map<int, std::vector<double> >::iterator mitr_cell=cellDSBDimensionMap.begin();mitr_cell!=cellDSBDimensionMap.end();mitr_cell++)
    {
       for (int i=0;i<mitr_cell->second.size();i++)
       {
           file_DSBDimension<<mitr_cell->first<<","<<(mitr_cell->second)[i]<<endl;
           
       }
       
    }
   
    file_DSBDimension.close();
    
    G4String outputDSB_Gy_GbpFileName=filename+"_DNADamage_Gy_Gbp"+".csv"; // write ratio_DSB_Gy_Gbp and ratio_SSB_Gy_Gbp out to file
   
    ofstream file_DNADamage;
    file_DNADamage.open(outputDSB_Gy_GbpFileName);
    file_DNADamage<<"cellID"<<","<<"ratio_DSB_Gy_Gbp"<<","<<"ratio_DSB_Gy_Gbp_Std"<<","<<"ratio_SSB_Gy_Gbp"<<","<<"ratio_SSB_Gy_Gbp_Std"<<endl;
   
    for (std::map<int, double>::iterator mitr_cell=ratio_SSB_Gy_Gbp_MeanMap.begin();mitr_cell!=ratio_SSB_Gy_Gbp_MeanMap.end();mitr_cell++)
    {
       file_DNADamage<<mitr_cell->first<<","<< ratio_DSB_Gy_Gbp_MeanMap[mitr_cell->first]<<","<<ratio_DSB_Gy_Gbp_StdMap[mitr_cell->first]\
       <<","<<ratio_SSB_Gy_Gbp_MeanMap[mitr_cell->first]<<","<<ratio_SSB_Gy_Gbp_StdMap[mitr_cell->first]<<endl;
       
    }
    file_DNADamage.close();

}

void G4RunDNADamageAnalysis::ClearRunData()
{
    SSBTallyMap.clear();
    DSBTallyMap.clear();
    cellDSBMeanMap.clear();
    cellDSBStdMap.clear();
    cellSSBMeanMap.clear();
    cellSSBStdMap.clear();
    cellDSBDimensionMap.clear();
    ratio_DSB_Gy_Gbp_MeanMap.clear();
    ratio_DSB_Gy_Gbp_StdMap.clear();
    ratio_SSB_Gy_Gbp_MeanMap.clear();
    ratio_SSB_Gy_Gbp_StdMap.clear();
    cellDoseMeanMap.clear();
}


