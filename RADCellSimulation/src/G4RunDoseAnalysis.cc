#include "G4RunDoseAnalysis.hh"

#include "CellDoseAnalysis.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>

EventDoseTallyMap G4RunDoseAnalysis::doseTallyMap; // declare an instance of this static variable


void G4RunDoseAnalysis::CollectEventDoseTallyMap(EventDoseTallyMap doseTallyMapPerEvent)
{
    doseTallyMap.insert(doseTallyMapPerEvent.begin(),doseTallyMapPerEvent.end());

}

void G4RunDoseAnalysis::ImportTotalParticleNumber(int num)
{
    totalParticleNumber = num;

}

void G4RunDoseAnalysis::CellDoseTally(const DetectorConstruction* detector)
{
    CellDoseAnalysis doseAnalysis;
//     cout<<"The size of doseTallyMap in G4RunDataAnalysis is "<<doseTallyMap.size()<<endl;
    doseAnalysis.ImportDoseTallyMapForWholeRun(doseTallyMap);
    doseAnalysis.ImportTotalParticleNumber(totalParticleNumber);
    doseAnalysis.DoseTally();
    cellDoseMeanMap=doseAnalysis.GetTotalCellDoseMeanMap();
    cellDoseStdMap=doseAnalysis.GetTotalCellDoseStdMap();
    nucleusDoseMeanMap=doseAnalysis.GetTotalNucleusDoseMeanMap();
    nucleusDoseStdMap=doseAnalysis.GetTotalNucleusDoseStdMap();
    
    std::map<int, double> cellMassMap=detector->GetCellMassMap(); // get cell mass map, so as to calculate cell dose, mass unit is kilogram
    std::map<int, double> nucleusMassMap=detector->GetNucleusMassMap();
    
    for (std::map<int, double>::iterator mitr_cell=cellDoseMeanMap.begin();mitr_cell!=cellDoseMeanMap.end();mitr_cell++)
    {
       // Since the doseMap contains the total energy deposited inside cell, it is not real dose, so here we need to calcuate real dose 
       // by using equation dose=edep/mass 
       
        cellDoseMeanMap[mitr_cell->first]=cellDoseMeanMap[mitr_cell->first]*eV/joule; // change edep unit from eV to J
        cellDoseMeanMap[mitr_cell->first]=cellDoseMeanMap[mitr_cell->first]/(cellMassMap[mitr_cell->first]); // get dose, unit in Gy
    }
    
    for (std::map<int, double>::iterator mitr_cell=cellDoseStdMap.begin();mitr_cell!=cellDoseStdMap.end();mitr_cell++)
    {
        // get the std of dose (unit in Gy)
        cellDoseStdMap[mitr_cell->first] = cellDoseStdMap[mitr_cell->first]*eV/joule/(cellMassMap[mitr_cell->first]);
    }
       
    for (std::map<int, double>::iterator mitr_cell=nucleusDoseMeanMap.begin();mitr_cell!=nucleusDoseMeanMap.end();mitr_cell++)
    {
        // Since the doseMap contains the total energy deposited inside cell, it is not real dose, so here we need to calcuate real dose 
       // by using equation dose=edep/mass 
        nucleusDoseMeanMap[mitr_cell->first]=nucleusDoseMeanMap[mitr_cell->first]*eV/joule; // change edep unit from eV to J
        nucleusDoseMeanMap[mitr_cell->first]=nucleusDoseMeanMap[mitr_cell->first]/(nucleusMassMap[mitr_cell->first]);
        
    }
    
    for (std::map<int, double>::iterator mitr_cell=nucleusDoseStdMap.begin();mitr_cell!=nucleusDoseStdMap.end();mitr_cell++)
    {
        nucleusDoseStdMap[mitr_cell->first]=nucleusDoseStdMap[mitr_cell->first]*eV/joule/(nucleusMassMap[mitr_cell->first]);
    }
    

}

void G4RunDoseAnalysis::RadiationTransportDoseOutPut(string filename)
{
    G4String outputDoseFileName=filename+"_dose"+".csv";// write dose tally information out to a file
    ofstream file_Dose;
    file_Dose.open(outputDoseFileName);
    file_Dose<<"cellID"<<","<<"totalCellDose"<<","<<"statisticalError"<<","<<"nucleusDose"<<","<<"statisticalError"<<endl; // file title
        
    for (std::map<int, double>::iterator mitr_cell=cellDoseMeanMap.begin();mitr_cell!=cellDoseMeanMap.end();mitr_cell++)
    {
        file_Dose<<mitr_cell->first<<","<<cellDoseMeanMap[mitr_cell->first]<<","<<cellDoseStdMap[mitr_cell->first]<<","\
        <<nucleusDoseMeanMap[mitr_cell->first]<<","<<nucleusDoseStdMap[mitr_cell->first]<<endl;
        
    }
    file_Dose.close();

}

void G4RunDoseAnalysis::ClearRunData()
{
    doseTallyMap.clear();
    cellDoseMeanMap.clear();
    cellDoseStdMap.clear();
    nucleusDoseMeanMap.clear();
    nucleusDoseStdMap.clear();
}

