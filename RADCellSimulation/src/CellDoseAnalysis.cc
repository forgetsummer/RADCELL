#include "CellDoseAnalysis.hh"
#include "CellSimulationStatisticalErrorAnalysis.hh"
#include "projectDataTypeGlobals.hh"
#include "EnergyDepositionEventInfo.hh"
#include "G4EventDoseAnalysis.hh"


#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

void CellDoseAnalysis::ReadEnergyDepositionFile(std::string filename) // this is the function for post anal
{
    ifstream file (filename.c_str());//open edep file
    std::string line;
    int line_num=0;
    
    
    EnergyDepositionEventMap edepEventMap;

    while (std::getline(file,line))
    {
        
        if (line_num==0) // read the first line
        {
            std::istringstream s(line);
            std::string field;
            int element_num=0;
            while(getline(s,field,','))
            {
                if (element_num==1)
                {
                    istringstream buffer(field);
                    buffer>>totalParticleNumber; // get total event number;
                    
                }
                element_num=element_num+1;
                
            }
            
        }

      
        if (line_num>1) // just starting read file from third line to end
        {
            std::istringstream s(line);
            std::string field;
            int element_num=0;
            
            int eventID=0;
            int cellID=0;
            double X=0;
            double Y=0;
            double Z=0;
            double edep=0;
            std::string organelle;
            

            while (getline(s,field,','))
            {
                if (element_num==0)
                {
                    istringstream buffer(field);
                    buffer>>eventID;
                }
                if (element_num==1)
                {
                    istringstream buffer(field);
                    buffer>>cellID;
                }
                if (element_num==2)
                {
                    istringstream buffer(field);
                    buffer>>X;
                }
                if (element_num==3)
                {
                    istringstream buffer(field);
                    buffer>>Y;
                }
                if(element_num==4)
                {
                    istringstream buffer(field);
                    buffer>>Z;
                }
                if (element_num==5)
                {
                    istringstream buffer(field);
                    buffer>>edep;
                }
                if (element_num==6)
                {
                    istringstream buffer(field);
                    buffer>>organelle;
                }
                
                element_num=element_num+1;
            }// end second while

            
            if (edepEventMap.find(eventID)==edepEventMap.end()) // which means we find a new eventID
            {
                // when we find a new eventID, then we start to process the edep points belonging to old eventID
                G4EventDoseAnalysis eventDataAnalysis;
                EventDoseTallyMap doseTallyMapForCurrentEvent=eventDataAnalysis.EventDoseTally(edepEventMap); // get doseTallyMap for current event
                doseTallyMapForWholeRun.insert(doseTallyMapForCurrentEvent.begin(),doseTallyMapForCurrentEvent.end()); // insert current event tally information to whole run tally information
                edepEventMap.clear(); // WARNING clear the map after processing current event
 
            }

            EnergyDepositionEventInfo edepPoint;
            edepPoint.eventID=eventID;
            edepPoint.cellID=cellID;
            edepPoint.pX=X;
            edepPoint.pY=Y;
            edepPoint.pZ=Z;
            edepPoint.edep=edep;
            edepPoint.affectedCellOrganelle=organelle;
            edepEventMap[eventID].push_back(edepPoint);

        } // end if 
         line_num=line_num+1; 
    }// end first while
    
    G4EventDoseAnalysis eventDataAnalysis;
    EventDoseTallyMap doseTallyMapForCurrentEvent=eventDataAnalysis.EventDoseTally(edepEventMap); // get doseTallyMap for current event
    doseTallyMapForWholeRun.insert(doseTallyMapForCurrentEvent.begin(),doseTallyMapForCurrentEvent.end());
    edepEventMap.clear();// WARNING clear the map after processing current event
    
    file.close();


}



void CellDoseAnalysis::DoseTally()
{
    std::map<int, std::map<int, double> >event_cell_TotalCellDose;
    std::map<int, std::map<int, double> >event_cell_NucleusDose;
    std::map<int, std::map<int, double> >event_cell_CytoplasmDose;
   
    for (EventDoseTallyMap::iterator mitr_event=doseTallyMapForWholeRun.begin();mitr_event!=doseTallyMapForWholeRun.end();mitr_event++)
    {
        for (CellDoseTallyMap::iterator mitr_cell=mitr_event->second.begin();mitr_cell!=mitr_event->second.end();mitr_cell++)
        {
            double sum_dose=0;
            for (OrganelleDoseMap::iterator mitr_organelle=mitr_cell->second.begin();mitr_organelle!=mitr_cell->second.end();mitr_organelle++)
            {
                sum_dose = sum_dose +mitr_organelle->second;
                if (mitr_organelle->first=="Nucleus")
                {
                    event_cell_NucleusDose[mitr_event->first][mitr_cell->first]=event_cell_NucleusDose[mitr_event->first][mitr_cell->first]+mitr_organelle->second;
                }
                if (mitr_organelle->first == "Cytoplasma")
                {
                    event_cell_CytoplasmDose[mitr_event->first][mitr_cell->first]=event_cell_CytoplasmDose[mitr_event->first][mitr_cell->first]+mitr_organelle->second;
                }
            }
            
            event_cell_TotalCellDose[mitr_event->first][mitr_cell->first] =  event_cell_TotalCellDose[mitr_event->first][mitr_cell->first]+sum_dose;
 
        }
    }
    
    cout<<"The size of event_cell_TotalCellDose is "<<event_cell_TotalCellDose.size()<<endl;
    
    CellSimulationStatisticalErrorAnalysis myAnalysis;
    myAnalysis.SimulationStatisticalAnalysis<double>(event_cell_TotalCellDose,totalParticleNumber);
    cellDoseMeanMap=myAnalysis.GetSimulationMean();
    cellDoseStdMap=myAnalysis.GetSimulationStd();
    CellSimulationStatisticalErrorAnalysis myAnalysis1;
    cout<<"the size of event_cell_NucleusDose is "<<event_cell_NucleusDose.size()<<endl;
    myAnalysis1.SimulationStatisticalAnalysis<double>(event_cell_NucleusDose,totalParticleNumber);
    nucleusDoseMeanMap=myAnalysis1.GetSimulationMean();
    nucleusDoseStdMap=myAnalysis1.GetSimulationStd();
    
    CellSimulationStatisticalErrorAnalysis myAnalysis2;
    myAnalysis2.SimulationStatisticalAnalysis<double>(event_cell_CytoplasmDose,totalParticleNumber);
    
//     cout<<"The cellDoseMeanMap is "<<cellDoseMeanMap[0]<<endl;
//     cout<<"The cellDoseStdMap is "<<cellDoseStdMap[0]<<endl;
//     cout<<"The nucleusDoseMeanMap is "<<nucleusDoseMeanMap[0]<<endl;
//     cout<<"The nucleusDoseStdMap is "<<nucleusDoseStdMap[0]<<endl;
//     

}

void CellDoseAnalysis::WriteCellDoseAnalysisFile(string doseTallyFileName)
{
    // this function is for writing the dose tally information to harddrive, so as to be used later for cell state analysis.
    // the writing format will be as :
    // cellID, totalCellDose,+-error, nucleusDose,+-errror;
    // it will be output as csv file
    
    ofstream file;
    file.open(doseTallyFileName.c_str());
    file<<"cellID"<<","<<"totalCellDose"<<","<<"statisticalError"<<","<<"nucleusDose"<<","<<"statisticalError"<<endl; // file title
    
    for (std::map<int, double>::iterator mitr_cell=cellDoseMeanMap.begin();mitr_cell!=cellDoseMeanMap.end();mitr_cell++)
    {
        file<<mitr_cell->first<<","<<cellDoseMeanMap[mitr_cell->first]<<","<<cellDoseStdMap[mitr_cell->first]<<","<<nucleusDoseMeanMap[mitr_cell->first]<<","<<nucleusDoseStdMap[mitr_cell->first]<<endl;
        
    }
    

}

void CellDoseAnalysis::ImportDoseTallyMapForWholeRun(EventDoseTallyMap doseTallyMap)
{
    doseTallyMapForWholeRun=doseTallyMap;

}

void CellDoseAnalysis::ImportTotalParticleNumber(int num)
{
    totalParticleNumber = num;

}


map< int, double > CellDoseAnalysis::GetTotalCellDoseMeanMap()
{
    return cellDoseMeanMap;
}
map< int, double > CellDoseAnalysis::GetTotalCellDoseStdMap()
{
    return cellDoseStdMap;
}
map< int, double > CellDoseAnalysis::GetTotalNucleusDoseMeanMap()
{
    return nucleusDoseMeanMap;
}
map< int, double > CellDoseAnalysis::GetTotalNucleusDoseStdMap()
{
    return nucleusDoseStdMap;
}

