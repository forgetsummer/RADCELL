#include "CellDNADamageAnalysis.hh"
#include "projectDataTypeGlobals.hh"
#include "G4EventDNADamageAnalysis.hh"
#include "CellSimulationStatisticalErrorAnalysis.hh"



#include <iostream>
#include <fstream>
#include <sstream>

void CellDNADamageAnalysis::ReadEnergyDepositionFile(string filename)
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
                    buffer>>totalParticleNumber; // get total particle number
                    
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
            }

            
            if (edepEventMap.find(eventID)==edepEventMap.end()) // which means we find a new eventID
            {
                // when we find a new eventID, then we start to process the edep points belonging to old eventID
                //EventDNADamageTally(EnergyDepositionEventMap edepEventMap,double SSBProb,double SSBEmin, double SSBEmax,double eps, int MinPts, double p)
                
                EventDNADamagesTallyInfo currentEventDNADamageTallyInfo; 
                G4EventDNADamageAnalysis eventDataAnalysis;
                currentEventDNADamageTallyInfo=eventDataAnalysis.EventDNADamageTally(edepEventMap,0.0143,10.79,63.4215,3.2,2,2);
                edepEventMap.clear();//WARNING clear the map after processing current event
                
                SSBTallyMap.insert(currentEventDNADamageTallyInfo.eventSSBMap.begin(),currentEventDNADamageTallyInfo.eventSSBMap.end());
                DSBTallyMap.insert(currentEventDNADamageTallyInfo.eventDSBMap.begin(),currentEventDNADamageTallyInfo.eventDSBMap.end());


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

        }
         line_num=line_num+1; 
    }
    
    EventDNADamagesTallyInfo currentEventDNADamageTallyInfo; 
    G4EventDNADamageAnalysis eventDataAnalysis;
    currentEventDNADamageTallyInfo=eventDataAnalysis.EventDNADamageTally(edepEventMap,0.0143,10.79,63.4215,3.2,2,2);
    edepEventMap.clear();//WARNING clear the map after processing current event
    
    SSBTallyMap.insert(currentEventDNADamageTallyInfo.eventSSBMap.begin(),currentEventDNADamageTallyInfo.eventSSBMap.end());
    DSBTallyMap.insert(currentEventDNADamageTallyInfo.eventDSBMap.begin(),currentEventDNADamageTallyInfo.eventDSBMap.end());
    N=currentEventDNADamageTallyInfo.N;

    file.close();


}

void CellDNADamageAnalysis::CellDNADamageTally()
{
    
    std::map<int, std::map<int, double> > SSB_event_cell_Map;// store the SSB number in this 2D map, first key is eventID, second key is cellID
    std::map<int, std::map<int, double> > DSB_event_cell_Map;
    
    for (EventDNASSBMap::iterator mitr_event=SSBTallyMap.begin();mitr_event!=SSBTallyMap.end();mitr_event++)
    {
        for (CellDNASSBMap::iterator mitr_cell=mitr_event->second.begin();mitr_cell!=mitr_event->second.end();mitr_cell++)
        {
            (SSB_event_cell_Map[mitr_event->first])[mitr_cell->first]=(SSBTallyMap[mitr_event->first])[mitr_cell->first].size()/double(N);
        }
    }
    
    for (EventDNADSBMap::iterator mitr_event=DSBTallyMap.begin();mitr_event!=DSBTallyMap.end();mitr_event++)
    {
        for (CellDNADSBMap::iterator mitr_cell=mitr_event->second.begin();mitr_cell!=mitr_event->second.end();mitr_cell++)
        {
            (DSB_event_cell_Map[mitr_event->first])[mitr_cell->first]=(DSBTallyMap[mitr_event->first])[mitr_cell->first].size()/double(N);
        }
    }
    
    
    
    cout<<"the size of DSB_event_cell_Map is "<<DSB_event_cell_Map.size()<<endl;
    cout<<"the size of SSB_event_cell_Map is "<<SSB_event_cell_Map.size()<<endl;

    
//     for ( std::map<int, std::map<int, double> >::iterator mitr_event=DSB_event_cell_Map.begin();mitr_event!=DSB_event_cell_Map.end();mitr_event++)
//     {
//         cout<<"the keys of DSB_event_cell_Map is "<<mitr_event->first<<endl;
//     }


    CellSimulationStatisticalErrorAnalysis DSBAnalysis;
    DSBAnalysis.SimulationStatisticalAnalysis<double>(DSB_event_cell_Map,totalParticleNumber);
    
    cellDSBMeanMap=DSBAnalysis.GetSimulationMean();
    cellDSBStdMap=DSBAnalysis.GetSimulationStd();
    
//     cout<<"The size of cellDSBMeanMap in CellDNADamageAnalysis is  "<<cellDSBMeanMap.size()<<endl;
    
//     for (std::map<int,double>::iterator mitr_cell=cellDSBMeanMap.begin();mitr_cell!=cellDSBMeanMap.end();mitr_cell++)
//     {
//         cout<<"The cellIDs in cellDSBMeanMap are "<<mitr_cell->first<<" The DSB number in this cell is "<<mitr_cell->second<<endl;
//     }
//     
    CellSimulationStatisticalErrorAnalysis SSBAnalysis;
    SSBAnalysis.SimulationStatisticalAnalysis<double>(SSB_event_cell_Map,totalParticleNumber);
    
    cellSSBMeanMap=SSBAnalysis.GetSimulationMean();
    cellSSBStdMap=SSBAnalysis.GetSimulationStd();
    
//     cout<<"the SSB mean of cell is "<<cellSSBMeanMap[0]<<endl;
//     cout<<"the SSB std of cell is "<<cellSSBStdMap[0]<<endl;
// 
//     cout<<"the DSB mean of cell is "<<cellDSBMeanMap[0]<<endl;
//     cout<<"the DSB std of cell is "<<cellDSBStdMap[0]<<endl;
                
    
//     cout<<"The size of cellSSBMeanMap in CellDNADamageAnalysis is "<<cellSSBMeanMap.size()<<endl;
//     
//     for (std::map<int,double>::iterator mitr_cell=cellSSBMeanMap.begin();mitr_cell!=cellSSBMeanMap.end();mitr_cell++)
//     {
//         cout<<"The cellIDs in cellSSBMeanMap are "<<mitr_cell->first<<" The SSB number in this cell is "<<mitr_cell->second<<endl;
//     }

}
map< int, double > CellDNADamageAnalysis::GetCellDSBNumberMeanMap()
{
    return cellDSBMeanMap;
}
map< int, double > CellDNADamageAnalysis::GetCellDSBNumberStdMap()
{
    return cellDSBStdMap;
}
map< int, double > CellDNADamageAnalysis::GetCellSSBNumberMeamMap()
{
    return cellSSBMeanMap;
}
map< int, double > CellDNADamageAnalysis::GetCellSSBNumberStdMap()
{
    return cellSSBStdMap;
}

void CellDNADamageAnalysis::WriteDNADamageTallyFile(string filename)
{
    ofstream file1;
    std::string  SSBFileName=filename+"_SSB.csv";
    std::string  DSBFileName=filename+"_DSB.csv";
    file1.open(SSBFileName.c_str());//open file for wrting SSB information of nucleus
    ofstream file2;
    file2.open(DSBFileName.c_str());// open file for writing DSB information of nuclues
    file1<<"Event ID"<<","<<"Cell ID"<<","<<"DNA Label"<<","<<"x"<<","<<"y"<<","<<"z"<<","<<"edep"<<endl;// write titles 
    for (EventDNASSBMap::iterator mitr_event=SSBTallyMap.begin();mitr_event!=SSBTallyMap.end();mitr_event++)
    {
            for (CellDNASSBMap::iterator mitr_cell=(mitr_event->second).begin();mitr_cell!=(mitr_event->second).end();mitr_cell++)
            {
                    for (int i=0; i<(mitr_cell->second).size();i++)
                    {
                            file1<<mitr_event->first<<","<<mitr_cell->first<<","<<((mitr_cell->second)[i]).DNALabel<<" , "<<((mitr_cell->second)[i]).edepPoint.x<<" , "
                            <<((mitr_cell->second)[i]).edepPoint.y<<" , "<<((mitr_cell->second)[i]).edepPoint.z<<" , "
                            <<((mitr_cell->second)[i]).edepPoint.edep<<endl;
                    }
            }
    }
    
    file2<<"Event ID"<<","<<"Cell ID"<<","<<"DSB Type"<<","<<"xC"<<","<<"yC"<<","<<"zC"<<","<<"DSBTotalEnergy"<<","<<"DSBDimension"<<endl;
    
    for (EventDNADSBMap::iterator mitr_event=DSBTallyMap.begin();mitr_event!=DSBTallyMap.end();mitr_event++)
    {
            for (CellDNADSBMap::iterator mitr_cell=(mitr_event->second).begin();mitr_cell!=(mitr_event->second).end();mitr_cell++)
            {
                    for (int i=0;i<(mitr_cell->second).size();i++)
                    {
                            file2<<mitr_event->first<<","<<mitr_cell->first<<","<<((mitr_cell->second)[i]).DSBType<<" , "<<((mitr_cell->second)[i]).xC<<" , "<<((mitr_cell->second)[i]).yC
                            <<" , "<<((mitr_cell->second)[i]).zC<<" , "<<((mitr_cell->second)[i]).DSBTotalEnergy
                            <<" , "<<((mitr_cell->second)[i]).DSBDimension<<endl;
                    }
            }
    }

}

map< int, vector< double > > CellDNADamageAnalysis::GetCellDSBDimensionMap()
{
    std::map<int, std::vector<double> > CellDSBDimensionMap;
    
    for (EventDNADSBMap::iterator mitr_event=DSBTallyMap.begin();mitr_event!=DSBTallyMap.end();mitr_event++)
    {
        for (CellDNADSBMap::iterator mitr_cell=mitr_event->second.begin();mitr_cell!=mitr_event->second.end();mitr_cell++)
        {
            for (int i=0;i<mitr_cell->second.size();i++)
            {
                CellDSBDimensionMap[mitr_cell->first].push_back(((DSBTallyMap[mitr_event->first])[mitr_cell->first])[i].DSBDimension);
            }
            
        }
    }

    
    return CellDSBDimensionMap;

}

void CellDNADamageAnalysis::ImportTotalParticleNumber(int num)
{
    totalParticleNumber = num;

}

void CellDNADamageAnalysis::ImportDNADamagesTallyMapForWholeRun(EventDNASSBMap ssbMap, EventDNADSBMap dsbMap, int num)
{
    SSBTallyMap = ssbMap;
    DSBTallyMap = dsbMap;
    N=num;

}


