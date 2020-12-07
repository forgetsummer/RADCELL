#include "EventAction.hh"

#include "G4EventDoseAnalysis.hh"
#include "G4EventDNADamageAnalysis.hh"
#include "projectDataTypeGlobals.hh"
#include <fstream>
#include "G4RunDoseAnalysis.hh"
#include "G4RunDNADamageAnalysis.hh"
#include <Cell.hh>
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "ArgumentInterpreter.hh"


EventAction::EventAction()
{

}

EventAction::~EventAction()
{

}

void EventAction::StoreEdepEvent(EnergyDepositionEventInfo edepEventInfo)
{
    int eventID=edepEventInfo.eventID;// get cellID
    (edepEventMap[eventID]).push_back(edepEventInfo);

}


void EventAction::BeginOfEventAction(const G4Event* event)
{
G4UserEventAction::BeginOfEventAction(event);
}

void EventAction::EndOfEventAction(const G4Event* event)
{
    cout<<"Now the processing event ID is: "<<event->GetEventID()<<endl;
    
    if (event->GetEventID()==0)
    {
            G4String filename = "FromEventRunAction.csv";
    ofstream file1;
    file1.open(filename);
    for (EnergyDepositionEventMap::iterator mitr=edepEventMap.begin();mitr!=edepEventMap.end();mitr++ )
    {
    //cerr<<"The size of vector is "<<mitr->second.size()<<endl;
        for (edepEventVector::iterator vitr=mitr->second.begin();vitr!=mitr->second.end();vitr++)
        {
            file1<<(*vitr).eventID<<","<<(*vitr).cellID<<","<<(*vitr).pX<<","<<(*vitr).pY<<","<<(*vitr).pZ<<","<<(*vitr).edep<<","<<(*vitr).affectedCellOrganelle<<endl;
        }
    }
    file1.close();
    }

    G4String simulationMode = ArgumentInterpreter::GetSimulationMode();
    if (simulationMode=="out")
    {
        bool postProcessMarker=ArgumentInterpreter::GetPostProcessMarker();
        if (postProcessMarker)
        {
            G4String fileName = ArgumentInterpreter::GetOutPutFileName();
            G4String outputEdepFileName=fileName+"_edep"+".csv";// each worker will have each own writting file
            ofstream file;
            file.open(outputEdepFileName,std::ios_base::app);
            for (EnergyDepositionEventMap::iterator mitr=edepEventMap.begin();mitr!=edepEventMap.end();mitr++ )
            {
            //cerr<<"The size of vector is "<<mitr->second.size()<<endl;
                for (edepEventVector::iterator vitr=mitr->second.begin();vitr!=mitr->second.end();vitr++)
                {
                    file<<(*vitr).eventID<<","<<(*vitr).cellID<<","<<(*vitr).pX<<","<<(*vitr).pY<<","<<(*vitr).pZ<<","<<(*vitr).edep<<","<<(*vitr).affectedCellOrganelle<<endl;
                }
            }
            file.close();   
        }

    }

    CellDoseTally();
    CellDNADamageTally();
    edepEventMap.clear();// WARNING free memory after data processing

}

void EventAction::CellDoseTally()
{
    G4EventDoseAnalysis eventDoseAnalysis;

    EventDoseTallyMap currentEventDoseTallyMap;
    currentEventDoseTallyMap = eventDoseAnalysis.EventDoseTally(edepEventMap);
    G4RunDoseAnalysis::CollectEventDoseTallyMap(currentEventDoseTallyMap);


}
void EventAction::CellDNADamageTally()
{
    G4EventDNADamageAnalysis eventDNAAnalysis;
    EventDNADamagesTallyInfo currentEventDNADamagesTallyInfo;
    currentEventDNADamagesTallyInfo=eventDNAAnalysis.EventDNADamageTally(edepEventMap,0.0143,10.79,63.4215,3.2,2,2);
    G4RunDNADamageAnalysis::CollectEventDNADamageTallyMap(currentEventDNADamagesTallyInfo);

}

