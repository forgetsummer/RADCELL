#include "G4EventDoseAnalysis.hh"

EventDoseTallyMap G4EventDoseAnalysis::EventDoseTally(EnergyDepositionEventMap edepEventMap)
{
    EventDoseTallyMap doseTallyMap;
    for (EnergyDepositionEventMap::iterator mitr_event=edepEventMap.begin();mitr_event!=edepEventMap.end();mitr_event++)
    {
        for (int i=0;i<mitr_event->second.size();i++)
        {
            if (mitr_event->second[i].cellID!=-1) // Here we neglect processing ECM dose, only considering cell organelle
            {
                doseTallyMap[mitr_event->first][mitr_event->second[i].cellID][mitr_event->second[i].affectedCellOrganelle]=\
                doseTallyMap[mitr_event->first][mitr_event->second[i].cellID][mitr_event->second[i].affectedCellOrganelle]+mitr_event->second[i].edep;    
            }
            
        }
    }
    return doseTallyMap;


}
