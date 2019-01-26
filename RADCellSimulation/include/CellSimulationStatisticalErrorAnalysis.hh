#ifndef CellSimulationStatisticalErrorAnalysis_h

#define CellSimulationStatisticalErrorAnalysis_h 1
#include <vector>
#include<math.h>
#include<iostream>
using namespace std;
class CellSimulationStatisticalErrorAnalysis
{
public:
    template <class T>
    void SimulationStatisticalAnalysis(std::map<int, std::map<int, T> > valueMap, int totalEventNumber);
    std::map<int, double> GetSimulationMean(){return meanValue;}; // function returning mean value of simulation for cell
    std::map<int, double> GetSimulationStd(){return stdValue;}; // function returning std value of simulation for cell
private:
    std::map<int, double> meanValue;
    std::map<int, double> stdValue;
    
};


template <class T>
void CellSimulationStatisticalErrorAnalysis::SimulationStatisticalAnalysis(std::map<int, std::map<int, T> > valueMap, int totalEventNumber)

{
    std::map<int, T> sum_cell_value;
    std::map<int, T> sumSquared_cell_value;
    std::map<int, double>meanSquaredValue;
    double confidence_level_z;
    for (typename std::map<int, std::map<int, T> >::iterator mitr_event=valueMap.begin();mitr_event!=valueMap.end();mitr_event++) //loop eventID
    {
        for (typename std::map<int, T>::iterator mitr_cell=(mitr_event->second).begin();mitr_cell!=(mitr_event->second).end();mitr_cell++) // loop cellID
        {
//             if (mitr_cell->first==0)
//             {
//             cout<<"The tally information by event "<<mitr_event->first<<" to cell "<<mitr_cell->first<<" is "<<valueMap[mitr_event->first][mitr_cell->first]<<endl;
//                 
//             }
            sum_cell_value[mitr_cell->first]=sum_cell_value[mitr_cell->first]+(valueMap[mitr_event->first])[mitr_cell->first];
            sumSquared_cell_value[mitr_cell->first]=sumSquared_cell_value[mitr_cell->first]+((valueMap[mitr_event->first])[mitr_cell->first])*((valueMap[mitr_event->first])[mitr_cell->first]);
        }
    }
    for (typename std::map<int, T>::iterator mitr_cell=sum_cell_value.begin();mitr_cell!=sum_cell_value.end();mitr_cell++)
    {
        meanValue[mitr_cell->first]= mitr_cell->second/double(totalEventNumber);
        
    }
    confidence_level_z=1.96;
    for (typename std::map<int, T>::iterator mitr_cell=sumSquared_cell_value.begin();mitr_cell!=sumSquared_cell_value.end();mitr_cell++)
    {
        meanSquaredValue[mitr_cell->first]=mitr_cell->second/double(totalEventNumber);
        stdValue[mitr_cell->first]=sqrt(meanSquaredValue[mitr_cell->first]-(meanValue[mitr_cell->first])*(meanValue[mitr_cell->first]))*confidence_level_z/sqrt(totalEventNumber);
        
    }
    
}

#endif