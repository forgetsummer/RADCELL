#include "ReadRadiationTransportInfo.hh"

#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;


void ReadRadiationTransportInfo::ReadDoseTallyOutPut(string doseFileName)
{
    ifstream file (doseFileName.c_str());//open edep file
    std::string line;
    int line_num=0;
    
    while (std::getline(file,line))
    {
    
        if (line_num>=1) // just starting read file from third line to end
        {
            std::istringstream s(line);
            std::string field;
            int element_num=0;
            int cellID;
            double totalDose;
            double doseStd;
            
            while (getline(s,field,','))
            {
                if (element_num==0)
                {
                    istringstream buffer(field);
                    buffer>>cellID; // get cellID
                }
                        
                if (element_num==1)
                {
                    istringstream buffer(field);
                    buffer>>totalDose; // get cell dose
                }
                    
                if (element_num==2)
                {
                    istringstream buffer (field);
                    buffer>>doseStd; // get dose std
                }
                
                element_num=element_num+1;
            }
            
            cellDoseMap[cellID]=totalDose;
            cellDoseStdMap[cellID]=doseStd;
        }
        
        line_num=line_num+1;
    }

}

    
void ReadRadiationTransportInfo::ReadDNADamageTallyOutPut(string DNAFileName)
{
    ifstream file (DNAFileName.c_str());//open edep file
    std::string line;
    int line_num=0; 
    
    while (std::getline(file,line))
    {
    
        if (line_num>=1) // just starting read file from third line to end
        {
            std::istringstream s(line);
            std::string field;
            int element_num=0;
            int cellID;
            double DSB;
            double DSBStd;
            
            while (getline(s,field,','))
            {
                if (element_num==0)
                {
                    istringstream buffer(field);
                    buffer>>cellID; // get cellID
                }
                        
                if (element_num==3)
                {
                    istringstream buffer(field);
                    buffer>>DSB; // get DSB number
                }
                    
                if (element_num==4)
                {
                    istringstream buffer (field);
                    buffer>>DSBStd; // get DSB std;
                }
                
                element_num=element_num+1;
            }
            
            cellDSBMap[cellID]=DSB;
            cellDSBStdMap[cellID]=DSBStd;

        }
        
        line_num=line_num+1;

    }
    
}

void ReadRadiationTransportInfo::ReadSingleCellDoseAsReference(string doseFileName)
{
    ifstream file (doseFileName.c_str());//open edep file
    std::string line;
    int line_num=0;
    
    while (std::getline(file,line))
    {
    
        if (line_num>=1) // just starting read file from third line to end
        {
            std::istringstream s(line);
            std::string field;
            int element_num=0;
            int cellID;
            double totalDose;
            double doseStd;
            
            while (getline(s,field,','))
            {
                if (element_num==0)
                {
                    istringstream buffer(field);
                    buffer>>cellID; // get cellID
                }
                        
                if (element_num==1)
                {
                    istringstream buffer(field);
                    buffer>>totalDose; // get cell dose
                }
                    
                if (element_num==2)
                {
                    istringstream buffer (field);
                    buffer>>doseStd; // get dose std
                }
                
                element_num=element_num+1;
            }
            
            singleCellDose = totalDose;
            singleCellDoseStd = doseStd;
        }
        
        line_num=line_num+1;
    }


}

void ReadRadiationTransportInfo::ReadSingleCellDNADamageAsReference(string DNAFileName)
{
    ifstream file (DNAFileName.c_str());//open edep file
    std::string line;
    int line_num=0; 
    
    while (std::getline(file,line))
    {
    
        if (line_num>=1) // just starting read file from third line to end
        {
            std::istringstream s(line);
            std::string field;
            int element_num=0;
            int cellID;
            double DSB;
            double DSBStd;
            
            while (getline(s,field,','))
            {
                if (element_num==0)
                {
                    istringstream buffer(field);
                    buffer>>cellID; // get cellID
                }
                        
                if (element_num==3)
                {
                    istringstream buffer(field);
                    buffer>>DSB; // get DSB number
                }
                    
                if (element_num==4)
                {
                    istringstream buffer (field);
                    buffer>>DSBStd; // get DSB std;
                }
                
                element_num=element_num+1;
            }
            
            singleCellDSB = DSB;
            singleCellDSBStd = DSBStd;

        }
        
        line_num=line_num+1;

    }
}


void ReadRadiationTransportInfo::GetCellDoseFractionMap()
{
    for ( std::map<int, double>::iterator mitr_cell = cellDoseMap.begin();mitr_cell!=cellDoseMap.end();mitr_cell++)
    {
        if (referenceDoseType =="mean")
        {
            double mean_dose = GetMeanCellDose();
            cellDoseFractionMap[mitr_cell->first] = mitr_cell->second/mean_dose;
        }
        if (referenceDoseType == "min")
        {
            double min_dose = GetMinimumCellDose();
            cellDoseFractionMap[mitr_cell->first] = mitr_cell->second/min_dose;
        }
        if (referenceDoseType == "max")
        {
            double max_dose = GetMaximumCellDose();
            cellDoseFractionMap[mitr_cell->first] = mitr_cell->second/max_dose;
        }
    }

}



double ReadRadiationTransportInfo::GetMCDoseOfCell(int cellID)
{
    // this function is getting the cell dose, the algorithm is searching the cell dose map by the input cellID
    // if there is cellID, then just return the dose, if there is no, then just return dose as zero.
    
    if (cellDoseMap.find(cellID)==cellDoseMap.end()) // if cellID is not found
    {
        return 0;
    }
    else
    {
        return cellDoseMap.at(cellID);
    }

}

double ReadRadiationTransportInfo::GetMCDoseStdOfCell(int cellID)
{
    if (cellDoseStdMap.find(cellID)==cellDoseStdMap.end())
    {
        return 0;
    }
    else
    {
        return cellDoseStdMap.at(cellID);
    }

}


double ReadRadiationTransportInfo::GetMCDSBOfCell(int cellID)
{
    if (cellDSBMap.find(cellID)==cellDSBMap.end())
    {
        return 0;
    }
    else
    {
        return cellDSBMap.at(cellID);
    }

}
double ReadRadiationTransportInfo::GetMCDSBStdOfCell(int cellID)
{
    if (cellDSBStdMap.find(cellID)==cellDSBStdMap.end())
    {
        return 0;
    }
    else
    {
        return cellDSBStdMap.at(cellID);
    }

}

void ReadRadiationTransportInfo::GetAbsoluteResultsByCOIMethod(double pDose, string directOrNot,string refDoseType)
{
    COIDose = pDose;
    COIMethod = true;
    FluenceMethod = false;
    referenceDoseType = refDoseType;
    if (directOrNot == "direct")
    {
        directResults = true;
    }
    if (directOrNot == "indirect")
    {
        indirectResults = true;
        GetCellDoseFractionMap();
    }
}
void ReadRadiationTransportInfo::GetAbsoluteResultsByFluenceMethod(double pFluence, string directOrNot)
{
    totalFluence = pFluence;
    FluenceMethod = true;
    COIMethod = false;
    if (directOrNot == "direct")
    {
        directResults = true;
    }
    if (directOrNot == "indirect")
    {
        indirectResults = true;
    }
}

double ReadRadiationTransportInfo::GetAbsDoseOfCell(int cellID)
{
    if (cellDoseMap.find(cellID)==cellDoseMap.end()) // if cellID is not found
    {
        return 0;
    }
    else
    {
        double absDose;
        if (directResults == true)
        {
            if (COIMethod == true)// if using COI method 
            {
                double min_dose = GetMinimumCellDose();
                double mean_dose = GetMeanCellDose();
                double max_dose = GetMaximumCellDose();
    //             cout<<"min_dose is "<<min_dose<<endl;
    //             cout<<"mean_dose is "<<mean_dose<<endl;
    //             cout<<"MC dose is "<<GetMCDoseOfCell(cellID)<<endl;
                if (referenceDoseType == "mean")
                {
                    absDose = GetMCDoseOfCell(cellID)/mean_dose*COIDose;
                }
                if (referenceDoseType == "min")
                {
                    absDose = GetMCDoseOfCell(cellID)/min_dose*COIDose;
                }
                if (referenceDoseType == "max")
                {
                    absDose = GetMCDoseOfCell(cellID)/max_dose*COIDose;
                }
                
            }
            if (FluenceMethod == true) // if using fluence method
            {
                absDose = GetMCDoseOfCell(cellID)*totalFluence;
            }  
        }
        if (indirectResults == true)
        {
//             cout<<"using indirectResults"<<endl;
            if (COIMethod == true)
            {
                absDose = cellDoseFractionMap.at(cellID)*COIDose;
            }

            if (FluenceMethod == true)
            {
                absDose = cellDoseFractionMap.at(cellID)*singleCellDose*totalFluence;
            }
        }

        return absDose;
    }
}

double ReadRadiationTransportInfo::GetAbsDoseStdOfCell(int cellID)
{
    if (cellDoseStdMap.find(cellID)==cellDoseStdMap.end())
    {
        return 0;
    }
    else
    {
        double absDoseStd;
        if (directResults == true)
        {
            if (COIMethod == true)
            {
                double min_dose = GetMinimumCellDose();
                double mean_dose = GetMeanCellDose();
                double max_dose = GetMaximumCellDose();
                if (referenceDoseType == "mean")
                {
                    absDoseStd = GetMCDoseStdOfCell(cellID)*COIDose/mean_dose;
                }
                if (referenceDoseType == "min")
                {
                    absDoseStd = GetMCDoseStdOfCell(cellID)*COIDose/min_dose;
                }
                if (referenceDoseType == "max")
                {
                    absDoseStd = GetMCDoseStdOfCell(cellID)*COIDose/max_dose;
                }
                
            }
            if (FluenceMethod == true)
            {
                absDoseStd = GetMCDoseStdOfCell(cellID)*totalFluence;
            }
            
        }
        if (indirectResults == true)
        {
            if (COIMethod == true)
            {
//                 cout<<"singleCellDose "<<singleCellDose<<endl;
                absDoseStd = cellDoseFractionMap.at(cellID)*COIDose/singleCellDose*singleCellDoseStd;
            }
            if (FluenceMethod == true)
            {
                absDoseStd = cellDoseFractionMap.at(cellID)*singleCellDoseStd*totalFluence;
            }
        }

        return  absDoseStd;
    }

}

double ReadRadiationTransportInfo::GetAbsDSBOfCell(int cellID)
{
    double absDSB;
    if (directResults == true)
    {
        if (cellDSBMap.find(cellID)==cellDSBMap.end())// when using direct method, just check whether it contains DSB information
        {
            absDSB = 0;
        }
        else
        {
            if (COIMethod == true)
            {
                double min_dose = GetMinimumCellDose();
                double mean_dose = GetMeanCellDose();
                double max_dose = GetMaximumCellDose();
                double N_particle;
                if (referenceDoseType == "mean")
                {
                    N_particle = COIDose/mean_dose;
                }
                if (referenceDoseType == "min")
                {
                    N_particle = COIDose/min_dose;
                }
                if (referenceDoseType == "max")
                {
                    N_particle = COIDose/max_dose;
                }
                absDSB = GetMCDSBOfCell(cellID)*N_particle;
            }
            if (FluenceMethod == true)
            {
                absDSB = GetMCDSBOfCell(cellID)*totalFluence;
            }
        }
        
    }
    if (indirectResults == true)
    {
        if (cellDoseMap.find(cellID)==cellDoseMap.end()) // when use indirect method, then check whether it contains dose information
        {
            absDSB = 0;
        }
        else
        {
            if (COIMethod == true)
            {
                absDSB = cellDoseFractionMap.at(cellID)*COIDose/singleCellDose*singleCellDSB;
//                 cout<<"absDSB is "<<absDSB<< endl;
            }
            if (FluenceMethod == true)
            {
                absDSB = cellDoseFractionMap.at(cellID)*singleCellDSB*totalFluence;
            }
        }
        
    }
    return absDSB;
}
double ReadRadiationTransportInfo::GetAbsDSBStdOfCell(int cellID)
{
    double absDSBStd;
    if (directResults == true)
    {
        if (cellDSBMap.find(cellID)==cellDSBMap.end())
        {
            absDSBStd = 0;
        }
        else
        {
            if (COIMethod == true)
            {
                double min_dose = GetMinimumCellDose();
                double mean_dose = GetMeanCellDose();
                double max_dose = GetMaximumCellDose();
                double N_particle;
                if (referenceDoseType == "mean")
                {
                    N_particle = COIDose/mean_dose;
                }
                if (referenceDoseType == "min")
                {
                    N_particle = COIDose/min_dose;
                }
                if (referenceDoseType == "max")
                {
                    N_particle = COIDose/max_dose;
                }
                absDSBStd = GetMCDSBStdOfCell(cellID)*N_particle;
            }
            if (FluenceMethod == true)
            {
                absDSBStd = GetMCDSBStdOfCell(cellID)*totalFluence;
            }
        }
    }
    if (indirectResults == true)
    {
        if (cellDoseMap.find(cellID)==cellDoseMap.end())
        {
            absDSBStd = 0;
        }
        else
        {
            if (COIMethod == true)
            {
                absDSBStd = cellDoseFractionMap.at(cellID)*COIDose/singleCellDose*singleCellDSBStd;
            }
            if (FluenceMethod == true)
            {
                absDSBStd = cellDoseFractionMap.at(cellID)*singleCellDSBStd*totalFluence;
            }
        }

    }
    return absDSBStd;
}


double ReadRadiationTransportInfo::GetMinimumCellDose()
{
    double min_dose;
    int cellID_minDose = 0;
    std::map<int, double>::iterator mitr_cell = cellDoseMap.begin();
    
    min_dose = mitr_cell->second;
    for (std::map<int, double>::iterator mitr_cell=cellDoseMap.begin();mitr_cell!=cellDoseMap.end();mitr_cell++)
    {
        if (mitr_cell->second<=min_dose)
        {
            min_dose = mitr_cell->second;
            cellID_minDose = mitr_cell->first;
        }
    }
    
    return min_dose; // return the smallest cell dose 
}

double ReadRadiationTransportInfo::GetMaximumCellDose()
{
    double maxCellDose;
    int cellID_maxDose = 0;
    std::map<int, double>::iterator mitr_cell = cellDoseMap.begin();
    maxCellDose = mitr_cell->second;
    for (std::map<int, double>::iterator mitr_cell=cellDoseMap.begin();mitr_cell!=cellDoseMap.end();mitr_cell++)
    {
        if (mitr_cell->second>=maxCellDose)
        {
            maxCellDose = mitr_cell->second;
            cellID_maxDose = mitr_cell->first;
        }
        
    }
    return maxCellDose; // return the maximum cell dose
}

double ReadRadiationTransportInfo::GetMeanCellDose()
{
    double meanCellDose;
    double sumDose = 0;
    for (std::map<int, double>::iterator mitr_cell=cellDoseMap.begin();mitr_cell!=cellDoseMap.end();mitr_cell++)
    {
        sumDose = sumDose + mitr_cell->second;
    }
    meanCellDose = sumDose/cellDoseMap.size();
    return meanCellDose;
}


double ReadRadiationTransportInfo::GetDoseFractionOfCell(int cellID)
{
    GetCellDoseFractionMap();
    if (cellDoseFractionMap.find(cellID)==cellDoseFractionMap.end()) // if cellID is not found
    {
        return 0;
    }
    else
    {
        return cellDoseFractionMap.at(cellID);
    }
    
}



