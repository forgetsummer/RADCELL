// This class is for reading the radiation transport results so to do the further simulation for cell biology
// Two types of information will be read, the first one is dose information for cell, 
// The second type of information is DNA damage information for cell

#ifndef ReadRadiationTransportInfo_h 
#define ReadRadiationTransportInfo_h 1



#include <map>
#include <string>
using namespace std;

class ReadRadiationTransportInfo
{
public:
    void ReadDoseTallyOutPut(string doseFileName);
    void ReadDNADamageTallyOutPut(string DNAFileName);
    void ReadSingleCellDoseAsReference(string doseFileName);
    void ReadSingleCellDNADamageAsReference(string DNAFileName);
    void GetAbsoluteResultsByCOIMethod(double pDose, string directOrNot,string refDoseType);// using cell of interest method to get absolute results
    void GetAbsoluteResultsByFluenceMethod(double pFluence, string directorNot); // using fluence method to get absolute results
    double GetMCDoseOfCell(int cellID);// 
    double GetMCDoseStdOfCell(int cellID);
    double GetMCDSBOfCell(int cellID);
    double GetMCDSBStdOfCell(int cellID);
    double GetAbsDoseOfCell(int cellID);
    double GetAbsDoseStdOfCell(int cellID);
    double GetAbsDSBOfCell(int cellID);
    double GetAbsDSBStdOfCell(int cellID);
    double GetDoseFractionOfCell(int cellID);
//     ReadRadiationTransportInfo()
//     {
//         directResults = false;
//         indirectResults = false;
//         referenceDoseType = "max";// set up the referenceDose, could be min_dose or mean_dose
//         DNABasePairNum = 6; // for human genome, 6Gbp DNA base pair
//     }
    
private:
    double COIDose;
    double totalFluence;
    bool directResults  = false ;
    bool indirectResults = false;
    string referenceDoseType = "max";// set up the referenceDose, could be min_dose or mean_dose
    double DNABasePairNum = 6 ; // for human genome, 6Gbp DNA base pair
    bool COIMethod= false;// if use COI method , this marker is true
    bool FluenceMethod = false; // if use fluence method, this marker is true
    double GetMinimumCellDose();
    double GetMeanCellDose();
    double GetMaximumCellDose();
    void GetCellDoseFractionMap();
    std::map<int, double> cellDoseMap;
    std::map<int, double> cellDoseStdMap;
    std::map<int, double> cellDSBMap;
    std::map<int, double> cellDSBStdMap;
    std::map<int, double> cellDoseFractionMap;
    double singleCellDose;
    double singleCellDoseStd;
    double singleCellDSB;
    double singleCellDSBStd;
   
    
};





#endif
