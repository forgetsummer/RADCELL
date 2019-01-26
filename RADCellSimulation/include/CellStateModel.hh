#ifndef CellStateModel_h
#define CellStateModel_h 1
#include <map>
#include <vector>
#include <string>
using namespace std;

class CellStateModel
{
public:
    void CellStateModelParameterSetup(double alphaVal, double betaVal, double E1, double E2, double E3,double sigma,double fG1,double fS, double fG2, double fM, double f1, double lambda1, double f2, double lambda2);// function for setting up the simulation parameters for cell state model
    void CellInformationSetup(int cType, double tG1, double sigmaG1, double tS, double sigmaS,\
        double tG2, double sigmaG2, double tM, double sigmaM, int freeReceptor
    );// function for basic cell information set up, setup based on cell type 
    void TissueGeometryInitialization(double xDim, double yDim, double zDim,double gridSize);
    void CellPositionInitialization(int cellID, double cx, double cy, double cz);
    void CellPhaseInitialization(int cellID, int cellType, string cphase);// call this function for each cell to initialize cell phase
    void CellPhaseInitializationRandom(int cellID,int cellType); // To initialize cell phase randomly
    void CellStateInitialization(int cellID, string cState);// function initialize cell state, S1,S2,S3
    void CellPhaseUpdate(int cellID, int cellType, bool proliferationState, double deltaT, int frequency);// function updating the cell phase with time going on
    void CellStateUpdate(int cellID, int cellType, int DSBNum, double integralConcentration, double deltaT, int frequency); // function updating the cell state, deltaT in unit of seconds
    std::map<int, std::string> GetCellPhase();
    std::map<int, double> GetCellPhaseDuration();
    std::map<int, double> GetCellAge();
    std::map<int, std::string> GetCellState();
    std::map<int, double> GetCellStateDuration();
    std::map<int, double> GetCellStateAge();
    std::map<int, int> GetCellAncestryID();
    std::map<int, double> GetCellPositionX();
    std::map<int, double> GetCellPositionY();
    std::map<int, double> GetCellPositionZ();
    void SetUpContactInhibition(bool considerOrNot); // determine whether considering contact inhibition or not, set up true when consider it

    

private:
    double d;// grid size of cell home
    int NX;
    int NY;
    int NZ;
    double alpha;// consant for translating DSB numbers to state energy, E = alpha*DSBs
    double beta;// constant for translating integral of concentration to state energy, E = beta*C
    double sigma;
    double E1;
    double E2;
    double E3;
    double fG1;
    double fS;
    double fG2;
    double fM;
    double f1;
    double f2;
    double lambda1;
    double lambda2;
    void CellPhaseTransition(int cellID, int cellType, double deltaT, int frequency);//function updating the cell phase with time going on
    void CellStateTransition(int cellID, double E1,double E21, double E3,double increaseTime, double ObservationTime,bool transitionType);
    double GaussianSampling(double mu,double sigma);  
    double GaussianSampling2Pi(double mu, double sigma);
    double GaussianCDF(double x, double mu, double sigma);
    double InstantaneousStateJumpProb(double p_sp,double increaseTime,double ObservationTime);
    double DelayedStateJumpProb(double p_sp,double increaseTime,double ObservationTime);
    bool considerContactInhibition = true;
    
    struct PhaseInfo // a struct containing cell age information
    {
        std::string phase; // age phase, {G1, S, G2, M}
        double age; // time cell has been staying at this phase
        double duration; // expected time the cell will stay at this phase, determined by GaussianSampling
        double telomereFraction;
        int ancestryID;// the ancestry id for the daughter cell, daughter cell will inherite mother ancestry ID
    };
    struct StateInfo
    {
        std::string state;// state , {S1,S2,S3}
        double age;
        double duration;
    };

    struct CycleInfo
    {
        double mTG1; // mean time for staying at G1
        double sigmaTG1; // sigma for G1
        double mTS;
        double sigmaTS;
        double mTG2;
        double sigmaTG2;
        double mTM;
        double sigmaTM;
    };
    struct PositionInfo
    {
         int i;
         int j;
         int k;
    };
    std::map<std::string,PositionInfo>CheckCellContactInhibitionCondition(int ci, int cj, int ck);
    std::vector<std::vector<std::vector<int> > > cellHomeResidence;
    std::map<int,PositionInfo > cellPositionMap;//key is cell id
    std::map<int, PhaseInfo> cellPhaseMap; // key is cell id
    std::map<int, StateInfo> cellStateMap; // map for storing the cell state, {S1, S2, S3} 
    std::map<int, double> cellStateObservationTimeMap;// key is cell id, map for storing the total time step during the observation time [0,T]
    std::map<int, std::string> cellPossibleJumpToStateMap;//key is cell id
    std::map<int, double> cellAccumulatedJumpProbMap;// key is cell id
    std::map<int, int> freeReceptorMap; // key is cell type
    std::map<int, CycleInfo> cellCycleMap; // key is cell type
    std::map<int, bool> cellProliferativeMap;// key is cell id, value is proliferative state, true or not 
    std::map<int, double> cellExternalPerturbationEnergyMap;// map for storing the external perturbation energy, key is cell id

};

#endif