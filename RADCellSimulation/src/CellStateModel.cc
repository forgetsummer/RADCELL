#include "CellStateModel.hh"
#include <stdlib.h>
#include <cmath>
#include "Randomize.hh"  // uisng random number
#include <iostream>
#include <stdlib.h>

void CellStateModel::CellInformationSetup(int cType, double tG1, double sigmaG1, double tS, double sigmaS, double tG2, double sigmaG2, double tM, double sigmaM, int freeReceptor)
{
    cellCycleMap[cType].mTG1 = tG1;
    cellCycleMap[cType].sigmaTG1 = sigmaG1;
    cellCycleMap[cType].mTS = tS;
    cellCycleMap[cType].sigmaTS = sigmaS;
    cellCycleMap[cType].mTG2 = tG2;
    cellCycleMap[cType].sigmaTG2 = sigmaG2;
    cellCycleMap[cType].mTM = tM;
    cellCycleMap[cType].sigmaTM = sigmaM;
    freeReceptorMap[cType]= freeReceptor;
    
}

double CellStateModel::GaussianSampling(double mu, double sigma)
{
    double Xf;
    double pi = 3.1415926535897;
    Xf =mu+ sigma*sqrt(-2*log(G4UniformRand()))*cos(pi/2.0*G4UniformRand());
    return Xf;
    
}
double CellStateModel::GaussianSampling2Pi(double mu, double sigma)
{
    double Xf;
    double pi = 3.1415926535897;
    Xf =mu+ sigma*sqrt(-2*log(G4UniformRand()))*cos(2.0*pi*G4UniformRand());
    return Xf;

}


double CellStateModel::GaussianCDF(double x, double mu, double sigma)
{
    double x_stdNORM=(x-mu)/sigma;
    double p_CDF;
    p_CDF = 0.5*(1+erf(x_stdNORM/sqrt(2)));
    return p_CDF;// return the probability of N(x,mu, sigma)

}

void CellStateModel::TissueGeometryInitialization(double xDim, double yDim, double zDim, double gridSize)
{
    d = gridSize;
    double cellHomeSizeX = gridSize;
    double cellHomeSizeY = gridSize;
    double cellHomeSizeZ = gridSize;
    int N_X=ceil(xDim/cellHomeSizeX);
    int N_Y=ceil(yDim/cellHomeSizeY);
    int N_Z=ceil(zDim/cellHomeSizeZ);
    if (N_Z==0)// in case of  ceil(0)=0
    {
        N_Z=1;
    }
    if (N_Y==0)
    {
        N_Y=1;
    }
    if (N_X==0)
    {
        N_X=1;
    }

    NX = N_X;
    NY = N_Y;
    NZ = N_Z;
    
    for (int i=0;i<N_X;i++)
    {
        cellHomeResidence.push_back(std::vector<std::vector<int> >());
        for (int j=0;j<N_Y;j++)
        {
            (cellHomeResidence[i]).push_back(std::vector<int>());
            for (int k=0;k<N_Z;k++)
            {
                (cellHomeResidence[i][j]).push_back(0);// intialize cellHomeResidence
            }
        }
    }

}


void CellStateModel::CellPositionInitialization(int cellID, double cx, double cy, double cz)
{
//     cout<<"cx ="<<cx<<endl;
//     cout<<"cy ="<<cy<<endl;
//     cout<<"cz ="<<cz<<endl;
    double c_i = floor(cx/d)-1;
    double c_j = floor(cy/d)-1;
    double c_k = floor(cz/d)-1;            
    if (c_i==-1)
    {
        c_i=0;
    }
    if (c_j==-1)
    {
        c_j=0;
    }
    if (c_k==-1)
    {
        c_k=0;
    }
    cellHomeResidence[c_i][c_j][c_k] = 1; // this cell home was occupied by the cell
    PositionInfo cellPosition;
    cellPosition.i = c_i;
    cellPosition.j = c_j;
    cellPosition.k = c_k;
    cellPositionMap[cellID] = cellPosition;

}


void CellStateModel::CellPhaseInitialization(int cellID,int cellType, string cphase)
{
    PhaseInfo cellPhase;
    if (cellCycleMap.find(cellType)==cellCycleMap.end())
    {
        cout<<endl<<"ERROR: "<<"THE CELL TYPE IS NOT EXISTING!"<<endl;
        exit(1);
    }

    
    if ((cphase!="G1") && ( cphase!="S") && (cphase!="G2")&& (cphase!="M") && (cphase!="G0"))
    {
        cout<<cphase<<endl;
        cout<<endl<<"ERROR: "<<"THIS PHASE DOES NOT EXIST!"<<endl;
        exit(1);
    }
    else
    {
        cellPhase.phase = cphase;
        cellPhase.age = 0;
        cellPhase.telomereFraction=1;// initially telomereFraction is 1
        cellPhase.ancestryID = cellID;
        if (cphase == "G1")
        {
            cellPhase.duration = GaussianSampling(cellCycleMap.at(cellType).mTG1,cellCycleMap.at(cellType).sigmaTG1);
        }
        if (cphase =="S")
        {
            cellPhase.duration = GaussianSampling(cellCycleMap.at(cellType).mTS,cellCycleMap.at(cellType).sigmaTS);
        }
        if (cphase == "G2")
        {
            cellPhase.duration = GaussianSampling(cellCycleMap.at(cellType).mTG2,cellCycleMap.at(cellType).sigmaTG2);
        }
        if (cphase == "M")
        {
            cellPhase.duration = GaussianSampling(cellCycleMap.at(cellType).mTM,cellCycleMap.at(cellType).sigmaTM);
        }
        if (cphase == "G0")
        {
            cellPhase.duration = 1E+18;// just initially give a very large number which means that cell will not proliferate
        }
    }
    
    cellPhaseMap[cellID] = cellPhase;
    
}


void CellStateModel::CellPhaseInitializationRandom(int cellID, int cellType)
{
    double eta = G4UniformRand();
    string phase;
    if (eta<=0.2)
    {
        phase = "G1";
    }
    if (0.2<eta && eta<=0.4)
    {
        phase = "G2";
    }
    if (0.4< eta && eta<=0.6)
    {
        phase = "S";
    }
    if (0.6<eta && eta <=0.8)
    {
        phase = "M";
    }
    if (0.8<eta && eta<=1.0)
    {
        phase = "G0";
    }
    PhaseInfo cellPhase;
    cellPhase.age = 0;
    cellPhase.phase = phase;
    cellPhase.telomereFraction = 1;
    cellPhase.ancestryID = cellID;
    
    if (phase == "G1")
    {
//         cout<<"meam time G1 "<<cellCycleMap.at(cellType).mTG1<<" sigma time "<<cellCycleMap.at(cellType).sigmaTG1<<endl;
        cellPhase.duration = GaussianSampling(cellCycleMap.at(cellType).mTG1,cellCycleMap.at(cellType).sigmaTG1);
    }
    if (phase =="S")
    {
//         cout<<"meam time S "<<cellCycleMap.at(cellType).mTS<<" sigma time "<<cellCycleMap.at(cellType).sigmaTS<<endl;
        cellPhase.duration = GaussianSampling(cellCycleMap.at(cellType).mTS,cellCycleMap.at(cellType).sigmaTS);
    }
    if (phase == "G2")
    {
//         cout<<"meam time G2 "<<cellCycleMap.at(cellType).mTG2<<" sigma time "<<cellCycleMap.at(cellType).sigmaTG2<<endl;
        cellPhase.duration = GaussianSampling(cellCycleMap.at(cellType).mTG2,cellCycleMap.at(cellType).sigmaTG2);
    }
    if (phase == "M")
    {
//         cout<<"meam time M "<<cellCycleMap.at(cellType).mTM<<" sigma time "<<cellCycleMap.at(cellType).sigmaTM<<endl;
        cellPhase.duration = GaussianSampling(cellCycleMap.at(cellType).mTM,cellCycleMap.at(cellType).sigmaTM);
    }
    if (phase == "G0")
    {
        cellPhase.duration = 1E+18;
    }
    
    cellPhaseMap[cellID] = cellPhase;

}

void CellStateModel::CellPhaseTransition(int cellID, int cellType, double deltaT, int frequency)
{
//     double increaseTime = deltaT*frequency;// change time unit from seconds to hours
    double increaseTime = deltaT*frequency/3600;// change time unit from seconds to hours
    if (cellProliferativeMap[cellID]==true) // if cell staying at proliferative state
    {
        if (cellPhaseMap[cellID].phase == "G0")
        {
//             cout<<"processing G0"<<endl;
            cellPhaseMap[cellID].phase = "G1"; // it goes to first phase when it becomes proliferative cell
            cellPhaseMap[cellID].age = 0;
            cellPhaseMap[cellID].duration = GaussianSampling(cellCycleMap.at(cellType).mTG1,cellCycleMap.at(cellType).sigmaTG1);
            return;
        }
        if (cellPhaseMap[cellID].phase == "G1")
        {
//             cout<<"processing G1"<<endl;
            cellPhaseMap[cellID].age = cellPhaseMap[cellID].age+increaseTime;// age increase at each time step
            if (cellPhaseMap[cellID].age>cellPhaseMap[cellID].duration)
            {
                cellPhaseMap[cellID].phase = "S";
                cellPhaseMap[cellID].age = 0;
                cellPhaseMap[cellID].duration = GaussianSampling(cellCycleMap.at(cellType).mTS,cellCycleMap.at(cellType).sigmaTS);
            }
            return;
        }
        if (cellPhaseMap[cellID].phase == "S")
        {
//             cout<<"processing S"<<endl;
            cellPhaseMap[cellID].age = cellPhaseMap[cellID].age+increaseTime;// age increase at each time step
            if (cellPhaseMap[cellID].age>cellPhaseMap[cellID].duration)
            {
                cellPhaseMap[cellID].phase = "G2";  
                cellPhaseMap[cellID].age = 0;
                cellPhaseMap[cellID].duration = GaussianSampling(cellCycleMap.at(cellType).mTG2,cellCycleMap.at(cellType).sigmaTG2);
            }
            return;
        }
        if (cellPhaseMap[cellID].phase == "G2")
        {
//             cout<<"processing G2"<<endl;
            cellPhaseMap[cellID].age = cellPhaseMap[cellID].age+increaseTime;// age increase at each time step 
            if (cellPhaseMap[cellID].age>cellPhaseMap[cellID].duration)
            {
                cellPhaseMap[cellID].phase = "M"; 
                cellPhaseMap[cellID].age = 0;
                cellPhaseMap[cellID].duration = GaussianSampling(cellCycleMap.at(cellType).mTM,cellCycleMap.at(cellType).sigmaTM);
            }
            return;
        }
        if (cellPhaseMap[cellID].phase == "M")
        {
//             cout<<"processing M"<<endl;
//             cout<<"The M Duration of Cell:"<<cellID<<" "<< cellPhaseMap[cellID].duration<<endl;
//             cout<<"the telomereFraction is "<<cellPhaseMap[cellID].telomereFraction<<endl;
            cellPhaseMap[cellID].age = cellPhaseMap[cellID].age+increaseTime;// age increase at each time step  
            bool allowableInCellAge = cellPhaseMap[cellID].age>cellPhaseMap[cellID].duration;
            bool allowableInCellTelomereLength = cellPhaseMap[cellID].telomereFraction>8.88E-16;
            
            std::map<std::string,PositionInfo> possibleCellHomeForMitosis;// key is the possible cell home location
            int i = cellPositionMap[cellID].i;
            int j = cellPositionMap[cellID].j;
            int k = cellPositionMap[cellID].k;
//             cout<<"cellID is "<<cellID<<" The i j k are "<<i<<" , "<<j<<" , "<<k<<endl;
            bool allowableInSpace;
//             considerContactInhibition = true; // if consider contact inhibition set this as true, otherwise false
//             bool considerContactInhibition = false; // if consider contact inhibition set this as true, otherwise false
            if (considerContactInhibition)
            {
                possibleCellHomeForMitosis = CheckCellContactInhibitionCondition(i,j,k);
                 allowableInSpace = possibleCellHomeForMitosis.size()>0;
            }
            else
            {
                possibleCellHomeForMitosis["original"] = cellPositionMap.at(cellID);// This is for case without conisidering contact inhibition 
                allowableInSpace = 1;
            }
            
//             cout<<"cellID is "<<cellID<<" allowableInSpace is "<<allowableInSpace<<endl;

            if (allowableInCellAge && allowableInCellTelomereLength && allowableInSpace) // if meeting those three condition then mitosis
            {
                int cellID_max = 0;
                for (std::map<int, PhaseInfo>::iterator mitr_cell1 = cellPhaseMap.begin();mitr_cell1!= cellPhaseMap.end();mitr_cell1++)
                {
                    if (cellID_max<mitr_cell1->first)
                    {
                        cellID_max = mitr_cell1->first;
                    }
                }

                cellPhaseMap[cellID_max+1].phase = "G0";// creat two new daughter cells, staying at G0 phase
                cellPhaseMap[cellID_max+1].age = 0;
                cellPhaseMap[cellID_max+1].duration = 1E+18;
                cellPhaseMap[cellID_max+1].telomereFraction = 0.5*cellPhaseMap[cellID].telomereFraction;// telomereFraction reduce half
                cellPhaseMap[cellID_max+1].ancestryID = cellPhaseMap[cellID].ancestryID;// inherete mother ancestryID
                cellStateMap[cellID_max+1].state = "S1";// assign a healthy state for new born cell
                cellStateMap[cellID_max+1].age = 0;
                cellStateMap[cellID_max+1].duration = 1E+18;
                cellPositionMap[cellID_max+1] = cellPositionMap[cellID]; // the first daughter cell takes the mother's place
                cellPhaseMap[cellID_max+2].phase = "G0";// creat two new daughter cells, staying at G0 phase
                cellPhaseMap[cellID_max+2].age = 0;
                cellPhaseMap[cellID_max+2].duration = 1E+18; 
                cellPhaseMap[cellID_max+2].telomereFraction = 0.5*cellPhaseMap[cellID].telomereFraction;
                cellPhaseMap[cellID_max+2].ancestryID = cellPhaseMap[cellID].ancestryID;// inherete mother ancestryID
                cellStateMap[cellID_max+2].state = "S1";// assign a healthy state for new born cell
                cellStateMap[cellID_max+2].age = 0;
                cellStateMap[cellID_max+2].duration = 1E+18;
                int ii  = floor(G4UniformRand()*possibleCellHomeForMitosis.size());
                int pi = 0;
                if (pi ==ii)
                {
                    cellPositionMap[cellID_max+2] = possibleCellHomeForMitosis[possibleCellHomeForMitosis.begin()->first];
                    int ci = possibleCellHomeForMitosis.begin()->second.i;
                    int cj = possibleCellHomeForMitosis.begin()->second.j;
                    int ck = possibleCellHomeForMitosis.begin()->second.k;
                    cellHomeResidence[ci][cj][ck] = 1;
                }
                else
                {
                    for (std::map<std::string, PositionInfo>::iterator mitr = possibleCellHomeForMitosis.begin(); mitr!= possibleCellHomeForMitosis.end();mitr++)
                    {
                        pi = pi+1;
                        if (pi==ii+1)
                        {
                            cellPositionMap[cellID_max+2] = mitr->second; //randomly select a possible position for the second daughter 
//                             cout<<"The possible room for daughter cell is "<<mitr->first<<endl;
//                             cout<<"dimension is "<<NX<<" "<<NY<<" "<<NZ<<endl;
//                             cout<<"i j k for this room is "<<mitr->second.i<<" "<<mitr->second.j<<" "<<mitr->second.k<<endl;
                            cellHomeResidence[mitr->second.i][mitr->second.j][mitr->second.k] =1; // this empty place is occupied now
                        }
                    }
                }
//                 cout<<"the mitosis cell id is "<<cellID<<endl;
//                 cout<<"the cell id before delete old cell are "<<endl;
//                 for (std::map<int, PhaseInfo>::iterator mitr_cell=cellPhaseMap.begin();mitr_cell!=cellPhaseMap.end();mitr_cell++)
//                 {
//                     cout<<"cell id is "<<mitr_cell->first<<endl;
//                 }
                cellPhaseMap.erase(cellID); // delete old cell in cell phase map
                cellStateMap.erase(cellID);//  delete old cell in cell state map
                cellPositionMap.erase(cellID);// delete old cell in cell position map
//                 cout<<"after delete old cell"<<endl;
//                 cout<<"the cell id after delete old cell are "<<endl;
//                 for (std::map<int, PhaseInfo>::iterator mitr_cell=cellPhaseMap.begin();mitr_cell!=cellPhaseMap.end();mitr_cell++)
//                 {
//                     cout<<"cell id is "<<mitr_cell->first<<endl;
//                 }

            }
            return;
        }
    }
}
map< string, CellStateModel::PositionInfo > CellStateModel::CheckCellContactInhibitionCondition(int i, int j, int k)
{
    std::map<std::string,PositionInfo> possibleCellHomeForMitosis;
    if (NX>1&&NY>1&&NZ==1)// if only in X-Y plane
    {
        if (i==0)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
        }
        if (i==NX-1)
        {
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (0<i && i<NX-1)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (j==0)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            
        }

        if (j==NY-1)
        {
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
        }
                                                    
        if (0<j && j<NY-1)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
            
        }
        
    }

    if (NY>1&&NZ>1&&NX==1)//if in Y-Z plane
    {
        if (j==0)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            
        }

        if (j==NY-1)
        {
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
        }
                                                    
        if (0<j && j<NY-1)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
            
        }
        if (k==0)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }
        }
        if (k==NZ-1)
        {
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        }
        if (0<k && k<NZ-1)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }  
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        } 
        
    }
    if (NX>1&&NZ>1&&NY==1)// if in X-Z plane
    {
            if (i==0)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
        }
        if (i==NX-1)
        {
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (0<i && i<NX-1)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (k==0)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }
        }
        if (k==NZ-1)
        {
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        }
        if (0<k && k<NZ-1)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }  
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        } 
        
    }
    if (NX>1&&NY>1&&NZ>1)// 3D case
    {
        if (i==0)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
        }
        if (i==NX-1)
        {
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (0<i && i<NX-1)
        {
            if (cellHomeResidence[i+1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i+1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["right"] = possiblePosition; 
            }
            if (cellHomeResidence[i-1][j][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i-1;
                possiblePosition.j = j;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["left"] = possiblePosition;
            }
        }
        if (j==0)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            
        }

        if (j==NY-1)
        {
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
        }
                                                    
        if (0<j && j<NY-1)
        {
            if (cellHomeResidence[i][j+1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j+1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["front"] = possiblePosition;
            }
            if (cellHomeResidence[i][j-1][k]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j-1;
                possiblePosition.k = k;
                possibleCellHomeForMitosis["back"] = possiblePosition;
            }
            
        }
        if (k==0)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }
        }
        if (k==NZ-1)
        {
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        }
        if (0<k && k<NZ-1)
        {
            if (cellHomeResidence[i][j][k+1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k+1;
                possibleCellHomeForMitosis["up"] = possiblePosition;
            }  
            if (cellHomeResidence[i][j][k-1]==0)
            {
                PositionInfo possiblePosition;
                possiblePosition.i = i;
                possiblePosition.j = j;
                possiblePosition.k = k-1;
                possibleCellHomeForMitosis["down"] = possiblePosition;
            }
        } 
        
    }
    return possibleCellHomeForMitosis;
    
}

void CellStateModel::CellPhaseUpdate(int cellID, int cellType, bool proliferationState, double deltaT, int frequency)
{
    bool proliferativeFromCellState;
    if(cellStateMap[cellID].state=="S1")// only when the cell is in healthy state, then it has the condtion to proliferate
    {
        proliferativeFromCellState = true;
    }
    else
    {
        proliferativeFromCellState = false;
    }
    
//     cellProliferativeMap[cellID] =(proliferationState);//when the nutrient condtion and healthy state all true, then it goes to proliferation
    cellProliferativeMap[cellID] =(proliferationState&&proliferativeFromCellState);//when the nutrient condtion and healthy state all true, then it goes to proliferation
    CellPhaseTransition(cellID,cellType,deltaT,frequency);
//     cout<<"CELL: "<<cellID<<" : "<<"AGE IS: "<<cellPhaseMap[cellID].age<<" PHASE IS: "<<cellPhaseMap[cellID].phase<<endl;

}

void CellStateModel::CellStateInitialization(int cellID, string cState)
{
    StateInfo cellState;
    cellState.state = cState;
    cellState.age =0;
    if ((cState!="S1") && ( cState!="S21")&&(cState!="S22") && (cState!="S3"))
    {
        cout<<cState<<endl;
        cout<<endl<<"ERROR: "<<"THIS STATE DOES NOT EXIST!"<<endl;
        cout<<"THE POSSIBLE STATES ARE: S1, S21, S22, S3"<<endl;
        exit(1);
    }
    
    if(cState=="S1")
    {
        cellState.duration=1E+18;
    }
    if(cState=="S21")
    {
        double E20= GaussianSampling2Pi(E2,sigma);// Sample a E20 of S2 state
        while (E20<0)
        {
            E20= GaussianSampling2Pi(E2,sigma);// Sample a E20 of S2 state
        }
        double pi = 3.1415926535897;
        double arrestedDuration = E20/sqrt(2*pi)/sigma*(f1/lambda1+f2/lambda2);
        cellState.duration = arrestedDuration;
    }
    if(cState=="S22")
    {
        cellState.duration = 1E+18;
    }
    if(cState=="S3")
    {
        cellState.duration = 1E+18;
    }

    cellStateMap[cellID] = cellState;
    cellExternalPerturbationEnergyMap[cellID] = 0;// at the beginning, there is no external perturbation energy
    cellStateObservationTimeMap[cellID] = 0;// at the begining, the total observation time is 0
    cellAccumulatedJumpProbMap[cellID] = 0;
    cellPossibleJumpToStateMap[cellID] = "S21";

}

void CellStateModel::CellStateUpdate(int cellID,int cellType, int DSBNum, double integralConcentration, double deltaT, int frequency)
{
    double increaseTime = deltaT*frequency/3600; // increase the time step for cell state transition, change time unit from second  to hour
    double ObservationTime = 0;
    ObservationTime = cellCycleMap[cellType].mTG1 +  cellCycleMap[cellType].mTG2 + cellCycleMap[cellType].mTS+ cellCycleMap[cellType].mTM;
    cellExternalPerturbationEnergyMap[cellID] = alpha*DSBNum + beta*integralConcentration;// calcuate the external perturbation energy for current time step
//     cout<<"The DSBNum is "<<DSBNum<<" , "<<"The bystander signal concentration is "<<integralConcentration<<endl;
//     cout<< "The external perturbation energy is "<<cellExternalPerturbationEnergyMap[cellID]<<endl;
    bool transitionType;
    if (DSBNum>0) // here we simply take radiation exposure as the agenet inducing Instantaneous cell state jump
    {
        transitionType = true;
    }
    else
    {
        transitionType = false;
    }
    if(cellPhaseMap[cellID].phase=="G1")//if it is G1 phase cell
    {
        double E21 = E2; // mean state energy of first degenerate state of S2
        double E22 = E3-sigma; // mean state energy of second degenerate state of S2 
        CellStateTransition(cellID,E1,E21,E3,increaseTime,ObservationTime,transitionType);
//         if(cellStateMap[cellID].state=="S1")// when the original state is S1
//         {
//             double ES1ToS2Transient = -fabs(E1+cellExternalPerturbationEnergyMap.at(cellID)-E21)/2/sigma;
//             double ES1ToS2Permanent = -fabs(E1+cellExternalPerturbationEnergyMap.at(cellID)-E22)/2/sigma;
//             double ES1ToS3 = -fabs(E1+cellExternalPerturbationEnergyMap.at(cellID)-E3)/2/sigma;
//             double PS1ToS2Transient = 2*GaussianCDF(ES1ToS2Transient,0,1);
//             double PS1ToS2Permanent = 2*GaussianCDF(ES1ToS2Permanent,0,1);
//             double PS1ToS3= 2*GaussianCDF(ES1ToS3,0,1);
//             if (E1+cellExternalPerturbationEnergyMap.at(cellID)-E3>=0)// boundary condition 
//             {
//                 PS1ToS3 = 1;
//             }
//             else
//             {
//                 PS1ToS3 = 2*GaussianCDF(ES1ToS3,0,1);
//             }
//             
//             double PPS1TOS2Transient = PS1ToS2Transient/(PS1ToS2Transient + PS1ToS2Permanent + PS1ToS3);
//             double PPS1TOS2Permanet = PS1ToS2Permanent/(PS1ToS2Transient + PS1ToS2Permanent + PS1ToS3);
//             double PPS1TOS3 = PS1ToS3/(PS1ToS2Transient + PS1ToS2Permanent + PS1ToS3);
//             
//             std::string possibleJumpState;
//             double rand = G4UniformRand();
//             
// 
//             if (rand <= PPS1TOS2Transient)
//             {
//                 possibleJumpState = "S21";
//             }
//             else if (rand<= PPS1TOS2Transient + PPS1TOS2Permanet)
//             {
//                 possibleJumpState = "S22";
//             }
//             else 
//             {
//                 possibleJumpState = "S3";
//             }
// //             cout<<"The possible jump state is "<<possibleJumpState<<endl;
//             if (possibleJumpState == "S1")// if it stays at old state
//             {
//                 cellStateMap[cellID].age = cellStateMap[cellID].age +increaseTime;// if not jump then stay at old state, age increase
//             }
//             if(possibleJumpState == "S21")// if it is transient arrest
//             {
// 
//                 double eta = G4UniformRand();
//                 double p ;
//                 if (cellExternalPerturbationEnergyMap[cellID]>0)
//                 {
//                     p = InstantaneousStatincreaseTimeeJumpProb(PS1ToS2Transient,increaseTime,ObservationTime);
//                 }
//                 else 
//                 {
//                     p = DelayedStateJumpProb(PS1ToS2Transient,increaseTime,ObservationTime);
//                 }
//                 if (eta <p)
//                 {                
//                     cellStateMap[cellID].state = "S21";// jumping to S21 state
//                     cellStateMap[cellID].age = 0;
//                     double E20= GaussianSampling2Pi(E21,sigma);// Sample a E20 of S2 state
//                     double pi = 3.1415926535897;
//                     double arrestedDuration = E20/sqrt(2*pi)/sigma*(f1/lambda1+f2/lambda2);
//                     cellStateMap[cellID].duration = arrestedDuration;
//                 }
//                 else
//                 {
//                     cellStateMap[cellID].age = cellStateMap[cellID].age +increaseTime;// if not jump then stay at old state, age increase
//                 }
//             }
//             if(possibleJumpState == "S22") // if it is permanent arrest
//             {
//                 double eta = G4UniformRand();
//                 double p ;
//                 if (cellExternalPerturbationEnergyMap[cellID]>0)
//                 {
//                     p = InstantaneousStateJumpProb(PS1ToS2Permanent,increaseTime,ObservationTime);
//                 }
//                 else 
//                 {
//                     p = DelayedStateJumpProb(PS1ToS2Permanent,increaseTime,ObservationTime);
//                 }
//                 if (eta<p)
//                 {
//                     cellStateMap[cellID].state = "S22";
//                     cellStateMap[cellID].age = 0;
//                     cellStateMap[cellIDincreaseTime].duration = 1E+18;
//                 }
//                 else
//                 {
//                     cellStateMap[cellID].age = cellStateMap[cellID].age +increaseTime;
//                 }
//             }
//             if(possibleJumpState == "S3") // if it is death state
//             {
//                 double eta = G4UniformRand();
//                 double p ;
//                 if (cellExternalPerturbationEnergyMap[cellID]>0)
//                 {
//                     p = InstantaneousStateJumpProb(PS1ToS3,increaseTime,ObservationTime);
//                 }
//                 else 
//                 {
//                     p = DelayedStateJumpProb(PS1ToS3,increaseTime,ObservationTime);
//                 }
// 
//                 if(eta<p)
//                 {
//                     cellStateMap[cellID].state = "S3";
//                     cellStateMap[cellID].age = 0;
//                     cellStateMap[cellID].duration = 1E+18;
//                 }
//                 else
//                 {
//                     cellStateMap[cellID].age = cellStateMap[cellID].age +increaseTime;// if not jump then stay at old state, age increase               
//                 }
//             }
// 
//             return;// after finish processing the possibility of state S1, return
//         }
//         
//         if(cellStateMap[cellID].state=="S21")// when the original state is S21
//         {
//             double ES21ToS3 = -fabs(E21+cellExternalPerturbationEnergyMap.at(cellID)-E3)/2/sigma;
//             double PS21ToS3=2*GaussianCDF(ES21ToS3,0,1);
//             if (E21+cellExternalPerturbationEnergyMap.at(cellID)-E3>=0)// boundary condition
//             {
//                 PS21ToS3 = 1;
//             }
//             elseDSB number big than zero?
//             {
//                 PS21ToS3 = 2*GaussianCDF(ES21ToS3, 0,1);
//             }
//             
//             double eta = G4UniformRand();
//             double p ;
//             if (cellExternalPerturbationEnergyMap[cellID]>0)
//             {
//                 p = InstantaneousStateJumpProb(PS21ToS3,increaseTime,ObservationTime);
//             }
//             else 
//             {
//                 p = DelayedStateJumpProb(PS21ToS3,increaseTime,ObservationTime);
//             }
//             if(eta< p)
//             {increaseTime
//                 cellStateMap[cellID].state = "S3";
//                 cellStateMap[cellID].age = 0;
//                 cellStateMap[cellID].duration = 1E+18;
//             }
//             else
//             {
//                 cellStateMap[cellID].age = cellStateMap[cellID].age +increaseTime;// if not jump then stay at old state, age increase     
//                 if(cellStateMap[cellID].age>=cellStateMap[cellID].duration)
//                 {
//                     cellStateMap[cellID].state = "S1"; // after arrest, cell jumps back to health state S1 from transient arrest state S21
//                     cellStateMap[cellID].age = 0;
//                     cellStateMap[cellID].duration = 1E+18;
//                 }
//             }
//             return;// after finish processing probability of state S21, return
//         }
//         
//         if(cellStateMap[cellID].state=="S22") // when the original state is S22
//         {
//             cellStateMap[cellID].age = cellStateMap[cellID].age+increaseTime; // when it is in permanent arrest, it is senesence state
//             return;
//             
//         }
//         if(cellStateMap[cellID].state=="S3")
//         {
//             cellStateMap[cellID].age = cellStateMap[cellID].age+increaseTime;
//             return;
//         }
        
    }
    if(cellPhaseMap[cellID].phase=="S")// if it is S phase cell
    {
        double E1_S=0;
        double E3_S = sqrt(E3*E3-8*sigma*sigma*log(fS));// the E3 for S phase cell 
        double E21_S = E3_S-sqrt((E3-E2)*(E3-E2)-8*sigma*sigma*log(fS)); // mean state energy of first degenerate state of S2
//         cout<<"Mean state energies are "<<" E3="<<E3_S<<" E2="<<E21_S<<endl;
        CellStateTransition(cellID,E1_S,E21_S,E3_S,increaseTime,ObservationTime,transitionType);
    
    }
    if(cellPhaseMap[cellID].phase=="G2")
    {
        double E1_G2=0;
        double E3_G2 = sqrt(E3*E3-8*sigma*sigma*log(fG2));// the E3 for S phase cell  
        double E21_G2 = E3_G2-sqrt((E3-E2)*(E3-E2)-8*sigma*sigma*log(fG2)); // mean state energy of first degenerate state of S2, it is different compared to G1
//         cout<<"Mean state energies are "<<" E3="<<E3_G2<<" E2="<<E21_G2<<endl;
        CellStateTransition(cellID,E1_G2,E21_G2,E3_G2,increaseTime,ObservationTime,transitionType);
        
    }
    if(cellPhaseMap[cellID].phase=="M")
    {
        double E1_M=0;
        double E3_M = sqrt(E3*E3-8*sigma*sigma*log(fM));// the E3 for S phase cell  
        double E21_M = E3_M-sqrt((E3-E2)*(E3-E2)-8*sigma*sigma*log(fM)); // mean state energy of first degenerate state of S2, it is different compared to G1
//         cout<<"Mean state energies are "<<" E3="<<E3_M<<" E2="<<E21_M<<endl;
        CellStateTransition(cellID,E1_M,E21_M,E3_M,increaseTime,ObservationTime,transitionType);
        
    }
    if(cellPhaseMap[cellID].phase == "G0")// when the cell phase is G0 , here take its radiation sensitivity as G1
    {
        double E1_G0=0;
        double E21_G0=E2;
        double E3_G0=E3;
        CellStateTransition(cellID,E1_G0,E21_G0,E3_G0,increaseTime,ObservationTime,transitionType);

    }
   
}
void CellStateModel::CellStateTransition(int cellID, double E1,double E21, double E3,double increaseTime, double ObservationTime,bool transitionType)
{
    if(cellStateMap[cellID].state=="S1")// when the original state is S1
    {

        double ES1ToS2Transient = -fabs(E1+cellExternalPerturbationEnergyMap.at(cellID)-E21)/2/sigma;
        double PS1ToS2Transient = 2*GaussianCDF(ES1ToS2Transient,0,1);
        double ES1ToS3 = -fabs(E1+cellExternalPerturbationEnergyMap.at(cellID)-E3)/2/sigma;
        double PS1ToS3 = 2*GaussianCDF(ES1ToS3,0,1); 
        
        if (E1+cellExternalPerturbationEnergyMap.at(cellID)-E3>=0)// boundary condition
        {
            PS1ToS3 = 1;
        }
        else
        {
            PS1ToS3 = 2*GaussianCDF(ES1ToS3,0,1); 
        }
        double PPS1TOS2Transient = PS1ToS2Transient/(PS1ToS2Transient + PS1ToS3);
        double PPS1TOS3 = PS1ToS3/(PS1ToS2Transient+ PS1ToS3);
        
        std::string possibleJumpState;
        double rand = G4UniformRand();
        
        if (rand <= PPS1TOS2Transient)
        {
            possibleJumpState = "S21";
        }
        else 
        {
            possibleJumpState = "S3";
        }
        
        if(possibleJumpState=="S21")// if it is transient arrest
        {
//             cout<<" try to jump transient state"<<endl;
            double eta = G4UniformRand();
            double p ;
            if (transitionType==true)// if transitionType is 1, then it is Instantaneous jump model, otherwise delayed model
            {
                p = InstantaneousStateJumpProb(PS1ToS2Transient,increaseTime,ObservationTime);
            }
            else 
            {
                p = DelayedStateJumpProb(PS1ToS2Transient,increaseTime,ObservationTime);
            }
            if(eta<p)
            {
                cellStateMap[cellID].state = "S21";// jumping to S21 state
                cellStateMap[cellID].age = 0;
                double E20= GaussianSampling2Pi(E21,sigma);// Sample a E20 of S2 state
                double pi = 3.1415926535897;
                double arrestedDuration = E20/sqrt(2*pi)/sigma*(f1/lambda1+f2/lambda2);
                cellStateMap[cellID].duration = arrestedDuration;
            }
            else
            {
                cellStateMap[cellID].age = cellStateMap[cellID].age +increaseTime;// if not jump then stay at old state, age increase
            }
        }

        if( possibleJumpState == "S3") // if it is death state
        {
//             cout<<"try to jump to death state"<<endl;
            double eta = G4UniformRand();
            double p ;
            if (transitionType==1)
            {
                p = InstantaneousStateJumpProb(PS1ToS3,increaseTime,ObservationTime);
            }
            else 
            {
                p = DelayedStateJumpProb(PS1ToS3,increaseTime,ObservationTime);
            }
            if(eta<p)
            {
                cellStateMap[cellID].state = "S3";
                cellStateMap[cellID].age = 0;
                cellStateMap[cellID].duration = 1E+18;
            }
            else
            {
                cellStateMap[cellID].age = cellStateMap[cellID].age +increaseTime;// if not jump then stay at old state, age increase               
            }
        }
        return;
    }
    
    if(cellStateMap[cellID].state=="S21")// when the original state is S21
    {
        double ES21ToS3 = -fabs(E21+cellExternalPerturbationEnergyMap.at(cellID)-E3)/2/sigma;
        double PS21ToS3 = 2*GaussianCDF(ES21ToS3, 0,1);
        if (E21+cellExternalPerturbationEnergyMap.at(cellID)-E3>=0)// boundary condition
        {
            PS21ToS3 = 1;
        }
        else
        {
            PS21ToS3 = 2*GaussianCDF(ES21ToS3, 0,1);
        }
        double eta = G4UniformRand();
        double p ;
        if (transitionType==1)
        {
            p = InstantaneousStateJumpProb(PS21ToS3,increaseTime,ObservationTime);
        }
        else 
        {
            p = DelayedStateJumpProb(PS21ToS3,increaseTime,ObservationTime);
        }
        if(eta<p)
        {
            cellStateMap[cellID].state = "S3";
            cellStateMap[cellID].age = 0;
            cellStateMap[cellID].duration = 1E+18;
        }
        else
        {
            cellStateMap[cellID].age = cellStateMap[cellID].age +increaseTime;// if not jump then stay at old state, age increase     
            if(cellStateMap[cellID].age>=cellStateMap[cellID].duration)
            {
                cellStateMap[cellID].state = "S1"; // after arrest, cell jumps back to health state S1 from transient arrest state S21
                cellStateMap[cellID].age = 0;
                cellStateMap[cellID].duration = 1E+18;
            }
        }
        return;
    }
    if(cellStateMap[cellID].state=="S3")
    {
        cellStateMap[cellID].age = cellStateMap[cellID].age+increaseTime;
        return;
    }
        

}

double CellStateModel::InstantaneousStateJumpProb(double p_sp, double increaseTime, double ObservationTime)
{
    //     double p = (1-pow((1-p_sp),increaseTime/ObservationTime));// probability in one time step
    double p = p_sp;
    return p;

}
double CellStateModel::DelayedStateJumpProb(double p_sp, double increaseTime, double ObservationTime)
{
        double p = (1-pow((1-p_sp),increaseTime/ObservationTime));// probability in one time step
//     double p = p_sp;
    return p;

}


void CellStateModel::CellStateModelParameterSetup(double alphaVal, double betaVal, double E1Val,\
double E2Val, double E3Val,double sigmaVal, double fG1Val,double fSVal, double fG2Val, double fMVal,\
double f1Val,double lambda1Val, double f2Val, double lambda2Val)
{
    alpha = alphaVal; // set up the alpha,for deltaE = alpha*DSB_num
    beta = betaVal;// set up for the beta, for deltaE = beta* integralConcentration
    E1 = E1Val;// mean state energy of S1 state
    E2 = E2Val;// mean state energy of S2 state
    E3 = E3Val;// mean state energy of S3 state
    sigma = sigmaVal;
    fG1 = fG1Val;
    fS = fSVal;
    fG2 = fG2Val;
    fM = fMVal;
    f1=f1Val;
    f2= f2Val;
    lambda1= lambda1Val;
    lambda2 = lambda2Val;

}



map< int, string > CellStateModel::GetCellPhase()
{
    std::map<int, std::string> phaseMap;
    for (std::map<int,PhaseInfo>::iterator mitr_cell=cellPhaseMap.begin();mitr_cell!=cellPhaseMap.end();mitr_cell++)
    {
        phaseMap[mitr_cell->first] = mitr_cell->second.phase;
    }
    return phaseMap;

}

map< int, double > CellStateModel::GetCellPhaseDuration()
{
    std::map<int, double> durationMap;
    for (std::map<int,PhaseInfo>::iterator mitr_cell=cellPhaseMap.begin();mitr_cell!=cellPhaseMap.end();mitr_cell++)
    {
        durationMap[mitr_cell->first] = mitr_cell->second.duration;
    }
    return durationMap;

}
map< int, double > CellStateModel::GetCellAge()
{
    std::map<int, double> ageMap;
    for (std::map<int,PhaseInfo>::iterator mitr_cell=cellPhaseMap.begin();mitr_cell!=cellPhaseMap.end();mitr_cell++)
    {
        ageMap[mitr_cell->first] = mitr_cell->second.age;
    }
    return ageMap;

}

map< int, string > CellStateModel::GetCellState()
{
    std::map<int, std::string> stateMap;
    for(std::map<int, StateInfo>::iterator mitr_cell=cellStateMap.begin();mitr_cell!=cellStateMap.end();mitr_cell++)
    {
        stateMap[mitr_cell->first] = mitr_cell->second.state;
    }

    return stateMap;
}
map< int, double > CellStateModel::GetCellStateDuration()
{
    std::map<int, double> durationMap;
    for(std::map<int, StateInfo>::iterator mitr_cell=cellStateMap.begin();mitr_cell!=cellStateMap.end();mitr_cell++)
    {
        durationMap[mitr_cell->first] = mitr_cell->second.duration;
    }    
    return durationMap;
}

map< int, double > CellStateModel::GetCellStateAge()
{
    std::map<int, double> ageMap;
    for(std::map<int, StateInfo>::iterator mitr_cell=cellStateMap.begin();mitr_cell!=cellStateMap.end();mitr_cell++)
    {
        ageMap[mitr_cell->first] = mitr_cell->second.age;
    }     
    return ageMap;
}

map< int, double > CellStateModel::GetCellPositionX()
{
    std::map<int, double> XPositionMap;
    for (std::map<int, PositionInfo>::iterator mitr_cell = cellPositionMap.begin();mitr_cell!=cellPositionMap.end();mitr_cell++)
    {
//         XPositionMap[mitr_cell->first] = d/2+mitr_cell->second.i*d;
        XPositionMap[mitr_cell->first] = mitr_cell->second.i*d;// when take the grid points as center of cell
    }
    return XPositionMap;

}
map< int, double > CellStateModel::GetCellPositionY()
{
    std::map<int, double> YPositionMap;
    for (std::map<int, PositionInfo>::iterator mitr_cell = cellPositionMap.begin();mitr_cell!=cellPositionMap.end();mitr_cell++)
    {
//         YPositionMap[mitr_cell->first] = d/2+mitr_cell->second.j*d;
        YPositionMap[mitr_cell->first] = mitr_cell->second.j*d;//when take the grid points as center of cell
    }
    return YPositionMap;
}

map< int, double > CellStateModel::GetCellPositionZ()
{
    std::map<int, double> ZPositionMap;
    for (std::map<int, PositionInfo>::iterator mitr_cell = cellPositionMap.begin();mitr_cell!=cellPositionMap.end();mitr_cell++)
    {
//         ZPositionMap[mitr_cell->first] = d/2 + mitr_cell->second.k*d;
        ZPositionMap[mitr_cell->first] =  mitr_cell->second.k*d;//when take the grid points as center of cell
    }
    return ZPositionMap;
}

map< int, int > CellStateModel::GetCellAncestryID()
{
    std::map<int, int> cellAncestryIDMap;
    for (std::map<int, PhaseInfo>::iterator mitr_cell = cellPhaseMap.begin();mitr_cell!=cellPhaseMap.end();mitr_cell++)
    {
        cellAncestryIDMap[mitr_cell->first] = mitr_cell->second.ancestryID;
    }
    return cellAncestryIDMap;

}

void CellStateModel::SetUpContactInhibition(bool considerOrNot)
{
    considerContactInhibition = considerOrNot;

}
