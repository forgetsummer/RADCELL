#include "ReactionDiffusionSimulation.hh"
#include "DiffusionReactionSolver.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>

void DiffusionReactionSolver::SetVerbose(int verboseLevel)
{
    verbose = verboseLevel;

}

void DiffusionReactionSolver::CellStateUpdateInitialize()
{
    cX.clear();
    cY.clear();
    cZ.clear();
    cellState.clear();// clear old cell state information 
    secreteRate.clear(); // clear old cell signalling rate information 

}

void DiffusionReactionSolver::CellStateUpdate(int cellID, double x, double y, double z, int state,double fmu)// need to loop all the cells
{

    cX[cellID]=x;
    cY[cellID]=y;
    cZ[cellID]=z;
    cellState[cellID]=state; // update cell state information 
    secreteRate[cellID] = fmu; // update the cell siganlling rate
    
}


void DiffusionReactionSolver::DiffusionReactionInitialization(double dX, double dY, double dZ, double gridSize, double DiffC, double reactionRate, int receptorNum, double delta_T)
{
    dimX=dX;
    dimY=dY;
    dimZ=dZ;
    d=gridSize;
    D=DiffC;
    r=reactionRate;
    Rt=receptorNum;
    dt=delta_T;

    N_X = ceil(dX/d);
    N_Y = ceil(dY/d);
    N_Z = ceil(dZ/d);

    if (N_Z==0)
    {
        N_Z =1;
    }
    if (N_Y==0)
    {
        N_Y =1;
    }
    if (N_X==0)
    {
        N_X =1;
    }

    for (int i=0;i<N_X;i++ )
    {
        initialC.push_back(std::vector<std::vector<double> >());
        for (int j=0;j<N_Y;j++)
        {
            initialC[i].push_back(std::vector<double>());
            for (int k=0;k<N_Z;k++)
            {
                initialC[i][j].push_back(0); // initial concentration
            }
        }
    }
    for (int i=0;i<N_X;i++ )
    {
        integralConcentration.push_back(std::vector<std::vector<double> >());
        for (int j=0;j<N_Y;j++)
        {
            integralConcentration[i].push_back(std::vector<double>());
            for (int k=0;k<N_Z;k++)
            {
                integralConcentration[i][j].push_back(0); // initial concentration
            }
        }
    }

}

void DiffusionReactionSolver::DiffusionReactionCalculation(double deltaT_diffusion, int frequency)
{
    double totalTime;
    totalTime = deltaT_diffusion*frequency;
    ReactionDiffusionSimulation react;
    react.SetVerbose(verbose);
    react.MeshGeometry(dimX,dimY,dimZ,d);// mesh the geometry for diffusion reaction calculation
    react.CellStateInitialization(cX,cY,cZ,cellState,secreteRate); // initialize the cell state
    react.DiffusionReactionCalculation(initialC,D,r,Rt,dt,totalTime); // call the diffusion solver 
    concentration = react.GetConcentration();
    std::vector<std::vector<std::vector<double> > > integralC;
    integralC = react.GetIntegralConcentration();
    for (int i=0;i<N_X;i++ )
    {
        for (int j=0;j<N_Y;j++)
        {
            for (int k=0;k<N_Z;k++)
            {
                integralConcentration[i][j][k] = integralConcentration[i][j][k]+ integralC[i][j][k]; // initial concentration
            }
        }
    }

    initialC = concentration; // update initialC using the calculated concentration

}


double DiffusionReactionSolver::GetConcentration(double x, double y, double z)
{
    double c_i = floor(x/d)-1;
    double c_j = floor(y/d)-1;
    double c_k = floor(z/d)-1;            
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
//     cout<<"x="<<x<<" , "<<"y="<<y<<" , "<<"z="<<z<<endl;
//     cout<<endl<<"i="<<c_i<<" , "<<"j="<<c_j<<" , "<<"k="<<c_k<<endl;
//     cout<<endl<<"N_X="<<N_X<<" , "<<"N_Y="<<N_Y<<" , "<<"N_Z="<<N_Z<<endl;
//     exit(0);
    
//     cout<<"THE INDEX OF POINT IS "<<c_i<<","<<c_j<<","<<c_k<<endl;
    if (c_i<0||c_i>(N_X-1)||c_j<0||c_j>(N_Y-1)||c_k<0||c_k>(N_Z-1))// check the point whether it is within the system dimension
    {
        cout<<endl<<"i="<<c_i<<" , "<<"j="<<c_j<<" , "<<"k="<<c_k<<endl;
        cout<<endl<<"N_X="<<N_X<<" , "<<"N_Y="<<N_Y<<" , "<<"N_Z="<<N_Z<<endl;
        cout<<"ERROR:THE POINT IS OUT OFF THE SYSTEM DIMENSION !"<<endl;
        exit(1);
    }
    return concentration[c_i][c_j][c_k];

}

double DiffusionReactionSolver::GetIntegralConcentration(double x, double y, double z)
{ 
    double c_i = floor(x/d)-1;
    double c_j = floor(y/d)-1;
    double c_k = floor(z/d)-1;            
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
    
//     cout<<"THE INDEX OF POINT IS "<<c_i<<","<<c_j<<","<<c_k<<endl;
    if (c_i<0||c_i>(N_X-1)||c_j<0||c_j>(N_Y-1)||c_k<0||c_k>(N_Z-1))// check the point whether it is within the system dimension
    {
        cout<<"ERROR:THE POINT IS OUT OFF THE SYSTEM DIMENSION !"<<endl;
        exit(1);
    }
    return integralConcentration[c_i][c_j][c_k];
}


void DiffusionReactionSolver::WriteConcentrationToFile(string path,int fileID)
{
    cout<<endl<<"WRITING CALCULATED CONCENTRATION TO FILE:"<<endl;
    stringstream ss;
    ss<<fileID;
    string filename;
    string filename1;
    filename = path+ "concentration-"+ss.str()+".csv";
    filename1= path+ "integralConcentration-"+ss.str()+".csv";
    ofstream file;
    ofstream file1;
    file.open(filename.c_str());
//     file1.open(filename1.c_str());
    for (int i=0;i<concentration.size();i++)
    {
        for (int j=0;j<concentration[i].size();j++)
        {
            for (int k=0;k<concentration[i][j].size();k++)
            {
//                 cout<<i<<" "<<j<<" "<<k<<" : "<<Concentration[nn][i][j][k]<<endl;
                file<<i<<","<<j<<","<<k<<","<<concentration[i][j][k]<<endl;
//                 file1<<i<<","<<j<<","<<k<<","<<integralConcentration[i][j][k]<<endl;
            }
        }
        
    } 

}


