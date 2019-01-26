#include "ReactionDiffusionSimulation.hh"

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

void ReactionDiffusionSimulation::SetVerbose(int verboseLevel)
{
    verbose = verboseLevel;
}

void ReactionDiffusionSimulation::MeshGeometry(double dX, double dY, double dZ, double gridSize)
{
    xDim = dX;
    yDim = dY;
    zDim = dZ;
    d = gridSize;

    int N_X = ceil(dX/d);// Total grid cell number in x dimention
    int N_Y = ceil(dY/d);// Total grid cell number in y dimension
    int N_Z = ceil(dZ/d);// Total grid cell number in z dimension
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
    gridNumX = N_X;
    gridNumY = N_Y;
    gridNumZ = N_Z;


    for (int i=0;i<N_X;i++)// loop grids in z dimension
    {        
        gridIndex.push_back(std::vector<std::vector<IJKVec> >());
        for (int j=0;j<N_Y;j++) //loop grids in y dimension
        {
            gridIndex[i].push_back(std::vector<IJKVec>());
            {
                for (int k=0;k<N_Z;k++) // loop grids in z dimension
                {
                    IJKVec index;
                    index.i = i;
                    index.j = j;
                    index.k = k;
                    gridIndex[i][j].push_back(index);  
                }
            }
        }
    }

/*
    for (int i = 0;i<N_X;i++)
    {
        for (int j =0;j<N_Y;j++)
        {
            for (int k=0; k<N_Z;k++)
            {
                cout<<"The grid information is as "<<gridIndex[i][j][k].i<<","<<gridIndex[i][j][k].j<<","<<gridIndex[i][j][k].k<<endl;
            }
        }
    }*/
    
    if (verbose==1)// if verbose is 1, then output these resluts, otherwise, do not output
    {
        cout<<"Dim number is "<<N_X <<","<<N_Y<<","<<N_Z<<endl;  
        cout<<"The size of gridIndex in X dimension is "<<gridIndex.size()<<endl;
        cout<<"The size of gridIndex in Y dimension is "<<gridIndex[0].size()<<endl;
        cout<<"The size of gridIndex in Z dimension is "<<gridIndex[0][0].size()<<endl;
    }

}

void ReactionDiffusionSimulation::CellStateInitialization(map< int, double > cX, map< int, double > cY, map< int, double > cZ, map< int, int > cState, map< int, double > secRMap)
{
    cellState = cState; // set up the cell state map
    secreteRateMap = secRMap; // set up the cell singalling rate map
    for (std::map<int,double>::iterator mitr_cell = cX.begin();mitr_cell!=cX.end();mitr_cell++) // loop all the biological cells
    {
        double c_i = floor(mitr_cell->second/d);
        double c_j = floor(cY[mitr_cell->first]/d);
        double c_k = floor(cZ[mitr_cell->first]/d);
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
        IJKVec c_index;
        c_index.i = c_i;
        c_index.j = c_j;
        c_index.k = c_k;
//         cout<<"The cell index is "<<c_index.i<<" "<<c_index.j<<" "<<c_index.k<<endl;
        cellIndex[mitr_cell->first] = c_index; // get the biological cell's index in cellular automaton
        
    }

}


pair< map< int, vector< vector< vector< double > > > >, vector< vector< vector< double > > > > ReactionDiffusionSimulation::DiffusionReactionCalculation3D(vector< vector< vector< double > > >& initialC, double DiffC, double reactionRate, int Rt, double delta_T, double T)
{
    cout<<endl<<"RUNNING 3D DIFFUSION-REACTION CALCULATION"<<endl;
    double D = DiffC; // diffusion coefficient, unit is mm^2/s
    std::map<int, double> mu = secreteRateMap; // secretion rate, unit is moleclues/second
    double r = reactionRate; // reaction rate coefficient, unit in 1/moleclues/second
    
    double deltaT_max;
    deltaT_max = pow(d,2)/D/6; // The maximum deltaT  constraint: 6D*deltaT <= d^2 , the units used here is micrometer
    if (delta_T>deltaT_max)
    {
        deltaT = deltaT_max;
    }
    else
    {
        deltaT = delta_T;
    }
    
    int N =ceil(T/deltaT); // calcuate the total number of time steps
    if (verbose==1)
    {
        cout<<"Total time length is: "<<T <<" Time step deltaT is: "<<deltaT <<endl; 
        cout<<"Total time step number N is: "<<N<<endl; 
    }

    
    std::map<int, std::vector<std::vector<std::vector<double> > > > Concentration; // A map containing the chemical Concentration at all grid cell
    Concentration[0] = initialC; // initial Concentration
    std::vector<std::vector<std::vector<double> > > integralC;
    
    std::vector<std::vector<std::vector<int> > > initialActiveReceptor;
    std::vector<std::vector<std::vector<double> > >  secretionPerStep;

    for (int i=0;i<gridIndex.size();i++ )
    {
        secretionPerStep.push_back(std::vector<std::vector<double> >());
        for (int j=0;j<gridIndex[i].size();j++)
        {
            secretionPerStep[i].push_back(std::vector<double>());
            for (int k=0;k<gridIndex[i][j].size();k++)
            {
                secretionPerStep[i][j].push_back(0); // initial secretionPerStep
            }
        }
    }
    
    for (std::map<int, IJKVec>::iterator mitr_cell=cellIndex.begin();mitr_cell!=cellIndex.end();mitr_cell++)
    {
        if(cellState[mitr_cell->first]==2 ||cellState[mitr_cell->first]==3) // if cell state in arrest state or dead state it emits signal
        {
            int c_i = mitr_cell->second.i;
            int c_j = mitr_cell->second.j;
            int c_k = mitr_cell->second.k;
//             cout<<"The singalling cell index is "<<c_i<<" "<<c_j<<" "<<c_k<<endl;
            secretionPerStep[c_i][c_j][c_k] = mu[mitr_cell->first]*delta_T; // update secretionPerStep by looping all the biological cells per deltaT
        }
    }
         
    for (int i=0;i<gridIndex.size();i++ )
    {
        initialActiveReceptor.push_back(std::vector<std::vector<int> >());
        for (int j=0;j<gridIndex[i].size();j++)
        {
            initialActiveReceptor[i].push_back(std::vector<int>());
            for (int k=0;k<gridIndex[i][j].size();k++)
            {
                initialActiveReceptor[i][j].push_back(0);
            }
        }
    }

    for (std::map<int, IJKVec>::iterator mitr_cell=cellIndex.begin();mitr_cell!=cellIndex.end();mitr_cell++)
    {
        if (cellState[mitr_cell->first] == 1 || cellState[mitr_cell->first] == 2) // only healthy state and arrest state cell reacts with signals
        {
            int c_i = mitr_cell->second.i;
            int c_j = mitr_cell->second.j;
            int c_k = mitr_cell->second.k; 
            initialActiveReceptor[c_i][c_j][c_k] = Rt;

        }
    }
       
    std::vector<std::vector<std::vector<double> > >  CI;// initial Concentration at one time step
    std::vector<std::vector<std::vector<double> > >  CF; // final Concentration at one time step
    std::vector<std::vector<std::vector<int> > > RI;
    std::vector<std::vector<std::vector<int> > > RF;
    
    for (int i=0;i<gridIndex.size();i++ )
    {
        RF.push_back(std::vector<std::vector<int> >());
        for (int j=0;j<gridIndex[i].size();j++)
        {
            RF[i].push_back(std::vector<int>());
            for (int k=0;k<gridIndex[i][j].size();k++)
            {
                RF[i][j].push_back(0);
            }
        }
    }
    
    for (int i=0;i<gridIndex.size();i++ )
    {
        CF.push_back(std::vector<std::vector<double> >());
        for (int j=0;j<gridIndex[i].size();j++)
        {
            CF[i].push_back(std::vector<double>());
            for (int k=0;k<gridIndex[i][j].size();k++)
            {
                CF[i][j].push_back(0);
            }
        }

    }
    for (int i=0;i<gridIndex.size();i++)
    {
        integralC.push_back(std::vector<std::vector<double> >());
        for (int j=0;j<gridIndex[i].size();j++)
        {
            integralC[i].push_back(std::vector<double>());
            for (int k=0;k<gridIndex[i][j].size();k++)
            {
                integralC[i][j].push_back(0);
            }
        }
    }

    RI = initialActiveReceptor;
    CI = initialC; 
    if (verbose ==1)
    {
        cout<<"Before calculation"<<endl;
    }
      
    for (int n=0;n<=N;n++) // loop all the time steps
    {
        if (verbose==1)
        {
            cout<<"Get to time step: "<<n<<endl;     
        }

        for (int i=1;i<gridIndex.size()-1;i++)
        {
            for (int j=1;j<gridIndex[i].size()-1;j++)
            {
                for (int k=1;k<gridIndex[i][j].size()-1;k++)
                {
//                     cout<<"get here for k "<<i<<" "<<j<<" "<<k<<endl;

                    CF[i][j][k] = (1-6*D*deltaT/pow(d,2))*CI[i][j][k]+\
                    D*deltaT/pow(d,2)*(CI[i+1][j][k]+CI[i-1][j][k]+CI[i][j+1][k]+CI[i][j-1][k]+CI[i][j][k+1]+CI[i][j][k-1])\
                    -r*RI[i][j][k]*CI[i][j][k]*deltaT+secretionPerStep[i][j][k]; // update concentration
                                        
                    RF[i][j][k] = RI[i][j][k]-r*RI[i][j][k]*CI[i][j][k]*deltaT; // update free receptor number 
                    integralC[i][j][k] = integralC[i][j][k]+ r*RI[i][j][k]*CI[i][j][k]*deltaT; // update integralC
                } 
            }
        }
                
        for (int i=1;i<gridIndex.size()-1;i++)
        {
            for (int j=1;j<gridIndex[i].size()-1;j++)
            {
                for (int k=1;k<gridIndex[i][j].size()-1;k++)
                {
                    RI[i][j][k] = RF[i][j][k];
                }
            }
        } 
        
        for (int j=0;j<gridIndex[0].size();j++)
        {
            for (int k=0;k<gridIndex[0][j].size();k++)
            {
                RI[0][j][k] = RI[0][j][k]- r*RI[0][j][k]*CI[0][j][k]*deltaT;   //boundary condition for free receptor   
                integralC[0][j][k] = integralC[0][j][k]+r*RI[0][j][k]*CI[0][j][k]*deltaT;
            }
        }
        
        for (int i=0;i<gridIndex.size();i++)
        {
            for (int k=0;k<gridIndex[i][0].size();k++)
            {
                RI[i][0][k] = RI[i][0][k]- r*RI[i][0][k]*CI[i][0][k]*deltaT;  //boundary condition for free receptor  
                integralC[i][0][k] = integralC[i][0][k]+r*RI[i][0][k]*CI[i][0][k]*deltaT; 
            }
        }
        
        for (int i=0;i<gridIndex.size();i++)
        {
            for (int j=0;j<gridIndex[i].size();j++)
            {
                RI[i][j][0] = RI[i][j][0]- r*RI[i][j][0]*CI[i][j][0]*deltaT;    //boundary condition for free receptor  
                integralC[i][j][0] = integralC[i][j][0]+ r*RI[i][j][0]*CI[i][j][0]*deltaT; 
            }
        }
        

        int Im= gridIndex.size()-1;     
        for (int j=0;j<gridIndex[Im].size();j++)
        {
            for (int k=0;k<gridIndex[Im][j].size();k++)
            {
                RI[Im][j][k] = RI[Im][j][k]- r*RI[Im][j][k]*CI[Im][j][k]*deltaT;  //boundary condition for free receptor 
                integralC[Im][j][k] =  integralC[Im][j][k]+r*RI[Im][j][k]*CI[Im][j][k]*deltaT; 
            }
        }
        
        for (int i=0;i<gridIndex.size();i++)
        {
            int Jm = gridIndex[i].size()-1;
            for (int k=0;k<gridIndex[i][Jm].size();k++)
            {
                RI[i][Jm][k] = RI[i][Jm][k]- r*RI[i][Jm][k]*CI[i][Jm][k]*deltaT;   //boundary condition for free receptor 
                integralC[i][Jm][k] = integralC[i][Jm][k]+r*RI[i][Jm][k]*CI[i][Jm][k]*deltaT; 
            }
        }
        
        for (int i=0;i<gridIndex.size();i++)
        {
            for (int j=0;j<gridIndex[i].size();j++)
            {
                int Km = gridIndex[i][j].size()-1;
                RI[i][j][Km] = RI[i][j][Km]- r*RI[i][j][Km]*CI[i][j][Km]*deltaT;    //boundary condition for free receptor 
                integralC[i][j][Km] = integralC[i][j][Km]+r*RI[i][j][Km]*CI[i][j][Km]*deltaT; 
            }
        }
        
        for (int i=1;i<gridIndex.size()-1;i++)
        {
            for (int j=1;j<gridIndex[i].size()-1;j++)
            {
                for (int k=1;k<gridIndex[i][j].size()-1;k++)
                {
                    CI[i][j][k] = CF[i][j][k]; // update initial concentration for next time step calculation
                }
            }
        } 

        for (int j=0;j<gridIndex[0].size();j++)
        {
            for (int k=0;k<gridIndex[0][j].size();k++)
            {
                CI[0][j][k] = CI[1][j][k];//boundary condition for concentration
                
            }
        }
        
        for (int i=0;i<gridIndex.size();i++)
        {
            for (int k=0;k<gridIndex[i][0].size();k++)
            {
                CI[i][0][k] = CI[i][1][k];//boundary condition for concentration
            }
        }
        
        for (int i=0;i<gridIndex.size();i++)
        {
            for (int j=0;j<gridIndex[i].size();j++)
            {
                CI[i][j][0] = CI[i][j][1];//boundary condition for concentration
            } 
        }
        
        
        for (int j=0;j<gridIndex[Im].size();j++)
        {
            for (int k=0;k<gridIndex[Im][j].size();k++)
            {
                CI[Im][j][k] = CI[Im-1][j][k];//boundary condition for concentration
            }
        }
        
        for (int i=0;i<gridIndex.size();i++)
        {
            int Jm = gridIndex[i].size()-1;
            for (int k=0;k<gridIndex[i][Jm].size();k++)
            {
                CI[i][Jm][k] = CI[i][Jm-1][k];//boundary condition for concentration
            }
        }
        
        for (int i=0;i<gridIndex.size();i++)
        {
            for (int j=0;j<gridIndex[i].size();j++)
            {
                int Km = gridIndex[i][j].size()-1;
                CI[i][j][Km] = CI[i][j][Km-1];//boundary condition for concentration
            }
        }
        
        Concentration[n] = CI;

    }
    cout<<endl<<"FINISH 3D DIFFUSION REACTION CALCULATION"<<endl;
    std::pair<std::map<int, std::vector<std::vector<std::vector<double> > > >,std::vector<std::vector<std::vector<double> > > > diffusionResults;
    diffusionResults.first = Concentration;
    diffusionResults.second = integralC;
    return diffusionResults;

}


pair< map< int, vector< vector< double > > >, vector< vector< double > > > ReactionDiffusionSimulation::DiffusionReactionCalculation2D(vector< vector< double > >& initialC, double DiffC, double reactionRate, int Rt, double delta_T, double T)
{
    cout<<endl<<"RUNNING 2D DIFFUSION-REACTION CALCULATION"<<endl;
    double D = DiffC; // diffusion coefficient, unit is mm^2/s
    std::map<int, double> mu = secreteRateMap;// secretion rate, unit is moleclues/second
    double r = reactionRate; // reaction rate coefficient, unit in 1/moleclues/second
    
    double deltaT_max;
    deltaT_max = pow(d,2)/D/4; // The maximum deltaT  constraint: 4D*deltaT <= d^2 , the units used here is micrometer
    if (delta_T>deltaT_max)
    {
        deltaT = deltaT_max;
    }
    else
    {
        deltaT = delta_T;
    }
    
    int N =ceil(T/deltaT); // calcuate the total number of time steps
    if(verbose==1)
    {
        cout<<"Total time length is: "<<T <<" Time step deltaT is: "<<deltaT <<endl; 
        cout<<"Total time step number N is: "<<N<<endl; 
    }

    std::map<int, std::vector<std::vector<double> > > Concentration2D;
    Concentration2D[0] = initialC; // initial Concentration
    
    std::vector<std::vector<int> >initialActiveReceptor;
    std::vector<std::vector<double> > secretionPerStep;
    std::vector<std::vector<double> > integralC;
    
    for (int i=0;i<gridIndex2D.size();i++)
    {
        secretionPerStep.push_back(std::vector<double>());
        for (int j=0;j<gridIndex2D[i].size();j++)
        {
            secretionPerStep[i].push_back(0);// initial secretionPerStep
        }
    }
    
    for (std::map<int, IJVec>::iterator mitr_cell=cellIndex2D.begin();mitr_cell!=cellIndex2D.end();mitr_cell++)
    {
        if(cellState[mitr_cell->first]==2 ||cellState[mitr_cell->first]==3) // if cell state in arrest state or dead state it emits signal
        {
            int c_i = mitr_cell->second.i;
            int c_j = mitr_cell->second.j;
//             cout<<"The singalling cell index is: "<<c_i<<" "<<c_j<<endl;
            secretionPerStep[c_i][c_j] = mu[mitr_cell->first]*delta_T; // update secretionPerStep by looping all the biological cells per deltaT
        }
        
    }
    
    for (int i=0;i<gridIndex2D.size();i++)
    {
        initialActiveReceptor.push_back(std::vector<int>());
        for (int j=0;j<gridIndex2D[i].size();j++)
        {
            initialActiveReceptor[i].push_back(0);
        }
    }
    
    for (std::map<int, IJVec>::iterator mitr_cell=cellIndex2D.begin();mitr_cell!=cellIndex2D.end();mitr_cell++)
    {
         if (cellState[mitr_cell->first] == 1 || cellState[mitr_cell->first] == 2) // only healthy state and arrest state cell reacts with signals
         {
            int c_i = mitr_cell->second.i;
            int c_j = mitr_cell->second.j;
            initialActiveReceptor[c_i][c_j] = Rt;     
         }   
    }
    
    std::vector<std::vector<double> >CI;// initial Concentration at one time step
    std::vector<std::vector<double> >CF;// final Concentration at one time step
    
    std::vector<std::vector<int> > RI;
    std::vector<std::vector<int> > RF;
    
    for (int i=0;i<gridIndex2D.size();i++)
    {
        RF.push_back(std::vector<int>());
        for (int j=0;j<gridIndex2D[i].size();j++)
        {
            RF[i].push_back(0);
        }
    }
    
    for (int i=0;i<gridIndex2D.size();i++)
    {
        CF.push_back(std::vector<double>());
        for (int j=0;j<gridIndex2D[i].size();j++)
        {
            CF[i].push_back(0);
        }
    }
    
    for (int i=0;i<gridIndex2D.size();i++)
    {
        integralC.push_back(std::vector<double>());
        for (int j=0;j<gridIndex2D[i].size();j++)
        {
            integralC[i].push_back(0);
        }
    }

    RI = initialActiveReceptor;
    CI = initialC; 
    
    if(verbose==1)
    {
        cout<<"Begin calculation:"<<endl;
        cout<<"size of initialC is "<<initialC.size() <<" "<<initialC[0].size()<<endl; 
    }
    for (int n=1;n<=N;n++)// loop all the time steps
    {
//         cout<<"Get to time step: "<<n<<endl;
        for (int i=1;i<gridIndex2D.size()-1;i++)
        {
            for (int j=1;j<gridIndex2D[i].size()-1;j++)
            {          
                CF[i][j] = (1-4*D*deltaT/pow(d,2))*CI[i][j]+D*deltaT/pow(d,2)*(CI[i+1][j]+CI[i-1][j]+CI[i][j+1]+CI[i][j-1])\
                -r*RI[i][j]*CI[i][j]*deltaT+secretionPerStep[i][j]; // update concentration
//                 cout<<"seretionPerStep[i][j]"<< i<<","<<j<<" "<<secretionPerStep[i][j]<<endl;
//                 cout<<"CF[i][j]="<<CF[i][j]<<endl;
                
                RF[i][j] = RI[i][j]-r*RI[i][j]*CI[i][j]*deltaT; // update free receptor number, comment this line for model verification using manufactured solution method
                integralC[i][j] = integralC[i][j]+r*RI[i][j]*CI[i][j]*deltaT;
            }
        }
        
        for (int i=1;i<gridIndex2D.size()-1;i++)
        {
            for (int j=1;j<gridIndex2D[i].size()-1;j++)
            {
                RI[i][j] = RF[i][j]; //update total receptor number here
            }
        }
        
        for (int j=0;j<gridIndex2D.size();j++)
        {
            RI[0][j] = RI[0][j]-r*RI[0][j]*CI[0][j]*deltaT;
        }
         
        for (int i=0;i<gridIndex2D.size();i++)
        {
            RI[i][0] = RI[i][0]-r*RI[i][0]*CI[i][0]*deltaT;
        }
        
        int Im = gridIndex2D.size()-1;
        for (int j= 0;j<gridIndex2D[Im].size();j++)
        {
            RI[Im][j] = RI[Im][j]-r*RI[Im][j]*CI[Im][j]*deltaT;
        }
       
        for (int i =0;i<gridIndex2D.size();i++)
        {
            int Jm = gridIndex2D[i].size()-1;
            RI[i][Jm] = RI[i][Jm]-r*RI[i][Jm]*CI[i][Jm]*deltaT;

        }
    
        for (int i=1;i<gridIndex2D.size()-1;i++)
        {
            for (int j=1;j<gridIndex2D[i].size()-1;j++)
            {
                CI[i][j] = CF[i][j]; // update concentration here
            }
        }
        
        for (int j=0;j<gridIndex2D.size();j++)
        {
            CI[0][j] = CI[1][j];
        }
         
        for (int i=0;i<gridIndex2D.size();i++)
        {
            CI[i][0] = CI[i][1];
        }
        
        for (int j= 0;j<gridIndex2D[Im].size();j++)
        {
            CI[Im][j] = CI[Im-1][j];
        }
       
        for (int i =0;i<gridIndex2D.size();i++)
        {
            int Jm = gridIndex2D[i].size()-1;
            CI[i][Jm] = CI[i][Jm-1];
        }
        
        Concentration2D[n] = CI;

    }
    cout<<endl<<"FINISH 2D DIFFUSION REACTION CALCULATION"<<endl;
    std::pair<std::map<int, std::vector<std::vector<double> > >, std::vector<std::vector<double> > > diffusionResults;
    diffusionResults.first = Concentration2D;
    diffusionResults.second = integralC;
    return diffusionResults;

}

void ReactionDiffusionSimulation::DiffusionReactionCalculation(vector< vector< vector< double > > >& initialC, double DiffC, double reactionRate, int Rt, double delta_T, double T)
{
    if (gridIndex.size()>1 && gridIndex[0].size()>1 && gridIndex[0][0].size()>1) // 3D diffusion-reaction 
    {
        std::pair<std::map<int, std::vector<std::vector<std::vector<double> > > >,std::vector<std::vector<std::vector<double> > > > diffusionResults;
        diffusionResults = DiffusionReactionCalculation3D(initialC,DiffC,reactionRate,Rt,delta_T,T);
        concentration = diffusionResults.first;
        integralConcentration = diffusionResults.second;
    }
    
    if (gridIndex.size()>1 && gridIndex[0].size()>1 && gridIndex[0][0].size()==1) // 2D diffusion-reaction in X, Y plane
    { 
        std::vector<std::vector<double> > C_0;
        for (int i=0;i<initialC.size();i++)
        {
            C_0.push_back(std::vector<double>());
            for (int j=0;j<initialC[i].size();j++)
            {
                for (int k=0;k<initialC[i][j].size();k++)
                {
                    C_0[i].push_back(0);
                }
            }
            
        }
        
        for (int i=0;i<initialC.size();i++)
        {
            for (int j=0;j<initialC[i].size();j++)
            {
                for (int k=0;k<initialC[i][j].size();k++)
                {
                    C_0[i][j] = initialC[i][j][k];
                }
            }
        }
      
        for (int i=0;i<gridIndex.size();i++)
        {
            gridIndex2D.push_back(std::vector<IJVec>());
            for (int j=0;j<gridIndex[i].size();j++)
            {
                for (int k=0;k<gridIndex[i][j].size();k++)
                {       
                    IJVec index;  
                    index.i = gridIndex[i][j][k].i; 
                    index.j = gridIndex[i][j][k].j; 
                    gridIndex2D[i].push_back(index);
                }
            }
        }
        
        for (std::map<int, IJKVec>::iterator mitr_cell=cellIndex.begin();mitr_cell!=cellIndex.end();mitr_cell++)
        {
            cellIndex2D[mitr_cell->first].i = mitr_cell->second.i;
            cellIndex2D[mitr_cell->first].j = mitr_cell->second.j;
        }
        std::map<int, std::vector<std::vector<double> > > C2D;
        std::pair<std::map<int, std::vector<std::vector<double> > >, std::vector<std::vector<double> > > diffusionResults;
        diffusionResults = DiffusionReactionCalculation2D(C_0,DiffC,reactionRate,Rt,delta_T,T);
        C2D = diffusionResults.first;

        std::vector<std::vector<std::vector<double> > > C;
  
        for (int i=0;i<gridIndex.size();i++)
        {                  
            C.push_back(std::vector<std::vector<double> >());
            for (int j=0;j<gridIndex[i].size();j++)
            {   
                C[i].push_back(std::vector<double>());
                for (int k=0;k<gridIndex[i][j].size();k++)
                {
                    C[i][j].push_back(0);
                }
            }
        }
             
        for (std::map<int, std::vector<std::vector<double> > >::iterator mitr_time=C2D.begin();\
            mitr_time!=C2D.end();mitr_time++
        )
        {
            concentration[mitr_time->first] = C; // initialize concentration 
        }
        
        for (std::map<int, std::vector<std::vector<double> > >::iterator mitr_time=C2D.begin();\
            mitr_time!=C2D.end();mitr_time++
        )
        {
            for (int i=0;i<mitr_time->second.size();i++)
            {
                for (int j=0;j<mitr_time->second[i].size();j++)
                {
                    concentration[mitr_time->first][i][j][0] = mitr_time->second[i][j];// assigning values for concentration
                }
            }
        }
               
        std::vector<std::vector<double> > integralC = diffusionResults.second;
        for (int i=0;i<gridIndex.size();i++)
        {                  
            integralConcentration.push_back(std::vector<std::vector<double> >());
            for (int j=0;j<gridIndex[i].size();j++)
            {   
                integralConcentration[i].push_back(std::vector<double>());
                for (int k=0;k<gridIndex[i][j].size();k++)
                {
                    integralConcentration[i][j].push_back(0);
                }
            }
        }

        for (int i =0;i<integralC.size();i++)
        {
            for (int j=0;j<integralC[i].size();j++)
            {       
                integralConcentration[i][j][0] = integralC[i][j]; 
            }
        }

    }
    
    if (gridIndex.size()==1 && gridIndex[0].size()>1 && gridIndex[0][0].size()>1)// 2D, Y, Z plane
    {
        std::vector<std::vector<double> > C_0;
        for (int i=0;i<initialC.size();i++)
        {
            for (int j=0;j<initialC[i].size();j++)
            {
                C_0.push_back(std::vector<double>());
                for (int k=0;k<initialC[i][j].size();k++)
                {
                    C_0[j].push_back(0);
                }
            }  
        }
    
        for (int i=0;i<initialC.size();i++)
        {
            for (int j=0;j<initialC[i].size();j++)
            {
                for (int k=0;k<initialC[i][j].size();k++)
                {
                    C_0[j][k] = initialC[i][j][k];
                }
            }
        }
        
        for (int i=0;i<gridIndex.size();i++)
        {
            for (int j=0;j<gridIndex[i].size();j++)
            {
                gridIndex2D.push_back(std::vector<IJVec>());
                for (int k=0;k<gridIndex[i][j].size();k++)
                {
                    IJVec index;
                    index.i = gridIndex[i][j][k].j;
                    index.j = gridIndex[i][j][k].k; 
                    gridIndex2D[j].push_back(index);
                }
            }
        }
        
        for (std::map<int, IJKVec>::iterator mitr_cell=cellIndex.begin();mitr_cell!=cellIndex.end();mitr_cell++)
        {
            cellIndex2D[mitr_cell->first].i = mitr_cell->second.j;
            cellIndex2D[mitr_cell->first].j = mitr_cell->second.k;
        }
        std::map<int, std::vector<std::vector<double> > > C2D;
        std::pair<std::map<int, std::vector<std::vector<double> > >, std::vector<std::vector<double> > > diffusionResults;
        diffusionResults = DiffusionReactionCalculation2D(C_0,DiffC,reactionRate,Rt,delta_T,T);
        C2D = diffusionResults.first;
        
        std::vector<std::vector<std::vector<double> > > C;
  
        for (int i=0;i<gridIndex.size();i++)
        {                  
            C.push_back(std::vector<std::vector<double> >());
            for (int j=0;j<gridIndex[i].size();j++)
            {   
                C[i].push_back(std::vector<double>());
                for (int k=0;k<gridIndex[i][j].size();k++)
                {
                    C[i][j].push_back(0);
                }
            }
        }
             
        for (std::map<int, std::vector<std::vector<double> > >::iterator mitr_time=C2D.begin();\
            mitr_time!=C2D.end();mitr_time++
        )
        {
            concentration[mitr_time->first] = C;
        }
        
        
        for (std::map<int, std::vector<std::vector<double> > >::iterator mitr_time=C2D.begin();\
            mitr_time!=C2D.end();mitr_time++
        )
        {
            for (int i=0;i<mitr_time->second.size();i++)
            {
                for (int j=0;j<mitr_time->second[i].size();j++)
                {
                    concentration[mitr_time->first][0][i][j] = mitr_time->second[i][j];
                }
            }
        }
        
        std::vector<std::vector<double> > integralC = diffusionResults.second;
        for (int i=0;i<gridIndex.size();i++)
        {                  
            integralConcentration.push_back(std::vector<std::vector<double> >());
            for (int j=0;j<gridIndex[i].size();j++)
            {   
                integralConcentration[i].push_back(std::vector<double>());
                for (int k=0;k<gridIndex[i][j].size();k++)
                {
                    integralConcentration[i][j].push_back(0);
                }
            }
        }
        for (int i =0;i<integralC.size();i++)
        {
            for (int j=0;j<integralC[i].size();j++)
            {
                integralConcentration[0][i][j] = integralC[i][j];
            }
        }

    }
    
    if (gridIndex.size()>1 && gridIndex[0].size()==1 && gridIndex[0][0].size()>1)//2D X, Z plane
    {
        std::vector<std::vector<double> > C_0;
        for (int i=0;i<initialC.size();i++)
        {
            C_0.push_back(std::vector<double>());
            for (int j=0;j<initialC[i].size();j++)
            {
                for (int k=0;k<initialC[i][j].size();k++)
                {
                    C_0[i].push_back(0);
                }
            }  
        }

        for (int i=0;i<initialC.size();i++)
        {
            for (int j=0;j<initialC[0].size();j++)
            {
                for (int k=0;k<initialC[i][j].size();k++)
                {
                    C_0[i][k] = initialC[i][j][k];
                }
            }
        }

        for (int i=0;i<gridIndex.size();i++)
        {
            gridIndex2D.push_back(std::vector<IJVec>());
            for (int j=0;j<gridIndex[i].size();j++)
            {
                for (int k=0;k<gridIndex[i][j].size();k++)
                {
                    IJVec index;
                    index.i = gridIndex[i][j][k].i;
                    index.j = gridIndex[i][j][k].k; 
                    gridIndex2D[i].push_back(index);
                }
            }
        }
        
        for (std::map<int, IJKVec>::iterator mitr_cell=cellIndex.begin();mitr_cell!=cellIndex.end();mitr_cell++)
        {
            cellIndex2D[mitr_cell->first].i = mitr_cell->second.i;
            cellIndex2D[mitr_cell->first].j = mitr_cell->second.k;
        }
        
        std::map<int, std::vector<std::vector<double> > > C2D;
        std::pair<std::map<int, std::vector<std::vector<double> > >, std::vector<std::vector<double> > > diffusionResults;
        diffusionResults = DiffusionReactionCalculation2D(C_0,DiffC,reactionRate,Rt,delta_T,T);
        C2D = diffusionResults.first;
        
        std::vector<std::vector<std::vector<double> > > C;
        for (int i=0;i<gridIndex.size();i++)
        {                  
            C.push_back(std::vector<std::vector<double> >());
            for (int j=0;j<gridIndex[i].size();j++)
            {   
                C[i].push_back(std::vector<double>());
                for (int k=0;k<gridIndex[i][j].size();k++)
                {
                    C[i][j].push_back(0);
                }
            }
        }
             
        for (std::map<int, std::vector<std::vector<double> > >::iterator mitr_time=C2D.begin();\
            mitr_time!=C2D.end();mitr_time++
        )
        {
            concentration[mitr_time->first] = C;
        }
   
        for (std::map<int, std::vector<std::vector<double> > >::iterator mitr_time=C2D.begin();\
            mitr_time!=C2D.end();mitr_time++
        )
        {
            for (int i=0;i<mitr_time->second.size();i++)
            {
                for (int j=0;j<mitr_time->second[i].size();j++)
                {
                    concentration[mitr_time->first][i][0][j] = mitr_time->second[i][j];
                }
            }
        }
        
        std::vector<std::vector<double> > integralC = diffusionResults.second;
        for (int i=0;i<gridIndex.size();i++)
        {                  
            integralConcentration.push_back(std::vector<std::vector<double> >());
            for (int j=0;j<gridIndex[i].size();j++)
            {   
                integralConcentration[i].push_back(std::vector<double>());
                for (int k=0;k<gridIndex[i][j].size();k++)
                {
                    integralConcentration[i][j].push_back(0);
                }
            }
        }
        for (int i =0;i<integralC.size();i++)
        {
            for (int j=0;j<integralC[i].size();j++)
            {
                integralConcentration[i][0][j] = integralC[i][j];
            }
        }
    }

}


void ReactionDiffusionSimulation::WriteConcentrationProfile(string path)
{
    cout<<endl<<"WRTING CALCULATED CONCENTRATION TO FILE:"<<endl;
    for (std::map<int, std::vector<std::vector<std::vector<double> > > >::iterator mitr_time = concentration.begin();\
        mitr_time!=concentration.end(); mitr_time++)// Loop all the time steps
    {      
        stringstream ss;
        ss<<mitr_time->first;
        string filename;
        filename = "concentration-"+ss.str()+".csv";
        string filePath = path + filename;
        cout<<"filePath is "<<filePath<<endl;
        ofstream file;
        file.open(filePath.c_str());
        for (int i=0;i<gridIndex.size();i++)
        {
            for (int j=0;j<gridIndex[i].size();j++)
            {
                for (int k=0;k<gridIndex[i][j].size();k++)
                {
    //                 cout<<i<<" "<<j<<" "<<k<<" : "<<Concentration[nn][i][j][k]<<endl;
//                     file<<i<<","<<j<<","<<k<<","<<concentration[mitr_time->first][i][j][k]<<endl;
                }
            }
        }  
    }

}

vector< vector< vector< double > > > ReactionDiffusionSimulation::GetConcentration()
{
    std::map<int, std::vector<std::vector<std::vector<double> > > >::reverse_iterator rmitr_time=concentration.rbegin();
    return rmitr_time->second; // return the last element which is the concentration at the final time step
}

vector< vector< vector< double > > > ReactionDiffusionSimulation::GetIntegralConcentration()
{
    return integralConcentration;
}


int ReactionDiffusionSimulation::GetMeshDIMX()
{
    return gridNumX;
}
int ReactionDiffusionSimulation::GetMeshDIMY()
{
    return gridNumY;
}
int ReactionDiffusionSimulation::GetMeshDIMZ()
{
    return gridNumZ;
}