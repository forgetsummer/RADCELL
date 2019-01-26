#include "MonteCarloErrorAnalysis.hh"
#include <math.h>
#include<iostream>
using namespace std;
void MonteCarloErrorAnalysis::ImportSimulationResults(vector< vector< double > > result)
{
    simResult = result;
}

vector< double > MonteCarloErrorAnalysis::GetMeanOfSimulationResults()
{
//     cout<<"size of simResult is "<<simResult.size()<<endl;
//     for (int i=0;i<simResult.size();i++)
//     {
//         for (int j=0;j<simResult[i].size();j++)
//         {
//             cout<<"simResult"<<i<<","<<j<<" = "<<simResult[i][j]<<endl;
//         }
//     }
//     cout<<"ssize "<<simResult[0].size()<<endl;
    std::vector<double> sum;
   
    for (int i=0;i<simResult[0].size();i++)
    {
        sum.push_back(0);
    }    
    for (int i=0;i<simResult.size();i++) // loop all the simulation histories
    {
        for (int j=0;j<simResult[i].size();j++)
        {
            sum[j] = sum[j] + simResult[i][j];
        }
    }

    std::vector<double> mean;
    for (int i=0;i<simResult[0].size();i++)
    {
        mean.push_back(0);
    }
    
    for (int i =0;i<sum.size();i++)
    {
        mean[i] = sum[i]/simResult.size();
    }
    return mean;
}

vector< double > MonteCarloErrorAnalysis::GetStdOfSimulatioinResults()
{
    std::vector<double> sum;
    for (int i=0;i<simResult[0].size();i++)
    {
        sum.push_back(0);
    }
    
    std::vector<double> mean = GetMeanOfSimulationResults();
    for (int j=0;j<simResult[0].size();j++)
    {
        for (int i=0;i<simResult.size();i++)
        {
            sum[j] = sum[j] + (simResult[i][j]-mean[j])*(simResult[i][j]-mean[j]);
        }
    }
    std::vector<double>std;
    std::vector<double> statisticalStd;
    for (int i =0;i<mean.size();i++)
    {
        std.push_back(0);
        statisticalStd.push_back(0);
    }
    for (int i=0;i<mean.size();i++)
    {
        std[i] = sqrt(sum[i]/(simResult.size()-1));
        statisticalStd[i] = std[i]*1.96/sqrt(simResult.size());
    }
    return statisticalStd;
//     return std;
}






