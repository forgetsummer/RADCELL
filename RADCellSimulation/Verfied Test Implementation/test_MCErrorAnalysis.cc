#include "RADCellSimulation.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include "G4RegionStore.hh"
#include "CellLayoutInitializer.hh"
#include "CellDoseAnalysis.hh"
#include "CellDNADamageAnalysis.hh"
#include <fstream>
#include <time.h>
#include "ReadRadiationTransportInfo.hh"
#include "ReactionDiffusionSimulation.hh"
#include "DiffusionReactionSolver.hh"
#include "CellStateModel.hh"
#include "RadiationSourceSettings.hh"
#include "DataProcessUsingPython.hh"
#include "math.h"
#include <sstream>
#include <stdlib.h>
#include "MonteCarloErrorAnalysis.hh"
double Rand2()
{
    int idum=12;// could change idum in here, here we just set it as 12
    const unsigned long IM1=2147483563,IM2=2147483399;
    const unsigned long IA1=40014,IA2=40692,IQ1=53668,IQ2=52774;
    const unsigned long IR1=12211,IR2=3791,NTAB=32,IMM1=IM1-1;
    const unsigned long NDIV=1+IMM1/NTAB;
    const double EPS=3.0e-16,AM=1.0/IM1,RNMX=(1.0-EPS);
    static int iy=0,idum2=314159269;
    static vector<int> iv(NTAB);
    int j,k;
    double temp;
    
    if ( idum <=0 )
        {
        idum=(idum ==0 ? 1 : -idum);
        idum2=idum;
        for ( j=NTAB+7; j>=0; j--) 
                {
            k=idum/IQ1;
            idum=IA1*(idum-k*IQ1)-k*IR1;
            if (idum < 0) idum+=IM1;
            if (j < NTAB) iv[j] = idum;
        }
        iy=iv[0];
    }
    k=idum/IQ1;
    idum=IA1*(idum-k*IQ1)-k*IR1;
    if (idum < 0) idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 +=IM2;
    j = iy/NDIV;
    iy = iv[j]-idum2;
    iv[j] = idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy)>RNMX ) return RNMX;
    else return temp;
}
int main(int argc,char** argv)
{
    
    MonteCarloErrorAnalysis testAnalysis;
    std::vector<std::vector<double> > data;
    int m=10000;
    int n=30;
    for (int i=0;i<m;i++)
    {
        data.push_back(std::vector<double>());
        for (int j=0;j<n;j++)
        {
            data[i].push_back(Rand2());
        }
    }
    
    for (int i= 0;i<data.size();i++)
    {
        for (int j=0;j<data[i].size();j++)
        {
            cout<<"data"<<i<<","<<j<<" = "<<data[i][j]<<endl;
        }
    }
    
    testAnalysis.ImportSimulationResults(data);
    std::vector<double> mean;
    mean = testAnalysis.GetMeanOfSimulationResults();
    for (int i=0;i<mean.size();i++)
    {
        cout<<"mean = "<<mean[i]<<endl;
    }
    std::vector<double> std;
    std = testAnalysis.GetStdOfSimulatioinResults();
    for (int i=0;i<mean.size();i++)
    {
        cout<<"std = "<<std[i]<<endl;
    }
    
    return 0;
}