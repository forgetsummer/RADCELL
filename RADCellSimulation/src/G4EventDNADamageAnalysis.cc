#include "G4EventDNADamageAnalysis.hh"

#include "DBSCAN.hh"
#include <iostream>
#include<fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include "Randomize.hh"  // uisng random number


EventDNADamagesTallyInfo G4EventDNADamageAnalysis::EventDNADamageTally(EnergyDepositionEventMap edepEventMap, double SSBProb, double SSBEmin, double SSBEmax, double eps, int MinPts, double p)
{
    EventEdepMap edepMap;
    
    for (EnergyDepositionEventMap::iterator mitr_event=edepEventMap.begin();mitr_event!=edepEventMap.end();mitr_event++)
    {

        for (int i=0;i<mitr_event->second.size();i++)
        {
            EnergyDepositionPoint edepPoint;
            double edep=((mitr_event->second)[i]).edep;
            double X=((mitr_event->second)[i]).pX;
            double Y=((mitr_event->second)[i]).pY;
            double Z= ((mitr_event->second)[i]).pZ;
            string affectedCellOrganelle=((mitr_event->second)[i]).affectedCellOrganelle;
            edepPoint.x=X;
            edepPoint.y=Y;
            edepPoint.z=Z;
            edepPoint.edep=edep;
            affectedCellOrganelle=mitr_event->second[i].affectedCellOrganelle;
            if (affectedCellOrganelle=="Nucleus")
            {
                edepMap[mitr_event->first][mitr_event->second[i].cellID].push_back(edepPoint); //only store the energy deposition in nuclues 
            }
            
        }

    }
    
  
        for (EventEdepMap::iterator mitr_event=edepMap.begin();mitr_event!=edepMap.end();mitr_event++)
        {
            if (mitr_event->first==0) // print the edep  point information of first event in cell nucleus
            {
                ofstream file;

                file.open("eventID_0_edepInfoInNucleus.csv");
                for (CellEdepMap::iterator mitr_cell=mitr_event->second.begin();mitr_cell!=mitr_event->second.end();mitr_cell++)
                {
                    for (int i=0;i<mitr_cell->second.size();i++)
                    {
                        file<<"EventID"<<","<<mitr_event->first<<","<< mitr_cell->second[i].x<<","<<mitr_cell->second[i].y<<","<<mitr_cell->second[i].z<<","<<mitr_cell->second[i].edep<<endl;
                    }
                }
            }

    }

    
    int N=1000; // the number for determining the times of processing the edep points data
//     int N=10; // this reprocessing number was used to do DSB quantification test
    double probOfInSensitiveRegion=SSBProb;
    
//     for (EventEdepMap::iterator mitr_event=edepMap.begin();mitr_event!=edepMap.end();mitr_event++)
//     {
//         cout<<"The eventID now in G4EventDNADamageAnalysis is :"<<mitr_event->first<<endl;
//         for (CellEdepMap::iterator mitr_cell=mitr_event->second.begin();mitr_cell!=mitr_event->second.end();mitr_cell++)
//         {
//             cout<<"The cellID for this event is: "<<mitr_cell->first<<endl;
//             cout<<"And the edep points number for this cell is: "<<mitr_cell->second.size()<<endl;
//         }
//     }
    
    for (int i=0;i<N;i++)
    {
        EventDNASSBMap potentialSSBTallyMap;
        potentialSSBTallyMap=EventSSBTally(edepMap,SSBEmin, SSBEmax);

        EventDSBTally(potentialSSBTallyMap,probOfInSensitiveRegion,eps, MinPts,p);
   
    }
    
    

    EventDNADamagesTallyInfo currentDNADamageTallyInfo;
    
    currentDNADamageTallyInfo.eventDSBMap=DSBTallyMap;
    currentDNADamageTallyInfo.eventSSBMap=SSBTallyMap;
    currentDNADamageTallyInfo.N=N;
    
    return currentDNADamageTallyInfo;
    

}


EventDNASSBMap G4EventDNADamageAnalysis::EventSSBTally(EventEdepMap edepMap, double SSBEmin, double SSBEmax)
{
    EventDNASSBMap potentialSSBTallyMap;// a map storing the potential SSB which is made by the edep points whose energy is meeting the energy requirment for SSB
    for (EventEdepMap::iterator mitr_event=edepMap.begin();mitr_event!=edepMap.end();mitr_event++)// looping the eventIDs
    {
        for(CellEdepMap::iterator mitr_cell=(mitr_event->second).begin();mitr_cell!=(mitr_event->second).end();mitr_cell++)// looping the cellIDs
        {
            for (int i=0;i<(mitr_cell->second).size();i++)// looping the vector storing the edep points in one cell
            {

                double p= (((mitr_cell->second)[i]).edep-SSBEmin)/(SSBEmax-SSBEmin);
                
                    if (G4UniformRand()<p)
//                 if (((mitr_cell->second)[i]).edep>10.79)
                    {
                        int DNALabel;
                        SSB possibleSSB;
                        if (G4UniformRand()<0.5)
                        {
                            DNALabel=0;
                        }
                        else
                        {
                            DNALabel=1;   
                        }

                        possibleSSB.edepPoint=((edepMap[mitr_event->first])[mitr_cell->first])[i];
                        possibleSSB.DNALabel=DNALabel;

                        (( potentialSSBTallyMap[mitr_event->first])[mitr_cell->first]).push_back(possibleSSB);

                        
                    }

            }
                
        }
    }
    
    return potentialSSBTallyMap;

}

void G4EventDNADamageAnalysis::EventDSBTally(EventDNASSBMap& potentialSSBTallyMap, double probOfInSensitiveRegion, double eps, int MinPts, double p)
{
    for (EventDNASSBMap::iterator mitr_event=potentialSSBTallyMap.begin();mitr_event!=potentialSSBTallyMap.end();mitr_event++)
    {
//         cout<<"Now processing eventID "<<mitr_event->first<<endl;
        for (CellDNASSBMap::iterator mitr_cell=(mitr_event->second).begin();mitr_cell!=(mitr_event->second).end();mitr_cell++)
        {
            DBSCAN DSB_DBSCAN; // declare a DBSCAN object for DSB clustering
            std::vector< std::vector<double> > pointPosition; //  this is the SSB points for clustring to get DSBs 
            for (int i=0; i<(mitr_cell->second).size();i++)// loop the vector storing SSB objects
            {
                std::vector<double> position;
                double x =((mitr_cell->second)[i]).edepPoint.x;
                double y = ((mitr_cell->second)[i]).edepPoint.y;
                double z = ((mitr_cell->second)[i]).edepPoint.z;
                position.push_back(x);
                position.push_back(y);
                position.push_back(z);
                pointPosition.push_back(position);

            }
            DSB_DBSCAN.ImportPoints(pointPosition);
            DSB_DBSCAN.DBSCANClustering(eps,MinPts,p);// eps in 3.2nm, MinPts is 2, using Euclidean distance for clustering
            std::map<int, std::vector<int> > SSBClusterMap;
            SSBClusterMap=DSB_DBSCAN.GetClusterMap();
            SSBClusterMap.erase(0);// delte the noise points, since the DSB is formed for SSB clusters , WARNING, be careful in here, still do not know whether should be delete it or not, but it seems that delete it will be better

            if (G4UniformRand()<probOfInSensitiveRegion) // if the cluster is in sensitive region, then it has the probablity of inducing DSB, otherwise none.
            {
//                 cout<<"enter test whether it is in sensitive region"<<endl;
                DSBVector DSBVec=GetDSBCluster(SSBClusterMap,mitr_cell->second);
                for (int i=0;i<DSBVec.size();i++)
                {
                    ((DSBTallyMap[mitr_event->first])[mitr_cell->first]).push_back(DSBVec[i]); //push_back DSB 
                }
//                 cout<<"The size of DSBTallyMap in  EventDSBTally is "<<DSBTallyMap.size()<<endl;
                for (int i=0;i<mitr_cell->second.size();i++)
                {
                    ((SSBTallyMap[mitr_event->first])[mitr_cell->first]).push_back((mitr_cell->second)[i]); // push_back SSB
                }
//                 cout<<"The size of SSBTallyMap in EventDSBTally is "<<SSBTallyMap.size()<<endl;
                
            }
                       
        }
    }

}

DSBVector G4EventDNADamageAnalysis::GetDSBCluster(map< int, vector< int > >& SSBClusterMap, SSBVector& SSBVec)
{
    DSBVector DSBVec;
    for (std::map<int, std::vector<int> >::iterator mitr_SSBCluster=SSBClusterMap.begin();mitr_SSBCluster!=SSBClusterMap.end();mitr_SSBCluster++)
    {
            int strandSum=0;
            for (int i=0; i<(mitr_SSBCluster->second).size();i++)
            {
                    int II=(mitr_SSBCluster->second)[i];//get the point index in SSBVector
                    strandSum=strandSum+(SSBVec[II]).DNALabel;
            }

            if (strandSum>0 && strandSum<(mitr_SSBCluster->second).size())
            {
                DSB theDSB;
                
                double DSBTotalEnergy=0;
                double tX=0;
                double tY=0;
                double tZ=0;
                double xC;
                double yC;
                double zC;
                double DSBDimension;
                int DSBType;
                for (int i=0; i<(mitr_SSBCluster->second).size();i++)
                {
                        int II=(mitr_SSBCluster->second)[i];
                        DSBTotalEnergy=DSBTotalEnergy+(SSBVec[II]).edepPoint.edep;
                        tX=tX+(SSBVec[II]).edepPoint.x;
                        tY=tY+(SSBVec[II]).edepPoint.y;
                        tZ=tZ+(SSBVec[II]).edepPoint.z;
                }
                xC=tX/(mitr_SSBCluster->second).size();
                yC=tY/(mitr_SSBCluster->second).size();
                zC=tZ/(mitr_SSBCluster->second).size();
//                 cout<<"the cluster size is "<<mitr_SSBCluster->second.size()<<endl;
//                 cout<<"the strandSum is "<<strandSum<<endl;
                DSBDimension=GetDSBDimension(mitr_SSBCluster->second,SSBVec);
                if ((mitr_SSBCluster->second).size()==2)
                {
                        DSBType=1;//it is simple DSB
                        
                }
                if ((mitr_SSBCluster->second).size()>2)
                {
                        DSBType=2; // it is complex DSB
                }
                theDSB.DSBTotalEnergy=DSBTotalEnergy;
                theDSB.xC=xC;
                theDSB.yC=yC;
                theDSB.zC=zC;
                theDSB.DSBDimension=DSBDimension;
                theDSB.DSBType=DSBType;
                DSBVec.push_back(theDSB);
            }   


    }
    return DSBVec;

}

double G4EventDNADamageAnalysis::GetDSBDimension(vector< int >& SSBPts, SSBVector& SSBVec)
{
    double maxDistance=0;
    std::vector<double> distance;
    for (int i=0; i<SSBPts.size();i++)
    {
        for (int j=0; j<SSBPts.size();j++)
        { 
            if (i!=j)
            {
                double d=0;
                double sum=0;
                double dx=0;
                double dy=0;
                double dz=0;
                int II=SSBPts[i];
                int JJ=SSBPts[j];
//                 cout<<"the ssb lable for i "<<SSBVec[II].DNALabel<<endl;
//                 cout<<"the ssb label for j "<<SSBVec[JJ].DNALabel<<endl;
                dx=((SSBVec[II]).edepPoint.x)-((SSBVec[JJ]).edepPoint.x);
                dy=((SSBVec[II]).edepPoint.y)-((SSBVec[JJ]).edepPoint.y);
                dz=((SSBVec[II]).edepPoint.z)-((SSBVec[JJ]).edepPoint.z);
                dx=pow(dx,2);
                dy=pow(dy,2);
                dz=pow(dz,2);
                sum=dx+dy+dz;
                d=pow(sum,0.5);
//                 cout<<"distance between i and j "<<d<<endl;
                distance.push_back(d);  
            }
            

        }
            
    }
    for (int k=0;k<distance.size();k++)
    {
        if (maxDistance<distance[k])
        {
            maxDistance=distance[k];
        }
    }
//     cout<<"the dsb distance is "<<maxDistance<<endl;
    return maxDistance;

}

double G4EventDNADamageAnalysis::Rand2()
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

