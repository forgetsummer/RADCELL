#include "DBSCAN.hh"
#include <iostream>
#include <math.h>
using namespace std;
void DBSCAN::testfunction()
{
    cout<<"function for testing class"<<endl;
/* 	cout<<"The database for clustering is as below:"<<endl;
    for (DataBaseMap::iterator mitr=dataBase.begin();mitr!=dataBase.end();mitr++ )
    {
            for (int i=0;i<(mitr->second).size();i++)
            {
                    cout<<(mitr->second)[i]<<",";
            }
            cout<<endl;
    } */
    cout<<"The cluster information is as below:"<<endl;
    cout<<"The number of clusters is "<<clusterID<<endl;
    std::map< int, int > number;
    for (std::map<int, int>::iterator mitr_cluster=clusterIDMap.begin();mitr_cluster!=clusterIDMap.end();mitr_cluster++)
    {
            number[mitr_cluster->second]=number[mitr_cluster->second]+1;
    }
    cout<<"The number for each cluster is: "<<endl;
    for (int i=0;i<number.size();i++)
    {
            cout<<"The number of points in cluster "<<i<<" is "<<number[i]<<endl;
            
    }

}
std::map<int, std::vector<int > > DBSCAN::GetClusterMap()
{
    std::map<int, std::vector<int> > clusterMap;
    for (std::map<int, int>::iterator mitr_cluster=clusterIDMap.begin();mitr_cluster!=clusterIDMap.end();mitr_cluster++)
    {
            //number[mitr_cluster->second]=number[mitr_cluster->second]+1;
            (clusterMap[mitr_cluster->second]).push_back(mitr_cluster->first);
    }
    return clusterMap;
    
}
void DBSCAN::ImportPoints(std::vector<std::vector<double> >& points) // function for importing database for clustering
{
    for (int i=0;i<points.size();i++) // loop line
    {
            for (int j=0;j<(points[i]).size();j++)// loop coloumn
            {
                    (dataBase[i]).push_back((points[i])[j]);
            }
    }

}
void DBSCAN::MarkClusterID(int index,int clID)
{
    clusterIDMap[index]=clID;
}
void DBSCAN::MarkVisitLabel(int index,int visit)
{
    visitLabelMap[index]=visit;
}
void DBSCAN::DBSCANClustering(double eps, int MinPts,double p)
{
    clusterID=0;// set clusterID as 0 initially
    for (std::map<int, int>::iterator mitr_cluster=clusterIDMap.begin();mitr_cluster!=clusterIDMap.end();mitr_cluster++)
    {
            MarkClusterID(mitr_cluster->first,0); // the initial cluster ID is 0, which means it is noise point
    }
    for (std::map<int, int>::iterator mitr_visit=visitLabelMap.begin();mitr_visit!=visitLabelMap.end();mitr_visit++)
    {
            MarkVisitLabel(mitr_visit->first, 0); // the initial visit label is 0,which means the points are not visited
    }
    for (DataBaseMap::iterator mitr_Pt=dataBase.begin();mitr_Pt!=dataBase.end();mitr_Pt++)
    {
            if (visitLabelMap[mitr_Pt->first]==0)
            {
                    MarkVisitLabel(mitr_Pt->first,1); // mark P as visited
                    DataBaseMap neighborPts;
                    neighborPts=regionQuery(mitr_Pt->first,dataBase,eps,p);
                    if (neighborPts.size()<MinPts)
                    {
                            MarkClusterID(mitr_Pt->first,0);// mark P as noise
                    }
                    else
                    {
                            clusterID=clusterID+1;
                            expandCluster(mitr_Pt->first,neighborPts,eps,MinPts,p);
                    }	
            }
    }
            
}
void DBSCAN::expandCluster(int index, DataBaseMap& neighborPts,double eps,int MinPts,double p)
{
    MarkClusterID(index,clusterID);// add P to cluster C
    bool expand_able=true;
    while (expand_able)
    {
            int initalSize=0;
            initalSize=neighborPts.size();
            for (DataBaseMap::iterator mitr_neighb=neighborPts.begin();mitr_neighb!=neighborPts.end();mitr_neighb++)
            {	
                    if (visitLabelMap[mitr_neighb->first]==0)// check new points in this neighborPts set
                    {
                            MarkVisitLabel(mitr_neighb->first,1);//mark p' as visited
                            DataBaseMap neighborPts_add;
                            neighborPts_add=regionQuery(mitr_neighb->first,dataBase,eps,p);
                            if (neighborPts_add.size()>=MinPts)
                            {
                                    for (DataBaseMap::iterator mitr_neighb_add=neighborPts_add.begin();mitr_neighb_add!=neighborPts_add.end();mitr_neighb_add++)
                                    {
                                            neighborPts[mitr_neighb_add->first]=mitr_neighb_add->second; // update neighborPts
                                    }
                            }
                    }
                    if (clusterIDMap[mitr_neighb->first]==0)// if P' is not yet member of any cluster
                    {
                            MarkClusterID(mitr_neighb->first,clusterID);// add p' to cluster with ID as clusterID
                    }
            }
            if (neighborPts.size()>initalSize)
            {
                    expand_able=true;
            }
            else
            {
                    expand_able=false;
            }
    }

}
DBSCAN::DataBaseMap DBSCAN::regionQuery(int index,DataBaseMap& DataBase,double eps,double p)
{

    DataBaseMap neighborPts;
    for  (DataBaseMap::iterator mitr=DataBase.begin();mitr!=DataBase.end();mitr++)
    {
            double distance=0;
            distance=Minkowskidistance(DataBase[index],mitr->second,p);// calulate Euclidean distance
            if (distance<eps)
            {
                    neighborPts[mitr->first]=mitr->second; // add this points to neighborPts 
            }	
    }
    return neighborPts; // return neighborPts
}
double DBSCAN::Minkowskidistance(std::vector<double>& P1,std::vector<double>& P2, double p )// function to calculate Minkowski distance between two points
{
    double sum=0;
    for (int i=0;i<P1.size();i++)
    {
            sum=sum+pow((fabs(P1[i]-P2[i])),p);
    }
    return pow(sum,1/p);
} 
