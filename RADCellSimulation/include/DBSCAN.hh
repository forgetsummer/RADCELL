#ifndef DBSCAN_h
#define DBSCAN_h 1
#include <map>
#include <vector>
class DBSCAN
{
    public:
        void ImportPoints(std::vector<std::vector<double> >& points);// position datas stored as a matrix
        void DBSCANClustering(double eps, int MinPts,double p);// eps is the minimum radius of cluster, MinPts is minimum points of cluster, p is order for Minkowski distance 
        void testfunction();
        int GetClusterNumber(){return clusterID;}// get cluster number
        int testNumber(){return 1;}
        std::map<int, std::vector<int > > GetClusterMap();// get cluster map, key is clusterID, element is points index in each cluster
        
        
    private:
        typedef std::vector<double> Coordinate;
        typedef std::map<int, Coordinate> DataBaseMap;
        void MarkClusterID(int index,int clID);
        void MarkVisitLabel(int index,int visit);
        double Minkowskidistance(std::vector<double>& P1,std::vector<double>& P2, double p );// function to calculate Minkowski distance between two points
        DataBaseMap regionQuery(int index,DataBaseMap& D,double eps,double p);// return a cluster region
        void expandCluster(int index, DataBaseMap& neighborPts, double eps,int MinPts,double p);
        DataBaseMap dataBase; // a map named DataMap to store all the points for clustering, key is point index
        std::map<int, int>clusterIDMap;//a map named ClusterIDMap to store cluster ID of all points, key is clusterID
        std::map<int, int>visitLabelMap;// a map to store vist label for all points, key is point index
        int clusterID;

	   
};
#endif