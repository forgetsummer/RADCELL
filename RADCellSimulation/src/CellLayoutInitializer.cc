#include "CellLayoutInitializer.hh"
#include <iostream>
#include <math.h> // in case of using floor() function 
#include <vector>
#include "Randomize.hh"  // uisng random number
#include <G4Types.hh>

using namespace std;

CellLayoutInitializer::CellLayoutInitializer()
{
    cellHomeSizeX=1;
    cellHomeSizeY=1;
    cellHomeSizeZ=1;
}
CellLayoutInitializer::~CellLayoutInitializer()
{
    
}
void CellLayoutInitializer::RectangularSlab(double xDim, double yDim, double zDim,int cellNumber)
{
    int N_X=ceil(xDim/cellHomeSizeX);
    int N_Y=ceil(yDim/cellHomeSizeY);
    int N_Z=ceil(zDim/cellHomeSizeZ);// calculating the total number of cell home
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
    std::vector<std::vector<std::vector<int> > > cellHomeResidence;
    int effectiveCellNumber;

    for (int k=0;k<N_Z;k++)
    {
        cellHomeResidence.push_back(std::vector<std::vector<int> >());
        for (int i=0;i<N_X;i++)
        {
            (cellHomeResidence[k]).push_back(std::vector<int>());
            for (int j=0;j<N_Y;j++)
            {
                (cellHomeResidence[k][i]).push_back(0);// intialize cellHomeResidence
            }
        }
    }

    if (cellNumber<=N_X*N_Y*N_Z)
    {
        effectiveCellNumber=cellNumber;
    }
    else
    {
        effectiveCellNumber=N_X*N_Y*N_Z;
    }
    
//     cout<<"the nx*ny*nz is "<<N_X*N_Y*N_Z<<endl;
    
//     cout<<"the effectiveCellNumber now is "<<effectiveCellNumber<<endl;
    for (int i=0; i<effectiveCellNumber;i++)
    {
        Index tempIndex;
        tempIndex.i=0;
        tempIndex.j=0;
        tempIndex.k=0;
        positionIndex.push_back(tempIndex); //initialize positionIndex
        cellPositionX.push_back(0);
        cellPositionY.push_back(0);
        cellPositionZ.push_back(0);
    }
    
    
//     cout<<"the size of cellPositionX now is "<<cellPositionX.size()<<endl;

    int I=0;
    while (I<effectiveCellNumber)
    { 
        int it;
        int jt;
        int kt;
        if (N_X==1)
        {
            it=0;
        }
        else
        {
            it=G4UniformRand()*N_X;
        }
        if (N_Y==1)
        {
            jt=0;
        }
        else
        {
           jt=G4UniformRand()*N_Y;
        }
        if (N_Z==1)
        {
            kt=0;
        }
        else
        {
           kt=G4UniformRand()*N_Z;
        }
        
        if (cellHomeResidence[kt][it][jt]==0)
        {
            cellHomeResidence[kt][it][jt]=1;
            Index tempIndex;
            tempIndex.i=it;
            tempIndex.j=jt;
            tempIndex.k=kt;
            positionIndex[I]=tempIndex; // asign position index for positionIndex vector
            I=I+1;  
        } 

    }
    for (int n=0; n<effectiveCellNumber;n++)
    {
        double x;
        double y;
        double z;
        if (xDim==0) // when xDim is 0, just set up x as 0, this is 2D case for y-z plane
        {
            x = 0;
        }
        else 
        {
            x=-1*xDim/2+((positionIndex[n]).i)*cellHomeSizeX+cellHomeSizeX/2;
        }
        if (yDim ==0)//this is 2D case for x-z plane
        {
            y = 0;
        }
        else 
        {
            y=-1*yDim/2+((positionIndex[n]).j)*cellHomeSizeY+cellHomeSizeY/2;
        }
        if (zDim == 0)//this is 2D case for x-y plane
        {
            z = 0;
        }
        else
        {
            z=-1*zDim/2+((positionIndex[n]).k)*cellHomeSizeZ+cellHomeSizeZ/2;
        }
        cellPositionX[n]=x; // asign x for cellPositionX
        cellPositionY[n]=y; // asign y for cellPositionY
        cellPositionZ[n]=z; // asign z for cellPositionZ
    }
    
}

void CellLayoutInitializer::Blob(double dDim, int seededCellNumber)
{
    double xDim;
    double yDim;
    double zDim;
    xDim=dDim;
    yDim=dDim;
    zDim=dDim;
    
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
    std::vector<std::vector<std::vector<int> > > cellHomeResidence;
    
    for (int k=0;k<N_Z;k++)
    {
        cellHomeResidence.push_back(std::vector<std::vector<int> >());
        for (int i=0;i<N_X;i++)
        {
            (cellHomeResidence[k]).push_back(std::vector<int>());
            for (int j=0;j<N_Y;j++)
            {
                (cellHomeResidence[k][i]).push_back(0);// intialize cellHomeResidence
            }
        }
    }

    int I=0;
    int totalCellInsideRegion = 0;
    for (int i =0;i<N_X;i++)
    {
        for (int j =0;j<N_Y;j++) 
        {
            for (int k =0;k<N_Z;k++)
            {
                double x;
                double y;
                double z;
                x=-1*xDim/2+i*cellHomeSizeX+cellHomeSizeX/2;
                y=-1*yDim/2+j*cellHomeSizeY+cellHomeSizeY/2;
                z=-z*zDim/2+k*cellHomeSizeZ+cellHomeSizeZ/2;
                if (x*x+y*y+z*z<0.25*dDim*dDim)
                {
                    totalCellInsideRegion = totalCellInsideRegion + 1;     
                }
            } 
        }
    }
    
    int effectiveCellNumber;
    if (seededCellNumber<totalCellInsideRegion)
    {
        effectiveCellNumber = seededCellNumber;
    }
    else
    {
        effectiveCellNumber = totalCellInsideRegion;
    }
    
    for (int i=0; i<effectiveCellNumber;i++)
    {
        Index tempIndex;
        tempIndex.i=0;
        tempIndex.j=0;
        tempIndex.k=0;
        positionIndex.push_back(tempIndex); //initialize positionIndex
        cellPositionX.push_back(0);
        cellPositionY.push_back(0);
        cellPositionZ.push_back(0);
    }
    
    int count=0;
    while (count<effectiveCellNumber)
    {
        int it;
        int jt;
        int kt;
        if (N_X==1)
        {
            it=0;
        }
        else
        {
            it=G4UniformRand()*N_X;
        }
        if (N_Y==1)
        {
            jt=0;
        }
        else
        {
           jt=G4UniformRand()*N_Y;
        }
        if (N_Z==1)
        {
            kt=0;
        }
        else
        {
           kt=G4UniformRand()*N_Z;
        }
        
        double x;
        double y;
        double z;
        x=-1*xDim/2+it*cellHomeSizeX+cellHomeSizeX/2;
        y=-1*yDim/2+jt*cellHomeSizeY+cellHomeSizeY/2;
        z=-1*zDim/2+kt*cellHomeSizeZ+cellHomeSizeZ/2;
        if ((cellHomeResidence[kt][it][jt]==0)&&(x*x+y*y+z*z<0.25*dDim*dDim))
        {
            cellHomeResidence[kt][it][jt]=1;
            Index tempIndex;
            tempIndex.i=it;
            tempIndex.j=jt;
            tempIndex.k=kt;
            positionIndex[count]=tempIndex; // asign posiiton index for positionIndex vector
            count=count+1;
        } 
        
    }

    for (int n=0; n<count;n++)
    {
        double x;
        double y;
        double z;
        x=-1*xDim/2+((positionIndex[n]).i)*cellHomeSizeX+cellHomeSizeX/2;
        y=-1*yDim/2+((positionIndex[n]).j)*cellHomeSizeY+cellHomeSizeY/2;
        z=-1*zDim/2+((positionIndex[n]).k)*cellHomeSizeZ+cellHomeSizeZ/2;
        cellPositionX[n]=x; // asign x for cellPositionX
        cellPositionY[n]=y; // asign y for cellPositionY
        cellPositionZ[n]=z; // asign z for cellPositionZ
    }

}


void CellLayoutInitializer::Blob2D(double dDim, int seededCellNumber)
{
    double xDim;
    double yDim;
    double zDim;
    xDim=dDim;
    yDim=dDim;
    zDim=0;
    
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
    std::vector<std::vector<std::vector<int> > > cellHomeResidence;
    
    for (int k=0;k<N_Z;k++)
    {
        cellHomeResidence.push_back(std::vector<std::vector<int> >());
        for (int i=0;i<N_X;i++)
        {
            (cellHomeResidence[k]).push_back(std::vector<int>());
            for (int j=0;j<N_Y;j++)
            {
                (cellHomeResidence[k][i]).push_back(0);// intialize cellHomeResidence
            }
        }
    }
    
    int I=0;
    int totalCellInsideRegion = 0;
    for (int i =0;i<N_X;i++)
    {
        for (int j =0;j<N_Y;j++) 
        {
            double x;
            double y;
            double z;
            x=-1*xDim/2+i*cellHomeSizeX+cellHomeSizeX/2;
            y=-1*yDim/2+j*cellHomeSizeY+cellHomeSizeY/2;
            if (x*x+y*y<0.25*dDim*dDim)
            {
                totalCellInsideRegion = totalCellInsideRegion + 1; 
            }
            
        }
    }
    
    int effectiveCellNumber;
    if (seededCellNumber<totalCellInsideRegion)
    {
        effectiveCellNumber = seededCellNumber;
    }
    else
    {
        effectiveCellNumber = totalCellInsideRegion;
    }
    
    for (int i=0; i<effectiveCellNumber;i++)
    {
        Index tempIndex;
        tempIndex.i=0;
        tempIndex.j=0;
        tempIndex.k=0;
        positionIndex.push_back(tempIndex); //initialize positionIndex
        cellPositionX.push_back(0);
        cellPositionY.push_back(0);
        cellPositionZ.push_back(0);
    }
    
    int count=0;
    while (count<effectiveCellNumber)
    {
        int it;
        int jt;
        int kt;
        if (N_X==1)
        {
            it=0;
        }
        else
        {
            it=G4UniformRand()*N_X;
        }
        if (N_Y==1)
        {
            jt=0;
        }
        else
        {
           jt=G4UniformRand()*N_Y;
        }
        if (N_Z==1)
        {
            kt=0;
        }
        else
        {
           kt=G4UniformRand()*N_Z;
        }
        
        double x;
        double y;
        double z;
        x=-1*xDim/2+it*cellHomeSizeX+cellHomeSizeX/2;
        y=-1*yDim/2+jt*cellHomeSizeY+cellHomeSizeY/2;
        z=0; // since this 2D case, so just set up z as 0, this will give a cicle in x-y plane
        if ((cellHomeResidence[kt][it][jt]==0)&&(x*x+y*y<0.25*dDim*dDim))
        {
            cellHomeResidence[kt][it][jt]=1;
            Index tempIndex;
            tempIndex.i=it;
            tempIndex.j=jt;
            tempIndex.k=kt;
            positionIndex[count]=tempIndex; // asign posiiton index for positionIndex vector
            count=count+1;
        } 
        
    }

    for (int n=0; n<count;n++)
    {
        double x;
        double y;
        double z;
        x=-1*xDim/2+((positionIndex[n]).i)*cellHomeSizeX+cellHomeSizeX/2;
        y=-1*yDim/2+((positionIndex[n]).j)*cellHomeSizeY+cellHomeSizeY/2;
        z=0;// since it is 2D case, so here just give circle in x-y plane
        cellPositionX[n]=x; // asign x for cellPositionX
        cellPositionY[n]=y; // asign y for cellPositionY
        cellPositionZ[n]=z; // asign z for cellPositionZ
    }

}


void CellLayoutInitializer::SetCellHomeParamter(double sizeX, double sizeY, double sizeZ)
{
    cellHomeSizeX=sizeX;
    cellHomeSizeY=sizeY;
    cellHomeSizeZ=sizeZ;
}

void CellLayoutInitializer::Heart(double dDim, int seededCellNumber)
{
    double xDim;
    double yDim;
    double zDim;
    xDim=dDim;
    yDim=dDim;
    zDim=dDim;
    double d = dDim/2*0.7; // 0.7 here is a scaling factor in case of the heart shape is larger than the world dimension
    
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
    std::vector<std::vector<std::vector<int> > > cellHomeResidence;
    
    for (int k=0;k<N_Z;k++)
    {
        cellHomeResidence.push_back(std::vector<std::vector<int> >());
        for (int i=0;i<N_X;i++)
        {
            (cellHomeResidence[k]).push_back(std::vector<int>());
            for (int j=0;j<N_Y;j++)
            {
                (cellHomeResidence[k][i]).push_back(0);// intialize cellHomeResidence
            }
        }
    }
    
    int I=0;
    int totalCellInsideRegion = 0;
    for (int i =0;i<N_X;i++)
    {
        for (int j =0;j<N_Y;j++) 
        {
            for (int k =0;k<N_Z;k++)
            {
                double x;
                double y;
                double z;
                x=-1*xDim/2+i*cellHomeSizeX+cellHomeSizeX/2;
                y=-1*yDim/2+j*cellHomeSizeY+cellHomeSizeY/2;
                z=-1*zDim/2+k*cellHomeSizeZ+cellHomeSizeZ/2;
                if (pow((x*x+y*y-d*d),3)-pow(x,2)*pow(y,3)<0)
                {
                    totalCellInsideRegion = totalCellInsideRegion + 1; 
                }
            }

            
        }
    }
    
    int effectiveCellNumber;
    if (seededCellNumber<totalCellInsideRegion)
    {
        effectiveCellNumber = seededCellNumber;
    }
    else
    {
        effectiveCellNumber = totalCellInsideRegion;
    }
    
    for (int i=0; i<effectiveCellNumber;i++)
    {
        Index tempIndex;
        tempIndex.i=0;
        tempIndex.j=0;
        tempIndex.k=0;
        positionIndex.push_back(tempIndex); //initialize positionIndex
        cellPositionX.push_back(0);
        cellPositionY.push_back(0);
        cellPositionZ.push_back(0);
    }
    
    int count=0;
    while (count<effectiveCellNumber)
    {
        int it;
        int jt;
        int kt;
        if (N_X==1)
        {
            it=0;
        }
        else
        {
            it=G4UniformRand()*N_X;
        }
        if (N_Y==1)
        {
            jt=0;
        }
        else
        {
           jt=G4UniformRand()*N_Y;
        }
        if (N_Z==1)
        {
            kt=0;
        }
        else
        {
           kt=G4UniformRand()*N_Z;
        }
        
        double x;
        double y;
        double z;
        x=-1*xDim/2+it*cellHomeSizeX+cellHomeSizeX/2;
        y=-1*yDim/2+jt*cellHomeSizeY+cellHomeSizeY/2;
        z=-1*zDim/2+kt*cellHomeSizeZ+cellHomeSizeZ/2;

        if (cellHomeResidence[kt][it][jt]==0&&(pow((x*x+y*y-d*d),3)-pow(x,2)*pow(y,3)<0))
//         if ((cellHomeResidence[kt][it][jt]==0)&&(x*x+y*y<0.25*dDim*dDim))
        {
            cellHomeResidence[kt][it][jt]=1;
            Index tempIndex;
            tempIndex.i=it;
            tempIndex.j=jt;
            tempIndex.k=kt;
            positionIndex[count]=tempIndex; // asign posiiton index for positionIndex vector
            count=count+1;
        } 
        
    }

    for (int n=0; n<count;n++)
    {
        double x;
        double y;
        double z;
        x=-1*xDim/2+((positionIndex[n]).i)*cellHomeSizeX+cellHomeSizeX/2;
        y=-1*yDim/2+((positionIndex[n]).j)*cellHomeSizeY+cellHomeSizeY/2;
        z=-1*zDim/2+((positionIndex[n]).k)*cellHomeSizeZ+cellHomeSizeZ/2;
        cellPositionX[n]=x; // asign x for cellPositionX
        cellPositionY[n]=y; // asign y for cellPositionY
        cellPositionZ[n]=z; // asign z for cellPositionZ
    }

}

void CellLayoutInitializer::Heart2D(double dDim, int seededCellNumber)
{
    double xDim;
    double yDim;
    double zDim;
    xDim=dDim;
    yDim=dDim;
    zDim=0;
    double d = dDim/2*0.7; // 0.7 here is a scaling factor in case of the heart shape is larger than the world dimension
    
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
    std::vector<std::vector<std::vector<int> > > cellHomeResidence;
    
    for (int k=0;k<N_Z;k++)
    {
        cellHomeResidence.push_back(std::vector<std::vector<int> >());
        for (int i=0;i<N_X;i++)
        {
            (cellHomeResidence[k]).push_back(std::vector<int>());
            for (int j=0;j<N_Y;j++)
            {
                (cellHomeResidence[k][i]).push_back(0);// intialize cellHomeResidence
            }
        }
    }
    
    int I=0;
    int totalCellInsideRegion = 0;
    for (int i =0;i<N_X;i++)
    {
        for (int j =0;j<N_Y;j++) 
        {
            double x;
            double y;
            double z;
            x=-1*xDim/2+i*cellHomeSizeX+cellHomeSizeX/2;
            y=-1*yDim/2+j*cellHomeSizeY+cellHomeSizeY/2;
            if (pow((x*x+y*y-d*d),3)-pow(x,2)*pow(y,3)<0)
            {
                totalCellInsideRegion = totalCellInsideRegion + 1; 
            }
            
        }
    }
    
    int effectiveCellNumber;
    if (seededCellNumber<totalCellInsideRegion)
    {
        effectiveCellNumber = seededCellNumber;
    }
    else
    {
        effectiveCellNumber = totalCellInsideRegion;
    }
    
    for (int i=0; i<effectiveCellNumber;i++)
    {
        Index tempIndex;
        tempIndex.i=0;
        tempIndex.j=0;
        tempIndex.k=0;
        positionIndex.push_back(tempIndex); //initialize positionIndex
        cellPositionX.push_back(0);
        cellPositionY.push_back(0);
        cellPositionZ.push_back(0);
    }
    
    int count=0;
    while (count<effectiveCellNumber)
    {
        int it;
        int jt;
        int kt;
        if (N_X==1)
        {
            it=0;
        }
        else
        {
            it=G4UniformRand()*N_X;
        }
        if (N_Y==1)
        {
            jt=0;
        }
        else
        {
           jt=G4UniformRand()*N_Y;
        }
        if (N_Z==1)
        {
            kt=0;
        }
        else
        {
           kt=G4UniformRand()*N_Z;
        }
        
        double x;
        double y;
        double z;
        x=-1*xDim/2+it*cellHomeSizeX+cellHomeSizeX/2;
        y=-1*yDim/2+jt*cellHomeSizeY+cellHomeSizeY/2;
        z=0; // since this 2D case, so just set up z as 0, this will give a cicle in x-y plane

        if (cellHomeResidence[kt][it][jt]==0&&(pow((x*x+y*y-d*d),3)-pow(x,2)*pow(y,3)<0))
//         if ((cellHomeResidence[kt][it][jt]==0)&&(x*x+y*y<0.25*dDim*dDim))
        {
            cellHomeResidence[kt][it][jt]=1;
            Index tempIndex;
            tempIndex.i=it;
            tempIndex.j=jt;
            tempIndex.k=kt;
            positionIndex[count]=tempIndex; // asign posiiton index for positionIndex vector
            count=count+1;
        } 
        
    }

    for (int n=0; n<count;n++)
    {
        double x;
        double y;
        double z;
        x=-1*xDim/2+((positionIndex[n]).i)*cellHomeSizeX+cellHomeSizeX/2;
        y=-1*yDim/2+((positionIndex[n]).j)*cellHomeSizeY+cellHomeSizeY/2;
        z=0;// since it is 2D case, so here just give circle in x-y plane
        cellPositionX[n]=x; // asign x for cellPositionX
        cellPositionY[n]=y; // asign y for cellPositionY
        cellPositionZ[n]=z; // asign z for cellPositionZ
    }
}