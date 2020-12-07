import sys
import radcell
import os
import itertools

from collections import OrderedDict

print 'Start calling subprocess RADCellSimulation here'
print 'The directory now is: ', os.getcwd()
cellInformationPath = os.getcwd() +"/cellInformation.csv"
print 'the cellInformationPath now is ', cellInformationPath
print 'check files is: ',os.path.isfile(cellInformationPath)

try:
    ifile = open(cellInformationPath,"r")
    line=ifile.readlines()
    cellInfoList=[]
    cellInfoDict=OrderedDict()
    for i in range(len(line)):
        a=line[i].split(",")
        if i==0:
            dX=float(a[0])
            dY=float(a[1])
            dZ=float(a[2])
        else:
            cellInfoList=[float(a[1]),float(a[2]),float(a[3]),int(a[4])]
            cellInfoDict[a[0]]=cellInfoList

    ifile.close()
except IOError:
    fp = open(cellInformationPath,"w")


sim=radcell.RADCellSimulation()

sim.SetCellWorld(dX,dY,dZ) # change unit here in mm
print 'the dimension of tissue is ',dX, dY, dZ

# cell color and nucleus color could be as :
# Gray, Grey, Black, Green, Red, Blue, Magenta, Yellow

cellType_list=[]
cellType_list.append("Epithelia0")
cellType_list.append("Epithelia1")
cellType_list.append("Epithelia2")
cellType_list.append("Epithelia3")


sim.CreateCell(cellType_list[0],"Cytoplasma Nucleus","Sphere","Green Green")
sim.CreateCell(cellType_list[1],"Cytoplasma Nucleus","Sphere","Blue Green")
sim.CreateCell(cellType_list[2],"Cytoplasma Nucleus","Sphere","Red Green")
sim.CreateCell(cellType_list[3],"Cytoplasma Nucleus","Sphere","Yellow Green")
sim.SetCellSimulationParameter(1,"Nucleus")


for key, value in cellInfoDict.items():
    sim.CellConstruction(int(key),cellType_list[value[3]-1],value[0],value[1],value[2],10,5)


# #     print 'the location of cell',int(key), 'is ',cell_x,cell_y,cell_z
# #     sim.CellConstruction(int(key),"Epithelia0",cell_x,cell_y,cell_z,2,1)
    
#     if int(key)==1 :
#         print 'the location of cell',int(key), 'is ',cell_x,cell_y,cell_z
#         sim.CellConstruction(int(key),"Epithelia0",cell_x,cell_y,cell_z,2,1)
#     if int(key)==2:
#         print 'the location of cell',int(key), 'is ',cell_x,cell_y,cell_z
#         sim.CellConstruction(int(key),"Epithelia1",cell_x,cell_y,cell_z,2,1))
#     if int(key)==3:
#         print 'the location of cell',int(key), 'is ',cell_x,cell_y,cell_z
#         sim.CellConstruction(int(key),"Epithelia2",cell_x,cell_y,cell_z,2,1)
    
    
print 'Before starting calling radiation energy deposition caclculation function \n'


sim.RADCellSimulationInitializePyWrapper(len(sys.argv), sys.argv)

print 'the argument is ', sys.argv

runMode=sys.argv[1]
radiationSource = sys.argv[2]+".in"

sim.EnergyDistributionCalculation(runMode,radiationSource,False)

