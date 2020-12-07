from PySteppables import *
from PySteppablesExamples import MitosisSteppableBase
import CompuCell
import sys
from random import uniform
import math
import random
import radcell
import RadiationTransportModule as RTM
import scipy.stats as stats
import os

class VolumeParamSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1,):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.fieldNameVEGF2 = 'VEGF2'
        self.fieldNameGlucose = 'Glucose'
        self.numberOfTumorCellList = []  
        self.numberOfVascularCellList = []
        self.numberOfNeovascularCellList = []
        self.numberOfNecroticCellList = []
        self.pDose = 8
        self.simulationID = 'single_MRT_'+str(self.pDose)+'Gy_5_hyperfraction'
#         self.simulationID = 'single_MRT_20Gy_2fractions'
        self.cellNumberFileOut = open(os.path.dirname(__file__)+"/"+ self.simulationID+'_cellNumber.csv','w')
    def start(self):
        self.pW=self.addNewPlotWindow(_title='Cell Number VS Time',_xAxisTitle='MonteCarlo Step (MCS)',_yAxisTitle='Cell Number', _xScaleType='linear',_yScaleType='linear')
        self.pW.addPlot('Tumor Cells',_style='Dots',_color='green',_size=5)
        self.pW.addPlot('Vascular Cells',_style='Dots',_color='red',_size=5)
        self.pW.addPlot('Neovascular Cells',_style='Dots',_color='yellow',_size=5)
        self.pW.addPlot('Necrotic Cells',_style='Dots',_color='blue',_size=5)

        for cell in self.cellList:
            if cell.type==self.VASCULAR or cell.type==self.NEOVASCULAR:
            #due to pressue from chemotaxis to vegf1, cell.volume is smaller that cell.target volume
            #in this simulation the offset is about 10 voxels.
                cell.targetVolume=64.0+10.0
                cell.lambdaVolume=20.0
            else:
                cell.targetVolume=32.0
                cell.lambdaVolume=20.0
                
        DCF = 10
        dX = self.dim.x*DCF/float(1000) #change unit into mm
        dY = self.dim.y*DCF/float(1000)
        dZ = self.dim.z*DCF/float(1000)
       

        radSource = RTM.RadiationSource('testInputSource')
#         radSource.PlaneSource_MRT(0,0,dZ/2.0,dX/2.0,dY/20.0,'e-',1000,0)
        radSource.PlaneSource(0,0,dZ/2.0,dX/2.0,dY/20.0,'e-',1000,1000)
        radTransport = RTM.RadiationTransport(self.dim,self.cellList,DCF)# importing tissue and cell geometry from CC3D
#         radTransport.Run("gui",radSource.GetSourceID()) # check geometry and radiation source 
        


    def step(self,mcs):
        fieldVEGF2=CompuCell.getConcentrationField(self.simulator,self.fieldNameVEGF2)
        fieldGlucose=CompuCell.getConcentrationField(self.simulator,self.fieldNameGlucose)
        
#         if mcs == 10:
#         if mcs == 5000:
#         if mcs == 7000:
#         if mcs == 9000:
#         if mcs == 11000:
#         if mcs == 14000:
#         if mcs == 1000 or mcs == 5000 or mcs == 7000 or mcs == 9000 or mcs == 11000 or mcs==14000:
#         if mcs == 12000 or mcs == 16000: # hyper-fraction scheme
        if mcs == 12000 or mcs==13000 or mcs == 14000 or mcs == 15000 or mcs == 16000:
            DCF = 10
            dX = self.dim.x*DCF/float(1000) #change unit into mm
            dY = self.dim.y*DCF/float(1000)
            dZ = self.dim.z*DCF/float(1000)
           

#             radSource = RTM.RadiationSource('testInputSource')
#             radSource.PlaneSource(0,0,dZ/4.0,dX/100.0,dY/100.0,'e-',1000,0)

#             radTransport = RTM.RadiationTransport(self.dim,self.cellList,DCF)# importing tissue and cell geometry from CC3D
#             radTransport.Run("gui",radSource.GetSourceID()) # check geometry and radiation source 
            
            radSource = RTM.RadiationSource('testInputSource')
            radSource.PlaneSource_MRT(0,0,dZ/2.0,dX/2.0,dY/20.0,'e-',1000,1000)
#             radSource.PlaneSource(0,0,dZ/2.0,dX/2.0,dY/20.0,'e-',1000,500)
            radTransport = RTM.RadiationTransport(self.dim,self.cellList,DCF)# importing tissue and cell geometry from CC3D
            runMode = 'gui'+' ' + self.simulationID
            print 'the runMode is: ', runMode
#             radTransport.Run("gui",radSource.GetSourceID()) # check geometry and radiation source 
#             radSource.PlaneSource(0,0,dZ/2.0,dX/2.0,dY/20.0,'e-',1000,5000)
            radTransport.Run(runMode,radSource.GetSourceID()) # check geometry and radiation source
            simResult = RTM.CellSimulationResult(self.simulationID)
            simResult.ReadTransportTally(self.pDose)

            for cell in self.cellList:

                alpha = 0.25 # alpha value for quantifying effective energy
                beta  = 0.03 # beta value for quantifying effective energy
                E1 = 0 # mean state energy of S1 state
                E3= 48.47 # mean state energy of S2 state
                sigma = 6.96 # sigma of state energy distribution
                
                cellDose, cellDSB, normalized_dose,normalized_DSB= simResult.GetCellDose(cell.id)
#                 print 'The dose for cell: ',cell.id, ' is: ',cellDose
#                 print 'The DSB for cell: ', cell.id, 'is: ' , cellDSB
#                 print 'The normalized dose for cell: ', cell.id, 'is: ', normalized_dose
#                 print 'The normalized DSB for cell: ', cell.id, 'is: ', normalized_DSB
                
                delta_E = alpha*cellDSB # here we just consider the direct effect
                x = -abs(E1+delta_E-E3)/2/sigma
                print 'delta energy is ',x
                p_total = 2*stats.norm.cdf(x,0,1)
                if E1+delta_E-E3>=0:
                    p_total = 1# boundary condition
                print 'the total jumping probability is: ',p_total
                jumpingProbabilityPerMCS = p_total 
                if cellDose!=0:
                    if random.random()<jumpingProbabilityPerMCS:
                        cell.type = self.NECROTIC # if so, cell dies due to radiation killing
        
        numberOfTumorCells = 0
        numberOfVascularCells = 0
        numberOfNeovascularCells = 0
        numberOfNecroticCells = 0
        for cell in self.cellList:
            #print cell.volume
            #NeoVascular
            if cell.type == self.VASCULAR:
                numberOfVascularCells  = numberOfVascularCells  + 1
            if cell.type == self.NEOVASCULAR:
                numberOfNeovascularCells = numberOfNeovascularCells + 1
                totalArea = 0
                # pt=CompuCell.Point3D()
                # pt.x=int(round(cell.xCM/max(float(cell.volume),0.001)))
                # pt.y=int(round(cell.yCM/max(float(cell.volume),0.001)))
                # pt.z=int(round(cell.zCM/max(float(cell.volume),0.001)))
                
                # VEGFConcentration=fieldVEGF2.get(pt)
                
                VEGFConcentration=fieldVEGF2[int(round(cell.xCOM)),int(round(cell.yCOM)),int(round(cell.zCOM))]
                
                # cellNeighborList=CellNeighborListAuto(self.nTrackerPlugin,cell)
                cellNeighborList=self.getCellNeighbors(cell)
                for neighborSurfaceData in cellNeighborList:
                    #Check to ensure cell neighbor is not medium
                    if neighborSurfaceData.neighborAddress:
                        if neighborSurfaceData.neighborAddress.type == self.VASCULAR or neighborSurfaceData.neighborAddress.type == self.NEOVASCULAR:                            
                            #sum up common surface area of cell with its neighbors
                            totalArea+=neighborSurfaceData.commonSurfaceArea 
                            #print "  commonSurfaceArea:",neighborSurfaceData.commonSurfaceArea
                #print totalArea        
                if totalArea < 45:
                    #Growth rate equation
                    
                    cell.targetVolume+=2.0*VEGFConcentration/(0.01 + VEGFConcentration)
                    print "totalArea", totalArea,"cell growth rat& voxel/MCS \\e: ", 2.0*VEGFConcentration/(0.01 + VEGFConcentration),"cell Volume: ", cell.volume
         
            #Proliferating Cells
            if cell.type == self.PROLIFERATING:
                
                numberOfTumorCells = numberOfTumorCells + 1 # tallying the tumor cells
                # pt=CompuCell.Point3D()
                # pt.x=int(round(cell.xCM/max(float(cell.volume),0.001)))
                # pt.y=int(round(cell.yCM/max(float(cell.volume),0.001)))
                # pt.z=int(round(cell.zCM/max(float(cell.volume),0.001)))
                # GlucoseConcentration=fieldGlucose.get(pt)
                
                GlucoseConcentration=fieldGlucose[int(round(cell.xCOM)),int(round(cell.yCOM)),int(round(cell.zCOM))]
                
                # Proliferating Cells become Necrotic when GlucoseConcentration is low
                if  GlucoseConcentration < 0.001 and mcs>1000:
                    cell.type = self.NECROTIC
                    #set growth rate equation -- fastest cell cycle is 24hours or 1440 mcs--- 32voxels/1440mcs= 0.022 voxel/mcs
                cell.targetVolume+=0.022*GlucoseConcentration/(0.05 + GlucoseConcentration)
                #print "growth rate: ", 0.044*GlucoseConcentration/(0.05 + GlucoseConcentration), "GlucoseConcentration", GlucoseConcentration

            #Necrotic Cells
            if cell.type == self.NECROTIC:
                #sNecrotic Cells shrink at a constant rate
                numberOfNecroticCells = numberOfNecroticCells + 1
                cell.targetVolume-=0.1
        
        self.numberOfTumorCellList.append(numberOfTumorCells)
        self.numberOfVascularCellList.append(numberOfVascularCells)
        self.numberOfNeovascularCellList.append(numberOfNeovascularCells)
        self.numberOfNecroticCellList.append(numberOfNecroticCells)
        
        self.pW.addDataPoint("Tumor Cells", mcs,numberOfTumorCells)
        self.pW.addDataPoint("Vascular Cells", mcs,numberOfVascularCells)
        self.pW.addDataPoint("Neovascular Cells",mcs,numberOfNeovascularCells)
        self.pW.addDataPoint("Necrotic Cells",mcs,numberOfNecroticCells)
 
    
    def finish(self):
        print 'the number of tumor  cells',self.numberOfTumorCellList
        for index,i in enumerate( self.numberOfTumorCellList):
            vascularNum = self.numberOfVascularCellList[index]
            neovascularNum = self.numberOfNeovascularCellList[index]
            necroticNum = self.numberOfNecroticCellList[index]
            self.cellNumberFileOut.write('{0:d},{1:d},{2:d},{3:d},{4:d}\n'.format(index,i,vascularNum,neovascularNum,necroticNum)) # save the tumor cell number for each MCS
        self.cellNumberFileOut.close()

        
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
     
    def step(self,mcs):
        
        cells_to_divide=[]
          
        for cell in self.cellList:
            if cell.type == self.PROLIFERATING and cell.volume>64:
                cells_to_divide.append(cell)
            if cell.type== self.NEOVASCULAR and cell.volume>128:
                cells_to_divide.append(cell)

                     
        for cell in cells_to_divide:

            self.divideCellRandomOrientation(cell)
            
    def updateAttributes(self):    
        self.parentCell.targetVolume /= 2.0 # reducing parent target volume                 
        self.cloneParent2Child()            
              
                