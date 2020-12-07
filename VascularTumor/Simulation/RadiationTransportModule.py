
from PySteppables import *
import CompuCell
import sys
import os
from scipy import interpolate
from subprocess import call
import subprocess
from collections import OrderedDict
import csv
import radcell
import math
import random
import numpy as np
import radcell

class WriteCellInformation():
    def __init__(self, dim,cellDict,_DCF):
        self.dim = []
        self.cellDict = OrderedDict()
        self.dim = dim
        self.cellDict = cellDict
        self.DCF = _DCF
        
    def Write(self, path):
        ofile = open (path,"w")
        TDF = 1000.0 # tissue dimension factor, tissue dimension unit in Geant4 is mm, in CC3D is um
        dX = self.dim[0][0]*self.DCF/TDF # dimension in X
        dY = self.dim[0][1]*self.DCF/TDF # dimension in Y
        dZ = self.dim[0][2]*self.DCF/TDF # dimension in Z
 
        writer = csv.writer(ofile)
        writer.writerow([dX,dY,dZ])
        
        for key, value in self.cellDict.items():
            cell_x=value[0]*self.DCF/TDF-dX/2.0  # this is changing coordinate system between CC3D and Geant4
            cell_y=value[1]*self.DCF/TDF-dY/2.0
            cell_z=value[2]*self.DCF/TDF-dZ/2.0
            cellType = value[3]
            writer.writerow([key,cell_x,cell_y,cell_z,cellType])
        ofile.close()
      

class RadiationTransport(SteppableBasePy):
    def __init__(self,dim,cell_List,_DCF):
        self.cellList = cell_List
        self.ddim = dim
        self.DCF = _DCF  # dimension conversion factor, convert pixel in CC3D to phyical length in Geant4
    def Run(self,runMode,radiationSource):
        cellInfoDict=OrderedDict()
        cellInfoList=[]  
        dim=[]
        dim.append((self.ddim.x,self.ddim.y,self.ddim.z))
        
        for cell in self.cellList:
            cellInfoList=[cell.xCOM,cell.yCOM,cell.zCOM,cell.type]
            cellInfoDict[str(cell.id)]=cellInfoList
        
        writeCellInfo = WriteCellInformation(dim,cellInfoDict,self.DCF) 
        cellInformationPath = os.path.dirname(__file__) +"/cellInformation.csv"
        writeCellInfo.Write(cellInformationPath)
        child = subprocess.Popen(["python","RADCellSimulation.py",runMode,radiationSource],cwd=os.path.dirname(__file__))#check the geometry and radiation source
        child.wait()

class RadiationSource(SteppableBasePy):
    def __init__(self,sourceName):
        self.sourceName = sourceName
    def GetSourceID(self):
        return self.sourceName

    def PointSource(self,cX,cY,cZ,radiation,energy,beamNum):
        ofile = open (os.path.dirname(__file__)+"/"+self.sourceName+ ".in","w") # creating the radiation source file
        ofile.write('/control/verbose 0 \n')
        ofile.write('/tracking/verbose 0 \n')
        ofile.write('/gps/verbose 0 \n')
        ofile.write('/gps/particle'+ ' '+radiation+'\n') # could be e-,proton, alpha, etc...
        ofile.write('/gps/pos/type Point \n')
        ofile.write('/gps/pos/centre'+' '+ str(cX)+' '+str(cY)+ ' '+ str(cZ)+ ' mm'+'\n')
        ofile.write('/gps/ang/rot2 0 0 -1 \n')
        ofile.write('/gps/ene/mono ' + str(energy)+' '+'keV'+ '\n')
        ofile.write('/run/beamOn' + ' ' +str(beamNum)+ ' '+ '\n')
        ofile.close() 
        
    def PlaneSource(self,cX,cY,cZ,pX,pY,radiation,energy,beamNum):
        ofile = open (os.path.dirname(__file__)+"/"+self.sourceName+ ".in","w") # creating the radiation source file
        ofile.write('/control/verbose 0 \n')
        ofile.write('/tracking/verbose 0 \n')
        ofile.write('/gps/verbose 0 \n')
        ofile.write('/gps/particle'+ ' '+ radiation+'\n') # could be e-,proton, alpha, etc...
        ofile.write('/gps/pos/type Plane \n')
        ofile.write('/gps/pos/shape Square \n')
        ofile.write('/gps/pos/centre'+' '+ str(cX)+' '+str(cY)+ ' '+ str(cZ)+ ' mm'+'\n')
        ofile.write('/gps/pos/halfx '+str(pX)+' mm'+'\n')
        ofile.write('/gps/pos/halfy '+str(pY)+' mm'+'\n')
        ofile.write('/gps/ang/rot2 0 0 -1 \n')
        ofile.write('/gps/ene/mono ' + str(energy)+' '+'keV'+ '\n')
        ofile.write('/run/beamOn' + ' ' +str(beamNum)+ ' '+ '\n')
        ofile.close() 
        
    def PlaneSource_MRT(self,cX,cY,cZ,pX,pY,radiation,energy,beamNum):
        ofile = open (os.path.dirname(__file__)+"/"+self.sourceName+ ".in","w") # creating the radiation source file
        ofile.write('/control/verbose 0 \n')
        ofile.write('/tracking/verbose 0 \n')
        ofile.write('/gps/verbose 0 \n')
        ofile.write('/gps/source/intensity 1\n')
        ofile.write('/gps/particle'+ ' '+ radiation+'\n') # could be e-,proton, alpha, etc...
        ofile.write('/gps/pos/type Plane \n')
        ofile.write('/gps/pos/shape Square \n')
        ofile.write('/gps/pos/centre'+' '+ str(cX)+' '+str(cY)+ ' '+ str(cZ)+ ' mm'+'\n')
        ofile.write('/gps/pos/halfx '+str(pX)+' mm'+'\n')
        ofile.write('/gps/pos/halfy '+str(pY)+' mm'+'\n')
        ofile.write('/gps/ang/rot2 0 0 -1 \n')
        ofile.write('/gps/ene/mono ' + str(energy)+' '+'keV'+ '\n')
        ofile.write('/gps/source/add 1\n') # 1 added source
        ofile.write('/gps/particle'+ ' '+ radiation+'\n') # could be e-,proton, alpha, etc...
        ofile.write('/gps/pos/type Plane \n')
        ofile.write('/gps/pos/shape Square \n')
        ofile.write('/gps/pos/centre'+' '+ str(cX)+' '+str(cY+0.1)+ ' '+ str(cZ)+ ' mm'+'\n')
        ofile.write('/gps/pos/halfx '+str(pX)+' mm'+'\n')
        ofile.write('/gps/pos/halfy '+str(pY) +' mm'+'\n')
        ofile.write('/gps/ang/rot2 0 0 -1 \n')
        ofile.write('/gps/ene/mono ' + str(energy)+' '+'keV'+ '\n')
        ofile.write('/gps/source/add 1\n') # 2 added source
        ofile.write('/gps/particle'+ ' '+ radiation+'\n') # could be e-,proton, alpha, etc...
        ofile.write('/gps/pos/type Plane \n')
        ofile.write('/gps/pos/shape Square \n')
        ofile.write('/gps/pos/centre'+' '+ str(cX)+' '+str(cY-0.1)+ ' '+ str(cZ)+ ' mm'+'\n')
        ofile.write('/gps/pos/halfx '+str(pX)+' mm'+'\n')
        ofile.write('/gps/pos/halfy '+str(pY) +' mm'+'\n')
        ofile.write('/gps/ang/rot2 0 0 -1 \n')
        ofile.write('/gps/ene/mono ' + str(energy)+' '+'keV'+ '\n')
        ofile.write('/gps/source/add 1\n') # 3 added source
        ofile.write('/gps/particle'+ ' '+ radiation+'\n') # could be e-,proton, alpha, etc...
        ofile.write('/gps/pos/type Plane \n')
        ofile.write('/gps/pos/shape Square \n')
        ofile.write('/gps/pos/centre'+' '+ str(cX)+' '+str(cY-0.2)+ ' '+ str(cZ)+ ' mm'+'\n')
        ofile.write('/gps/pos/halfx '+str(pX)+' mm'+'\n')
        ofile.write('/gps/pos/halfy '+str(pY) +' mm'+'\n')
        ofile.write('/gps/ang/rot2 0 0 -1 \n')
        ofile.write('/gps/ene/mono ' + str(energy)+' '+'keV'+ '\n')
        ofile.write('/gps/source/add 1\n') # 4 added source
        ofile.write('/gps/particle'+ ' '+ radiation+'\n') # could be e-,proton, alpha, etc...
        ofile.write('/gps/pos/type Plane \n')
        ofile.write('/gps/pos/shape Square \n')
        ofile.write('/gps/pos/centre'+' '+ str(cX)+' '+str(cY+0.2)+ ' '+ str(cZ)+ ' mm'+'\n')
        ofile.write('/gps/pos/halfx '+str(pX)+' mm'+'\n')
        ofile.write('/gps/pos/halfy '+str(pY) +' mm'+'\n')
        ofile.write('/gps/ang/rot2 0 0 -1 \n')
        ofile.write('/gps/ene/mono ' + str(energy)+' '+'keV'+ '\n')
        
        ofile.write('/run/beamOn' + ' ' +str(beamNum)+ ' '+ '\n')
        ofile.close() 
        


class CellSimulationResult():
        def __init__(self,_fileName):
            self.fileName = _fileName
            self.singleCellFileName = 'Mono_Electron_Eng_1PNum_100000'
            self.readResult = radcell.ReadRadiationTransportInfo()
 
        def ReadTransportTally(self,_prescribed_dose):
            prescribed_dose = _prescribed_dose
            self.readResult.ReadDoseTallyOutPut(os.path.dirname(__file__)+"/"+str(self.fileName)+"_dose.csv") # read normalized cell dose file
            self.readResult.ReadDNADamageTallyOutPut(os.path.dirname(__file__)+"/"+str(self.fileName)+"_DNADamage.csv") # read normalized cell DSB file
            singleCellDoseFileName = self.singleCellFileName+"_dose.csv"
            singleCellDSBFileName = self.singleCellFileName+"_DNADamage.csv"
            print 'single cell file name is ', singleCellDoseFileName, singleCellDSBFileName 
            self.readResult.ReadSingleCellDoseAsReference(os.path.dirname(__file__)+"/"+singleCellDoseFileName)
            self.readResult.ReadSingleCellDNADamageAsReference(os.path.dirname(__file__)+"/"+singleCellDSBFileName)
            self.readResult.GetAbsoluteResultsByCOIMethod(prescribed_dose, "indirect","mean")
        def GetCellDose(self,cellID):
            electron_energy = 1
            particle_flux = 3.38E9 # this is the particle flux of 1MeV electron in water for calculating the absolute dose
            absolute_dose = self.readResult.GetAbsDoseOfCell(cellID)
            absolute_DSB = self.readResult.GetAbsDSBOfCell(cellID)
            cellDose = self.readResult.GetMCDoseOfCell(cellID)
            cellDSB = self.readResult.GetMCDSBOfCell(cellID)
            
            return absolute_dose, absolute_DSB, cellDose, cellDSB
            

        
      