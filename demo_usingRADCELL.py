import radcell
import sys

argc = len(sys.argv)
argv = sys.argv
argv_vec = radcell.vectorstring()
for i in range(argc):
    argv_vec.push_back(argv[i])
    
mytest=radcell.testCellPosition()
N=300
mytest.setPosition(N)
mySim=radcell.RADCellSimulation()
mySim.SetCellWorld(1.1,1.1,1.1)
mySim.SetCellSimulationParameter(1,"Nucleus")
mySim.CreateCell("Epithelia","Cytoplasma Nucleus","Sphere","Red Green")
for i in range (N):
    mySim.CellConstruction(i,"Epithelia",mytest.passPositionX(i),mytest.passPositionY(i),mytest.passPositionZ(i),20,5)
mySim.EnergyDistributionCalculationPyWrapper(argc,argv_vec)
