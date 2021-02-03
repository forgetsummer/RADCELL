# RADCELL
RADCELL is a radiation transport simulation module developed for conducting radiation transport simulation in cells. It can simulate the cell dose and cell DNA damages, such as single-strand breaks (SSBs) and double-strand breaks (DSBs). The RADCELL is developed based on the microdosimetry example of Geant4.  

The primary function of RADCELL is to calculate the radiation dose to cell organelles, and DNA damages to cells. We propose a three-dimensional (3D) cellular compartment model, which incorporates two cellular compartments including nucleus and cytoplasm. We use the sphere to approximate the cell shape which substantially simplifies the complexity of cell geometry but there is no significant big accuracy penalty. The nucleus is modeled as a sphere which is located at the center of the cell. The size of the cell and the nucleus can be customized according to the biological cell which will be studied during the simulation. 
During the radiation transport simulation, the energy deposition information in each cell will be collected, and the information, such as the index for indicating the affected cells, cellular dose, and affected cell organelles, etc., could be used to quantify the cell dose and DNA damages. 

**RADCellSimulation**

The RADCellSimulation is the class used to control the whole process of radiation transport simulation of cell. It takes care of building the geometry of cell culture, cell visualization, radiation track structure visualization, cell geometry update, etc. The CreateCell() method is used to create a cell line. The SetCellWorld() method is used to build cell culture geometry. These methods are used in the initial phase when we describe cell and tissue geometry in the Geant4 simulation. The UpdateGeometryInitialize() method and UpdateGeometryFinalize() methods are used to update the cell geometry since the cell geometry may change due to cell proliferation or cell death. The EnergyDistributionCalculation() method is used to invoke the Geant4 transport kernel to implement the radiation transport simulation of cells.

**Cell**

The Cell class is used to describe one type of cell line. Some attributes are used to describe the characteristics of cell. For instance, cellType indicates the type of cell, cellOrganelle is a vector that stores the cell organelles, and cellShape indicates the shape of cell. The CellConstruct() method is used to construct a cell line. This class is not to be confused with CompuCell3D CellG class.  RADCELL’s Cell class is a representation of the CellG class in the RADCELL module, and thus every CellG CompuCell3D object will have matching RADCELL Cell object. 
 
 **DetectorConstruction**
 
The DetectorConstruction class is used to build simulation geometry. It serves the function of translating the geometry of cell culture to the corresponding geometry in radiation transport simulation. The SetCellWorld() method is used to define the boundaries of simulation geometry. In Geant4, a world is used to define the boundaries of geometry. Usually, we can use a box to define the world where the whole cell culture locates inside. The cells in the cell culture may belong to different cell lines, so before seeding the cell into cell culture, the cell objects of cell lines should be created using the method CreateCell(). The information of cell lines will be stored in a vector generalCellMap. The SetCellPlacement() method is used to seed the cell in cell culture. This method will create an actual physical volume to represent the cells seeded in the cell culture. The SetCellSimulationParameter() method is used to specify the simulation parameters for simulation, such as setting the production cuts of secondary particles and specifying the target logical volume which will apply Geant4-DNA crosssection data. The GetEdepCellInformation()  method is used to determine which cell and to which cell organelle the energy deposition point belongs.

**DBSCAN**

The DBSCAN is a class used to implement the DBSCAN clustering algorithm [4]. It serves as a clustering kernel for clustering two-dimensional and three-dimensional data. In RADCellSimulation, the energy deposition points will be clustered using the DBSCAN algorithm.

**EventAction**

In Geant4, there are two big units of simulation, i.e., Run and Event. An event is the basic unit of simulation in Geant4. The event corresponds to the interactions of one particle from the radiation source with matter (in our case cell). One event contains the whole transport process of one primary particle. The run is the largest unit of simulation. One run consists of a sequence of events. Thus a run corresponds to simulating a time interval in which many events of radiation-matter interactions take place [5]. To control or to keep track of what happens during the event, Geant4 introduces the EventAction class (class directly derived from the base class G4UserEventAction in Geant4 toolkit). The EventAction is the class for defining actions of event. This class has two virtual methods which are invoked by G4EventManager for each event, and they are BeginOfEventAction() and EndOfEventAction(). The EndOfEventAction() method is invoked at the very end of event processing. It is typically used for a simple analysis of the processed event. For scoring simulation results, tally is the process of scoring the parameters of interest; for example, keeping track of energy depositions. In radiation transport simulation, the tally of cell dose and DNA damage for each event will be invoked in this method. In RADCELL, we use CellDoseTally() method to get the cell dose tally of the event and  the CellDNADamageTally() to get the cell DNA damage tally of the event.

**RunAction**

To keep track of what happens during multiple Geant4 events i.e., during multiple events or radiation-matter interaction events, RADCELL uses RunAction class which is directly derived from the base class G4UserRunAction in Geant4 toolkit. This class has two key methods, BeginOfRunAction() and EndOfRunAction(), defined in the base class G4UserRunAction. The EndOfRunAction() is invoked at the very end of a run, and we use it for storing/printing histogram or manipulating run summary. For example, we implemented final tally of cell dose and DNA damage in this method. For final tallies, Geant4 uses G4RunDoseAnalysis class to do the final tally of cell dose. In G4RunDoseAnalysis, the CollectEventDoseTallyMap()method will collect the dose tally results from all the events, then the CellDoseTally() method will do the final tally of cell dose and statistical error analysis. Similarly, G4RunDNADamageAnalysis class is created to do the final tally of DNA damage and statistical error analysis. 
In G4RunDNADamageAnalysis, the CollectEventDNADamageTallyMap() method will collect the DNA damage tally results from all the events, then the DNADamageTally() method will do the final tally of cell DNA damage. When the run is finished, one object of G4RunDoseAnalysis will be instantiated inside the EndOfRunAction() method to get the final dose tally results of a run. Similarly, one object of G4RunDNADamageAnalysis will be instantiated inside the EndOfRunAction() method to get the final DNA damage tally of run. The calculation process of dose tally and DNA damage tally is shown in Figure 3.
 

**PhysicsList**

The way Geant4 simulates interaction of radiation with matter is by keeping track of various physical mechanisms in which the radiation will occur. Specifying which interactions are relevant for our simulation takes place in the PhysicsList class which is derived, in large part, from the microdosimetry example included in the Geant4 toolkit. This class activates a full suite of Geant4-DNA specific physics included in Geant4. Physics processes and models in Geant4 are enabled by region which is known in program as G4Region. In G4Region, Geant4-DNA physics are enabled.

**PrimaryGeneratorAction**

Once we defined cells in the Geant4 simulation along with all the infrastructure necessary to keep track of radiation-induced effects, we need to specify radiation source. We use so called general physical source (GPS) within radiation transport solver to describe properties of the radiation source. The GPS allows users to define the source according to the efficient built-in accommodations for source location distributions, energy distributions, angular distributions, and other features.

**SteppingAction**

PhysicsList defines many types of interactions that may happen when a single event of radiation-matter interaction happens. Geant4 will simulate those types of interaction sequentially, and it uses SteppingAction class to allow users to keep track of what happens when a single type of interaction takes place. The SteppingAction is derived, in large part, from the microdosimetry example included in the Geant4 toolkit. SteppingAction inherits the UserSteppingAction() method from the base class G4UserSteppingAction. When the UserSteppingAction() is invoked, we keep track of the information of energy deposition and store it in the EventAction object.

 

**RADCELL Implementation**

We first instantiate RADCellSimulation object and this object’s methods to build the simulation geometry and run the transport simulation. The RADCellSimulation supports two run modes: gui mode and detailed mode. When gui mode is chosen, RADCellSimulation will only run radiation transport, and the information of energy deposition points will not be collected and processed during the simulation. Typically, we use this mode to check whether the cell geometry is right or not. When a detailed mode is chosen, the RADCellSimulation will collect and process the information of energy deposition points during the simulation. The detailed mode should be chosen if we want to get the cell dose and DNA damage tallies.
 
For radiation transport simulation of multicellular system, firstly, the computerized cells are created in the model, then conducting radiation transport calculation to obtain the simulation results, such as cellular dose and double-strand breaks. It is worth noting that the geometry in Geant4 should be updated if the multicellular system changed due to mitosis or cell death. 


**About installation and using RADCELL**
Currently, RADCELL only supports Linux system. 
1. Install Geant4 in your local computer
2. Compile RADCellSimulation using Geant4
3. After installing RADCellSimulation, add the PYTHONPATH and LD_LIBRARY_PATH for RADCELL in bashrc file.
   export PYTHONPATH=<path of RADCELL installation>/bin:$PYTHONPATH 
   export LD_LIBRARY_PATH=<path of RADCELL installation>/bin:$LD_LIBRARY_PATH 
 4. After setting up the environemtnal variable as step 3, run demo_usingRADCELL.py, to see whether RADCELL path is set up correctly
 5. Install CompuCell3D in your local computer (https://compucell3d.org/)
 6. For runnig the example of simulating vascular tumor reponse,  in CompuCell3D, open project VascularTumor. Simply open file VascularTumor.cc3d. Then run the model.
 
 
 **Citing RADCELL**
 
RADCELL is described in a technical paper, available in https://iopscience.iop.org/article/10.1088/1361-6560/abd4f9. If you would like the cite the article, you may use the following citation:

Liu R, Higley KA, Swat MH, Chaplain M, Powathil G, Glazier JA. Development of a coupled simulation toolkit for computational radiation biology based on Geant4 and CompuCell3D. Phys Med Biol. 2020 Dec 18. doi: 10.1088/1361-6560/abd4f9. Epub ahead of print. PMID: 33339019.


