
// This program is largely amended from the microdosimetry simulation from Geant4-DNA project
// This is program is written by Ruirui Liu, at Department of Nuclear Engineering and Radiation Health physics,
// Oregon State University
// December 10,2015 
// reference paper is :
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// $ID$
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "G4Run.hh"
#include "TrackingAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "Analysis.hh"
#include "G4Threading.hh"
#include <fstream>
#include "G4CsvAnalysisManager.hh"
#include <sstream>
#include "ArgumentInterpreter.hh"
#include "CellDNADamageAnalysis.hh"
#include "CellDoseAnalysis.hh"
#include <math.h>
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunDoseAnalysis.hh"
#include "G4RunDNADamageAnalysis.hh"

using namespace std;

void PrintNParticles(std::map<const G4ParticleDefinition*, int>& container);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction() : G4UserRunAction(),
      fpTrackingAction(0), fInitialized(0), fDebug(false)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void RunAction::BeginOfRunAction(const G4Run* run)
{
  // In this example, we considered that the same class was
  // used for both master and worker threads.
  // However, in case the run action is long,
  // for better code review, this practice is not recommanded.
  //
  // Please note, in the example provided with the Geant4 X beta version,
  // this RunAction class were not used by the master thread.

  totalNumberOfEvent=run->GetNumberOfEventToBeProcessed();// get the total number of event to be processed, which is same as the N in /beamOn N
  
  
  bool sequential = (G4RunManager::GetRunManager()->GetRunManagerType() == 
                     G4RunManager::sequentialRM);

 
  
  if(isMaster && sequential == false ) // if it is false, then it is the master mode, so use BeginMaster(run)
  // WARNING : in sequential mode, isMaster == true
  {
    BeginMaster(run);
  }
  else BeginWorker(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{

  bool sequential = (G4RunManager::GetRunManager()->GetRunManagerType() == 
                     G4RunManager::sequentialRM);

  if(isMaster && sequential == false)
  {
    EndMaster(run);
  }
  else
  {
    EndWorker(run);
  }
   
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginMaster(const G4Run* run)
{
  if(fDebug)
  {
    bool sequential = (G4RunManager::GetRunManager()->GetRunManagerType() == 
                       G4RunManager::sequentialRM);
    G4cout << "===================================" << G4endl;
    if(!sequential)
      G4cout << "================ RunAction::BeginMaster" << G4endl;
    PrintRunInfo(run);
    G4cout << "===================================" << G4endl;
  }
  
  G4String simulationMode = ArgumentInterpreter::GetSimulationMode();
  if (simulationMode=="out")
  {
      bool postProcessMarker=ArgumentInterpreter::GetPostProcessMarker();
      if (postProcessMarker)
      {
        G4String fileName = ArgumentInterpreter::GetOutPutFileName();
        G4String outputEdepFileName=fileName+"_edep"+".csv";// each worker will have each own writting file
        ofstream file;
        std::remove(outputEdepFileName); // first, remove the old file if the file exists
        file.open(outputEdepFileName); // then open the new file for writing
        file<<"TotalEventNumber"<<","<<totalNumberOfEvent<<endl; // write totalNumberOfEvent to the file in first line
        file<<"Event ID"<<","<<"Cell ID"<<","<<"X/nm"<<","<<"Y/nm"<<","<<"Z/nm"<<","<<"Edep/eV"<<","<<"Cell Organelle"<<endl;// write titles
        file.close();

       }

  }

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginWorker(const G4Run* run)
{
  if (fDebug)
  {
    G4cout << "===================================" << G4endl;
    G4cout << "================ RunAction::BeginWorker" << G4endl;
    PrintRunInfo(run);
    G4cout << "===================================" << G4endl;
  }
  if(fInitialized == false) InitializeWorker(run);

  
    G4String simulationMode = ArgumentInterpreter::GetSimulationMode();
  if (simulationMode=="out")
  {
      bool postProcessMarker=ArgumentInterpreter::GetPostProcessMarker();
      if (postProcessMarker)
      {
        CreateHistogram();  // creat histogram
        
        G4String fileName = ArgumentInterpreter::GetOutPutFileName();
        G4String outputEdepFileName=fileName+"_edep"+".csv";// each worker will have each own writting file
        ofstream file;
        std::remove(outputEdepFileName); // first, remove the old file if the file exists
        file.open(outputEdepFileName);
        file<<"TotalEventNumber"<<","<<totalNumberOfEvent<<endl; // write totalNumberOfEvent to the file in first line
        file<<"Event ID"<<","<<"Cell ID"<<","<<"X/nm"<<","<<"Y/nm"<<","<<"Z/nm"<<","<<"Edep/eV"<<","<<"Cell Organelle"<<endl;// write titles
        file.close();

       }

  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndMaster(const G4Run*)
{
    // The code below is largely similar as the final results analysis done for local run in the function of EndWorker()
    // The results obtained below are  for the multithreading mode, so the global results were obtained by combing all the local run results.
    cout<<endl<<"The global results were as below: "<<endl<<endl;

    G4String fileName;
    fileName = ArgumentInterpreter::GetOutPutFileName();

    /////////////////////////////////////////////////////////////////////////////
    /// Dose analysis
    //////////////////////////
    
    G4RunDoseAnalysis runDoseAnalysis;
    runDoseAnalysis.ImportTotalParticleNumber(totalNumberOfEvent);
    
    const DetectorConstruction* detector= static_cast< const DetectorConstruction*> \
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    
    runDoseAnalysis.CellDoseTally(detector);
    
    std::map<int, double> cellDoseMeanMap;// map storing the mean cell total dose, the key is cellID, value is total dose of cell
    std::map<int, double> cellDoseStdMap; // map storing the statistical error of mean cell total dose, the key is cellID
    std::map<int, double> nucleusDoseMeanMap;
    std::map<int, double> nucleusDoseStdMap;
    
    cellDoseMeanMap = runDoseAnalysis.GetCellDoseMeanMap();
    cellDoseStdMap = runDoseAnalysis.GetCellDoseStdMap();
    nucleusDoseMeanMap = runDoseAnalysis.GetNucleusDoseMeanMap();
    nucleusDoseStdMap = runDoseAnalysis.GetNucleusDoseStdMap();
    
    
    runDoseAnalysis.RadiationTransportDoseOutPut(fileName); // write out the dose tally information to file
        

    /////////////////////////////////////////////////////////////////////////////////////////
        
    G4RunDNADamageAnalysis runDNADamageAnalysis;
    
    runDNADamageAnalysis.ImportTotalParticleNumber(totalNumberOfEvent);
    runDNADamageAnalysis.CellDNADamageTally(cellDoseMeanMap,cellDoseStdMap);
    
    std::map<int, double> cellSSBMeanMap;
    std::map<int, double> cellSSBStdMap;
    std::map<int, double> cellDSBMeanMap;
    std::map<int, double> cellDSBStdMap;
    
    cellSSBMeanMap = runDNADamageAnalysis.GetCellSSBMeanMap();
    cellSSBStdMap = runDNADamageAnalysis.GetCellSSBStdMap();
    cellDSBMeanMap = runDNADamageAnalysis.GetCellDSBMeanMap();
    cellDSBStdMap = runDNADamageAnalysis.GetCellDSBStdMap();

    runDNADamageAnalysis.RadiationTransportDNADamageOutPut(fileName); // write out the DNA damage tally information
    
    cout<<endl<<"The radiation transport simulation results are as below:"<<endl;
    cout<<endl;
    cout<<"Cell Dose Information: "<<endl<<endl;
    
    for (std::map<int, double>::iterator mitr_cell=cellDoseMeanMap.begin();mitr_cell!=cellDoseMeanMap.end();mitr_cell++)
    {
        int space=6;
        int space1=10;
        int accuracy=5;
        cout<<"mean dose "<<cellDoseMeanMap[mitr_cell->first]<<endl;
        
        cout<<"cellID: "<<std::setw(space)<<left<<mitr_cell->first<<" "<<"TotalDose: "<<std::setw(space1)<<std::scientific\
        <<std::setprecision(accuracy)<<cellDoseMeanMap[mitr_cell->first]<<" +- "<<std::setw(space1)<<std::setprecision(accuracy)<<\
        cellDoseStdMap[mitr_cell->first]<<" Gy"<<" "<<" RError: "<<std::fixed<<std::setw(space1)<<\
        cellDoseStdMap[mitr_cell->first]/cellDoseMeanMap[mitr_cell->first]*100<<" %  "<<\
        "NucleusDose: "<<std::scientific<<std::setw(space1)<<std::setprecision(accuracy)<<\
        nucleusDoseMeanMap[mitr_cell->first]<<" +- "<<std::setw(space1)<<std::setprecision(accuracy)<<nucleusDoseStdMap[mitr_cell->first]<<" Gy"<<\
        " RError: "<<std::fixed<<std::setw(space1)<<nucleusDoseStdMap[mitr_cell->first]/nucleusDoseMeanMap[mitr_cell->first]*100<<" %"<<endl;
        
    }
    cout<<endl;
    cout<<"DNA Damage Information: "<<endl<<endl;
    for (std::map<int, double>::iterator mitr_cell=cellSSBMeanMap.begin();mitr_cell!=cellSSBMeanMap.end();mitr_cell++)
    {
        int space=6;
        int space1=8;
        int accuracy=5;
        cout<<"cellID: "<<std::setw(space)<<left<<mitr_cell->first<<" "<<"SSB: "<<std::setw(space1)<<std::scientific\
        <<std::setprecision(accuracy)<<cellSSBMeanMap[mitr_cell->first]<<" +- "<<std::setw(space1)<<std::setprecision(accuracy)<<\
        cellSSBStdMap[mitr_cell->first]<<" RError: "<<std::fixed<<std::setw(space1)<<cellSSBStdMap[mitr_cell->first]/cellSSBMeanMap[mitr_cell->first]\
        *100<<" %  "\
        <<" DSB: "<<std::setw(space1)<<std::setprecision(accuracy)<<std::scientific<<\
        cellDSBMeanMap[mitr_cell->first]<<" +- "<<std::setw(space1)<<std::setprecision(accuracy)<<cellDSBStdMap[mitr_cell->first]\
        <<" RError: "<<std::fixed<<std::setw(space1)<<cellDSBStdMap[mitr_cell->first]/cellDSBMeanMap[mitr_cell->first]*100<<" %"<<endl;
    }
    
    /////After run, clear the data in whole run
        
    runDoseAnalysis.ClearRunData();
    runDNADamageAnalysis.ClearRunData();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndWorker(const G4Run* run)
{
  if(fDebug)
  {
    G4cout << "===================================" << G4endl;
    G4cout << "================ RunAction::EndWorker" << G4endl;
    PrintRunInfo(run);
    G4cout << "===================================" << G4endl;
  }

  int nofEvents = run->GetNumberOfEvent();
  if ( nofEvents == 0 )
  {
    if(fDebug)
    {
      G4cout << "================ NO EVENTS TREATED IN THIS RUN ==> Exit"
             << G4endl;
    }
    return;
  }
  


//   ///////////////
//   // Write Histo
//   //
//   WriteHistogram();
// 
//   ///////////////
//   // Complete cleanup
//   //
//   delete G4AnalysisManager::Instanc++ static function undefined reference toce();

  ///////////////
  // Printouts
  //
  std::map<const G4ParticleDefinition*, int>&
  particlesCreatedInWorld = fpTrackingAction->GetNParticlesCreatedInWorld();

  G4cout << "Number and type of particles created outside region \"Target\" :"
         << G4endl;

  PrintNParticles(particlesCreatedInWorld);

  G4cout << "_______________________" << G4endl;

  std::map<const G4ParticleDefinition*, int>&
  particlesCreatedInTarget = fpTrackingAction->GetNParticlesCreatedInTarget();

  G4cout << "Number and type of particles created in region \"Target\" :"
         << G4endl;

  PrintNParticles(particlesCreatedInTarget);
  
  /////////////////////////////////////////////////////// Write information to external file
  
    bool sequential = (G4RunManager::GetRunManager()->GetRunManagerType() == 
                     G4RunManager::sequentialRM);
    if(isMaster && sequential == true)// When the code is implemented in sequential way do this, otherwise donot do this
    {
        
        G4String simulationMode = ArgumentInterpreter::GetSimulationMode();
        if (simulationMode == "out")
        {
            G4String fileName;
            fileName = ArgumentInterpreter::GetOutPutFileName();
    
            bool postProcessMarker=ArgumentInterpreter::GetPostProcessMarker();
            if (postProcessMarker)
            {
                    ///////////////
                    // Write Histo
                    //
                    WriteHistogram();

                    ///////////////
                    // Complete cleanup
                    //
                    delete G4AnalysisManager::Instance();
            }
            
            G4RunDoseAnalysis runDoseAnalysis;
            runDoseAnalysis.ImportTotalParticleNumber(totalNumberOfEvent);
            
            const DetectorConstruction* detector= static_cast< const DetectorConstruction*> \
            (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
            
            runDoseAnalysis.CellDoseTally(detector);
            
            std::map<int, double> cellDoseMeanMap;// map storing the mean cell total dose, the key is cellID, value is total dose of cell
            std::map<int, double> cellDoseStdMap; // map storing the statistical error of mean cell total dose, the key is cellID
            std::map<int, double> nucleusDoseMeanMap;
            std::map<int, double> nucleusDoseStdMap;
            
            cellDoseMeanMap = runDoseAnalysis.GetCellDoseMeanMap();
            cellDoseStdMap = runDoseAnalysis.GetCellDoseStdMap();
            nucleusDoseMeanMap = runDoseAnalysis.GetNucleusDoseMeanMap();
            nucleusDoseStdMap = runDoseAnalysis.GetNucleusDoseStdMap();
            
            
            runDoseAnalysis.RadiationTransportDoseOutPut(fileName); // write out the dose tally information to file
            

        /////////////////////////////////////////////////////////////////////////////////////////
            
            G4RunDNADamageAnalysis runDNADamageAnalysis;
            
            runDNADamageAnalysis.ImportTotalParticleNumber(totalNumberOfEvent);
            runDNADamageAnalysis.CellDNADamageTally(cellDoseMeanMap,cellDoseStdMap);
            
            std::map<int, double> cellSSBMeanMap;
            std::map<int, double> cellSSBStdMap;
            std::map<int, double> cellDSBMeanMap;
            std::map<int, double> cellDSBStdMap;
            
            cellSSBMeanMap = runDNADamageAnalysis.GetCellSSBMeanMap();
            cellSSBStdMap = runDNADamageAnalysis.GetCellSSBStdMap();
            cellDSBMeanMap = runDNADamageAnalysis.GetCellDSBMeanMap();
            cellDSBStdMap = runDNADamageAnalysis.GetCellDSBStdMap();

            runDNADamageAnalysis.RadiationTransportDNADamageOutPut(fileName); // write out the DNA damage tally information
            
            cout<<"The radiation transport simulation results are as below:"<<endl;
            cout<<endl;
            cout<<"Cell Dose Information: "<<endl<<endl;
            
            for (std::map<int, double>::iterator mitr_cell=cellDoseMeanMap.begin();mitr_cell!=cellDoseMeanMap.end();mitr_cell++)
            {
                int space=6;
                int space1=10;
                int accuracy=5;
                
                cout<<"cellID: "<<std::setw(space)<<left<<mitr_cell->first<<" "<<"TotalDose: "<<std::setw(space1)<<std::scientific\
                <<std::setprecision(accuracy)<<cellDoseMeanMap[mitr_cell->first]<<" +- "<<std::setw(space1)<<std::setprecision(accuracy)<<\
                cellDoseStdMap[mitr_cell->first]<<" Gy"<<" "<<" RError: "<<std::fixed<<std::setw(space1)<<\
                cellDoseStdMap[mitr_cell->first]/cellDoseMeanMap[mitr_cell->first]*100<<" %  "<<\
                "NucleusDose: "<<std::scientific<<std::setw(space1)<<std::setprecision(accuracy)<<\
                nucleusDoseMeanMap[mitr_cell->first]<<" +- "<<std::setw(space1)<<std::setprecision(accuracy)<<nucleusDoseStdMap[mitr_cell->first]<<" Gy"<<\
                " RError: "<<std::fixed<<std::setw(space1)<<nucleusDoseStdMap[mitr_cell->first]/nucleusDoseMeanMap[mitr_cell->first]*100<<" %"<<endl;
                
            }
            cout<<endl;
            cout<<"DNA Damage Information: "<<endl<<endl;
            for (std::map<int, double>::iterator mitr_cell=cellSSBMeanMap.begin();mitr_cell!=cellSSBMeanMap.end();mitr_cell++)
            {
                int space=6;
                int space1=8;
                int accuracy=5;
                cout<<"cellID: "<<std::setw(space)<<left<<mitr_cell->first<<" "<<"SSB: "<<std::setw(space1)<<std::scientific\
                <<std::setprecision(accuracy)<<cellSSBMeanMap[mitr_cell->first]<<" +- "<<std::setw(space1)<<std::setprecision(accuracy)<<\
                cellSSBStdMap[mitr_cell->first]<<" RError: "<<std::fixed<<std::setw(space1)<<cellSSBStdMap[mitr_cell->first]/cellSSBMeanMap[mitr_cell->first]\
                *100<<" %  "\
                <<" DSB: "<<std::setw(space1)<<std::setprecision(accuracy)<<std::scientific<<\
                cellDSBMeanMap[mitr_cell->first]<<" +- "<<std::setw(space1)<<std::setprecision(accuracy)<<cellDSBStdMap[mitr_cell->first]\
                <<" RError: "<<std::fixed<<std::setw(space1)<<cellDSBStdMap[mitr_cell->first]/cellDSBMeanMap[mitr_cell->first]*100<<" %"<<endl;
            }
            /////After run, clear the data in whole run

            runDoseAnalysis.ClearRunData();
            runDNADamageAnalysis.ClearRunData();

        }     
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::InitializeWorker(const G4Run*)
{
  RunInitManager::Instance()->Initialize();

  if (fpTrackingAction == 0)
  {
    fpTrackingAction = (TrackingAction*) G4RunManager::GetRunManager()->
        GetUserTrackingAction();

    if(fpTrackingAction == 0 && isMaster == false)
    {
      G4ExceptionDescription exDescrption ;
      exDescrption << "fpTrackingAction is a null pointer. "
          "Has it been correctly initialized ?";
      G4Exception("RunAction::BeginOfRunAction",
          "RunAction001",FatalException, exDescrption);
    }
  }

  fInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....com..

void RunAction::CreateHistogram()
{
     
  // Book histograms, ntuple

  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in Analysis.hh
    
    G4String simulationMode=ArgumentInterpreter::GetSimulationMode();


    
   if (simulationMode=="gui") return;

  G4cout << "##### Create analysis manager " << "  " << this << G4endl;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//  if(!analysisManager->IsActive()) {return; }

  G4cout << "Using " << analysisManager->GetType() <<
      " analysis manager" << G4endl;
      
      
  analysisManager->SetVerboseLevel(1);


    G4String fileName = ArgumentInterpreter::GetOutPutFileName();
    analysisManager->OpenFile(fileName);
  

  // Creating ntuple

  analysisManager->CreateNtuple("microdosimetry", "physics");
  analysisManager->CreateNtupleDColumn("flagParticle");
  analysisManager->CreateNtupleDColumn("flagProcess");
  analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->CreateNtupleDColumn("totalEnergyDeposit");
  analysisManager->CreateNtupleDColumn("stepLength");
  analysisManager->CreateNtupleDColumn("kineticEnergyDifference");
  analysisManager->CreateNtupleDColumn("cellID");
    
  analysisManager->FinishNtuple();
  
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::WriteHistogram()
{
    G4String simulationMode = ArgumentInterpreter::GetSimulationMode();

   if (simulationMode=="gui") return;

  // print histogram statistics
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
//  if(!analysisManager->IsActive()) {return; }

  // save histograms
  //
  analysisManager->Write();
  analysisManager->CloseFile();
  
  
  if(fDebug)
  {
    G4cout << "================ ROOT FILES HAVE BEEN WRITTEN"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::PrintRunInfo(const G4Run* run)
{
  G4cout << "================ Run is = "
         << run->GetRunID() << G4endl;
  G4cout << "================ Run type is = "
         << G4RunManager::GetRunManager()->GetRunManagerType() << G4endl;
  G4cout << "================ Event processed = "
         << run->GetNumberOfEventToBeProcessed() << G4endl;
  G4cout << "================ Nevent = "
         << run->GetNumberOfEvent() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrintNParticles(std::map<const G4ParticleDefinition*, int>& container)
{
  std::map<const G4ParticleDefinition*, int>::iterator it;
  for(it = container.begin() ; it != container.end(); it ++)
  {
    G4cout << "N " << it->first->GetParticleName() << " : "
        << it->second << G4endl;
  }
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
