#include "BACDetectorConstruction_aerogel3_mppc4.hh"
#include "BACPrimaryGeneratorAction.hh"
#include "BACRunAction.hh"
#include "BACEventAction.hh"
#include "BACStackingAction.hh"
#include "BACAnalysisManager.hh"
#include "AeroSD.hh"
#include "AeroHit.hh"
#include "MPPCSD.hh"
#include "MPPCHit.hh"

//#ifdef G4Multithreded
//#include "G4MTRunManager.hh"
//#else
#include "G4RunManager.hh"
//#endif

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4DeexPrecoParameters.hh"
#include "G4NuclearLevelData.hh"
//#include "G4RunManagerFactory.hh"

#include "G4VisExecutive.hh"
#include "G4ScoringManager.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  CLHEP::HepRandom::setTheSeed((unsigned)time(NULL));

  G4String histname;

  G4String par1_put;
  G4String par2_put;
  G4String par3_put;
  G4String par4_put;

  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {

    //std::cout<<"Please enter the version"<<std::endl;
    ui = new G4UIExecutive(argc, argv);
    histname = "geant4_test.root";
    par1_put ="0";
    par2_put = "0";
    par3_put = "780";
    par4_put = "1";
  }



  
  

  else if(argc>=3){
    histname = argv[2];
    par1_put = argv[3];
    par2_put = argv[4];
    par3_put = argv[5];
    par4_put = argv[6];
  }

  G4int geo_version = stod(par4_put);
  

  // Optionally: choose a different Random engine...
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
  
  // Construct the default run manager
  //
  G4RunManager* runManager = new G4RunManager;
  
  
  //    G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  // Set mandatory initialization classes
  //
  // Detector construction
  if(geo_version ==1)runManager->SetUserInitialization(new BACDetectorConstruction_aerogel3_mppc4());







  

  // Physics list
  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  G4VModularPhysicsList* physicsList = new QGSP_BERT;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  //opticalPhysics->SetScintillationYieldFactor(1.);
  //opticalPhysics->SetScintillationExcitationRatio(0.);
  //opticalPhysics->SetWLSTimeProfile("delta");
  opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  //opticalPhysics->Configure(kCerenkov,false);
  //opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);
  opticalPhysics->Configure(kScintillation,false);
  physicsList->RegisterPhysics(opticalPhysics);

  runManager->SetUserInitialization(physicsList);

  BACAnalysisManager *anaMan = new BACAnalysisManager(histname);
  //BACPrimaryGeneratorAction *priGen = new BACPrimaryGeneratorAction();
  BACPrimaryGeneratorAction *priGen = new BACPrimaryGeneratorAction(par1_put,par2_put,par3_put);
  BACRunAction *runAction = new BACRunAction(anaMan);
  BACEventAction *eventAction = new BACEventAction(anaMan);
  BACStackingAction *stackAction = new BACStackingAction();

  runManager->SetUserAction(priGen);
  runManager->SetUserAction(runAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(stackAction);
    
  // User action initialization
  //runManager->SetUserInitialization(new BACActionInitialization());
  runManager->Initialize();

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);

  }
  else { 
    // interactive mode
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;

  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !


  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
