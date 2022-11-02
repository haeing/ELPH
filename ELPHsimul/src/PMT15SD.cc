#include "PMT15SD.hh"
#include "PMT15Hit.hh"

#include "G4HCofThisEvent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TouchableHandle.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include "TGraph.h"

PMT15SD::PMT15SD(G4String name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("Pmt15Collection");
  SetQETable();
}

PMT15SD::~PMT15SD()
{}

void PMT15SD::Initialize(G4HCofThisEvent *HCTE)
{
  static int HCID = -1;
  Pmt15Collection = new PMT15HitsCollection(GetName(), collectionName[0]);

  if (HCID<0) HCID = GetCollectionID(0);

  HCTE->AddHitsCollection(HCID, Pmt15Collection);


  

}

G4bool PMT15SD::ProcessHits(G4Step *astep, G4TouchableHistory *ROhist)
{
  
  const G4StepPoint* preStepPoint = astep-> GetPreStepPoint();
  G4Track* atrack = astep->GetTrack();
  G4int pid = astep->GetTrack()->GetDefinition()-> GetPDGEncoding();
  G4ThreeVector worldPos =preStepPoint->GetPosition();
  G4ThreeVector pos = preStepPoint->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(worldPos);

  G4double energy = atrack->GetTotalEnergy();
  G4ThreeVector hitmom = atrack->GetMomentum()/CLHEP::eV;
  G4double E_p = hitmom.mag();

  const G4double h = 6.628e-34;
  const G4double c = 3.0e+8;
  G4double wavelength = ((h*c)/(energy*1.6e-13))*(1e+9); //nm

  G4int copynumber = preStepPoint->GetTouchableHandle()->GetCopyNumber();

  
  
  
  auto hitTime = astep->GetPreStepPoint()->GetGlobalTime();


  atrack->SetTrackStatus(fStopAndKill);

  G4double random = G4UniformRand();

  G4double qe_value = QETable->Eval(E_p);
  if (qe_value < 0 || qe_value > 1) qe_value = 0;

  if(random<qe_value){
    PMT15Hit* ahit = new PMT15Hit(pos,worldPos, hitTime,pid,wavelength, copynumber);
    Pmt15Collection->insert(ahit);
   }
  return true;
      

}


void PMT15SD::EndOfEvent(G4HCofThisEvent *HCTE)
{
  Pmt15Collection->PrintAllHits();

}

void PMT15SD::SetQETable(){
  const int num_qe = 14;

  //Original
  //G4double pmt15_ep1[num_qe] = {1.82,1.88,1.91,1.97,2.07,2.10,2.25,2.48,2.76,2.92,3.35,3.76,4.28,4.43};
  //G4double pmt15_effi[num_qe] = {0.001,0.0025,0.005,0.01,0.025,0.05,0.1,0.16,0.22,0.25,0.29,0.25,0.1,0.025};

  //UV
  G4double pmt15_ep1[num_qe] = {1.82,1.88,1.91,1.97,2.07,2.10,2.25,2.48,2.76,2.92,3.35,3.76,5.64,6.53};
  G4double pmt15_effi[num_qe] = {0.001,0.0025,0.005,0.01,0.025,0.05,0.1,0.16,0.22,0.25,0.29,0.25,0.1,0.05};
  

  QETable = new TGraph(num_qe,pmt15_ep1,pmt15_effi);
  //QETable->SetName("PMT15_eff");
  

  
}
  
