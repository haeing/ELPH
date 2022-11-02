#include "PMT28SD.hh"
#include "PMT28Hit.hh"

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

PMT28SD::PMT28SD(G4String name)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("Pmt28Collection");
  SetQETable();
}

PMT28SD::~PMT28SD()
{}

void PMT28SD::Initialize(G4HCofThisEvent *HCTE)
{
  static int HCID = -1;
  Pmt28Collection = new PMT28HitsCollection(GetName(), collectionName[0]);

  if (HCID<0) HCID = GetCollectionID(0);

  HCTE->AddHitsCollection(HCID, Pmt28Collection);


  

}

G4bool PMT28SD::ProcessHits(G4Step *astep, G4TouchableHistory *ROhist)
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
    PMT28Hit* ahit = new PMT28Hit(pos,worldPos, hitTime,pid,wavelength, copynumber);
    Pmt28Collection->insert(ahit);
    }
  return true;
      

}


void PMT28SD::EndOfEvent(G4HCofThisEvent *HCTE)
{
  Pmt28Collection->PrintAllHits();

}

void PMT28SD::SetQETable(){
  const int num_qe = 10;

  
  G4double pmt28_ep1[num_qe] = {1.91,1.97,2.07,2.18,2.48,3.02,3.54,3.81,4.20,4.35};
  G4double pmt28_effi[num_qe] = {0.01,0.02,0.05,0.1,0.2,0.3,0.36,0.3,0.2,0.1};

  QETable = new TGraph(num_qe,pmt28_ep1,pmt28_effi);
  //QETable->SetName("PMT28_eff");
  

  
}
  
