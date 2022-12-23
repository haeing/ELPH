#include "AirSD.hh"
#include "AirHit.hh"
#include "typeinfo"

#include <TMath.h>
#include <TRandom.h>

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
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessType.hh"
#include "G4StepStatus.hh"

AirSD::AirSD(G4String name)
  : G4VSensitiveDetector(name)//, fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert("AirCollection"); //have to modify
}

AirSD::~AirSD()
{}

void AirSD::Initialize(G4HCofThisEvent *HCTE)
{
  static int HCID = -1;
  AirCollection = new AirHitsCollection(SensitiveDetectorName, collectionName[0]);

  if(HCID<0) HCID = GetCollectionID(0);
  
  HCTE->AddHitsCollection(HCID, AirCollection);
  
}

G4bool AirSD::ProcessHits(G4Step *astep, G4TouchableHistory *ROhist)
{
  const G4StepPoint* preStepPoint = astep-> GetPreStepPoint();
  const G4StepPoint* postStepPoint = astep-> GetPostStepPoint();
  G4Track* atrack = astep->GetTrack();
  G4int pid = astep->GetTrack()->GetDefinition()->GetPDGEncoding();
  G4ThreeVector pos =preStepPoint->GetPosition();

  //Angle
  G4double angle;
  G4ThreeVector mom;
  mom = postStepPoint -> GetMomentum();
  angle = TMath::ACos(mom.z()/mom.mag())*180/TMath::Pi();

  G4double edep = astep->GetTotalEnergyDeposit();

      

  G4int copynumber = preStepPoint->GetTouchableHandle()->GetCopyNumber();

  
  G4double time = astep->GetPreStepPoint()->GetGlobalTime();
  G4double time_re = gRandom -> Gaus(time,0.1);
  G4double tof = time+time_re;
				   
  


  AirHit* ahit = new AirHit(pos,tof,pid,angle,copynumber,edep);

  AirCollection->insert(ahit);

  return true;  
}

void AirSD::EndOfEvent(G4HCofThisEvent* HCTE)
{

  AirCollection->PrintAllHits();


}

