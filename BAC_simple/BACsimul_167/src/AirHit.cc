#include "AirHit.hh"
#include "BACDetectorConstruction.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4VisAttributes.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "G4AttDefStore.hh"

#include <iomanip>

G4Allocator<AirHit> AirHitAllocator;

AirHit::AirHit()
  :xyz(0.,0.,0.),tof(0.)
{}

AirHit::AirHit(G4ThreeVector& axyz, G4double t)
  :xyz(axyz),tof(t)
{}

AirHit::~AirHit()
{}

AirHit::AirHit(G4ThreeVector &axyz, G4double t, G4int pid, G4double angle,G4int cn,G4double en_depo)
  :xyz(axyz), tof(t),particleID(pid),angleget(angle),copynum(cn),energydeposit(en_depo)
{}
