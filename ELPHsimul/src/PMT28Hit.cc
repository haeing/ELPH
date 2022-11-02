#include "PMT28Hit.hh"
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

G4Allocator<PMT28Hit> PMT28HitAllocator;

PMT28Hit::PMT28Hit()
  :xyz(0.,0.,0.),tof(0.)
{}

PMT28Hit::PMT28Hit(G4ThreeVector& axyz, G4double t)
  :xyz(axyz),tof(t)
{}

PMT28Hit::~PMT28Hit()
{}

PMT28Hit::PMT28Hit(G4ThreeVector &axyz, G4ThreeVector &wxyz, G4double t, G4int pid, G4double Wavelength, G4int cn)
  :xyz(axyz), worldxyz(wxyz), tof(t),particleID(pid), wavelengthMP(Wavelength), copynum(cn)
{}

/*
const PMT28Hit& PMT28Hit::operator=(const PMT28Hit &right)
{
  fTime = right.fTime;
  //fId = right.fId;
  fPhotons = right.fPhotons;
  fPos = right.fPos;
  return *this;
}



const std::map<G4String,G4AttDef>* PMT28Hit::GetAttDefs() const
{
  G4bool isNew;
  auto store = G4AttDefStore::GetInstance("PMT28Hit",isNew);

  if (isNew){
    //(*store)["ID"] = G4AttDef("ID","ID","Physics","","G4int");
    (*store)["Time"] = G4AttDef("Time","Time","Physics","ns","G4double");
    (*store)["Hits"] = G4AttDef("Hit","Photon Hit","Physics","G4BestUnit","G4double");
    (*store)["Pos"] = G4AttDef("Pos","Position","Physics","mm""G4ThreeVector");
  }
  return store;
}

std::vector<G4AttValue>* PMT28Hit::CreateAttValues() const
{
  auto values = new std::vector<G4AttValue>;
  values ->push_back(G4AttValue("Energy",G4BestUnit(fEdep,"Energy"),""));
  values ->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));
  //values ->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(fId),""));
  values->push_back(G4AttValue("Hits",G4BestUnit(fPhotons,"Hits"),""));
  values ->push_back(G4AttValue("Pos",G4BestUnit(fPos,"Length"),""));
  return values;
}

*/



