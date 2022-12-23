#ifndef AirSD_h
#define AirSD_h 1

#include "G4VSensitiveDetector.hh"
#include "AirHit.hh"
//#include "TGraph.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class AirSD : public G4VSensitiveDetector
{
public:
  AirHitsCollection *AirCollection;

public:
  AirSD(G4String name);
  //AirSD(G4String name, ParamMan*);
  virtual ~AirSD();

  virtual void Initialize( G4HCofThisEvent *HCTE);
  virtual G4bool ProcessHits(G4Step* astep, G4TouchableHistory *ROhist);
  virtual void EndOfEvent(G4HCofThisEvent *HCTE);
};

#endif
