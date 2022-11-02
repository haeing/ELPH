#ifndef PMT28SD_h
#define PMT28SD_h 1

#include "G4VSensitiveDetector.hh"
#include "PMT28Hit.hh"
#include "TGraph.h"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class PMT28SD : public G4VSensitiveDetector
{
public:
  PMT28SD(G4String name);
  //PMT28SD(G4String name, ParamMan*);
  virtual ~PMT28SD();

  virtual void Initialize( G4HCofThisEvent *HCTE );
  virtual G4bool ProcessHits(G4Step* astep, G4TouchableHistory *ROhist);
  virtual void EndOfEvent(G4HCofThisEvent *HCTE);

  void SetQETable();
  TGraph* QETable;
  

private:
  PMT28HitsCollection *Pmt28Collection;

};

#endif
