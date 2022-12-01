#ifndef BACStackingAction_H
#define BACStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4HCofThisEvent;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BACStackingAction : public G4UserStackingAction
{
  public:
    BACStackingAction();
    virtual ~BACStackingAction();

  public:
    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
    virtual void NewStage();
    virtual void PrepareNewEvent();


  private:
    G4int fScintillationCounter;
    G4int fCerenkovCounter;
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
