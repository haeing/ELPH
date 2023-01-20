#ifndef BACPrimaryGeneratorAction_hh
#define BACPrimaryGeneratorAction_hh

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <TRandom3.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>

class G4ParticleGun;
class G4Event;
class G4Box;

class BACPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  BACPrimaryGeneratorAction(const G4String &parameter1, const G4String &parameter2, const G4String &parameter3);
  virtual ~BACPrimaryGeneratorAction();
  virtual void GeneratePrimaries(G4Event *anEvent);
  virtual void GenerateBeamKaonMBr(G4Event *anEvent, G4ThreeVector D, G4ThreeVector P,G4String particle);
  virtual void ReadBeamProfile(G4ThreeVector &X, G4ThreeVector &P);
  const G4ParticleGun* GetParticleGun() const {return fParticleGun;}



private:

  G4double pa1;
  G4double pa2;
  G4double pa3;
  
  G4ParticleGun* fParticleGun;
  G4ParticleTable* particleTable;
  Int_t bp_file_ndata;
  Int_t bp_nAccess=0;
  const double mass_kaonm = 0.493677;
  const double mass_pim = 0.139570;
  G4String particle = "pion";
  G4double energy;

  G4String parameter1;
  G4String parameter2;
  G4String parameter3;

};

#endif