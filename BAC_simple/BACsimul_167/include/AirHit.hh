#ifndef AirHit_hh
#define AirHit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"

class G4AttDef;
class G4AttValue;

class AirHit : public G4VHit
{
private:
  G4ThreeVector xyz;
  G4int particleID;
  G4double tof;
  G4double angleget;
  G4double copynum;
  G4double energydeposit;

  
public:
  AirHit();
  AirHit(G4ThreeVector& axyz, G4double t);
  AirHit(G4ThreeVector& axyz, G4double t, G4int pid, G4double angle,G4int cn,G4double en_depo);
  virtual ~AirHit();

  AirHit(const AirHit& right);
  const AirHit& operator=(const AirHit &right);
  //G4bool operator==(const AirHit &right) const;

  void *operator new(size_t);
  void operator delete(void *aHit);

  //virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  //virtual std::vector<G4AttValue>* CreateAttValues() const;

  const G4ThreeVector& GetPosition() const {return xyz;}
  G4double GetTOF() const {return tof;}
  
  G4int GetParticleID() const {return particleID;}
  G4double GetAngle() const {return angleget;}
  G4int GetCopyNum() const {return copynum;}
  G4double GetEnergyDeposit() const{return energydeposit;}
  //G4int GetNum() const {return count_ce;}

  //void SetTime(G4double dt) {fTime = dt;}


  //void IncPhotonCount() {fPhotons++;}
  //G4double GetPhotonCount() {return fPhotons;}
  //void ClearPhotonCount() {fPhotons=0.0;}

  //void SetPos(G4ThreeVector xyz) {fPos = xyz;}
  //G4ThreeVector GetPos() const {return fPos;}
  
  
};

inline AirHit::AirHit(const AirHit& right)
  : G4VHit()
{
  xyz = right.xyz;
  tof = right.tof;
  particleID = right.particleID;
  angleget = right.angleget;
  copynum = right.copynum;
  energydeposit = right.energydeposit;

}

inline const AirHit& AirHit::operator=
(const AirHit& right)
{
  xyz = right.xyz;
  tof = right.tof;
  particleID = right.particleID;
  angleget = right.angleget;
  copynum = right.copynum;
  energydeposit = right.energydeposit;

  return *this;
}

//using AirHitsCollection = G4THitsCollection<AirHit>;
typedef G4THitsCollection<AirHit> AirHitsCollection;
//extern G4ThreadLocal G4Allocator<AirHit>* AirHitAllocator;
extern G4Allocator<AirHit> AirHitAllocator;

inline void* AirHit::operator new(size_t)
{
  //return (void*)AirHitAllocator-> MallocSingle();
  return static_cast<void*>(AirHitAllocator.MallocSingle());

}

inline void AirHit::operator delete(void* aHit)
{
  //AirHitAllocator->FreeSingle((AirHit*) aHit);
  AirHitAllocator.FreeSingle(static_cast<AirHit*> (aHit));
}

#endif
  
