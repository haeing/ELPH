#ifndef PMT28Hit_hh
#define PMT28Hit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"

class G4AttDef;
class G4AttValue;

class PMT28Hit : public G4VHit
{
public:
  PMT28Hit();
  PMT28Hit(G4ThreeVector& axyz, G4double t);
  PMT28Hit(G4ThreeVector& axyz, G4ThreeVector& wxyz, G4double t, G4int pid, G4double Wavelength, G4int cn);
  virtual ~PMT28Hit();

  //copy constructor & assignment operator
  PMT28Hit(const PMT28Hit& right);
  const PMT28Hit& operator=(const PMT28Hit &right);
  //G4bool operator==(const PMT28Hit &right) const;

  //new/delete operators
  void *operator new(size_t);
  void operator delete(void *aHit);

  /*
  virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  virtual std::vector<G4AttValue>* CreateAttValues() const;
  */

  const G4ThreeVector& GetPosition() const { return xyz; }
  const G4ThreeVector& GetWorldPosition() const { return worldxyz; }
  G4double GetTOF() const { return tof; }
  //void SetEdep(G4double de) {fEdep = de;}
  //void AddEdep(G4double de) {fEdep += de;}
  //G4double GetEdep() const {return fEdep;}
  G4double GetWavelength() const {return wavelengthMP; }

  G4int GetParticleID() const {return particleID;}

  //void SetTime(G4double dt) {fTime = dt;}
  //G4double GetTime() const {return fTime;}

  void IncPhotonCount() {fPhotons++;}
  G4double GetPhotonCount() {return fPhotons;}
  void ClearPhotonCount() {fPhotons=0.0;}

  G4int GetCopyNum() const {return copynum;}




  //void SetPos(G4ThreeVector xyz) {fPos = xyz;}
  //G4ThreeVector GetPos() const {return fPos;}
  
private:
  G4ThreeVector xyz;
  G4ThreeVector worldxyz;
  G4double tof;
  //G4int fId;
  G4int particleID;
  G4double wavelengthMP;
  //G4double fEdep;
  //G4double fTime;
  G4double fPhotons;
  //G4ThreeVector fPos;
  G4double copynum;

};

inline PMT28Hit::PMT28Hit(const PMT28Hit& right)
  : G4VHit()
{
  xyz = right.xyz;
  worldxyz = right.worldxyz;
  tof = right.tof;
  wavelengthMP = right.wavelengthMP;
  copynum = right.copynum;
}

inline const PMT28Hit& PMT28Hit::operator=
(const PMT28Hit& right)
{
  xyz = right.xyz;
  worldxyz = right.worldxyz;
  tof = right.tof;
  wavelengthMP = right.wavelengthMP;
  copynum = right.copynum;
  return *this;
}

typedef G4THitsCollection<PMT28Hit> PMT28HitsCollection;
extern G4Allocator<PMT28Hit> PMT28HitAllocator;


/*
using PMT28HitsCollection = G4THitsCollection<PMT28Hit>;

extern G4ThreadLocal G4Allocator<PMT28Hit>* PMT28HitAllocator;
*/
inline void* PMT28Hit::operator new(size_t)
{
  return static_cast<void*>(PMT28HitAllocator.MallocSingle());
}

inline void PMT28Hit::operator delete(void* aHit)
{
  PMT28HitAllocator.FreeSingle(static_cast<PMT28Hit*> (aHit));
}

#endif
  
