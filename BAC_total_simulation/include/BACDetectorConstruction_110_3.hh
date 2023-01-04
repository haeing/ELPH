#ifndef BACDetectorConstruction_110_3_hh
#define BACDetectorConstruction_110_3_hh

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4PVPlacement;
class G4VisAttributes;


class BACDetectorConstruction_110_3 : public G4VUserDetectorConstruction
{
public:
  BACDetectorConstruction_110_3();
  virtual ~BACDetectorConstruction_110_3();

  virtual void ConstructSDandField();

  virtual G4VPhysicalVolume* Construct();

private:


  
  G4VPhysicalVolume* physWorld;
  G4VPhysicalVolume* physDetect;
  G4LogicalVolume* Aero1LW;
  G4LogicalVolume* Aero2LW;
  G4LogicalVolume* Aero3LW;
  G4LogicalVolume* BlackLW;
  G4LogicalVolume* trdworldLW;
  G4LogicalVolume* TrdLW;
  G4LogicalVolume* UpReflLW;
  G4LogicalVolume* DownReflLW;
  G4LogicalVolume* MPPCLW;
  G4LogicalVolume* WinstonLW;
  G4LogicalVolume* CCPCLW;
  G4LogicalVolume* UpLW;
  G4LogicalVolume* BottomLW;
  G4LogicalVolume* DetLW;
  G4LogicalVolume* ACLW;
  G4LogicalVolume* MPPCLW1;
  G4LogicalVolume* MPPC1LW;
  G4LogicalVolume* BehindLW;
  G4LogicalVolume* FrameLW;
  


  G4LogicalVolume* Part1LW;
  G4LogicalVolume* Part2LW;
  G4LogicalVolume* CheckLW;
  G4LogicalVolume* ReflectLW;
  G4LogicalVolume* Reflect1LW;
  G4LogicalVolume* ReflectBLW;
  G4LogicalVolume* DetectLW;
  G4LogicalVolume* SideLW;
  G4LogicalVolume* HolderLW;



  std::vector<G4VisAttributes*> fVisAttributes;

  G4int version;
  //const int version = stoi(version_put);

  //G4String version_put;
  G4String version_in;
  G4String num_aero;
  //G4String parameter1;
  //G4String parameter2;
  
  //BACLens *Lens;


  
};

#endif
