#ifndef BACDetectorConstruction_aerogel3_mppc4_hh
#define BACDetectorConstruction_aerogel3_mppc4_hh

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSensitiveDetector;
class G4PVPlacement;
class G4VisAttributes;


class BACDetectorConstruction_aerogel3_mppc4 : public G4VUserDetectorConstruction
{
public:

  BACDetectorConstruction_aerogel3_mppc4();
  virtual ~BACDetectorConstruction_aerogel3_mppc4();

  virtual void ConstructSDandField();

  virtual G4VPhysicalVolume* Construct();

private:



  


  
  G4VPhysicalVolume* physWorld;
  G4VPhysicalVolume* physDetect;


  G4LogicalVolume* Aero1LW;
  G4LogicalVolume* Aero2LW;
  G4LogicalVolume* Aero3LW;
  G4LogicalVolume* HolderLW;
  G4LogicalVolume* BehindLW;
  G4LogicalVolume* Ae_sideLW;
  G4LogicalVolume* Behind_filmLW;
  G4LogicalVolume* BottomLW;
  G4LogicalVolume* ReflectLW;
  G4LogicalVolume* FilmLW;
  G4LogicalVolume* SideLW;
  G4LogicalVolume* MPPCLW;



  std::vector<G4VisAttributes*> fVisAttributes;



  
};

#endif
