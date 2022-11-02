#include "BACDetectorConstruction.hh"
#include "AeroSD.hh"
#include "MPPCSD.hh"
#include "PMT15SD.hh"
#include "PMT28SD.hh"


#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4ExtrudedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4UnionSolid.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trap.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"
#include "TString.h"
#include "TMath.h"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4RotationMatrix.hh"

double theta = 62*TMath::Pi()/180;
G4double r_enter = 12*mm;
//double theta = 45*TMath::Pi()/180;
//G4double r_enter = 8*mm;

double sint = TMath::Sin(theta);
double cost = TMath::Cos(theta);

double f(G4double r, G4double z){
  return pow(r*cost+z*sint,2)+2*r_enter*pow(1+sint,2)*r-2*r_enter*cost*pow(2+sint,2)*z-pow(r_enter,2)*(1+sint)*(3+sint);
}


BACDetectorConstruction::BACDetectorConstruction(const G4String &version_put,const G4String &num_aerogel, const G4String &par1,const G4String &par2,const G4String &par3)
  : G4VUserDetectorConstruction(),version_in(version_put), num_aero(num_aerogel),parameter1(par1), parameter2(par2), parameter3(par3)
{
}

BACDetectorConstruction::~BACDetectorConstruction()
{
  for (auto visAttributes: fVisAttributes){
    delete visAttributes;
  }
}


G4VPhysicalVolume* BACDetectorConstruction::Construct()
{

  version=stoi(version_in);
  G4int num_aero_i = stoi(num_aero);

  

  if(version==1){
    reflect_part_length_d = stod(parameter1);
    light_guide_length_d = stod(parameter2);
    middle_length_d = stod(parameter3);
  }
  
  G4NistManager* nist = G4NistManager::Instance();
  G4bool checkOverlaps = true;

  if(version==2){
    theta1 = stod(parameter1)*degree;
    theta2 = stod(parameter2)*degree;
    theta3 = stod(parameter3)*degree;
  }

  if(version==3){
    p = stod(parameter1);
    ref_z = stod(parameter2);
    ref_theta = stod(parameter3)*degree;
  }

  if(version==4){
    p = 36;
    Dpartz = stod(parameter1)*cm;
    ref_theta = stod(parameter2)*degree;
    mppc_theta = stod(parameter3)*degree;

  }
  


  //Material---------------------------------------------------------------

  G4Material* world_mat = nist-> FindOrBuildMaterial("G4_AIR");
  G4Material* Al = nist-> FindOrBuildMaterial("G4_Al");

  G4Element* C = new G4Element("Carbon","C",6.,12.011*g/mole);
  G4Element* H = new G4Element("Hydrogen","H",1.,1.00794*g/mole);
  G4Element* O = new G4Element("Oxygen","O",15.,15.9994*g/mole);
  //G4Element* Al = new G4Element("Aluminum","Al",13.,26.981538*g/mole);
  G4Element* B = new G4Element("Boron","B",5.,10.811*g/mole);
  G4Element* Na = new G4Element("Na","Na",11.,23.0*g/mole);
  G4Element* Si = new G4Element("Silicon","Si",14.,28.0844*g/mole);
  G4Element* Cl = new G4Element("Chlorine","Cl",17.,35.453*g/mole);
  G4Element* K = new G4Element("Potassium","K",19.,39.093*g/mole);
  G4Element* F = new G4Element("Fluorine","F",9.,18.998*g/mole);
  

  G4Material *Aerogel = new G4Material("Aerogel",0.2000*g/cm3,2);
  Aerogel->AddElement(Si,1);
  Aerogel->AddElement(O,2);

  G4Material* blacksheet = new G4Material("blacksheet",0.95*g/cm3,2);
  blacksheet->AddElement(C,1);
  blacksheet->AddElement(H,2);

  G4Material *Mylar = new G4Material("Mylar", 1.39*g/cm3, 3);
  Mylar->AddElement(C, 5);
  Mylar->AddElement(H, 4);
  Mylar->AddElement(O, 2);

  G4Material* Tyvek = new G4Material("Tyvek",0.38*g/cm3,2);
  Tyvek->AddElement(C,1);
  Tyvek->AddElement(H,2);

  G4Material* Epoxi = new G4Material("Epoxi",1.1*g/cm3,4);
  Epoxi->AddElement(C, 21);
  Epoxi->AddElement(H, 25);
  Epoxi->AddElement(O,  5);
  Epoxi->AddElement(Cl, 1);


  G4Material* Acrylic = new G4Material("Acrylic", 1.19*g/cm3, 3);
  Acrylic->AddElement(C, 5);
  Acrylic->AddElement(H, 8);
  Acrylic->AddElement(O, 2);

  G4Material* Teflon = new G4Material("Teflon",2.2*g/cm3,2);
  Teflon->AddElement(C,2);
  Teflon->AddElement(F,4);
  



  //Property--------------------------------------------------------------

  //air property---------------------------------------------------------

   
  G4MaterialPropertiesTable* prop_air = new G4MaterialPropertiesTable();
  G4double air_ep[] = {1.3*eV,7.*eV};
  G4double air_rindex[] = {1.0,1.0};
  prop_air->AddProperty("RINDEX",air_ep,air_rindex,2)->SetSpline(true);
  world_mat->SetMaterialPropertiesTable(prop_air);


  //Aerogel Property------------------------------------------------------
  G4MaterialPropertiesTable* prop_aerogel = new G4MaterialPropertiesTable();

  

  G4double aerogel_ep[] = {1.3*eV,7.*eV};
  //G4double aerogel_abs[] = {20*mm,20*mm};
  G4double aerogel_abs[] = {50*mm,50*mm};
  //G4double aerogel_abs[] = {300*mm,300*mm};
  //G4double aerogel_rindex[]={1.10,1.10};
  G4double aerogel_rindex[]={1.167,1.1670};
  //G4double aerogel_rindex[]={1.05,1.05};
  //G4double aerogel_ray[] = {6.16*pow(10,10),6.16*pow(10,10)};

  assert(sizeof(aerogel_ep_abs)==sizeof(aerogel_abs));
  //const G4int num_aerogel = sizeof(aerogel_ep_abs)/sizeof(G4double);
  
  prop_aerogel->AddProperty("RINDEX",aerogel_ep,aerogel_rindex,2)->SetSpline(true);
  prop_aerogel->AddProperty("ABSLENGTH",aerogel_ep,aerogel_abs,2)->SetSpline(true);
  //prop_aerogel->AddProperty("RAYLEIGH",aerogel_ep,aerogel_ray,2)->SetSpline(true);
  //prop_aerogel->AddProperty("FASTCOMPONENT",scin_ep1,scin_fast,numentries_scin1)->SetSpline(true);
  //prop_aerogel->AddProperty("SLOWCOMPONENT",scin_ep1,scin_fast,numentries_scin1)->SetSpline(true);
  //prop_aerogel->AddConstProperty("SCINTILLATIONYIELD",10000./MeV);

  prop_aerogel->AddConstProperty("RESOLUTIONSCALE",1.0);
  prop_aerogel->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  prop_aerogel->AddConstProperty("SLOWTIMECONSTANT",2.8*ns);
  prop_aerogel->AddConstProperty("YIELDRATIO",0.8);
  Aerogel->SetMaterialPropertiesTable(prop_aerogel);


  //Black sheet property-------------------------------------------------
  G4MaterialPropertiesTable* prop_bs = new G4MaterialPropertiesTable();

  G4double bs_ep[] = {1.3*eV,7.*eV};
  G4double bs_abs[] = {1.0e-9*cm,1.0e-9*cm};
  G4double bs_rindex[]={1.6,1.6};
  
  prop_bs->AddProperty("RINDEX",bs_ep,bs_rindex,2)->SetSpline(true);
  prop_bs->AddProperty("ABSLENGTH",bs_ep,bs_abs,2)->SetSpline(true);
  blacksheet->SetMaterialPropertiesTable(prop_bs);


  //Mylar property
  G4double mylar_ep[]={1.4*eV,1.48*eV,1.52*eV,1.56*eV,1.6*eV,1.7*eV,1.8*eV,1.9*eV,2*eV,2.2*eV,2.4*eV,2.6*eV,2.8*eV,3*eV,3.4*eV,3.8*eV,4*eV,5*eV,6*eV,7*eV};

  G4double mylar_real[]={2.2802,2.6945,2.7668,2.7675,2.6154,2.1606,1.8301,1.5724,1.366,1.0728,0.8734,0.7278,0.6079,0.52135,0.39877,0.31474,0.28003,0.18137,0.12677,0.094236};
  G4double mylar_ima[]={8.1134,8.1878,8.2573,8.3866,8.4914,8.3565,8.0601,7.7354,7.4052,6.7839,6.2418,5.7781,5.3676,5.0008,4.3957,3.9165,3.7081,2.9029,2.3563,1.9519};

  G4double mylar_ep1[] = {1.3*eV,7.*eV};

  assert (sizeof(mylar_ep) == sizeof(mylar_real));
  assert (sizeof(mylar_ep) == sizeof(mylar_ima));
  
  const G4int numentries_mylar = sizeof(mylar_ep)/sizeof(G4double);
  G4double mylar_abs[] = {1.0e-9*cm,1.0e-9*cm};
  G4MaterialPropertiesTable* prop_mylar = new G4MaterialPropertiesTable();
  prop_mylar->AddProperty("REALRINDEX",mylar_ep,mylar_real,numentries_mylar)->SetSpline(true);
  prop_mylar->AddProperty("IMAGINARYRINDEX",mylar_ep,mylar_ima,numentries_mylar)->SetSpline(true);
  prop_mylar->AddProperty("ABSLENGTH",mylar_ep1,mylar_abs,2)->SetSpline(true);
  Mylar->SetMaterialPropertiesTable(prop_mylar);

  //Tyvek property
  G4MaterialPropertiesTable* prop_tyvek = new G4MaterialPropertiesTable();
  
  G4double tyvek_rindex[2]={1.5,1.5};
  G4double tyvek_ep[] = {1.6*eV,7.*eV};
  G4double tyvek_abs[] = {1.0e-9*cm,1.0e-9*cm};

  prop_tyvek->AddProperty("RINDEX",tyvek_ep,tyvek_rindex,2)->SetSpline(true);
  prop_tyvek->AddProperty("ABSLENGTH",tyvek_ep,tyvek_abs,2)->SetSpline(true);
  Tyvek->SetMaterialPropertiesTable(prop_tyvek);

  //Tyvek property
  G4MaterialPropertiesTable* prop_teflon = new G4MaterialPropertiesTable();
  
  G4double teflon_rindex[2]={1.35,1.35};
  G4double teflon_ep[] = {1.6*eV,7.*eV};
  G4double teflon_abs[] = {1.0e-9*cm,1.0e-9*cm};

  prop_teflon->AddProperty("RINDEX",teflon_ep,teflon_rindex,2)->SetSpline(true);
  prop_teflon->AddProperty("ABSLENGTH",teflon_ep,teflon_abs,2)->SetSpline(true);
  Teflon->SetMaterialPropertiesTable(prop_teflon);

  
  //MPPC property
  G4MaterialPropertiesTable* prop_mppc = new G4MaterialPropertiesTable();
  
  G4double mppc_rindex[2]={1.,1.};
  G4double mppc_ep[] = {1.6*eV,7.*eV};
  G4double mppc_abs[] = {1.0*cm,1.0*cm};
  //G4double mppc_abs[] = {1.0*cm,1.0*cm};

  prop_mppc->AddProperty("RINDEX",mppc_ep,mppc_rindex,2)->SetSpline(true);
  prop_mppc->AddProperty("ABSLENGTH",mppc_ep,mppc_abs,2)->SetSpline(true);
  Al->SetMaterialPropertiesTable(prop_mppc);

  //MPPC surface (Epoxi) property
  G4MaterialPropertiesTable* prop_epoxi = new G4MaterialPropertiesTable();
 
  G4double epoxi_rindex[2]={1.5,1.5};
  G4double epoxi_ep[] = {1.3*eV,7.*eV};
  G4double epoxi_abs[] = {1.0*cm,1.0*cm};


  prop_epoxi->AddProperty("RINDEX",epoxi_ep,epoxi_rindex,2)->SetSpline(true);
  prop_epoxi->AddProperty("ABSLENGTH",epoxi_ep,epoxi_abs,2)->SetSpline(true);
  Epoxi->SetMaterialPropertiesTable(prop_epoxi);


  //////////////////////////////////////////////////////////////////
//               ACRYLIC Optical properties
//////////////////////////////////////////////////////////////////

// Refractive index

  G4double lambda_min = 200*nm;
  G4double lambda_max = 700*nm;

  const G4int NENTRIES = 11 ;
  G4double LAMBDA_ACRYLIC[NENTRIES] ;


  G4double RINDEX_ACRYLIC[NENTRIES] ;
  G4double ENERGY_ACRYLIC[NENTRIES] ;

// Parameterization for refractive index of High Grade PMMA 

  G4double bParam[4] = {1760.7010,-1.3687,2.4388e-3,-1.5178e-6} ; 
  
  for(G4int i=0;i<NENTRIES; i++){
 
    LAMBDA_ACRYLIC[i] = lambda_min + i*(lambda_max-lambda_min)/float(NENTRIES-1) ;
    RINDEX_ACRYLIC[i] = 0.0 ;

    for (G4int jj=0 ; jj<4 ; jj++)
    {
      RINDEX_ACRYLIC[i] +=  (bParam[jj]/1000.0)*std::pow(LAMBDA_ACRYLIC[i]/nm,jj) ; 
    }

    ENERGY_ACRYLIC[i] =   CLHEP::h_Planck*CLHEP::c_light/LAMBDA_ACRYLIC[i] ;  // Convert from wavelength to energy ;
//  G4cout << ENERGY_ACRYLIC[i]/eV << " " << LAMBDA_ACRYLIC[i]/nm << " " << RINDEX_ACRYLIC[i] << G4endl ;

  }

  G4MaterialPropertiesTable *MPT_Acrylic = new G4MaterialPropertiesTable();
  MPT_Acrylic->AddProperty("RINDEX", ENERGY_ACRYLIC, RINDEX_ACRYLIC, NENTRIES);


// Absorption
  const G4int NENT = 25 ;
  G4double LAMBDAABS[NENT] = 
  {
    100.0,
    246.528671, 260.605103, 263.853516, 266.019104, 268.726105,    
    271.433136, 273.598724, 276.305725, 279.554138, 300.127380,    
    320.159241, 340.191101, 360.764343, 381.337585, 399.745239,    
    421.401276, 440.891724, 460.382172, 480.414001, 500.987274,    
    520.477722, 540.509583, 559.458618,
    700.0    
  } ;

  G4double ABS[NENT] =   // Transmission (in %) of  3mm thick PMMA 
  { 
    0.0000000,
    0.0000000,  5.295952,  9.657321, 19.937695, 29.283491, 
    39.252335, 48.598133, 58.255451, 65.109039, 79.439247,
    85.669785, 89.719627, 91.277260, 91.588783, 91.900307,
    91.588783, 91.277260, 91.277260, 91.588783, 91.588783,
    91.900307, 91.900307, 91.588783,
    91.5
  } ;


  MPT_Acrylic->AddProperty("ABSLENGTH", new G4MaterialPropertyVector()) ;
  for(G4int i=0;i<NENT; i++){
    G4double energy    = CLHEP::h_Planck*CLHEP::c_light/(LAMBDAABS[i]*nm) ;
    G4double abslength ;

    if (ABS[i] <= 0.0) {
      abslength = 1.0/kInfinity ;
    }
    else {
      abslength = -3.0*mm/(G4Log(ABS[i]/100.0)) ;
    }

    MPT_Acrylic->AddEntry("ABSLENGTH", energy, abslength);

  }

  Acrylic->SetMaterialPropertiesTable(MPT_Acrylic);
  

  
  //Geometry---------------------------------------------------------------

  //World------------------------------------------------------------------
  G4double world_size = 1*m;
  G4Box* solidWorld = new G4Box("World",world_size, world_size, world_size); 
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  physWorld = new G4PVPlacement(0,G4ThreeVector(), logicWorld, "World",0,false,0,checkOverlaps);


  G4double top_ae = 90*mm;
  G4double bottom_ae = 117*mm;
  G4double height_ae = 123*mm;
  G4double length_ae = 40*mm;


  if(version==1){
  //Size-----------------------------------------------------------------------
  G4double extra = 10*mm;
  G4double Aerox = 125.0 *mm+extra;
  G4double Aeroy = 125.0 *mm+extra;
  //G4double Aeroz = 12.0 *mm*num_aero_i;
  G4double Aeroz = 4*cm;

  
  G4double Aerox_real = 125.0 *mm;
  G4double Aeroy_real = 125.0 *mm;

  //G4double Aerox_real = 50.0 *mm;
  //G4double Aeroy_real = 50.0 *mm;
  G4double Aeroz_real = 12.0 *mm*num_aero_i;

  G4double reflect_thick = 0.3*mm;
  G4double mppc_thick = 1*mm;
  G4double air_thin = 0.2*mm;

  //G4double empty_part1_z = 12.0*num_aero_i*mm;
  G4double empty_part1_z = 0;
  G4double empty_part2_z = reflect_part_length_d*mm;
  //G4double empty_part2_z = Aeroy+air_thin*2;

  //original
  G4double trd_dxa = 2.4*cm;    //-z position x length
  //G4double trd_dxa = 4.8*cm;    //-z position x length
  G4double trd_dxb = Aerox*0.5+air_thin;
  //G4double trd_dxb = Aerox+air_thin*2;
  G4double trd_dya = 2.4*cm;
  G4double trd_dyb = empty_part2_z;                        
  G4double trd_dz  = light_guide_length_d*mm;


  //G4int numRZ = 33;


  G4double extra_trap = -10*mm;
  G4double top_trap =12.0*mm*num_aero_i+extra_trap;
  G4double bottom_trap = empty_part2_z;
  G4double middle_trap = middle_length_d*mm;
  G4double height_trap = Aeroy+air_thin*2;
  G4double length_trap = Aerox+air_thin*2;


  /*
  G4double top_ae = 90*mm;
  G4double bottom_ae = 117*mm;
  G4double height_ae = 123*mm;
  G4double length_ae = 40*mm;
  */



  //Part1-----------------------------------------------------------------------------
  G4Box* part1_cover = new G4Box("part1_cover",Aerox/2+air_thin+reflect_thick,Aeroy/2+air_thin+reflect_thick,(empty_part1_z+air_thin+reflect_thick)/2);
  G4Box* part1_hole = new G4Box("part1_hole",Aerox/2+air_thin,Aeroy/2+air_thin,(empty_part1_z+air_thin)/2);
  G4SubtractionSolid* Part1 = new G4SubtractionSolid("Part1",part1_cover,part1_hole,0,G4ThreeVector(0,0,+reflect_thick/2));

  Part1LW = new G4LogicalVolume(Part1,Mylar,"Part1");
  new G4PVPlacement(0,G4ThreeVector(0,0,-(empty_part1_z+air_thin+reflect_thick)/2),Part1LW,"Part1",logicWorld,false,0,checkOverlaps);
    

  
  //Part2 - Aerogel-------------------------------------------------------------------
  //G4Box* Aero = new G4Box("Aero",Aerox_real/2,Aeroy_real/2,Aeroz_real/2);
  G4Trap* Aero = new G4Trap("Aero",length_ae/2, 0*degree, 0*degree, height_ae/2, top_ae/2, bottom_ae/2, 0*degree, height_ae/2, top_ae/2, bottom_ae/2, 0*degree);
  AeroLW = new G4LogicalVolume(Aero,Aerogel,"Aero");
  new G4PVPlacement(0,G4ThreeVector(0,0*mm,3*mm+Aeroz_real/2),AeroLW,"Aero",logicWorld,false,0,checkOverlaps);


  // virtual plane to check angle
  /*
  G4Box* AC = new G4Box("AC",Aerox_real/2,Aeroy_real/2,5*mm);
  ACLW = new G4LogicalVolume(AC,world_mat, "AC");
  new G4PVPlacement(0,G4ThreeVector(0,0,3*mm+Aeroz_real+5*mm),ACLW,"AC",logicWorld,false,0,checkOverlaps);
  */



  //Part2---------------------------------------------------------------------------

    

  //new
  G4Box* part2_cover = new G4Box("part2_cover",Aerox/2+air_thin+reflect_thick,Aeroy/2+air_thin+reflect_thick+trd_dz,(empty_part2_z+reflect_thick)/2);


  std::vector<G4TwoVector> polygon(5);
  polygon[0].set(-0.5*bottom_trap,0.5*height_trap);
  polygon[1].set(-0.5*bottom_trap,-height_trap*0.5);
  polygon[2].set(0.5*bottom_trap,-height_trap*0.5);
  polygon[3].set(0.5*bottom_trap-middle_trap,0);
  polygon[4].set(0.5*bottom_trap,0.5*height_trap);


  G4TwoVector offsetA(0.,0.), offsetB(0.,0.);
  G4double scaleA=1., scaleB=1.;
  G4ExtrudedSolid* part2_hole1 = new G4ExtrudedSolid("part2_hole1",polygon,length_trap/2,offsetA,scaleA, offsetB, scaleB);




  G4double anx = 0;
  G4double any = -90*degree;
  G4double anz = 180*degree;
    
  G4double sina = TMath::Sin(anx);
  G4double cosa = TMath::Cos(anx);
  G4double sinb = TMath::Sin(any);
  G4double cosb = TMath::Cos(any);
  G4double sinc = TMath::Sin(anz);
  G4double cosc = TMath::Cos(anz);
  G4RotationMatrix *rot = new G4RotationMatrix(G4ThreeVector(cosb*cosc,sina*sinb*cosc-cosa*sinc,cosa*sinb*cosc+sina*sinc),
					       G4ThreeVector(cosb*sinc,sina*sinb*sinc+cosa*cosc,cosa*sinb*sinc-sina*cosc),
					       G4ThreeVector(-sinb,sina*sinb,cosa*cosb));

  G4SubtractionSolid* part2_cover_second = new G4SubtractionSolid("part2_cover_second",part2_cover,part2_hole1,rot,G4ThreeVector(0,0,-reflect_thick*0.5));
  G4RotationMatrix *rotY = new G4RotationMatrix();
  rotY->rotateY(90*degree);


  G4Trd* trd_hole =   new G4Trd("trd_hole",0.5*trd_dxa, 0.5*trd_dxb,0.5*trd_dya, 0.5*trd_dyb, 0.5*(trd_dz+reflect_thick));
  G4RotationMatrix *rotX = new G4RotationMatrix();
  rotX->rotateX(270*degree);
  G4RotationMatrix *rotX90 = new G4RotationMatrix();
  rotX90->rotateX(90*degree);


  //Two MPPC per each side
  
  //light guide bottom part
  G4SubtractionSolid* Part2_1 = new G4SubtractionSolid("Part2",part2_cover_second,trd_hole,rotX90,G4ThreeVector(-Aerox*0.25-air_thin*0.5,-(Aeroy/2+air_thin+reflect_thick/2+trd_dz/2),-reflect_thick*0.5));
  G4SubtractionSolid* Part2_2 = new G4SubtractionSolid("Part2",Part2_1,trd_hole,rotX90,G4ThreeVector(Aerox*0.25+air_thin*0.5,-(Aeroy/2+air_thin+reflect_thick/2+trd_dz/2),-reflect_thick*0.5));
  
  //light guide top part
  G4SubtractionSolid* Part2_3 = new G4SubtractionSolid("Part2",Part2_2,trd_hole,rotX,G4ThreeVector(-Aerox*0.25-air_thin*0.5,Aeroy/2+air_thin+reflect_thick/2+trd_dz/2,-reflect_thick*0.5));
  G4SubtractionSolid* Part2 = new G4SubtractionSolid("Part2",Part2_3,trd_hole,rotX,G4ThreeVector(Aerox*0.25+air_thin*0.5,Aeroy/2+air_thin+reflect_thick/2+trd_dz/2,-reflect_thick*0.5));


  
  //One MPPC per each side
  /*    
    G4SubtractionSolid* Part2_1 = new G4SubtractionSolid("Part2",part2_cover_second,trd_hole,rotX90,G4ThreeVector(0,-(Aeroy/2+air_thin+reflect_thick/2+trd_dz/2),-reflect_thick*0.5));
    G4SubtractionSolid* Part2 = new G4SubtractionSolid("Part2",Part2_1,trd_hole,rotX,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick/2+trd_dz/2,-reflect_thick*0.5));
  */



  G4Box* Block = new G4Box("Block",trd_dxb,trd_dz*0.5,1*mm);
  G4LogicalVolume* BlockLW =new G4LogicalVolume(Block,Mylar,"Block");

  

  Part2LW = new G4LogicalVolume(Part2,Mylar,"Part2");

  new G4PVPlacement(0,G4ThreeVector(0,0,(empty_part2_z+reflect_thick)/2),Part2LW,"Part2",logicWorld,false,0,checkOverlaps);
    

  //MPPC---------------------------------------------------------------------------
  //G4Box* MPPC = new G4Box("MPPC",8*cm,8*cm,mppc_thick/2);
  G4Box* MPPC = new G4Box("MPPC",12*mm,12*mm,mppc_thick/2);
  //G4Tubs* PMT15 = new G4Tubs("PMT15",0*mm,9*mm,mppc_thick/20,0,2*TMath::Pi());
  //G4Tubs* PMT28 = new G4Tubs("PMT15",0*mm,14*mm,mppc_thick/20,0,2*TMath::Pi());
  //MPPCLW = new G4LogicalVolume(MPPC,Al,"MPPC");
  MPPCLW = new G4LogicalVolume(MPPC,Epoxi,"MPPC");
  //PMT15LW = new G4LogicalVolume(PMT15,Epoxi,"PMT15");
  //PMT28LW = new G4LogicalVolume(PMT28,Epoxi,"PMT28");
  //new G4PVPlacement(rotX,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2,empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,0,checkOverlaps);

  /*
  new G4PVPlacement(rotX,G4ThreeVector(0,Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2,empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,0,checkOverlaps);
  new G4PVPlacement(rotX,G4ThreeVector(0,-(Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2),empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,0,checkOverlaps);


  */

  new G4PVPlacement(rotX,G4ThreeVector(-trd_dxb/2,Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2,empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,1,checkOverlaps);
  new G4PVPlacement(rotX,G4ThreeVector(trd_dxb/2,Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2,empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,2,checkOverlaps);
  new G4PVPlacement(rotX,G4ThreeVector(-trd_dxb/2,-(Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2),empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,4,checkOverlaps);
  new G4PVPlacement(rotX,G4ThreeVector(trd_dxb/2,-(Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2),empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,3,checkOverlaps); //setting number to be same with the comsic test
  


  /*
  new G4PVPlacement(rotX,G4ThreeVector(-1.2*cm,Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2,empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,1,checkOverlaps);
  new G4PVPlacement(rotX,G4ThreeVector(1.2*cm,Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2,empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,2,checkOverlaps);
  new G4PVPlacement(rotX,G4ThreeVector(-1.2*cm,-(Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2),empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,3,checkOverlaps);
  new G4PVPlacement(rotX,G4ThreeVector(1.2*cm,-(Aeroy/2+air_thin+reflect_thick+trd_dz+mppc_thick/2),empty_part2_z/2),MPPCLW,"MPPC",logicWorld,false,4,checkOverlaps);
  */

  
  }


  if(version==2){

    //size
    G4double Aerox_real = 125.0 *mm;
    G4double Aeroy_real = 125.0 *mm;
    G4double Aeroz_real = 12.0 *mm*num_aero_i;

    G4double extrax = 20*mm;
    G4double extray = 10*mm;
    G4double Aerox = 125.0 *mm+extrax;
    G4double Aeroy = 125.0 *mm+extray;

    

    //G4int numRZ = 30;
    G4int numRZ = 100;
    G4double win_thick = 1*mm;
    G4double mppc_thick = 1*mm;
    G4double re_thick = 1*mm;




    //Aerogel 
    G4Box* Aero = new G4Box("Aero",Aerox_real/2,Aeroy_real/2,Aeroz_real/2);
    AeroLW = new G4LogicalVolume(Aero,Aerogel,"Aero");
    new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,0*mm),AeroLW,"Aero",logicWorld,false,0,checkOverlaps);

    //Winston calculation
    
    G4double r_in[numRZ];
    G4double z[numRZ];
    G4double r_out[numRZ];
    
    for(int i=0;i<numRZ;i++){
      z[i] = i*mm;
    }
    /*
    r_in[0] = r_enter;
    r_out[0] = r_in[0]+win_thick;
    for(int i=1;i<numRZ;i++){
      G4double x = r_in[i-1];
      while(fabs(f(x,z[i]))>0.1){
	x+=0.0001;
      }
      r_in[i] = x;
      r_out[i] = x+win_thick;
    }
    */


    //Reflect
    /*
    G4Box* Reflect = new G4Box("Reflect",15*cm,4*cm,1*mm);
    G4Box* ReflectB = new G4Box("ReflectB",15*cm,8*cm,1*mm);
    ReflectLW = new G4LogicalVolume(Reflect,Mylar,"Reflect");
    ReflectBLW = new G4LogicalVolume(ReflectB,Mylar,"ReflectB");
    G4RotationMatrix *rotXT = new G4RotationMatrix();
    G4RotationMatrix *rotXB = new G4RotationMatrix();
    G4RotationMatrix *rotXR = new G4RotationMatrix();
    rotXT->rotateX(-10*degree);
    rotXB->rotateX(-45*degree);
    //rotXR->rotateX(50*degree);
    rotXR->rotateX(60*degree);
    new G4PVPlacement(rotXT,G4ThreeVector(0,4*cm,100*mm),ReflectLW,"Reflect",logicWorld,false,0,checkOverlaps);
    new G4PVPlacement(rotXB,G4ThreeVector(0,-2*cm,75*mm),ReflectBLW,"ReflectB",logicWorld,false,0,checkOverlaps);
    new G4PVPlacement(rotXR,G4ThreeVector(0,9.5*cm,90*mm),ReflectLW,"Reflect",logicWorld,false,0,checkOverlaps);
    */


    r_out[numRZ-1]=15*mm;
    //Combining reflectors
    std::vector<G4TwoVector> poly(8);
    poly[0].set(0,-Aeroy/2);
    poly[1].set(re_thick,-Aeroy/2);
    poly[2].set(Aeroy/(2*TMath::Tan(theta1))+re_thick,0);
    poly[3].set(Aeroy/2*(1/TMath::Tan(theta1)+1/TMath::Tan(theta2))+re_thick,Aeroy/2);
    poly[4].set(2*r_out[numRZ-1]*TMath::Cos(theta3)+re_thick,2*r_out[numRZ-1]*TMath::Sin(theta3)+Aeroy/2);
    poly[5].set(2*r_out[numRZ-1]*TMath::Cos(theta3),2*r_out[numRZ-1]*TMath::Sin(theta3)+Aeroy/2);
    poly[6].set(Aeroy/2*(1/TMath::Tan(theta1)+1/TMath::Tan(theta2)),Aeroy/2);
    poly[7].set(Aeroy/(2*TMath::Tan(theta1)),0);
    for(int i=0;i<8;i++){
      std::cout<<i<<" : "<<poly[i]<<std::endl;
    }


    G4TwoVector offsetA(0.,0.), offsetB(0.,0.);
    G4double scaleA=1., scaleB=1.;
    G4ExtrudedSolid* Reflect = new G4ExtrudedSolid("Reflect",poly,Aerox/2,offsetA,scaleA, offsetB, scaleB);
    ReflectLW = new G4LogicalVolume(Reflect,Mylar,"Reflect");
    G4RotationMatrix *rotY = new G4RotationMatrix();
    rotY->rotateY(+90*degree);
    new G4PVPlacement(rotY,G4ThreeVector(0,0,Aeroz_real/2),ReflectLW,"Reflect",logicWorld,false,0,checkOverlaps);

    //side reflector
    G4Box* Side = new G4Box("Side",1*mm,30*cm,30*cm);
    SideLW = new G4LogicalVolume(Side,Mylar,"Side");
    new G4PVPlacement(0,G4ThreeVector(-72.5*mm,0,0),SideLW,"Side",logicWorld,false,0,checkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(72.5*mm,0,0),SideLW,"Side",logicWorld,false,0,checkOverlaps);

     
    
      

    //Detection part

    G4RotationMatrix *rotd = new G4RotationMatrix();
    rotd->rotateX(270*degree+theta3);
    
    //G4Box* Detect = new G4Box("Detect",r_out[numRZ-1], r_out[numRZ-1], (z[numRZ-1]+mppc_thick)/2);
    G4Box* Detect = new G4Box("Detect",Aerox/2, r_out[numRZ-1], (z[numRZ-1]+mppc_thick)/2);
    DetectLW = new G4LogicalVolume(Detect, world_mat, "Detect");


    //Winston cone part
    /*
    G4Polycone* Winston = new G4Polycone("Winston",0,TMath::Pi()*2,numRZ,z,r_in,r_out);
    //G4Polycone* Winston = new G4Polycone("Winston",0,TMath::Pi()*2,numRZ_real,z_real,r_in_real,r_out_real);
    WinstonLW = new G4LogicalVolume(Winston,Mylar,"Winston");
    */
    //new G4PVPlacement(rotX90,G4ThreeVector(Aerox_real/4,Aeroy_real/2+z[numRZ-1]+3*cm,Aeroz_real/2),WinstonLW,"Winston",logicWorld,false,0,checkOverlaps);
    //new G4PVPlacement(rotX90,G4ThreeVector(-Aerox_real/4,Aeroy_real/2+z[numRZ-1]+3*cm,Aeroz_real/2),WinstonLW,"Winston",logicWorld,false,0,checkOverlaps);

    //new G4PVPlacement(0,G4ThreeVector(0,0,-z[numRZ-1]/2),WinstonLW,"Winston",DetectLW,false,0,checkOverlaps);
    //new G4PVPlacement(rotX90,G4ThreeVector(0,Aeroy_real/2+z[numRZ],0),WinstonLW,"Winston",logicWorld,false,0,checkOverlaps);


    //CCPC
    /*
    std::vector<G4TwoVector> polygon_in(4);
    polygon_in[0].set(-12*mm,12*mm);
    polygon_in[1].set(12*mm,12*mm);
    polygon_in[2].set(12*mm,-12*mm);
    polygon_in[3].set(-12*mm,-12*mm);

    std::vector<G4TwoVector> polygon_out(4);
    polygon_out[0].set(-12*mm-win_thick,12*mm+win_thick);
    polygon_out[1].set(12*mm+win_thick,12*mm+win_thick);
    polygon_out[2].set(12*mm+win_thick,-12*mm-win_thick);
    polygon_out[3].set(-12*mm-win_thick,-12*mm-win_thick);


    //std::vector<G4ExtrudedSolid::ZSection> zsections(numRZ);
    std::vector<G4ExtrudedSolid::ZSection> zsections_in;
    G4ExtrudedSolid::ZSection z_zero_in(0*mm, {0,0},1);
    zsections_in.push_back(z_zero_in);

    std::vector<G4ExtrudedSolid::ZSection> zsections_out;
    G4ExtrudedSolid::ZSection z_zero_out(0*mm, {0,0},1);
    zsections_out.push_back(z_zero_out);

    

    for(int i=1;i<numRZ;i++){
      G4ExtrudedSolid::ZSection z_i_in(i*mm,{0,0},r_in[i]/12*mm);
      zsections_in.push_back(z_i_in);

      G4ExtrudedSolid::ZSection z_i_out(i*mm,{0,0},r_out[i]/(12*mm+win_thick));
      zsections_out.push_back(z_i_out);
    }
    G4ExtrudedSolid* CCPC_in = new G4ExtrudedSolid("CCPCin",polygon_in,zsections_in);
    G4ExtrudedSolid* CCPC_out = new G4ExtrudedSolid("CCPCout",polygon_out,zsections_out);

    G4SubtractionSolid* CCPC = new G4SubtractionSolid("CCPC",CCPC_out,CCPC_in,0,G4ThreeVector(0,0,0));
    
    CCPCLW = new G4LogicalVolume(CCPC,Mylar,"CCPC");
    new G4PVPlacement(0,G4ThreeVector(0,0,-z[numRZ-1]/2),CCPCLW,"CCPC",DetectLW,false,0,checkOverlaps);
    */
    

    

    //MPPC---------------------------------------------------------------------------
    //G4Tubs* MPPC = new G4Tubs("MPPC",0*mm,1.2*cm,mppc_thick/2,0,2*TMath::Pi());
    //G4Box* MPPC = new G4Box("MPPC",1.2*cm,1.2*cm,mppc_thick/2);
    //G4Box* MPPC = new G4Box("Detect",Aerox/2,r_out[numRZ-1],mppc_thick/2);
    G4Box* MPPC = new G4Box("Detect",12*mm,12*mm,mppc_thick/2);
    MPPCLW = new G4LogicalVolume(MPPC,Epoxi,"MPPC");
    //new G4PVPlacement(0,G4ThreeVector(0,0,-z[numRZ-1]/2),MPPCLW,"MPPC",DetectLW,false,0,checkOverlaps);
    //new G4PVPlacement(0,G4ThreeVector(0,0,z[numRZ-1]/2),MPPCLW,"MPPC",DetectLW,false,1,checkOverlaps);
    G4double ml = 32;
    for(int i=0;i<4;i++){
      new G4PVPlacement(0,G4ThreeVector((-ml*1.5)*mm+(29*i)*mm,0,z[numRZ-1]/2),MPPCLW,"MPPC",DetectLW,false,5-i,checkOverlaps);
    }




    //Calculation of position
    G4double cross = TMath::Sqrt(r_out[numRZ-1]*r_out[numRZ-1]+pow((z[numRZ-1]+mppc_thick)/2,2));
    G4double theta_l = (180*TMath::ATan((z[numRZ-1]+mppc_thick)/(2*r_out[numRZ-1]))/TMath::Pi());
    G4double ly = Aeroy/2+cross*TMath::Sin(theta_l*degree+theta3);
    G4double lz = Aeroz_real/2+cross*TMath::Cos(theta_l*degree+theta3);

    

    
    //new G4PVPlacement(rotd,G4ThreeVector(Aerox_real/4,ly,lz), DetectLW, "Detect",logicWorld,false,0,checkOverlaps);
    //new G4PVPlacement(rotd,G4ThreeVector(-Aerox_real/4,ly,lz), DetectLW, "Detect",logicWorld,false,0,checkOverlaps);

    new G4PVPlacement(rotd,G4ThreeVector(0,ly,lz), DetectLW, "Detect",logicWorld,false,0,checkOverlaps);




  }


    if(version==3){

    //size
     
    G4double Aerox_real = 145.0 *mm;
    G4double Aeroy_real = 113.0 *mm;
    //G4double Aeroz_real = 12.0 *mm*num_aero_i;
    //G4double Aeroz_real = 12.0 *mm*3;
    G4double Aeroz_real = 33*mm;

    G4double extrax = 20*mm;
    G4double extray = 10*mm;
    G4double Aerox = 125.0 *mm+extrax;
    G4double Aeroy = 125.0 *mm+extray;

    

    //G4int numRZ = 30;
    G4int numRZ = 150;
    G4double win_thick = 1*mm;
    G4double mppc_thick = 1*mm;
    G4double re_thick = 1*mm;
    G4double mppc_place = 70*mm;




    //Aerogel 
    G4Box* Aero = new G4Box("Aero",Aerox_real/2,Aeroy_real/2,Aeroz_real/2);
    AeroLW = new G4LogicalVolume(Aero,Aerogel,"Aero");
    new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,0*mm),AeroLW,"Aero",logicWorld,false,0,checkOverlaps);


    //Reflector
    G4Box* Reflect = new G4Box("Reflect",Aerox_real/2,Aeroy_real/2,1*mm);
    ReflectLW = new G4LogicalVolume(Reflect,Teflon,"Reflect");

    G4Box* Reflect1 = new G4Box("Reflect1",Aeroz_real/2,Aeroy_real/2,1*mm);
    Reflect1LW = new G4LogicalVolume(Reflect1,Teflon,"Reflect1");

    G4RotationMatrix *rotR = new G4RotationMatrix();
    rotR->rotateY(90*degree);
    
    new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,Aeroz_real/2+1*mm),ReflectLW,"Reflect",logicWorld,false,0,checkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,-Aeroz_real/2-1*mm),ReflectLW,"Reflect",logicWorld,false,0,checkOverlaps);

    new G4PVPlacement(rotR,G4ThreeVector(Aerox_real/2+1*mm,0*mm),Reflect1LW,"Reflect1",logicWorld,false,0,checkOverlaps);
    new G4PVPlacement(rotR,G4ThreeVector(-Aerox_real/2-1*mm,0*mm),Reflect1LW,"Reflect1",logicWorld,false,0,checkOverlaps);

    
    

    
    
    


   
    //MPPC---------------------------------------------------------------------------

    G4Tubs* PMT15 = new G4Tubs("PMT15",0*mm,12.5*mm,mppc_thick,0,2*TMath::Pi());
    PMT15LW = new G4LogicalVolume(PMT15,Epoxi,"PMT15");
    
    //G4Box* MPPC = new G4Box("MPPC",12*mm,12*mm,mppc_thick/2);
    
    //G4Box* MPPC = new G4Box("MPPC",12*mm,12*mm,mppc_thick/2);
    //MPPCLW = new G4LogicalVolume(MPPC,Epoxi,"MPPC");
    G4RotationMatrix *rotM = new G4RotationMatrix();
    //rotM->rotateX(90*degree+theta3);
    rotM->rotateX(90*degree);
    new G4PVPlacement(rotM,G4ThreeVector(-41.5*mm,Aeroy_real/2+mppc_thick,0),PMT15LW,"PMT15",logicWorld,false,1,checkOverlaps);
    new G4PVPlacement(rotM,G4ThreeVector(0*mm,Aeroy_real/2+mppc_thick,0),PMT15LW,"PMT15",logicWorld,false,2,checkOverlaps);
    new G4PVPlacement(rotM,G4ThreeVector(41.5*mm,Aeroy_real/2+mppc_thick,0),PMT15LW,"PMT15",logicWorld,false,3,checkOverlaps);
    new G4PVPlacement(rotM,G4ThreeVector(-41.5*mm,-Aeroy_real/2-mppc_thick,0),PMT15LW,"PMT15",logicWorld,false,4,checkOverlaps);
    new G4PVPlacement(rotM,G4ThreeVector(0*mm,-Aeroy_real/2-mppc_thick,0),PMT15LW,"PMT15",logicWorld,false,5,checkOverlaps);
    new G4PVPlacement(rotM,G4ThreeVector(41.5*mm,-Aeroy_real/2-mppc_thick,0),PMT15LW,"PMT15",logicWorld,false,6,checkOverlaps);

    //Al frame
    G4Box* Frame_be = new G4Box("Frame_be",Aerox_real/2,mppc_thick/2,Aeroz_real/2);
    G4Tubs* PMThole = new G4Tubs("PMThole",0*mm,13*mm,mppc_thick,0,2*TMath::Pi());
    
    G4SubtractionSolid* Frame_1 = new G4SubtractionSolid("Frame_1",Frame_be,PMThole,rotM,G4ThreeVector(-41.5*mm,0,0));
    G4SubtractionSolid* Frame_2 = new G4SubtractionSolid("Frame_2",Frame_1,PMThole,rotM,G4ThreeVector(0,0,0));
    G4SubtractionSolid* Frame = new G4SubtractionSolid("Frame",Frame_2,PMThole,rotM,G4ThreeVector(41.5*mm,0,0));

    FrameLW = new G4LogicalVolume(Frame,Al,"Frame");
    new G4PVPlacement(0,G4ThreeVector(0,Aeroy_real/2+mppc_thick/2,0),FrameLW,"Frame",logicWorld,false,0,checkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0,-Aeroy_real/2-mppc_thick/2,0),FrameLW,"Frame",logicWorld,false,0,checkOverlaps);

    
    
    /*
    for(int i=0;i<5;i++){
      new G4PVPlacement(rotM,G4ThreeVector((29*(i-2))*mm,Aeroy/2+mppc_place/2*TMath::Sin(theta3),Aeroz_real/2+mppc_place/2*TMath::Cos(theta3)),MPPCLW,"MPPC",logicWorld,false,i+1,checkOverlaps);
    }
    */

    /*
    for(int i=0;i<11;i++){
      new G4PVPlacement(rotM,G4ThreeVector((12*(i-5))*mm,Aeroy/2+mppc_place/2*TMath::Sin(theta3),Aeroz_real/2+mppc_place/2*TMath::Cos(theta3)),MPPCLW,"MPPC",logicWorld,false,i+1,checkOverlaps);
    }
    */

    //new G4PVPlacement(rotM,G4ThreeVector(0,Aeroy/2+mppc_place/2*TMath::Sin(theta3),Aeroz_real/2+mppc_place/2*TMath::Cos(theta3)),MPPCLW,"MPPC",logicWorld,false,0,checkOverlaps);
    //new G4PVPlacement(rotM,G4ThreeVector(0,Aeroy/2,Aeroz_real/2+60*mm),MPPCLW,"MPPC",logicWorld,false,1,checkOverlaps);
    //new G4PVPlacement(0,G4ThreeVector(0,0,5*cm-mppc_thick/2),MPPCLW,"MPPC",DetectLW,false,0,checkOverlaps);

  }





    
    if(version==4){

    //size

    G4double Aerox_real = 105.0 *mm;
    G4double Aeroy_real = 58.0 *mm;
    //G4double Aeroz_real = 12.0 *mm*num_aero_i;
    G4double Aeroz_real = 3.1*cm;
    //ELPH
    //G4double Aeroz_real = 3*cm;


      /*
      G4double Aerox_real = 58.0 *mm;
      G4double Aeroy_real = 105.0 *mm;
      G4double Aeroz_real = 3.1*cm;
      */

    G4double extrax = 20*mm;
    G4double extray = 10*mm;
    G4double Aerox = 125.0 *mm+extrax;
    G4double Aeroy = 125.0 *mm+extray;

    

    //G4int numRZ = 30;
    G4int numRZ = 150;
    G4double win_thick = 1*mm;
    G4double mppc_thick = 1*mm;
    G4double re_thick = 1*mm;
    G4double mppc_place = 30*mm;





    //lower position
    //G4double lower = 23*mm;
    G4double lower = 0*mm;
    //Aerogel
    G4Box* thin = new G4Box("thin",Aerox_real/2,Aeroy_real/2,0.5*mm/2);
    G4Box* Aero1 = new G4Box("Aero1",Aerox_real/2,Aeroy_real/2,Aeroz_real/2);
    G4SubtractionSolid* Aero2 = new G4SubtractionSolid("Aero2",Aero1,thin,0,G4ThreeVector(0,0,-5*mm-0.25*mm));
    G4SubtractionSolid* Aero = new G4SubtractionSolid("Aero",Aero2,thin,0,G4ThreeVector(0,0,5*mm+0.25*mm));
    //G4Trap* Aero = new G4Trap("Aero",length_ae/2, 0*degree, 0*degree, height_ae/2, top_ae/2, bottom_ae/2, 0*degree, height_ae/2, top_ae/2, bottom_ae/2, 0*degree);
    AeroLW = new G4LogicalVolume(Aero,Aerogel,"Aero");
    G4RotationMatrix *rotAero = new G4RotationMatrix();
    //rotAero->rotateZ(180*degree);
    //new G4PVPlacement(rotAero,G4ThreeVector(0*mm,0*mm,5*mm),AeroLW,"Aero",logicWorld,false,0,checkOverlaps);

    //ELPH
    new G4PVPlacement(0,G4ThreeVector(0*mm,10*mm-lower,5*mm),AeroLW,"Aero",logicWorld,false,0,checkOverlaps);


    //ELPH reflector behind the aerogel
    G4Box* Behind = new G4Box("Behind",Aerox_real/2,Aeroy_real/2,1*mm);
    BehindLW = new G4LogicalVolume(Behind,Mylar,"Behind");
    //new G4PVPlacement(0,G4ThreeVector(0,10*mm-lower,-Aeroz_real/2-1*mm),BehindLW,"Behind",logicWorld,false,0,checkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0,10*mm-lower,-Aeroz_real/2-1*mm+5*mm),BehindLW,"Behind",logicWorld,false,0,checkOverlaps);

    //ELPH reflector bottom

    G4Box *Bottom = new G4Box("Bottom",Aerox_real/2,1*mm,Aeroz_real/2);
    BottomLW = new G4LogicalVolume(Bottom,Mylar,"Bottom");
    new G4PVPlacement(0,G4ThreeVector(0,-Aeroy_real/2+10*mm-3.5*mm-1*mm,5*mm),BottomLW,"Bottom",logicWorld,false,0,checkOverlaps);
    //new G4PVPlacement(0,G4ThreeVector(0,Aeroy/2,0),BottomLW,"Bottom",logicWorld,false,0,checkOverlaps);

    //Aerogel Holder
    G4Box *Box = new G4Box("Box",Aerox_real/2+3.5*mm,Aeroy_real/2+3.5*mm,Aeroz_real/2+3.5*mm/2);
    G4SubtractionSolid* Box1 = new G4SubtractionSolid("Box1",Box,Aero,0,G4ThreeVector(0,0,-3.5*mm/2));
    G4Box *Boxwin = new G4Box("Boxwin",Aerox_real/2-3.5*mm,Aeroy_real/2-3.5*mm,Aeroz_real/2+3.5*mm/2);
    G4SubtractionSolid* Box2 = new G4SubtractionSolid("Box2",Box1,Boxwin,0,G4ThreeVector(0,0,3.5*mm/2));
    G4Box *Boxwin2 = new G4Box("Boxwin2",Aerox_real,Aeroy_real/2-3.5*mm,Aeroz_real/2-3.5*mm/2);
    G4SubtractionSolid* Box3 = new G4SubtractionSolid("Box3",Box2,Boxwin2,0,G4ThreeVector(0,0,0));
    G4Box *Boxwin3 = new G4Box("Boxwin3",Aerox_real/2-3.5*mm,Aeroy_real,Aeroz_real/2-3.5*mm/2);
    G4SubtractionSolid* Holder = new G4SubtractionSolid("Holder",Box3,Boxwin3,0,G4ThreeVector(0,0,0));
    
    HolderLW = new G4LogicalVolume(Holder,blacksheet,"Holder");
    new G4PVPlacement(0,G4ThreeVector(0*mm,10*mm-lower,5*mm+3.5*mm/2),HolderLW,"Holder",logicWorld,false,0,checkOverlaps);
    



    //Parabola
    G4double x[numRZ];
    G4double y[numRZ];
    G4double x_out[numRZ];
    //G4double p = mppc_place/2;
    for(int i=0;i<numRZ;i++){
      x[i] = i;
      x_out[i] = (i+1);
      y[i] = x[i]*x[i]/(4*p);
      
    }

    std::vector<G4TwoVector> poly(2*numRZ);
    for(int i=0;i<numRZ;i++){
      poly[i].set(x[i]*mm,-y[i]*mm);
      poly[2*numRZ-1-i].set(x_out[i]*mm,-y[i]*mm);
    }

    
    G4TwoVector offsetA(0.,0.), offsetB(0.,0.);
    G4double scaleA=1., scaleB=1.;
    G4ExtrudedSolid* Reflect = new G4ExtrudedSolid("Reflect",poly,Aerox/2,offsetA,scaleA, offsetB, scaleB);
    ReflectLW = new G4LogicalVolume(Reflect,Mylar,"Reflect");
    G4RotationMatrix *rotY = new G4RotationMatrix();
    rotY->rotateY(+90*degree);
    rotY->rotateZ(ref_theta);

    G4double pary = Aeroy/2+mppc_place/2*TMath::Sin(theta3)+5*cm;
    G4double parz = Aeroz_real/2+mppc_place*TMath::Cos(theta3)+5*cm;

    new G4PVPlacement(rotY,G4ThreeVector(0,pary,parz),ReflectLW,"Reflect",logicWorld,false,0,checkOverlaps);

    



    //side reflector
    G4Box* Side = new G4Box("Side",1*mm,12*cm,5*cm);
    SideLW = new G4LogicalVolume(Side,Mylar,"Side");
    //new G4PVPlacement(0,G4ThreeVector(-82.5*mm,0,16.5*cm),SideLW,"Side",logicWorld,false,0,checkOverlaps);
    //new G4PVPlacement(0,G4ThreeVector(+82.5*mm,0,16.5*cm),SideLW,"Side",logicWorld,false,0,checkOverlaps);

    new G4PVPlacement(0,G4ThreeVector(-82.5*mm,0,7*cm),SideLW,"Side",logicWorld,false,0,checkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(+82.5*mm,0,7*cm),SideLW,"Side",logicWorld,false,0,checkOverlaps);

    //for efficiency of upper part
    G4Box* Up = new G4Box("Up",Aerox/2,15*mm,1*mm);
    UpLW = new G4LogicalVolume(Up,Mylar,"Up");
    G4RotationMatrix *rotU = new G4RotationMatrix();
    rotU->rotateX(0*degree);
    //new G4PVPlacement(rotU,G4ThreeVector(0,pary-mppc_place/2*TMath::Sin(ref_theta),parz-0.5*mm),UpLW,"Up",logicWorld,false,0,checkOverlaps);
    
     
    
      

   
    //MPPC---------------------------------------------------------------------------


    



  
    G4RotationMatrix *rotM = new G4RotationMatrix();
    rotM->rotateX(90*degree+mppc_theta);

    /*
    G4double Dparty = 8*cm;
    //G4double Dpartz = 3*cm;
    G4double mthick = 2*mm;

    //Calculation of position
    G4double cross = TMath::Sqrt(pow(Dparty/2,2)+pow(Dpartz/2,2));
    G4double theta_l = (180*TMath::ATan((Dpartz+mppc_thick)/Dparty)/TMath::Pi());
    G4double ly = Aeroy/2+cross*TMath::Sin(theta_l*degree+mppc_theta);
    G4double lz = Aeroz_real/2+cross*TMath::Cos(theta_l*degree+mppc_theta);

    G4Box* Det = new G4Box("Det",Aerox/4+mthick/2,Dparty/2+mthick/2,Dpartz/2+mppc_thick/2);
    DetLW = new G4LogicalVolume(Det,world_mat,"Det");
    new G4PVPlacement(rotM,G4ThreeVector(0,ly,lz), DetLW, "Det",logicWorld,false,0,checkOverlaps);


    G4double trd_dxb = 2.9*cm*4;

    G4double trd_dxa = Aerox;

    G4double trd_dyb = 2.7*cm;
    G4double trd_dya = Dparty;                        
    G4double trd_dz  = Dpartz;
    
    G4Trd* Trd_in =   new G4Trd("Trd_in",0.5*trd_dxa, 0.5*trd_dxb,0.5*trd_dya, 0.5*trd_dyb, 0.5*trd_dz);
    G4Trd* Trd_out =   new G4Trd("Trd_out",0.5*trd_dxa+mthick/2, 0.5*trd_dxb+mthick/2,0.5*trd_dya+mthick/2, 0.5*trd_dyb+mthick/2, 0.5*trd_dz);
    G4SubtractionSolid* Trd = new G4SubtractionSolid("Trd",Trd_out,Trd_in,0,G4ThreeVector(0,0,0));
    
    TrdLW = new G4LogicalVolume(Trd,Mylar,"Trd");
    new G4PVPlacement(0,G4ThreeVector(0,0,-mppc_thick/2),TrdLW,"Trd",DetLW,false,0,checkOverlaps);
    */

    G4Box* MPPC = new G4Box("MPPC",12*mm,12*mm,mppc_thick/2);

    //SNP
    //G4Box* MPPC = new G4Box("MPPC",100*mm,150*mm,mppc_thick/2);
    G4Box* MPPC1 = new G4Box("MPPC1",9*mm,12*mm,mppc_thick/2);
    
    //G4Box* MPPC = new G4Box("MPPC",1.5*mm,1.5*mm,mppc_thick/2);
    MPPCLW = new G4LogicalVolume(MPPC,Epoxi,"MPPC");
    MPPC1LW = new G4LogicalVolume(MPPC1,Epoxi,"MPPC1");

    /*
    G4Tubs* PMT15 = new G4Tubs("PMT15",0*mm,7.5*mm,mppc_thick/20,0,2*TMath::Pi());
    G4Tubs* PMT28 = new G4Tubs("PMT15",0*mm,14*mm,mppc_thick/20,0,2*TMath::Pi());
    PMT15LW = new G4LogicalVolume(PMT15,Epoxi,"PMT15");
    PMT28LW = new G4LogicalVolume(PMT28,Epoxi,"PMT28");
    */

    

    /*
    
    for(int i=0;i<4;i++){
      new G4PVPlacement(0,G4ThreeVector(-43.5*mm+(29*i)*mm,0,Dpartz/2),MPPCLW,"MPPC",DetLW,false,i+1,checkOverlaps);
    }
    */


    

    //originally 29
    G4double ml = 29;
    //mppc_theta = 20*degree;


    G4double factor = 0*mm;

    //E72 Real!!

    /*
     for(int i=0;i<4;i++){
      new G4PVPlacement(rotM,G4ThreeVector(-(ml*1.5)*mm+(ml*i)*mm,Aeroy/2+mppc_place*1.5*TMath::Sin(mppc_theta)+factor*TMath::Sin(mppc_theta),Aeroz_real/2+mppc_place*1.5*TMath::Cos(mppc_theta)+factor*TMath::Cos(mppc_theta)),MPPCLW,"MPPC",logicWorld,false,i+1,checkOverlaps);
    }
    */





      



    //SNP
    //factor = 3*cm;
    //new G4PVPlacement(rotM,G4ThreeVector(0*mm,Aeroy/2+mppc_place*1.5*TMath::Sin(mppc_theta)+factor*TMath::Sin(mppc_theta),Aeroz_real/2+mppc_place*1.5*TMath::Cos(mppc_theta)-factor*TMath::Cos(mppc_theta)),MPPCLW,"MPPC",logicWorld,false,1,checkOverlaps);
    

    //ELPH condition!!!!!!!
     
    for(int i=0;i<4;i++){
      new G4PVPlacement(rotM,G4ThreeVector(-(ml*1.5)*mm+(ml*i)*mm,Aeroy/2+mppc_place*1.5*TMath::Sin(mppc_theta)+3*cm,Aeroz_real/2+mppc_place*1.5*TMath::Cos(mppc_theta)-0.8*cm),MPPCLW,"MPPC",logicWorld,false,i+1,checkOverlaps);
    }


    
    for(int i=0;i<4;i++){
      /*
      if(i==0||i==3){
	new G4PVPlacement(rotM,G4ThreeVector(-(ml*1.5)*mm+(ml*i)*mm,Aeroy/2+mppc_place*1.5*TMath::Sin(mppc_theta),Aeroz_real/2+mppc_place*1.5*TMath::Cos(mppc_theta)),MPPCLW,"MPPC",logicWorld,false,i+1,checkOverlaps);
      }
      else if(i==1||i==2){
	new G4PVPlacement(rotM,G4ThreeVector(-(ml*1.5)*mm+(ml*i)*mm+3*mm,Aeroy/2+mppc_place*1.5*TMath::Sin(mppc_theta),Aeroz_real/2+mppc_place*1.5*TMath::Cos(mppc_theta)),MPPC1LW,"MPPC1",logicWorld,false,i+1,checkOverlaps);
      }
}
    }


   
    new G4PVPlacement(rotM,G4ThreeVector(-(ml*1.5)*mm+(ml*3)*mm,Aeroy/2+mppc_place*1.5*TMath::Sin(mppc_theta),Aeroz_real/2+mppc_place*1.5*TMath::Cos(mppc_theta)),PMT28LW,"PMT28",logicWorld,false,3+1,checkOverlaps);
      */
    }


    

    
    /*
    for(int i=0;i<5;i++){
      for(int j=0;j<4;j++){
	for(int k=0;k<4;k++){
	  new G4PVPlacement(0,G4ThreeVector(((29*(i-2))-9+6*j)*mm,(-9+6*k)*mm,0),MPPCLW,"MPPC",DetLW,false,i*16+j*4+k+1,checkOverlaps);
	}
      }
    }
    */




    /*
    for(int i=0;i<11;i++){
      new G4PVPlacement(rotM,G4ThreeVector((12*(i-5))*mm,Aeroy/2+mppc_place*1.5*TMath::Sin(mppc_theta),Aeroz_real/2+mppc_place*1.5*TMath::Cos(mppc_theta)),MPPCLW,"MPPC",logicWorld,false,i+1,checkOverlaps);
    }
    */


  }


    




    

  



  //visattributes------------------------------------------------
  auto visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  visAttributes -> SetVisibility(false);
  logicWorld->SetVisAttributes(visAttributes);
  //if(version==2||version==3)DetectLW->SetVisAttributes(visAttributes);
  //if(version==4)DetLW->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  
  visAttributes = new G4VisAttributes(G4Color::Blue()); 
  AeroLW->SetVisAttributes(visAttributes);
  if(version==3)PMT15LW->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Color::Red()); 
  MPPCLW->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Color::Brown()); 
  ReflectLW->SetVisAttributes(visAttributes);
  SideLW->SetVisAttributes(visAttributes);
  BehindLW->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  

  visAttributes = new G4VisAttributes(G4Colour(0.5,0.5,0.0));
  fVisAttributes.push_back(visAttributes);





  //Surface--------------------------------------------------------------
  //black sheet surface------------------
  
  G4OpticalSurface* surface_bs = new G4OpticalSurface("surface_bs");
  surface_bs->SetType(dielectric_dielectric);
  surface_bs->SetModel(unified);
  surface_bs->SetFinish(polishedtyvekair);

  G4MaterialPropertiesTable* sp_bs = new G4MaterialPropertiesTable();

  G4double bs_effi[] = {0.0,0.0};
  G4double bs_reflec[]= {0.05,0.05};
  G4double bs_specularLobe[] = {0.3,0.3};
  G4double bs_specularSpike[] = {0.2,0.2};
  G4double bs_backScatter[] = {0,0};

  sp_bs->AddProperty("EFFICIENCY",bs_ep,bs_effi,2)->SetSpline(true);
  sp_bs->AddProperty("REFLECTIVITY",bs_ep,bs_reflec,2)->SetSpline(true);
  sp_bs->AddProperty("SPECULARLOBECONSTANT",bs_ep,bs_specularLobe,2)->SetSpline(true);
  sp_bs->AddProperty("SPECULARSPIKECONSTANT",bs_ep,bs_specularSpike,2)->SetSpline(true);
  sp_bs->AddProperty("BACKSCATTERCONSTANT",bs_ep,bs_backScatter,2)->SetSpline(true);
  surface_bs->SetMaterialPropertiesTable(sp_bs);

  new G4LogicalSkinSurface("bs_surface",HolderLW,surface_bs);

  
  
  

  //mylar_al surface-------------------
  G4OpticalSurface* surface_mylar = new G4OpticalSurface("surface_mylar");
  surface_mylar->SetType(dielectric_metal);
  surface_mylar->SetFinish(polished);    
  surface_mylar->SetModel(unified);
  
  G4MaterialPropertiesTable* sp_mylar = new G4MaterialPropertiesTable();
  G4double mylar_reflec[] = {0.80,0.80};  //for metal, reflectivity is calculated using rindex, they use polarization, angle, energy
  G4double mylar_effi[] = {0.0,0.0};
  G4double mylar_specularLobe[] = {0.85,0.85};
  G4double mylar_specularSpike[]={0.87,0.87};
  G4double mylar_backScatter[] = {0,0};
  sp_mylar->AddProperty("EFFICIENCY",mylar_ep1,mylar_effi,2)->SetSpline(true);
  sp_mylar->AddProperty("REFLECTIVITY",air_ep,mylar_reflec,2)->SetSpline(true);
  sp_mylar->AddProperty("SPECULARLOBECONSTANT",mylar_ep1,mylar_specularLobe,2)->SetSpline(true);
  sp_mylar->AddProperty("SPECULARSPIKECONSTANT",mylar_ep1,mylar_specularSpike,2)->SetSpline(true);
  sp_mylar->AddProperty("BACKSCATTERCONSTANT",mylar_ep1,mylar_backScatter,2)->SetSpline(true);
  surface_mylar->SetMaterialPropertiesTable(sp_mylar);
  
  //new G4LogicalSkinSurface("mylar_surface",TrdLW,surface_mylar);

  new G4LogicalSkinSurface("mylar_surface",Part1LW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",Part2LW,surface_mylar);

  //221012 ichikawa-san's
  new G4LogicalSkinSurface("mylar_surface",ReflectLW,surface_mylar);
  
  new G4LogicalSkinSurface("mylar_surface",ReflectBLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",WinstonLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",CCPCLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",SideLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",UpLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",BottomLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",TrdLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",Reflect1LW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",BehindLW,surface_mylar);
  new G4LogicalSkinSurface("mylar_surface",FrameLW,surface_mylar);




  
  //tyvek surface--------------------------------------------------
  G4OpticalSurface* surface_tyvek = new G4OpticalSurface("surface_tyvek");
  surface_tyvek->SetType(dielectric_dielectric);
  surface_tyvek->SetModel(unified);
  surface_tyvek->SetFinish(groundfrontpainted);

  //teflon surface--------------------------------------------------
  G4OpticalSurface* surface_teflon = new G4OpticalSurface("surface_teflon");
  surface_teflon->SetType(dielectric_dielectric);
  surface_teflon->SetModel(unified);
  surface_teflon->SetFinish(groundfrontpainted);

  G4MaterialPropertiesTable* sp_teflon = new G4MaterialPropertiesTable();

  G4double teflon_reflec[] = {0.99,0.99};  //for metal, reflectivity is calculated using rindex, they use polarization, angle, energy
  G4double teflon_effi[] = {0.0,0.0};

  //sp_teflon->AddProperty("EFFICIENCY",air_ep,teflon_effi,2)->SetSpline(true);
  sp_teflon->AddProperty("REFLECTIVITY",air_ep,teflon_reflec,2)->SetSpline(true);

  surface_teflon->SetMaterialPropertiesTable(sp_teflon);
  

  //new G4LogicalSkinSurface("teflon_surface",ReflectLW,surface_teflon);
  //new G4LogicalSkinSurface("teflon_surface",Reflect1LW,surface_teflon);

  



  //MPPC surface--------------------------------------------------
  //G4OpticalSurface* surface_mppc = new G4OpticalSurface("surface_mppc",glisur, ground, dielectric_metal, polished);


  /*
  G4OpticalSurface* surface_mppc = new G4OpticalSurface("surface_mppc");
  //surface_mppc->SetType(dielectric_metal);
  surface_mppc->SetType(dielectric_dielectric);
  surface_mppc->SetFinish(polished);    
  surface_mppc->SetModel(unified);
  G4MaterialPropertiesTable* sp_mppc = new G4MaterialPropertiesTable();
  G4double mppc_reflec[]={0.0,0.0};
  G4double mppc_effi[]={1.0,1.0};
  
  //G4double mppc_ep1[] = {1.38*eV,1.43*eV,1.47*eV,1.51*eV,1.56*eV,1.61*eV,1.66*eV,1.7*eV,1.74*eV,1.79*eV,1.84*eV,1.88*eV,1.93*eV,1.97*eV,2*eV,2.06*eV,2.1*eV,2.15*eV,2.19*eV,2.24*eV,2.3*eV,2.33*eV,2.4*eV,2.47*eV,2.57*eV,2.7*eV,2.85*eV,2.96*eV,3.05*eV,3.15*eV,3.22*eV,3.28*eV,3.35*eV,3.41*eV,3.5*eV,3.57*eV,3.65*eV,3.67*eV,3.71*eV,3.73*eV,3.77*eV,3.79*eV,3.83*eV,3.9*eV};
  //G4double mppc_effi[] = {0.035,0.048,0.057,0.07,0.085,0.098,0.113,0.126,0.142,0.158,0.172,0.191,0.206,0.226,0.243,0.258,0.276,0.294,0.308,0.326,0.344,0.356,0.371,0.385,0.395,0.399,0.392,0.376,0.360,0.342,0.321,0.300,0.278,0.251,0.228,0.201,0.175,0.141,0.120,0.098,0.079,0.059,0.039,0.021};



  //assert (sizeof(mppc_ep1) == sizeof(mppc_effi));
  //const G4int numentries_mppc = sizeof(mppc_ep1)/sizeof(G4double);
    
  sp_mppc->AddProperty("REFLECTIVITY", mppc_ep, mppc_reflec,2);
  sp_mppc->AddProperty("EFFICIENCY", mppc_ep, mppc_effi,2);
  surface_mppc->SetMaterialPropertiesTable(sp_mppc);

  new G4LogicalSkinSurface("mppc_surface",MPPCLW,surface_mppc);
  */  


 


  return physWorld;
}		    


void BACDetectorConstruction::ConstructSDandField()
{
  
  auto aeroSD = new AeroSD("aeroSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(aeroSD);
  AeroLW->SetSensitiveDetector(aeroSD);
  //ACLW->SetSensitiveDetector(aeroSD);





  auto mppcSD = new MPPCSD("mppcSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(mppcSD);
  MPPCLW->SetSensitiveDetector(mppcSD);


  //MPPCLW1->SetSensitiveDetector(mppcSD);

  auto pmt15SD = new PMT15SD("pmt15SD");
  G4SDManager::GetSDMpointer()->AddNewDetector(pmt15SD);
  //PMT15LW->SetSensitiveDetector(pmt15SD);
  /*
  MPPCLW->SetSensitiveDetector(pmt15SD);
  //MPPC1LW->SetSensitiveDetector(pmt15SD);
  */

  /*
  auto pmt28SD = new PMT28SD("pmt28SD");
  G4SDManager::GetSDMpointer()->AddNewDetector(pmt28SD);
  PMT28LW->SetSensitiveDetector(pmt28SD);
  */





}



  
