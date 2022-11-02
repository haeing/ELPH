#ifndef BACAnalysisManager_h
#define BACAnalysisManager_h

#include "globals.hh"
#include "G4ThreeVector.hh"

#include <cmath>
#include <vector>

class G4Run;
class G4Event;
class TFile;
class TTree;
class EvtVector4R;

const int num_evtgen = 100000;
const int num_mppchit = 100000;
const int num_aerohit = 100000;
const int num_pmt15hit = 100000;
const int num_pmt28hit = 100000;

class BACAnalysisManager
{
public:
  void BeginOfRun(const G4Run*);
  void EndOfRun(const G4Run*);
  void BeginOfEvent(const G4Event*);
  void EndOfEvent(const G4Event*);

private:
  G4String outfile;
  TFile *hfile;
  TTree *tree;

public:
  BACAnalysisManager(const G4String &histname);
  virtual ~BACAnalysisManager();
  BACAnalysisManager(const BACAnalysisManager&);
  BACAnalysisManager& operator=(const BACAnalysisManager&);

  void SetEvtGen(int j, int partnum);
  void Terminate(void) const;
  void SaveFile(void) const;

private:
  G4bool fActive_;
  int event;

  int nEvt;
  int evtid[num_evtgen];
  //int evtpid[num_evtgen];
  int evtpid;
  double evtposx;
  double evtposy;
  double evtposz;
  int evtnumce;
  

  int nhMppc;
  //int mppcmulti;
  int mppcpid[num_mppchit];
  double mppctime[num_mppchit];
  double mppcposx[num_mppchit];
  double mppcposy[num_mppchit];
  double mppcposz[num_mppchit];
  double mppcwavelength[num_mppchit];
  int mppcnum[num_mppchit];

  int nhAero;
  double aeroangle[num_aerohit];
  int aeropid[num_aerohit];
  double aerotime[num_aerohit];
  double aeroposx[num_aerohit];
  double aeroposy[num_aerohit];
  double aeroposz[num_aerohit];

  int nhPmt15;
  //int mppcmulti;
  int pmt15pid[num_pmt15hit];
  double pmt15time[num_pmt15hit];
  double pmt15posx[num_pmt15hit];
  double pmt15posy[num_pmt15hit];
  double pmt15posz[num_pmt15hit];
  double pmt15wavelength[num_pmt15hit];
  int pmt15num[num_pmt15hit];

  int nhPmt28;
  //int mppcmulti;
  int pmt28pid[num_pmt28hit];
  double pmt28time[num_pmt28hit];
  double pmt28posx[num_pmt28hit];
  double pmt28posy[num_pmt28hit];
  double pmt28posz[num_pmt28hit];
  double pmt28wavelength[num_pmt28hit];
  int pmt28num[num_pmt28hit];

  
  
};

#endif
