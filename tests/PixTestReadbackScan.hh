#ifndef PIXTESTREADBACKSCAN_H
#define PIXTESTREADBACKSCAN_H

#include "PixTest.hh"
#include "PHCalibration.hh"

#include <TH1.h>
#include <TH2.h>
#include <TProfile2D.h>

#include <TTree.h>

class DLLEXPORT PixTestReadbackScan: public PixTest {
public:
  PixTestReadbackScan(PixSetup *, std::string);
  PixTestReadbackScan();
  virtual ~PixTestReadbackScan();
  void init(); 
  void setToolTips();
  void bookHist(std::string);  
  void runCommand(std::string command);
  virtual bool setParameter(std::string parName, std::string sval);
  bool setTrgFrequency(uint8_t TrgTkDel);
  void pgToDefault();
  void setHistos();
  void ProcessData(uint16_t numevents = 1000);
  void FinalCleaning();
  uint8_t daqReadback(std::string dac, uint8_t vana, int8_t parReadback);
  uint8_t daqReadback(std::string dac, double vana, int8_t parReadback);
  void CalibrateIa();
  void CalibrateVana();
  void CalibrateVd();
  void CalibrateVa();
  void getCalibratedVbg();
  void doTest();

private:

  void stop();


  uint8_t  fParReadback;
  uint16_t fParPeriod;
  uint16_t fParStretch; 
  bool     fParFillTree;
  uint16_t fParTriggerFrequency;
  bool	   fParResetROC;
  
  bool     fPhCalOK;
  PHCalibration fPhCal;
  bool	   fParOutOfRange;
  bool     fDaq_loop;

  double fPar0VdCal;
  double fPar1VdCal;
  
  std::vector<std::pair<std::string, uint8_t> > fPg_setup;

  std::vector<TH2D*> fHits;
  std::vector<TProfile2D*> fPhmap;
  std::vector<TH1D*> fPh;
  std::vector<TH1D*> fQ;
  std::vector<TProfile2D*> fQmap;

  ClassDef(PixTestReadbackScan, 1)

};
#endif
