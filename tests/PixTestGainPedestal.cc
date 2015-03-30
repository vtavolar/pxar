#include <stdlib.h>     /* atof, atoi */
#include <algorithm>    // std::find
#include <iostream>
#include <fstream>

#include <TH1.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TStopwatch.h>

#include "PixTestGainPedestal.hh"
#include "PixUtil.hh"
#include "log.h"


using namespace std;
using namespace pxar;

ClassImp(PixTestGainPedestal)

// ----------------------------------------------------------------------
PixTestGainPedestal::PixTestGainPedestal(PixSetup *a, std::string name) : PixTest(a, name), 
  fParNtrig(-1), fParShowFits(0), fParExtended(0), fParDumpHists(0)  {
  PixTest::init();
  init(); 
}


//----------------------------------------------------------
PixTestGainPedestal::PixTestGainPedestal() : PixTest() {
  //  LOG(logDEBUG) << "PixTestGainPedestal ctor()";
}

// ----------------------------------------------------------------------
bool PixTestGainPedestal::setParameter(string parName, string sval) {
  bool found(false);
  std::transform(parName.begin(), parName.end(), parName.begin(), ::tolower);
  for (unsigned int i = 0; i < fParameters.size(); ++i) {
    if (fParameters[i].first == parName) {
      found = true; 
      sval.erase(remove(sval.begin(), sval.end(), ' '), sval.end());
      if (!parName.compare("ntrig")) {
	fParNtrig = atoi(sval.c_str()); 
      }
      if (!parName.compare("showfits")) {
	PixUtil::replaceAll(sval, "checkbox(", "");
	PixUtil::replaceAll(sval, ")", "");
	fParShowFits = atoi(sval.c_str()); 
      }
      if (!parName.compare("extended")) {
	PixUtil::replaceAll(sval, "checkbox(", "");
	PixUtil::replaceAll(sval, ")", "");
	fParExtended = atoi(sval.c_str()); 
      }
      if (!parName.compare("dumphists")) {
	PixUtil::replaceAll(sval, "checkbox(", "");
	PixUtil::replaceAll(sval, ")", "");
	fParDumpHists = atoi(sval.c_str()); 
      }
      setToolTips();
      break;
    }
  }
  
  return found; 
}


// ----------------------------------------------------------------------
void PixTestGainPedestal::setToolTips() {
  fTestTip    = string(Form("measure and fit pulseheight vs VCAL (combining low- and high-range)\n")); 
  fSummaryTip = string("all ROCs are displayed side-by-side. Note the orientation:")
    + string("\nthe canvas bottom corresponds to the narrow module side with the cable")
    ;
}

// ----------------------------------------------------------------------
void PixTestGainPedestal::init() {

  setToolTips(); 

  fDirectory = gFile->GetDirectory(fName.c_str()); 
  if (!fDirectory) {
    fDirectory = gFile->mkdir(fName.c_str()); 
  } 
  fDirectory->cd(); 

}

// ----------------------------------------------------------------------
void PixTestGainPedestal::bookHist(string /*name*/) {
  fDirectory->cd(); 
  //  fHistList.clear();

}


//----------------------------------------------------------
PixTestGainPedestal::~PixTestGainPedestal() {
  LOG(logDEBUG) << "PixTestGainPedestal dtor";
}


// ----------------------------------------------------------------------
void PixTestGainPedestal::doTest() {

  TStopwatch t;

  fDirectory->cd();
  PixTest::update(); 
  bigBanner(Form("PixTestGainPedestal::doTest() ntrig = %d", fParNtrig));

  measure();
  fit();
  saveGainPedestalParameters();

  int seconds = t.RealTime(); 
  LOG(logINFO) << "PixTestGainPedestal::doTest() done, duration: " << seconds << " seconds";
}


// ----------------------------------------------------------------------
void PixTestGainPedestal::fullTest() {

  TStopwatch t;

  fDirectory->cd();
  PixTest::update(); 
  bigBanner(Form("PixTestGainPedestal::fullTest() ntrig = %d", fParNtrig));

  //  fParDumpHists = 1; 

  measure();
  fit();
  saveGainPedestalParameters();

  int seconds = t.RealTime(); 
  LOG(logINFO) << "PixTestGainPedestal::doTest() done, duration: " << seconds << " seconds";
}


// ----------------------------------------------------------------------
void PixTestGainPedestal::runCommand(string command) {
  std::transform(command.begin(), command.end(), command.begin(), ::tolower);
  LOG(logDEBUG) << "running command: " << command;
  if (!command.compare("measure")) {
    measure(); 
    return;
  }
  if (!command.compare("fit")) {
    fit(); 
    return;
  }
  if (!command.compare("save")) {
    saveGainPedestalParameters(); 
    return;
  }
  return;
}



// ----------------------------------------------------------------------
void PixTestGainPedestal::measure() {
  uint16_t FLAGS = FLAG_FORCE_MASKED;
  LOG(logDEBUG) << " using FLAGS = "  << (int)FLAGS; 


  fLpoints.clear();
  fLpoints.push_back(50); 
  fLpoints.push_back(100); 
  fLpoints.push_back(150); 
  fLpoints.push_back(200); 
  fLpoints.push_back(250); 

  fHpoints.clear();
  if (1 == fParExtended) {
    fHpoints.push_back(10); //new:  70
    fHpoints.push_back(17); //new: 119
    fHpoints.push_back(24); //new: 168
  }
  fHpoints.push_back(30); 
  fHpoints.push_back(50); 
  fHpoints.push_back(70); 
  //  fHpoints.push_back(90); 
  if (1 == fParExtended) {
    fHpoints.push_back(120); //new
  }
  //  fHpoints.push_back(200); 


  cacheDacs();
 
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 

  shist256 *pshistBlock  = new (fPixSetup->fPxarMemory) shist256[16*52*80]; 
  shist256 *ph;
  
  int idx(0);
  fHists.clear();
  for (unsigned int iroc = 0; iroc < rocIds.size(); ++iroc) {
    for (unsigned int ic = 0; ic < 52; ++ic) {
      for (unsigned int ir = 0; ir < 80; ++ir) {
	idx = PixUtil::rcr2idx(iroc, ic, ir); 
	ph = pshistBlock + idx;
	fHists.push_back(ph); 
      }
    }
  }

  fApi->_dut->testAllPixels(true);
  fApi->_dut->maskAllPixels(false);

  // -- first low range 
  fApi->setDAC("ctrlreg", 0);

  vector<pair<uint8_t, vector<pixel> > > rresult, lresult, hresult; 
  for (unsigned int i = 0; i < fLpoints.size(); ++i) {
    LOG(logINFO) << "scanning low vcal = " << fLpoints[i];
    int cnt(0); 
    bool done = false;
    while (!done){
      try {
	rresult = fApi->getPulseheightVsDAC("vcal", fLpoints[i], fLpoints[i], FLAGS, fParNtrig);
	copy(rresult.begin(), rresult.end(), back_inserter(lresult)); 
	done = true; // got our data successfully
      }
      catch(pxar::DataMissingEvent &e){
	LOG(logCRITICAL) << "problem with readout: "<< e.what() << " missing " << e.numberMissing << " events"; 
	++cnt;
	if (e.numberMissing > 10) done = true; 
      } catch(pxarException &e) {
	LOG(logCRITICAL) << "pXar execption: "<< e.what(); 
	++cnt;
      }
      done = (cnt>2) || done;
    }
  }

  // -- and high range
  fApi->setDAC("ctrlreg", 4);
  for (unsigned int i = 0; i < fHpoints.size(); ++i) {
    LOG(logINFO) << "scanning high vcal = " << fHpoints[i] << " (= " << 7*fHpoints[i] << " in low range)";
    int cnt(0); 
    bool done = false;
    while (!done){
      try {
	rresult = fApi->getPulseheightVsDAC("vcal", fHpoints[i], fHpoints[i], FLAGS, fParNtrig);
	copy(rresult.begin(), rresult.end(), back_inserter(hresult)); 
	done = true; // got our data successfully
      }
      catch(pxar::DataMissingEvent &e){
	LOG(logCRITICAL) << "problem with readout: "<< e.what() << " missing " << e.numberMissing << " events"; 
	++cnt;
	if (e.numberMissing > 10) done = true; 
      } catch(pxarException &e) {
	LOG(logCRITICAL) << "pXar execption: "<< e.what(); 
	++cnt;
      }
      done = (cnt>2) || done;
    }
  }

  for (unsigned int i = 0; i < lresult.size(); ++i) {
    int dac = lresult[i].first; 
    int dacbin(-1); 
    for (unsigned int v = 0; v < fLpoints.size(); ++v) {
      if (fLpoints[v] == dac) {
	dacbin = v; 
	break;
      }
    }
    vector<pixel> vpix = lresult[i].second;
    for (unsigned int ipx = 0; ipx < vpix.size(); ++ipx) {
      int roc = vpix[ipx].roc();
      int ic = vpix[ipx].column();
      int ir = vpix[ipx].row();

      double val = vpix[ipx].value();
      int idx = PixUtil::rcr2idx(getIdxFromId(roc), ic, ir);
      if (idx > -1) fHists[idx]->fill(dacbin+1, val);

    } 
  } 

  for (unsigned int i = 0; i < hresult.size(); ++i) {
    int dac = hresult[i].first; 
    int dacbin(-1); 
    for (unsigned int v = 0; v < fHpoints.size(); ++v) {
      if (fHpoints[v] == dac) {
	dacbin = 100 + v; 
	break;
      }
    }
    vector<pixel> vpix = hresult[i].second;
    for (unsigned int ipx = 0; ipx < vpix.size(); ++ipx) {
      int roc = vpix[ipx].roc();
      int ic = vpix[ipx].column();
      int ir = vpix[ipx].row();
      double val = vpix[ipx].value();
      int idx = PixUtil::rcr2idx(getIdxFromId(roc), ic, ir);
      if (idx > -1) fHists[idx]->fill(dacbin+1, val);

    } 
  } 

  printHistograms();

  PixTest::update(); 
  restoreDacs();
  LOG(logINFO) << "PixTestGainPedestal::measure() done ";
}




// ----------------------------------------------------------------------
void PixTestGainPedestal::fit() {
  PixTest::update(); 
  fDirectory->cd();

  string name = Form("gainPedestal_c%d_r%d_C%d", 0, 0, 0); 
  TH1D *h1 = bookTH1D(name, name, 1800, 0., 1800.);
  h1->SetMinimum(0);
  h1->SetMaximum(260.); 
  h1->SetNdivisions(506);
  h1->SetMarkerStyle(20);
  h1->SetMarkerSize(1.);

  TF1 *f(0); 

  //  vector<vector<gainPedestalParameters> > v;
  vector<vector<gainPedestalParameters> > v;
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
  gainPedestalParameters a; 
  //  a.p0 = a.p1 = a.p2 = a.p3 = 0.;
  a.p0 = a.p1 = a.p2  = 0.;
  TH1D* h(0);
  vector<TH1D*> p2list; 
  double fracErr(0.05); 
  for (unsigned int i = 0; i < rocIds.size(); ++i) {
    LOG(logDEBUG) << "Create hist " << Form("gainPedestalP1_C%d", rocIds[i]); 
    h = bookTH1D(Form("gainPedestalP1_C%d", i), Form("gainPedestalP1_C%d", rocIds[i]), 100, 0., 2.); 
    setTitles(h, "p1", "Entries / Bin"); 
    p2list.push_back(h); 
    vector<gainPedestalParameters> vroc;
    for (unsigned j = 0; j < 4160; ++j) {
      vroc.push_back(a); 
    }
    v.push_back(vroc); 
  }
  
  int iroc(0), ic(0), ir(0); 
  for (unsigned int i = 0; i < fHists.size(); ++i) {

    h1->Reset();
    for (int ib = 0; ib < static_cast<int>(fLpoints.size()); ++ib) {
      h1->SetBinContent(fLpoints[ib]+1, fHists[i]->get(ib+1));
      h1->SetBinError(fLpoints[ib]+1, 1.); 
    }
    for (int ib = 0; ib < static_cast<int>(fHpoints.size()); ++ib) {
      h1->SetBinContent(7*fHpoints[ib]+1, fHists[i]->get(100+ib+1));
      h1->SetBinError(7*fHpoints[ib]+1, 1.); 
    }
    
    //f = fPIF->gpTanH(h1); 
    f = fPIF->gpPol2(h1);


    if (h1->Integral() < 1) continue;
    PixUtil::idx2rcr(i, iroc, ic, ir);
    if (fParShowFits) {
      TH1D *hc = (TH1D*)h1->Clone(Form("gainPedestal_c%d_r%d_C%d", ic, ir, iroc));
      hc->SetTitle(Form("gainPedestal_c%d_r%d_C%d", ic, ir, iroc)); 
      string hcname = hc->GetName();
      LOG(logDEBUG) << hcname; 
      hc->Fit(f, "r");
      fHistList.push_back(hc); 
      PixTest::update(); 
    } else {
      //      cout << Form("gainPedestal_c%d_r%d_C%d", ic, ir, iroc) << endl;
      if (fParDumpHists) {
	h1->SetTitle(Form("gainPedestal_c%d_r%d_C%d", ic, ir, iroc)); 
	h1->SetName(Form("gainPedestal_c%d_r%d_C%d", ic, ir, iroc)); 
      }
      h1->Fit(f, "rq");
      if (fParDumpHists) {
	h1->SetDirectory(fDirectory); 
	h1->Write();
      }
    }
    int idx = ic*80 + ir; 
    v[iroc][idx].p0 = f->GetParameter(0); 
    v[iroc][idx].p1 = f->GetParameter(1); 
    v[iroc][idx].p2 = f->GetParameter(2); 
    //    v[iroc][idx].p3 = f->GetParameter(3); 

    p2list[getIdxFromId(iroc)]->Fill(f->GetParameter(2)); 
  }

  fPixSetup->getConfigParameters()->setGainPedestalParameters(v);

  copy(p2list.begin(), p2list.end(), back_inserter(fHistList));
  h = (TH1D*)(fHistList.back());
  h->Draw();

  string p2MeanString(""), p2RmsString(""); 
  for (unsigned int i = 0; i < p2list.size(); ++i) {
    p2MeanString += Form(" %5.3f", p2list[i]->GetMean()); 
    p2RmsString += Form(" %5.3f", p2list[i]->GetRMS()); 
  }

  fDisplayedHist = find(fHistList.begin(), fHistList.end(), h);
  PixTest::update(); 

  LOG(logINFO) << "PixTestGainPedestal::fit() done"; 
  LOG(logINFO) << "p2 mean: " << p2MeanString; 
  LOG(logINFO) << "p2 RMS:  " << p2RmsString; 
}



// ----------------------------------------------------------------------
void PixTestGainPedestal::saveGainPedestalParameters() {
  fPixSetup->getConfigParameters()->writeGainPedestalParameters();
}


// ----------------------------------------------------------------------
void PixTestGainPedestal::printHistograms() {

  ofstream OutputFile;
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
  int nRocs = static_cast<int>(rocIds.size()); 

  for (int iroc = 0; iroc < nRocs; ++iroc) {

    OutputFile.open(Form("%s/%s_C%d.dat", fPixSetup->getConfigParameters()->getDirectory().c_str(), 
			 fPixSetup->getConfigParameters()->getGainPedestalFileName().c_str(), 
			 iroc));

    OutputFile << "Pulse heights for the following Vcal values:" << endl;
    OutputFile << "Low range:  ";
    for (unsigned int i = 0; i < fLpoints.size(); ++i) OutputFile << fLpoints[i] << " "; 
    OutputFile << endl;
    OutputFile << "High range:  ";
    for (unsigned int i = 0; i < fHpoints.size(); ++i) OutputFile << fHpoints[i] << " "; 
    OutputFile << endl;
    OutputFile << endl;

    int roc(0), ic(0), ir(0); 
    string line(""); 
    for (int i = iroc*4160; i < (iroc+1)*4160; ++i) {
      PixUtil::idx2rcr(i, roc, ic, ir);
      if (roc != iroc) {
	LOG(logDEBUG) << "BIG CONFUSION?! iroc = " << iroc << " roc = " << roc;
      }

      line.clear();
      for (int ib = 0; ib < static_cast<int>(fLpoints.size()); ++ib) {
	line += Form(" %3d", static_cast<int>(fHists[i]->get(ib+1))); 
      }
      
      for (int ib = 0; ib < static_cast<int>(fHpoints.size()); ++ib) {
	line += Form(" %3d", static_cast<int>(fHists[i]->get(100+ib+1))); 
      }
      
      line += Form("    Pix %2d %2d", ic, ir); 
      OutputFile << line << endl;
    }

    OutputFile.close();
  }
}
