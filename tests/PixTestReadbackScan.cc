#include <stdlib.h>     /* atof, atoi */
#include <algorithm>    // std::find
#include <iostream>
#include "PixTestReadbackScan.hh"
#include "log.h"
#include "helper.h"
#include "timer.h"


using namespace std;
using namespace pxar;

ClassImp(PixTestReadbackScan)

// ----------------------------------------------------------------------
PixTestReadbackScan::PixTestReadbackScan(PixSetup *a, std::string name) : PixTest(a, name), fParStretch(0), fParTriggerFrequency(0), fParResetROC(0) {
  PixTest::init();
  init(); 
  LOG(logDEBUG) << "PixTestReadbackScan ctor(PixSetup &a, string, TGTab *)";
  fTree = 0; 

  fPhCal.setPHParameters(fPixSetup->getConfigParameters()->getGainPedestalParameters());
  fPhCalOK = fPhCal.initialized();
}


//----------------------------------------------------------
PixTestReadbackScan::PixTestReadbackScan() : PixTest() {
  LOG(logDEBUG) << "PixTestReadbackScan ctor()";
  fTree = 0; 
}

//----------------------------------------------------------
PixTestReadbackScan::~PixTestReadbackScan() {
	LOG(logDEBUG) << "PixTestReadbackScan dtor, saving tree ... ";
	fDirectory->cd();
	if (fTree && fParFillTree) fTree->Write();
}

// ----------------------------------------------------------------------
void PixTestReadbackScan::init() {
  LOG(logDEBUG) << "PixTestReadbackScan::init()";

  setToolTips();
  fDirectory = gFile->GetDirectory(fName.c_str()); 
  if (!fDirectory) {
    fDirectory = gFile->mkdir(fName.c_str()); 
  } 
  fDirectory->cd(); 
}

// ----------------------------------------------------------------------
void PixTestReadbackScan::setToolTips() {
  fTestTip    = string("Run DAQ - data from each run will be added to the same histogram.") ;
  fSummaryTip = string("to be implemented") ;
  fStopTip    = string("Stop DAQ and save data.");
}

// ----------------------------------------------------------------------
void PixTestReadbackScan::bookHist(string name) {
	fDirectory->cd();
	LOG(logDEBUG) << "nothing done with " << name;
}

// ----------------------------------------------------------------------
void PixTestReadbackScan::stop(){
	// Interrupt the test 
	fDaq_loop = false;
	LOG(logINFO) << "Stop pressed. Ending test.";
}

// ----------------------------------------------------------------------
void PixTestReadbackScan::runCommand(std::string command) {

	if (command == "stop")
		stop();
	else
		LOG(logINFO) << "Command " << command << " not implemented.";
}

// ----------------------------------------------------------------------
bool PixTestReadbackScan::setParameter(string parName, string sval) {
	bool found(false);
	fParOutOfRange = false;
	std::transform(parName.begin(), parName.end(), parName.begin(), ::tolower);
	for (unsigned int i = 0; i < fParameters.size(); ++i) {
		if (fParameters[i].first == parName) {
			found = true;
			if (!parName.compare("readback")) {
				fParReadback = atoi(sval.c_str());
				setToolTips();
			}
			if (!parName.compare("clockstretch")) {
				fParStretch = atoi(sval.c_str());
				setToolTips();
			}
			if (!parName.compare("filltree")) {
				fParFillTree = !(atoi(sval.c_str()) == 0);
				setToolTips();
			}
			if (!parName.compare("trgfrequency(khz)")){   // trigger frequency in kHz.
				fParTriggerFrequency = atoi(sval.c_str());
				if (fParTriggerFrequency == 0) {
					LOG(logWARNING) << "PixTestReadbackScan::setParameter() trgfrequency must be different from zero";
					found = false; fParOutOfRange = true;
				}
			}
		}
	}
	return found;
}

//----------------------------------------------------------
bool PixTestReadbackScan::setTrgFrequency(uint8_t TrgTkDel){

	int nDel = 0;
	uint8_t trgtkdel= TrgTkDel;
	double period_ns = 1 / (double)fParTriggerFrequency * 1000000; // trigger frequency in kHz.
	fParPeriod = (uint16_t)period_ns / 25;
	uint16_t ClkDelays = fParPeriod - trgtkdel;
	
	//add right delay between triggers:
	
	if (fParResetROC) {       //by default not reset (already done before daqstart)
		fPg_setup.push_back(make_pair("resetroc", 15));
		ClkDelays -= 15;
		nDel++;
	}

	while (ClkDelays>255){
		fPg_setup.push_back(make_pair("delay", 255));
		ClkDelays = ClkDelays - 255;
		nDel ++;
	}
	fPg_setup.push_back(make_pair("delay", ClkDelays));

	//then send trigger and token:
	fPg_setup.push_back(make_pair("trg", trgtkdel));
	fPg_setup.push_back(make_pair("tok", 0));
	
	fParPeriod = fParPeriod + 4 + nDel; //to align to the new pg minimum (1 additional clk cycle per PG call);

	return true;
}

// ----------------------------------------------------------------------
void PixTestReadbackScan::pgToDefault() {
  fPg_setup.clear();
  LOG(logDEBUG) << "PixTestPattern::PG_Setup clean";
  
  fPg_setup = fPixSetup->getConfigParameters()->getTbPgSettings();
  fApi->setPatternGenerator(fPg_setup);
  LOG(logINFO) << "PixTestPattern::       pg_setup set to default.";
}

// ----------------------------------------------------------------------
void PixTestReadbackScan::setHistos(){
	
	if (fParFillTree) bookTree();
	fHits.clear(); fPhmap.clear(); fPh.clear(); fQmap.clear(); fQ.clear();

	std::vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs();
	TH1D *h1(0);
	TH2D *h2(0);
	TProfile2D *p2(0);
	for (unsigned int iroc = 0; iroc < rocIds.size(); ++iroc){
		h2 = bookTH2D(Form("hits_C%d", rocIds[iroc]), Form("hits_C%d", rocIds[iroc]), 52, 0., 52., 80, 0., 80.);
		h2->SetMinimum(0.);
		h2->SetDirectory(fDirectory);
		setTitles(h2, "col", "row");
		fHistOptions.insert(make_pair(h2, "colz"));
		fHits.push_back(h2);

		p2 = bookTProfile2D(Form("phMap_C%d", rocIds[iroc]), Form("phMap_C%d", rocIds[iroc]), 52, 0., 52., 80, 0., 80.);
		p2->SetMinimum(0.);
		p2->SetDirectory(fDirectory);
		setTitles(p2, "col", "row");
		fHistOptions.insert(make_pair(p2, "colz"));
		fPhmap.push_back(p2);

		h1 = bookTH1D(Form("ph_C%d", rocIds[iroc]), Form("ph_C%d", rocIds[iroc]), 256, 0., 256.);
		h1->SetMinimum(0.);
		h1->SetDirectory(fDirectory);
		setTitles(h1, "ADC", "Entries/bin");
		fPh.push_back(h1);

		p2 = bookTProfile2D(Form("qMap_C%d", rocIds[iroc]), Form("qMap_C%d", rocIds[iroc]), 52, 0., 52., 80, 0., 80.);
		p2->SetMinimum(0.);
		p2->SetDirectory(fDirectory);
		setTitles(p2, "col", "row");
		fHistOptions.insert(make_pair(p2, "colz"));
		fQmap.push_back(p2);

		h1 = bookTH1D(Form("q_C%d", rocIds[iroc]), Form("q_C%d", rocIds[iroc]), 200, 0., 1000.);
		h1->SetMinimum(0.);
		h1->SetDirectory(fDirectory);
		setTitles(h1, "Q [Vcal]", "Entries/bin");
		fQ.push_back(h1);
	}
}


// ----------------------------------------------------------------------
void PixTestReadbackScan::ProcessData(uint16_t numevents){

	LOG(logDEBUG) << "Getting Event Buffer";
	std::vector<pxar::Event> daqdat;

	if (numevents > 0) {
		for (unsigned int i = 0; i < numevents; i++) {
			pxar::Event evt = fApi->daqGetEvent();
			//Check if event is empty?
			if (evt.pixels.size() > 0)
				daqdat.push_back(evt);
		}
	}
	else
		daqdat = fApi->daqGetEventBuffer();

	LOG(logDEBUG) << "Processing Data: " << daqdat.size() << " events.";

	int pixCnt(0);
	int idx(-1);
	uint16_t q;
	vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs();
	for (std::vector<pxar::Event>::iterator it = daqdat.begin(); it != daqdat.end(); ++it) {
		pixCnt += it->pixels.size();

		if (fParFillTree) {
			fTreeEvent.header = it->header;
			fTreeEvent.dac = 0;
			fTreeEvent.trailer = it->trailer;
			fTreeEvent.npix = it->pixels.size();
		}

		for (unsigned int ipix = 0; ipix < it->pixels.size(); ++ipix) {
			idx = getIdxFromId(it->pixels[ipix].roc());
			if(idx == -1) {
				LOG(logWARNING) << "PixTestReadbackScan::ProcessData() wrong 'idx' value --> return";
				return;    			
			}
			fHits[idx]->Fill(it->pixels[ipix].column(), it->pixels[ipix].row());
			fPhmap[idx]->Fill(it->pixels[ipix].column(), it->pixels[ipix].row(), it->pixels[ipix].value());
			fPh[idx]->Fill(it->pixels[ipix].value());

			if (fPhCalOK) {
				q = static_cast<uint16_t>(fPhCal.vcal(it->pixels[ipix].roc(), it->pixels[ipix].column(),	
								      it->pixels[ipix].row(), it->pixels[ipix].value()));
			}
			else {
				q = 0;
			}
			fQ[idx]->Fill(q);
			fQmap[idx]->Fill(it->pixels[ipix].column(), it->pixels[ipix].row(), q);
				if (fParFillTree) {
				fTreeEvent.proc[ipix] = it->pixels[ipix].roc();
				fTreeEvent.pcol[ipix] = it->pixels[ipix].column();
				fTreeEvent.prow[ipix] = it->pixels[ipix].row();
				fTreeEvent.pval[ipix] = it->pixels[ipix].value();
				fTreeEvent.pq[ipix] = q;
			}
		}
		if (fParFillTree) fTree->Fill();
	}

  	//to draw the hitsmap as 'online' check.
	TH2D* h2 = (TH2D*)(fHits.back());
	h2->Draw(getHistOption(h2).c_str());
	fDisplayedHist = find(fHistList.begin(), fHistList.end(), h2);
	PixTest::update();

	LOG(logINFO) << Form("events read: %6ld, pixels seen: %3d, hist entries: %4d",
	                 daqdat.size(), pixCnt,	static_cast<int>(fHits[0]->GetEntries()));	
}

// ----------------------------------------------------------------------
void PixTestReadbackScan::FinalCleaning() {

	// Reset the pg_setup to default value.
	pgToDefault();
	//clean local variables:
	fPg_setup.clear();
}

// ----------------------------------------------------------------------
void PixTestReadbackScan::doTest() {
  LOG(logINFO) << "PixTestReadbackScan::doTest() start.";
  int readback=0;
 
  fParReadback=10;
 
 // TH1D* h1= new TH1D("rbScan","rbScan", 255, 0., 255.);
 // for(uint8_t vana=0; vana<255; vana++){
 //   readback=daqReadback("vana", vana, fParReadback);
 //   h1->Fill(vana, readback);
 //   
 // }

  CalibrateVd();
  //  CalibrateVa();
  getCalibratedVbg();
 //fHistList.push_back(h1);
  for (list<TH1*>::iterator il = fHistList.begin(); il != fHistList.end(); ++il) {
    (*il)->Draw((getHistOption(*il)).c_str()); 
  }


 LOG(logINFO) << "PixTestReadbackScan::doTest() done";

}


void PixTestReadbackScan::CalibrateIa(){

  //readback DAC set to 12 (i.e. Ia)
  fParReadback=12;

  int readback=0;

  TH1D* h_rbIa = new TH1D("rbIa","rbIa", 256, 0., 256.);
  TH1D* h_tbIa = new TH1D("tbIa","tbIa", 256, 0., 256.);
  double tbIa = 0.;
  vector<double > rbIa;
  
  for(uint8_t vana=0; vana<255; vana++){
    readback=daqReadback("vana", vana, fParReadback);
    rbIa.push_back(readback);
    fApi->setDAC("vana", vana);
    tbIa = fApi->getTBia()*1E3;
    h_rbIa->Fill(vana, readback);
    h_tbIa->Fill(vana, tbIa);
  }

  double rb_vanaMax=0.;
  double tb_vanaMax=0.;
  double diff=0.;
//  for(int ibin=1; ibin<=h_rbIa->GetNbinsX()-5; ibin++){
//    diff = h_rbIa->GetBinContent(ibin) -  h_rbIa->GetBinContent(ibin+5);
//    if(diff==0){
//      rb_vanaMax=h_rbIa->GetBinCenter(ibin+5);
//      break;
//    }
//  }
//
// for(int ibin=1; ibin<=h_tbIa->GetNbinsX()-5; ibin++){
//    diff = h_tbIa->GetBinContent(ibin) -  h_tbIa->GetBinContent(ibin+5);
//    if(diff==0){
//      tb_vanaMax=h_tbIa->GetBinCenter(ibin+5);
//      break;
//    }
//  }


 rb_vanaMax = h_rbIa->GetBinCenter(h_rbIa->FindFirstBinAbove(254));


LOG(logDEBUG)<<"Vana max for fit:"<<endl<<"rb: "<<rb_vanaMax<<endl<<"tb :"<<tb_vanaMax;

  TF1* frb = new TF1("lin_rb", "[0] + x*[1]", 0, rb_vanaMax);
  TF1* ftb = new TF1("lin_ftb", "[0] + x*[1]", 0, 255);

  h_rbIa->Fit(frb, "W", "", 0., rb_vanaMax);
  h_tbIa->Fit(ftb);

  LOG(logDEBUG)<<"Number of points for rb fit "<<frb->GetNumberFitPoints();

  double rbpar0=0., rbpar1=0., tbpar0=0., tbpar1=0.;
  rbpar0=frb->GetParameter(0);
  rbpar1=frb->GetParameter(1);
  tbpar0=ftb->GetParameter(0);
  tbpar1=ftb->GetParameter(1);

  TH1D* h_rbIaCal = new TH1D("rbIa","rbIa", 256, 0., 256.);
  for(int vana=0; vana<256; vana++){
    h_rbIaCal->Fill(vana, ((tbpar1/rbpar1)*(rbIa[vana]-rbpar0)+tbpar0));
  }

  h_rbIaCal->SetLineColor(kBlue);
  fHistOptions.insert(make_pair(h_rbIaCal,"same"));

  fHistList.push_back(h_rbIa);
  fHistList.push_back(h_rbIaCal);
  fHistList.push_back(h_tbIa);
//  h_rbIaCal->Draw();
//  h_tbIa->Draw("same");


}


void PixTestReadbackScan::CalibrateVana(){

  //readback DAC set to 11 (i.e. Vana)
  fParReadback=11;

  int readback=0;

  TH1D* h_rbVana = new TH1D("rbVana","rbVana", 256, 0., 256.);
  TH1D* h_dacVana = new TH1D("dacVana","dacVana", 256, 0., 256.);
  vector<double > rbVana;
  
  for(uint8_t vana=0; vana<255; vana++){
    readback=daqReadback("vana", vana, fParReadback);
    rbVana.push_back(readback);
    h_rbVana->Fill(vana, readback);
    h_dacVana->Fill(vana, vana);
  }

  double rb_vanaMax=0.;

  rb_vanaMax = h_rbVana->GetBinCenter(h_rbVana->FindFirstBinAbove(254));


  LOG(logDEBUG)<<"Vana max for fit:"<<endl<<"rb: "<<rb_vanaMax;

  TF1* frb = new TF1("lin_rb", "[0] + x*[1]", 0, rb_vanaMax);
  TF1* fdac = new TF1("lin_fdac", "[0] + x*[1]", 0, 255);

  h_rbVana->Fit(frb, "W", "", 0., rb_vanaMax);
  h_dacVana->Fit(fdac);

  LOG(logDEBUG)<<"Number of points for rb fit "<<frb->GetNumberFitPoints();

  double rbpar0=0., rbpar1=0., dacpar0=0., dacpar1=0.;
  rbpar0=frb->GetParameter(0);
  rbpar1=frb->GetParameter(1);
  dacpar0=fdac->GetParameter(0);
  dacpar1=fdac->GetParameter(1);

  TH1D* h_rbVanaCal = new TH1D("rbVana","rbVana", 256, 0., 256.);
  for(int vana=0; vana<256; vana++){
    h_rbVanaCal->Fill(vana, ((dacpar1/rbpar1)*(rbVana[vana]-rbpar0)+dacpar0));
  }

  h_rbVanaCal->SetLineColor(kBlue);
  fHistOptions.insert(make_pair(h_rbVanaCal,"same"));

  fHistList.push_back(h_rbVana);
  fHistList.push_back(h_rbVanaCal);
  fHistList.push_back(h_dacVana);
//  h_rbVanaCal->Draw();
//  h_dacVana->Draw("same");


}

void PixTestReadbackScan::CalibrateVd(){

  //readback DAC set to 8 (i.e. Vd)
  fParReadback=8;

  int readback=0;

  TH1D* h_rbVd = new TH1D("rbVd","rbVd", 500, 0., 5.);
  TH1D* h_dacVd = new TH1D("dacVd","dacVd", 500, 0., 5.);
  vector<double > rbVd;
  double Vd;

  for(int iVd=0; iVd<45; iVd++){
    LOG(logDEBUG)<<"/****:::::: CALIBRATE VD FUNCTION :::::****/";
    Vd = 2.1 + iVd*0.02;
    LOG(logDEBUG)<<"Digital voltage will be set to: "<<Vd;
    readback=daqReadback("vd", Vd, fParReadback);
    LOG(logDEBUG)<<"Voltage "<<Vd<<", readback "<<readback;
    rbVd.push_back(readback);
    h_rbVd->Fill(Vd, readback);
    h_dacVd->Fill(Vd, Vd);
  }

  // double rb_VdMax=0.;
  
  //rb_VdMax = h_rbVd->GetBinCenter(h_rbVd->FindFirstBinAbove(254));
  
  
  //  LOG(logDEBUG)<<"Vd max for fit:"<<endl<<"rb: "<<rb_VdMax;
  
  TF1* frb = new TF1("lin_vd", "[0] + x*[1]");
  

  h_rbVd->Fit(frb, "W", "");

  LOG(logDEBUG)<<"Number of points for rb fit "<<frb->GetNumberFitPoints();

  //  double rbpar0=0., rbpar1=0.;
  fPar0VdCal=frb->GetParameter(0);
  fPar1VdCal=frb->GetParameter(1);

  fHistList.push_back(h_rbVd);
  //fHistList.push_back(h_rbVdCal);
  fHistList.push_back(h_dacVd);
//  h_rbVdCal->Draw();
//  h_dacVd->Draw("same");


}


void PixTestReadbackScan::getCalibratedVbg(){

  //readback DAC set to 11 (i.e. Vbg)
  fParReadback=11;

  int readback=0;

  TH1D* h_rbVbg = new TH1D("rbVbg","rbVbg", 500, 0., 5.);
  TH1D* h_dacVbg = new TH1D("dacVbg","dacVbg", 500, 0., 5.);
  vector<double > rbVbg;
  double Vd;

  int n_meas=0.;
  double avReadback=0.;

  for(int i=0; i<10; i++){
    LOG(logDEBUG)<<"/****:::::: CALIBRATE VBG FUNCTION :::::****/";
    Vd = 2.5;
    LOG(logDEBUG)<<"Digital voltage will be set to: "<<Vd;
    readback = daqReadback("vd", Vd, fParReadback);
    if (0!=readback){
      avReadback+=(double)readback;
      n_meas++;
    }
    LOG(logDEBUG)<<"Voltage "<<Vd<<", average readback "<<(double)readback/(i+1);
    rbVbg.push_back(readback);
    h_rbVbg->Fill(Vd, readback);
    h_dacVbg->Fill(Vd, Vd);
  }


  avReadback/=n_meas;
 

  double calVbg=0;
  //0.5 needed because Vbg rb has twice the sensitivity Vd and Va have
  calVbg=0.5*(avReadback-fPar0VdCal)/fPar1VdCal;
  LOG(logDEBUG)<<"/*/*/*/*::: Calibrated Vbg = "<<calVbg<<" :::*/*/*/*/";
  //  double rbpar0=0., rbpar1=0.;
  // fPar0VdCal=frb->GetParameter(0);
  // fPar1VdCal=frb->GetParameter(1);

  fHistList.push_back(h_rbVbg);
  //fHistList.push_back(h_rbVdCal);
  fHistList.push_back(h_dacVbg);
//  h_rbVdCal->Draw();
//  h_dacVd->Draw("same");


}





void PixTestReadbackScan::CalibrateVa(){

  //readback DAC set to 9 (i.e. Va)
  fParReadback=11;

  int readback=0;

  TH1D* h_rbVa = new TH1D("rbVa","rbVa", 500, 0., 5.);
  TH1D* h_dacVa = new TH1D("dacVa","dacVa", 500, 0., 5.);
  vector<double > rbVa;
  double Va;

  for(int iVa=0; iVa<35; iVa++){
    LOG(logDEBUG)<<"/****:::::: CALIBRATE VA FUNCTION :::::****/";
    Va = 1.5 + iVa*0.02;
    LOG(logDEBUG)<<"Digital voltage will be set to: "<<Va;
    readback=daqReadback("va", Va, fParReadback);
    LOG(logDEBUG)<<"Voltage "<<Va<<", readback "<<readback;
    rbVa.push_back(readback);
    h_rbVa->Fill(Va, readback);
    h_dacVa->Fill(Va, Va);
  }

  //double rb_VaMax=0.;
  //
  //rb_VaMax = h_rbVa->GetBinCenter(h_rbVa->FindFirstBinAbove(254));
  //
  //
  //LOG(logDEBUG)<<"Va max for fit:"<<endl<<"rb: "<<rb_VaMax;
  //
  //TF1* frb = new TF1("lin_rb", "[0] + x*[1]", 0, rb_VaMax);
  //TF1* fdac = new TF1("lin_fdac", "[0] + x*[1]", 0, 255);

  //  h_rbVa->Fit(frb, "W", "", 0., rb_VaMax);
  //  h_dacVa->Fit(fdac);

  //LOG(logDEBUG)<<"Number of points for rb fit "<<frb->GetNumberFitPoints();

  //double rbpar0=0., rbpar1=0., dacpar0=0., dacpar1=0.;
  //rbpar0=frb->GetParameter(0);
  //rbpar1=frb->GetParameter(1);
  //dacpar0=fdac->GetParameter(0);
  //dacpar1=fdac->GetParameter(1);
  //
  //TH1D* h_rbVaCal = new TH1D("rbVa","rbVa", 256, 0., 256.);
  //for(int Va=0; Va<256; Va++){
  //  h_rbVaCal->Fill(Va, ((dacpar1/rbpar1)*(rbVa[Va]-rbpar0)+dacpar0));
  //}
  //
  //h_rbVaCal->SetLineColor(kBlue);
  //fHistOptions.insert(make_pair(h_rbVaCal,"same"));

  fHistList.push_back(h_rbVa);
  //fHistList.push_back(h_rbVaCal);
  fHistList.push_back(h_dacVa);
//  h_rbVaCal->Draw();
//  h_dacVa->Draw("same");


}




uint8_t PixTestReadbackScan::daqReadback(string dac, double vana, int8_t parReadback){


  PixTest::update();
  fDirectory->cd();
  fPg_setup.clear();

  if (!dac.compare("vana")){
    LOG(logDEBUG)<<"Wrong daqReadback function called!!!";
    //    fApi->setDAC(dac.c_str, (uint8_t)vana);
  }
  else if (!dac.compare("vd")){
    vector<pair<string,double > > powerset;
    powerset.push_back(make_pair("ia", 1.19));
    powerset.push_back(make_pair("id", 1.10));
    powerset.push_back(make_pair("va", fApi->getTBva()));
    powerset.push_back(make_pair("vd", vana));
    fApi->setTestboardPower(powerset);
    //    fApi->_hal->setTBvd(vana);
  }
else if (!dac.compare("va")){
    vector<pair<string,double > > powerset;
    powerset.push_back(make_pair("ia", 1.19));
    powerset.push_back(make_pair("id", 1.10));
    powerset.push_back(make_pair("va", vana));
    powerset.push_back(make_pair("vd", fApi->getTBvd()));
    fApi->setTestboardPower(powerset);
    //    fApi->_hal->setTBvd(vana);
  }

  fApi->setDAC("readback", parReadback);


  //Immediately stop if parameters not in range	
  if (fParOutOfRange) return 255;
  
 

  //Set the ClockStretch
  fApi->setClockStretch(0, 0, fParStretch); //Stretch after trigger, 0 delay
   
  //Set the histograms:
  if(fHistList.size() == 0) setHistos();  //to book histo only for the first 'doTest' (or after Clear).

  //To print on shell the number of masked pixels per ROC:
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs();
  LOG(logINFO) << "PixTestReadbackScan::Number of masked pixels:";
  for (unsigned int iroc = 0; iroc < rocIds.size(); ++iroc) {
	  LOG(logINFO) << "PixTestReadbackScan::    ROC " << static_cast<int>(iroc) << ": " << fApi->_dut->getNMaskedPixels(static_cast<int>(iroc));
  }  

  
  // Start the DAQ:
  //::::::::::::::::::::::::::::::::

  //:::Setting register to read back a given quantity::::://
  

  //First send only a RES:
  fPg_setup.push_back(make_pair("resetroc", 0));     // PG_RESR b001000 
  uint16_t period = 28;

  //Set the pattern generator:
  fApi->setPatternGenerator(fPg_setup);

  fApi->daqStart();

  //Send only one trigger to reset:
  fApi->daqTrigger(1, period);
  LOG(logINFO) << "PixTestReadbackScan::RES sent once ";

  fApi->daqStop();

  fPg_setup.clear();
  LOG(logINFO) << "PixTestReadbackScan::PG_Setup clean";

  //Set the pattern wrt the trigger frequency:
  LOG(logINFO) << "PG set to have trigger frequency = " << fParTriggerFrequency << " kHz";
  if (!setTrgFrequency(20)){
	  FinalCleaning();
	  return 255;
  }

  //Set pattern generator:
  fApi->setPatternGenerator(fPg_setup);

  fDaq_loop = true;

  //Start the DAQ:
  fApi->daqStart();

  int  Ntrig=32;
  //Send the triggers:
  fApi->daqTrigger(Ntrig, fParPeriod);
  gSystem->ProcessEvents();
  ProcessData(0);
 
  fApi->daqStop(); 

  std::vector<std::vector<uint16_t> > rb;
  rb = fApi->daqGetReadback();
  uint8_t rb_val=0;


  for(int i=0; i<rb.size(); i++){
    for(int j=0; j<rb[i].size(); j++){
      LOG(logDEBUG)<<"Readback values for vana = "<<(int)vana<<" : "<<(int)(rb[i][j]&0xff);
      rb_val=(rb[i][j]&0xff);
    }
  }
  

  
  //::::::::::::::::::::::::::::::
  //DAQ - THE END.

  //to draw and save histograms
 // TH1D *h1(0);
 // TH2D *h2(0);
 // TProfile2D *p2(0);
 // copy(fQ.begin(), fQ.end(), back_inserter(fHistList));
 // copy(fQmap.begin(), fQmap.end(), back_inserter(fHistList));
 // copy(fPh.begin(), fPh.end(), back_inserter(fHistList));
 // copy(fPhmap.begin(), fPhmap.end(), back_inserter(fHistList));
 // copy(fHits.begin(), fHits.end(), back_inserter(fHistList));
 // for (list<TH1*>::iterator il = fHistList.begin(); il != fHistList.end(); ++il) {
 //  	(*il)->Draw((getHistOption(*il)).c_str()); 
 // }
 // fDisplayedHist = find(fHistList.begin(), fHistList.end(), h1);
 // fDisplayedHist = find(fHistList.begin(), fHistList.end(), p2);
 // fDisplayedHist = find(fHistList.begin(), fHistList.end(), h1);
 // fDisplayedHist = find(fHistList.begin(), fHistList.end(), p2);
 // fDisplayedHist = find(fHistList.begin(), fHistList.end(), h2);
 // PixTest::update();

  FinalCleaning();
  fApi->setClockStretch(0, 0, 0); //No Stretch after trigger, 0 delay
  return rb_val;
 }


uint8_t PixTestReadbackScan::daqReadback(string dac, uint8_t vana, int8_t parReadback){


  PixTest::update();
  fDirectory->cd();
  fPg_setup.clear();

  if (!dac.compare("vana")){
    fApi->setDAC(dac.c_str(), vana);
  }
  else if (!dac.compare("vd")){
    LOG(logDEBUG)<<"Wrong daqReadback function called!!!";
    //fApi->_hal->setTBvd(vana);
  }

  fApi->setDAC("readback", parReadback);


  //Immediately stop if parameters not in range	
  if (fParOutOfRange) return 255;
  
 

  //Set the ClockStretch
  fApi->setClockStretch(0, 0, fParStretch); //Stretch after trigger, 0 delay
   
  //Set the histograms:
  if(fHistList.size() == 0) setHistos();  //to book histo only for the first 'doTest' (or after Clear).

  //To print on shell the number of masked pixels per ROC:
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs();
  LOG(logINFO) << "PixTestReadbackScan::Number of masked pixels:";
  for (unsigned int iroc = 0; iroc < rocIds.size(); ++iroc) {
	  LOG(logINFO) << "PixTestReadbackScan::    ROC " << static_cast<int>(iroc) << ": " << fApi->_dut->getNMaskedPixels(static_cast<int>(iroc));
  }  

  
  // Start the DAQ:
  //::::::::::::::::::::::::::::::::

  //:::Setting register to read back a given quantity::::://
  

  //First send only a RES:
  fPg_setup.push_back(make_pair("resetroc", 0));     // PG_RESR b001000 
  uint16_t period = 28;

  //Set the pattern generator:
  fApi->setPatternGenerator(fPg_setup);

  fApi->daqStart();

  //Send only one trigger to reset:
  fApi->daqTrigger(1, period);
  LOG(logINFO) << "PixTestReadbackScan::RES sent once ";

  fApi->daqStop();

  fPg_setup.clear();
  LOG(logINFO) << "PixTestReadbackScan::PG_Setup clean";

  //Set the pattern wrt the trigger frequency:
  LOG(logINFO) << "PG set to have trigger frequency = " << fParTriggerFrequency << " kHz";
  if (!setTrgFrequency(20)){
	  FinalCleaning();
	  return 255;
  }

  //Set pattern generator:
  fApi->setPatternGenerator(fPg_setup);

  fDaq_loop = true;

  //Start the DAQ:
  fApi->daqStart();

  int  Ntrig=32;
  //Send the triggers:
  fApi->daqTrigger(Ntrig, fParPeriod);
  gSystem->ProcessEvents();
  ProcessData(0);
 
  fApi->daqStop(); 

  std::vector<std::vector<uint16_t> > rb;
  rb = fApi->daqGetReadback();
  uint8_t rb_val=0;


  for(int i=0; i<rb.size(); i++){
    for(int j=0; j<rb[i].size(); j++){
      LOG(logDEBUG)<<"Readback values for vana = "<<(int)vana<<" : "<<(int)(rb[i][j]&0xff);
      rb_val=(rb[i][j]&0xff);
    }
  }
  

  
  //::::::::::::::::::::::::::::::
  //DAQ - THE END.

  //to draw and save histograms
 // TH1D *h1(0);
 // TH2D *h2(0);
 // TProfile2D *p2(0);
 // copy(fQ.begin(), fQ.end(), back_inserter(fHistList));
 // copy(fQmap.begin(), fQmap.end(), back_inserter(fHistList));
 // copy(fPh.begin(), fPh.end(), back_inserter(fHistList));
 // copy(fPhmap.begin(), fPhmap.end(), back_inserter(fHistList));
 // copy(fHits.begin(), fHits.end(), back_inserter(fHistList));
 // for (list<TH1*>::iterator il = fHistList.begin(); il != fHistList.end(); ++il) {
 //  	(*il)->Draw((getHistOption(*il)).c_str()); 
 // }
 // fDisplayedHist = find(fHistList.begin(), fHistList.end(), h1);
 // fDisplayedHist = find(fHistList.begin(), fHistList.end(), p2);
 // fDisplayedHist = find(fHistList.begin(), fHistList.end(), h1);
 // fDisplayedHist = find(fHistList.begin(), fHistList.end(), p2);
 // fDisplayedHist = find(fHistList.begin(), fHistList.end(), h2);
 // PixTest::update();

  FinalCleaning();
  fApi->setClockStretch(0, 0, 0); //No Stretch after trigger, 0 delay
  return rb_val;
 }
