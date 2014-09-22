#include <stdlib.h>  // atof, atoi
#include <algorithm> // std::find
#include <sstream>   // parsing

#include "PixTestPhOptimization.hh"
#include "log.h"

using namespace std;
using namespace pxar;

ClassImp(PixTestPhOptimization)

PixTestPhOptimization::PixTestPhOptimization() {}

PixTestPhOptimization::PixTestPhOptimization( PixSetup *a, std::string name ) :  PixTest(a, name), fParNtrig(-1), fParDAC("nada"), fParDacVal(100),   fFlagSinglePix(true), fSafetyMargin(20) {
  PixTest::init();
  init();
}


bool PixTestPhOptimization::setParameter(string parName, string sval) {
  bool found(false);
  string str1, str2;
  string::size_type s1;
  int pixc, pixr;
  std::transform(parName.begin(), parName.end(), parName.begin(), ::tolower);
  for (unsigned int i = 0; i < fParameters.size(); ++i) {
    if (!fParameters[i].first.compare(parName)) {
      found = true;
      sval.erase(remove(sval.begin(), sval.end(), ' '), sval.end());
      if (!parName.compare("ntrig")) {
	setTestParameter("ntrig", sval); 
	fParNtrig = atoi( sval.c_str() );
	LOG(logDEBUG) << "  setting fParNtrig  ->" << fParNtrig
		      << "<- from sval = " << sval;
      }
      if (!parName.compare("safetymargin")) {
	fSafetyMargin = atoi( sval.c_str() );
	LOG(logDEBUG) << "  setting fSafetyMargin  ->" << fSafetyMargin
		      << "<- from sval = " << sval;
      }
      if (!parName.compare("singlepix")) {
	fFlagSinglePix = atoi( sval.c_str() );
	LOG(logDEBUG) << "  setting fFlagSinglePix  ->" << fFlagSinglePix
		      << "<- from sval = " << sval;
      }
      if (!parName.compare("dac")) {
	setTestParameter("dac", sval); 
	fParDAC = sval;
	LOG(logDEBUG) << "  setting fParDAC  ->" << fParDAC
		      << "<- from sval = " << sval;
      }

      if (!parName.compare("dacval")) {
	setTestParameter("dacval", sval); 
	fParDacVal = atoi(sval.c_str());
	LOG(logDEBUG) << "  setting fParDacVal  ->" << fParDacVal
		      << "<- from sval = " << sval;
      }
      
      if (!parName.compare("pix")) {
        s1 = sval.find(",");
        if (string::npos != s1) {
	  str1 = sval.substr(0, s1);
	  pixc = atoi(str1.c_str());
	  str2 = sval.substr(s1+1);
	  pixr = atoi(str2.c_str());
	  fPIX.push_back(make_pair(pixc, pixr));
	  addSelectedPixels(sval); 
	  LOG(logDEBUG) << "  adding to FPIX ->" << pixc << "/" << pixr << " fPIX.size() = " << fPIX.size() ;
	} else {
	  clearSelectedPixels();
	  LOG(logDEBUG) << "  clear fPIX: " << fPIX.size(); 
	}
      }
      break;
    }
  }
  return found;
}

void PixTestPhOptimization::init() {
  fDirectory = gFile->GetDirectory(fName.c_str());
  if(!fDirectory) {
    fDirectory = gFile->mkdir(fName.c_str());
  }
  fDirectory->cd();
}

void PixTestPhOptimization::bookHist(string /*name*/) {}

PixTestPhOptimization::~PixTestPhOptimization() {}

void PixTestPhOptimization::doTest() {

  cacheDacs();
  bigBanner(Form("PixTestPhOptimization::doTest() Ntrig = %d, singlePix = %d", fParNtrig, (fFlagSinglePix?1:0)));
  fDirectory->cd();
  PixTest::update();

  TH1D *h1(0); 
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
  map<string, TH1D*> hists; 
  string name, title;

  //looking for inefficient pixels, so that they can be avoided
  std::vector<std::pair<uint8_t, pair<int,int> > > badPixels;
  BlacklistPixels(badPixels, 10);

  //vcal threshold map in order to choose the low-vcal value the PH will be sampled at
  fApi->_dut->testAllPixels(true);
  fApi->_dut->maskAllPixels(false);
  fApi->setDAC("ctrlreg",4);
  std::vector<pxar::pixel> thrmap;
  int cnt(0); 
  bool done(false);
  while (!done) {
    try {
      // Scanning the full VCal range, so no need to specify bounds:
      thrmap = fApi->getThresholdMap("vcal", FLAG_RISING_EDGE, 10);
      done = true;
    } catch(pxarException &e) {
      LOG(logCRITICAL) << "pXar execption: "<< e.what(); 
      ++cnt;
    }
    done = (cnt>5) || done;
  }

  int minthr=0;
  LOG(logDEBUG) << "thr map size "<<thrmap.size()<<endl;

  //flag allows to choose between PhOpt on single pixel (1) or on the whole roc (0)
  pair<int, pxar::pixel> maxpixel;
  pair<int, pxar::pixel> minpixel;
  vector<pair<int, pxar::pixel > > maxpixels;
  vector<pair<int, pxar::pixel > > minpixels;

  
  for (int iroc = 0; iroc < rocIds.size(); ++iroc){
    if(fFlagSinglePix){
      LOG(logDEBUG)<<"**********Ph range will be optimised on a single random pixel***********";
      pxar::pixel randomPix;
      randomPix= *(RandomPixel(badPixels, rocIds[iroc]));
      LOG(logDEBUG)<<"In doTest(), randomCol "<<(int)randomPix.column<<", randomRow "<<(int)randomPix.row<<", pixel "<<randomPix;
      maxpixel.first = rocIds[iroc]; 
      maxpixel.second.roc_id = rocIds[iroc];
      maxpixel.second.column = randomPix.column;
      maxpixel.second.row    = randomPix.row;
      minpixel.first = rocIds[iroc]; 
      minpixel.second.roc_id = rocIds[iroc];
      minpixel.second.column = randomPix.column;
      minpixel.second.row    = randomPix.row;
      LOG(logDEBUG)<<"random pixel: "<<maxpixel.second<<", "<<minpixel.second<<"is not on the blacklist";
      maxpixels.push_back(maxpixel);
      minpixels.push_back(minpixel);
      //retrieving info from the vcal thr map for THIS random pixel
    }
    else{
      LOG(logDEBUG)<<"**********Ph range will be optimised on the whole ROC***********";
      //getting highest ph pixel
      GetMaxPhPixel(maxpixels, badPixels);
      // getting pixel showing the largest vcal threshold (i.e., all other pixels are responding)
      GetMinPixel(minpixels, thrmap, badPixels);
    }
  }
  /*
  for(std::vector<pxar::pixel>::iterator thrit = thrmap.begin(); thrit != thrmap.end(); thrit++){
    if(thrit->column == minpixel.column && thrit->row == minpixel.row){
      minthr=static_cast<int>(thrit->getValue());
    }
  }

  fApi->_dut->testAllPixels(false);
  fApi->_dut->maskAllPixels(true);
  fApi->_dut->testPixel(maxpixel.column,maxpixel.row,true);
  fApi->_dut->maskPixel(maxpixel.column,maxpixel.row,false);
  fApi->setDAC("vcal",255);
  fApi->setDAC("ctrlreg",4);
  //scanning through offset and scale for max pixel (or randpixel)
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > dacdac_max;

  cnt = 0; 
  done = false;
  while (!done) {
    try {
      dacdac_max = fApi->getPulseheightVsDACDAC("phoffset",0,255,"phscale",0,255,0,10);
      done = true;
    } catch(pxarException &e) {
      LOG(logCRITICAL) << "pXar execption: "<< e.what(); 
      ++cnt;
    }
    done = (cnt>5) || done;
  }

  fApi->_dut->testAllPixels(false);
  fApi->_dut->maskAllPixels(true);
  fApi->_dut->testPixel(minpixel.column,minpixel.row,true);
  fApi->_dut->maskPixel(minpixel.column,minpixel.row,false);
  minthr = (minthr==255) ? (minthr) : (minthr + 5); 
  LOG(logDEBUG)<<"minthr, i.e. vcal value where PH is sampled on the low edge, is: "<<minthr;
  //minthr = 50;
  fApi->setDAC("ctrlreg",4);
  fApi->setDAC("vcal",minthr);
  //scanning through offset and scale for min pixel (or same randpixel)
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > dacdac_min;
  cnt = 0; 
  done = false;
  while (!done) {
    try {
      dacdac_min = fApi->getPulseheightVsDACDAC("phoffset",0,255,"phscale",0,255,0,10);
      done = true;
    } catch(pxarException &e) {
      LOG(logCRITICAL) << "pXar execption: "<< e.what(); 
      ++cnt;
    }
    done = (cnt>5) || done;
  }

  //search for optimal dac values in 3 steps
  //1. shrinking the PH to be completely inside the ADC range, adjusting phscale
  int ps_opt = 999, po_opt = 999;
  po_opt=200;
  ps_opt = InsideRangePH(po_opt, dacdac_max, dacdac_min);
  //check for opt failing
  if(ps_opt==999){
    LOG(logDEBUG)<<"PH optimization failed"<<endl<<"Please run PreTest or try PhOptimization on a random pixel";
    return;
  }
  //2. centring PH curve adjusting phoffset
  po_opt = CentrePhRange(po_opt, ps_opt, dacdac_max, dacdac_min);
  //3. stretching curve adjusting phscale
  ps_opt = StretchPH(po_opt, ps_opt, dacdac_max, dacdac_min);


  fApi->setDAC("ctrlreg",4);
  fApi->setDAC("phscale",ps_opt);
  fApi->setDAC("phoffset",po_opt);
  //draw PH curve for max and min pixel
  name  = Form("PH_c%d_r%d_C%d", maxpixel.column, maxpixel.row, 0);
  title = Form("PH_c%d_r%d_C%d, phscale = %d, phoffset = %d, maxpixel", maxpixel.column, maxpixel.row, 0, ps_opt, po_opt);
  h1 = bookTH1D(name, name, 256, 0., 256.);
  vector<pair<uint8_t, vector<pixel> > > results;
  fApi->_dut->testAllPixels(false);
  fApi->_dut->maskAllPixels(true);
  fApi->_dut->testPixel(maxpixel.column, maxpixel.row, true);
  fApi->_dut->maskPixel(maxpixel.column, maxpixel.row, false);
  cnt = 0; 
  done = false;
  while (!done) {
    try {
      results = fApi->getPulseheightVsDAC("vcal", 0, 255, FLAG_FORCE_MASKED, 10);
      done = true;
    } catch(pxarException &e) {
      LOG(logCRITICAL) << "pXar execption: "<< e.what(); 
      ++cnt;
    }
    done = (cnt>5) || done;
  }

  for (unsigned int iroc = 0; iroc < rocIds.size(); ++iroc){
    for (unsigned int i = 0; i < results.size(); ++i) {
      pair<uint8_t, vector<pixel> > v = results[i];
      int idac = v.first; 
      vector<pixel> vpix = v.second;
      for (unsigned int ipix = 0; ipix < vpix.size(); ++ipix) {
	if (vpix[ipix].roc_id == rocIds[iroc]) {
	  h1->Fill(idac, vpix[ipix].getValue());
	}
      }
    }
  }
  h1->SetMinimum(0);
  setTitles(h1, title.c_str(), "average PH");
  hists.insert(make_pair(name, h1));
  fHistList.push_back(h1);  

  results.clear();
  
  name  = Form("PH_c%d_r%d_C%d", minpixel.column, minpixel.row, 0);
  title = Form("PH_c%d_r%d_C%d, phscale = %d, phoffset = %d, minpixel", minpixel.column, minpixel.row, 0, ps_opt, po_opt);
  h1 = bookTH1D(name, name, 256, 0., 256.);
  fApi->_dut->testAllPixels(false);
  fApi->_dut->maskAllPixels(true);
  fApi->_dut->testPixel(minpixel.column, minpixel.row, true);
  fApi->_dut->maskPixel(minpixel.column, minpixel.row, false);
  cnt = 0; 
  done = false;
  while (!done) {
    try {
      results = fApi->getPulseheightVsDAC("vcal", 0, 255, FLAG_FORCE_MASKED, 10);
      done = true;
    } catch(pxarException &e) {
      LOG(logCRITICAL) << "pXar execption: "<< e.what(); 
      ++cnt;
    }
    done = (cnt>5) || done;
  }

  for (unsigned int iroc = 0; iroc < rocIds.size(); ++iroc){
    for (unsigned int i = 0; i < results.size(); ++i) {
      pair<uint8_t, vector<pixel> > v = results[i];
      int idac = v.first; 
      vector<pixel> vpix = v.second;
      for (unsigned int ipix = 0; ipix < vpix.size(); ++ipix) {
	if (vpix[ipix].roc_id == rocIds[iroc]) {
	  h1->Fill(idac, vpix[ipix].getValue());
	}
      }
    }
  }
  h1->SetMinimum(0);
  setTitles(h1, title.c_str(), "average PH");
  hists.insert(make_pair(name, h1));
  fHistList.push_back(h1);  
  */
  /*

  LOG(logDEBUG)<<"PH map for min vcal will now be done";
  std::vector<pxar::pixel> phMapMinVcal;
  fApi->_dut->testAllPixels(true);
  fApi->_dut->maskAllPixels(false);
  fApi->setDAC("ctrlreg",4);
  fApi->setDAC("vcal",minthr);
  fApi->setDAC("phscale",ps_opt);
  fApi->setDAC("phoffset",po_opt);
  phMapMinVcal = fApi->getPulseheightMap(0, 10);
  TH2D *h2(0); 
  name  = Form("PhMapMinVcal_SM%d", fSafetyMargin);
  h2 = bookTH2D(name, "PhMap at lower edge vcal value", 80, 0., 80., 52, 0., 52.);
  for(std::vector<pxar::pixel>::iterator px = phMapMinVcal.begin(); px != phMapMinVcal.end(); px++) {
    h2->Fill(px->row, px->column, px->getValue());
  }
  for(std::vector<std::pair<int, int> >::iterator bad_it = badPixels.begin(); bad_it != badPixels.end(); bad_it++){
	h2->SetBinContent(h2->GetXaxis()->FindBin(bad_it->second), h2->GetYaxis()->FindBin(bad_it->first), -1);
  }
  //  hists.insert(make_pair("PhMapMinVcal", h2));
  fHistList.push_back(h2);
  fHistOptions.insert(make_pair(h2, "colz"));

  LOG(logDEBUG)<<"PH map for max vcal will now be done";
  std::vector<pxar::pixel> phMapMaxVcal;
  fApi->setDAC("ctrlreg",4);
  fApi->setDAC("vcal",250);
  phMapMaxVcal = fApi->getPulseheightMap(0, 10);
  TH2D *h3(0); 
  name  = Form("PhMapMaxVcal_SM%d", fSafetyMargin);
  h3 = bookTH2D(name, "PhMap at upper edge vcal value", 80, 0., 80., 52, 0., 52.);
  for(std::vector<pxar::pixel>::iterator px = phMapMaxVcal.begin(); px != phMapMaxVcal.end(); px++) {
    h3->Fill(px->row, px->column, px->getValue());
  }
  for(std::vector<std::pair<int, int> >::iterator bad_it = badPixels.begin(); bad_it != badPixels.end(); bad_it++){
	h3->SetBinContent(h3->GetXaxis()->FindBin(bad_it->second), h3->GetYaxis()->FindBin(bad_it->first), -1);
  }
  fHistList.push_back(h3);
  fHistOptions.insert(make_pair(h3, "colz"));
  */

//  LOG(logDEBUG)<<"Just before entering DynamicRange()";
//  DynamicRange();
/*
  for (list<TH1*>::iterator il = fHistList.begin(); il != fHistList.end(); ++il) {
    (*il)->Draw(getHistOption(*il).c_str());
    PixTest::update();
  }
  fDisplayedHist = find(fHistList.begin(), fHistList.end(), h1);
  restoreDacs(); 

  // -- FIXME: should be ROC specific! Must come after restoreDacs()!
  fApi->setDAC("phscale",ps_opt);
  fApi->setDAC("phoffset",po_opt);
  saveDacs();

  // -- print summary information
  string psString(""), poString(""); 
  for (unsigned int i = 0; i < rocIds.size(); ++i) {
    psString += Form(" %3d", fApi->_dut->getDAC(rocIds[i], "phscale"));
    poString += Form(" %3d", fApi->_dut->getDAC(rocIds[i], "phoffset"));
  }
  LOG(logINFO) << "PixTestPhOptimization::doTest() done";
  LOG(logINFO) << "PH scale (per ROC):  " << psString;
  LOG(logINFO) << "PH offset (per ROC): " << poString;
*/
}

void PixTestPhOptimization::DynamicRange(){
  LOG(logDEBUG)<<"Welcome to DynamicRange subtest";

  fApi->setDAC("ctrlreg", 4); 
  
  TH1D *h1= new TH1D("h1", "h1", 255, 0., 255.);
  TH2D *h2(0);
  string name;
  double range;
  name  = Form("PhMapDynamicRange_SM%d", fSafetyMargin);
  h2 = bookTH2D(name, "Map of dynamic range for optimized PH", 80, 0., 80., 52, 0., 52.);
  vector<pair<uint8_t, vector<pixel> > > result;
  fApi->_dut->testAllPixels(false);
  fApi->_dut->maskAllPixels(true);
  for (unsigned int ir = 0; ir <80 ; ++ir) {
    for (unsigned int ic = 0; ic < 52; ++ic) {
      LOG(logDEBUG)<<"Preparing to have PHvsVcal for pix "<<ic<<","<<ir;
      fApi->_dut->testPixel(ic, ir, true);
      fApi->_dut->maskPixel(ic, ir, false);
      result = fApi->getPulseheightVsDAC("vcal", 0, 255, FLAG_FORCE_MASKED, 10);
      LOG(logDEBUG)<<"PHvsVcal for pix "<<ic<<","<<ir<<" acquired";
      LOG(logDEBUG)<<"Size of PHvsVcal vector is "<<result.size();
      for (int ires=0; ires< result.size(); ires++){
	pair<uint8_t, vector<pixel> > v = result[ires];
	int idac = v.first;
	vector<pixel> vpix = v.second;
	LOG(logDEBUG)<<"Size of pixel vector is "<<vpix.size();
	for(int ipix=0; ipix<vpix.size(); ipix++){
	  h1->Fill(idac, vpix[0].getValue());
	}
      }
      range =  h1->GetBinContent(h1->GetMaximumBin()) - h1->GetMinimum(0.1);
      h2->Fill(ir, ic, range);

      h1->Reset();
      fApi->_dut->testAllPixels(false);
      fApi->_dut->maskAllPixels(true);
    }
  }

  fHistList.push_back(h2);
  fHistOptions.insert(make_pair(h2, "colz"));

}


void PixTestPhOptimization::BlacklistPixels(std::vector<std::pair<uint8_t, pair<int, int> > > &badPixels, int aliveTrig){
  //makes a list of inefficient pixels, to be avoided during optimization
  fApi->_dut->testAllPixels(true);
  fApi->_dut->maskAllPixels(false);

  vector<uint8_t> vVcal = getDacs("vcal"); 
  vector<uint8_t> vCreg = getDacs("ctrlreg"); 

  vector<TH2D*> testEff = efficiencyMaps("PixelAlive", aliveTrig);
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
  std::pair<uint8_t, pair<int, int> > badPix;
  Double_t eff=0.;
  for(uint8_t rocid = 0; rocid<rocIds.size(); rocid++){
    for(int r=0; r<80; r++){
      for(int c=0; c<52; c++){
	eff = testEff[rocIds[rocid]]->GetBinContent( testEff[rocIds[rocid]]->FindFixBin((double)c + 0.5, (double)r+0.5) );
	if(eff<aliveTrig){
	  LOG(logDEBUG)<<"Pixel ["<<rocIds[rocid]<<", "<<c<<", "<<r<<"] has eff "<<eff<<"/"<<aliveTrig;
	  badPix.first = rocIds[rocid];
	  badPix.second.first = c;
	  badPix.second.second = r;
	  LOG(logDEBUG)<<"bad Pixel found and blacklisted: ["<<badPix.first<<", "<<badPix.second.first<<", "<<badPix.second.second;
	  (badPixels).push_back(badPix);
	}
      }
    }
  }
  setDacs("vcal", vVcal); 
  setDacs("ctrlreg", vCreg); 
  LOG(logDEBUG)<<"Number of bad pixels found: "<<badPixels.size();
}


pxar::pixel* PixTestPhOptimization::RandomPixel(std::vector<std::pair<uint8_t, pair<int, int> > > &badPixels, uint8_t iroc){
  //Returns a random pixel, taking care it is not on the blacklist
  fApi->setDAC("ctrlreg",4);
  bool isPixGood=true;
  pxar::pixel *randPixel= new pixel();
  srand(int(time(NULL)));
  int random_col=-1, random_row=-1;
  do{
    random_col = rand() % 52;
    random_row = rand() % 80;
    LOG(logDEBUG)<<"random pixel: ["<<iroc<<", "<<random_col<<", "<<random_row<<"]";
    isPixGood=true;
    for(std::vector<std::pair<uint8_t, pair<int, int> > >::iterator bad_it = badPixels.begin(); bad_it != badPixels.end(); bad_it++){
      if(bad_it->first == iroc && bad_it->second.first == random_col && bad_it->second.second == random_row){
	isPixGood=false;
      }
    }
    LOG(logDEBUG)<<"is the random pixel good? "<<isPixGood;
  }while(!isPixGood);
  randPixel->roc_id = iroc;
  randPixel->column = random_col;
  randPixel->row    = random_row;
  LOG(logDEBUG)<<"In RandomPixel(), rocId "<<iroc<<", randomCol "<<(int)randPixel->column<<", randomRow "<<(int)randPixel->row<<", pixel "<<randPixel;
  return randPixel;
}


void PixTestPhOptimization::GetMaxPhPixel(vector<pair <int, pxar::pixel > > &maxpixels,   std::vector<std::pair<uint8_t, pair<int,int> > >  &badPixels){
  //looks for the pixel with the highest Ph at vcal = 255, taking care the pixels are not already saturating (ph=255)
    fApi->_dut->testAllPixels(true);
    fApi->_dut->maskAllPixels(false);
    vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
    bool isPixGood=true;
    int maxph = 255;
    /*********changed to test module w v2 chips***** put it back at 255 later********/
    int init_phScale =25;
    int flag_maxPh=0;
    pair<int, pxar::pixel> maxpixel;
    maxpixel.second.setValue(0);
    std::vector<pxar::pixel> result;
    while(maxph>254 && flag_maxPh<52){
      result.clear();
      fApi->setDAC("phscale", init_phScale);
      fApi->setDAC("vcal",255);
      fApi->setDAC("ctrlreg",4);
      fApi->setDAC("phoffset",200);  
      int cnt(0); 
      bool done(false);
      while (!done) {
	try {
	  result = fApi->getPulseheightMap(0, 10);
	  done = true;
	} catch(pxarException &e) {
	  LOG(logCRITICAL) << "pXar execption: "<< e.what(); 
	  ++cnt;
	}
	done = (cnt>5) || done;
      }
      
      maxph=0;
      LOG(logDEBUG) << "result size "<<result.size()<<endl;
      for(std::vector<pxar::pixel>::iterator px = result.begin(); px != result.end(); px++) {
	isPixGood=true;
	for(std::vector<std::pair<uint8_t, pair<int, int> > >::iterator bad_it = badPixels.begin(); bad_it != badPixels.end(); bad_it++){
	  if(bad_it->second.first == px->column && bad_it->second.second == px->row && bad_it->first == px->roc_id){
	    isPixGood=false;
	  }
	}
	if(isPixGood && px->getValue() > maxph){
	  maxph = px->getValue();
	}
      }
      //should have flag for v2 or v2.1
      /*********changed to test module w v2 chips***** put it back at -=5 later********/
      init_phScale+=5;
      flag_maxPh++;
    }
    
    
      // Look for pixel with max. pulse height:

    for(int iroc=0; iroc< rocIds.size(); iroc++){
      maxph=0;
      for(std::vector<pxar::pixel>::iterator px = result.begin(); px != result.end(); px++) {
	isPixGood=true;
	if(px->getValue() > maxpixel.second.getValue() && px->roc_id == rocIds[iroc]){
	  for(std::vector<std::pair<uint8_t, pair<int, int> > >::iterator bad_it = badPixels.begin(); bad_it != badPixels.end(); bad_it++){
	    if(bad_it->second.first == px->column && bad_it->second.second == px->row && bad_it->first == px->roc_id){
	      isPixGood=false;
	    }
	  }
	  if(isPixGood){
	    maxpixel = make_pair(rocIds[iroc],*px);
	    maxph = px->getValue();
	    }
	}
      }
      maxpixels.push_back(maxpixel);
    }
    LOG(logDEBUG) << "maxPh " << maxph <<" for ROC "<<maxpixel.first<<" on pixel "<<maxpixel.second << endl ;      

}




void PixTestPhOptimization::GetMinPixel(vector<pair<int, pxar::pixel> > &minpixels, std::vector<pxar::pixel> &thrmap,   std::vector<std::pair<uint8_t, pair<int,int> > > &badPixels){
  //the minimum pixel is the pixel showing the largest vcal threshold: it is the smallest vcal we can probe the Ph for all pixels, and this pix has the smallest ph for this vcal
  int minthr=0;
  pxar::pixel minpixel;
  bool  isPixGood=true;
  fApi->_dut->testAllPixels(true);
  fApi->_dut->maskAllPixels(false);  
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
  for(int iroc=0; iroc<rocIds.size(); iroc++){
    for(std::vector<pxar::pixel>::iterator thrit = thrmap.begin(); thrit != thrmap.end(); thrit++){
      if(thrit->roc_id != rocIds[iroc]){
	continue;
      }
      isPixGood=true;
      for(std::vector<std::pair<uint8_t, pair<int, int> > >::iterator bad_it = badPixels.begin(); bad_it != badPixels.end(); bad_it++){
	if(bad_it->second.first == thrit->column && bad_it->second.second == thrit->row && bad_it->first ==  thrit->roc_id){
	  isPixGood=false;
	  break;
	}
      }
      if(thrit->getValue() > minthr && isPixGood) {
	minpixel = *thrit;
	minthr = static_cast<int>(minpixel.getValue());
	
      }
    }
    minpixels.push_back(make_pair(minthr, minpixel));
  }
  LOG(logDEBUG) << "min vcal thr " << minthr << "found for pixel "<<minpixel<<endl ;
}


int PixTestPhOptimization::InsideRangePH(int po_opt,  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_max,   std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_min){
  //adjusting phscale so that the PH curve is fully inside the ADC range
  int ps_opt = 999;
  int maxPh(0);
  int minPh(0);
  bool lowEd=false, upEd=false;
  int upEd_dist=255, lowEd_dist=255;
  int safetyMargin = 50;
  int dist = 255;
  int bestDist = 255;
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_max = dacdac_max.begin();
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_min = dacdac_min.begin();
  //or two for cycles??
  LOG(logDEBUG) << "dacdac at max vcal has size "<<dacdac_max.size()<<endl;
  LOG(logDEBUG) << "dacdac at min vcal has size "<<dacdac_min.size()<<endl;
  while(dacit_max != dacdac_max.end() || dacit_min != dacdac_min.end()){
    if(dacit_max->first == po_opt && dacit_min->first == po_opt && dacit_min->second.second.size() && dacit_max->second.second.size()) {
      maxPh=dacit_max->second.second[0].getValue();
      minPh=dacit_min->second.second[0].getValue();
      lowEd = (minPh > safetyMargin);
      upEd = (maxPh < 255 - safetyMargin);
      upEd_dist = abs(maxPh - (255 - safetyMargin));
      lowEd_dist = abs(minPh - safetyMargin);
      dist = (upEd_dist > lowEd_dist ) ? (upEd_dist) : (lowEd_dist);
      if(dist < bestDist && upEd && lowEd){
	ps_opt = dacit_max->second.first;
	bestDist=dist;
      }
    }
    dacit_max++;
    dacit_min++;
  }
  LOG(logDEBUG)<<"opt step 1: po fixed to"<<po_opt<<" and scale adjusted to "<<ps_opt<<", with distance "<<bestDist;
  return ps_opt;
}



int PixTestPhOptimization::CentrePhRange(int po_opt_in, int ps_opt,  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_max,   std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_min){
  //centring PH curve adjusting phoffset   
  int po_opt_out = po_opt_in;
  int maxPh(0);
  int minPh(0);
  int dist = 255;
  int bestDist = 255;
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_max = dacdac_max.begin();
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_min = dacdac_min.begin();
  //or two for cycles??
  while(dacit_max != dacdac_max.end() || dacit_min != dacdac_min.end()){
    if(dacit_max->second.first == ps_opt && dacit_min->second.first == ps_opt && dacit_min->second.second.size() && dacit_max->second.second.size()) {
      maxPh=dacit_max->second.second.at(0).getValue();
      minPh=dacit_min->second.second.at(0).getValue();
      dist = abs(minPh - (255 - maxPh));
      if (dist < bestDist){
	po_opt_out = dacit_max->first;
	bestDist = dist;
      } 
    }
    dacit_max++;
    dacit_min++;
  }
  LOG(logDEBUG)<<"opt centring step: po "<<po_opt_out<<" and scale "<<ps_opt<<", with distance "<<bestDist;
  return po_opt_out;
}



int PixTestPhOptimization::StretchPH(int po_opt, int ps_opt_in,  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_max,   std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_min){
  //stretching PH curve to exploit the full ADC range, adjusting phscale             
  int ps_opt_out = ps_opt_in;
  int maxPh(0);
  int minPh(0);
  bool lowEd=false, upEd=false;
  int upEd_dist=255, lowEd_dist=255;
  int safetyMargin = fSafetyMargin;
  LOG(logDEBUG)<<"safety margin for stretching set to "<<fSafetyMargin;
  int dist = 255;
  int bestDist = 255;
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_max = dacdac_max.begin();
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_min = dacdac_min.begin();
  while(dacit_max != dacdac_max.end() || dacit_min != dacdac_min.end()){
    if(dacit_max->first == po_opt && dacit_min->first == po_opt && dacit_min->second.second.size() && dacit_max->second.second.size()) {
      maxPh=dacit_max->second.second.at(0).getValue();
      minPh=dacit_min->second.second.at(0).getValue();
      lowEd = (minPh > safetyMargin);
      upEd = (maxPh < 255 - safetyMargin);
      upEd_dist = abs(maxPh - (255 - safetyMargin));
      lowEd_dist = abs(minPh - safetyMargin);
      dist = (upEd_dist < lowEd_dist ) ? (upEd_dist) : (lowEd_dist);
      if(dist<bestDist && lowEd && upEd){
	ps_opt_out = dacit_max->second.first;
	bestDist=dist;
      }
    }
    dacit_max++;
    dacit_min++;
  }
  LOG(logDEBUG)<<"opt final step: po fixed to"<<po_opt<<" and scale adjusted to "<<ps_opt_out<<", with distance "<<bestDist;
  return ps_opt_out;
}
