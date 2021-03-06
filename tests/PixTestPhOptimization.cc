#include <stdlib.h>  // atof, atoi
#include <algorithm> // std::find
#include <sstream>   // parsing

#include <TStopwatch.h>

#include "PixTestPhOptimization.hh"
#include "log.h"

using namespace std;
using namespace pxar;

ClassImp(PixTestPhOptimization)

PixTestPhOptimization::PixTestPhOptimization() {}

PixTestPhOptimization::PixTestPhOptimization( PixSetup *a, std::string name ) :  PixTest(a, name), fParNtrig(-1), fParDAC("nada"), fParDacVal(100),   fFlagSinglePix(true), fSafetyMarginUp(10), fSafetyMarginLow(15) {
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
      if (!parName.compare("safetymarginup")) {
	fSafetyMarginUp = atoi( sval.c_str() );
	LOG(logDEBUG) << "  setting fSafetyMarginUp  ->" << fSafetyMarginUp
		      << "<- from sval = " << sval;
      }
      if (!parName.compare("safetymarginlow")) {
	fSafetyMarginLow = atoi( sval.c_str() );
	LOG(logDEBUG) << "  setting fSafetyMarginLow  ->" << fSafetyMarginLow
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

  TStopwatch t;

  cacheDacs();
  bigBanner(Form("PixTestPhOptimization::doTest() Ntrig = %d, singlePix = %d", fParNtrig, (fFlagSinglePix?1:0)));
  fDirectory->cd();
  PixTest::update();

  TH1D *h1(0); 
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
//  LOG(logDEBUG)<<"Enabled ROCs vector has size: "<<rocIds.size();
//  LOG(logDEBUG)<<"ROC "<<(int)rocIds[0]<<" is enabled";
//  for(unsigned int iroc=0; iroc < rocIds.size(); iroc++){
//    LOG(logDEBUG)<<"ROC "<<(int)rocIds[iroc]<<" is enabled";
//  }
  string name, title;

  //looking for inefficient pixels, so that they can be avoided
  std::vector<std::pair<uint8_t, pair<int,int> > > badPixels;
  BlacklistPixels(badPixels, 10);

  //Set from trim file  minimum vcal where PH is sampled for minphpix search
  SetMinThr();

  //flag allows to choose between PhOpt on single pixel (1) or on the whole roc (0)
  pair<int, pxar::pixel> maxpixel;
  pair<int, pxar::pixel> minpixel;
  map<int, pxar::pixel> maxpixels;
  map<int, pxar::pixel> minpixels;
  map<int, int> minVcal;
  
  if(fFlagSinglePix){
    for (unsigned int iroc = 0; iroc < rocIds.size(); ++iroc){
      LOG(logDEBUG)<<"**********Ph range will be optimised on a single random pixel***********";
      pxar::pixel randomPix;
      randomPix= *(RandomPixel(badPixels, rocIds[iroc]));
      LOG(logDEBUG)<<"In doTest(), randomCol "<<(int)randomPix.column()<<", randomRow "<<(int)randomPix.row()<<", pixel "<<randomPix;
      maxpixel.first = rocIds[iroc]; 
      maxpixel.second.setRoc(rocIds[iroc]);
      maxpixel.second.setColumn(randomPix.column());
      maxpixel.second.setRow(randomPix.row());
      minpixel.first = iroc; 
      minpixel.second.setRoc(rocIds[iroc]);
      minpixel.second.setColumn(randomPix.column());
      minpixel.second.setRow(randomPix.row());
      LOG(logDEBUG)<<"random pixel: "<<maxpixel.second<<", "<<minpixel.second<<"is not on the blacklist";
      maxpixels.insert(maxpixel);
      minpixels.insert(minpixel);
    }
  }
  else{
    LOG(logDEBUG)<<"**********Ph range will be optimised on the whole ROC***********";
    //getting highest ph pixel
    GetMaxPhPixel(maxpixels, badPixels);
    //getting min ph pixel and finding its vcal threshold
    GetMinPhPixel(minpixels, minVcal, badPixels);
  }

  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    LOG(logDEBUG)<<"vcal min "<<minVcal[roc_it]<<" on ROC"<<(int)rocIds[roc_it];
  }

  //scan phoffset and phscale for max and min ph pixels
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > dacdac_max;
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > dacdac_min;
  MaxPhVsDacDac(dacdac_max, maxpixels);
  MinPhVsDacDac(dacdac_min, minpixels, minVcal);

  //search for optimal dac values in 3 steps
  //1. shrinking the PH to be completely inside the ADC range, adjusting phscale
  map<uint8_t, int> ps_opt, po_opt;
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    po_opt[rocIds[roc_it]] = 120;
  }
  ps_opt = InsideRangePH(po_opt, dacdac_max, dacdac_min);
  //check for opt failure
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    if(ps_opt[rocIds[roc_it]]==999){
      LOG(logDEBUG)<<"PH optimization failed on ROC "<<(int)rocIds[roc_it]<<endl<<"Please run PreTest or try PhOptimization on a random pixel";
    }
  }
  //2. centring PH curve adjusting phoffset
  po_opt = CentrePhRange(po_opt, ps_opt, dacdac_max, dacdac_min);
  
  //3. stretching curve adjusting phscale
  ps_opt = StretchPH(po_opt, ps_opt, dacdac_max, dacdac_min);
  
  //set optimized dacs and save
  restoreDacs();
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    fApi->setDAC("phscale",ps_opt[rocIds[roc_it]], rocIds[roc_it] );
    fApi->setDAC("phoffset",po_opt[rocIds[roc_it]], rocIds[roc_it]);
  }
  saveDacs();

  cacheDacs(); 
  //draw figures of merit of optimization
  DrawPhMaps(minVcal, badPixels);
  DrawPhCurves(maxpixels, minpixels, po_opt, ps_opt);
  restoreDacs();

  for (list<TH1*>::iterator il = fHistList.begin(); il != fHistList.end(); ++il) {
    (*il)->Draw(getHistOption(*il).c_str());
    PixTest::update();
  }
  fDisplayedHist = find(fHistList.begin(), fHistList.end(), h1);

  // -- print summary information
  string psString(""), poString(""); 
  for (unsigned int i = 0; i < rocIds.size(); ++i) {
    psString += Form(" %3d", fApi->_dut->getDAC(rocIds[i], "phscale"));
    poString += Form(" %3d", fApi->_dut->getDAC(rocIds[i], "phoffset"));
  }

  int seconds = t.RealTime(); 
  LOG(logINFO) << "PixTestPhOptimization::doTest() done, duration: " << seconds << " seconds";
  LOG(logINFO) << "PH scale (per ROC):  " << psString;
  LOG(logINFO) << "PH offset (per ROC): " << poString;
  
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
	eff = testEff[rocid]->GetBinContent( testEff[rocid]->FindFixBin((double)c + 0.5, (double)r+0.5) );
	if(eff<aliveTrig){
	  LOG(logDEBUG)<<"Pixel ["<<(int)rocIds[rocid]<<", "<<(int)c<<", "<<(int)r<<"] has eff "<<eff<<"/"<<aliveTrig;
	  badPix.first = rocIds[rocid];
	  badPix.second.first = c;
	  badPix.second.second = r;
	  LOG(logDEBUG)<<"bad Pixel found and blacklisted: ["<<(int)badPix.first<<", "<<(int)badPix.second.first<<", "<<(int)badPix.second.second<<"]";
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
  randPixel->setRoc(iroc);
  randPixel->setColumn(random_col);
  randPixel->setRow(random_row);
  LOG(logDEBUG)<<"In RandomPixel(), rocId "<<iroc<<", randomCol "<<(int)randPixel->column()<<", randomRow "<<(int)randPixel->row()<<", pixel "<<randPixel;
  return randPixel;
}

void PixTestPhOptimization::GetMaxPhPixel(map<int, pxar::pixel > &maxpixels,   std::vector<std::pair<uint8_t, pair<int,int> > >  &badPixels){
  //looks for the pixel with the highest Ph at vcal = 255, taking care the pixels are not already saturating (ph=255)
    fApi->_dut->testAllPixels(true);
    fApi->_dut->maskAllPixels(false);
    vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
    bool isPixGood=true;
    int maxph = 255;
    fApi->setDAC("phoffset", 200);
    int init_phScale =200;
    int flag_maxPh=0;
    pair<int, pxar::pixel> maxpixel;
    maxpixel.second.setValue(0);
    std::vector<pxar::pixel> result;
    while((maxph>254 || maxph==0) && flag_maxPh<52){
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
      //check that the pixel showing highest PH on the module is not reaching 255
      for(std::vector<pxar::pixel>::iterator px = result.begin(); px != result.end(); px++) {
	isPixGood=true;
	for(std::vector<std::pair<uint8_t, pair<int, int> > >::iterator bad_it = badPixels.begin(); bad_it != badPixels.end(); bad_it++){
	  if(bad_it->second.first == px->column() && bad_it->second.second == px->row() && bad_it->first == px->roc()){
	    isPixGood=false;
	  }
	}
	if(isPixGood && px->value() > maxph){
	  maxph = px->value();
	}
      }
      //should have flag for v2 or v2.1
      init_phScale+=5;
      flag_maxPh++;
    }
    // Look for pixel with max. pulse height on every ROC:
    for(unsigned int iroc=0; iroc< rocIds.size(); iroc++){
      maxph=0;
      for(std::vector<pxar::pixel>::iterator px = result.begin(); px != result.end(); px++) {
	isPixGood=true;
	if(px->value() > maxph && px->roc() == rocIds[iroc]){
	  for(std::vector<std::pair<uint8_t, pair<int, int> > >::iterator bad_it = badPixels.begin(); bad_it != badPixels.end(); bad_it++){
	    if(bad_it->second.first == px->column() && bad_it->second.second == px->row() && bad_it->first == px->roc()){
	      isPixGood=false;
	      break;
	    }
	  }
	  if(isPixGood){
	    maxpixel = make_pair(iroc,*px);
	    maxph = px->value();
	    }
	}
      }
      maxpixels.insert(maxpixel);
      LOG(logDEBUG) << "maxPh " << maxph <<" for ROC "<<maxpixel.first<<" on pixel "<<maxpixel.second << endl ;      
    }
}

void PixTestPhOptimization::GetMinPhPixel(map<int, pxar::pixel > &minpixels, map<int, int> &minVcal, std::vector<std::pair<uint8_t, pair<int,int> > >  &badPixels){
  //looks for the pixel with the lowest Ph at vcal = 1.2*trimValue, taking care the pixels are correclty responding (ph>0)
  //finds the min vcal at which the minphpixel can be sampled
  fApi->_dut->testAllPixels(true);
  fApi->_dut->maskAllPixels(false);
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
  bool isPixGood=true;
  int minph = 0;
  int init_phScale = 100;
  int flag_minPh=0;
  pair<int, pxar::pixel> minpixel;
  minpixel.second.setValue(0);
  std::vector<pxar::pixel> result;
  while(minph<1 && flag_minPh<52){
    result.clear();
    fApi->setDAC("phscale", init_phScale);
    fApi->setDAC("ctrlreg",0);
    //    fApi->setDAC("vcal",fMinThr*1.2);
    fApi->setDAC("vcal",130);
    fApi->setDAC("phoffset",150);  
    int cnt(0); 
    bool done(false);
    int size = 0;
    while (!(done && size !=0)) {
      try {
	result = fApi->getPulseheightMap(0, 10);
	size = result.size();
	++cnt;
	done = true;
      } catch(pxarException &e) {
	LOG(logCRITICAL) << "pXar execption: "<< e.what(); 
	++cnt;
	}
      done = (cnt>5) || done;
    }
    minph=255;
    LOG(logDEBUG) << "result size "<<result.size()<<endl;
    //check that the pixel showing lowest PH above 0 on the module
    for(std::vector<pxar::pixel>::iterator px = result.begin(); px != result.end(); px++) {
      isPixGood=true;
      for(std::vector<std::pair<uint8_t, pair<int, int> > >::iterator bad_it = badPixels.begin(); bad_it != badPixels.end(); bad_it++){
	if(bad_it->second.first == px->column() && bad_it->second.second == px->row() && bad_it->first == px->roc()){
	  isPixGood=false;
	  break;
	  }
      }
      if(isPixGood && px->value() < minph){
	minph = px->value();
      }
    }
    //should have flag for v2 or v2.1
    init_phScale+=5;
    flag_minPh++;
  }
  // Look for pixel with min pulse height on every ROC:
  for(unsigned int iroc=0; iroc< rocIds.size(); iroc++){
    minph=255;
    for(std::vector<pxar::pixel>::iterator px = result.begin(); px != result.end(); px++) {
      isPixGood=true;
      if(px->value() < minph && px->roc() == rocIds[iroc]){
	for(std::vector<std::pair<uint8_t, pair<int, int> > >::iterator bad_it = badPixels.begin(); bad_it != badPixels.end(); bad_it++){
	  if(bad_it->second.first == px->column() && bad_it->second.second == px->row() && bad_it->first == px->roc()){
	    isPixGood=false;
	    break;
	  }
	}
	if(isPixGood){
	  minpixel = make_pair(iroc,*px);
	  minph = px->value();
	}
      }
    }
    minpixels.insert(minpixel);
    LOG(logDEBUG) << "minPh " << minph <<" for ROC "<<minpixel.first<<" on pixel "<<minpixel.second << endl ;      
  }
  
  //finds min vcal
  int cnt(0); 
  bool done(false);
  TH1D *h1(0); 
  vector<pair<uint8_t, vector<pixel> > > results;
  pair<int, int> vcalmin;
  int vcalthr = 0;
  h1 = bookTH1D("h1", "h1", 256, 0., 256.);
  unsigned int NRocs = rocIds.size();
  for(unsigned int roc_it = 0; roc_it < NRocs; roc_it++){
    for(unsigned int roc_kt = 0; roc_kt < NRocs; roc_kt++){
      fApi->_dut->setROCEnable(roc_kt, true);
    }
    fApi->_dut->testAllPixels(false);
    fApi->_dut->maskAllPixels(true);
    //    fApi->_dut->setROCEnable(roc_it, true);
    //    fApi->_dut->testPixel(minpixels[rocIds[roc_it]].column(), minpixels[rocIds[roc_it]].row(), true, rocIds[roc_it]);
    //    fApi->_dut->maskPixel(minpixels[rocIds[roc_it]].column(), minpixels[rocIds[roc_it]].row(), false, rocIds[roc_it]);
    fApi->_dut->testPixel(minpixels[rocIds[roc_it]].column(), minpixels[rocIds[roc_it]].row(), true);
    fApi->_dut->maskPixel(minpixels[rocIds[roc_it]].column(), minpixels[rocIds[roc_it]].row(), false);
    LOG(logDEBUG)<<"enabling pixels "<<(int)minpixels[roc_it].column()<<", "<<(int)minpixels[roc_it].row()<<", "<<(int)minpixels[roc_it].roc()<<" "<<(int)roc_it;
    for(unsigned int roc_jt = 0; roc_jt < NRocs; roc_jt++){
      if(roc_jt!= roc_it){
	fApi->_dut->setROCEnable(roc_jt, false);
      }
    }
    //    fApi->_dut->testPixel(minpixels[roc_it].column(), minpixels[roc_it].row(), true, getIdFromIdx((int)roc_it));
    //fApi->_dut->maskPixel(minpixels[roc_it].column(), minpixels[roc_it].row(), false, getIdFromIdx((int)roc_it));
    fApi->_dut->info();
    cnt = 0; 
    done = false;
    while (!done) {
      try {
	//	LOG(logDEBUG)<<"trying phVsDAC for roc "<<(int)roc_it<<" "<<(int)rocIds[roc_it];
	results = fApi->getPulseheightVsDAC("vcal", 0, 255, FLAG_FORCE_MASKED, 10);
	//LOG(logDEBUG)<<"worked? phVsDAC for roc "<<(int)roc_it<<" "<<(int)rocIds[roc_it];
	done = true;
      } catch(pxarException &e) {
	LOG(logCRITICAL) << "pXar execption: "<< e.what(); 
	++cnt;
      }
      done = (cnt>5) || done;
    }
    
    LOG(logDEBUG)<<"size of results "<<results.size();
    for (unsigned int i = 0; i < results.size(); ++i) {
      //      LOG(logDEBUG)<<"analyzing postion "<<int(i)<<" of results";
      pair<uint8_t, vector<pixel> > v = results[i];
      int idac = v.first; 
      // LOG(logDEBUG)<<"idac is "<<(int)idac;
      vector<pixel> vpix = v.second;
      //LOG(logDEBUG)<<"vpix size is "<<vpix.size();
      for (unsigned int ipix = 0; ipix < vpix.size(); ++ipix) {
	//	LOG(logDEBUG)<<"vcalmin loop: ipix = "<<ipix<<", dac = "<<idac<<", ph = "<<vpix[ipix].value();
	h1->Fill(idac, vpix[ipix].value());
      }
    }
    vcalthr = static_cast<int>( h1->GetBinCenter( h1->FindFirstBinAbove(1.) ) );
    vcalmin = make_pair(roc_it, vcalthr);
    minVcal.insert(vcalmin);
    h1->Reset();
  }
  for(unsigned int roc_kt = 0; roc_kt < NRocs; roc_kt++){
    fApi->_dut->setROCEnable(roc_kt, false);
  }
  for(unsigned int roc_kt = 0; roc_kt < NRocs; roc_kt++){
    fApi->_dut->setROCEnable(roc_kt, true);
  }
}

map<uint8_t, int> PixTestPhOptimization::InsideRangePH(map<uint8_t,int> &po_opt,  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_max,   std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_min){
  //adjusting phscale so that the PH curve is fully inside the ADC range
  map<uint8_t, int> ps_opt;
  int maxPh(0);
  int minPh(0);
  bool lowEd=false, upEd=false;
  int upEd_dist=255, lowEd_dist=255;
  int safetyMargin = 40;
  int dist = 255;
  map<uint8_t, int> bestDist;
  LOG(logDEBUG) << "dacdac at max vcal has size "<<dacdac_max.size()<<endl;
  LOG(logDEBUG) << "dacdac at min vcal has size "<<dacdac_min.size()<<endl;
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    bestDist[rocIds[roc_it]] = 255;
    LOG(logDEBUG)<<"Bestdist at roc_it "<<roc_it<<" initialized with "<<bestDist[roc_it]<<" "<<bestDist[rocIds[roc_it]];
    ps_opt[rocIds[roc_it]] = 999;
  }
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_max = dacdac_max.begin();
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_min = dacdac_min.begin();
//  for(dacit_max = dacdac_max.begin(); dacit_max != dacdac_max.end(); dacit_max++){
//    for(int pixit = 0; pixit < dacit_max->second.second.size(); pixit++){
//      LOG(logDEBUG)<<"dacit_max: pixel "<<(int)dacit_max->second.second[pixit].roc()<<", "<<(int)dacit_max->second.second[pixit].column()<<", "<<(int)dacit_max->second.second[pixit].row()<<", ph offset "<<(int)dacit_max->first<<", ph scale "<<(int)dacit_max->second.first<<", max ph "<<(int)dacit_max->second.second[pixit].value();
//      LOG(logDEBUG)<<"dacit_min: pixel "<<(int)dacit_min->second.second[pixit].roc()<<", "<<(int)dacit_min->second.second[pixit].column()<<", "<<(int)dacit_min->second.second[pixit].row()<<", ph offset "<<(int)dacit_min->first<<", ph scale "<<(int)dacit_min->second.first<<", min ph "<<(int)dacit_min->second.second[pixit].value();
//    }
//    dacit_min++;
//  }
  int pixsize_max=0;
  int pixsize_min=0;
  LOG(logDEBUG)<<"InsideRange() subtest";
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    for(dacit_max = dacdac_max.begin(); dacit_max != dacdac_max.end(); dacit_max++){
      if(dacit_max->first == po_opt[rocIds[roc_it]]){
	for(dacit_min = dacdac_min.begin(); dacit_min != dacdac_min.end(); dacit_min++)
	  if(dacit_min->first == po_opt[rocIds[roc_it]] && dacit_min->second.first == dacit_max->second.first){
    	    pixsize_max = dacit_max->second.second.size();
	    pixsize_min = dacit_min->second.second.size();
    	    for(int pix=0; pix < pixsize_max; pix++){
	      if((dacit_max->second.second[pix].roc()!=rocIds[roc_it] || dacit_min->second.second[pix].roc()!=rocIds[roc_it]) && dacit_max->second.second[pix].roc() != dacit_min->second.second[pix].roc()){
		//	LOG(logDEBUG)<<"####//// this time roc ids DO NOT match "<<(int)dacit_min->second.second[pix].roc()<<" "<<(int)dacit_max->second.second[pix].roc()<<" "<<(int)rocIds[roc_it]<<" "<<(int)roc_it<<" ////####";
		continue;
	      }
	      //LOG(logDEBUG)<<"####//// this time roc ids match "<<(int)rocIds[roc_it]<<" "<<(int)roc_it<<" ////####";
	      maxPh=dacit_max->second.second[pix].value();
	      minPh=dacit_min->second.second[pix].value();
	      if(dacit_max->second.second[pix].roc() != dacit_min->second.second[pix].roc()){
		LOG(logDEBUG) << "InsideRangePH: ROC ids do not correspond";
	      }
	      lowEd = (minPh > safetyMargin);
	      upEd = (maxPh < 255 - safetyMargin);
	      lowEd_dist = abs(minPh - safetyMargin);
	      upEd_dist = abs(maxPh - (255 - safetyMargin));
	      dist = (upEd_dist > lowEd_dist ) ? (upEd_dist) : (lowEd_dist);
	      if(dist < bestDist[dacit_max->second.second[pix].roc()] && upEd && lowEd){
		LOG(logDEBUG)<<"New distance "<<dist<<" is smaller than previous bestDist "<<bestDist[dacit_max->second.second[pix].roc()]<<" and edges are ok, so... ";
		ps_opt[dacit_max->second.second[pix].roc()] = dacit_max->second.first;
		bestDist[dacit_max->second.second[pix].roc()]=dist;
		LOG(logDEBUG)<<"... new bestDist is "<<bestDist[dacit_max->second.second[pix].roc()]<<" for ps_opt = "<<ps_opt[dacit_max->second.second[pix].roc()];
	      }
	    }
	  }
      }
    }
  }
  
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    LOG(logDEBUG)<<"opt step 1: po fixed to"<<po_opt[rocIds[roc_it]]<<" and scale adjusted to "<<ps_opt[rocIds[roc_it]]<<" for ROC "<<(int)rocIds[roc_it]<<", with distance "<<bestDist[rocIds[roc_it]];
  }
  return ps_opt;
}

map<uint8_t, int> PixTestPhOptimization::CentrePhRange(map<uint8_t, int> &po_opt, map<uint8_t, int> &ps_opt,  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_max,   std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_min){
  //centring PH curve adjusting phoffset   
  LOG(logDEBUG)<<"Welcome to CentrePhRange()"; 
  int maxPh(0);
  int minPh(0);
  int dist = 255;
  map<uint8_t, int> bestDist;
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    bestDist[rocIds[roc_it]] = 255;
  }
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_max = dacdac_max.begin();
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_min = dacdac_min.begin();
  int pixsize_max=0;
  int pixsize_min=0;
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    dist = 255;
    for(dacit_max = dacdac_max.begin(); dacit_max != dacdac_max.end(); dacit_max++){
      if(dacit_max->second.first != ps_opt[rocIds[roc_it]])continue;
      for(dacit_min = dacdac_min.begin(); dacit_min != dacdac_min.end(); dacit_min++){
	if(dacit_min->second.first == ps_opt[rocIds[roc_it]] && dacit_min->first == dacit_max->first){
	  pixsize_max = dacit_max->second.second.size();
	  pixsize_min = dacit_min->second.second.size();
	  for(int pix=0; pix < pixsize_max; pix++){
	    if(dacit_max->second.second[pix].roc()!=rocIds[roc_it] || dacit_min->second.second[pix].roc()!=rocIds[roc_it]) continue;
	    maxPh=dacit_max->second.second[pix].value();
	    minPh=dacit_min->second.second[pix].value();
	    if(dacit_max->second.second[pix].roc() != dacit_min->second.second[pix].roc()){
	      LOG(logDEBUG) << "CentrePhRange: ROC ids do not correspond";
	    }
	    dist = abs(minPh - (255 - maxPh));
	    if (dist < bestDist[dacit_max->second.second[pix].roc()]){
	      po_opt[dacit_max->second.second[pix].roc()] = dacit_max->first;
	      bestDist[dacit_max->second.second[pix].roc()] = dist;
	    } 
	  }
	}
      }
    }
  }
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    LOG(logDEBUG)<<"opt centring step: po "<<po_opt[rocIds[roc_it]]<<" and scale "<<ps_opt[rocIds[roc_it]]<<", with distance "<<bestDist[rocIds[roc_it]]<<" on ROC "<<(int)rocIds[roc_it];
  }
  return po_opt;
}

map<uint8_t, int> PixTestPhOptimization::StretchPH(map<uint8_t, int> &po_opt, map<uint8_t, int> &ps_opt,  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_max,   std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_min){
  //stretching PH curve to exploit the full ADC range, adjusting phscale             
  int maxPh(0);
  int minPh(0);
  bool lowEd=false, upEd=false;
  int upEd_dist=255, lowEd_dist=255;
  int safetyMarginUp = fSafetyMarginUp;
  int safetyMarginLow = fSafetyMarginLow;
  LOG(logDEBUG)<<"safety margin for stretching set to "<<fSafetyMarginLow<<" (lower edge) and "<<fSafetyMarginUp<<"(upper edge)";
  int dist = 255;
  map<uint8_t, int> bestDist;
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs();
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    bestDist[rocIds[roc_it]] = 255;
  }
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_max = dacdac_max.begin();
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_min = dacdac_min.begin();
  while(dacit_max != dacdac_max.end() || dacit_min != dacdac_min.end()){
    for(unsigned int pix=0; pix < dacit_min->second.second.size() && pix < dacit_max->second.second.size(); pix++){
      if(dacit_max->first == po_opt[dacit_max->second.second[pix].roc()] && dacit_min->first == po_opt[dacit_min->second.second[pix].roc()]){
	maxPh=dacit_max->second.second[pix].value();
	minPh=dacit_min->second.second[pix].value();
	if(dacit_max->second.second[pix].roc() != dacit_min->second.second[pix].roc()){
	  LOG(logDEBUG) << "CentrePhRange: ROC ids do not correspond";
	}
	lowEd = (minPh > safetyMarginLow);
	upEd = (maxPh < 255 - safetyMarginUp);
	upEd_dist = abs(maxPh - (255 - safetyMarginUp));
	lowEd_dist = abs(minPh - safetyMarginLow);
	dist = (upEd_dist < lowEd_dist ) ? (upEd_dist) : (lowEd_dist);
	if(dist < bestDist[dacit_max->second.second[pix].roc()] && lowEd && upEd){
	  ps_opt[dacit_max->second.second[pix].roc()] = dacit_max->second.first;
	  bestDist[dacit_max->second.second[pix].roc()]=dist;
	}
      }
    }
    dacit_max++;
    dacit_min++;
  }
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    LOG(logDEBUG)<<"opt final step: po fixed to"<<po_opt[rocIds[roc_it]]<<" and scale adjusted to "<<ps_opt[rocIds[roc_it]]<<", with distance "<<bestDist[rocIds[roc_it]]<<" on ROC "<<(int)rocIds[roc_it];
  }
  return ps_opt;
}

void PixTestPhOptimization::DrawPhMaps(std::map<int, int> &minVcal, std::vector<std::pair<uint8_t, std::pair<int, int> > > &badPixels){
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs();
  string name, title;
  //draw PH maps and extract validation distributions from them
    fApi->_dut->testAllPixels(true);
    fApi->_dut->maskAllPixels(false);
    
    std::vector<pxar::pixel> result_map;
    map<int, TH2D* > h2_PhMaps;
    map<int, TH1D* > h1_PhMaps;
    TH2D* h2_PhMap(0);
    TH1D* h1_PhMap(0);
    fApi->setDAC("ctrlreg",4);
    fApi->setDAC("vcal",255);
    //pulseheight map at vcal=255
    result_map = fApi->getPulseheightMap(0,10);
    //unpacking data from map and filling one histo per ROC
    for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
      name  = Form("PH_mapHiVcal_C%d", rocIds[roc_it]);
      title = Form("PH_mapHiVcal_C%d", rocIds[roc_it]);
      h2_PhMap = bookTH2D(name, name, 52, 0., 52., 80, 0., 80.);
      h2_PhMaps.insert(make_pair(rocIds[roc_it], h2_PhMap));
      fHistList.push_back(h2_PhMaps[rocIds[roc_it]]);  
      fHistOptions.insert( make_pair(h2_PhMaps[rocIds[roc_it]], "colz")  ); 
    }
    for (unsigned int i = 0; i < result_map.size(); ++i) {
      h2_PhMaps[ (int) result_map[i].roc() ]->Fill( result_map[i].column(), result_map[i].row(), result_map[i].value());
    }
    //adjust z axis range and extract ph distribution
    for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
      h2_PhMaps[ (int) rocIds[roc_it] ]->GetZaxis()->SetRangeUser(h2_PhMaps[ (int) rocIds[roc_it] ]->GetMinimum(), 255. );
      h1_PhMap = distribution( h2_PhMaps[ (int) rocIds[roc_it] ], 255, 0., 255.);
      h1_PhMaps.insert( make_pair(rocIds[roc_it], h1_PhMap) );
      fHistList.push_back(h1_PhMaps[rocIds[roc_it]]);  
    }
    
    //PH map for lower vcal sampling point
    for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
      fApi->setDAC("ctrlreg",0);
      fApi->setDAC("vcal",minVcal[roc_it]+10, rocIds[roc_it] );
    }
    map<int, TH2D* > h2_PhMapsMin;
    map<int, TH1D* > h1_PhMapsMin;
    //phmap
    result_map = fApi->getPulseheightMap(0,10);
    //unpacking data from map and filling one histo per ROC
    for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
      name  = Form("PH_mapLowVcal_C%d", rocIds[roc_it]);
      title = Form("PH_mapLowVcal_C%d", rocIds[roc_it]);
      h2_PhMap = bookTH2D(name, name, 52, 0., 52., 80, 0., 80.);
      h2_PhMapsMin.insert(make_pair(rocIds[roc_it], h2_PhMap));
      fHistList.push_back(h2_PhMapsMin[rocIds[roc_it]]);  
      fHistOptions.insert( make_pair(h2_PhMapsMin[rocIds[roc_it]], "colz")  ); 
    }
    for (unsigned int i = 0; i < result_map.size(); ++i) {
      h2_PhMapsMin[ (int) result_map[i].roc() ]->Fill( result_map[i].column(), result_map[i].row(), result_map[i].value());
    }
    //removing blacklisted pixels (bincontent = 0) from histos
    for(std::vector<std::pair<uint8_t, pair<int, int> > >::iterator bad_it = badPixels.begin(); bad_it != badPixels.end(); bad_it++){
      h2_PhMaps[ bad_it->first] ->SetBinContent( h2_PhMaps[ bad_it->first] ->FindFixBin(bad_it->second.first, bad_it->second.second), 0);
      h2_PhMapsMin[ bad_it->first] ->SetBinContent( h2_PhMapsMin[ bad_it->first] ->FindFixBin(bad_it->second.first, bad_it->second.second), 0);
    }
    //adjust z axis range and extract ph distribution
    for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
      name  = Form("PH_distr_LowVcal_C%d", rocIds[roc_it]);
      title = Form("PH_distr_LowVcal_C%d", rocIds[roc_it]);
      h2_PhMapsMin[ (int) rocIds[roc_it] ]->GetZaxis()->SetRangeUser(0., h2_PhMapsMin[ (int) rocIds[roc_it] ]->GetMaximum() );
      h1_PhMap = distribution( h2_PhMapsMin[ (int) rocIds[roc_it] ], 255, 0., 255.);
      h1_PhMapsMin.insert( make_pair(rocIds[roc_it], h1_PhMap) );
      fHistList.push_back(h1_PhMapsMin[rocIds[roc_it]]);  
    }
}

void PixTestPhOptimization::DrawPhCurves(map<int, pxar::pixel > &maxpixels, map<int, pxar::pixel > &minpixels, std::map<uint8_t, int> &po_opt, std::map<uint8_t, int> &ps_opt){
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs();
  string name, title;
  TH1D *h1(0); 
    for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
      fApi->setDAC("ctrlreg",4);
    }
  vector<pair<uint8_t, vector<pixel> > > results;
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
  //draw PH curve for max and min pixel on every ROC
    results.clear();
    name  = Form("PH_c%d_r%d_C%d", maxpixels[roc_it].column(), maxpixels[roc_it].row(), rocIds[roc_it]);
    title = Form("PH_c%d_r%d_C%d, phscale = %d, phoffset = %d, maxpixel", maxpixels[roc_it].column(), maxpixels[roc_it].row(), rocIds[roc_it], ps_opt[rocIds[roc_it]], po_opt[rocIds[roc_it]]);
    h1 = bookTH1D(name, name, 256, 0., 256.);
    fApi->_dut->testAllPixels(false);
    fApi->_dut->maskAllPixels(true);
    fApi->_dut->testPixel(maxpixels[roc_it].column(), maxpixels[roc_it].row(), true, roc_it);
    fApi->_dut->maskPixel(maxpixels[roc_it].column(), maxpixels[roc_it].row(), false, roc_it);
    int  cnt = 0; 
    bool done = false;
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
    
    for (unsigned int i = 0; i < results.size(); ++i) {
      pair<uint8_t, vector<pixel> > v = results[i];
      int idac = v.first; 
      vector<pixel> vpix = v.second;
      for (unsigned int ipix = 0; ipix < vpix.size(); ++ipix) {
	h1->Fill(idac, vpix[ipix].value());
      }
    }
    h1->SetMinimum(0);
    setTitles(h1, title.c_str(), "average PH");
    fHistList.push_back(h1);  
    
    results.clear();
    
    name  = Form("PH_c%d_r%d_C%d", minpixels[roc_it].column(), minpixels[roc_it].row(), rocIds[roc_it]);
    title = Form("PH_c%d_r%d_C%d, phscale = %d, phoffset = %d, minpixel", minpixels[roc_it].column(), minpixels[roc_it].row(), rocIds[roc_it], ps_opt[rocIds[roc_it]], po_opt[rocIds[roc_it]]);
    h1 = bookTH1D(name, name, 256, 0., 256.);
    fApi->_dut->testAllPixels(false);
    fApi->_dut->maskAllPixels(true);
    fApi->_dut->testPixel(minpixels[roc_it].column(), minpixels[roc_it].row(), true, roc_it);
    fApi->_dut->maskPixel(minpixels[roc_it].column(), minpixels[roc_it].row(), false, roc_it);
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
    
    for (unsigned int i = 0; i < results.size(); ++i) {
      pair<uint8_t, vector<pixel> > v = results[i];
      int idac = v.first; 
      vector<pixel> vpix = v.second;
      for (unsigned int ipix = 0; ipix < vpix.size(); ++ipix) {
	h1->Fill(idac, vpix[ipix].value());
      }
    }
    h1->SetMinimum(0);
    setTitles(h1, title.c_str(), "average PH");
    fHistList.push_back(h1);
  }
}

void PixTestPhOptimization::MaxPhVsDacDac(std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_max, map<int, pxar::pixel> maxpixels){
  
  fApi->_dut->testAllPixels(false);
  fApi->_dut->maskAllPixels(true);
  for(std::map<int, pxar::pixel>::iterator maxp_it = maxpixels.begin(); maxp_it != maxpixels.end(); maxp_it++){
    fApi->_dut->testPixel(maxp_it->second.column(),maxp_it->second.row(),true, getIdxFromId(maxp_it->second.roc()));
    fApi->_dut->maskPixel(maxp_it->second.column(),maxp_it->second.row(),false, getIdxFromId(maxp_it->second.roc()));
  } 
  fApi->setDAC("vcal",255);
  fApi->setDAC("ctrlreg",4);
  
  //scanning through offset and scale for max pixel (or randpixel)
  int cnt = 0; 
  bool done = false;
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
  //std::map<uint8_t, std::map<std::pair<uint8_t, uint8_t>, pxar::pixel > >  dacdacmax_map;
  //std::map<std::pair<uint8_t ,uint8_t >, pxar::pixel > tempMap;
  for(  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_max = dacdac_max.begin(); dacit_max < dacdac_max.end(); dacit_max+=10){
    LOG(logDEBUG)<<"size "<<(int)(dacdac_max.end() - dacdac_max.begin())<<"pos "<<(int)(dacdac_max.end() - dacit_max)<<" sizepix "<<(int)dacit_max->second.second.size();
 //   for(int ipix=0; ipix < dacit_max->second.second.size(); ipix++ ){
 //     LOG(logDEBUG)<<"roc "<< (int)dacit_max->second.second[ipix].roc() <<" col "<<(int)dacit_max->second.second[ipix].column() << " row "<<(int)dacit_max->second.second[ipix].row();
 //   }
  // for(int ipix=0; ipix < dacit_max->second.second.size(); ipix++ ){
  //   tempMap[ make_pair(dacit_max->first, dacit_max->second.first) ] = dacit_max->second.second.at(ipix);
  //   dacdacmax_map[dacit_max->second.second[ipix].roc()] = tempMap;
  //   tempMap.clear();
  //  }
  }
  
  
}

void PixTestPhOptimization::MinPhVsDacDac(std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_min, map<int, pxar::pixel> minpixels, std::map<int, int> &minVcal){
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs();
  fApi->_dut->testAllPixels(false);
  fApi->_dut->maskAllPixels(true);
  unsigned int NRocs = rocIds.size();
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > dacdac_min_part;

  for(std::map<int, pxar::pixel>::iterator minp_it = minpixels.begin(); minp_it != minpixels.end(); minp_it++){
    for(unsigned int roc_kt = 0; roc_kt < NRocs; roc_kt++){
      fApi->_dut->setROCEnable(roc_kt, false);
    }
    fApi->_dut->setROCEnable(minp_it->first, true);
    fApi->_dut->testPixel(minp_it->second.column(),minp_it->second.row(),true);
    fApi->_dut->maskPixel(minp_it->second.column(),minp_it->second.row(),false);
  
    fApi->setDAC("ctrlreg",0);
    fApi->setDAC("vcal",minVcal[minp_it->first]+10, rocIds[minp_it->first]);
    
  //scanning through offset and scale for min pixel (or same randpixel)
    int cnt = 0; 
    int done = false;
    while (!done) {
      try {
	dacdac_min_part = fApi->getPulseheightVsDACDAC("phoffset",0,255,"phscale",0,255,0,10);
	done = true;
      } catch(pxarException &e) {
	LOG(logCRITICAL) << "pXar execption: "<< e.what(); 
	++cnt;
    }
      done = (cnt>5) || done;
    }
    int dacdacsize_old = dacdac_min.size();
    dacdac_min.resize(dacdac_min.size() + dacdac_min_part.size());
    for(int idacdacpart = 0; idacdacpart< dacdac_min_part.size(); idacdacpart++){
      dacdac_min[dacdacsize_old + idacdacpart] = dacdac_min_part[idacdacpart];
    }
    dacdac_min_part.clear();
  }

  //  std::map<uint8_t, std::map<std::pair<uint8_t, uint8_t>, pxar::pixel > >  dacdacmin_map;
  //std::map<std::pair<uint8_t ,uint8_t >, pxar::pixel > tempMap;
  for(  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_min = dacdac_min.begin(); dacit_min < dacdac_min.end(); dacit_min+=1000){
    LOG(logDEBUG)<<"dacdac_min: size "<<(int)(dacdac_min.end() - dacdac_min.begin())<<"pos "<<(int)(dacdac_min.end() - dacit_min)<<" sizepix "<<(int)dacit_min->second.second.size();
    for(int ipix=0; ipix < dacit_min->second.second.size(); ipix++ ){
      LOG(logDEBUG)<<"roc "<< (int)dacit_min->second.second[ipix].roc() <<" col "<<(int)dacit_min->second.second[ipix].column() << " row "<<(int)dacit_min->second.second[ipix].row();
    }
//   for(int ipix=0; ipix < dacit_min->second.second.size(); ipix++ ){
//     tempMap[ make_pair(dacit_min->first, dacit_min->second.first) ] = dacit_min->second.second.at(ipix);
//     dacdacmin_map[dacit_min->second.second[ipix].roc()] = tempMap;
//     tempMap.clear();
//    }
  }

}

void PixTestPhOptimization::SetMinThr(){
  fMinThr=0;
  string trimfile = fPixSetup->getConfigParameters()->getTrimParameterFileName() + fPixSetup->getConfigParameters()->getTrimVcalSufix();
  trimfile.erase(0,14); //erases 'trimParamerers' from the name of the file
  if(0!=(trimfile.compare(""))){
    fMinThr = atoi(trimfile.c_str());
  }
  else{
    LOG(logINFO)<<"***::: The test requires a TRIMMED module, but no TrimParameterFile is loaded :::***";
    LOG(logINFO)<<"Vcal lower sample point will be set to 40";
    fMinThr=40;
  }
}


/*map<uint8_t, int> PixTestPhOptimization::InsideRangePH_new(map<uint8_t,int> &po_opt,  std::map<uint8_t, std::map<std::pair<uint8_t, uint8_t>, pxar::pixel > >  dacdacmax_map, std::map<uint8_t, std::map<std::pair<uint8_t, uint8_t>, pxar::pixel > >  dacdacmin_map){
  //adjusting phscale so that the PH curve is fully inside the ADC range
  map<uint8_t, int> ps_opt;
  int maxPh(0);
  int minPh(0);
  bool lowEd=false, upEd=false;
  int upEd_dist=255, lowEd_dist=255;
  int safetyMargin = 50;
  int dist = 255;
  map<uint8_t, int> bestDist;
  LOG(logDEBUG) << "dacdac at max vcal has size "<<dacdac_max.size()<<endl;
  LOG(logDEBUG) << "dacdac at min vcal has size "<<dacdac_min.size()<<endl;
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs(); 
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    bestDist[rocIds[roc_it]] = 255;
    LOG(logDEBUG)<<"Bestdist at roc_it "<<roc_it<<" initialized with "<<bestDist[roc_it]<<" "<<bestDist[rocIds[roc_it]];
    ps_opt[rocIds[roc_it]] = 999;
  }
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_max = dacdac_max.begin();
  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > >::iterator dacit_min = dacdac_min.begin();
//  for(dacit_max = dacdac_max.begin(); dacit_max != dacdac_max.end(); dacit_max++){
//    for(int pixit = 0; pixit < dacit_max->second.second.size(); pixit++){
//      LOG(logDEBUG)<<"dacit_max: pixel "<<(int)dacit_max->second.second[pixit].roc()<<", "<<(int)dacit_max->second.second[pixit].column()<<", "<<(int)dacit_max->second.second[pixit].row()<<", ph offset "<<(int)dacit_max->first<<", ph scale "<<(int)dacit_max->second.first<<", max ph "<<(int)dacit_max->second.second[pixit].value();
//      LOG(logDEBUG)<<"dacit_min: pixel "<<(int)dacit_min->second.second[pixit].roc()<<", "<<(int)dacit_min->second.second[pixit].column()<<", "<<(int)dacit_min->second.second[pixit].row()<<", ph offset "<<(int)dacit_min->first<<", ph scale "<<(int)dacit_min->second.first<<", min ph "<<(int)dacit_min->second.second[pixit].value();
//    }
//    dacit_min++;
//  }
  int pixsize_max=0;
  int pixsize_min=0;
  LOG(logDEBUG)<<"InsideRange() subtest";
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    for(dacit_max = dacdac_max.begin(); dacit_max != dacdac_max.end(); dacit_max++){
      if(dacit_max->first == po_opt[rocIds[roc_it]]){
	for(dacit_min = dacdac_min.begin(); dacit_min != dacdac_min.end(); dacit_min++)
	  if(dacit_min->first == po_opt[rocIds[roc_it]] && dacit_min->second.first == dacit_max->second.first){
    	    pixsize_max = dacit_max->second.second.size();
	    pixsize_min = dacit_min->second.second.size();
    	    for(int pix=0; pix < pixsize_max; pix++){
	      if(dacit_max->second.second[pix].roc()!=rocIds[roc_it] || dacit_min->second.second[pix].roc()!=rocIds[roc_it]){
		LOG(logDEBUG)<<"####//// this time roc ids DO NOT match "<<(int)dacit_min->second.second[pix].roc()<<" "<<(int)dacit_max->second.second[pix].roc()<<" "<<(int)rocIds[roc_it]<<" "<<(int)roc_it<<" ////####";
		continue;
	      }
	      LOG(logDEBUG)<<"####//// this time roc ids match "<<(int)rocIds[roc_it]<<" "<<(int)roc_it<<" ////####";
	      maxPh=dacit_max->second.second[pix].value();
	      minPh=dacit_min->second.second[pix].value();
	      if(dacit_max->second.second[pix].roc() != dacit_min->second.second[pix].roc()){
		LOG(logDEBUG) << "InsideRangePH: ROC ids do not correspond";
	      }
	      lowEd = (minPh > safetyMargin);
	      upEd = (maxPh < 255 - safetyMargin);
	      lowEd_dist = abs(minPh - safetyMargin);
	      upEd_dist = abs(maxPh - (255 - safetyMargin));
	      dist = (upEd_dist > lowEd_dist ) ? (upEd_dist) : (lowEd_dist);
	      if(dist < bestDist[dacit_max->second.second[pix].roc()] && upEd && lowEd){
		LOG(logDEBUG)<<"New distance "<<dist<<" is smaller than previous bestDist "<<bestDist[dacit_max->second.second[pix].roc()]<<" and edges are ok, so... ";
		ps_opt[dacit_max->second.second[pix].roc()] = dacit_max->second.first;
		bestDist[dacit_max->second.second[pix].roc()]=dist;
		LOG(logDEBUG)<<"... new bestDist is "<<bestDist[dacit_max->second.second[pix].roc()]<<" for ps_opt = "<<ps_opt[dacit_max->second.second[pix].roc()];
	      }
	    }
	  }
      }
    }
  }
  
  for(unsigned int roc_it = 0; roc_it < rocIds.size(); roc_it++){
    LOG(logDEBUG)<<"opt step 1: po fixed to"<<po_opt[rocIds[roc_it]]<<" and scale adjusted to "<<ps_opt[rocIds[roc_it]]<<" for ROC "<<(int)rocIds[roc_it]<<", with distance "<<bestDist[rocIds[roc_it]];
  }
  return ps_opt;
}*/
