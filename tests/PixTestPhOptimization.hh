#ifndef PixTestPhOptimization_H
#define PixTestPhOptimization_H

#include "PixTest.hh"


class DLLEXPORT PixTestPhOptimization: public PixTest {
public:
  PixTestPhOptimization(PixSetup *, std::string);
  PixTestPhOptimization();
  virtual ~PixTestPhOptimization();
  virtual bool setParameter(std::string parName, std::string sval); 
  void init(); 
  void bookHist(std::string); 
  void BlacklistPixels(std::vector<std::pair<uint8_t, std::pair<int, int> > > &badPixels, int aliveTrig);  
  pxar::pixel* RandomPixel(std::vector<std::pair<uint8_t, std::pair<int, int> > > &badPixels, uint8_t iroc);
  void GetMaxPhPixel(std::vector<std::pair<int, pxar::pixel> > &maxpixel, std::vector<std::pair<uint8_t, std::pair<int, int> > > &badPixels);
  void GetMinPixel(std::vector<std::pair<int, pxar::pixel> > &minpixel, std::vector<pxar::pixel> &thrmap, std::vector<std::pair<uint8_t, std::pair<int, int> > > &badPixels);
  int InsideRangePH(int po_opt,  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_max,   std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_min);
  int CentrePhRange(int po_opt, int ps_opt,  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_max,   std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_min);
  int StretchPH(int po_opt, int ps_opt_in,  std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_max,   std::vector< std::pair<uint8_t, std::pair<uint8_t, std::vector<pxar::pixel> > > > &dacdac_min);
  void DynamicRange();
  void doTest(); 

private:

  int     fParNtrig; 
  std::string fParDAC; 
  int     fParDacVal;
  bool fFlagSinglePix;
  int fSafetyMargin;

  ClassDef(PixTestPhOptimization, 1)

};
#endif
