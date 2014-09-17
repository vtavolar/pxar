#include "TH2F.h"
#include "TFile.h"
#include "TH1F.h"
#include "iostream"
#include "string"

void PhOptAnalysis(){

  TH1F* h_max = new TH1F("h_max", "h_max", 80, 0., 80.);
  TH1F* h_min = new TH1F("h_min", "h_min", 80, 0., 80.);
  TH2D* h2_min = new TH2D();
  TH2D* h2_max = new TH2D();
    

  for(int i=0; i<4; i++){
    string filename = Form("pxarWR%d.root", i+1); 
    std::cout<<filename<<std::endl;
    TFile* f_in = new TFile(filename.c_str(), "READ"); 
    std::cout<<"Reading file "<<filename<<" for input"<<std::endl; 
    f_in->ls();
    f_in->cd("PhOptimization");
    gDirectory->pwd();
    f_in->ls();
    for(int sm = 0; sm<14; sm++){ 
      string histoname = Form("PhOptimization/PhMapMinVcal_SM%d_V0", 10+sm*5);
      std::cout<<"histoname: "<<histoname<<std::endl;
      h2_min = (TH2D*)f_in->Get(histoname.c_str());
      int nbins;
      nbins=h2_min->GetSize();
      std::cout << "Number of bins h2_min: "<<nbins<<std::endl;
      int cnt=0;
      double val=0;
      for(int j=1; j<=nbins; j++){
	val = h2_min->GetBinContent(j);
	if(val==0){
	  cnt++;
	}
      }
      h_min->Fill(10+sm*5, cnt);



      histoname = Form("PhOptimization/PhMapMaxVcal_SM%d_V0", 10+sm*5);
      h2_max = (TH2D*)f_in->Get(histoname.c_str());
      nbins=h2_max->GetSize();
      std::cout << "Number of bins h2_max: "<<nbins<<std::endl;
      cnt=0;
      val=0;
      for(int j=1; j<=nbins; j++){
	val = h2_max->GetBinContent(j);
	if( val==255){
	  cnt++;
	}
      }
      h_max->Fill(10+sm*5, cnt);

      
    }
    f_in->Close();
  }

  h_min->Scale(0.2);
  h_max->Scale(0.2);
  
  h_min->SetMarkerStyle(20);
  h_min->SetMarkerColor(kBlue);

  h_max->SetMarkerStyle(20);
  h_max->SetMarkerColor(kRed);


  TCanvas* c1 = new TCanvas("c1", "", 1);
  h_min->Draw("pe");
  TCanvas* c2 = new TCanvas("c2", "", 1);
  h_max->Draw("pe");

  h_min->Print("all");
  h_max->Print("all");

  c1->SaveAs("minPhOutOfRange_WR.jpg");
  c2->SaveAs("maxPhOutOfRange_WR.jpg");

}
