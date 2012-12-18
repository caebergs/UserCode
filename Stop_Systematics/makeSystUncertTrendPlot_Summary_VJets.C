#include <TROOT.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

#include <math.h>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;


void makeSystUncertTrendPlot_Summary_VJets(const int UseCase){

  gStyle->SetOptFit(0);

  const int NbOfUseCase = 4;
  std::string teffname[NbOfUseCase] = {"LowDM", "IntDM1", "IntDM2", "HighDM"};
// Use case : 
  // 0 : low DM
  // 1 : Int DM 1
  // 2 : Int DM 2
  // 3 : High DM
//  const int UseCase = 3;
  const UInt_t channel = 3 ;
  std::string label = "channel" ;
  label += (channel==1?"1":(channel==2?"2":(channel==3?"3":"X")));

  const int nsyst = 4;

  double rel_syst_uncert[NbOfUseCase][nsyst][2];
  rel_syst_uncert[0][0][0]=-0.1424080; // low DM, jes -
  rel_syst_uncert[0][0][1]= 0.0482035; // low DM, jes +
  rel_syst_uncert[1][0][0]=-0.0011384; // int DM1, jes -
  rel_syst_uncert[1][0][1]= 0.0238592; // int DM1, jes +
  rel_syst_uncert[2][0][0]= 0.0232405; // int DM2, jes -
  rel_syst_uncert[2][0][1]= 0.0075060; // int DM2, jes +
  rel_syst_uncert[3][0][0]= 0.0137523; // high DM, jes -
  rel_syst_uncert[3][0][1]= 0.0116015; // high DM, jes +

  rel_syst_uncert[0][1][0]=-0.122193; // low DM, matching -
  rel_syst_uncert[0][1][1]=-0.057742; // low DM, matching +
  rel_syst_uncert[1][1][0]=-0.291199; // int DM1, matching -
  rel_syst_uncert[1][1][1]=-0.201351; // int DM1, matching +
  rel_syst_uncert[2][1][0]=-0.267173; // int DM2, matching -
  rel_syst_uncert[2][1][1]=-0.150011; // int DM2, matching +
  rel_syst_uncert[3][1][0]=-0.268341; // high DM, matching -
  rel_syst_uncert[3][1][1]=-0.149720; // high DM, matching +

  rel_syst_uncert[0][2][0]=-0.130539; // low DM,  scale -
  rel_syst_uncert[0][2][1]=-0.183195; // low DM,  scale +
  rel_syst_uncert[1][2][0]=-0.027962; // int DM1, scale -
  rel_syst_uncert[1][2][1]=-0.272830; // int DM1, scale +
  rel_syst_uncert[2][2][0]=-0.096214; // int DM2, scale -
  rel_syst_uncert[2][2][1]=-0.317220; // int DM2, scale +
  rel_syst_uncert[3][2][0]=-0.174006; // high DM, scale -
  rel_syst_uncert[3][2][1]=-0.413494; // high DM, scale +

  rel_syst_uncert[0][3][0]=-0.130539; // low DM,  5D Rew.
  rel_syst_uncert[0][3][1]=-0.183195; // low DM,  5D Rew.
  rel_syst_uncert[1][3][0]=-0.027962; // int DM1, 5D Rew.
  rel_syst_uncert[1][3][1]=-0.272830; // int DM1, 5D Rew.
  rel_syst_uncert[2][3][0]=-0.096214; // int DM2, 5D Rew.
  rel_syst_uncert[2][3][1]=-0.317220; // int DM2, 5D Rew.
  rel_syst_uncert[3][3][0]=-0.174006; // high DM, 5D Rew.
  rel_syst_uncert[3][3][1]=-0.413494; // high DM, 5D Rew.

  const double RV[NbOfUseCase] = {0.001789,0.001605,0.007944,0.010535};

  string binlabel[nsyst*2-1+1] = {"JES (-)","JES (+)","ME-PS matching (-)","ME-PS matching (+)","Scale Q^{2} (-)","Scale Q^{2} (+)","5D Rew.","Tot. syst."};
  char name[100];
  
  TFile *outputfile = new TFile((std::string()+"RV_trendPlotSummary_"+label+"_"+teffname[UseCase]+"_Output.root").c_str(),"RECREATE");
  TH1F *h  = new TH1F("h","",nsyst*2-1+1,0,nsyst*2+1);
  TH1F *h1 = new TH1F("h1","",nsyst*2-1+1,0,nsyst*2+1);
 
  TAxis *Xaxis = h->GetXaxis();
  TAxis *Yaxis = h->GetYaxis();

  double CombUncer_low  = 0, CombUncer_error_low   = 0;
  double CombUncer_high = 0, CombUncer_error_high  = 0;

  for (int i=0;i<nsyst;i++){
    Xaxis->SetBinLabel((i*2)+1,binlabel[(i*2)].c_str());
    if (i==3) {
      if (rel_syst_uncert[UseCase][i][0]<0) rel_syst_uncert[UseCase][i][0] *= -1.;
      if (rel_syst_uncert[UseCase][i][1]>0) rel_syst_uncert[UseCase][i][1] *= -1.;
    } else {
      Xaxis->SetBinLabel((i*2)+2,binlabel[(i*2)+1].c_str());
    }
    
    if(rel_syst_uncert[UseCase][i][0]>0 && rel_syst_uncert[UseCase][i][1]<0){
      h->SetBinContent((i*2)+1,rel_syst_uncert[UseCase][i][0]);
      if (i==3) {
        h1->SetBinContent((i*2)+1,rel_syst_uncert[UseCase][i][1]);
      } else {
        h1->SetBinContent((i*2)+2,rel_syst_uncert[UseCase][i][1]);
      }
    //h->SetBinError(i+1,XXX);
      cout<<"Bin+ "<<i<<", Content : "<<rel_syst_uncert[UseCase][i][0]<<endl;
      cout<<"Bin- "<<i<<", Content : "<<rel_syst_uncert[UseCase][i][1]<<endl;
      CombUncer_high += pow(rel_syst_uncert[UseCase][i][0],2);
      CombUncer_low  += pow(rel_syst_uncert[UseCase][i][1],2);
    }
    else if(rel_syst_uncert[UseCase][i][0]<0 && rel_syst_uncert[UseCase][i][1]>0){
      h->SetBinContent((i*2)+1,rel_syst_uncert[UseCase][i][1]);
      h1->SetBinContent((i*2)+2,rel_syst_uncert[UseCase][i][0]);
      //h->SetBinError(i+1,XXX);
      cout<<"Bin+ "<<i<<", Content : "<<rel_syst_uncert[UseCase][i][1]<<endl;
      cout<<"Bin- "<<i<<", Content : "<<rel_syst_uncert[UseCase][i][0]<<endl;
      CombUncer_high += pow(rel_syst_uncert[UseCase][i][1],2);
      CombUncer_low  += pow(rel_syst_uncert[UseCase][i][0],2);
    }
    else if(rel_syst_uncert[UseCase][i][0]>0 && rel_syst_uncert[UseCase][i][1]>0){
      h->SetBinContent((i*2)+1,rel_syst_uncert[UseCase][i][0]);
      h->SetBinContent((i*2)+2,rel_syst_uncert[UseCase][i][1]);
      if(rel_syst_uncert[UseCase][i][0]>rel_syst_uncert[UseCase][i][1]) CombUncer_high += pow(rel_syst_uncert[UseCase][i][0],2);
      else                                                              CombUncer_high += pow(rel_syst_uncert[UseCase][i][1],2);
    }
    else if(rel_syst_uncert[UseCase][i][0]<0 && rel_syst_uncert[UseCase][i][1]<0){
      h1->SetBinContent((i*2)+1,rel_syst_uncert[UseCase][i][0]);
      h1->SetBinContent((i*2)+2,rel_syst_uncert[UseCase][i][1]);
      if(rel_syst_uncert[UseCase][i][0]>rel_syst_uncert[UseCase][i][1]) CombUncer_low += pow(rel_syst_uncert[UseCase][i][1],2);
      else                                                              CombUncer_low += pow(rel_syst_uncert[UseCase][i][0],2);
    }
    else cout<<"Am I missing something?"<<endl;
  }

  // Systematic uncertainties
  CombUncer_high = (CombUncer_high>0? sqrt(CombUncer_high) : 0);
  cout<<"Rel. Comb. (+) syst. uncert : "<<CombUncer_high<<endl;
  h->SetBinContent(nsyst*2-1+1, CombUncer_high);
  // Adding statistic uncertainties
/*
  cout<<"Rel. stat. uncert. : "<<yerr_high[0]/y[0]<<endl;
  float TotUncer_high = sqrt(pow(CombUncer_high,2)+pow(yerr_high[0]/y[0],2));
  cout<<"Rel. Comb. (+) (stat.+syst.) uncert : "<<TotUncer_high<<endl;
*/

  CombUncer_low = (CombUncer_low >0? sqrt(CombUncer_low)  : 0);
  cout<<"Rel. Comb. (-) syst. uncert : "<<CombUncer_low<<endl; //"+/-"<<endl;
  h1->SetBinContent(nsyst*2-1+1,-CombUncer_low);
  // Adding statistic uncertainties
/*
  cout<<"Rel. stat. uncert. : "<<yerr_high[0]/y[0]<<endl;
  float TotUncer_low = sqrt(pow(CombUncer_low,2)+pow(yerr_low[0]/y[0],2));
  cout<<"Rel. Comb. (-) (stat.+syst.) uncert : "<<TotUncer_low<<endl;
*/
  cout<<"R_{V} = "<<RV[UseCase]<<"^{+"<<CombUncer_high*RV[UseCase]<<"}_{-"<<CombUncer_low*RV[UseCase]<<"}"<<endl;
  
  h->GetXaxis()->SetBinLabel(nsyst*2-1+1,binlabel[nsyst*2-1].c_str());
  h->GetXaxis()->LabelsOption("d");
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitle("(R_{V}-R_{V}^{Nominal})/R_{V}^{Nominal}");
  h->GetYaxis()->SetNdivisions(505);
  //h->GetYaxis()->SetRangeUser(0.0004,0.00119);

  gStyle->SetHistMinimumZero();

  TCanvas *c = new TCanvas("c","",600,600);
  c->cd();
  h->SetFillColor(kBlue);
  h->GetYaxis()->SetRangeUser(-0.5,0.5);
  h->Draw("Bar");

  h1->SetFillColor(kRed);
  h1->GetYaxis()->SetRangeUser(h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  h1->Draw("Barsame");

  TLatex *tex = new TLatex(0.14,0.995,"CMS Preliminary, 5.0 fb^{-1} at #sqrt{s} = 7 TeV");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();

  c->Print((std::string()+"RV_TrendPlotSummary_"+label+"_"+teffname[UseCase]+".eps").c_str());
  c->Write((std::string()+"RV_TrendPlotSummary_"+label+"_"+teffname[UseCase]).c_str());
//  printf("Closing file\n");
  outputfile->Close();

}
