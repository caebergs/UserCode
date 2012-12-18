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


void makeSystUncertTrendPlot_JES_VJets(const int UseCase, const int channel){

  gStyle->SetOptFit(0);

  const int NbOfUseCase = 4;
  // Use case : 
  // 0 : low DM
  // 1 : Int DM 1
  // 2 : Int DM 2
  // 3 : High DM
//  const int UseCase = 0;

  //const int NbOfFiles = 3;
  //TFile *inputfiles[NbOfFiles][2][2];
  TFile *inputfiles[2][2];
  std::string nompath = "$HOME/AnalysisCode/GentStopAnalysis/UGentCode/Leg3/v1/Samples_21102012/";

  inputfiles[0][0] = TFile::Open((nompath+"Wjets_mu.root").c_str());
  inputfiles[0][1] = TFile::Open((nompath+"Zjets_mu.root").c_str());
  inputfiles[1][0] = TFile::Open((nompath+"Wjets_el.root").c_str());
  inputfiles[1][1] = TFile::Open((nompath+"Zjets_el.root").c_str());

  const int bin = 20;
  const int nsyst = 2;
  double x[nsyst+1][bin]={{0}}, xerr[nsyst+1][bin]={{0}}, y[nsyst+1][bin]={{0}}, yerr_low[nsyst+1][bin]={{0}}, yerr_high[nsyst+1][bin]={{0}};
  string binlabel[nsyst] = {"JES (-)","JES (+)"};//,"Total syst."};
  std::string label = "JES" ;
  char name[100];
  const double NNcut[NbOfUseCase] = {0.9,0.95,0.95,0.95};
  
  string teffname[NbOfUseCase] = {"Eff_LowDM_bin", "Eff_IntDM1_bin", "Eff_IntDM2_bin", "Eff_HighDM_bin"};
  string effname = "";
  cout << "Use case : "<<teffname[UseCase]<<endl;

  TEfficiency *tEff1, *tEff2, *tEff3, *tEff4;
  TEfficiency *tDummy;
  TList *pList = new TList();
  TGraphAsymmErrors *tg = 0;

  for(int i=0;i<nsyst+1;i++){

    effname = teffname[UseCase];

    if(i==1) effname += "_m";
    else if(i==2) effname += "_p";

    // W+jets
    tEff1 = (TEfficiency*)inputfiles[0][0]->Get(effname.c_str());
    //tEff2 = ( inputfiles[i][1] ? (TEfficiency*)inputfiles[i][1]->Get(effname.c_str()) : 0 );
    tEff2 = (TEfficiency*)inputfiles[1][0]->Get(effname.c_str());
    
    if(tEff1==0 || tEff2==0) continue;
    tEff1->SetStatisticOption(TEfficiency::kBBayesian);
    tEff2->SetStatisticOption(TEfficiency::kBBayesian);

    //tEff1->Add(*tEff2);
    
    // Z+jets
    tEff3 = (TEfficiency*)inputfiles[0][1]->Get(effname.c_str());
   //tEff2 = ( inputfiles[i][1] ? (TEfficiency*)inputfiles[i][1]->Get(effname.c_str()) : 0 );
    tEff4 = (TEfficiency*)inputfiles[1][1]->Get(effname.c_str());
    
    if(tEff3==0 || tEff4==0) continue;
    tEff3->SetStatisticOption(TEfficiency::kBBayesian);
    tEff4->SetStatisticOption(TEfficiency::kBBayesian);

    std::vector<Double_t>  w_vjets;

    if ((channel&1)!=0) {
      pList->Add((TEfficiency*)tEff1);
      w_vjets.push_back(18934./tEff1->GetTotalHistogram()->GetEntries());
    }
    if ((channel&2)!=0) {
      pList->Add((TEfficiency*)tEff2);
      w_vjets.push_back(14970./tEff2->GetTotalHistogram()->GetEntries());
    }
    //tEff3->Add(*tEff4);
    if ((channel&1)!=0) {
      pList->Add((TEfficiency*)tEff3);
      w_vjets.push_back(1400./tEff3->GetTotalHistogram()->GetEntries());
    }
    if ((channel&2)!=0) {
      pList->Add((TEfficiency*)tEff4);
      w_vjets.push_back(2223./tEff4->GetTotalHistogram()->GetEntries());
    }
    
    // weights
    
    tg = (TGraphAsymmErrors*)tDummy->Combine(pList,"",w_vjets.size(),&w_vjets[0]);//,"",2);
    for(int j=0;j<bin;j++){
      tg->GetPoint(j,x[i][j],y[i][j]);
      //x[i][j] -= tg->GetErrorXlow(j);
      xerr[i][j] = tg->GetErrorXlow(j);      
      cout<<"NN cut : "<<x[i][j]-tg->GetErrorXlow(j)<<" / Eff = "<<y[i][j]<<endl;
      yerr_low[i][j]  = tg->GetErrorYlow(j);
      yerr_high[i][j] = tg->GetErrorYhigh(j);
    }
    tEff1->Delete();
    tEff2->Delete();
    tEff3->Delete();
    tEff4->Delete();
    tg->Delete();
    pList->Clear();
  }

  /**************************************************************/
  // systematic uncertainty trend
  /**************************************************************/
/*
  TH1F *h  = new TH1F("h","",bin,0,bin);
  TH1F *h1 = new TH1F("h1","",bin,0,bin);
  TAxis *Xaxis = h->GetXaxis();
  TAxis *Yaxis = h->GetYaxis();
*/
  TFile *outputfile = new TFile((std::string()+"RV_trendPlot_"+label+"_"+teffname[UseCase]+"_Output.root").c_str(),"RECREATE");
  TH1F *h[nsyst];
  TGraphErrors *tgE[nsyst];
  TMultiGraph *mg = new TMultiGraph;
  mg->SetTitle("title;NN;(R_{V}-R_{V}^{Nominal})/R_{V}^{Nominal}");

  double BinContent[bin] = {0};
  double BinError[bin]   = {0};
  //double CombUncer_low  = 0, CombUncer_error_low   = 0;
  //double CombUncer_high = 0, CombUncer_error_high  = 0;
  for (int i=1;i<nsyst+1;i++){
    h[i-1] = new TH1F(binlabel[i-1].c_str(),"",20,0,1);
    for(int j=0;j<bin;j++){
      BinContent[j] = (y[i][j]-y[0][j])/y[0][j];
      h[i-1]->SetBinContent(j+1,BinContent[j]);
      cout<<"x["<<i<<"]["<<j<<"] = "<<x[i][j]<<" / BinContent["<<j<<"] = "<<BinContent[j]<<endl;
//      if(j>10) cout<<BinContent[j]<<endl;
      BinError[j] = (y[i][j]/y[0][j])*sqrt(pow(yerr_high[i][j]/y[i][j],2)+pow(yerr_high[0][j]/y[0][j],2));
    }
    tgE[i-1] = new TGraphErrors(bin,x[i],BinContent,xerr[i],BinError);
    tgE[i-1]->SetLineColor(i);
    tgE[i-1]->SetLineWidth(2);
  }
//  gStyle->SetHistMinimumZero();
/*
  //h[0]->Smooth(2);
  h[0]->SetFillColor(kRed);
  h[1]->SetFillColor(kRed);
  h[0]->Draw();
  h[1]->Draw("same");
*/

  TF1 *fit1 = new TF1("fit1","pol1",0.2,0.8);
  fit1->SetLineColor(kRed);
  fit1->SetLineWidth(3);
  TF1 *fit2 = new TF1("fit2","pol1",0.2,0.8);
  fit2->SetLineColor(kBlack);
  fit2->SetLineWidth(3);

  tgE[0]->SetFillColor(kRed);
  tgE[0]->SetFillStyle(3004);
  tgE[0]->SetMarkerStyle(23);
  tgE[0]->Fit(fit1,"EMR");//,"",0.2,0.8);

  tgE[1]->SetFillColor(kBlue);
  tgE[1]->SetFillStyle(3005);
  tgE[1]->SetMarkerStyle(22);
  tgE[1]->Fit(fit2,"EMR");//,"",0.2,0.8);

  TCanvas *c = new TCanvas("c","",600,600);
  c->cd();
  mg->Add(tgE[0]);
  mg->Add(tgE[1]);
  mg->Draw("AP2");
//  mg->Draw("P0");
  
  cout<<"Extrapolated uncert. ("<<binlabel[0]<<") at "<<NNcut[UseCase]<<" = "<<fit1->Eval(NNcut[UseCase])<<endl;
  cout<<"Extrapolated uncert. ("<<binlabel[1]<<") at "<<NNcut[UseCase]<<" = "<<fit2->Eval(NNcut[UseCase])<<endl;
/*
//    Xaxis->SetBinLabel(j,binlabel[j-1].c_str());

    if(BinContent>0){
      h->SetBinContent(j,BinContent);
      h->SetBinError(j,BinError);
      cout<<"Bin+ "<<j<<", Content : "<<BinContent<<" / Error : "<<BinError<<endl;//" / Syst² : "<<pow(BinContent,2)-pow(BinError,2)<<endl;
      CombUncer_high += pow(BinContent,2);//-pow(BinError,2);
      CombUncer_error_high += pow(BinError,2);
    }
    else{
      h1->SetBinContent(j,BinContent);
      h1->SetBinError(j,BinError);
      cout<<"Bin- "<<j<<", Content : "<<BinContent<<" / Error : "<<BinError<<endl;//" / Syst² : "<<pow(BinContent,2)-pow(BinError,2)<<endl;
      CombUncer_low       += pow(BinContent,2);//-pow(BinError,2);
      CombUncer_error_low += pow(BinError,2);
    }
  }

  // Systematic uncertainties
  CombUncer_high = (CombUncer_high>0? sqrt(CombUncer_high) : 0);
  CombUncer_error_high = sqrt(CombUncer_error_high); 
  cout<<"Rel. Comb. (+) syst. uncert : "<<CombUncer_high<<"+/-"<<CombUncer_error_high<<endl;
  h->SetBinContent(nsyst, CombUncer_high);
  // Adding statistic uncertainties
  cout<<"Rel. stat. uncert. : "<<yerr_high[0]/y[0]<<endl;
  float TotUncer_high = sqrt(pow(CombUncer_high,2)+pow(yerr_high[0]/y[0],2));
  cout<<"Rel. Comb. (+) (stat.+syst.) uncert : "<<TotUncer_high<<endl;

  CombUncer_low = (CombUncer_low >0? sqrt(CombUncer_low)  : 0);
  CombUncer_error_low = sqrt(CombUncer_error_low);
  cout<<"Rel. Comb. (-) syst. uncert : "<<CombUncer_low<<"+/-"<<CombUncer_error_low<<endl;
  h1->SetBinContent(nsyst,-CombUncer_low);
  // Adding statistic uncertainties
  cout<<"Rel. stat. uncert. : "<<yerr_high[0]/y[0]<<endl;
  float TotUncer_low = sqrt(pow(CombUncer_low,2)+pow(yerr_low[0]/y[0],2));
  cout<<"Rel. Comb. (-) (stat.+syst.) uncert : "<<TotUncer_low<<endl;

  cout<<"*********** Results *****************"<<endl;
  cout<<std::fixed<<setprecision(6);
  cout<<"$R_{V}$ & $"<<y[0]<<"\\pm"<<yerr_low[0]<<"^{+"<<y[0]*CombUncer_high<<"}_{-"<<y[0]*CombUncer_low<<"}$"<<endl;
  cout<<"$R_{V}$ & $"<<y[0]<<"\\pm"<<yerr_high[0]<<"^{+"<<y[0]*CombUncer_high<<"}_{-"<<y[0]*CombUncer_low<<"}$"<<endl;
  cout<<"*************************************"<<endl;
  
  h->GetXaxis()->SetBinLabel(nsyst,binlabel[nsyst-1].c_str());
  h->GetXaxis()->LabelsOption("d");
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitle("(R_{V}-R_{V}^{Nominal})/R_{V}^{Nominal}");
  h->GetYaxis()->SetNdivisions(505);
  //h->GetYaxis()->SetRangeUser(0.0004,0.00119);

  gStyle->SetHistMinimumZero();

  h->SetFillColor(kBlue);
  h->GetYaxis()->SetRangeUser(-0.5,0.5);
  h->Draw("E1Bar");

  h1->SetFillColor(kRed);
  h1->GetYaxis()->SetRangeUser(h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  h1->Draw("E1Barsame");
*/
  TLegend *leg = new TLegend(0.17,0.75,0.65,0.92,NULL,"brNDC");
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->AddEntry(tgE[0],binlabel[0].c_str(),"fp");
  leg->AddEntry(tgE[1],binlabel[1].c_str(),"fp");

  leg->Draw("same");

  TLatex *tex = new TLatex(0.14,0.995,"CMS Preliminary, 5.0 fb^{-1} at #sqrt{s} = 7 TeV");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();



  c->Print((std::string()+"RV_TrendPlot_"+label+"_"+teffname[UseCase]+".eps").c_str());
  c->Write((std::string()+"RV_TrendPlot_"+label+"_"+teffname[UseCase]).c_str());

  outputfile->Close();
  c->Close();
  tgE[0]->Delete();
  tgE[1]->Delete();
  fit1->Delete();
  fit2->Delete();
  

}
