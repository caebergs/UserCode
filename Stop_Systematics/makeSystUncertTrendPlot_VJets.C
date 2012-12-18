#include <TROOT.h>
#include <TAxis.h>
#include <TCanvas.h>
//#include <TColor.h>
#include <TEfficiency.h>
#include <TF1.h>
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
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


void makeSystUncertTrendPlot_VJets(const Int_t UseCase, UInt_t uncertMode, const UInt_t channel){

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  const int NbOfUseCase = 4;
  // Use case : 
  // 0 : low DM
  // 1 : Int DM 1
  // 2 : Int DM 2
  // 3 : High DM
  
  //const int NbOfFiles = 3;
  //TFile *inputfiles[NbOfFiles][2][2];
  const int bin = 20;
  const int nsyst = 2; //Number of variants for a given systematics
  const int nVjetsSamples = 2;
  const int nNominal = (uncertMode==2 ? 2 : 1);
  const int nSyNo = nsyst + nNominal +(uncertMode==2?2:0) ; //For the combined (nominal and syst)
  
//  UInt_t channel=3 ;
  std::string path = "/user/caebergs/VJetEstimation/testChain/samplesSyst/";
  std::string nompath = "/user/caebergs/Leg3/Leg3Output_5Dreweigh/";
  TFile ***nominalinputfiles = new TFile**[2];
  for (int i=0; i<2; i++) {
    nominalinputfiles[i] = new TFile*[nVjetsSamples];
    for (int j=0; j<nVjetsSamples; j++) {
      nominalinputfiles[i][j] = NULL;
    }
  }
  TFile ***inputfiles = new TFile**[nsyst*2] ;
  for (int i=0; i<nsyst*2; i++) {
    inputfiles[i] = new TFile*[nVjetsSamples];
    for (int j=0; j<nVjetsSamples; j++) {
      inputfiles[i][j] = NULL;
    }
  }
  
//  std::string path = "$HOME/AnalysisCode/GentStopAnalysis/UGentCode/Leg3/v1/Samples_21102012/";

  nominalinputfiles[0][0] = TFile::Open((nompath+"Wjets_mu.root").c_str());
  nominalinputfiles[0][1] = TFile::Open((nompath+"Zjets_mu.root").c_str());
  nominalinputfiles[1][0] = TFile::Open((nompath+"Wjets_el.root").c_str());
  nominalinputfiles[1][1] = TFile::Open((nompath+"Zjets_el.root").c_str());
//  inputfiles[1][0] = TFile::Open((nompath+"Wjets_el.root").c_str());
//  inputfiles[1][1] = TFile::Open((nompath+"Zjets_el.root").c_str());

  inputfiles[0][0] = TFile::Open((uncertMode==0 ? (path+"Wjets-matchingDown_mu.root").c_str() :
                                  (uncertMode==1 ? (path+"Wjets-scaleDown_mu.root").c_str() : 
                                   (uncertMode==2 ? (nompath+"Wjets_mu.root").c_str(): ""))));
  inputfiles[0][1] = TFile::Open((nompath+"Zjets_mu.root").c_str());
  inputfiles[1][0] = TFile::Open((uncertMode==0 ? (path+"Wjets-matchingDown_el.root").c_str() :
                                  (uncertMode==1 ? (path+"Wjets-scaleDown_el.root").c_str() : 
                                   (uncertMode==2 ? (nompath+"Wjets_el.root").c_str():""))));
  inputfiles[1][1] = TFile::Open((nompath+"Zjets_el.root").c_str());

  inputfiles[2][0] = TFile::Open((uncertMode==0 ? (path+"Wjets-matchingUp_mu.root").c_str() :
                                  (uncertMode==1 ? (path+"Wjets-scaleUp_mu.root").c_str() :
                                   (uncertMode==2 ? (nompath+"Wjets_mu.root").c_str():""))));
  inputfiles[2][1] = TFile::Open((nompath+"Zjets_mu.root").c_str());
  inputfiles[3][0] = TFile::Open((uncertMode==0 ? (path+"Wjets-matchingUp_el.root").c_str() :
                                  (uncertMode==1 ? (path+"Wjets-scaleUp_el.root").c_str() : 
                                   (uncertMode==2 ? (nompath+"Wjets_el.root").c_str():""))));
  inputfiles[3][1] = TFile::Open((nompath+"Zjets_el.root").c_str());

    double x[nSyNo][bin], xerr[nSyNo][bin], y[nSyNo][bin], yerr_low[nSyNo][bin], yerr_high[nSyNo][bin];
  for (int i=0; i<nSyNo; i++) {
    for (int j=0; j<bin; j++) {
      x[i][j] = 0. ;
      xerr[i][j] = 0. ;
      y[i][j] = 0. ;
      yerr_low[i][j] = 0. ;
      yerr_high[i][j] = 0. ;
    }
  }
  
  //  string binlabel[nsyst] = {"ME-PS matching (-)","ME-PS matching (+)"};//,"Total syst."};
//  string binlabel[2*nsyst] = { "ME-PS matching (-)", "ME-PS matching (+)", "Scale Q^{2} (-)", "Scale Q^{2} (+)", "5D Rew. (muon)", "5D Rew. (elec.)" } ;
  std::string *binlabel = new std::string[nsyst+(uncertMode==2?1:0)];
  binlabel[0] = (uncertMode==0 ? "ME-PS matching (-)":
                 (uncertMode==1 ? "Scale Q^{2} (-)":
                  (uncertMode==2 ? "5D Rew. (muon)":"")));
  binlabel[1] = (uncertMode==0 ? "ME-PS matching (+)":
                 (uncertMode==1 ? "Scale Q^{2} (+)":
                  (uncertMode==2 ? "5D Rew. (elec.)":"")));//,"Total syst."};
  if (uncertMode==2) {
    binlabel[2] = "5D Rew. (comb.)" ;
  }
  std::string label = (uncertMode==0 ? "Matching":
                                  (uncertMode==1 ? "Scale":
                                   (uncertMode==2 ? "5DRew":"")));
  label+="_channel";
  label+=(channel==1?"1":(channel==2?"2":(channel==3?"3":"X")));
  for (int i=0; i<nsyst; i++) {
    cout << binlabel[i] << endl;
  }
//  char name[100];
  const double NNcut[NbOfUseCase] = {0.9,0.95,0.95,0.95};
  
  string NNsuffix[NbOfUseCase] = {"","_IM1","_IM2","_HM"};
  string teffname[NbOfUseCase] = {"Eff_LowDM_bin", "Eff_IntDM1_bin", "Eff_IntDM2_bin", "Eff_HighDM_bin"};
  string effname = "";
  cout << "Use case : "<<teffname[UseCase]<<endl;

    TEfficiency *tEff1 = NULL, *tEff2 = NULL, *tEff3 = NULL, *tEff4 = NULL;
TList *pList = new TList();
TGraphAsymmErrors *tg = NULL;
  

  std::vector<Double_t> w_v; // For 5D reweighting combination
//  printf("Before looping on files\n");
  for(int i=0;i<nNominal+nsyst;i++){
//    printf(" %d\n", i);
    if ((uncertMode==2)&&(i>0/*Only once for the weights*/)) {
      continue ;
    }
//    printf(" %d passed\n", i);
    effname = teffname[UseCase];
    w_v.clear();
    TFile ***fich = /*(TFile***)*/ (i<nNominal ? nominalinputfiles : inputfiles) ;
//    printf("fich loaded\n");
//    printf(" __> %s\n", fich[0+nVjetsSamples*(i<nNominal?i:i-nNominal)][0]->GetName());
//    if(i==1) effname += "_m";
//    else if(i==2) effname += "_p";
    // W+jets
//    printf(" -->> %d %d \n", 0+nVjetsSamples*(i<nNominal?i:i-nNominal), 0/*, fich[0+2*(i<nNominal?i:i-nNominal)][0]->GetName()*/);
  
    tEff1 = (TEfficiency*)fich[0+nVjetsSamples*(i<nNominal?i:i-nNominal)][0]->Get(effname.c_str());
    //tEff2 = ( inputfiles[i][1] ? (TEfficiency*)inputfiles[i][1]->Get(effname.c_str()) : 0 );
    tEff2 = (TEfficiency*)fich[1+nVjetsSamples*(i<nNominal?i:i-nNominal)][0]->Get(effname.c_str());
    
    //tEff1->Add(*tEff2);
    
    // Z+jets
    tEff3 = (TEfficiency*)fich[0+nVjetsSamples*(i<nNominal?i:i-nNominal)][1]->Get(effname.c_str());
    //tEff2 = ( inputfiles[i][1] ? (TEfficiency*)inputfiles[i][1]->Get(effname.c_str()) : 0 );
    tEff4 = (TEfficiency*)fich[1+nVjetsSamples*(i<nNominal?i:i-nNominal)][1]->Get(effname.c_str());
    
    if(tEff1==NULL || tEff2==NULL) { cout<<"Couldn't find the TEFF 1&2"<<endl; continue; }
    tEff1->SetStatisticOption(TEfficiency::kBBayesian);
    tEff2->SetStatisticOption(TEfficiency::kBBayesian);
    if(tEff3==NULL || tEff4==NULL) { cout<<"Couldn't find the TEFF 1&2"<<endl; continue; }
    tEff3->SetStatisticOption(TEfficiency::kBBayesian);
    tEff4->SetStatisticOption(TEfficiency::kBBayesian);
    
//    printf("Efficiencies loaded\n");
    
    //tEff3->Add(*tEff4);
    
    // weights
    std::vector<Double_t> w_vjets; // For TEfficiency combination
    if ((channel&1)!=0) {
      pList->Add((TEfficiency*)tEff1);
      w_vjets.push_back(18934/tEff1->GetTotalHistogram()->GetEntries());
      w_v.push_back(18934);
    }
    if ((channel&2)!=0) {
      pList->Add((TEfficiency*)tEff2);
      w_vjets.push_back(14970/tEff2->GetTotalHistogram()->GetEntries());
      w_v.push_back(14970);
    }

    //tEff3->Add(*tEff4);
    if ((channel&1)!=0) {
      pList->Add((TEfficiency*)tEff3);
      w_vjets.push_back(1400/tEff3->GetTotalHistogram()->GetEntries());
      w_v.push_back(1400);
    }
    if ((channel&2)!=0) {
      pList->Add((TEfficiency*)tEff4);
      w_vjets.push_back(2223/tEff4->GetTotalHistogram()->GetEntries());
      w_v.push_back(2223);
    }
//    printf("TGraph\n");
    
    tg = (TGraphAsymmErrors*) TEfficiency::Combine(pList,"mode",w_vjets.size(),&w_vjets[0]);//,"",2);
    for(int j=0;j<bin;j++){
      tg->GetPoint(j,x[i][j],y[i][j]);
      //x[i][j] -= tg->GetErrorXlow(j);
      xerr[i][j] = tg->GetErrorXlow(j);      
      cout<<"NN cut : "<<x[i][j]-tg->GetErrorXlow(j)<<" / Eff = "<<y[i][j]<<endl;
      yerr_low[i][j]  = tg->GetErrorYlow(j);
      yerr_high[i][j] = tg->GetErrorYhigh(j);
    }
//    printf("End\n");
    
    tEff1->Delete();
    tEff2->Delete();
    tEff3->Delete();
    tEff4->Delete();
    tg->Delete();
    pList->Clear();
  }
//  printf("At files Mid-loop\n");
  if (uncertMode==2) {
    for( int j=0 ; j<bin ; j++ ) {
    Double_t xPoint = 0. + j *0.05 ;
    Double_t xError = 0.025 ;
    
    float Num_mu_nominal = 0.;
    float Num_el_nominal = 0.;
    float Num_nominal    = 0.;
    float Num_mu         = 0.;
    float Num_el         = 0.;
    float Num            = 0.;
    
    float Den_mu_nominal = 0.;
    float Den_el_nominal = 0.;
    float Den_nominal    = 0.;
    float Den_mu         = 0.;
    float Den_el         = 0.;
    float Den            = 0.;
    
    std::vector<UInt_t> idxVect ;
    idxVect.push_back(0);
    idxVect.push_back(1);
    
    Int_t w_idx = 0;
    for (std::vector<UInt_t>::iterator idx=idxVect.begin() ; idx!=idxVect.end() ; idx++) {
/*
      printf("%d %d %d\n", *idx, w_v.size(), w_idx);
      printf(" %s \n", nominalinputfiles[0][*idx]->GetName());
      printf(" %s \n", nominalinputfiles[1][*idx]->GetName());
      printf(" %s \n", inputfiles[0][*idx]->GetName());
      printf(" %s \n", inputfiles[1][*idx]->GetName());
*/
      TH1D *Rv_mu_nominal = (TH1D*)nominalinputfiles[0+2*0][*idx]->Get(("MVA"+NNsuffix[UseCase]+"_Mu_4j").c_str());
      TH1D *Rv_el_nominal = (TH1D*)nominalinputfiles[1+2*0][*idx]->Get(("MVA"+NNsuffix[UseCase]+"_El_4j").c_str());
      TH1D *Rv_mu = (TH1D*)inputfiles[0+2*0][*idx]->Get(("MVA"+NNsuffix[UseCase]+"_Mu_4j_ShapeCorrected").c_str());
      TH1D *Rv_el = (TH1D*)inputfiles[1+2*0][*idx]->Get(("MVA"+NNsuffix[UseCase]+"_El_4j_ShapeCorrected").c_str());
//      printf(" %lu %lu %lu %lu \n", (unsigned long) Rv_mu_nominal, (unsigned long) Rv_el_nominal, (unsigned long) Rv_mu, (unsigned long) Rv_el);
      //  int CutBin = (UseCase==0 ? 233 : 237);
      
      if ((channel&1)!=0) {
        Num_mu_nominal += Rv_mu_nominal->Integral(/*CutBin*/Rv_mu_nominal->FindBin(xPoint),Rv_mu_nominal->GetNbinsX()+1) *(w_v[w_idx]/Rv_mu_nominal->Integral(0,-1));
        Num_mu         += Rv_mu->Integral(/*CutBin*/Rv_mu->FindBin(xPoint),Rv_mu->GetNbinsX()+1) *(w_v[w_idx]/Rv_mu->Integral(0,-1));
        Den_mu_nominal += Rv_mu_nominal->Integral(0,Rv_mu_nominal->GetNbinsX()+1) *(w_v[w_idx]/Rv_mu_nominal->Integral(0,-1));
        Den_mu         += Rv_mu->Integral(0,Rv_mu->GetNbinsX()+1) *(w_v[w_idx]/Rv_mu->Integral(0,-1));
        w_idx++;
      }
      if ((channel&2)!=0) {
        Num_el_nominal += Rv_el_nominal->Integral(/*CutBin*/Rv_el_nominal->FindBin(xPoint),Rv_el_nominal->GetNbinsX()+1) *(w_v[w_idx]/Rv_el_nominal->Integral(0,-1));
        Num_el         += Rv_el->Integral(/*CutBin*/Rv_el->FindBin(xPoint),Rv_el->GetNbinsX()+1) *(w_v[w_idx]/Rv_el->Integral(0,-1));
        Den_el_nominal += Rv_el_nominal->Integral(0,Rv_el_nominal->GetNbinsX()+1) *(w_v[w_idx]/Rv_el_nominal->Integral(0,-1));
        Den_el         += Rv_el->Integral(0,Rv_el->GetNbinsX()+1) *(w_v[w_idx]/Rv_el->Integral(0,-1));
        w_idx++;
      }
      Num_nominal    = Num_mu_nominal+Num_el_nominal;
      Num            = Num_mu+Num_el;
      
      Den_nominal    = Den_mu_nominal+Den_el_nominal;
      Den            = Den_mu+Den_el;
    }
    
    TList* tList = new TList();
    TH1D passe("passe", "passe", 6, 0., 6.);
    TH1D total("total", "total", 6, 0., 6.);

    passe.SetBinContent(1, Num_nominal);
    total.SetBinContent(1, Den_nominal);
    passe.SetBinContent(2, Num);
    total.SetBinContent(2, Den);
    passe.SetBinContent(3, Num_mu_nominal);
    total.SetBinContent(3, Den_mu_nominal);
    passe.SetBinContent(4, Num_mu);
    total.SetBinContent(4, Den_mu);
    passe.SetBinContent(5, Num_el_nominal);
    total.SetBinContent(5, Den_el_nominal);
    passe.SetBinContent(6, Num_el);
    total.SetBinContent(6, Den_el);
    TEfficiency eff(passe, total);
    tList->Add(&eff);
    
    TGraphAsymmErrors *tgr = (TGraphAsymmErrors*) TEfficiency::Combine(tList,"",1,NULL);//,"",2);
    Double_t dummy = 0.;
    Double_t nominal = -1.;
    tgr->GetPoint(2,dummy,y[0][j]);
    tgr->GetPoint(4,dummy,y[1][j]);
    tgr->GetPoint(3,dummy,y[2][j]);
    tgr->GetPoint(5,dummy,y[3][j]);
    
    tgr->GetPoint(0,dummy,y[4][j]);
    tgr->GetPoint(1,dummy,y[5][j]);

    x[0][j] = x[1][j] = x[2][j] = x[3][j]     = x[4][j] = x[5][j] = xPoint ;
    xerr[0][j] = xerr[1][j] = xerr[2][j] = xerr[3][j]   = xerr[4][j] = xerr[5][j] = xError ;
//    y[3][j] = y[0][j] * y[3][j] / nominal ;

    //x[i][j] -= tg->GetErrorXlow(j);
    //xerr[3][j] = tg->GetErrorXlow(j);      
    cout<<"NN cut : "<<x[3][j]-xerr[3][j]<<" / Eff = "<<y[3][j]<<endl;
    yerr_low[0][j]  = tgr->GetErrorYlow(2);
    yerr_high[0][j] = tgr->GetErrorYhigh(2);
    yerr_low[1][j]  = tgr->GetErrorYlow(4);
    yerr_high[1][j] = tgr->GetErrorYhigh(4);
    yerr_low[2][j]  = tgr->GetErrorYlow(3);
    yerr_high[2][j] = tgr->GetErrorYhigh(3);
    yerr_low[3][j]  = tgr->GetErrorYlow(5);
    yerr_high[3][j] = tgr->GetErrorYhigh(5);

    yerr_low[4][j]  = tgr->GetErrorYlow(0);
    yerr_high[4][j] = tgr->GetErrorYhigh(0);
    yerr_low[5][j]  = tgr->GetErrorYlow(1);
    yerr_high[5][j] = tgr->GetErrorYhigh(1);

//    int k=0;
//    printf("%f %f\n", yerr_high[][j], yerr_high[][j]);
    //    yerr_low[3][j]  = tgr->GetErrorYlow(1)   * y[3][j] / nominal ;
//    yerr_high[3][j] = tgr->GetErrorYhigh(1)  * y[3][j] / nominal ;
    tList->Delete();
    tgr->Delete();
    
    }
  }
  for (int i=0; i<nNominal+2*nsyst; i++) {
    TFile ***fich = /*(TFile***)*/ (i<nNominal ? nominalinputfiles : inputfiles) ;
    for (int j=0; j<nVjetsSamples; j++) {
      fich[(i<nNominal ? i : i-nNominal)][j]->Close();
    }
  }
  
  printf("End of loop on files\n");
    

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
  printf("Output file opened\n");
  
  TH1F *h[nsyst];
  TGraphErrors *tgE[nsyst];
  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("title;NN;(R_{V}-R_{V}^{Nominal})/R_{V}^{Nominal}");
  printf("Fit result array\n");
//  TList *FitResults = new TList();// nsyst ;
  TFitResultPtr **FitResults = new TFitResultPtr*[nsyst] ;

  double BinContent[bin];
  double BinError[bin];
  //double CombUncer_low  = 0, CombUncer_error_low   = 0;
  //double CombUncer_high = 0, CombUncer_error_high  = 0;
  for (int i=nNominal;i<nSyNo;i++){
    printf("Loop on nsyst : %d\n", i);
    if (uncertMode==2 && i==nSyNo-1) {
      continue;
    }
    h[i-nNominal] = new TH1F(binlabel[i-nNominal].c_str(),"",20,0,1);
    h[i-nNominal]->SetDirectory(NULL);
    for(int j=0;j<bin;j++){
      if (uncertMode==2) {
        if (i==nSyNo-2) {
          BinContent[j] = (y[i+1][j]-y[i][j])/y[i][j];
          BinError[j] = (y[i+1][j]/y[i][j])*sqrt(pow(yerr_high[i+1][j]/y[i+1][j],2)+pow(yerr_high[i][j]/y[i][j],2));
        } else {
          BinContent[j] = (y[i][j]-y[i-nNominal][j])/y[i-nNominal][j];
          BinError[j] = (y[i][j]/y[i-nNominal][j])*sqrt(pow(yerr_high[i][j]/y[i][j],2)+pow(yerr_high[i-nNominal][j]/y[i-nNominal][j],2));
        }
      } else {
        BinContent[j] = (y[i][j]-y[0][j])/y[0][j];
        BinError[j] = (y[i][j]/y[0][j])*sqrt(pow(yerr_high[i][j]/y[i][j],2)+pow(yerr_high[0][j]/y[0][j],2));
      }
      h[i-nNominal]->SetBinContent(j+1,BinContent[j]);
      cout<<"x["<<i<<"]["<<j<<"] = "<<x[i][j]<<" / BinContent["<<j<<"] = "<<BinContent[j]<<endl;
//      if(j>10) cout<<BinContent[j]<<endl;
    }
    tgE[i-nNominal] = new TGraphErrors(bin,x[i],BinContent,xerr[i],BinError);
//    tgE[i-nNominal]->SetDirectory(NULL);
    tgE[i-nNominal]->SetLineColor(i-nNominal);
    tgE[i-nNominal]->SetLineWidth(2);
  }
//  printf("Bin content and error set\n");
  
//  gStyle->SetHistMinimumZero();
/*
  //h[0]->Smooth(2);
  h[0]->SetFillColor(kRed);
  h[1]->SetFillColor(kRed);
  h[0]->Draw();
  h[1]->Draw("same");
*/
//  printf("Defining fits\n");
  
  TF1 **fits = new TF1*[nsyst+(uncertMode==2?1:0)];
//  TList* fits = new TList(); //nsyst
  
  for (int i=0; i<nsyst+(uncertMode==2?1:0); i++) {
    char name[50];
    char function[50];
    sprintf(name, "fit%d", 1+i);
    sprintf(function, "pol1");
    TF1* tmp = new TF1(name,function,0.2,0.8);
    tmp->SetLineColor((i==0 ? kRed : kBlack));
    tmp->SetLineWidth(3);
    fits[i] = tmp;
//    fits->AddAt(tmp, i);
}
//  printf("Fitting\n");
  for (int i=0; i<nsyst+(uncertMode==2?1:0); i++) {
    tgE[i]->SetFillColor((i==0 ? kRed : kBlue));
    tgE[i]->SetFillStyle(3004+i);
    tgE[i]->SetMarkerStyle(22);
//TFitResultPtr ptr = 
    FitResults[i] = new TFitResultPtr (tgE[i]->Fit( fits[i],"SEMR"));
//    FitResults[i] = new TFitResultPtr (tgE[i]->Fit((TF1*) fits->At(i),"SEMR"));
//    FitResults->AddAt(new TFitResultPtr (tgE[i]->Fit((TF1*) fits->At(i),"SEMR")), i) ;
    //    new ( FitResults[i]) TFitResultPtr(tgE[i]->Fit(fits[i],"SEMR"));
//    .Copy(    *(FitResults[i]));//,"",0.2,0.8);
  }
  
  TCanvas *c = new TCanvas("c","",600,600);
  c->cd();
//  printf("Drawing\n");
  for (int i=0; i<nsyst; i++) {
    mg->Add(tgE[i]);
  }
  mg->SetMaximum(1.);
  mg->SetMinimum(-1.);
    mg->Draw("AP2");
  //mg->Draw("P0");
  
  for (UInt_t i=0; i<nsyst+(uncertMode==2?1:0); i++) {
    Double_t fittedVal = (fits[i])->Eval(NNcut[UseCase]) ;
    if ((uncertMode==2) && (i==3) && false) {
      fittedVal = sqrt( pow((fits[0])->Eval(NNcut[UseCase]), 2.) + pow((fits[1])->Eval(NNcut[UseCase]), 2.) ) ;
    }
    if ((uncertMode==2) && ((i+1)!=channel)) {
      cout << "--> extr unc ("<<binlabel[i]<<") at "<<NNcut[UseCase]<<" = "<< fittedVal <<endl;
      continue;
    }
    cout<<"Extrapolated uncert. ("<<binlabel[i]<<") at "<<NNcut[UseCase]<<" = "<< fittedVal <<endl;
//    cout<<"Extrapolated uncert. ("<<binlabel[i]<<") at "<<NNcut[UseCase]<<" = "<<((TF1*)fits->At(i))->Eval(NNcut[UseCase])<<endl;
  }
  
  /*
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  assert(fitter != 0);
  fitter->PrintResults(0,0);
*/
  //double * covMatrix = fitter->GetCovarianceMatrixElement(0,1);

  double covMatrixElem_fits[nsyst];
  double  y_err[nsyst];
  for (int i=0; i<nsyst+(uncertMode==2?1:0); i++) {
//    covMatrixElem_fits[i] = (*((TFitResultPtr*) FitResults->At(i)))->GetCovarianceMatrix()(0,1);
    covMatrixElem_fits[i] = (*((TFitResultPtr*) FitResults[i]))->GetCovarianceMatrix()(0,1);
    cout<<"Fit "<<i<<" : Element 0/1 : "<<covMatrixElem_fits[i]<<endl;
    y_err[i] = sqrt(pow(NNcut[UseCase]*fits[i]->GetParError(1),2)+pow(fits[i]->GetParError(0),2)+2*NNcut[UseCase]*covMatrixElem_fits[i]);
//    y_err[i] = sqrt(pow(NNcut[UseCase]*((TF1*)fits->At(i))->GetParError(1),2)+pow(((TF1*)fits->At(i))->GetParError(0),2)+2*NNcut[UseCase]*covMatrixElem_fits[i]);
}

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
//  printf("Legend\n");
    TLegend *leg = new TLegend(0.17,0.80,0.65,0.92,NULL,"brNDC");
  leg->SetTextFont(62);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);

  for (int i=0; i<nsyst; i++) {
    leg->AddEntry(tgE[i],binlabel[i].c_str(),"f");
  }
  
  leg->Draw("same");

  TLatex *tex = new TLatex(0.14,0.995,"CMS Preliminary, 5.0 fb^{-1} at #sqrt{s} = 7 TeV");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();

//  c->Update();
  //  c->SetGridy();
  c->Print((std::string()+"RV_TrendPlot_"+label+"_"+teffname[UseCase]+".eps").c_str());
  c->Write((std::string()+"RV_TrendPlot_"+label+"_"+teffname[UseCase]).c_str());
//  printf("Closing file\n");
  outputfile->Close();
//  c->SetDirectory(NULL);
//  c->Close();
//  c->Delete();
  
//  printf("Freeing the memory\n");
  for (int i=0; i<nsyst; i++) {
    printf("_\n");
    delete FitResults[i];
    printf("a\n");
    fits[i]->Delete() ;
    printf("b\n");
    tgE[i]->Delete() ;
    printf("c\n");
//    h[i]->Delete() ;
//    printf("d\n");

  }

  printf("1\n");
  delete[] fits ;
//  printf("2\n");
//  delete[] tgE ;
//  printf("3\n");
  delete[] FitResults ;
  printf("4\n");

  delete[] binlabel ;
}

int main(const int argc, const char* argv[]) 
{
  makeSystUncertTrendPlot_VJets(0, 1, 3);
}

    
