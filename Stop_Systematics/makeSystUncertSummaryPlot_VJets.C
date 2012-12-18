#include <TROOT.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

#include <math.h>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;


void makeSystUncertSummaryPlot_VJets(UInt_t UseCase, UInt_t channel){
  printf("Starting the program\n");

  gStyle->SetOptFit(0);
  printf("Style options set\n");

  const int NbOfUseCase = 4;
  // Use case : 
  // 0 : low DM
  // 1 : Int DM 1
  // 2 : Int DM 2
  // 3 : High DM
  //const int UseCase = 0;

  //const int NbOfFiles = 3;
  //TFile *inputfiles[NbOfFiles][2][2];
  TFile *inputfiles[2][2];

  std::string nompath = "$HOME/AnalysisCode/GentStopAnalysis/UGentCode/Leg3/v1/Samples_15112012/";
  inputfiles[0][0] = TFile::Open((nompath+"Wjets_mu.root").c_str());
  inputfiles[0][1] = TFile::Open((nompath+"Wjets_el.root").c_str());
  inputfiles[1][0] = TFile::Open((nompath+"Zjets_mu.root").c_str());
  inputfiles[1][1] = TFile::Open((nompath+"Zjets_el.root").c_str());

  int bin = 0;
  const int npoints = 4;
  double x[npoints]={0}, y[npoints]={0}, yerr_low[npoints]={0}, yerr_high[npoints]={0};
  string binlabel[npoints] = {"JES (-)","JES (+)", "5D rew.","Total syst."};
  char name[100];
  
  std::string comment = "";
  std::string chanName[4] = { "", "SemiMuon", "SemiElectron", "Combined" };
  string NNsuffix[NbOfUseCase] = {"","_IM1","_IM2","_HM"};
  string teffname[NbOfUseCase] = {"Eff_LowDM", "Eff_IntDM1", "Eff_IntDM2", "Eff_HighDM"};
  string effname = "";
  cout << "Use case : "<<teffname[UseCase]<<endl;

  TEfficiency *tEff1, *tEff2, *tEff3, *tEff4;
  TList *pList = new TList();
  TGraphAsymmErrors *tg = 0;
  printf("vector\n");
  std::vector<Double_t> w_v; // For 5D reweighting combination
  printf("Before looping on the files") ;

  for(int i=0;i<npoints-1;i++){
    w_v.clear();
    
    effname = teffname[UseCase];

    if(i==1) effname += "_m";
    else if(i==2) effname += "_p";

    // weights
    std::vector<Double_t> w_vjets; // For TEfficiency combination
    
    // W+jets
    tEff1 = (TEfficiency*)inputfiles[0][0]->Get(effname.c_str());
    tEff1->SetStatisticOption(TEfficiency::kBBayesian);
    //tEff2 = ( inputfiles[i][1] ? (TEfficiency*)inputfiles[i][1]->Get(effname.c_str()) : 0 );
    tEff2 = (TEfficiency*)inputfiles[0][1]->Get(effname.c_str());
    tEff2->SetStatisticOption(TEfficiency::kBBayesian);

    if(tEff1==0 || tEff2==0) continue;

    
    // Z+jets
    tEff3 = (TEfficiency*)inputfiles[1][0]->Get(effname.c_str());
    tEff3->SetStatisticOption(TEfficiency::kBBayesian);
    //tEff2 = ( inputfiles[i][1] ? (TEfficiency*)inputfiles[i][1]->Get(effname.c_str()) : 0 );
    tEff4 = (TEfficiency*)inputfiles[1][1]->Get(effname.c_str());
    tEff4->SetStatisticOption(TEfficiency::kBBayesian);

    if(tEff3==0 || tEff4==0) continue;

    //tEff1->Add(*tEff2);
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
    
    
    tg = (TGraphAsymmErrors*)TEfficiency::Combine(pList,"mode",w_vjets.size(),&w_vjets[0]);//,"",2);
/*    
    tEff1->Add(*tEff3);
    y[i]         = tEff1->GetEfficiency(bin);
    yerr_low[i]  = tEff1->GetEfficiencyErrorLow(bin);
    yerr_high[i] = tEff1->GetEfficiencyErrorUp(bin);
*/
    tg->GetPoint(bin,x[i],y[i]);
    yerr_low[i]  = tg->GetErrorYlow(bin);
    yerr_high[i] = tg->GetErrorYhigh(bin);

    tEff1->Delete();
    tEff2->Delete();
    tEff3->Delete();
    tEff4->Delete();
    tg->Delete();
    pList->Clear();
  }



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
    
    TH1D *Rv_mu_nominal = (TH1D*)inputfiles[*idx][0]->Get(("MVA"+NNsuffix[UseCase]+"_Mu_4j").c_str());
    TH1D *Rv_mu = (TH1D*)inputfiles[*idx][0]->Get(("MVA"+NNsuffix[UseCase]+"_Mu_4j_ShapeCorrected").c_str());
    TH1D *Rv_el_nominal = (TH1D*)inputfiles[*idx][1]->Get(("MVA"+NNsuffix[UseCase]+"_El_4j").c_str());
    TH1D *Rv_el = (TH1D*)inputfiles[*idx][1]->Get(("MVA"+NNsuffix[UseCase]+"_El_4j_ShapeCorrected").c_str());
    
    int CutBin = (UseCase==0 ? 233 : 237);
    
    if ((channel&1)!=0) {
      Num_mu_nominal += Rv_mu_nominal->Integral(/*CutBin*/Rv_mu_nominal->FindBin((UseCase==0?0.9:0.95)),Rv_mu_nominal->GetNbinsX()+1) *(w_v[w_idx]/Rv_mu_nominal->Integral(0,-1));
      Num_mu         += Rv_mu->Integral(/*CutBin*/Rv_mu->FindBin((UseCase==0?0.9:0.95)),Rv_mu->GetNbinsX()+1) *(w_v[w_idx]/Rv_mu->Integral(0,-1));
      Den_mu_nominal += Rv_mu_nominal->Integral(0,Rv_mu_nominal->GetNbinsX()+1) *(w_v[w_idx]/Rv_mu_nominal->Integral(0,-1));
      Den_mu         += Rv_mu->Integral(0,Rv_mu->GetNbinsX()+1) *(w_v[w_idx]/Rv_mu->Integral(0,-1));
      w_idx++;
    }
    if ((channel&2)!=0) {
      Num_el_nominal += Rv_el_nominal->Integral(/*CutBin*/Rv_el_nominal->FindBin((UseCase==0?0.9:0.95)),Rv_el_nominal->GetNbinsX()+1) *(w_v[w_idx]/Rv_el_nominal->Integral(0,-1));
      Num_el         += Rv_el->Integral(/*CutBin*/Rv_el->FindBin((UseCase==0?0.9:0.95)),Rv_el->GetNbinsX()+1) *(w_v[w_idx]/Rv_el->Integral(0,-1));
      Den_el_nominal += Rv_el_nominal->Integral(0,Rv_el_nominal->GetNbinsX()+1) *(w_v[w_idx]/Rv_el_nominal->Integral(0,-1));
      Den_el         += Rv_el->Integral(0,Rv_el->GetNbinsX()+1) *(w_v[w_idx]/Rv_el->Integral(0,-1));
      w_idx++;
    }
    Num_nominal    = Num_mu_nominal+Num_el_nominal;
    Num            = Num_mu+Num_el;
    
    Den_nominal    = Den_mu_nominal+Den_el_nominal;
    Den            = Den_mu+Den_el;
  }
  
  cout<<"R_V for 5D reweighting"<<endl;
                                  
  cout<<"Rv_mu (nominal) ="<<Num_mu_nominal<<"/"<<Den_mu_nominal<<" = "<<Num_mu_nominal/Den_mu_nominal<<endl;
  cout<<"Rv_mu (reweighted) ="<<Num_mu<<"/"<<Den_mu<<" = "<<Num_mu/Den_mu<<endl;

  cout<<"Rv_el (nominal) ="<<Num_el_nominal<<"/"<<Den_el_nominal<<" = "<<Num_el_nominal/Den_el_nominal<<endl;
  cout<<"Rv_el (reweighted) ="<<Num_el<<"/"<<Den_el<<" = "<<Num_el/Den_el<<endl;

  cout<<"Rv (mu+el) (nominal) ="<<Num_nominal<<"/"<<Den_nominal<<" = " <<Num_nominal/Den_nominal<<endl;
  cout<<"Rv (mu+el) (reweighted) ="<<Num<<"/"<<Den<<" = " <<Num/Den<<endl;;
  
  float syst_Rv_mu = (Num_mu/Den_mu) - (Num_mu_nominal/Den_mu_nominal);
  cout<<"syst_Rv_mu = "<<syst_Rv_mu<<" / rel. uncert = "<<syst_Rv_mu/(Num_mu_nominal/Den_mu_nominal)*100<<"%"<<endl;
  float syst_Rv_el = (Num_el/Den_el) - (Num_el_nominal/Den_el_nominal);
  cout<<"syst_Rv_el = "<<syst_Rv_el<<" / rel. uncert = "<<syst_Rv_el/(Num_el_nominal/Den_el_nominal)*100<<"%"<<endl;

//  float syst_Rv_comb = sqrt(pow(w_ttjets[0]*syst_Rtt_mu,2)+pow(w_ttjets[1]*syst_Rtt_el,2));
  float syst_Rv_comb = sqrt(((channel&1)!=0?pow(syst_Rv_mu,2.):0.) + ((channel&2)!=0?pow(syst_Rv_el,2.):0.));
  cout<<"syst_Rv_comb = "<<syst_Rv_comb<<endl;

  if ( (((channel&1)==0)||(((channel&1)!=0)&&(syst_Rv_mu<=0.))) && (((channel&2)==0)||(((channel&2)!=0)&&(syst_Rv_el<=0.))) ) {
    syst_Rv_comb *= -1. ;
  } else if ( (((channel&1)==0)||(((channel&1)!=0)&&(syst_Rv_mu>=0.))) && (((channel&2)==0)||(((channel&2)!=0)&&(syst_Rv_el>=0.))) ) {
    syst_Rv_comb *= 1. ;
  } else {
    printf("\n\n\nAttention (error, failure, crash, et cetaera) !!!!!!! : quadratic sum of opposite directions for 5D reweighting\n\n\n");
    if (syst_Rv_mu+syst_Rv_el<0.) {
      syst_Rv_comb *= -1. ;
    }
//    exit(1);
  }


  //y[10] = Num/Den;
  y[3] = y[0]+syst_Rv_comb;
  cout<<"Relative syst. uncertainty on Rv(comb. channel) : "<<fabs((y[3]-y[0])/y[0])*100<<"%"<<endl;
  yerr_low[3] = 0.00; // TMP
  yerr_high[3] = 0.00; // TMP

  TFile *outputfile = new TFile((std::string()+"RV_"+chanName[channel]+"_"+teffname[UseCase]+"_"+comment+"_Output.root").c_str(),"RECREATE");


  /**************************************************************/
  // systematic Uncert. summary
  /**************************************************************/
  
  TH1F *h  = new TH1F("h","",npoints,0,npoints);
  TH1F *h1 = new TH1F("h1","",npoints,0,npoints);
  TAxis *Xaxis = h->GetXaxis();
  TAxis *Yaxis = h->GetYaxis();
  double BinContent = 0;
  double BinError   = 0;
  double CombUncer_low  = 0, CombUncer_error_low   = 0;
  double CombUncer_high = 0, CombUncer_error_high  = 0;
  for (int j=1;j<npoints;j++){

    BinContent = (y[j]-y[0])/y[0];
    BinError = (y[j]/y[0])*sqrt(pow(yerr_high[j]/y[j],2)+pow(yerr_high[0]/y[0],2));

    Xaxis->SetBinLabel(j,binlabel[j-1].c_str());

    if(j==3){
      h->SetBinContent(j,fabs(BinContent));
      h1->SetBinContent(j,-fabs(BinContent));
      CombUncer_high       += pow(BinContent,2);//-pow(BinError,2);
      CombUncer_error_high += pow(BinError,2);
      CombUncer_low        += pow(BinContent,2);//-pow(BinError,2);
      CombUncer_error_low  += pow(BinError,2);
    }
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
  h->SetBinContent(npoints, CombUncer_high);
  // Adding statistic uncertainties
  cout<<"Rel. stat. uncert. : "<<yerr_high[0]/y[0]<<endl;
  float TotUncer_high = sqrt(pow(CombUncer_high,2)+pow(yerr_high[0]/y[0],2));
  cout<<"Rel. Comb. (+) (stat.+syst.) uncert : "<<TotUncer_high<<endl;

  CombUncer_low = (CombUncer_low >0? sqrt(CombUncer_low)  : 0);
  CombUncer_error_low = sqrt(CombUncer_error_low);
  cout<<"Rel. Comb. (-) syst. uncert : "<<CombUncer_low<<"+/-"<<CombUncer_error_low<<endl;
  h1->SetBinContent(npoints,-CombUncer_low);
  // Adding statistic uncertainties
  cout<<"Rel. stat. uncert. : "<<yerr_high[0]/y[0]<<endl;
  float TotUncer_low = sqrt(pow(CombUncer_low,2)+pow(yerr_low[0]/y[0],2));
  cout<<"Rel. Comb. (-) (stat.+syst.) uncert : "<<TotUncer_low<<endl;

  cout<<"*********** Results *****************"<<endl;
  cout<<std::fixed<<setprecision(6);
  cout<<"$R_{V}$ & $"<<y[0]<<"\\pm"<<yerr_low[0]<<"^{+"<<y[0]*CombUncer_high<<"}_{-"<<y[0]*CombUncer_low<<"}$"<<endl;
  cout<<"$R_{V}$ & $"<<y[0]<<"\\pm"<<yerr_high[0]<<"^{+"<<y[0]*CombUncer_high<<"}_{-"<<y[0]*CombUncer_low<<"}$"<<endl;
  cout<<"*************************************"<<endl;
  
  h->GetXaxis()->SetBinLabel(npoints,binlabel[npoints-1].c_str());
  h->GetXaxis()->LabelsOption("d");
  h->GetXaxis()->SetLabelSize(0.05);
  //h->GetYaxis()->SetTitle("R_{TT}");
  h->GetYaxis()->SetTitle("(R_{V}-R_{V}^{Nominal})/R_{V}^{Nominal}");
  h->GetYaxis()->SetNdivisions(505);
  //h->GetYaxis()->SetRangeUser(0.0004,0.00119);


  TCanvas *c = new TCanvas("c","",600,600);
  c->cd();


  gStyle->SetHistMinimumZero();

  h->SetFillColor(kBlue);
  h->GetYaxis()->SetRangeUser(-0.5,0.5);
  h->Draw("E1Bar");

  h1->SetFillColor(kRed);
  h1->GetYaxis()->SetRangeUser(h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
  h1->Draw("E1Barsame");
  //h1->Draw("2Bar");
/*
  TLine *line = new TLine(0,1,npoints,1);
  line->SetHorizontal();
  line->SetLineColor(kRed);
  line->Draw("same");
*/

  TLatex *tex = new TLatex(0.14,0.995,"CMS Preliminary, 5.0 fb^{-1} at #sqrt{s} = 7 TeV");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();

  c->SetGridy();
  c->Print(("RV_UncertSummary_"+chanName[channel]+"_"+teffname[UseCase]+comment+".eps").c_str());
  c->Write(("RV_UncertSummary_"+chanName[channel]+"_"+teffname[UseCase]+comment).c_str());

  outputfile->Close();

/*
  tex = new TLatex(0.14,1.005,"L_{int} = 5.0 fb^{-1}, #sqrt{s} = 7 TeV");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();
*/
//  TGraphAsymmErrors *tg = new TGraphAsymmErrors(npoints,x,y,0,0,yerr_low,yerr_high);
//  tg->Draw("P");
//  tg->SetMarkerStyle(21);
  inputfiles[0][0]->Close();
  inputfiles[0][1]->Close();
  inputfiles[1][0]->Close();
  inputfiles[1][1]->Close();
}
