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


void makeSystUncertSummaryPlot(UInt_t UseCase, UInt_t channel){
    printf("Starting the program\n");

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  printf("Style options set\n");

  const int NbOfUseCase = 4;
  // Use case : 
  // 0 : low DM
  // 1 : Int DM 1
  // 2 : Int DM 2
  // 3 : High DM
  //const int UseCase = 0;
  string NNsuffix[4] = {"","_IM1","_IM2","_HM"};
  std::string chanName[4] = { "", "SemiMuon", "SemiElectron", "Combined" };
    
  const int NbOfFiles = 11;
  TFile *inputfiles[NbOfFiles][2];
  std::string path = "$HOME/AnalysisCode/GentStopAnalysis/UGentCode/Leg3/v1/Samples_21102012/";
  std::string nompath = "$HOME/AnalysisCode/GentStopAnalysis/UGentCode/Leg3/v1/Samples_15112012/";
//  string path = "$HOME/AnalysisCode/GentStopAnalysis/UGentCode/Leg3/v1/Samples_15112012/";
  
//  inputfiles[0][0] = TFile::Open((path+"TTbar_mu.root").c_str());
//  inputfiles[0][1] = TFile::Open((path+"TTbar_el.root").c_str());
  inputfiles[0][0] = TFile::Open((nompath+"TTbar_mu.root").c_str());
  inputfiles[0][1] = TFile::Open((nompath+"TTbar_el.root").c_str());
  inputfiles[1][0] = TFile::Open((nompath+"TTbar_mu.root").c_str());
  inputfiles[1][1] = TFile::Open((nompath+"TTbar_el.root").c_str());
  inputfiles[2][0] = TFile::Open((nompath+"TTbar_mu.root").c_str());
  inputfiles[2][1] = TFile::Open((nompath+"TTbar_el.root").c_str());
  
  inputfiles[3][0] = TFile::Open((path+"TTbar-Mtop169p5_mu.root").c_str());
  inputfiles[3][1] = TFile::Open((path+"TTbar-Mtop169p5_el.root").c_str());
  inputfiles[4][0] = TFile::Open((path+"TTbar-Mtop175p5_mu.root").c_str());
  inputfiles[4][1] = TFile::Open((path+"TTbar-Mtop175p5_el.root").c_str());
  inputfiles[5][0] = TFile::Open((path+"TTbar-MatchingDown_mu.root").c_str());
  inputfiles[5][1] = TFile::Open((path+"TTbar-MatchingDown_el.root").c_str());
  inputfiles[6][0] = TFile::Open((path+"TTbar-MatchingUp_mu.root").c_str());
  inputfiles[6][1] = TFile::Open((path+"TTbar-MatchingUp_el.root").c_str());
  inputfiles[7][0] = TFile::Open((path+"TTbar-ScaleDown_mu.root").c_str());
  inputfiles[7][1] = TFile::Open((path+"TTbar-ScaleDown_el.root").c_str());
  inputfiles[8][0] = TFile::Open((path+"TTbar-ScaleUp_mu.root").c_str());
  inputfiles[8][1] = TFile::Open((path+"TTbar-ScaleUp_el.root").c_str());
  inputfiles[9][0] = TFile::Open((nompath+"TTbar_mu.root").c_str());
  inputfiles[9][1] = TFile::Open((nompath+"TTbar_el.root").c_str());
  inputfiles[10][0]= TFile::Open((nompath+"TTbar_mu.root").c_str());
  inputfiles[10][1]= TFile::Open((nompath+"TTbar_el.root").c_str());

  int bin = 0;
  const int npoints = NbOfFiles;
  double x[npoints]={0}, y[npoints]={0}, yerr_low[npoints]={0}, yerr_high[npoints]={0};
  //string binlabel[npoints+1] = {"Nominal","JES (-)","JES (+)","m_{t}=173.3 (-1#sigma)","m_{t}=173.3 (+1#sigma)","Matching (-)","Matching (+)","Scale Q^{2} (-)","Scale Q^{2} (+)","SF veto","OF veto","Total syst."};
  string binlabel[npoints+1] = {"Nominal","JES (-)","JES (+)","m_{t}=173.3 (-1#sigma)","m_{t}=173.3 (+1#sigma)","Matching (-)","Matching (+)","Scale Q^{2} (-)","Scale Q^{2} (+)","Lept. SF","5D rew.","Total sys."};
  //char name[100];
  
  string teffname[NbOfUseCase] = {"Eff_LowDM", "Eff_IntDM1", "Eff_IntDM2", "Eff_HighDM"};
  string effname = "";
  string comment = "";

  TFile *outputfile = new TFile((std::string()+"RTT_"+chanName[channel]+"_"+teffname[UseCase]+"_"+comment+"_Output.root").c_str(),"RECREATE");

  TEfficiency *tEff1, *tEff2;
  TList *pList = new TList();
  TGraphAsymmErrors *tg = 0;
  std::vector<Double_t> w_ttjets(2,0.);
  printf("vector\n");
  w_ttjets[0] = 26500./(26500.+19666.);//((TH1D*)inputfiles[0][0]->Get("Entries"))->GetEntries();//165;
  w_ttjets[1] = 19666./(26500.+19666.);//((TH1D*)inputfiles[0][1]->Get("Entries"))->GetEntries();//165;
    printf("Before looping on the files") ;

  for(int i=0;i<npoints-1;i++){

    effname = teffname[UseCase];

    if(i==1) effname += "_m";
    else if(i==2) effname += "_p";
    else if(i==9) effname += "_LW";    
    //else if(i==9) effname += "_noSFveto";
    //else if(i==10)effname += "_noOFveto";
    
    tEff1 = (TEfficiency*)inputfiles[i][0]->Get(effname.c_str());
    tEff1->SetStatisticOption(TEfficiency::kBBayesian);
    //tEff2 = ( inputfiles[i][1] ? (TEfficiency*)inputfiles[i][1]->Get(effname.c_str()) : 0 );
    tEff2 = (TEfficiency*)inputfiles[i][1]->Get(effname.c_str());
    tEff2->SetStatisticOption(TEfficiency::kBBayesian);

    if(tEff1==0 || tEff2==0) continue;
/*
    tEff1->Add(*tEff2);
    y[i]         = tEff1->GetEfficiency(bin);
    yerr_low[i]  = tEff1->GetEfficiencyErrorLow(bin);
    yerr_high[i] = tEff1->GetEfficiencyErrorUp(bin);
*/
    if ((channel&1)!=0) {
      pList->Add((TEfficiency*)tEff1);
      //w_ttjets.push_back(26500./tEff1->GetTotalHistogram()->GetEntries());///((TH1D*)inputfiles[0][0]->Get("Entries"))->GetEntries();//165;
    }
    if ((channel&2)!=0) {
      pList->Add((TEfficiency*)tEff2);
      //w_ttjets.push_back(19666./tEff2->GetTotalHistogram()->GetEntries());///((TH1D*)inputfiles[0][1]->Get("Entries"))->GetEntries();//165;
    }
    tg = (TGraphAsymmErrors*)TEfficiency::Combine(pList,"mode",((channel+1)/2)*2);//,&w_ttjets[0]);
    tg->GetPoint(bin,x[i],y[i]);
    yerr_low[i]  = tg->GetErrorYlow(bin);
    yerr_high[i] = tg->GetErrorYhigh(bin);
    if(i==0){
      cout<<"Rtt(mu,nominal) = "<<tEff1->GetEfficiency(bin+1)<<" +/- "<<tEff1->GetEfficiencyErrorUp(bin+1)<<"/"<<tEff1->GetEfficiencyErrorLow(bin+1)<<endl;
      cout<<"Rtt(el,nominal) = "<<tEff2->GetEfficiency(bin+1)<<" +/- "<<tEff2->GetEfficiencyErrorUp(bin+1)<<"/"<<tEff2->GetEfficiencyErrorLow(bin+1)<<endl;
      cout<<"Rtt(comb,nominal) = "<<y[i]<<" +/- "<<yerr_high[i]<<"/"<<yerr_low[i]<<endl;
    }

    tEff1->Delete();
    tEff2->Delete();

    tg->Delete();
    pList->Clear();
  }
  
  int CutBin = (UseCase==0 ? 233 : 237);
  
  TH1D *Rtt_mu_nominal = (TH1D*)inputfiles[10][0]->Get(("MVA"+NNsuffix[UseCase]+"_Mu_4j").c_str());
  TH1D *Rtt_mu = (TH1D*)inputfiles[10][0]->Get(("MVA"+NNsuffix[UseCase]+"_Mu_4j_ShapeCorrected").c_str());
  TH1D *Rtt_el_nominal = (TH1D*)inputfiles[10][1]->Get(("MVA"+NNsuffix[UseCase]+"_El_4j").c_str());
  TH1D *Rtt_el = (TH1D*)inputfiles[10][1]->Get(("MVA"+NNsuffix[UseCase]+"_El_4j_ShapeCorrected").c_str());
	
  float Num_mu_nominal = Rtt_mu_nominal->Integral(CutBin,Rtt_mu_nominal->GetNbinsX()+1);
  float Num_el_nominal = Rtt_el_nominal->Integral(CutBin,Rtt_el_nominal->GetNbinsX()+1);
  float Num_nominal    = Num_mu_nominal+Num_el_nominal;
  float Num_mu         = Rtt_mu->Integral(CutBin,Rtt_mu->GetNbinsX()+1);
  float Num_el         = Rtt_el->Integral(CutBin,Rtt_el->GetNbinsX()+1);
  float Num            = Num_mu+Num_el;

  float Den_mu_nominal = Rtt_mu_nominal->Integral(0,Rtt_mu_nominal->GetNbinsX()+1);
  float Den_el_nominal = Rtt_el_nominal->Integral(0,Rtt_el_nominal->GetNbinsX()+1);
  float Den_nominal    = Den_mu_nominal+Den_el_nominal;
  float Den_mu         = Rtt_mu->Integral(0,Rtt_mu->GetNbinsX()+1);
  float Den_el         = Rtt_el->Integral(0,Rtt_el->GetNbinsX()+1);
  float Den            = Den_mu+Den_el;
  
  //Rtt_mu->Scale(w_ttjets[0]/Den_mu);
  //Rtt_el->Scale(w_ttjets[1]/Den_el);
  
  cout<<"Rtt_mu (nominal) ="<<Num_mu_nominal<<"/"<<Den_mu_nominal<<" = "<<Num_mu_nominal/Den_mu_nominal<<endl;
  cout<<"Rtt_mu (reweighted) ="<<Num_mu<<"/"<<Den_mu<<" = "<<Num_mu/Den_mu<<endl;

  cout<<"Rtt_el (nominal) ="<<Num_el_nominal<<"/"<<Den_el_nominal<<" = "<<Num_el_nominal/Den_el_nominal<<endl;
  cout<<"Rtt_el (reweighted) ="<<Num_el<<"/"<<Den_el<<" = "<<Num_el/Den_el<<endl;

  cout<<"Rtt (mu+el) (nominal) ="<<Num_nominal<<"/"<<Den_nominal<<" = " <<Num_nominal/Den_nominal<<endl;
  cout<<"Rtt (mu+el) (reweighted) ="<<Num<<"/"<<Den<<" = " <<Num/Den<<endl;;
  
  float syst_Rtt_mu = (Num_mu/Den_mu) - (Num_mu_nominal/Den_mu_nominal);
  cout<<"syst_Rtt_mu = "<<syst_Rtt_mu<<" / rel. uncert = "<<syst_Rtt_mu/(Num_mu_nominal/Den_mu_nominal)*100<<"%"<<endl;
  float syst_Rtt_el = (Num_el/Den_el) - (Num_el_nominal/Den_el_nominal);
  cout<<"syst_Rtt_el = "<<syst_Rtt_el<<" / rel. uncert = "<<syst_Rtt_el/(Num_el_nominal/Den_el_nominal)*100<<"%"<<endl;
  float syst_Rtt_comb = sqrt(((channel&1)!=0?pow(w_ttjets[0]*syst_Rtt_mu,2):0.) + ((channel&2)!=0?pow(w_ttjets[1]*syst_Rtt_el,2):0.));
  if ( (((channel&1)==0)||(((channel&1)!=0)&&(syst_Rtt_mu<=0.))) && (((channel&2)==0)||(((channel&2)!=0)&&(syst_Rtt_el<=0.))) ) {
    syst_Rtt_comb *= -1. ;
  } else if ( (((channel&1)==0)||(((channel&1)!=0)&&(syst_Rtt_mu>=0.))) && (((channel&2)==0)||(((channel&2)!=0)&&(syst_Rtt_el>=0.))) ) {
    syst_Rtt_comb *= 1. ;
  } else {
    printf("\n\n\nAttention (error, failure, crash, et cetaera) !!!!!!! : quadratic sum of opposite directions for 5D reweighting\n\n\n");
    if (syst_Rtt_mu+syst_Rtt_el<0.) {
      syst_Rtt_comb *= -1. ;
    }
//    exit(1);
  }

  cout<<"syst_Rtt_comb = "<<syst_Rtt_comb<<endl;

  //y[10] = Num/Den;
  y[10] = y[0]+syst_Rtt_comb;
  cout<<"Relative syst. uncertainty on Rtt(comb. channel) : "<<fabs((y[10]-y[0])/y[0])*100<<"%"<<endl;
  yerr_low[10] = 0.00; // TMP
  yerr_high[10] = 0.00; // TMP
  /****************************************/
  // Uncert. related to the top quark mass
  /****************************************/
  TCanvas *c = new TCanvas("c","",600,600);
  c->cd();

  TGraphErrors *mtop = new TGraphErrors(3);//"","",3,168,177);
  TF1 *f1 = new TF1("f1","[0]+[1]*x");

  mtop->SetPoint(1,172.5,0.);
  mtop->SetPointError(1,1.5,0.);
  mtop->SetPoint(0,169.5,(y[3]-y[0])/y[0]);
  mtop->SetPointError(0,1.5,(y[3]/y[0])*sqrt(pow(yerr_high[3]/y[3],2)+pow(yerr_high[0]/y[0],2)));
  mtop->SetPoint(2,175.5,(y[4]-y[0])/y[0]);
  mtop->SetPointError(2,1.5,(y[4]/y[0])*sqrt(pow(yerr_high[4]/y[4],2)+pow(yerr_high[0]/y[0],2)));
  
  mtop->Draw("AP");
  mtop->GetXaxis()->SetTitle("m_{top quark} [GeV/c^{2}]");
  mtop->GetYaxis()->SetTitle("(R_{TT}-R_{TT}^{Nominal})/R_{TT}^{Nominal}");

  mtop->GetXaxis()->SetRangeUser(168,177);
//  mtop->GetYaxis()->SetRangeUser(-0.3,0.1);

  f1->SetParameter(0,-6.);
  f1->SetParameter(1,0.03);

  //mtop->SetStats(false);
  mtop->Fit(f1,"LW");
  //mtop->GetXaxis()->CenterLabels();
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  assert(fitter != 0);
  fitter->PrintResults(0,0);
  //double * covMatrix = fitter->GetCovarianceMatrixElement(0,1);
  double covMatrixElem00 = fitter->GetCovarianceMatrixElement(0,0);
  double covMatrixElem01 = fitter->GetCovarianceMatrixElement(0,1);
  double covMatrixElem10 = fitter->GetCovarianceMatrixElement(1,0);
  double covMatrixElem11 = fitter->GetCovarianceMatrixElement(1,1);
  cout<<"Element 0/0 : "<<covMatrixElem00<<endl;
  cout<<"Element 0/1 : "<<covMatrixElem01<<endl;
  cout<<"Element 1/0 : "<<covMatrixElem10<<endl;
  cout<<"Element 1/1 : "<<covMatrixElem11<<endl;

  float y_min = mtop->GetYaxis()->GetXmin();

  float y_mean = f1->Eval(173.3);
  float y_m = f1->Eval(171.9);
  float y_p = f1->Eval(174.7);

  float y_m_err = sqrt(pow(171.9*f1->GetParError(1),2)+pow(f1->GetParError(0),2)+2*171.9*covMatrixElem01);
  float y_p_err = sqrt(pow(174.7*f1->GetParError(1),2)+pow(f1->GetParError(0),2)+2*174.7*covMatrixElem01);

  TLine *lmtop_m = new TLine(171.9,y_min,171.9,y_m);
  TLine *lmtop_p = new TLine(174.7,y_min,174.7,y_p);
  TLine *lerr_m  = new TLine(168,y_m,171.9,y_m);
  TLine *lerr_p  = new TLine(168,y_p,174.7,y_p);
  
  lmtop_m->SetLineStyle(3);
  lmtop_p->SetLineStyle(3);
  lerr_m->SetLineStyle(3);
  lerr_p->SetLineStyle(3);

  lmtop_m->Draw();
  lmtop_p->Draw();
  lerr_m->Draw();
  lerr_p->Draw();
  
  TLatex *latex = new TLatex(0.14,0.92,"LHC top quark mass = (173.3#pm1.4) GeV/c^{2}");
  latex->SetNDC();
  latex->SetTextAlign(13);
  latex->SetTextFont(42);
  latex->SetTextSize(0.045);
  latex->SetLineWidth(2);
  latex->Draw();

  TLatex *tex = new TLatex(/*0.65*/0.14,0.995,"CMS Preliminary, 5.0 fb^{-1} at #sqrt{s} = 7 TeV");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
  tex->Draw();

  c->Print(("Mtop_Uncert"+chanName[channel]+"_"+teffname[UseCase]+comment+".eps").c_str());
  c->Write(("Mtop_Uncert_"+chanName[channel]+"_"+teffname[UseCase]+comment).c_str());

  cout<<"Rel. Uncertainty on Rtt : ["<<y_m-y_mean<<","<<y_p-y_mean<<"]"<<endl;
  cout<<"Errors on y_m : "<<y_m_err<<endl;
  cout<<"Errors on y_p : "<<y_p_err<<endl;

  /**************************************************************/
  // systematic Uncert. summary
  /**************************************************************/
  
  TH1F *h  = new TH1F("h","",npoints,0,npoints);
  TH1F *h1 = new TH1F("h1","",npoints,0,npoints);
  TAxis *Xaxis = h->GetXaxis();
  //TAxis *Yaxis = h->GetYaxis();
  double BinContent = 0;
  double BinError   = 0;
  double CombUncer_low  = 0, CombUncer_error_low   = 0;
  double CombUncer_high = 0, CombUncer_error_high  = 0;
  for (int j=1;j<npoints;j++){

    if(j==3)      BinContent = y_m-y_mean;
    else if(j==4) BinContent = y_p-y_mean;
    //else if(j==5) BinContent = (y[0]-y[j])/y[0];
    else          BinContent = (y[j]-y[0])/y[0];

    if(j==3)      BinError = y_m_err;
    else if(j==4) BinError = y_p_err;
    else          BinError = (y[j]/y[0])*sqrt(pow(yerr_high[j]/y[j],2)+pow(yerr_high[0]/y[0],2));

    Xaxis->SetBinLabel(j,binlabel[j].c_str());

    if(j==10){
      h->SetBinContent(j,fabs(BinContent));
      h1->SetBinContent(j,-fabs(BinContent));
      CombUncer_high       += pow(BinContent,2);//-pow(BinError,2);
      CombUncer_error_high += pow(BinError,2);
      CombUncer_low        += pow(BinContent,2);//-pow(BinError,2);
      CombUncer_error_low  += pow(BinError,2);
    }
    else if(BinContent>0){
      h->SetBinContent(j,BinContent);
      h->SetBinError(j,BinError);
      cout<<"Bin+ "<<j<<", Content : "<<BinContent<<" / Error : "<<BinError<<" / Syst² : "<<pow(BinContent,2)-pow(BinError,2)<<endl;
      CombUncer_high += pow(BinContent,2);//-pow(BinError,2);
      CombUncer_error_high += pow(BinError,2);
    }
    else{
      h1->SetBinContent(j,BinContent);
      h1->SetBinError(j,BinError);
      cout<<"Bin- "<<j<<", Content : "<<BinContent<<" / Error : "<<BinError<<" / Syst² : "<<pow(BinContent,2)-pow(BinError,2)<<endl;
      CombUncer_low       += pow(BinContent,2);//-pow(BinError,2);
      CombUncer_error_low += pow(BinError,2);
    }
  }

  // Systematic uncertainties
  CombUncer_high = (CombUncer_high>0? sqrt(CombUncer_high) : 0);
  CombUncer_error_high = sqrt(CombUncer_error_high); 
  cout<<"Rel. Comb. (+) syst. uncert : "<<CombUncer_high<<"+/-"<<CombUncer_error_high<<endl;
  h->SetBinContent(npoints, CombUncer_high);
  //h->SetBinError(npoints+1, CombUncer_error_high);

  // Adding statistic uncertainties
  cout<<"Rel. stat. uncert. : "<<yerr_high[0]/y[0]<<endl;
  float TotUncer_high = sqrt(pow(CombUncer_high,2)+pow(yerr_high[0]/y[0],2));
  cout<<"Rel. Comb. (+) (stat.+syst.) uncert : "<<TotUncer_high<<endl;

  CombUncer_low = (CombUncer_low >0? sqrt(CombUncer_low)  : 0);
  CombUncer_error_low = sqrt(CombUncer_error_low);
  cout<<"Rel. Comb. (-) syst. uncert : "<<CombUncer_low<<"+/-"<<CombUncer_error_low<<endl;
  h1->SetBinContent(npoints,-CombUncer_low);
  //h1->SetBinError(npoints+1, CombUncer_error_low);

  // Adding statistic uncertainties
  cout<<"Rel. stat. uncert. : "<<yerr_high[0]/y[0]<<endl;
  float TotUncer_low = sqrt(pow(CombUncer_low,2)+pow(yerr_low[0]/y[0],2));
  cout<<"Rel. Comb. (-) (stat.+syst.) uncert : "<<TotUncer_low<<endl;

  cout<<"*********** Results *****************"<<endl;
  cout<<std::fixed<<setprecision(6);
  cout<<"$R_{TT}$ & $"<<y[0]<<"\\pm"<<yerr_low[0]<<"^{+"<<y[0]*CombUncer_high<<"}_{-"<<y[0]*CombUncer_low<<"}$"<<endl;
  cout<<"$R_{TT}$ & $"<<y[0]<<"\\pm"<<yerr_high[0]<<"^{+"<<y[0]*CombUncer_high<<"}_{-"<<y[0]*CombUncer_low<<"}$"<<endl;
  cout<<"*************************************"<<endl;

  h->GetXaxis()->SetBinLabel(npoints,binlabel[npoints].c_str());
  h->GetXaxis()->LabelsOption("d");
  h->GetXaxis()->SetLabelSize(0.04);
  //h->GetYaxis()->SetTitle("R_{TT}");
  h->GetYaxis()->SetTitle("(R_{TT}-R_{TT}^{Nominal})/R_{TT}^{Nominal}");
  h->GetYaxis()->SetNdivisions(505);
  //h->GetYaxis()->SetRangeUser(0.0004,0.00119);

  gStyle->SetHistMinimumZero();

  h->SetFillColor(kBlue);
  float max = h->GetMaximum()*1.1;
  float min = h1->GetMinimum()*1.1;
  float range = (fabs(min)<fabs(max) ? fabs(max) : fabs(min));
  h->GetYaxis()->SetRangeUser(-range,range);
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
/*
  TLatex *tex = new TLatex(0.14,0.995,"CMS Preliminary, 5.0 fb^{-1} at #sqrt{s} = 7 TeV");
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextFont(42);
  tex->SetTextSize(0.045);
  tex->SetLineWidth(2);
*/
  tex->Draw();
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
  c->SetGridy();
  c->Print(("RTT_UncertSummary_"+chanName[channel]+"_"+teffname[UseCase]+comment+".eps").c_str());
  c->Write(("RTT_UncertSummary_"+chanName[channel]+"_"+teffname[UseCase]+comment).c_str());
//  c->Close();
  
  outputfile->Close();
}
