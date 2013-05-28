#include "TStyle.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

  // Root headers
#include "TArrow.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TEntryList.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2I.h"
#include "THStack.h"
#include "TMath.h"
#include "TMarker.h"
#include "TPaveStats.h"

#include "TFitter.h"
  // RooFit headers
#include "RooAddition.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooFormula.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"

#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooProdPdf.h"
#include "RooSimultaneous.h"

#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooNLLVar.h"

#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooTable.h"
#include "Roo1DTable.h"
#include "TF1.h"

#include "TChainElement.h"

using namespace RooFit;
using namespace std;

  //#include "MistagFuncs.C"
#include "Leg3_8TeV/include/BScaleFac.h"
#include "VJetEstimation.h"
  //  #include "ScaleFactors_btagPOG_Moriond2013/SFlightFuncs_Moriond2013.C"
  //  #include "ScaleFactors_btagPOG_Moriond2013/SFb_codified.C"

/**
 Results from BTV-11-003, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG link to https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/eff_b_c-ttbar_payload.txt 
 In Range : 0<=|eta|<=2.4 , 30<=p_T<=200  (approximately)
 */
/** TCHE */
/** //7 TeV only
 Float_t eff_b_ttbar_FtCM(Float_t /*b-disc* / x) {
 return -3.67153247396e-07*x*x*x*x +  -2.81599797034e-05*x*x*x +  0.00293190163243*x*x +  -0.0849600849778*x +  0.928524440715 ;
 };
 Float_t eff_b_ttbar_FtCM_pluserr(Float_t /*b-disc* / x) {
 return 3.03337430722e-06*x*x*x*x + -0.000171604835897*x*x*x + 0.00474711667943*x*x + -0.0929933040514*x + 0.978347619293 ;
 };
 Float_t eff_b_ttbar_MC(Float_t /*b-disc* / x) {
 return 3.90732786802e-06*x*x*x*x +  -0.000239934437355*x*x*x +  0.00664986827287*x*x +  -0.112578996016*x +  1.00775721404 ;
 };
 Float_t eff_c_ttbar_MC(Float_t /*b-disc* / x) {
 return 0.343760640168*exp(-0.00315525164823*x*x*x + 0.0805427315196*x*x + -0.867625139194*x + 1.44815935164 ) ;
 };
 */
const Int_t ptBins_NUM = 16;
float ptBins[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};

  // BTV-11-004
const Int_t mistag_eta_binning_loose_NUM = 5;
Double_t mistag_eta_binning_loose[mistag_eta_binning_loose_NUM]  = { 0., 0.5, 1., 1.5, 2.4 };
const Int_t mistag_eta_binning_medium_NUM = 4;
Double_t mistag_eta_binning_medium[mistag_eta_binning_medium_NUM] = { 0., 0.8, 1.6, 2.4 };
const Int_t mistag_eta_binning_tight_NUM = 2;
Double_t mistag_eta_binning_tight[mistag_eta_binning_tight_NUM]   = { 0., 2.4 };
  // + [0 ; 2.4] ;



/*
 TF1** mistag_light_max = new TF1*[mistag_eta_binning_medium_NUM-1];
 for (Int_t i=1; i<mistag_eta_binning_medium_NUM; i++) {
 mistag_light_max[i-1] = NULL;
 mistag_light_max[i-1] = GetMistagmax("TCHE","M",mistag_eta_binning_medium[i-1], mistag_eta_binning_medium[i]);
 }
 TF1** mistag_light_min = new TF1*[mistag_eta_binning_medium_NUM-1];
 for (Int_t i=1; i<mistag_eta_binning_medium_NUM; i++) {
 mistag_light_min[i-1] = NULL;
 mistag_light_min[i-1] = GetMistagmin("TCHE","M",mistag_eta_binning_medium[i-1], mistag_eta_binning_medium[i]);
 }
 */

const char* formatValErr(Double_t val, Double_t err)
{
  Int_t b=(Int_t) floor(log10(fabs(err))) ;
  Int_t c=(Int_t) floor(log10(fabs(val))) ;
  Int_t di=0;
  if (c-b >0) {
    di=c-b;
  }
  char *format = (char*) calloc((Int_t)(17+100+floor(di)), sizeof(char)) ;
  if (err/pow(10.,b)<5.) {
    //printf(" %.3lg \\pm %lg ", val, err);
    sprintf(format, "%%.%dlg \\pm %%.2lg", di+2);
  } else {
    //printf(" %lg \\pm %lg ", val, err);
    sprintf(format, "%%.%dlg \\pm %%.1lg", di+1);
  }
  return format ;
}

void printfValErr(Double_t val, Double_t err) {
  printf(formatValErr(val, err), val, err);
}

void sprintfValErr(char* string, Double_t val, Double_t err) {
  sprintf(string, formatValErr(val, err), val, err);
}


Float_t WilsonScoreIntervalHigh(Float_t Non, Float_t Ntot)
{
	Double_t T = (Ntot>0 ? 1/Ntot : 0);
	Double_t p_hat = (Ntot>0 && Non>=0 && Ntot>=Non ? Non/Ntot : 0);
	Double_t Int_High = ((p_hat+(T/2))/(1+T))+(sqrt(p_hat*(1-p_hat)*T+pow(T/2,2))/(1+T));
	return Int_High;
}

Float_t WilsonScoreIntervalLow(Float_t Non, Float_t Ntot)
{
	Double_t T = (Ntot>0 ? 1/Ntot : 0);
	Double_t p_hat = (Ntot>0 && Non>=0 && Ntot>=Non ? Non/Ntot : 0);
	Double_t Int_Low = ((p_hat+(T/2))/(1+T))-(sqrt(p_hat*(1-p_hat)*T+pow(T/2,2))/(1+T));
	return Int_Low;
}

Float_t WilsonScoreIntervalMean(Float_t Non, Float_t Ntot)
{
	Double_t Err_High = WilsonScoreIntervalHigh(Non, Ntot)-Non;
	Double_t Err_Low  = Non-WilsonScoreIntervalLow(Non, Ntot);
	return (Err_High+Err_Low)/2;
}

/*
 Compute the average efficiency, with SF applied
 flavour : 0 (guds), 4 (c), 5 (b)
 SFerrMode : -1 (err -), 0 (nominal), +1 (err +)
 */
TEfficiency* average(TEfficiency *eff, Int_t flavour, std::string run_per, Int_t SFerrMode=0, Bool_t isOnMC=kFALSE) {
  TF1* sf = NULL;
    //  printf("%lf   ", eff->GetPassedHistogram()->Integral()); eff->GetPassedHistogram()->Print();
    //  printf("%lf   ", eff->GetTotalHistogram()->Integral());  eff->GetTotalHistogram()->Print();
  //  TAxis *xaxis = eff->GetTotalHistogram()->GetXaxis() ;
  //  TAxis *yaxis = eff->GetTotalHistogram()->GetYaxis() ;
  /*
   TEfficiency *effSF = new TEfficiency((std::string()+eff->GetName()+"_Avg_SFed").c_str(),(std::string()+eff->GetName()+"_SFed").c_str(),
   1, xaxis->GetBinLowEdge(1), xaxis->GetBinUpEdge(xaxis->GetNbins()),
   1, yaxis->GetBinLowEdge(1), yaxis->GetBinUpEdge(yaxis->GetNbins()));
   */
  TEfficiency *effSF = new TEfficiency((std::string()+eff->GetName()+"_Avg_SFed").c_str(),(std::string()+eff->GetName()+"_SFed").c_str(),
                                       1, 0., 1.);
  Double_t avgTotal = 0.;
  Double_t avgPassed = 0.;
  const TH1* passedHisto   = eff->GetPassedHistogram();
  const TH1* totHisto      = eff->GetTotalHistogram();
  Int_t nbinsx = totHisto->GetNbinsX(); //ptBins_NUM
  Int_t nbinsy = totHisto->GetNbinsY(); //mistag_eta_binning_medium_NUM
                                        //  printf("   ( %d , %d ) \n", nbinsx, nbinsy);
  for (Int_t iy=1; iy<=nbinsy; iy++) {
    if (flavour==4 || flavour==5) {
      sf = GetSFHeavy("mean","CSV","M", 0.0, 2.4, "ABCD");
    } else {
      if (SFerrMode == -1) {
        sf = GetSFLight("min","CSV","M", mistag_eta_binning_medium[iy-1], mistag_eta_binning_medium[iy], run_per);
      } else if (SFerrMode == +1) {
        sf = GetSFLight("max","CSV","M", mistag_eta_binning_medium[iy-1], mistag_eta_binning_medium[iy], run_per);
      } else {
        sf = GetSFLight("mean","CSV","M", mistag_eta_binning_medium[iy-1], mistag_eta_binning_medium[iy], run_per);
      } 
    }
    for (Int_t ix=1; ix<=nbinsx; ix++) {
      Double_t ptBin = (ptBins[ix-1]+ptBins[ix])/2. ;
      Double_t scalefac = sf->Eval(ptBin,0,0);
      if (flavour==4 || flavour==5) {
        if (SFerrMode == -1) {
          scalefac -= SFHeavyError(ptBin);
        } else if (SFerrMode == +1) {
          scalefac += SFHeavyError(ptBin);
        }
      }
      if (isOnMC) {
        scalefac = 1. ;
      }
      Int_t bin = eff->FindFixBin(ix-1,iy-1);
      Double_t nPassed = passedHisto->GetBinContent(bin);
        //      passedHisto->SetBinContent(bin, scalefac*nPassed);
      Double_t nTotal = totHisto->GetBinContent(bin);
        //      printf("Fill ( %lf -> %lf ) / %lf\n", nPassed, scalefac * nPassed, nTotal);
      avgPassed += scalefac * nPassed ;
      avgTotal  += nTotal ;
    }
  }
  TH1 *tmp_passed = (TH1*) (effSF->GetPassedHistogram())->Clone();
  tmp_passed->SetBinContent(1 /*tmp_passed->FindFixBin(1,1)*/, avgPassed);
  TH1 *tmp_total  = (TH1*) (effSF->GetTotalHistogram())->Clone();
  tmp_total->SetBinContent(1 /*tmp_total->FindFixBin(1,1)*/, avgTotal);
  /*
   printf("Return from histo : %lf / %lf\n",
   tmp_passed->GetBinContent(4),
   tmp_total->GetBinContent(4) );
   */
  effSF->SetTotalHistogram (*tmp_total ,"f");
  effSF->SetPassedHistogram(*tmp_passed,"f");
    //  effSF->SetPassedHistogram(passedHisto, "f");
  /*
   printf("Return from average : %lf / %lf = %lf\n",
   effSF->GetPassedHistogram()->GetBinContent(4),
   effSF->GetTotalHistogram()->GetBinContent(4),
   effSF->GetEfficiency(4) );
   */
/*
  TH2* ttp = (TH2*) eff->GetPassedHistogram()->Clone();
  ttp->Rebin2D(ttp->GetNbinsX(), ttp->GetNbinsY());
  TH2* ttt = (TH2*) eff->GetTotalHistogram()->Clone();
  ttt->Rebin2D(ttp->GetNbinsX(), ttp->GetNbinsY());
  effSF->SetTotalHistogram(*ttt, "f");
  effSF->SetPassedHistogram(*ttp, "f");
  */
  return effSF;
}

/*
 Apply the bin-per-bin SF
 flavour : 0 (guds), 4 (c), 5 (b)
 */
TEfficiency* applySF(TEfficiency *eff, Int_t flavour, std::string run_per, Int_t SFerrMode=0, Bool_t isOnMC=kFALSE) {
  TF1* sf = NULL;
  TEfficiency *effSF = (TEfficiency*) eff->Clone((std::string()+eff->GetName()+"_SFed").c_str());
  TH1* passedHisto   = (TH1*) effSF->GetPassedHistogram()->Clone();
  const TH1* totHisto      = effSF->GetTotalHistogram();
  Int_t nbinsx = totHisto->GetNbinsX(); //ptBins_NUM
  Int_t nbinsy = totHisto->GetNbinsY(); //mistag_eta_binning_medium_NUM
  for (Int_t iy=1; iy<=nbinsy; iy++) {
    if (flavour==4 || flavour==5) {
      sf = GetSFHeavy("mean","CSV","M", 0.0, 2.4, "ABCD");
    } else {
      if (SFerrMode == -1) {
        sf = GetSFLight("min","CSV","M", mistag_eta_binning_medium[iy-1], mistag_eta_binning_medium[iy], run_per);
      } else if (SFerrMode == +1) {
        sf = GetSFLight("max","CSV","M", mistag_eta_binning_medium[iy-1], mistag_eta_binning_medium[iy], run_per);
      } else {
        sf = GetSFLight("mean","CSV","M", mistag_eta_binning_medium[iy-1], mistag_eta_binning_medium[iy], run_per);
      }
    } 
    for (Int_t ix=1; ix<=nbinsx; ix++) {
      Double_t ptBin = (ptBins[ix-1]+ptBins[ix])/2. ;
      Double_t scalefac = sf->Eval(ptBin,0,0);
      if (flavour==4 || flavour==5) {
        if (SFerrMode == -1) {
          scalefac -= SFHeavyError(ptBin);
        } else if (SFerrMode == +1) {
          scalefac += SFHeavyError(ptBin);
        }
      }
      if (isOnMC) {
        scalefac = 1. ;
      }
      Int_t bin = eff->FindFixBin(ix-1,iy-1);
      Double_t nPassed = passedHisto->GetBinContent(bin);
      passedHisto->SetBinContent(bin, scalefac*nPassed);
    }
  }
  effSF->SetPassedHistogram(*passedHisto, "f");
  return effSF;
}


/*
 apply the SF, on the SF
 flavour : 0 (guds), 4 (c), 5 (b)
 */
TGraphAsymmErrors* tgWithSF(TGraphAsymmErrors *eff, Int_t flavour, std::string run_per, Bool_t isOnMC=kFALSE) {
  TF1* sf = NULL;
  TF1* sfmin = NULL;
  TF1* sfmax = NULL;
  TGraphAsymmErrors *effSF = (TGraphAsymmErrors*) eff->Clone((std::string()+eff->GetName()+"_SFed").c_str());
  Int_t nbinsx = ptBins_NUM ;
  Int_t nbinsy = mistag_eta_binning_medium_NUM ;
  for (Int_t iy=1; iy<=nbinsy; iy++) {
    if (flavour==4 || flavour==5) {
      sf = GetSFHeavy("mean","CSV","M", 0.0, 2.4, "ABCD");
    } else {
      sf = GetSFLight("mean","CSV","M", mistag_eta_binning_medium[iy-1], mistag_eta_binning_medium[iy], run_per);
      sfmin = GetSFLight("min","CSV","M", mistag_eta_binning_medium[iy-1], mistag_eta_binning_medium[iy], run_per);
      sfmax = GetSFLight("max","CSV","M", mistag_eta_binning_medium[iy-1], mistag_eta_binning_medium[iy], run_per);
    }
    for (Int_t ix=1; ix<=nbinsx; ix++) {
      Double_t ptBin = (ptBins[ix-1]+ptBins[ix])/2. ;
      Double_t scalefac = sf->Eval(ptBin,0,0);
      Double_t scalefac_min=0., scalefac_max=0.;
      if (flavour==1 || flavour==2) {
        Double_t sfErr = SFHeavyError(ptBin) ;
        scalefac_min = scalefac - sfErr ;
        scalefac_max = scalefac + sfErr ;
      } else {
        scalefac_min = sfmin->Eval(ptBin,0,0);
        scalefac_max = sfmax->Eval(ptBin,0,0);
      }
      
      Int_t bin = ix*(nbinsx+2) + iy ; //eff->FindFixBin(ix,iy);
      Double_t x=0., y=0. ;
      eff->GetPoint(bin, x, y);
      Double_t errLow  = eff->GetErrorYlow(bin);
      Double_t errHigh = eff->GetErrorYhigh(bin);
      if (isOnMC) {
        scalefac = 1. ;
        effSF->SetPoint(bin, x, y);
        effSF->SetPointError(bin, eff->GetErrorXlow(bin), eff->GetErrorXhigh(bin),
                             errLow,
                             errHigh);
      } else {
        effSF->SetPoint(bin, x, y*scalefac);
        effSF->SetPointError(bin, eff->GetErrorXlow(bin), eff->GetErrorXhigh(bin),
                             sqrt((errLow/y)*(errLow/y)  +(1-scalefac_min/scalefac)*(1-scalefac_min/scalefac))*errLow*scalefac,
                             sqrt((errHigh/y)*(errHigh/y)+(1-scalefac_max/scalefac)*(1-scalefac_max/scalefac))*errHigh*scalefac);
      }
    }
  }
  return effSF;
}


/**
 void FillTriggerPattern(TH1I* histoTrigger, std::string triggerPattern, double weight) 
 {
 TAxis *xAxis = histoTrigger->GetXaxis();
 bool found = false;
 for(int i=1; i<=histoTrigger->GetNbinsX(); i++) {
 std::string name = xAxis->GetBinLabel(i);
 if (name==triggerPattern) {
 histoTrigger->Fill(triggerPattern.c_str(),weight);
 found=true;
 continue;
 }
 }
 if (!found) {
 for(int i=1; i<=histoTrigger->GetNbinsX(); i++) {
 /*
 if (xAxis->GetBinLabel(i)==NULL) {
 printf("Bin label is NULL\n");
 } else {
 printf("Bin label is \"%s\"\n", xAxis->GetBinLabel(i));
 }
 * /
 if (std::string(xAxis->GetBinLabel(i))==std::string("")) {
 xAxis->SetBinLabel(i, triggerPattern.c_str());
 histoTrigger->Fill(triggerPattern.c_str(), weight);
 continue ;
 }
 if (i==histoTrigger->GetNbinsX()) {
 printf("Histogram \"%s\"  not big enough for trigger\n", histoTrigger->GetName());
 //        exit(1);
 }
 }
 }
 };
 */


Int_t StopSearches_VJetsBckgdEst_NoFiles_SevJets(std::string signalFile="", Double_t signal_xs=0.)
{
  /**  //Not provided for 2012 anymore
   TF1** mistag_light_mean = new TF1*[mistag_eta_binning_medium_NUM-1];
   for (Int_t i=1; i<mistag_eta_binning_medium_NUM; i++) {
   mistag_light_mean[i-1] = NULL;
   mistag_light_mean[i-1] = GetMistagmean("TCHE","M",mistag_eta_binning_medium[i-1], mistag_eta_binning_medium[i]);
   }
   
   TF1** mistag_light_mean_WPloose = new TF1*[mistag_eta_binning_loose_NUM-1];
   for (Int_t i=1; i<mistag_eta_binning_loose_NUM; i++) {
   mistag_light_mean_WPloose[i-1] = NULL;
   mistag_light_mean_WPloose[i-1] = GetMistagmean("TCHE","L",mistag_eta_binning_loose[i-1], mistag_eta_binning_loose[i]);
   }
   */
  /*
  printfValErr(-532.46467543, -132.); printf("\n");
  printfValErr(-532.46467543, -13.2); printf("\n");
  printfValErr(-532.46467543, -1.32); printf("\n");
  printfValErr(-532.46467543, -0.000000132); printf("\n");
  printfValErr(-53.246467543, -0.000000132); printf("\n");
  printfValErr(-5.3246467543, -0.000000132); printf("\n");
  printfValErr(-5.2309873209, -.0015121); printf("\n");
  printfValErr(5.2309873209, -.0015121); printf("\n");
  printfValErr(5.2309873209, .0015121); printf("\n");
  printfValErr(-5.002342389, .006243); printf("\n");
  printfValErr(-0.0000001, -5.3209873); printf("\n");
  printfValErr(5., -5.); printf("\n");
  printfValErr(-5., 5.); printf("\n");
  */
  TString runPeriod = "ABCD" ;
  
  TF1** sf_light_mean = new TF1*[mistag_eta_binning_medium_NUM-1];
  for (Int_t i=1; i<mistag_eta_binning_medium_NUM; i++) {
    sf_light_mean[i-1] = NULL;
    sf_light_mean[i-1] = GetSFLight("mean","CSV","M", mistag_eta_binning_medium[i-1], mistag_eta_binning_medium[i], runPeriod);
  }
  
  TF1** sf_light_mean_WPloose = new TF1*[mistag_eta_binning_loose_NUM-1];
  for (Int_t i=1; i<mistag_eta_binning_loose_NUM; i++) {
    sf_light_mean_WPloose[i-1] = NULL;
    sf_light_mean_WPloose[i-1] = GetSFLight("mean","CSV","L",mistag_eta_binning_loose[i-1], mistag_eta_binning_loose[i], runPeriod);
  }
  
  
  clock_t start = clock();
  
  cout << "***************************************************************" << endl;
  cout << " Beginning of the program for the V+Jets background estimation " << endl;
  cout << "***************************************************************" << endl;
  
  bool Wbb_from_WJets = true;
  bool doPE = false;
  bool isOnData = false;
  bool reloadHisto = false;
    //  bool isWbbFromWbbMeas = false; //else, take the measure cross-section ; not using the Flavour History Path
  bool semiMuon = true;
  bool semiElectron = true;

  Bool_t WnJets = kTRUE ;
  Bool_t DYnJets = kTRUE ;
  Bool_t TTJets_inclusive = kTRUE ;

  Bool_t estimateQCDApart = kFALSE ;
  Bool_t addRareProcesses = kFALSE ;
  Bool_t isNewNaming = kTRUE ;

  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////Configuration ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  const UInt_t NbOfPE = 10000;
  const UInt_t NbOfJetBins = (isNewNaming ? 4 : 3);
  const UInt_t MinNbOfJets = 4;
  const bool isLastBinExclusive = false;
  const UInt_t NbOfWP = 1;
  const UInt_t muonMask = (1<<0);
  const UInt_t electronMask = (1<<1);
  const UInt_t NbOfChannels = muonMask + electronMask;
  // To read input files, btag-binned plots
  const UInt_t NbOfBtagBins = (isNewNaming ? 5 : 1 );
  const UInt_t MinNbOfBtagJets = 0;
  const bool isLastBtagBinExclusive = false;

    //  std::string leg3Dir = "/afs/cern.ch/work/a/aocampor/public/RunABCD/";
    //  std::string leg3Dir = "/afs/cern.ch/work/a/aocampor/public/RunABCD6930/";
  //      std::string leg3Dir = "$HOME/temp/CONFIGS/" ;
      //    std::string leg3Dir = "$HOME/Leg3ProductionAndVJets_8TeV/Leg3_user_aocampor_public/" ;

  std::string leg3Dir = "$HOME/RunABCD_AllSel_DD/" ;

  //    std::string leg3Dir = "$HOME/Leg3ProductionAndVJets_8TeV/RunABCD_AllSel_v3/" ;
  //  std::string leg3Dir = "./RunABCD_AllSel_v2/" ;
//  std::string leg3Dir = "./RunABCD6930/";
    //    std::string leg3Dir = "/afs/cern.ch/work/a/aocampor/public/RunABCD6930/" ;
  std::string qcdShapeFile = "../Downloads/histos_separate(2).root" ;
  std::string pathCMSSM = "$HOME/Leg3ProductionAndVJets_8TeV/" ;
  std::string WbbShapeRootfile = ".root"; //External histo for Wbb from flavour history path
  std::string btagJetMult_prefix = "BtaggedJet_Multiplicity";
  std::string btagJetMultbb_prefix = "BtaggedJet_Multiplicity_bb";
  std::string btagJetMult_suffix = "" ;
  //  std::string btagJetMult_suffix = "_met_cut_60" ;

    //Output ROOT file
  std::string postfix = "_VJetsBckgdEst";
  std::string channelpostfix = "_";
  if (semiMuon) {
    channelpostfix = channelpostfix + "SemiMuon"; 
  }
  if (semiElectron) {
    channelpostfix = channelpostfix + "SemiElectron"; 
  }
  std::string comment = "_Constants_euds";
  if (isOnData) {
    comment = comment + "_Data" ;
  }
  std::string rootFileName ("StopSearches"+postfix+channelpostfix+comment+".root");
  std::string histoFileName = rootFileName ;
  if (reloadHisto) {
    rootFileName = "_" + rootFileName;
  }
  
    //  std::string albertoFileName(".root");
  Bool_t forceE3bq = kFALSE;
  Float_t fracForceE3bq = 0.01;
  
  
  /*
   LumiReWeighting LumiWeights;
   LumiWeights = reweight::LumiReWeighting("../pileup_MC_Fall11.root", "../pileup_2011Data_UpToRun180252.root", "pileup", "pileup");
   */
  
  /** TO BE ADAPTED **/
  
  
  
  
    //  TChain *vbjets = new TChain("AnaTree");
    //  vbjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/XXX_skim.root");
    //  TEntryList *entrylist_vbjets = 0;
    //  vector<double>   *vjtDiscri_pf_vbjets = 0;
  
    //  TBranch        *b_vjtDiscri_pf = 0;
  
    // Retrieve int. lumi from InputTree (?)
  double *IntLumi = new double[NbOfChannels];
  for (UInt_t c=0; c<NbOfChannels; c++) {
    IntLumi[c] = 0.;
  }
  for (UInt_t c=0; c<NbOfChannels; c++) {
    if (((c+1) & muonMask)!=0) {
      if (runPeriod.Contains("A")) {
        IntLumi[c] += 0.890608 * 1000. ;
      }
      if (runPeriod.Contains("B")) {
        IntLumi[c] += 4.429 * 1000. ;
      }
      if (runPeriod.Contains("C")) {
        IntLumi[c] += 6.892003 * 1000. ;
      }
      if (runPeriod.Contains("D")) {
        IntLumi[c] += 6.752 * 1000. ;
      }
    }
    if (((c+1) & electronMask)!=0) {
      if (runPeriod.Contains("A")) {
        IntLumi[c] += 0.890608 * 1000. ;
      }
      if (runPeriod.Contains("B")) {
        IntLumi[c] += 4.429 * 1000. ;
      }
      if (runPeriod.Contains("C")) {
        IntLumi[c] += 6.892003 * 1000. ;
      }
      if (runPeriod.Contains("D")) {
        IntLumi[c] += 6.752 * 1000. ;
      }
    }
  }
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    printf("IntLumi[%d] = %lf\n", c, IntLumi[c]);
  }
  
  
  UInt_t NbOfCategories = 5;
  std::vector<std::vector<TFile*> > files(NbOfCategories+3, std::vector<TFile*>());
  std::vector<std::vector<std::string> > samples_names(NbOfCategories+3, std::vector<std::string>());
  std::vector<std::vector<std::string> > samples_labelBQ(NbOfCategories+3, std::vector<std::string>());
  std::vector<std::vector<Double_t> > samples_weights(NbOfCategories+3, std::vector<Double_t>());
  
  std::vector<std::string> dataNames;
  std::vector<TFile*> datafiles;
  std::vector<Double_t> weight_data;
  std::vector<Double_t> weight_data_est;
  
    //  std::string CChannel[2] = {"_mu","_el"};
  Double_t Wn_SF = 37509. / 30400. ; // To be checked on 8 TeV data !!!!!!!!
  Double_t DYn_SF = 3503.71 / 2950.0 ; // To be checked on 8 TeV data !!!!!!!!

  Double_t entries = -1.;
  
  TFile *fBtagEff = new TFile((pathCMSSM+"Leg3_8TeV_DD/test/BTagEff_MCsamples.root").c_str(), "read");
  
  /*
   for (int i=0; (1<<i)<=channelConf+1; i++) {
   Double_t entries = -1.;
   Int_t channel=(1<<i)-1 ;
   if (((channel+1)&(channelConf+1))==0) {
   continue;
   }
   printf("channel : %d\n", channel) ;
   */
  
  /*
   DY1Jets
   DY2Jets
   DY3Jets
   DY4Jets
   DYJetsInc
   W1Jets
   W2Jets
   W3Jets
   W4Jets
   WJetsInc
   TTJetsSemiLep
   TTJetsHad
   TTJetsFullLep
   TTJetsInc
   TopS
   TopT
   TopTW
   TbarS
   TbarT
   TbarTW
   QCDEM20
   QCDEM30
   QCDEM170
   QCDEM80
   QCDEM250
   QCDEM350
   QCDMu15
   QCDMu20
   QCDMu30
   QCDMu50
   QCDMu80
   QCDMu120
   QCDMu170
   QCDMu300
   QCDMu470
   QCDMu600
   QCDMu800
   QCDMu1000
   QCDMuInc
   TTWJets
   TTWWJets
   TTZJets
   WW
   WZ
   ZZ
   WWWJets
   WWZJets
   WZZJets
   ZZZJets
   S200N50
   S250N100
   S300N150
   S300N50
   S350N100
   S350N200
   S400N150
   S400N50
   S450N100
   S450N200
   S500N150
   S500N50
   S550N100
   S550N200
   S600N150
   S650N200
   */  
  
    //TTbar
  {
    if (TTJets_inclusive) {
      samples_names[0].push_back("TTJets"/*+CChannel[channel]*/);
      samples_labelBQ[0].push_back("TTJetsInc"/*+CChannel[channel]*/);
      files[0].push_back(TFile::Open((leg3Dir+(isNewNaming==kTRUE?"TTJetsInc":"TTJets")/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[0].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[0].rbegin())->GetName(), ((TH1D*) (* files[0].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[0].push_back((entries==0.? 0. : 234. /entries )); //TTbar
    } else {
       samples_names[0].push_back("TTJets"/*+CChannel[channel]*/);
      samples_labelBQ[0].push_back("TTJetsSemiLep"/*+CChannel[channel]*/);
      files[0].push_back(TFile::Open((leg3Dir+"TTJetsSemiLep"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[0].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[0].rbegin())->GetName(), ((TH1D*) (* files[0].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[0].push_back((entries==0.? 0. : 234. /entries )); //TTbar

       samples_names[0].push_back("TTJets"/*+CChannel[channel]*/);
      samples_labelBQ[0].push_back("TTJetsFullLep"/*+CChannel[channel]*/);
      files[0].push_back(TFile::Open((leg3Dir+"TTJetsFullLep"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[0].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[0].rbegin())->GetName(), ((TH1D*) (* files[0].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[0].push_back((entries==0.? 0. : 234. /entries )); //TTbar

       samples_names[0].push_back("TTJets"/*+CChannel[channel]*/);
      samples_labelBQ[0].push_back("TTJetsHad"/*+CChannel[channel]*/);
      files[0].push_back(TFile::Open((leg3Dir+"TTJetsHad"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[0].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[0].rbegin())->GetName(), ((TH1D*) (* files[0].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[0].push_back((entries==0.? 0. : 234. /entries )); //TTbar
   }
  }
  
    //W+jets
  {
    if (WnJets) {
      
      // samples_names[1].push_back("W1Jets"/*+CChannel[channel]*/);
      // samples_labelBQ[1].push_back("W1Jets");
      // files[1].push_back(TFile::Open((leg3Dir+"W1Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      // entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      // printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      // samples_weights[1].push_back((entries==0.? 0. : Wn_SF*5400. /entries )); //W1jets
      
      samples_names[1].push_back("W2Jets"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("W2Jets");
      files[1].push_back(TFile::Open((leg3Dir+"W2Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : Wn_SF*1750. /entries )); //W2jets
      
      samples_names[1].push_back("W3Jets"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("W3Jets");
      files[1].push_back(TFile::Open((leg3Dir+"W3Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : Wn_SF*519. /entries )); //W3jets
      
      samples_names[1].push_back("W4Jets"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("W4Jets");
      files[1].push_back(TFile::Open((leg3Dir+"W4Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : Wn_SF*214. /entries )); //W4jets
    } else {
      samples_names[1].push_back("WJets"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("WJetsInc");
      files[1].push_back(TFile::Open((leg3Dir+"WJets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : 37509.0 /entries )); //Wjets
    }
   
    //Z+jets
    if (DYnJets) {
      // samples_names[1].push_back("DY1Jets"/*+CChannel[channel]*/);
      // samples_labelBQ[1].push_back("DY1Jets");
      // files[1].push_back(TFile::Open((leg3Dir+"DY1Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      // entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      // printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      // samples_weights[1].push_back((entries==0.? 0. : DYn_SF*561. /entries )); //DY1jets

      samples_names[1].push_back("DY2Jets"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("DY2Jets");
      files[1].push_back(TFile::Open((leg3Dir+"DY2Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : DYn_SF*181. /entries )); //W2jets
      
      samples_names[1].push_back("DY3Jets"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("DY3Jets");
      files[1].push_back(TFile::Open((leg3Dir+"DY3Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : DYn_SF*51.1 /entries )); //W3jets
      
      samples_names[1].push_back("DY4Jets"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("DY4Jets");
      files[1].push_back(TFile::Open((leg3Dir+"DY4Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : DYn_SF*23.04 /entries )); //W4jets 
    } else  {
      samples_names[1].push_back("DYJets"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("DYJetsInc");
      files[1].push_back(TFile::Open((leg3Dir+"DYJets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : 3503.71 /entries )); //Zjets
    }  
 
    //QCD
    if (estimateQCDApart == kFALSE) {
      //      if (channel==0) {
      samples_names[1].push_back("QCD_mu20"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("QCDMuInc");
      files[1].push_back(TFile::Open((leg3Dir+"QCDMu"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : 3.64e8 * 3.7e-4 /entries )); //QCD_mu20
      
      //      } else {
      samples_names[1].push_back("QCD_el_20_30"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("QCDEM20");
      files[1].push_back(TFile::Open((leg3Dir+"QCDEM20"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : 2.886E8 * 0.0101 /entries )); //QCD_el_20_30
      
      samples_names[1].push_back("QCD_el_30_80"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("QCDEM30");
      files[1].push_back(TFile::Open((leg3Dir+"QCDEM30"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : 7.433E7 * 0.0621 /entries )); //QCD_el_30_80
      
      samples_names[1].push_back("QCD_el_80_170"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("QCDEM80");
      files[1].push_back(TFile::Open((leg3Dir+"QCDEM80"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : 1191000.0 * 0.1539 /entries )); //QCD_el_80_170
      
      samples_names[1].push_back("QCD_el_170_250"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("QCDEM170");
      files[1].push_back(TFile::Open((leg3Dir+"QCDEM170"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : 30990.0 * 0.148 /entries )); //QCD_el_170_250
      
      samples_names[1].push_back("QCD_el_250_350"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("QCDEM250");
      files[1].push_back(TFile::Open((leg3Dir+"QCDEM250"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : 4250.0 * 0.131 /entries )); //QCD_el_250_350
      
      samples_names[1].push_back("QCD_el_350_"/*+CChannel[channel]*/);
      samples_labelBQ[1].push_back("QCDEM350");
      files[1].push_back(TFile::Open((leg3Dir+"QCDEM350"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[1].rbegin())->GetName(), ((TH1D*) (* files[1].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[1].push_back((entries==0.? 0. : 810.0 * 0.11 /entries )); //QCD_el_350_
      
      //      }
    }
  }
  
  //V-bb component
  {
    //W+jets
    if (WnJets) {
      // samples_names[2].push_back("W1Jets"/*+CChannel[channel]*/);
      // samples_labelBQ[2].push_back("W1Jets");
      // files[2].push_back(TFile::Open((leg3Dir+"W1Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      // entries = ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2);
      // printf("Entries (%s) : %lf\n", (* files[2].rbegin())->GetName(), ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2));
      // samples_weights[2].push_back((entries==0.? 0. : Wn_SF*5400. /entries )); //W1jets

      samples_names[2].push_back("W2Jets"/*+CChannel[channel]*/);
      samples_labelBQ[2].push_back("W2Jets");
      files[2].push_back(TFile::Open((leg3Dir+"W2Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[2].rbegin())->GetName(), ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[2].push_back((entries==0.? 0. : Wn_SF*1750. /entries )); //W2jets
      
      samples_names[2].push_back("W3Jets"/*+CChannel[channel]*/);
      samples_labelBQ[2].push_back("W3Jets");
      files[2].push_back(TFile::Open((leg3Dir+"W3Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[2].rbegin())->GetName(), ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[2].push_back((entries==0.? 0. : Wn_SF*519. /entries )); //W3jets
      
      samples_names[2].push_back("W4Jets"/*+CChannel[channel]*/);
      samples_labelBQ[2].push_back("W4Jets");
      files[2].push_back(TFile::Open((leg3Dir+"W4Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[2].rbegin())->GetName(), ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[2].push_back((entries==0.? 0. : Wn_SF*214. /entries )); //W4jets
    } else {
      samples_names[2].push_back("WJets"/*+CChannel[channel]*/);
      samples_labelBQ[2].push_back("WJetsInc");
      files[2].push_back(TFile::Open((leg3Dir+"WJets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[2].rbegin())->GetName(), ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[2].push_back((entries==0.? 0. : 37509.0 /entries )); //Wjets
      
    }
    //Z+jets
    if (DYnJets) {
      // samples_names[2].push_back("DY1Jets"/*+CChannel[channel]*/);
      // samples_labelBQ[2].push_back("DY1Jets");
      // files[2].push_back(TFile::Open((leg3Dir+"DY1Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      // entries = ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2);
      // printf("Entries (%s) : %lf\n", (* files[2].rbegin())->GetName(), ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2));
      // samples_weights[2].push_back((entries==0.? 0. : DYn_SF*561. /entries )); //DY1jets
      
      samples_names[2].push_back("DY2Jets"/*+CChannel[channel]*/);
      samples_labelBQ[2].push_back("DY2Jets");
      files[2].push_back(TFile::Open((leg3Dir+"DY2Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[2].rbegin())->GetName(), ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[2].push_back((entries==0.? 0. : DYn_SF*181. /entries )); //W2jets
      
      samples_names[2].push_back("DY3Jets"/*+CChannel[channel]*/);
      samples_labelBQ[2].push_back("DY3Jets");
      files[2].push_back(TFile::Open((leg3Dir+"DY3Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[2].rbegin())->GetName(), ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[2].push_back((entries==0.? 0. : DYn_SF*51.1 /entries )); //W3jets
      
      samples_names[2].push_back("DY4Jets"/*+CChannel[channel]*/);
      samples_labelBQ[2].push_back("DY4Jets");
      files[2].push_back(TFile::Open((leg3Dir+"DY4Jets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[2].rbegin())->GetName(), ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[2].push_back((entries==0.? 0. : DYn_SF*23.04 /entries )); //W4jets 
    } else  {
      samples_names[2].push_back("DYJets"/*+CChannel[channel]*/);
      samples_labelBQ[2].push_back("DYJetsInc");
      files[2].push_back(TFile::Open((leg3Dir+"DYJets"/*+CChannel[channel]*/+".root").c_str(), "READ"));
      entries = ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2);
      printf("Entries (%s) : %lf\n", (* files[2].rbegin())->GetName(), ((TH1D*) (* files[2].rbegin())->Get("Entries"))->GetBinContent(2));
      samples_weights[2].push_back((entries==0.? 0. : 3503.71 /entries )); //Zjets
    }
  }
  

    //Single Top
  {
    samples_names[3].push_back("TopS"/*+CChannel[channel]*/);
    samples_labelBQ[3].push_back("TopS");
    files[3].push_back(TFile::Open((leg3Dir+"TopS"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[3].rbegin())->GetName(), ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[3].push_back((entries==0.? 0. : 3.79 /entries )); //
    
    samples_names[3].push_back("TbarS"/*+CChannel[channel]*/);
    samples_labelBQ[3].push_back("TbarS");
    files[3].push_back(TFile::Open((leg3Dir+"TbarS"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[3].rbegin())->GetName(), ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[3].push_back((entries==0.? 0. : 1.76 /entries )); //
    
    samples_names[3].push_back("TopT"/*+CChannel[channel]*/);
    samples_labelBQ[3].push_back("TopT");
    files[3].push_back(TFile::Open((leg3Dir+"TopT"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[3].rbegin())->GetName(), ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[3].push_back((entries==0.? 0. : 56.4 /entries )); //
    
    samples_names[3].push_back("TbarT"/*+CChannel[channel]*/);
    samples_labelBQ[3].push_back("TbarT");
    files[3].push_back(TFile::Open((leg3Dir+"TbarT"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[3].rbegin())->GetName(), ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[3].push_back((entries==0.? 0. : 30.7 /entries )); //
    
    samples_names[3].push_back("TopTW"/*+CChannel[channel]*/);
    samples_labelBQ[3].push_back("TopTW");
    files[3].push_back(TFile::Open((leg3Dir+"TopTW"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[3].rbegin())->GetName(), ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[3].push_back((entries==0.? 0. : 11.1 /entries )); //
    
    samples_names[3].push_back("TbarTW"/*+CChannel[channel]*/);
    samples_labelBQ[3].push_back("TbarTW");
    files[3].push_back(TFile::Open((leg3Dir+"TbarTW"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[3].rbegin())->GetName(), ((TH1D*) (* files[3].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[3].push_back((entries==0.? 0. : 11.1 /entries )); //
  }
  
    //VV
  {
    samples_names[4].push_back("WW"/*+CChannel[channel]*/);
    samples_labelBQ[4].push_back("WW");
    files[4].push_back(TFile::Open((leg3Dir+"WW"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[4].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[4].rbegin())->GetName(), ((TH1D*) (* files[4].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[4].push_back((entries==0.? 0. : 56.7532 /entries )); //
    
    samples_names[4].push_back("WZ"/*+CChannel[channel]*/);
    samples_labelBQ[4].push_back("WZ");
    files[4].push_back(TFile::Open((leg3Dir+"WZ"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[4].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[4].rbegin())->GetName(), ((TH1D*) (* files[4].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[4].push_back((entries==0.? 0. : 33.85 /entries )); //
    
    samples_names[4].push_back("ZZ"/*+CChannel[channel]*/);
    samples_labelBQ[4].push_back("ZZ");
    files[4].push_back(TFile::Open((leg3Dir+"ZZ"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[4].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[4].rbegin())->GetName(), ((TH1D*) (* files[4].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[4].push_back((entries==0.? 0. : 8.297 /entries )); //
  }
  
      //Extra Category
  if (estimateQCDApart == kTRUE) {
      //      if (channel==0) {
    samples_names[NbOfCategories].push_back("QCD_mu20"/*+CChannel[channel]*/);
    samples_labelBQ[NbOfCategories].push_back("QCDMuInc");
    files[NbOfCategories].push_back(TFile::Open((leg3Dir+"QCDMu"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[NbOfCategories].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[NbOfCategories].push_back((entries==0.? 0. : 3.64e8 * 3.7e-4 /entries )); //QCD_mu20
    
      //      } else {
    samples_names[NbOfCategories].push_back("QCD_el_20_30"/*+CChannel[channel]*/);
    samples_labelBQ[NbOfCategories].push_back("QCDEM20");
    files[NbOfCategories].push_back(TFile::Open((leg3Dir+"QCDEM20"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[NbOfCategories].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[NbOfCategories].push_back((entries==0.? 0. : 2.886E8 * 0.0101 /entries )); //QCD_el_20_30
    
    samples_names[NbOfCategories].push_back("QCD_el_30_80"/*+CChannel[channel]*/);
    samples_labelBQ[NbOfCategories].push_back("QCDEM30");
    files[NbOfCategories].push_back(TFile::Open((leg3Dir+"QCDEM30"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[NbOfCategories].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[NbOfCategories].push_back((entries==0.? 0. : 7.433E7 * 0.0621 /entries )); //QCD_el_30_80
    
    samples_names[NbOfCategories].push_back("QCD_el_80_170"/*+CChannel[channel]*/);
    samples_labelBQ[NbOfCategories].push_back("QCDEM80");
    files[NbOfCategories].push_back(TFile::Open((leg3Dir+"QCDEM80"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[NbOfCategories].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[NbOfCategories].push_back((entries==0.? 0. : 1191000.0 * 0.1539 /entries )); //QCD_el_80_170
    
    samples_names[NbOfCategories].push_back("QCD_el_170_250"/*+CChannel[channel]*/);
    samples_labelBQ[NbOfCategories].push_back("QCDEM170");
    files[NbOfCategories].push_back(TFile::Open((leg3Dir+"QCDEM170"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[NbOfCategories].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[NbOfCategories].push_back((entries==0.? 0. : 30990.0 * 0.148 /entries )); //QCD_el_170_250
    
    samples_names[NbOfCategories].push_back("QCD_el_250_350"/*+CChannel[channel]*/);
    samples_labelBQ[NbOfCategories].push_back("QCDEM250");
    files[NbOfCategories].push_back(TFile::Open((leg3Dir+"QCDEM250"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[NbOfCategories].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[NbOfCategories].push_back((entries==0.? 0. : 4250.0 * 0.131 /entries )); //QCD_el_250_350
    
    samples_names[NbOfCategories].push_back("QCD_el_350_"/*+CChannel[channel]*/);
    samples_labelBQ[NbOfCategories].push_back("QCDEM350");
    files[NbOfCategories].push_back(TFile::Open((leg3Dir+"QCDEM350"/*+CChannel[channel]*/+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[NbOfCategories].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[NbOfCategories].push_back((entries==0.? 0. : 810.0 * 0.11 /entries )); //QCD_el_350_
    
      //      }
  }

  // Rare processes
  {
    if (addRareProcesses == kTRUE) 
      { // ZZZJets , WZZJets , WWZJets , WWWJets , TTZJets , TTWWJets , TTWJets
        //      if (channel==0) {
        samples_names[NbOfCategories+1].push_back("ZZZJets");
        samples_labelBQ[NbOfCategories+1].push_back("ZZZJets");
        files[NbOfCategories+1].push_back(TFile::Open((leg3Dir+"ZZZJets"+".root").c_str(), "READ"));
        entries = ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2);
        printf("Entries (%s) : %lf\n", (* files[NbOfCategories+1].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2));
        samples_weights[NbOfCategories+1].push_back((entries==0.? 0. : 0.004587 /entries )); //ZZZJets
        
        samples_names[NbOfCategories+1].push_back("WZZJets");
        samples_labelBQ[NbOfCategories+1].push_back("WZZJets");
        files[NbOfCategories+1].push_back(TFile::Open((leg3Dir+"WZZJets"+".root").c_str(), "READ"));
        entries = ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2);
        printf("Entries (%s) : %lf\n", (* files[NbOfCategories+1].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2));
        samples_weights[NbOfCategories+1].push_back((entries==0.? 0. : 0.01922 /entries )); //WZZJets
        
        samples_names[NbOfCategories+1].push_back("WWZJets");
        samples_labelBQ[NbOfCategories+1].push_back("WWZJets");
        files[NbOfCategories+1].push_back(TFile::Open((leg3Dir+"WWZJets"+".root").c_str(), "READ"));
        entries = ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2);
        printf("Entries (%s) : %lf\n", (* files[NbOfCategories+1].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2));
        samples_weights[NbOfCategories+1].push_back((entries==0.? 0. : 0.0633 /entries )); //WWZJets
        
        samples_names[NbOfCategories+1].push_back("WWWJets");
        samples_labelBQ[NbOfCategories+1].push_back("WWWJets");
        files[NbOfCategories+1].push_back(TFile::Open((leg3Dir+"WWWJets"+".root").c_str(), "READ"));
        entries = ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2);
        printf("Entries (%s) : %lf\n", (* files[NbOfCategories+1].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2));
        samples_weights[NbOfCategories+1].push_back((entries==0.? 0. : 0.08217 /entries )); //WWWJets
        
        samples_names[NbOfCategories+1].push_back("TTZJets");
        samples_labelBQ[NbOfCategories+1].push_back("TTZJets");
        files[NbOfCategories+1].push_back(TFile::Open((leg3Dir+"TTZJets"+".root").c_str(), "READ"));
        entries = ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2);
        printf("Entries (%s) : %lf\n", (* files[NbOfCategories+1].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2));
        samples_weights[NbOfCategories+1].push_back((entries==0.? 0. : 0.172 /entries )); //TTZJets
        
        samples_names[NbOfCategories+1].push_back("TTWWJets");
        samples_labelBQ[NbOfCategories+1].push_back("TTWWJets");
        files[NbOfCategories+1].push_back(TFile::Open((leg3Dir+"TTWWJets"+".root").c_str(), "READ"));
        entries = ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2);
        printf("Entries (%s) : %lf\n", (* files[NbOfCategories+1].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2));
        samples_weights[NbOfCategories+1].push_back((entries==0.? 0. : 0.002037 /entries )); //TTWWJets
        
        samples_names[NbOfCategories+1].push_back("TTWJets");
        samples_labelBQ[NbOfCategories+1].push_back("TTWJets");
        files[NbOfCategories+1].push_back(TFile::Open((leg3Dir+"TTWJets"+".root").c_str(), "READ"));
        entries = ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2);
        printf("Entries (%s) : %lf\n", (* files[NbOfCategories+1].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories+1].rbegin())->Get("Entries"))->GetBinContent(2));
        samples_weights[NbOfCategories+1].push_back((entries==0.? 0. : 0.2149 /entries )); //TTWJets
      }
  }

  if (signalFile.compare("")!=0) {
      //      if (channel==0) {
    samples_names[NbOfCategories+2].push_back(signalFile);
    samples_labelBQ[NbOfCategories+2].push_back(signalFile);
    files[NbOfCategories+2].push_back(TFile::Open((leg3Dir+signalFile+".root").c_str(), "READ"));
    entries = ((TH1D*) (* files[NbOfCategories+2].rbegin())->Get("Entries"))->GetBinContent(2);
    printf("Entries (%s) : %lf\n", (* files[NbOfCategories+2].rbegin())->GetName(), ((TH1D*) (* files[NbOfCategories+2].rbegin())->Get("Entries"))->GetBinContent(2));
    samples_weights[NbOfCategories+2].push_back((entries==0.? 0. : signal_xs /entries )); //QCD_mu20
    
  }

    // Data

  {
    
    dataNames.push_back("Single Mu");
    datafiles.push_back(TFile::Open((leg3Dir+"Muons.root").c_str(), "READ"));
    weight_data.push_back( 1 );
    
    dataNames.push_back("Single El");
    datafiles.push_back(TFile::Open((leg3Dir+"Electrons.root").c_str(), "READ"));
    weight_data.push_back( 1 );    

  }
  /** 
  {
   //      if (((channel+1)&1)!=0) {
   if (runPeriod.Contains("A")) {
   dataNames.push_back("MuHad_13Jul_RunA");
   datafiles.push_back(TFile::Open((leg3Dir+(isNewNaming ?"MuHad13Jul.root" : "MuHad_13Jul_RunA.root")).c_str(), "READ"));
   weight_data.push_back( 1 );
   
   //   dataNames.push_back("MuHad_06Aug_RunA");
   //   datafiles.push_back(TFile::Open((leg3Dir+(isNewNaming ?"MuHad06Aug.root":"MuHad_06Aug_RunA.root")).c_str(), "READ"));
   //   weight_data.push_back( 1 );
   }
   
   if (runPeriod.Contains("B")) {
   dataNames.push_back("Single Mu RunB");
   datafiles.push_back(TFile::Open((leg3Dir+(isNewNaming ?"SingleMuB.root":"SingleMu_13Jul_RunB.root")).c_str(), "READ"));
   weight_data.push_back( 1 );
   }
   
   if (runPeriod.Contains("C")) {
     // dataNames.push_back("Single Mu RunC_v1");
     //  datafiles.push_back(TFile::Open((leg3Dir+"SingleMu_PRv1_RunC.root").c_str(), "READ"));
     // weight_data.push_back( 1 );
     dataNames.push_back("Single Mu RunC_v2");
     datafiles.push_back(TFile::Open((leg3Dir+(isNewNaming ?"SingleMuC_v2.root":"SingleMu_PRv2_RunC.root")).c_str(), "READ"));
     weight_data.push_back( 1 );
     
     if(isNewNaming) {
       dataNames.push_back("Single Mu RunC_11Dec");
       datafiles.push_back(TFile::Open((leg3Dir+"SingleMuC_11Dec.root").c_str(), "READ"));
       weight_data.push_back( 1 );
       
       dataNames.push_back("Single Mu RunC_24Aug");
       datafiles.push_back(TFile::Open((leg3Dir+"SingleMuC_24Aug.root").c_str(), "READ"));
       weight_data.push_back( 1 );
     }
   }
   
   if (runPeriod.Contains("D")) {
   dataNames.push_back("Single Mu RunD");
   datafiles.push_back(TFile::Open((leg3Dir+(isNewNaming ?"SingleMuD.root":"SingleMu_RunD.root")).c_str(), "READ"));
   weight_data.push_back( 1 ); //
   }
   //      }
   //      if (((channel+1)&2)!=0) {
   if (runPeriod.Contains("A")) {
   dataNames.push_back("ElHad_13Jul_RunA");
   datafiles.push_back(TFile::Open((leg3Dir+(isNewNaming ?"EleHad13Jul.root":"EleHad_13Jul_RunA.root")).c_str(), "READ"));
   weight_data.push_back( 1 );
   
   dataNames.push_back("ElHad_06Aug_RunA");
   datafiles.push_back(TFile::Open((leg3Dir+(isNewNaming ?"EleHad06Aug.root":"EleHad_06Aug_RunA.root")).c_str(), "READ"));
   weight_data.push_back( 1 );
   }
   
   if (runPeriod.Contains("B")) {
   dataNames.push_back("Single El RunB");
   datafiles.push_back(TFile::Open((leg3Dir+(isNewNaming ?"SingleEleB.root":"SingleEl_13Jul_RunB.root")).c_str(), "READ"));
   weight_data.push_back( 1 );
   }
   
   if (runPeriod.Contains("C")) {
     // dataNames.push_back("Single El RunC_v1");
     // datafiles.push_back(TFile::Open((leg3Dir+"SingleEl_PRv1_RunC.root").c_str(), "READ"));
     // weight_data.push_back( 1 );
     dataNames.push_back("Single El RunC_v2");
     datafiles.push_back(TFile::Open((leg3Dir+(isNewNaming ?"SingleEleC_v2.root":"SingleEl_PRv2_RunC.root")).c_str(), "READ"));
     weight_data.push_back( 1 );
     
     if(isNewNaming) {
       dataNames.push_back("Single El RunC_11Dec");
       datafiles.push_back(TFile::Open((leg3Dir+"SingleEleC_11Dec.root").c_str(), "READ"));
       weight_data.push_back( 1 );
       
       dataNames.push_back("Single El RunC_24Aug");
       datafiles.push_back(TFile::Open((leg3Dir+"SingleEleC_24Aug.root").c_str(), "READ"));
       weight_data.push_back( 1 );
     }
   }
   
   if (runPeriod.Contains("D")) {
   dataNames.push_back("Single El RunD");
   datafiles.push_back(TFile::Open((leg3Dir+(isNewNaming ? "SingleEleD.root":"SingleEl_RunD.root")).c_str(), "READ"));
   weight_data.push_back( 1 ); //
   }
   //      }
   }
   */
  /* } */
  entries = -1. ;
  
  
  
  
  
  
  
  
  printf("Parameters jobs definition\n");
  
    // E v e n t   s e l e c t i o n   c u t s
    // ---------------------------------------
  
  
  
  /** B-tagging working points ; https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP#B_tagging_Operating_Points_for_3
   TrackCountingHighEff      TCHEL 1.7
   TrackCountingHighEff 	    TCHEM 3.3
   TrackCountingHighEff 	    TCHET 10.2 (not supported)
   CombinedSecondaryVertex   CSVL  0.244
   CombinedSecondaryVertex   CSVM  0.679
   CombinedSecondaryVertex   CSVT  0.898 
   **/
  float *btagCuts = new float[NbOfWP];
  btagCuts[0] = 1.7;
  if (NbOfWP>1)
    btagCuts[1] = 3.3;
  if (NbOfWP>2)
    btagCuts[2] = 6.;
  if (NbOfWP>3)
    btagCuts[3] = 4.7;
  
  char** channelLeg3Suffix = new char*[NbOfChannels];
  for (UInt_t i=0; i<NbOfChannels; i++) {
    channelLeg3Suffix[i] = new char[50];
    sprintf(channelLeg3Suffix[i], "_%s", (std::string()
                                          +(((i+1)&(1<<0))?"Mu":"")
                                          +(((i+1)&(1<<1))?"El":"")).c_str());
  }    
  
  char** channelSuffix = new char*[NbOfChannels];
  for (UInt_t i=0; i<NbOfChannels; i++) {
    channelSuffix[i] = new char[50];
    sprintf(channelSuffix[i], "_%s", (std::string()
                                      +(((i+1)&(1<<0))?"SemiMuon":"")
                                      +(((i+1)&(1<<1))?"SemiElectron":"")).c_str());
  }    
  
  
  char** jetSuffix = new char*[NbOfJetBins];
  char** jtSuf = new char*[NbOfJetBins];
  for (UInt_t i=0; i<NbOfJetBins; i++) {
    jetSuffix[i] = new char[50];
    jtSuf[i] = new char[5];
    if (i==NbOfJetBins-1 && !isLastBinExclusive) {
      sprintf(jetSuffix[i], "_%djInc", MinNbOfJets+i);
      sprintf(jtSuf[i], "_%dj", MinNbOfJets+i);
    } else {
      sprintf(jetSuffix[i], "_%djExc", MinNbOfJets+i);
      sprintf(jtSuf[i], "_%dj", MinNbOfJets+i);
    }
  }    
  
  char** wpSuffix = new char*[NbOfWP];
  char** wpSuf = new char*[NbOfWP];
  if (NbOfWP==1) {
    for (UInt_t i=0; i<NbOfWP; i++) {
      wpSuffix[0] = new char[50];
      wpSuf[0] = new char[5];
      sprintf(wpSuffix[i], "");
      sprintf(wpSuf[i], "_wp%d", 0);
    }    
  } else {
    for (UInt_t i=0; i<NbOfWP; i++) {
      wpSuffix[i] = new char[50];
      wpSuf[i] = new char[5];
      sprintf(wpSuffix[i], "_WP%d", i);
      sprintf(wpSuf[i], "_wp%d", i);
    }
  }
  
  char**nbtagSuffix = new char*[NbOfBtagBins];
  for (UInt_t i=0; i<NbOfBtagBins; i++) {
    nbtagSuffix[i] = new char[50];
    if (isNewNaming == kTRUE) {
      if (i==NbOfBtagBins-1 && !isLastBtagBinExclusive) {
        sprintf(nbtagSuffix[i], "_%dbInc", MinNbOfBtagJets+i);
      } else {
        sprintf(nbtagSuffix[i], "_%dbExc", MinNbOfBtagJets+i);
      }
    } else {
      sprintf(nbtagSuffix[i], "");
    }
  }  

    // Muon channel
    //  TCut muon = "(muPt[0]>0.)&&(TMath::Abs(muEta[0])<=25.)";
    //  TCut muon = "(muPt[0]>20.)&&(TMath::Abs(muEta[0])<=2.1)";
    // Electron channel
    //TCut electron = "(elPt[0]>0.)&&(elSCenergy[0]>=17.0)";
    // Lepton selection
    //  TCut lveto = "NmuIdIso>=0";
    //  TCut lveto = "NmuIdIso==1";
  
    // Channel selection
    //  TCut singlep = muon+lveto;
  
    // Jet selection
    //  TCut pfnjet[NbOfJetBins];
    //  pfnjet[0] = "Njt_pf==3";
    //  pfnjet[1] = "Njt_pf==4";
    //  pfnjet[2] = "Njt_pf==5";
    //  pfnjet[3] = "Njt_pf>=6";
    //TCut pf4jet = "Njt_pf>3";
    //TCut btag   = "Njtbtag_pf>0";
  
    // Selection cuts
  TCut cuts = "";//pfnjet[JetIdx];
  
  printf("Defining interesting observables\n");
  
    // S e t t i n g s    f o r   m a x i m u m   l i k e l i h o o d   e s t i m a t i o n
    // ------------------------------------------------------------------------------------
  
    // Probabilities to have X b-quarks in the tt-like final state **DERIVED FROM MC SIMULATIONS**
    //  float e0bq_[NbOfChannels][NbOfJetBins];// = {0.00873498,0.00638495,0.00734619};//{0.0247448,0.00873498,0.00638495,0.00734619};
    //  float e1bq_[NbOfChannels][NbOfJetBins];// = {0.179928,0.145577,0.134986};//{0.379931,0.179928,0.145577,0.134986};
    //  float e2bq_[NbOfChannels][NbOfJetBins];// = {0.808598,0.83759,0.836547};//{0.594488,0.808598,0.83759,0.836547};
  
    // Initial values for minimization
  /*
   float** init_Nttlike = new float*[NbOfChannels];
   float** init_Nvlike = new float*[NbOfChannels];
   for (int c=0; c<NbOfChannels; c++) {
   init_Nttlike[c] = new float[NbOfJetBins];//{0.,24500.,9850.,4015.};
   init_Nvlike[c] = new float[NbOfJetBins];
   if (NbOfJetBins>=1) {
   init_Nttlike[c][0] =  24500.;
   init_Nvlike[c][0] = 23400.;
   }
   if (NbOfJetBins>=2) {
   init_Nttlike[c][1] = 9850.;
   init_Nvlike[c][1] =  3900.;
   }
   if (NbOfJetBins>=3) {
   init_Nttlike[c][2] =   4015.;
   init_Nvlike[c][2] =  750.;//{0.,23400.,3900., 750.};
   }
   }
   */  
  float** init_Eb = new float*[NbOfChannels];
  float** init_Eudsc = new float*[NbOfChannels];
  float** init_Euds = new float*[NbOfChannels];
    // Mean and error values for constraints
  float** Eb_const_mean = new float*[NbOfChannels];
  float** Eb_const_error = new float*[NbOfChannels];
  float** Eudsc_const_mean = new float*[NbOfChannels];
  float** Eudsc_const_error = new float*[NbOfChannels];
  float** Euds_const_mean = new float*[NbOfChannels];
  float** Euds_const_error = new float*[NbOfChannels];
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    init_Eb[c] = new float[NbOfWP];
    init_Eudsc[c] = new float[NbOfWP];
    init_Euds[c] = new float[NbOfWP];
      // Mean and error values for constraints
    Eb_const_mean[c] = new float[NbOfWP];
    Eb_const_error[c] = new float[NbOfWP];
    Eudsc_const_mean[c] = new float[NbOfWP];
    Eudsc_const_error[c] = new float[NbOfWP];
    Euds_const_mean[c] = new float[NbOfWP];
    Euds_const_error[c] = new float[NbOfWP];
    
  }
  
  /*
   Double_t** number_events_mistagLight = new Double_t*[NbOfChannels];
   Double_t** number_jets_mistagLight = new Double_t*[NbOfChannels];
   Double_t** sum_events_mistagLight = new Double_t*[NbOfChannels];
   Double_t** sum_jets_mistagLight = new Double_t*[NbOfChannels];
   Double_t** number_events_SFLight = new Double_t*[NbOfChannels];
   Double_t** number_jets_SFLight = new Double_t*[NbOfChannels];
   Double_t** sum_events_SFLight = new Double_t*[NbOfChannels];
   Double_t** sum_jets_SFLight = new Double_t*[NbOfChannels];
   */
  Double_t** nEffUds_pass = new Double_t*[NbOfChannels];
  Double_t** nEffUds_tot  = new Double_t*[NbOfChannels];
  Double_t** nEffUds_pass_tmp = new Double_t*[NbOfChannels];
  Double_t** nEffUds_tot_tmp= new Double_t*[NbOfChannels];
  Double_t** nEffUdsc_pass = new Double_t*[NbOfChannels];
  Double_t** nEffUdsc_tot  = new Double_t*[NbOfChannels];
  
  Double_t** nEffb_tt_pass = new Double_t*[NbOfChannels];
  Double_t** nEffb_tt_tot = new Double_t*[NbOfChannels];
  Double_t** nEffb_w_pass = new Double_t*[NbOfChannels];
  Double_t** nEffb_w_tot = new Double_t*[NbOfChannels];
  
  for(UInt_t c=0; c<NbOfChannels; c++) {
    /*
     number_events_mistagLight[c] = new Double_t[NbOfWP];
     number_jets_mistagLight[c] = new Double_t[NbOfWP];
     sum_events_mistagLight[c] = new Double_t[NbOfWP];
     sum_jets_mistagLight[c] = new Double_t[NbOfWP];
     number_events_SFLight[c] = new Double_t[NbOfWP];
     number_jets_SFLight[c] = new Double_t[NbOfWP];
     sum_events_SFLight[c] = new Double_t[NbOfWP];
     sum_jets_SFLight[c] = new Double_t[NbOfWP];
     */  
    nEffUds_pass    [c] = new Double_t[NbOfWP];
    nEffUds_tot     [c] = new Double_t[NbOfWP];
    nEffUds_pass_tmp[c] = new Double_t[NbOfWP];
    nEffUds_tot_tmp [c] = new Double_t[NbOfWP];
    nEffUdsc_pass   [c] = new Double_t[NbOfWP];
    nEffUdsc_tot    [c] = new Double_t[NbOfWP];
    nEffb_tt_pass   [c] = new Double_t[NbOfWP];
    nEffb_tt_tot    [c] = new Double_t[NbOfWP];
    nEffb_w_pass    [c] = new Double_t[NbOfWP];
    nEffb_w_tot     [c] = new Double_t[NbOfWP];
    for (UInt_t i=0; i<NbOfWP; i++) {
      /*
       number_events_mistagLight[c][i] = 0.;
       number_jets_mistagLight[c][i] = 0.;
       sum_events_mistagLight[c][i] = 0.;
       sum_jets_mistagLight[c][i] = 0.;
       number_events_SFLight[c][i] = 0.;
       number_jets_SFLight[c][i] = 0.;
       sum_events_SFLight[c][i] = 0.;
       sum_jets_SFLight[c][i] = 0.;
       */    
      nEffUds_pass[c][i] = 0.;
      nEffUds_tot[c][i]  = 0.;
      nEffUds_pass_tmp[c][i] = 0.;
      nEffUds_tot_tmp[c][i]= 0.;
      nEffUdsc_pass[c][i] = 0.;
      nEffUdsc_tot[c][i]  = 0.;
      
      nEffb_tt_pass[c][i] = 0.;
      nEffb_tt_tot[c][i] = 0.;
      nEffb_w_pass[c][i] = 0.;
      nEffb_w_tot[c][i] = 0.;    
    }
  }
  
  printf("Defining (or reloading) histograms and useful MC information\n");
  
  /*
   TH1F*** eXbq__ = new TH1F**[NbOfChannels];
   for (UInt_t c=0; c<NbOfChannels ; c++) {
   eXbq__[c] = new TH1F*[NbOfJetBins];
   for (UInt_t i=0; i<NbOfJetBins; i++) {
   eXbq__[c][i] = NULL;
   }
   }
   */
  
  
  
  TH1I** histoNJ_data   = new TH1I*[NbOfChannels]; for (UInt_t c=0; c<NbOfChannels; c++) { histoNJ_data[c]   = new TH1I((std::string()+"JetMult_data"+channelSuffix[c]).c_str(),"JetMult_datas",50,0,50); }
  TH1I** histoNJ_ttjets = new TH1I*[NbOfChannels]; for (UInt_t c=0; c<NbOfChannels; c++) { histoNJ_ttjets[c] = new TH1I((std::string()+"JetMult_ttjets"+channelSuffix[c]).c_str(),"JetMult_ttjets",50,0,50); }
  TH1I** histoNJ_wjets  = new TH1I*[NbOfChannels]; for (UInt_t c=0; c<NbOfChannels; c++) { histoNJ_wjets[c]  = new TH1I((std::string()+"JetMult_wjets"+channelSuffix[c]).c_str(),"JetMult_wjets",50,0,50); }
  TH1I** histoNJ_zjets  = new TH1I*[NbOfChannels]; for (UInt_t c=0; c<NbOfChannels; c++) { histoNJ_zjets[c]  = new TH1I((std::string()+"JetMult_zjets"+channelSuffix[c]).c_str(),"JetMult_zjets",50,0,50); }
  TH1I** histoNJ_stjets = new TH1I*[NbOfChannels]; for (UInt_t c=0; c<NbOfChannels; c++) { histoNJ_stjets[c] = new TH1I((std::string()+"JetMult_stjets"+channelSuffix[c]).c_str(),"JetMult_stjets",50,0,50); }
  
    // Counters for the b-tagged jet multiplicity
    // Histograms for the b-tagged jet multiplicity
  TH1F*** hMC_Nominal_N = new TH1F**[NbOfChannels]; for (UInt_t c=0; c<NbOfChannels; c++) { hMC_Nominal_N[c] = new TH1F*[NbOfJetBins]; for (UInt_t i=0; i<NbOfJetBins; i++) {
    hMC_Nominal_N[c][i] = new TH1F((std::string()+"MC_Nominal_N"+channelSuffix[c]+jetSuffix[i]).c_str() ,"MC_Nominal_N;", 20, 0, 20);
  }
  }
  
  TH1F*** hMC_Nominal_E = new TH1F**[NbOfChannels]; for (UInt_t c=0; c<NbOfChannels; c++) { hMC_Nominal_E[c] = new TH1F*[NbOfWP]; for (UInt_t j=0; j<NbOfWP; j++) {
    hMC_Nominal_E[c][j] = new TH1F((std::string()+"MC_Nominal_eff"+channelSuffix[c]+wpSuffix[j]).c_str() ,"MC_Nominal_eff", 20, 0, 20);
  }
  }
  TH1D**** hNbtaggedJets_ttjets = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_ttjets[c] = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_ttjets[c][i] = new TH1D*[NbOfWP]; } }
  TH1D**** hNbtaggedJets_wjets  = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_wjets[c]  = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_wjets[c][i]  = new TH1D*[NbOfWP]; } }
  TH1D**** hNbtaggedJets_zjets  = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_zjets[c]  = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_zjets[c][i]  = new TH1D*[NbOfWP]; } }
  TH1D**** hNbtaggedJets_stjets = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_stjets[c] = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_stjets[c][i] = new TH1D*[NbOfWP]; } }
  TH1D**** hNbtaggedJets_vbjets = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_vbjets[c] = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_vbjets[c][i] = new TH1D*[NbOfWP]; } }
  TH1D**** hNbtaggedJets_qcd    = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_qcd[c]    = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_qcd[c][i] = new TH1D*[NbOfWP]; } }
  TH1D**** hNbtaggedJets_qcd_shape=new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_qcd_shape[c]=new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_qcd_shape[c][i]=new TH1D*[NbOfWP]; } }
  TH1D**** hNbtaggedJets_vvjets = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_vvjets[c] = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_vvjets[c][i] = new TH1D*[NbOfWP]; } }
  TH1D**** hWbb_In_WJets        = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hWbb_In_WJets[c]        = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hWbb_In_WJets[c][i]  = new TH1D*[NbOfWP]; } }
  TH1D**** hWcx_In_WJets        = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hWcx_In_WJets[c]        = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hWcx_In_WJets[c][i]  = new TH1D*[NbOfWP]; } }
  TH1D**** hWLFg_In_WJets       = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hWLFg_In_WJets[c]       = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hWLFg_In_WJets[c][i] = new TH1D*[NbOfWP]; } }
  TH1D**** hNbtaggedJets_rare   = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_rare[c]   = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_rare[c][i] = new TH1D*[NbOfWP]; } }
  TH1D**** hNbtaggedJets_signal = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_signal[c] = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_signal[c][i] = new TH1D*[NbOfWP]; } }
  TH1D**** hNbtaggedJets        = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets[c] = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) hNbtaggedJets[c][i] = new TH1D*[NbOfWP]; }
  for(UInt_t c=0; c<NbOfChannels; c++) { 
      //    NbtaggedJets[c] = 0;
      //    nB_Flavoured_Jets[c] = 0;
      //    nC_Flavoured_Jets[c] = 0;
      //    nLightAndGluon_Flavoured_Jets[c] = 0;
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      for (UInt_t j=0; j<NbOfWP; j++) {
        hNbtaggedJets[c][i][j] = NULL;
        hNbtaggedJets_ttjets[c][i][j] = NULL;
        hNbtaggedJets_wjets[c][i][j] = NULL;
        hNbtaggedJets_zjets[c][i][j] = NULL;
        hNbtaggedJets_stjets[c][i][j] = NULL;
        hNbtaggedJets_vbjets[c][i][j] = NULL;
        hNbtaggedJets_qcd[c][i][j] = NULL;
        hNbtaggedJets_qcd_shape[c][i][j] = NULL;
        hNbtaggedJets_vvjets[c][i][j] = NULL;
        hWbb_In_WJets[c][i][j] = NULL;
        hWcx_In_WJets[c][i][j] = NULL;
        hWLFg_In_WJets[c][i][j] = NULL;
        hNbtaggedJets_rare[c][i][j] = NULL;
        hNbtaggedJets_signal[c][i][j] = NULL;
      }
    }
  }
  
  TH1D**** hNbtaggedJets_ttlike = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_ttlike[c] = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) hNbtaggedJets_ttlike[c][i] = new TH1D*[NbOfWP]; }
  TH1D**** hNbtaggedJets_vlike  = new TH1D***[NbOfChannels]; for(UInt_t c=0; c<NbOfChannels; c++) { hNbtaggedJets_vlike[c] = new TH1D**[NbOfJetBins]; for(UInt_t i=0; i<NbOfJetBins ; i++) hNbtaggedJets_vlike[c][i] = new TH1D*[NbOfWP]; }
  for(UInt_t c=0; c<NbOfChannels; c++) { 
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      for (UInt_t j=0; j<NbOfWP; j++) {
        hNbtaggedJets_ttlike[c][i][j] = NULL; //new TH1D((std::string()+"hNbtaggedJets_ttlike"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),";B-tagged jet multiplicity",10,0,10);
        hNbtaggedJets_vlike[c][i][j]  = NULL; //new TH1D((std::string()+"hNbtaggedJets_vlike"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str() ,";B-tagged jet multiplicity",10,0,10);
      }
    }
  }
  
  /*
   TH2F *ebVMC_histo = NULL; 
   TH2F *ebTTMC_histo = NULL;
   TH2F *eudsMC_histo = NULL; 
   TH2F *eudscMC_histo = NULL;
   */
  cout << "***************************************************************" << endl;
  cout << "                 Computing the efficiencies                    " << endl;
  cout << "***************************************************************" << endl;

  TH1D*** eXbq__ = new TH1D**[NbOfChannels] ;
  TH1D*** btaggingEff_tt_ = new TH1D**[NbOfChannels] ;
  TH1D*** btaggingEff_vlike_ = new TH1D**[NbOfChannels] ;
  printf("\n\n !!!! Here for loading the efficiencies from MC (Leg2bis) and applying the SF !!!!\n\n");
  char** channelLabel = new char*[NbOfChannels];
  for (UInt_t i=0; i<NbOfChannels; i++) {
    channelLabel[i] = new char[50];
    sprintf(channelLabel[i], "%s", (std::string()
                                    +(((i+1)&(1<<0))?"_mu":"")
                                    +(((i+1)&(1<<1))?"":"")).c_str());
  }
  
  std::vector<std::vector<Double_t> > effTT(3,vector<Double_t>(NbOfChannels,0.));
  std::vector<std::vector<Double_t> > errTT(3,vector<Double_t>(NbOfChannels,0.));
  std::vector<std::vector<Double_t> > effV(3,vector<Double_t>(NbOfChannels,0.));
  std::vector<std::vector<Double_t> > errV(3,vector<Double_t>(NbOfChannels,0.));
  
  {
    std::string *flavourKind = new std::string[3] ;
    flavourKind[0] = "_guds" ;
    flavourKind[1] = "_c" ;
    flavourKind[2] = "_b" ;
    
    printf("\nFor TT-like category\n\n");
    for(UInt_t chan=0; chan<NbOfChannels; chan++) {
      printf("\nChannel %d : \n", chan);
      TList *effListNominal = new TList();
      TList *effListNominal_NoSF = new TList();
      TList *effListPlus   = new TList();
      TList *effListMinus  = new TList();
      std::vector<Double_t> weights;
      for (UInt_t fl=0; fl<3; fl++) {
        
        for (UInt_t c=0; ((UInt_t)(1<<c))<=chan+1; c++) {
          Int_t channel=(1<<c)-1 ;
          if (((channel+1)&(chan+1))!=0) {
            printf("channel : %d     chan : %d\n", channel, chan);
            Bool_t nominalListWasEmpty = kTRUE;
            if (effListNominal->GetSize()!=0) {
              nominalListWasEmpty = /*kFALSE*/ kTRUE; //Not the same channel 
            }
            for (UInt_t k=0; k<files[0].size(); k++) {
              TEfficiency *eff_tmp = (TEfficiency*) fBtagEff->Get((samples_labelBQ[0][k]+flavourKind[fl]+channelLabel[channel]).c_str());
              printf("%s ",(samples_labelBQ[0][k]+flavourKind[fl]+channelLabel[channel]).c_str());
              TEfficiency *eff0bin = NULL;
              TEfficiency *eff1bin = NULL;
              if (eff_tmp!=NULL) {
                //                eff0bin = average(eff_tmp, ((fl==2)?5:((fl==1)?4:0)), std::string(runPeriod.Data()), 0, kTRUE);
                eff0bin = average(eff_tmp, ((fl==2)?5:((fl==1)?4:0)), std::string(runPeriod.Data()), 0, kTRUE);
                eff1bin = average(eff_tmp, ((fl==2)?5:((fl==1)?4:0)), std::string(runPeriod.Data()), 0, kFALSE);
                  //              eff1bin->GetPassedHistogram()->Print();
                  //              eff1bin->GetTotalHistogram()->Print();
                if (nominalListWasEmpty == kTRUE) {
                  effListNominal_NoSF->Add((TEfficiency*) eff0bin);
                  effListNominal->Add((TEfficiency*) eff1bin);
                    //                  ((TEfficiency*) effListNominal->At(effListNominal->GetSize()-1))->GetPassedHistogram()->Print();
                    //                  ((TEfficiency*) effListNominal->At(effListNominal->GetSize()-1))->GetTotalHistogram()->Print();
                  weights.push_back( /* IntLumi[c]* */samples_weights[0][k]);
                } else {
                  TH1* tmp = (TH1*) ((TEfficiency*) effListNominal->At(k))->GetPassedHistogram()->Clone();
                  tmp->Add(eff1bin->GetPassedHistogram());
                  ((TEfficiency*) effListNominal->At(k))->SetPassedHistogram(*tmp, "f"); //exclusive OR
                    //                  ((TEfficiency*) effListNominal->At(k))->Add(* (TEfficiency*) eff1bin);
                  tmp = (TH1*) ((TEfficiency*) effListNominal_NoSF->At(k))->GetPassedHistogram()->Clone();
                  tmp->Add(eff0bin->GetPassedHistogram());
                  ((TEfficiency*) effListNominal_NoSF->At(k))->SetPassedHistogram(*tmp, "f");
                }
                
                eff1bin = average(eff_tmp, (fl==2?5:(fl==1?4:0)), std::string(runPeriod.Data()), -1, kFALSE);
                if (nominalListWasEmpty == kTRUE) {
                  effListMinus->Add((TEfficiency*) eff1bin);
                    //                weights.push_back( /* IntLumi[c]* */samples_weights[0][k]);
                } else {
                  TH1* tmp = (TH1*) ((TEfficiency*) effListMinus->At(k))->GetPassedHistogram()->Clone();
                  tmp->Add(eff1bin->GetPassedHistogram());
                  ((TEfficiency*) effListMinus->At(k))->SetPassedHistogram(*tmp, "f"); //exclusive OR
                }
                eff1bin = average(eff_tmp, (fl==2?5:(fl==1?4:0)), std::string(runPeriod.Data()), +1, kFALSE);
                if (nominalListWasEmpty == kTRUE) {
                  effListPlus->Add((TEfficiency*) eff1bin);
                    //                weights.push_back( /* IntLumi[c]* */samples_weights[0][k]);
                } else {
                  TH1* tmp = (TH1*) ((TEfficiency*) effListPlus->At(k))->GetPassedHistogram()->Clone();
                  tmp->Add(eff1bin->GetPassedHistogram());
                  ((TEfficiency*) effListPlus->At(k))->SetPassedHistogram(*tmp, "f"); //exclusive OR
                }
              } 
            } 
            /*
             printf("Nominal size : %d\n", effListNominal->GetSize());
             for (UInt_t lll=0; lll<effListNominal->GetSize(); lll++) {
             ((TEfficiency*) effListNominal->At(lll))->GetTotalHistogram()->Print();
             (dynamic_cast<TEfficiency*>(effListNominal->At(lll)))->GetTotalHistogram()->Print();
             ((TEfficiency*) effListNominal->At(lll))->GetPassedHistogram()->Print();
             }
             */
            printf("\n");
            for (UInt_t lll=0; lll<(UInt_t) effListNominal->GetSize(); lll++) {
              printf("%lf ", ((TEfficiency*) effListNominal->At(lll))->GetEfficiency(1));
            }
            printf("\n");
            for (UInt_t lll=0; lll<weights.size(); lll++) {
              printf("%lf ", weights[lll]);
            }
            printf("\n");
            UInt_t bin = 0 /*4*/;
              //            printf("a %d -> [0] %lu\n", effListNominal->GetSize(), (unsigned long) effListNominal->At(0));
            TGraphAsymmErrors* tgn = NULL;
            tgn = TEfficiency::Combine(effListNominal, "mode", weights.size(), & weights[0]);
            TGraphAsymmErrors* tgn_NoSF = TEfficiency::Combine(effListNominal_NoSF, "mode", weights.size(), & weights[0]);
              //            effListNominal->Clear();
              //            printf("b\n");
            TGraphAsymmErrors* tgp = NULL;
            tgp = TEfficiency::Combine(effListPlus, "mode", weights.size(), & weights[0]);
              //            printf("c\n");
            TGraphAsymmErrors* tgm = NULL;
            tgm = TEfficiency::Combine(effListMinus, "mode", weights.size(), & weights[0]);
            /*          TCanvas c("test");
             c.cd();
             tgn->Draw();
             c.Print((std::string()+"tgnomTT"+flavourKind[fl]+channelSuffix[chan]+".pdf").c_str());
             */
            Double_t x=0. , ynom=0. , yminus=0., yplus=0., ynom_NoSF=0. ;
            if (tgn != NULL) {
              tgn->GetPoint(bin, x,ynom);
              tgn_NoSF->GetPoint(bin, x,ynom_NoSF);
              Double_t ylow = tgn->GetErrorYlow(bin);
              Double_t yhigh = tgn->GetErrorYhigh(bin);
              effTT[fl][chan] = ynom;
              printf("effTT[%d][%d] = %lf // %lf\n", fl, chan, ynom, ynom_NoSF);
              tgm->GetPoint(bin, x,yminus);
              tgp->GetPoint(bin, x,yplus);
              Double_t errrelM = sqrt( pow((ynom-yminus)/ynom,2.) + pow(ylow,2.));
              Double_t errrelP = sqrt( pow((yplus-ynom)/ynom,2.) + pow(yhigh,2.));
              errTT[fl][chan] = ynom * (errrelM<errrelP ? errrelP : errrelM);
            }
          }
        }
        if (fl==1 || fl==2) { // To make the guds and c categories together, while keeping results for guds alone.
          printf(" --> Reset\n");
          weights.clear();
          effListNominal->Clear();
          effListNominal_NoSF->Clear();
          effListPlus->Clear();
          effListMinus->Clear();
        }        
      }
    }
    
    printf("\nFor V-like category\n\n");
    for(UInt_t chan=0; chan<NbOfChannels; chan++) {
      printf("\nChannel %d : \n", chan);
      TList *effListNominal = new TList();
      TList *effListNominal_NoSF = new TList();
      TList *effListPlus   = new TList();
      TList *effListMinus  = new TList();
      std::vector<Double_t> weights;
      for (UInt_t fl=0; fl<3; fl++) {

        for (UInt_t c=0; ((UInt_t)(1<<c))<=chan+1; c++) {
          Int_t channel=(1<<c)-1 ;
          if (((channel+1)&(chan+1))==0) {
            continue;
          }
          Bool_t nominalListWasEmpty = kTRUE;
          if (effListNominal->GetSize()!=0) {
            nominalListWasEmpty = /*kFALSE*/ kTRUE; //Not the same channel 
          }
          for (UInt_t k=0; k<files[1].size(); k++) {
            TEfficiency *eff_tmp = (TEfficiency*) fBtagEff->Get((samples_labelBQ[1][k]+flavourKind[fl]+channelLabel[channel]).c_str());
            printf("%s ",(samples_labelBQ[1][k]+flavourKind[fl]+channelLabel[channel]).c_str());
            TEfficiency *eff0bin = NULL;
            TEfficiency *eff1bin = NULL;
            if (eff_tmp!=NULL) {
              eff0bin = average(eff_tmp, ((fl==2)?5:((fl==1)?4:0)), std::string(runPeriod.Data()), 0, kTRUE);
              eff1bin = average(eff_tmp, (fl==2?5:(fl==1?4:0)), std::string(runPeriod.Data()), 0, kFALSE);
              if (nominalListWasEmpty == kTRUE) {
                effListNominal_NoSF->Add((TEfficiency*) eff0bin);
                effListNominal->Add((TEfficiency*) eff1bin);
                weights.push_back( /* IntLumi[c]* */samples_weights[1][k]);
              } else {
                TH1* tmp = (TH1*) ((TEfficiency*) effListNominal->At(k))->GetPassedHistogram()->Clone();
                tmp->Add(eff1bin->GetPassedHistogram());
                ((TEfficiency*) effListNominal->At(k))->SetPassedHistogram(*tmp, "f"); //exclusive OR
                tmp = (TH1*) ((TEfficiency*) effListNominal_NoSF->At(k))->GetPassedHistogram()->Clone();
                tmp->Add(eff0bin->GetPassedHistogram());
                ((TEfficiency*) effListNominal_NoSF->At(k))->SetPassedHistogram(*tmp, "f");
              }
              eff1bin = average(eff_tmp, (fl==2?5:(fl==1?4:0)), std::string(runPeriod.Data()), -1, kFALSE);
              if (nominalListWasEmpty == kTRUE) {
                effListMinus->Add((TEfficiency*) eff1bin);
                  //                weights.push_back( /* IntLumi[c]* */samples_weights[1][k]);
              } else {
                TH1* tmp = (TH1*) ((TEfficiency*) effListMinus->At(k))->GetPassedHistogram()->Clone();
                tmp->Add(eff1bin->GetPassedHistogram());
                ((TEfficiency*) effListMinus->At(k))->SetPassedHistogram(*tmp, "f"); //exclusive OR
              }
              eff1bin = average(eff_tmp, (fl==2?5:(fl==1?4:0)), std::string(runPeriod.Data()), +1, kFALSE);
              if (nominalListWasEmpty == kTRUE) {
                effListPlus->Add((TEfficiency*) eff1bin);
                  //                weights.push_back( /* IntLumi[c]* */samples_weights[1][k]);
              } else {
                TH1* tmp = (TH1*) ((TEfficiency*) effListPlus->At(k))->GetPassedHistogram()->Clone();
                tmp->Add(eff1bin->GetPassedHistogram());
                ((TEfficiency*) effListPlus->At(k))->SetPassedHistogram(*tmp, "f"); //exclusive OR
              }
            } 
          }
          printf("-> Names\n");
          for (UInt_t lll=0; lll<(UInt_t) effListNominal->GetSize(); lll++) {
            printf("%s\t", ((TEfficiency*) effListNominal->At(lll))->GetName());
          }
          printf("\n");
          printf("-> %d\n", effListNominal->GetSize());
          for (UInt_t lll=0; lll<(UInt_t) effListNominal->GetSize(); lll++) {
            printf("%lf (No SF : %lf)\t", ((TEfficiency*) effListNominal->At(lll))->GetEfficiency(1), ((TEfficiency*) effListNominal_NoSF->At(lll))->GetEfficiency(1));
          }
          printf("\n");
          printf("-> Numerators\n");
          for (UInt_t lll=0; lll<(UInt_t) effListNominal->GetSize(); lll++) {
            printf("%lf\t", ((TEfficiency*) effListNominal->At(lll))->GetPassedHistogram()->GetBinContent(1));
          }
          printf("\n");
           printf("-> Denominators\n");
          for (UInt_t lll=0; lll<(UInt_t) effListNominal->GetSize(); lll++) {
            printf("%lf\t", ((TEfficiency*) effListNominal->At(lll))->GetTotalHistogram()->GetBinContent(1));
          }
          printf("\n");
           for (UInt_t lll=0; lll<weights.size(); lll++) {
            printf("%lf\t", weights[lll]);
          }
          printf("\n");
          UInt_t bin = 0 /*5*/;
          TGraphAsymmErrors* tgn = TEfficiency::Combine(effListNominal, "mode", weights.size(), & weights[0]);
          TGraphAsymmErrors* tgn_NoSF = TEfficiency::Combine(effListNominal_NoSF, "mode", weights.size(), & weights[0]);
          TGraphAsymmErrors* tgp = TEfficiency::Combine(effListPlus, "mode", weights.size(), & weights[0]);
          TGraphAsymmErrors* tgm = TEfficiency::Combine(effListMinus, "mode", weights.size(), & weights[0]);
          /*          TCanvas c("test");
           c.cd();
           tgn->Draw();
           c.Print((std::string()+"tgnomV"+flavourKind[fl]+channelSuffix[chan]+".pdf").c_str());
           */
          Double_t x=0. , ynom=0. , yminus=0., yplus=0. , ynom_NoSF=0.;
          if (tgn != NULL) {
            tgn->GetPoint(bin, x,ynom);
            tgn_NoSF->GetPoint(bin, x,ynom_NoSF);
            Double_t ylow = tgn->GetErrorYlow(bin);
            Double_t yhigh = tgn->GetErrorYhigh(bin);
            effV[fl][chan] = ynom;
            printf("effV[%d][%d] = %lf // %lf\n", fl, chan, ynom, ynom_NoSF);
            tgm->GetPoint(bin, x,yminus);
            tgp->GetPoint(bin, x,yplus);
            Double_t errrelM = sqrt( pow((ynom-yminus)/ynom,2.) + pow(ylow,2.));
            Double_t errrelP = sqrt( pow((yplus-ynom)/ynom,2.) + pow(yhigh,2.));
            errV[fl][chan] = ynom * (errrelM<errrelP ? errrelP : errrelM);
          }
        }
        if (fl==1 || fl==2) { // To make the guds and c categories  non-splitted
          printf("\n");
          weights.clear();
          effListNominal->Clear();
          effListNominal_NoSF->Clear();
          effListPlus->Clear();
          effListMinus->Clear();
        }            
      }
    }
    
    for (UInt_t fl=0; fl<3; fl++) {
      for(UInt_t chan=0; chan<NbOfChannels; chan++) {
        printf("TT[%u][%u] : %lf \\pm %lf \n", fl, chan, effTT[fl][chan], errTT[fl][chan]);
      }
      printf("\n");
    }
    for (UInt_t fl=0; fl<3; fl++) {
      for(UInt_t chan=0; chan<NbOfChannels; chan++) {
        printf("V[%u][%u] : %lf \\pm %lf \n", fl, chan, effV[fl][chan], errV[fl][chan]);
      }
      printf("\n");
    }
    
    
    for(UInt_t c=0; c<NbOfChannels; c++) {
      eXbq__[c] = new TH1D*[NbOfJetBins];
      btaggingEff_tt_[c] = new TH1D*[NbOfWP];
      btaggingEff_vlike_[c] = new TH1D*[NbOfWP];
      for (UInt_t i=0; i<NbOfJetBins; i++) {
        for (UInt_t j=0; j<NbOfWP; j++) {
          if (i==0) {
            btaggingEff_tt_[c][j] = NULL;
            btaggingEff_vlike_[c][j] = NULL;
          }
          if (j==0) {
            eXbq__[c][i] = NULL;
          }
        }
      }
    }
    TFile *fi = NULL;
    if (estimateQCDApart) {
      fi = new TFile(qcdShapeFile.c_str(), "READ");
      if (fi==NULL) {
        printf("Unable to open file \"%s\" to read the shapes for QCD estimation\n", qcdShapeFile.c_str());
        exit(1);
      }
    }
    
    for(UInt_t chan=0; chan<NbOfChannels; chan++) {
      for (UInt_t c=0; ((UInt_t)(1<<c))<=chan+1; c++) {
        Int_t channel=(1<<c)-1 ;
        if (((channel+1)&(chan+1))==0) {
          continue;
        }
        printf("\nchan : %d  ;   c : %d \n", chan, c);
        
        for (UInt_t i=0; i<NbOfJetBins; i++) {
          for (UInt_t ib=0; ib<NbOfBtagBins; ib++) {
            for (UInt_t j=0; j<NbOfWP; j++) {
              printf("%s\n", (std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str());
              //ttlike
            printf(" TTlike");
            if (files[0].size() != 0) {
              if (hNbtaggedJets_ttlike[chan][i][j]==NULL) {
                hNbtaggedJets_ttlike[chan][i][j] = (TH1D*) files[0][0]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                hNbtaggedJets_ttlike[chan][i][j]->SetDirectory(NULL);
                hNbtaggedJets_ttlike[chan][i][j]->Reset();
              }
              for (UInt_t k=0; k<files[0].size(); k++) {
                hNbtaggedJets_ttlike[chan][i][j]->Add((TH1D*) files[0][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[0][k]);
              }
              if (j==0) {
                printf("(BQuarks)");
                /*
                { // "Should be value"
                  if (eXbq__[chan][i]==NULL) {
                    eXbq__[chan][i] = (TH1D*) hNbtaggedJets_ttlike[chan][i][0]->Clone();
                    eXbq__[chan][i]->SetDirectory(NULL);
                    eXbq__[chan][i]->Reset();
                  }
                  if (chan==0) {
                      //SemiMuon
                    Double_t a[7]={0,0,0,0, 0.022348, 0.016215, 0.013229 };
                    Double_t b[7]={0,0,0,0, 0.283191, 0.238223, 0.211845 };
                    Double_t c[7]={0,0,0,0, 0.686278, 0.725439, 0.732775 };
                    eXbq__[chan][i]->Fill(0., a[MinNbOfJets+i]); 
                    eXbq__[chan][i]->Fill(1., b[MinNbOfJets+i]); 
                    eXbq__[chan][i]->Fill(2., c[MinNbOfJets+i]); 
                    eXbq__[chan][i]->Fill(3., 1.-eXbq__[chan][i]->Integral(0,-1));
                  } else if (chan==1) {
                      //SemiElectron
                    Double_t a[7]={0,0,0,0, 0.023797, 0.017614, 0.013601 };
                    Double_t b[7]={0,0,0,0, 0.288904, 0.242580, 0.217566 };
                    Double_t c[7]={0,0,0,0, 0.679199, 0.719735, 0.727453 };
                    eXbq__[chan][i]->Fill(0., a[MinNbOfJets+i]); 
                    eXbq__[chan][i]->Fill(1., b[MinNbOfJets+i]); 
                    eXbq__[chan][i]->Fill(2., c[MinNbOfJets+i]); 
                    eXbq__[chan][i]->Fill(3., 1.-eXbq__[chan][i]->Integral(0,-1));
                  } else if (chan==2) {
                      //SemiMuonSemiElectron
                    Double_t a[7]={0,0,0,0, 0.022928, 0.016823, 0.013386 };
                    Double_t b[7]={0,0,0,0, 0.284861, 0.240341, 0.214254 };
                    Double_t c[7]={0,0,0,0, 0.684076, 0.722721, 0.730549 };
                    eXbq__[chan][i]->Fill(0., a[MinNbOfJets+i]); 
                    eXbq__[chan][i]->Fill(1., b[MinNbOfJets+i]); 
                    eXbq__[chan][i]->Fill(2., c[MinNbOfJets+i]); 
                    eXbq__[chan][i]->Fill(3., 1.-eXbq__[chan][i]->Integral(0,-1));
                  }
                }
                */
                 if (eXbq__[chan][i]==NULL) {
                 eXbq__[chan][i] = (TH1D*) files[0][0]->Get((std::string()+ (isNewNaming==kTRUE ? "BJet_Multiplicity" : "B_Quarks")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]).c_str())->Clone();
                 printf("%s\n",eXbq__[chan][i]->GetName());
                 eXbq__[chan][i]->SetDirectory(NULL);
                 eXbq__[chan][i]->Reset();
                 }
                 
                
                 for (UInt_t k=0; k<files[0].size(); k++) {
                   printf("%s\n",((TH1D*) files[0][k]->Get((std::string()+ (isNewNaming==kTRUE ? "BJet_Multiplicity" : "B_Quarks")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]).c_str()))->GetName());
                   eXbq__[chan][i]->Add( ((TH1D*) files[0][k]->Get((std::string()+ (isNewNaming==kTRUE ? "BJet_Multiplicity" : "B_Quarks")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]).c_str())), IntLumi[c]*samples_weights[0][k]); //eXbq
                   printf(" [%d][%d] %lf\n", chan, i, eXbq__[chan][i]->Integral(0,-1));
                 }
                
              }
              if (i==0) {
                /**
                 histo2d  =  new TH2D("guds_jets","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"guds_jets");
                 histovec2d.push_back(histomap2d);  
                 histo2d  =  new TH2D("guds_tagged_jets","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"guds_tagged_jets");
                 histovec2d.push_back(histomap2d);  
                 histo2d  =  new TH2D("c_jets","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"c_jets");
                 histovec2d.push_back(histomap2d);  
                 histo2d  =  new TH2D("c_tagged_jets","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"c_tagged_jets");
                 histovec2d.push_back(histomap2d);  
                 histo2d  =  new TH2D("b_jets","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"b_jets");
                 histovec2d.push_back(histomap2d);  
                 histo2d  =  new TH2D("b_tagged_jets","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"b_tagged_jets");
                 histovec2d.push_back(histomap2d);  
                 
                 histo2d  =  new TH2D("guds_jets_mu","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"guds_jets_mu");
                 histovec2d.push_back(histomap2d);  
                 histo2d  =  new TH2D("guds_tagged_jets_mu","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"guds_tagged_jets_mu");
                 histovec2d.push_back(histomap2d);  
                 histo2d  =  new TH2D("c_jets_mu","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"c_jets_mu");
                 histovec2d.push_back(histomap2d);  
                 histo2d  =  new TH2D("c_tagged_jets_mu","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"c_tagged_jets_mu");
                 histovec2d.push_back(histomap2d);  
                 histo2d  =  new TH2D("b_jets_mu","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"b_jets_mu");
                 histovec2d.push_back(histomap2d);  
                 histo2d  =  new TH2D("b_tagged_jets_mu","",16,0,16,3,0,3);
                 histomap2d = make_pair(histo2d,"b_tagged_jets_mu");
                 histovec2d.push_back(histomap2d);  
                 
                 */
                printf("(BTaggingEff)");
                if (btaggingEff_tt_[chan][j]==NULL) {
                  btaggingEff_tt_[chan][j] = (TH1D*) files[0][0]->Get((std::string()+"Btagging_CSVM_"+channelLeg3Suffix[c]+wpSuffix[j]).c_str());
                  if (btaggingEff_tt_[chan][j] != NULL) {
                    btaggingEff_tt_[chan][j] = (TH1D*) btaggingEff_tt_[chan][j]/*->Clone()*/;
                    btaggingEff_tt_[chan][j]->SetDirectory(NULL);
                    btaggingEff_tt_[chan][j]->Reset();
                  }
                }
                if (btaggingEff_tt_[chan][j] != NULL) {
                  for (UInt_t k=0; k<files[0].size(); k++) {
                    btaggingEff_tt_[chan][j]->Add((TH1D*) files[0][0]->Get((std::string()+"Btagging_CSVM_"+channelLeg3Suffix[c]+wpSuffix[j]).c_str()), IntLumi[c]*samples_weights[0][k]);
                  }
                }
              }
            }
            
              //V-like
            printf(" Vlike");
            if (files[1].size() != 0) {
              if (hNbtaggedJets_vlike[chan][i][j]==NULL) {
                hNbtaggedJets_vlike[chan][i][j] = (TH1D*) files[1][0]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                hNbtaggedJets_vlike[chan][i][j]->SetDirectory(NULL);
                hNbtaggedJets_vlike[chan][i][j]->Reset();
              }
              for (UInt_t k=0; k<files[1].size(); k++) {
                hNbtaggedJets_vlike[chan][i][j]->Add((TH1D*) files[1][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[1][k]);
                  //                (TH1D*) files[1][0]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+wpSuffix[j]+btagJetMult_suffix).c_str());
              }
              if (i==0) {
                printf("(BTaggingEff)");
                if (btaggingEff_vlike_[chan][j]==NULL) {
                  btaggingEff_vlike_[chan][j] = (TH1D*) files[1][0]->Get((std::string()+"Btagging_CSVM_"+channelLeg3Suffix[c]+wpSuffix[j]).c_str());
                  if (btaggingEff_vlike_[chan][j] != NULL) {
                    btaggingEff_vlike_[chan][j] = (TH1D*) btaggingEff_vlike_[chan][j]/*->Clone()*/;
                    btaggingEff_vlike_[chan][j]->SetDirectory(NULL);
                    btaggingEff_vlike_[chan][j]->Reset();
                  }
                }
                if (btaggingEff_vlike_[chan][j] != NULL) {
                  for (UInt_t k=0; k<files[1].size(); k++) {
                    if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
                    if ((((c+1)&muonMask)!=0) && samples_names[1][k].find("QCD_mu")!=0) continue ;
                    if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
                    if ((((c+1)&electronMask)!=0) && samples_names[1][k].find("QCD_el")!=0) continue ;
                    btaggingEff_vlike_[chan][j]->Add((TH1D*) files[1][0]->Get((std::string()+"Btagging_CSVM_"+channelLeg3Suffix[c]+wpSuffix[j]).c_str()), IntLumi[c]*samples_weights[1][k]);
                  }
                }
              }
            }
            
              //V-bb
            printf(" Vbb");
            if (files[2].size() != 0) {
              if (hNbtaggedJets_vbjets[chan][i][j]==NULL) {
                hNbtaggedJets_vbjets[chan][i][j] = (TH1D*) files[2][0]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMultbb_prefix : "B_Jet_Multiplicity_bb")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                hNbtaggedJets_vbjets[chan][i][j]->SetDirectory(NULL);
                hNbtaggedJets_vbjets[chan][i][j]->Reset();
              }
              for (UInt_t k=0; k<files[2].size(); k++) {
                hNbtaggedJets_vbjets[chan][i][j]->Add((TH1D*) files[2][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMultbb_prefix : "B_Jet_Multiplicity_bb")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[2][k]);
                hNbtaggedJets_vlike[chan][i][j]->Add((TH1D*) files[2][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMultbb_prefix : "B_Jet_Multiplicity_bb")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), -1. * IntLumi[c]*samples_weights[2][k]);
              }
            }
            
              //Single Top
            printf(" ST");
            if (files[3].size() != 0) {
              if (hNbtaggedJets_stjets[chan][i][j]==NULL) {
                hNbtaggedJets_stjets[chan][i][j] = (TH1D*) files[3][0]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                hNbtaggedJets_stjets[chan][i][j]->SetDirectory(NULL);
                hNbtaggedJets_stjets[chan][i][j]->Reset();
              }
              for (UInt_t k=0; k<files[3].size(); k++) {
                hNbtaggedJets_stjets[chan][i][j]->Add((TH1D*) files[3][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[3][k]);
              }
            }
            
              //VV
            printf(" VV");
            if (files[4].size() != 0) {
              if (hNbtaggedJets_vvjets[chan][i][j]==NULL) {
                hNbtaggedJets_vvjets[chan][i][j] = (TH1D*) files[4][0]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                hNbtaggedJets_vvjets[chan][i][j]->SetDirectory(NULL);
                hNbtaggedJets_vvjets[chan][i][j]->Reset();
              }
              
              for (UInt_t k=0; k<files[4].size(); k++) {
                hNbtaggedJets_vvjets[chan][i][j]->Add((TH1D*) files[4][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[4][k]);
              }
            }
            
            printf(" MCPerProcess");
            for (UInt_t f=0; f<NbOfCategories; f++) {
              for (UInt_t k=0; k<files[f].size(); k++) {
                if(samples_names[f][k].find("TTJets")==0) {
                  if (hNbtaggedJets_ttjets[chan][i][j]==NULL) {
                    hNbtaggedJets_ttjets[chan][i][j] = (TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                    hNbtaggedJets_ttjets[chan][i][j]->SetDirectory(NULL);
                    hNbtaggedJets_ttjets[chan][i][j]->Reset();
                  }
                  hNbtaggedJets_ttjets[chan][i][j]->Add((TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[f][k]);
                }
                if((WnJets==kFALSE && samples_names[f][k].find("WJets")==0)
                   || ((WnJets==kTRUE)&&(samples_names[f][k].find("W1Jets")==0 || samples_names[f][k].find("W2Jets")==0 || samples_names[f][k].find("W3Jets")==0 || samples_names[f][k].find("W4Jets")==0))) {
                  if (hNbtaggedJets_wjets[chan][i][j]==NULL) {
                    hNbtaggedJets_wjets[chan][i][j] = (TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                    hNbtaggedJets_wjets[chan][i][j]->SetDirectory(NULL);
                    hNbtaggedJets_wjets[chan][i][j]->Reset();
                  }
                  hNbtaggedJets_wjets[chan][i][j]->Add((TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[f][k]);
                }
                if((DYnJets==kFALSE && samples_names[f][k].find("DYJets")==0)
                   || ((DYnJets==kTRUE)&&(samples_names[f][k].find("DY1Jets")==0 || samples_names[f][k].find("DY2Jets")==0 || samples_names[f][k].find("DY3Jets")==0 || samples_names[f][k].find("DY4Jets")==0))) {
                  if (hNbtaggedJets_zjets[chan][i][j]==NULL) {
                    hNbtaggedJets_zjets[chan][i][j] = (TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                    hNbtaggedJets_zjets[chan][i][j]->SetDirectory(NULL);
                    hNbtaggedJets_zjets[chan][i][j]->Reset();
                  }
                  hNbtaggedJets_zjets[chan][i][j]->Add((TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[f][k]);
                }
                // Will only be filled by this code if the QCD samples are in the V-like category
                if(samples_names[f][k].find("QCD")==0) {
                  if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
                  if ((((c+1)&muonMask)!=0) && samples_names[1][k].find("QCD_mu")!=0) continue ;
                  if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
                  if ((((c+1)&electronMask)!=0) && samples_names[1][k].find("QCD_el")!=0) continue ;
                  if (hNbtaggedJets_qcd[chan][i][j]==NULL) {
                    hNbtaggedJets_qcd[chan][i][j] = (TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                    hNbtaggedJets_qcd[chan][i][j]->SetDirectory(NULL);
                    hNbtaggedJets_qcd[chan][i][j]->Reset();
                  }
                  hNbtaggedJets_qcd[chan][i][j]->Add((TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[f][k]);
                }
              }
            }
            {
              //Extra-categories
              //QCD
              UInt_t f = NbOfCategories ;
              for (UInt_t k=0; k<files[f].size(); k++) {
                if(samples_names[f][k].find("QCD")==0) {
                  if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
                  if ((((c+1)&muonMask)!=0) && samples_names[NbOfCategories][k].find("QCD_mu")!=0) continue ;
                  if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
                  if ((((c+1)&electronMask)!=0) && samples_names[NbOfCategories][k].find("QCD_el")!=0) continue ;
                  if (hNbtaggedJets_qcd[chan][i][j]==NULL) {
                    hNbtaggedJets_qcd[chan][i][j] = (TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                    hNbtaggedJets_qcd[chan][i][j]->SetDirectory(NULL);
                    hNbtaggedJets_qcd[chan][i][j]->Reset();
                  }
                  hNbtaggedJets_qcd[chan][i][j]->Add((TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[f][k]);
                }
              }

              f = NbOfCategories + 1 ;
              for (UInt_t k=0; k<files[f].size(); k++) {
                if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
                if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
                if (hNbtaggedJets_rare[chan][i][j]==NULL) {
                  hNbtaggedJets_rare[chan][i][j] = (TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                  hNbtaggedJets_rare[chan][i][j]->SetDirectory(NULL);
                  hNbtaggedJets_rare[chan][i][j]->Reset();
                }
                hNbtaggedJets_rare[chan][i][j]->Add((TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[f][k]);
              }

              f = NbOfCategories + 2 ;
              for (UInt_t k=0; k<files[f].size(); k++) {
                if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
                if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
                if (hNbtaggedJets_signal[chan][i][j]==NULL) {
                  hNbtaggedJets_signal[chan][i][j] = (TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                  hNbtaggedJets_signal[chan][i][j]->SetDirectory(NULL);
                  hNbtaggedJets_signal[chan][i][j]->Reset();
                }
                hNbtaggedJets_signal[chan][i][j]->Add((TH1D*) files[f][k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), IntLumi[c]*samples_weights[f][k]);
              }
            }
            
            // QCD shape from data isolation side-band
            if (estimateQCDApart == kTRUE)    {
              if (!semiMuon && (((c+1)&muonMask)!=0)) break ;
              //if ((((c+1)&muonMask)!=0) && samples_names[NbOfCategories][k].find("QCD_mu")!=0) continue ;
              if (!semiElectron && (((c+1)&electronMask)!=0)) break ;
              //if ((((c+1)&electronMask)!=0) && samples_names[NbOfCategories][k].find("QCD_el")!=0) continue ;
              if (hNbtaggedJets_qcd_shape[chan][i][j]==NULL) {
                //                hNbtaggedJets_qcd_shape[chan][i][j] = (TH1D*) f->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+wpSuffix[j]+btagJetMult_suffix).c_str())/*->Clone()*/;
                hNbtaggedJets_qcd_shape[chan][i][j] = (TH1D*) fi->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+"_4j_1b_iso_0p3_0p5"/*+jetSuffix[i]+wpSuffix[j]+btagJetMult_suffix*/).c_str())/*->Clone()*/;
                hNbtaggedJets_qcd_shape[chan][i][j]->SetDirectory(NULL);
                hNbtaggedJets_qcd_shape[chan][i][j]->Reset();
              }
              //    hNbtaggedJets_qcd_shape[chan][i][j]->Add(hNbtaggedJets_qcd[chan][i][j], 1.);
              //              hNbtaggedJets_qcd_shape[chan][i][j]->Print((std::string()+"echo"+channelLeg3Suffix[c]+jetSuffix[i]+wpSuffix[j]+btagJetMult_suffix+".pdf").c_str());
              hNbtaggedJets_qcd_shape[chan][i][j]->Add((TH1D*) fi->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+"_4j_1b_iso_0p3_0p5"/*+jetSuffix[i]+wpSuffix[j]+btagJetMult_suffix*/).c_str())/*, IntLumi[c]*samples_weights[f][k]*/);
              //              hNbtaggedJets_qcd_shape[chan][i][j]->Add((TH1D*) fi->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+"_4j_1b_iso_0p5_0p7"/*+jetSuffix[i]+wpSuffix[j]+btagJetMult_suffix*/).c_str())/*, IntLumi[c]*samples_weights[f][k]*/);
              //              hNbtaggedJets_qcd_shape[chan][i][j]->Add((TH1D*) fi->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+"_4j_1b_iso_0p7_1"/*+jetSuffix[i]+wpSuffix[j]+btagJetMult_suffix*/).c_str())/*, IntLumi[c]*samples_weights[f][k]*/);

            }
            
            
            
            /*              
             hNbtaggedJets_ttjets[chan][i][j]->Scale(IntLumi[c]);
             hNbtaggedJets_wjets[chan][i][j]->Scale(IntLumi[c]);
             hNbtaggedJets_zjets[chan][i][j]->Scale(IntLumi[c]);
             hNbtaggedJets_stjets[chan][i][j]->Scale(IntLumi[c]);
             hNbtaggedJets_vbjets[chan][i][j]->Scale(IntLumi[c]);
             hWbb_In_WJets[chan][i][j]->Scale(IntLumi[c]);
             hWcx_In_WJets[chan][i][j]->Scale(IntLumi[c]);
             hWLFg_In_WJets[chan][i][j]->Scale(IntLumi[c]);
             */
            
            
              //Data
            if (datafiles.size() != 0) {
              if (isOnData) {
                printf(" data");
                if (hNbtaggedJets[chan][i][j]==NULL) {
                  hNbtaggedJets[chan][i][j] = (TH1D*) datafiles[0]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str());
                  hNbtaggedJets[chan][i][j]->SetDirectory(NULL);
                  hNbtaggedJets[chan][i][j]->Reset();
                }
                for (UInt_t k=0; k<datafiles.size(); k++) {
                  printf("Testing for c=%d , k=%d\n", c, k);
                  if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
                  if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
                  if ((((c+1)&muonMask)!=0) && (dataNames[k].find("MuHad")!=0 && dataNames[k].find("Single Mu")!=0)) continue;
                  if ((((c+1)&electronMask)!=0) && (dataNames[k].find("ElHad")!=0 && dataNames[k].find("Single El")!=0)) continue;
                  hNbtaggedJets[chan][i][j]->Add((TH1D*) datafiles[k]->Get((std::string()+(isNewNaming==kTRUE ? btagJetMult_prefix : "B_Jet_Multiplicity")+channelLeg3Suffix[c]+jetSuffix[i]+nbtagSuffix[ib]+wpSuffix[j]+btagJetMult_suffix).c_str()), weight_data[k]);
                   printf("\ndata[%d][%d][%d] = %lf\n", chan, i, j, hNbtaggedJets[chan][i][j]->Integral(0,-1));
                   hNbtaggedJets[chan][i][j]->Print();
                   
                }
              }
            }
            printf("\n");
            
            }
          }
        }
      }
    }
    if (fi != NULL) {
      fi->Close();
    }
  }
  
  
  {
    printf("Scale factors and efficiencies stuff\n");
    
    /*
     init_Euds[0] = 0.15;
     init_Euds[1] = 0.05;
     init_Eudsc[0] = 0.2;
     init_Eudsc[1] = 0.05;
     for (int wp=0; wp<NbOfWP; wp++) {
     Euds_const_mean[wp] = init_Euds[wp];
     Euds_const_error[wp] = init_Euds[wp] * 0.01;
     Euds_const_error[wp] *= 20.;
     Eudsc_const_mean[wp] = init_Eudsc[wp] ;
     Eudsc_const_error[wp] = Eudsc_const_mean[wp] * .15    * 3;
     init_Eb[wp] = Eb_const_mean[wp] ;
     }
     */
    for (UInt_t c=0; c<NbOfChannels; c++) {
        //Begin of inserted
        // Applying the scale factor if from MC
        // Applying the value if from Data
      /**
       Results from BTV-11-003, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG link to https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/eff_b_c-ttbar_payload.txt 
       In Range : 0<=|eta|<=2.4 , 30<=p_T<=200  (approximately)
       */
      for (UInt_t j=0; j<NbOfWP; j++) {
        /*
         //"Should be" values from 7 TeV
        Eb_const_mean[c][j] = (c==0? 0.66675 : (c==1? 0.68411 : 0.67908)); //0.7
        Eb_const_error[c][j] = Eb_const_mean[c][j]*0.05*3;
        if (isOnData) {
          Eb_const_mean[c][j] *= (c==0? (0.62365/0.66675) : (c==1? (0.57379/0.68411) : (0.60889/0.67908)));//0.99
          Eb_const_error[c][j] *= 1.;
        }
         */
        /*
          //"Should be" values, Michael's analysis, 8 TeV
        Eb_const_mean[c][j] = (c==0? 0.698 : (c==1? 0.696 : 0.692)); //0.7
        Eb_const_error[c][j] = (c==0? 0.002 : (c==1? 0.002 : 0.002));
        if (isOnData) {
          Eb_const_mean[c][j] = (c==0? 0.696 : (c==1? 0.688 : 0.692)); //0.7
          Eb_const_error[c][j] = (c==0? 0.019 : (c==1? 0.028 : 0.028));
        }
        */
        /**/
       Eb_const_mean[c][j] = effTT[2][c];
        Eb_const_error[c][j] = errTT[2][c];
        /**/
      }
      for (UInt_t j=0; j<NbOfWP; j++) {
        
        /*
         //"Should be" values from 7 TeV
        Euds_const_mean[c][j]   = (c==0? 0.034026 : (c==1? 0.040233 : 0.039805));//0.03
        Euds_const_error[c][j]  = Euds_const_mean[c][j] * 0.15 *6 ;
        Eudsc_const_mean[c][j]  = (c==0? 0.046730 : (c==1? 0.045695 : 0.046730));//0.04
        Eudsc_const_error[c][j] = Eudsc_const_mean[c][j] * 0.15 *6 ;
        
        if (isOnData) {
          Euds_const_mean[c][j]   *= (c==0? 0.029720/0.034026 : (c==1? 0.032032/0.040233 : 0.030064/0.039805)) ;//0.99
          Euds_const_error[c][j]  *= 1. ;
          Eudsc_const_mean[c][j]  *= (c==0? 0.053578/0.046730 : (c==1? 0.051578/0.045695 : 0.053589/0.046730)) ;//0.99
          Eudsc_const_error[c][j] *= 1. ;
        }
         */
        /*
         //"Should be" values, Michael's analysis, 8 TeV
        Euds_const_mean[c][j]   = (c==0? 0.030 : (c==1? 0.029 : 0.0295)) ;//0.03
        Euds_const_error[c][j]  = (c==0? 0.003 : (c==1? 0.003 : 0.003)) ;
        Eudsc_const_mean[c][j]  = (c==0? 0.030 : (c==1? 0.029 : 0.0295)) ;//0.04
        Eudsc_const_error[c][j] = (c==0? 0.003 : (c==1? 0.003 : 0.003)) ;
        
        if (isOnData) {
          Euds_const_mean[c][j]   = (c==0? 0.032 : (c==1? 0.037 : 0.0345)) ;//0.99
          Euds_const_error[c][j]  = (c==0? 0.004 : (c==1? 0.004 : 0.004)) ;
          Eudsc_const_mean[c][j]  = (c==0? 0.032 : (c==1? 0.037 : 0.0345)) ;//0.99
          Eudsc_const_error[c][j] = (c==0? 0.004 : (c==1? 0.004 : 0.004)) ;
        }
        /*
        /**/
        Euds_const_mean[c][j] = effV[1][c]; //udsc on V
        Euds_const_error[c][j] = errV[1][c]; //udsc on V
        /**/
        Eudsc_const_mean[c][j] = effTT[1][c]; //udsc on TT
        Eudsc_const_error[c][j] = errTT[1][c]; //udsc on TT
        /**/
        init_Eb[c][j] = Eb_const_mean[c][j];
        init_Euds[c][j] = Euds_const_mean[c][j];
        init_Eudsc[c][j] = Eudsc_const_mean[c][j];
        
        
          //end of inserted
        /*
        Euds_const_error[c][j] = WilsonScoreIntervalHigh(nEffUds_pass[c][j]*IntLumi[c],nEffUds_tot[c][j]*IntLumi[c]) - init_Euds[c][j];
        if ( init_Euds[c][j] - WilsonScoreIntervalLow(nEffUds_pass[c][j]*IntLumi[c],nEffUds_tot[c][j]*IntLumi[c]) > Euds_const_error[c][j] )
          Euds_const_error[c][j] = init_Euds[c][j] - WilsonScoreIntervalLow(nEffUds_pass[c][j]*IntLumi[c],nEffUds_tot[c][j]*IntLumi[c]);
        Euds_const_error[c][j] = Euds_const_mean[c][j] * 0.15;
        
        
        Eudsc_const_mean[c][j] = init_Eudsc[c][j] ;
        Eudsc_const_error[c][j] = WilsonScoreIntervalHigh(nEffUdsc_pass[c][j]*IntLumi[c],nEffUdsc_tot[c][j]*IntLumi[c]) - init_Eudsc[c][j];
        if ( init_Eudsc[c][j] - WilsonScoreIntervalLow(nEffUdsc_pass[c][j],nEffUdsc_tot[c][j]) > Eudsc_const_error[c][j] )
          Eudsc_const_error[c][j] = init_Eudsc[c][j] - WilsonScoreIntervalLow(nEffUdsc_pass[c][j],nEffUdsc_tot[c][j]);
        Eudsc_const_error[c][j] = Eudsc_const_mean[c][j] * .15 ;
        */
        init_Eb[c][j] = Eb_const_mean[c][j] ;
        
      }
    }
  }
  
  
  
  
  std::vector<int>::iterator Idx;
  std::vector<int> FixedVarIdx;
  //FixedVarIdx.push_back(0); // Fix eb, b-tagging efficiency
  //FixedVarIdx.push_back(1); // Fix eudsc, mis-tagging efficiency for tt-like events
    //FixedVarIdx.push_back(2); // Fix euds, mis-tagging efficiency for v-like events
  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////// Analysis ////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    // E v e n t   s e l e c t i o n
    // -----------------------------
  
  float **Nttlike = new float*[NbOfChannels];
  double *NttlikeTot = new double[NbOfChannels];
  float **Nvlike = new float*[NbOfChannels];
  double *NvlikeTot = new double[NbOfChannels];
  for (UInt_t c=0; c<NbOfChannels; c++) {
    Nttlike[c] = new float[NbOfJetBins];
    NttlikeTot[c] = 0.;
    Nvlike[c] = new float[NbOfJetBins];
    NvlikeTot[c] = 0.;
  }
  
  printf("Ready to go through all the events (if not reloading histograms)\n");
  
    //cout<<" - data : Muon = "<<
    //        anaMu_data->Show(0);
    //<<endl;
  
  if (isOnData) {
    /*
     std::vector<double> discri = *(anaMu_data->vjtDiscri_pf);
     std::vector<double> pt = *(anaMu_data->vjtPt_pf);
     std::vector<double> eta = *(anaMu_data->vjtEta_pf);
     */
    /*
     if (btagCuts[wp] > 1.7) {
     for (Int_t j=1; j<mistag_eta_binning_medium_NUM; j++) {
     if (mistag_eta_binning_medium[j-1] <=eta[i] && eta[i]<=mistag_eta_binning_medium[j]) {
     Float_t localMistag = mistag_light_mean[j-1]->Eval(pt[i]);
     sum_jets_mistagLight[c][wp] += localMistag;
     number_jets_mistagLight[c][wp] += 1.;
     sum_events_mistagLight[c][wp] += localMistag/discri.size();
     number_events_mistagLight[c][wp] += 1./discri.size();
     Float_t localSF = sf_light_mean[j-1]->Eval(pt[i]);
     sum_jets_SFLight[c][wp] += localSF;
     number_jets_SFLight[c][wp] += 1.;
     sum_events_SFLight[c][wp] += localSF/discri.size();
     number_events_SFLight[c][wp] += 1./discri.size();
     break;
     }
     }
     } else {
     for (Int_t j=1; j<mistag_eta_binning_loose_NUM; j++) {
     if (mistag_eta_binning_loose[j-1] <=eta[i] && eta[i]<=mistag_eta_binning_loose[j]) {
     Float_t localMistag = mistag_light_mean_WPloose[j-1]->Eval(pt[i]);
     sum_jets_mistagLight[c][wp] += localMistag;
     number_jets_mistagLight[c][wp] += 1.;
     sum_events_mistagLight[c][wp] += localMistag/discri.size();
     number_events_mistagLight[c][wp] += 1./discri.size();
     Float_t localSF = sf_light_mean_WPloose[j-1]->Eval(pt[i]);
     sum_jets_SFLight[c][wp] += localSF;
     number_jets_SFLight[c][wp] += 1.;
     sum_events_SFLight[c][wp] += localSF/discri.size();
     number_events_SFLight[c][wp] += 1./discri.size();
     break;
     }
     }
     }
     */
  } else {
    
      // For tt+jets events 
    /*
     double selWeight= (IntLumi[muonMask-1] * preselEff_ttjets_NUM[muonMask-1]) / (preselEff_ttjets_DEN[muonMask-1]*n_Mu);
     Double_t MCweight = anaMu_ttjets->MCweight;
     std::vector<double> discri = *(anaMu_ttjets->vjtDiscri_pf);
     std::vector<double> flav = *(anaMu_ttjets->vjtFlav_pf);
     std::vector<double> pt = *(anaMu_ttjets->vjtPt_pf);
     std::vector<double> eta = *(anaMu_ttjets->vjtEta_pf);
     for (int c=0; c<NbOfChannels; c++){
     if ( ((c+1) &muonMask) == 0)
     continue;
     nB_Flavoured_Jets = 0;
     for (int wp=0; wp<NbOfWP; wp++) {
     int NbtaggedJets = 0;
     for(UInt_t i = 0; i<discri.size(); ++i) {
     if (std::fabs(flav[i])==5.) {
     if (wp==0) {
     nB_Flavoured_Jets++;
     }
     if(discri[i]>btagCuts[wp]) {
     nEffb_tt_pass[c][wp] += MCweight*lumiWeight;
     }
     nEffb_tt_tot[c][wp] += MCweight*lumiWeight;
     }
     if(pt.at(i)<30.) { printf("Jet looser for tt_jets\n"); }
     if(discri[i]>btagCuts[wp]) NbtaggedJets++;
     if (discri[i]>btagCuts[wp] && std::fabs(flav[i])!=5.) {
     nEffUdsc_pass[c][wp]+=MCweight*lumiWeight;
     }
     if (std::fabs(flav[i])!=5.) {
     nEffUdsc_tot[c][wp]+=MCweight*lumiWeight;
     }
     
     }
     hNbtaggedJets_ttjets[c][njt_idx][wp]->Fill(NbtaggedJets, MCweight * lumiWeight * selWeight);
     }
     histoNJ_ttjets[c]->Fill(njt, MCweight*lumiWeight * selWeight);
     
     //Njt_pf
     switch (nB_Flavoured_Jets) {
     case 0:
     case 1:
     case 2:
     eXbq__[c][njt_idx]->Fill(nB_Flavoured_Jets, MCweight*lumiWeight * selWeight);
     break;
     default:
     eXbq__[c][njt_idx]->Fill(3, MCweight*lumiWeight * selWeight);
     break;
     }
     }
     }
     */      
    
    /*
     for (int c=0; c<NbOfChannels; c++) {
     printf("\nFor channel #%d\n", c);
     for (int j=0; j<NbOfWP; j++) {
     printf("----->>>>>> b-tagging eff in tt+jets sample WP %d : %lf up:%lf down:%lf = %lf / %lf\n", j, nEffb_tt_pass[c][j] / nEffb_tt_tot[c][j],
     WilsonScoreIntervalHigh(nEffb_tt_pass[c][j],nEffb_tt_tot[c][j]),
     WilsonScoreIntervalLow(nEffb_tt_pass[c][j],nEffb_tt_tot[c][j]),
     nEffb_tt_pass[c][j], nEffb_tt_tot[c][j]);
     printf("----->>>>>> Mistag rate in tt+jets sample WP %d : %lf mean %lf up %lf down %lf = %lf / %lf\n", j, nEffUdsc_pass[c][j] / nEffUdsc_tot[c][j],
     WilsonScoreIntervalMean(nEffUdsc_pass[c][j],nEffUdsc_tot[c][j]),
     WilsonScoreIntervalHigh(nEffUdsc_pass[c][j],nEffUdsc_tot[c][j]),
     WilsonScoreIntervalLow(nEffUdsc_pass[c][j],nEffUdsc_tot[c][j]),
     nEffUdsc_pass[c][j], nEffUdsc_tot[c][j]);
     
     
     for (int i=0; i<NbOfJetBins; i++) {
     hNbtaggedJets_ttjets[c][i][j]->Scale(XS_ttjets    ); // Selection, on skimmed dataset, to correct the number of entries of the tree
     }
     }
     }
     cout<< endl;
     }
     */
    
      // For w+jets events 
    /*
     for (int c=0; c<NbOfChannels; c++) {
     for (int j=0; j<NbOfWP; j++) {
     nEffUds_pass_tmp[c][j] = 0.;
     nEffUds_tot_tmp[c][j] = 0.;
     }
     }
     if (semiMuon && triggerMu && anaMu_wjets->CutMu_iso0p1_3jets(-1)==0) {
     lumiWeight = numInteractions->GetBinContent(numInteractions->FindBin(anaMu_wjets->truenuminter));
     double selWeight= (IntLumi[muonMask-1] * preselEff_wjets_NUM[muonMask-1]) / (preselEff_wjets_DEN[muonMask-1]*n_Mu);
     for (int c=0; c<NbOfChannels; c++){
     if ( ((c+1) &muonMask) == 0)
     continue;
     for (int wp=0; wp<NbOfWP; wp++) {
     nB_Flavoured_Jets = 0;
     nC_Flavoured_Jets = 0;
     nLightAndGluon_Flavoured_Jets = 0;
     
     NbtaggedJets = 0;
     for(UInt_t i = 0; i<anaMu_wjets->vjtDiscri_pf->size(); ++i)
     {
     if(anaMu_wjets->vjtPt_pf->at(i)<30.) { printf("Jet looser for w_jets\n"); }
     if((*(anaMu_wjets->vjtDiscri_pf))[i]>btagCuts[wp]) NbtaggedJets++;
     if (std::fabs((*(anaMu_wjets->vjtFlav_pf))[i])==5.) {
     nB_Flavoured_Jets++;
     } else if (std::fabs((*(anaMu_wjets->vjtFlav_pf))[i])==4.) {
     nC_Flavoured_Jets++;
     } else {
     nLightAndGluon_Flavoured_Jets++;
     }
     }
     if (nB_Flavoured_Jets == 0) {
     for(UInt_t i = 0; i<anaMu_wjets->vjtDiscri_pf->size(); ++i)
     {
     if ((*(anaMu_wjets->vjtDiscri_pf))[i]>btagCuts[wp]) {
     nEffUds_pass_tmp[c][wp]+=anaMu_wjets->MCweight*lumiWeight;
     }
     nEffUds_tot_tmp[c][wp]+=anaMu_wjets->MCweight*lumiWeight;
     }
     } else {
     for(UInt_t i = 0; i<anaMu_wjets->vjtDiscri_pf->size(); ++i)
     {
     if (std::fabs((*(anaMu_wjets->vjtFlav_pf))[i])==5.) {
     if((*(anaMu_wjets->vjtDiscri_pf))[i]>btagCuts[wp]) {
     nEffb_w_pass[c][wp] += anaMu_wjets->MCweight*lumiWeight;
     }
     nEffb_w_tot[c][wp] += anaMu_wjets->MCweight*lumiWeight;
     }
     }
     }
     hNbtaggedJets_wjets[c][njt_idx][wp]->Fill(NbtaggedJets, anaMu_wjets->MCweight*lumiWeight * selWeight);
     if (nB_Flavoured_Jets>0) {
     hWbb_In_WJets[c][njt_idx][wp]->Fill(NbtaggedJets, anaMu_wjets->MCweight*lumiWeight * selWeight);
     } else if (nC_Flavoured_Jets>0) {
     hWcx_In_WJets[c][njt_idx][wp]->Fill(NbtaggedJets, anaMu_wjets->MCweight*lumiWeight * selWeight);
     } else {
     hWLFg_In_WJets[c][njt_idx][wp]->Fill(NbtaggedJets, anaMu_wjets->MCweight*lumiWeight * selWeight);
     }
     }
     */
    
    /*
     for (int c=0; c<NbOfChannels ; c++) {
     printf("\nFor channel #%d\n", c);
     for (int j=0; j<NbOfWP; j++) {
     printf("----->>>>>> b-tagging eff in W+jets sample WP %d : %lf up:%lf down:%lf = %lf / %lf\n", j, nEffb_w_pass[c][j] / nEffb_w_tot[c][j],
     WilsonScoreIntervalHigh(nEffb_w_pass[c][j],nEffb_w_tot[c][j]),
     WilsonScoreIntervalLow(nEffb_w_pass[c][j],nEffb_w_tot[c][j]),
     nEffb_w_pass[c][j], nEffb_w_tot[c][j]);
     printf("----->>>>>> Mistag rate in W+jets sample WP %d : %lf mean %lf up %lf down %lf = %lf / %lf\n",
     j, nEffUds_pass_tmp[c][j] / nEffUds_tot_tmp[c][j],
     WilsonScoreIntervalMean(nEffUds_pass_tmp[c][j],nEffUds_tot_tmp[c][j]),
     WilsonScoreIntervalHigh(nEffUds_pass_tmp[c][j],nEffUds_tot_tmp[c][j]),
     WilsonScoreIntervalLow(nEffUds_pass_tmp[c][j],nEffUds_tot_tmp[c][j]),
     nEffUds_pass_tmp[c][j], nEffUds_tot_tmp[c][j]);
     
     nEffUds_pass[c][j] += nEffUds_pass_tmp[c][j];
     nEffUds_tot[c][j]  += nEffUds_tot_tmp[c][j];
     
     for (int i=0; i<NbOfJetBins; i++) {
     hNbtaggedJets_wjets[c][i][j]->Scale(XS_wjets /*NormFact_wjets*IntLumi[c] * / ); // Selection, on skimmed dataset, to correct the number of entries of the tree
     }
     }
     }
     cout<< endl;
     }
     */
    
    /*
     for (int c=0; c<NbOfChannels ; c++) {
     for (int j=0; j<NbOfWP; j++) {
     for (int i=0; i<NbOfJetBins; i++) {
     Double_t normChangeWbb = 1.; //If from W+jets inclusive
     if (isWbbFromWbbMeas) {
     //      normChange = NormFact_vbjets/NormFact_wjets;
     //      (XS_vbjets *IntLumi[c]) / hWbb_In_WJets->GetIntegral();
     normChangeWbb = (XS_vbjets) / hWbb_In_WJets[c][i][j]->Integral();
     }
     if (reloadHisto) {
     hWbb_In_WJets[c][i][j]->Scale(1.*normChangeWbb);
     hWcx_In_WJets[c][i][j]->Scale(1.);
     hWLFg_In_WJets[c][i][j]->Scale(IntLumi[c]);
     } else {
     hWbb_In_WJets[c][i][j]->Scale(XS_wjets /* NormFact_wjets*IntLumi[c] * / *normChangeWbb);
     hWcx_In_WJets[c][i][j]->Scale(XS_wjets /* NormFact_wjets*IntLumi[c] * / );
     hWLFg_In_WJets[c][i][j]->Scale(XS_wjets /* NormFact_wjets*IntLumi[c] * / );
     }
     }
     }
     }
     */
    
    
      // For z+jets events 
    
    /*
     NbtaggedJets = 0;
     nB_Flavoured_Jets = 0;
     for(UInt_t i = 0; i<anaMu_zjets->vjtDiscri_pf->size(); ++i) {
     if(anaMu_zjets->vjtPt_pf->at(i)<30.) { printf("Jet looser for z_jets\n"); }
     if((*(anaMu_zjets->vjtDiscri_pf))[i]>btagCuts[wp]) NbtaggedJets++;
     if (std::fabs((*(anaMu_zjets->vjtFlav_pf))[i])==5.) {
     nB_Flavoured_Jets++;
     }
     }
     if (nB_Flavoured_Jets == 0) {
     for(UInt_t i = 0; i<anaMu_zjets->vjtDiscri_pf->size(); ++i) {
     if ((*(anaMu_zjets->vjtDiscri_pf))[i]>btagCuts[wp]) {
     nEffUds_pass_tmp[c][wp]+=anaMu_zjets->MCweight*lumiWeight;
     }
     nEffUds_tot_tmp[c][wp]+=anaMu_zjets->MCweight*lumiWeight;
     }
     }
     
     */
    
    /*
     for (int c=0; c<NbOfChannels; c++) {
     printf("\nFor channel #%d\n", c);
     for (int j=0; j<NbOfWP; j++) {
     printf("\n\n----->>>>>> Mistag rate in Z+jets sample WP %d : %lf mean %lf up %lf down %lf = %lf / %lf\n", j,
     nEffUds_pass_tmp[c][j] / nEffUds_tot_tmp[c][j],
     WilsonScoreIntervalMean(nEffUds_pass_tmp[c][j],nEffUds_tot_tmp[c][j]),
     WilsonScoreIntervalHigh(nEffUds_pass_tmp[c][j],nEffUds_tot_tmp[c][j]),
     WilsonScoreIntervalLow(nEffUds_pass_tmp[c][j],nEffUds_tot_tmp[c][j]),
     nEffUds_pass_tmp[j], nEffUds_tot_tmp[c][j]);
     nEffUds_pass[c][j] += nEffUds_pass_tmp[c][j];
     nEffUds_tot[c][j]  += nEffUds_tot_tmp[c][j];
     }
     for (int i=0; i<NbOfJetBins; i++) {
     for (int wp=0; wp<NbOfWP; wp++) {
     hNbtaggedJets_zjets[c][i][wp]->Scale(XS_zjets/*NormFact_zjets*IntLumi[c] * /    ); // Selection, on skimmed dataset, to correct the number of entries of the tree);
     }
     }
     }
     cout <<endl;
     }
     */
      // For st+jets events 
    /*
     printf("\n\n\n");
     for (int c=0; c<NbOfChannels; c++) {
     printf("\nFor channel #%d\n", c);
     for (int wp=0; wp<NbOfWP; wp++) {
     printf("----->>>>>> Mistag rate in V+jets sample (with NormFact) WP %d : %lf mean %lf up %lf down %lf = %lf / %lf\n", wp,
     nEffUds_pass[c][wp] / nEffUds_tot[c][wp],
     WilsonScoreIntervalMean(nEffUds_pass[c][wp]*IntLumi[c],nEffUds_tot[c][wp]*IntLumi[c]),
     WilsonScoreIntervalHigh(nEffUds_pass[c][wp]*IntLumi[c],nEffUds_tot[c][wp]*IntLumi[c]),
     WilsonScoreIntervalLow(nEffUds_pass[c][wp]*IntLumi[c],nEffUds_tot[c][wp]*IntLumi[c]),
     nEffUds_pass[c][wp]*IntLumi[c], nEffUds_tot[c][wp]*IntLumi[c]);
     }
     }
     */
  }
  
  cout << endl << endl << endl << endl << endl << endl << "********************************************" << endl;
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to perform the crawling in the event trees" << endl;
  cout << "********************************************" << endl;
  cout << " Starting the fit part now          " << endl;
  cout << "********************************************" << endl << endl << endl << endl << endl;
  
  
  
    //  pileUpRewightingFile->close();
  
  /** Njtbtag CANNOT BE USED AS IT INCLUDES NON-SELECTED JETS !!!!
   
   // Fill 2D histograms Nb b-tagged jets vs Nb jets
   ttjets->Draw("Njtbtag_pf>>hNbtaggedJets_ttjets(10,0,10)",cuts);
   wjets->Draw("Njtbtag_pf>>hNbtaggedJets_wjets(10,0,10)",cuts);
   zjets->Draw("Njtbtag_pf>>hNbtaggedJets_zjets(10,0,10)",cuts);
   stjets->Draw("Njtbtag_pf>>hNbtaggedJets_stjets(10,0,10)",cuts);
   //  vbjets->Draw("Njtbtag_pf>>hNbtaggedJets_vbjets(10,0,10)",cuts);
   
   // Retrieve the b-tagged jet multiplicity histograms
   hNbtaggedJets_ttjets = (TH1F*)gDirectory->Get("hNbtaggedJets_ttjets");
   hNbtaggedJets_wjets  = (TH1F*)gDirectory->Get("hNbtaggedJets_wjets");
   hNbtaggedJets_zjets  = (TH1F*)gDirectory->Get("hNbtaggedJets_zjets");
   hNbtaggedJets_stjets = (TH1F*)gDirectory->Get("hNbtaggedJets_stjets");
   //  hNbtaggedJets_vbjets = (TH1F*)gDirectory->Get("hNbtaggedJets_vbjets");
   
   // Set title axis for histograms
   hNbtaggedJets_ttjets->SetTitle(";B-tagged jet multiplicity");
   hNbtaggedJets_wjets->SetTitle(";B-tagged jet multiplicity");
   hNbtaggedJets_zjets->SetTitle(";B-tagged jet multiplicity");
   hNbtaggedJets_stjets->SetTitle(";B-tagged jet multiplicity");
   //  hNbtaggedJets_vbjets->SetTitle(";B-tagged jet multiplicity");
   */
  
    // Select processes to be considered in the estimation
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      /*
       for (int wp=0; wp<NbOfWP; wp++) {
       // - tt-like processes
       hNbtaggedJets_ttlike[c][i][wp]->Add(hNbtaggedJets_ttjets[c][i][wp]);
       // - v-like processes
       //  hNbtaggedJets_vlike[i]->Add(hNbtaggedJets_wjets[i]);
       hNbtaggedJets_vlike[c][i][wp]->Add(hWcx_In_WJets[c][i][wp]);
       hNbtaggedJets_vlike[c][i][wp]->Add(hWLFg_In_WJets[c][i][wp]);
       hNbtaggedJets_vlike[c][i][wp]->Add(hNbtaggedJets_zjets[c][i][wp]);
       
       if (hNbtaggedJets_qcd[c][i][wp]!=NULL) {
       for(int k=0; k<hNbtaggedJets_qcd[c][i][wp]->GetNbinsX(); k++) {
       //                          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_qcd[c][i][wp]);
       hNbtaggedJets_vlike[c][i][wp]->Fill(hNbtaggedJets_qcd[c][i][wp]->GetBinLowEdge(k), hNbtaggedJets_qcd[c][i][wp]->GetBinContent(k)/*0.* /);
       }
       }
       }
       */
        // Set true numbers of tt-like and v-like events
      Nttlike[c][i] = hNbtaggedJets_ttlike[c][i][0]->Integral();
      NttlikeTot[c] += Nttlike[c][i];
      Nvlike[c][i]  = hNbtaggedJets_vlike[c][i][0]->Integral();
      NvlikeTot[c] += Nvlike[c][i];
    }
  }
  
  
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    if (!Wbb_from_WJets) {
      UInt_t Nbins = 10 ;
      TFile *f = new TFile(WbbShapeRootfile.c_str(), "READ");
      if (!f->IsZombie()) {
        char name[100];
        for (UInt_t j=0; j<12; j++) {
          if ((j!=1)&&(j!=2)&&(j!=5)) {
            continue;
          }
          for (UInt_t i=0; i<NbOfJetBins; i++) {
            for (UInt_t wp=0; wp<NbOfWP; wp++) {
              sprintf(name, (std::string()+"HistoryFlavor--FlavorHistoryPath_%d"+channelSuffix[c]+jetSuffix[i]+wpSuffix[wp]).c_str(), j);
              TH1I *h = (TH1I*) f->Get(name);
              hNbtaggedJets_vbjets[c][i][wp]->Reset();
              for (UInt_t k=1; k<=Nbins; k++) {
                hNbtaggedJets_vbjets[c][i][wp]->Fill(k-1,h->GetBinContent(k));
              }
            }
          }
        }
      }
      f->Close();
    }
  }
  
  
  
  
    // Set numbers of background events and associated uncertainties
  float** Nstjets = new float*[NbOfChannels];
  float** Nstjets_uncert = new float*[NbOfChannels];
  float** Nvbjets = new float*[NbOfChannels];
  float** Nvbjets_uncert = new float*[NbOfChannels];
  float** Nqcd_mean = new float*[NbOfChannels];
  float** Nqcd_uncert = new float*[NbOfChannels];
  float** Nvvjets = new float*[NbOfChannels];
  float** Nvvjets_uncert = new float*[NbOfChannels];
  double* NstjetsTot = new double[NbOfChannels];
  double* NvbjetsTot = new double[NbOfChannels];
  double* NmultijetsTot = new double[NbOfChannels];
  double* NvvjetsTot = new double[NbOfChannels];
  for (UInt_t c=0; c<NbOfChannels; c++) {
    Nstjets[c] = new float[NbOfJetBins];
    Nstjets_uncert[c] = new float[NbOfJetBins];
    Nvbjets[c] = new float[NbOfJetBins];
    Nvbjets_uncert[c] = new float[NbOfJetBins];
    Nqcd_mean[c] = new float[NbOfJetBins];
    Nqcd_uncert[c] = new float[NbOfJetBins];
    Nvvjets[c] = new float[NbOfJetBins];
    Nvvjets_uncert[c] = new float[NbOfJetBins];
    NstjetsTot[c] = 0.;
    NvbjetsTot[c] = 0.;
    NmultijetsTot[c] = 0.;
    NvvjetsTot[c] = 0.;
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      Nstjets[c][i] = hNbtaggedJets_stjets[c][i][0]->Integral();
      NstjetsTot[c] += Nstjets[c][i];
      Nstjets_uncert[c][i] = 0.3; // 30% uncertainty on the st+jets XS ** TO BE ADAPTED **
      Nvbjets[c][i] = hNbtaggedJets_vbjets[c][i][0]->Integral();
      NvbjetsTot[c] += Nvbjets[c][i];
      Nvbjets_uncert[c][i] = 0.3; // 30% uncertainty on the vb+jets XS ** TO BE ADAPTED **
      
      Nvvjets[c][i] = 0.;
      if (hNbtaggedJets_vvjets[c][i][0]!=NULL) {
        Nvvjets[c][i] = hNbtaggedJets_vvjets[c][i][0]->Integral()/*0.*/;
      }
      Nvvjets_uncert[c][i] = 0.30;
      NvvjetsTot[c] += Nvvjets[c][i];
      Nqcd_mean[c][i] = 0.;
      if (hNbtaggedJets_qcd[c][i][0]!=NULL) {
        if (estimateQCDApart) {
          Nqcd_mean[c][i] = hNbtaggedJets_qcd[c][i][0]->Integral()/*0.*/;
        } else {
          Double_t isoSidebandSF = (c+1==muonMask  ? 4. : 1. ) ;
          Nqcd_mean[c][i] = hNbtaggedJets_qcd[c][i][0]->Integral() * isoSidebandSF ;
        }
      }
      Nqcd_uncert[c][i] = 1. ;
      NmultijetsTot[c] += Nqcd_mean[c][i];
    }
  }
  
  
    // Sum over all datasets
  if (!isOnData) {
    for (UInt_t c=0; c<NbOfChannels; c++) {
      for (UInt_t wp=0; wp<NbOfWP; wp++) {
        for (UInt_t i=0; i<NbOfJetBins; i++) {
          if (hNbtaggedJets[c][i][wp]==NULL) {
            hNbtaggedJets[c][i][wp] = (TH1D*) hNbtaggedJets_ttlike[c][i][wp]->Clone() ;
            hNbtaggedJets[c][i][wp]->SetDirectory(NULL) ;
            hNbtaggedJets[c][i][wp]->Reset() ;
          }
          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_ttlike[c][i][wp]);
          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_vlike[c][i][wp]);
          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_stjets[c][i][wp]);
          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_vbjets[c][i][wp]);
          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_vvjets[c][i][wp]);
          
          
          if (estimateQCDApart) {
            hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_qcd[c][i][wp], 1.);
          }
          if (addRareProcesses == kTRUE) {
              hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_rare[c][i][wp]);
          }
          
          if (signalFile.compare("")!=0) {
            hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_signal[c][i][wp]);
          }

            //          if (hNbtaggedJets_qcd[c][i][wp]!=NULL) {
            //            for(int k=0; k<hNbtaggedJets_qcd[c][i][wp]->GetNbinsX(); k++) {
            ////                          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_qcd[c][i][wp]);
            //              hNbtaggedJets[c][i][wp]->Fill(hNbtaggedJets_qcd[c][i][wp]->GetBinLowEdge(k), hNbtaggedJets_qcd[c][i][wp]->GetBinContent(k)/*0.*/);
            //            }
            //          }
          
            //          if (hNbtaggedJets_vvjets[c][i][wp]!=NULL) {
            //            for(int k=0; k<hNbtaggedJets_vvjets[c][i][wp]->GetNbinsX(); k++) {
            //                //                            hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_vvjets[c][i][wp]);
            //              hNbtaggedJets[c][i][wp]->Fill(hNbtaggedJets_vvjets[c][i][wp]->GetBinLowEdge(k), hNbtaggedJets_vvjets[c][i][wp]->GetBinContent(k)/*0.*/);
            //            }
            //          }
        }
      }
    }
  }
  
    // Total nb of events (rescaled to the required int. lumi.)
  Int_t** nExpected = new Int_t*[NbOfChannels];
  Int_t* nTot = new Int_t[NbOfChannels];
  for (UInt_t c=0; c<NbOfChannels; c++) {
    nExpected[c] = new Int_t[NbOfJetBins];
    nTot[c] = 0;
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      nExpected[c][i] = 0;
      for(UInt_t j=0; j<NbOfWP; j++) {
        if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
        if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
        nExpected[c][i] += (int) hNbtaggedJets[c][i][0]->Integral();
      }
      nTot[c] += nExpected[c][i];
    }
    
    for(UInt_t i=0; i<NbOfJetBins; i++) {
        //      printf("eXbq__[%d][%d] %lu\n", c, i, (unsigned long) eXbq__[c][i]);
        //      printf("Rescale [%d][%d] : %lf\n", c, i, eXbq__[c][i]->Integral(0,-1));
      if (forceE3bq==kTRUE) {
        eXbq__[c][i]->SetBinContent(eXbq__[c][i]->FindBin(3.), ((eXbq__[c][i]->Integral(0,-1))-eXbq__[c][i]->GetBinContent(eXbq__[c][i]->FindBin(3.))) * fracForceE3bq/(1.-fracForceE3bq));
      }
      eXbq__[c][i]->Scale(1./eXbq__[c][i]->Integral(0,-1));
    }
  }
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
    for (UInt_t wp=0; wp<NbOfWP; wp++) {
      for (UInt_t i=0; i<NbOfJetBins; i++) {
        for (UInt_t k=1; k<=((UInt_t) hNbtaggedJets_ttlike[c][i][wp]->GetNbinsX()); k++) {
          hNbtaggedJets_ttlike[c][i][wp]->SetBinError(k, TMath::Sqrt(hNbtaggedJets_ttlike[c][i][wp]->GetBinContent(k)));
        }
        for (UInt_t k=1; k<=((UInt_t) hNbtaggedJets_vlike[c][i][wp]->GetNbinsX()); k++) {
          hNbtaggedJets_vlike[c][i][wp]->SetBinError(k, TMath::Sqrt(hNbtaggedJets_vlike[c][i][wp]->GetBinContent(k)));
        }
        for (UInt_t k=1; k<=((UInt_t) hNbtaggedJets_stjets[c][i][wp]->GetNbinsX()); k++) {
          hNbtaggedJets_stjets[c][i][wp]->SetBinError(k, TMath::Sqrt(hNbtaggedJets_stjets[c][i][wp]->GetBinContent(k)));
        }
        for (UInt_t k=1; k<=((UInt_t) hNbtaggedJets_vbjets[c][i][wp]->GetNbinsX()); k++) {
          hNbtaggedJets_vbjets[c][i][wp]->SetBinError(k, TMath::Sqrt(hNbtaggedJets_vbjets[c][i][wp]->GetBinContent(k)));
        }
        for (UInt_t k=1; k<=((UInt_t) hNbtaggedJets[c][i][wp]->GetNbinsX()); k++) {
          hNbtaggedJets[c][i][wp]->SetBinError(k, TMath::Sqrt(hNbtaggedJets[c][i][wp]->GetBinContent(k)));
        }
      }
    }
  }
  
  {
    for(UInt_t c=0; c<NbOfChannels; c++) { 
      if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
      if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
      printf("\nOn MC (before fitting if you are impatient): \n");
      printf("For channel : %s\n", channelSuffix[c]);
      printf("\\begin{tabular}{r|ccccc}\n");
      for (UInt_t j=0; j<NbOfWP; j++) {
        if (NbOfWP!=1) {
          printf("\\multicolumn{6}{l}{For working point %s :}\\\\ \n", wpSuffix[j]);
          printf("\\hline");
        }
        for (UInt_t i=0; i<NbOfJetBins; i++) {
          std::string jts = jetSuffix[i];
          if (jts.find("_")!=std::string::npos) {
            jts=jts.substr(jts.find("_")+1,std::string::npos);
          }
          if (jts.find("Inc")!=std::string::npos) {
            jts = " \\geq "+jts.substr(0,jts.find("Inc"))+jts.substr(jts.find("Inc")+3,std::string::npos);
          } else if (jts.find("Exc")!=std::string::npos) {
            jts = " = "+jts.substr(0,jts.find("Exc"))+jts.substr(jts.find("Exc")+3,std::string::npos);
          }
          
          printf("\\multicolumn{6}{l}{For jet bin $ %s $ :}\\\\ \n", jts.c_str());
          printf("$\\#b-tags$ & $t\\bar{t}-like$ & $V-like$ & Single Top & $W-bb$ & $VV$ \\\\ \\hline \n");
        for (UInt_t k=0; k<6; k++) { //because up to 5 b-tags in the equations
            printf("$ %d $ & $ %lf $ & $ %lf $ & $ %lf $ & $ %lf $ & $ %lf $ \\\\ \n", k,
                   hNbtaggedJets_ttlike[c][i][j]->GetBinContent( hNbtaggedJets_ttlike[c][i][j]->FindBin(k)),
                   hNbtaggedJets_vlike[c][i][j]->GetBinContent(hNbtaggedJets_vlike[c][i][j]->FindBin(k)),
                   hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k)),
                   hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k)),
                   hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k)));
          
          }
        }
      }
      printf("\\end{tabular}\n");
    }
  }
  
 

 {
    cout<<"***************************************"<<endl;
    cout<<"***************************************"<<endl;
    for (UInt_t c=0; c<NbOfChannels; c++) {
      if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
      if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
      cout<<"- Number of entries "<<endl;
      cout<<"-- for tt-like processes : "<<endl;
      for (UInt_t i=0; i<NbOfJetBins; i++) {
        cout<<"--- tt-like (jet bin "<< i << ") = "<<hNbtaggedJets_ttjets[c][i][0]->Integral()<<endl;
      }
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NttlikeTot[c]<<endl;
      cout<<"---------------------------------"<<endl;
      cout<<"-- for v-like processes : "<<endl;
      /*
       for (UInt_t i=0; i<NbOfJetBins; i++) {
       cout<<"--- W_LFg+jets (jet bin "<< i << ") = "<<hWLFg_In_WJets[c][i][0]->Integral()<<endl;
       cout<<"--- W_cx+jets (jet bin "<< i << ") = "<<hWcx_In_WJets[c][i][0]->Integral()<<endl;
       }
       for (int i=0; i<NbOfJetBins; i++) {
       cout<<"--- Z+jets (jet bin "<< i << ") = "<<hNbtaggedJets_zjets[c][i][0]->Integral()<<endl;
       }
       */
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NvlikeTot[c]<<endl;
      cout<<"---------------------------------"<<endl;
      cout<<"-- for background processes : "<<endl;
      for (UInt_t i=0; i<NbOfJetBins; i++) {
        cout<<"--- st+jets (jet bin "<< i << ") = "<<hNbtaggedJets_stjets[c][i][0]->Integral()<<endl;
      }
        //    cout<<"--- st+jets = "<<hNbtaggedJets_stjets->GetEntries()<<"/ "<<NormFact_stjets<<endl;
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NstjetsTot[c]<<endl;
      for (UInt_t i=0; i<NbOfJetBins; i++) {
        cout<<"--- vb+jets (jet bin "<< i << ") = "<<hNbtaggedJets_vbjets[c][i][0]->Integral()<<endl;//"/ "<<NormFact_vbjets<<endl;
      }
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NvbjetsTot[c]<<endl;
      for (UInt_t i=0; i<NbOfJetBins; i++) {
        if (hNbtaggedJets_qcd[c][i][0]!=NULL) {
          cout<<"--- multijets (jet bin "<< i << ") = "<<hNbtaggedJets_qcd[c][i][0]->Integral()<<endl;
        }
      }
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NmultijetsTot[c]<<endl;
      for (UInt_t i=0; i<NbOfJetBins; i++) {
        if (hNbtaggedJets_vvjets[c][i][0]!=NULL) {
          cout<<"--- VV+jets (jet bin "<< i << ") = "<<hNbtaggedJets_vvjets[c][i][0]->Integral()<<endl;
        }
      }
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NvvjetsTot[c]<<endl;
      cout<<"---------------------------------"<<endl;
      if (isOnData) {
        cout<<"-- for Data : "<<endl;
        double sum = 0.;
        double sumMC = 0.;
        for (UInt_t i=0; i<NbOfJetBins; i++) {
          sum += hNbtaggedJets[c][i][0]->Integral();
          double tmpMCsum = hNbtaggedJets_ttlike[c][i][0]->Integral() + hNbtaggedJets_vlike[c][i][0]->Integral() + hNbtaggedJets_stjets[c][i][0]->Integral() + hNbtaggedJets_vbjets[c][i][0]->Integral();
          sumMC += tmpMCsum;
          cout<<"--- total (jet bin "<< i << ") = "<<hNbtaggedJets[c][i][0]->Integral() << "         and for MC : " << tmpMCsum <<endl;
        }
        cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<< sum /*nTot[c]/NbOfWP*/ << "          and for MC : " << sumMC << endl;
      } else {
        cout<<"-- for all processes : "<<endl;
        for (UInt_t i=0; i<NbOfJetBins; i++) {
          cout<<"--- total (jet bin "<< i << ") = "<<hNbtaggedJets[c][i][0]->Integral()<<endl;
        }
        cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<nTot[c]/NbOfWP<<endl;
      }
      cout<<"***************************************"<<endl;
    }
    cout<<"***************************************"<<endl;
  }
  
  
    // C r e a t e   o b s e r v a b l e
    // ---------------------------------
  printf("Efficiencies\n");
  
  
  RooRealVar ***eb = new RooRealVar**[NbOfChannels];
  RooRealVar ***eudsc = new RooRealVar**[NbOfChannels];
  RooRealVar ***euds = new RooRealVar**[NbOfChannels];
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    eb[c] = new RooRealVar*[NbOfWP];
    eudsc[c] = new RooRealVar*[NbOfWP];
    euds[c] = new RooRealVar*[NbOfWP];
    
    for (UInt_t j=0; j<NbOfWP; j++) {
      eb[c][j] = new RooRealVar(    (std::string()+"eb"+channelSuffix[c]+wpSuf[j]).c_str(),  "#epsilon_{b-tag} ",      init_Eb[c][j],0.0,1.0);
      Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 0 );
      if(Idx != FixedVarIdx.end()) eb[c][j]->setConstant(kTRUE) ;
      eudsc[c][j] = new RooRealVar((std::string()+"eudsc"+channelSuffix[c]+wpSuf[j]).c_str(),"#epsilon_{mis-tag}",    init_Eudsc[c][j],0.0,1.0);
      Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 1 );
      if(Idx != FixedVarIdx.end()) eudsc[c][j]->setConstant(kTRUE) ;
      euds[c][j] = new RooRealVar( (std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str(),  "#epsilon^{,}_{mis-tag}",init_Euds[c][j],0.0,1.0);
      Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 2 );
      if(Idx != FixedVarIdx.end()) euds[c][j]->setConstant(kTRUE) ;
    }
  }
  
    // Estimation parameters
    //  RooRealVar Ntt("Ntt","N_{t#bar{t}-like}",init_Nttlike[JetIdx],0.0,init_Nttlike[JetIdx]*3);
    //  RooRealVar Nv("Nv","N_{V-like}",init_Nvlike[JetIdx],0.0,init_Nvlike[JetIdx]*4);
  
    //RooConstVar e3bq("e3bq","e3bq",e3bq_[njets]);
  RooRealVar ***Ntt = new RooRealVar**[NbOfChannels];
  RooRealVar ***Nv  = new RooRealVar**[NbOfChannels];
  RooRealVar ***Nst = new RooRealVar**[NbOfChannels];
  RooRealVar ***Nvb = new RooRealVar**[NbOfChannels];
  RooRealVar ***Nqcd = new RooRealVar**[NbOfChannels];
  RooRealVar ***Nvv = new RooRealVar**[NbOfChannels];
  
  RooConstVar ***n  = new RooConstVar**[NbOfChannels];
  
  RooConstVar ***e0bq = new RooConstVar**[NbOfChannels];
  RooConstVar ***e1bq = new RooConstVar**[NbOfChannels];
  RooConstVar ***e2bq = new RooConstVar**[NbOfChannels];
  RooConstVar ***e3bq = new RooConstVar**[NbOfChannels];

  RooConstVar**** n0bjets_st = new RooConstVar***[NbOfChannels];
  RooConstVar**** n1bjets_st = new RooConstVar***[NbOfChannels];
  RooConstVar**** n2bjets_st = new RooConstVar***[NbOfChannels];
  RooConstVar**** n3bjets_st = new RooConstVar***[NbOfChannels];
  RooConstVar**** n4bjets_st = new RooConstVar***[NbOfChannels];
  RooConstVar**** n5bjets_st = new RooConstVar***[NbOfChannels];
  RooConstVar**** n0bjets_vb = new RooConstVar***[NbOfChannels];
  RooConstVar**** n1bjets_vb = new RooConstVar***[NbOfChannels];
  RooConstVar**** n2bjets_vb = new RooConstVar***[NbOfChannels];
  RooConstVar**** n3bjets_vb = new RooConstVar***[NbOfChannels];
  RooConstVar**** n4bjets_vb = new RooConstVar***[NbOfChannels];
  RooConstVar**** n5bjets_vb = new RooConstVar***[NbOfChannels];
  RooConstVar**** n0bjets_qcd = new RooConstVar***[NbOfChannels];
  RooConstVar**** n1bjets_qcd = new RooConstVar***[NbOfChannels];
  RooConstVar**** n2bjets_qcd = new RooConstVar***[NbOfChannels];
  RooConstVar**** n3bjets_qcd = new RooConstVar***[NbOfChannels];
  RooConstVar**** n4bjets_qcd = new RooConstVar***[NbOfChannels];
  RooConstVar**** n5bjets_qcd = new RooConstVar***[NbOfChannels];
  RooConstVar**** n0bjets_vv = new RooConstVar***[NbOfChannels];
  RooConstVar**** n1bjets_vv = new RooConstVar***[NbOfChannels];
  RooConstVar**** n2bjets_vv = new RooConstVar***[NbOfChannels];
  RooConstVar**** n3bjets_vv = new RooConstVar***[NbOfChannels];
  RooConstVar**** n4bjets_vv = new RooConstVar***[NbOfChannels];
  RooConstVar**** n5bjets_vv = new RooConstVar***[NbOfChannels];
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    
    Ntt[c] = new RooRealVar*[NbOfJetBins];
    Nv[c]  = new RooRealVar*[NbOfJetBins];
    Nst[c]  = new RooRealVar*[NbOfJetBins];
    Nvb[c]  = new RooRealVar*[NbOfJetBins];
    Nqcd[c]  = new RooRealVar*[NbOfJetBins];
    Nvv[c]  = new RooRealVar*[NbOfJetBins];
    n[c] = new RooConstVar*[NbOfJetBins];
    e0bq[c] = new RooConstVar*[NbOfJetBins];
    e1bq[c] = new RooConstVar*[NbOfJetBins];
    e2bq[c] = new RooConstVar*[NbOfJetBins];
    e3bq[c] = new RooConstVar*[NbOfJetBins];
    n0bjets_st[c] = new RooConstVar**[NbOfJetBins]; n1bjets_st[c] = new RooConstVar**[NbOfJetBins]; n2bjets_st[c] = new RooConstVar**[NbOfJetBins];
    n3bjets_st[c] = new RooConstVar**[NbOfJetBins]; n4bjets_st[c] = new RooConstVar**[NbOfJetBins]; n5bjets_st[c] = new RooConstVar**[NbOfJetBins];
    n0bjets_vb[c] = new RooConstVar**[NbOfJetBins]; n1bjets_vb[c] = new RooConstVar**[NbOfJetBins]; n2bjets_vb[c] = new RooConstVar**[NbOfJetBins];
    n3bjets_vb[c] = new RooConstVar**[NbOfJetBins]; n4bjets_vb[c] = new RooConstVar**[NbOfJetBins]; n5bjets_vb[c] = new RooConstVar**[NbOfJetBins];
    n0bjets_qcd[c] = new RooConstVar**[NbOfJetBins]; n1bjets_qcd[c] = new RooConstVar**[NbOfJetBins]; n2bjets_qcd[c] = new RooConstVar**[NbOfJetBins];
    n3bjets_qcd[c] = new RooConstVar**[NbOfJetBins]; n4bjets_qcd[c] = new RooConstVar**[NbOfJetBins]; n5bjets_qcd[c] = new RooConstVar**[NbOfJetBins];
    n0bjets_vv[c] = new RooConstVar**[NbOfJetBins]; n1bjets_vv[c] = new RooConstVar**[NbOfJetBins]; n2bjets_vv[c] = new RooConstVar**[NbOfJetBins];
    n3bjets_vv[c] = new RooConstVar**[NbOfJetBins]; n4bjets_vv[c] = new RooConstVar**[NbOfJetBins]; n5bjets_vv[c] = new RooConstVar**[NbOfJetBins];
    
    printf("Other variables\n");
    
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      
      Ntt[c][i] = new RooRealVar((std::string()+"Ntt"+channelSuffix[c]+jtSuf[i]).c_str(),(std::string()+"N_{t#bar{t}-like}"+jetSuffix[i]).c_str(),Nttlike[c][i],0.0,Nttlike[c][i]*3);
      Nv[c][i] = new RooRealVar((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str(),(std::string()+"N_{V-like}"+jetSuffix[i]).c_str(),Nvlike[c][i],0.0,Nvlike[c][i]*4);
      Nst[c][i] = new RooRealVar((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str(),(std::string()+"N_{single top}"+jetSuffix[i]).c_str(),Nstjets[c][i],0.0,Nstjets[c][i]*4);
      Nvb[c][i] = new RooRealVar((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str(),(std::string()+"N_{V+b-jets}"+jetSuffix[i]).c_str(),Nvbjets[c][i],0.0,Nvbjets[c][i]*4);
      Nqcd[c][i] = new RooRealVar((std::string()+"Nqcd"+channelSuffix[c]+jtSuf[i]).c_str(),(std::string()+"N_{multijets}"+jetSuffix[i]).c_str(),Nqcd_mean[c][i],0.0,Nqcd_mean[c][i]*4/*+100000*/);
      Nvv[c][i] = new RooRealVar((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str(),(std::string()+"N_{VV}"+jetSuffix[i]).c_str(),Nvvjets[c][i],0.0,Nvvjets[c][i]*4/*+100000*/);
      
      n[c][i] = new RooConstVar((std::string()+"n"+channelSuffix[c]+jtSuf[i]).c_str(),"number of selected jets",MinNbOfJets+i) ;
      
      e0bq[c][i] = new RooConstVar((std::string()+"e0bq"+channelSuffix[c]+jtSuf[i]).c_str(),(std::string()+"e0bq"+jetSuffix[i]).c_str(),eXbq__[c][i]->GetBinContent(1));
      e1bq[c][i] = new RooConstVar((std::string()+"e1bq"+channelSuffix[c]+jtSuf[i]).c_str(),(std::string()+"e1bq"+jetSuffix[i]).c_str(),eXbq__[c][i]->GetBinContent(2));
      e2bq[c][i] = new RooConstVar((std::string()+"e2bq"+channelSuffix[c]+jtSuf[i]).c_str(),(std::string()+"e2bq"+jetSuffix[i]).c_str(),eXbq__[c][i]->GetBinContent(3));
      e3bq[c][i] = new RooConstVar((std::string()+"e3bq"+channelSuffix[c]+jtSuf[i]).c_str(),(std::string()+"e3bq"+jetSuffix[i]).c_str(),eXbq__[c][i]->Integral(0,-1)-(e0bq[c][i]->getVal()+e1bq[c][i]->getVal()+e2bq[c][i]->getVal()));

      /*
       RooConstVar e0bq("e0bq","e0bq",e0bq_[JetIdx]);
       RooConstVar e1bq("e1bq","e1bq",e1bq_[JetIdx]);
       RooConstVar e2bq("e2bq","e2bq",e2bq_[JetIdx]);
       */
        //  RooConstVar e3bq("e3bq","e3bq",e3bq__->GetBinContent(JetIdx+1));
      
      n0bjets_st[c][i] = new RooConstVar*[NbOfWP]; n1bjets_st[c][i] = new RooConstVar*[NbOfWP]; n2bjets_st[c][i] = new RooConstVar*[NbOfWP];
      n3bjets_st[c][i] = new RooConstVar*[NbOfWP]; n4bjets_st[c][i] = new RooConstVar*[NbOfWP]; n5bjets_st[c][i] = new RooConstVar*[NbOfWP];
      n0bjets_vb[c][i] = new RooConstVar*[NbOfWP]; n1bjets_vb[c][i] = new RooConstVar*[NbOfWP]; n2bjets_vb[c][i] = new RooConstVar*[NbOfWP];
      n3bjets_vb[c][i] = new RooConstVar*[NbOfWP]; n4bjets_vb[c][i] = new RooConstVar*[NbOfWP]; n5bjets_vb[c][i] = new RooConstVar*[NbOfWP];
      n0bjets_qcd[c][i] = new RooConstVar*[NbOfWP]; n1bjets_qcd[c][i] = new RooConstVar*[NbOfWP]; n2bjets_qcd[c][i] = new RooConstVar*[NbOfWP];
      n3bjets_qcd[c][i] = new RooConstVar*[NbOfWP]; n4bjets_qcd[c][i] = new RooConstVar*[NbOfWP]; n5bjets_qcd[c][i] = new RooConstVar*[NbOfWP];
      n0bjets_vv[c][i] = new RooConstVar*[NbOfWP]; n1bjets_vv[c][i] = new RooConstVar*[NbOfWP]; n2bjets_vv[c][i] = new RooConstVar*[NbOfWP];
      n3bjets_vv[c][i] = new RooConstVar*[NbOfWP]; n4bjets_vv[c][i] = new RooConstVar*[NbOfWP]; n5bjets_vv[c][i] = new RooConstVar*[NbOfWP];
      for (UInt_t j=0; j<NbOfWP; j++) {
        n0bjets_st[c][i][j] = new RooConstVar((std::string()+"n0bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{0 b-jet}", hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(0.)));
        n1bjets_st[c][i][j] = new RooConstVar((std::string()+"n1bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{1 b-jet}", hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(1.)));
        n2bjets_st[c][i][j] = new RooConstVar((std::string()+"n2bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{2 b-jets}",hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(2.)));
        n3bjets_st[c][i][j] = new RooConstVar((std::string()+"n3bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{3 b-jets}",hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(3.)));
        n4bjets_st[c][i][j] = new RooConstVar((std::string()+"n4bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{4 b-jets}",hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(4.)));
        //n5bjets_st[c][i][j] = new RooConstVar((std::string()+"n5bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_stjets[c][i][j]->Integral(hNbtaggedJets_stjets[c][i][j]->FindBin(5.),-1));
        n5bjets_st[c][i][j] = new RooConstVar((std::string()+"n5bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(5.)/*,-1*/));

        n0bjets_vb[c][i][j] = new RooConstVar((std::string()+"n0bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{0 b-jet}", hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(0.)));
        n1bjets_vb[c][i][j] = new RooConstVar((std::string()+"n1bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{1 b-jet}", hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(1.)));
        n2bjets_vb[c][i][j] = new RooConstVar((std::string()+"n2bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{2 b-jets}",hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(2.)));
        n3bjets_vb[c][i][j] = new RooConstVar((std::string()+"n3bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{3 b-jets}",hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(3.)));
        n4bjets_vb[c][i][j] = new RooConstVar((std::string()+"n4bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{4 b-jets}",hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(4.)));
        //n5bjets_vb[c][i][j] = new RooConstVar((std::string()+"n5bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_vbjets[c][i][j]->Integral(hNbtaggedJets_vbjets[c][i][j]->FindBin(5.),-1));
        n5bjets_vb[c][i][j] = new RooConstVar((std::string()+"n5bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(5.)/*,-1*/));

        if (estimateQCDApart) {
          n0bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n0bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{0 b-jet}", hNbtaggedJets_qcd_shape[c][i][j]->GetBinContent(hNbtaggedJets_qcd_shape[c][i][j]->FindBin(0.)));
          n1bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n1bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{1 b-jet}", hNbtaggedJets_qcd_shape[c][i][j]->GetBinContent(hNbtaggedJets_qcd_shape[c][i][j]->FindBin(1.)));
          n2bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n2bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{2 b-jets}",hNbtaggedJets_qcd_shape[c][i][j]->GetBinContent(hNbtaggedJets_qcd_shape[c][i][j]->FindBin(2.)));
          n3bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n3bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{3 b-jets}",hNbtaggedJets_qcd_shape[c][i][j]->GetBinContent(hNbtaggedJets_qcd_shape[c][i][j]->FindBin(3.)));
          n4bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n4bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{4 b-jets}",hNbtaggedJets_qcd_shape[c][i][j]->GetBinContent(hNbtaggedJets_qcd_shape[c][i][j]->FindBin(4.)));
          //          n5bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n5bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_qcd_shape[c][i][j]->Integral(hNbtaggedJets_qcd_shape[c][i][j]->FindBin(5.),-1));
          n5bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n5bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_qcd_shape[c][i][j]->GetBinContent(hNbtaggedJets_qcd_shape[c][i][j]->FindBin(5.)/*,-1*/));
        } else {
          n0bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n0bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{0 b-jet}", hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(0.)));
          n1bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n1bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{1 b-jet}", hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(1.)));
          n2bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n2bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{2 b-jets}",hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(2.)));
          n3bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n3bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{3 b-jets}",hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(3.)));
          n4bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n4bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{4 b-jets}",hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(4.)));
          //          n5bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n5bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_qcd[c][i][j]->Integral(hNbtaggedJets_qcd[c][i][j]->FindBin(5.),-1));            
          n5bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n5bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(5.)/*,-1*/));            
        }
        
        n0bjets_vv[c][i][j] = new RooConstVar((std::string()+"n0bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{0 b-jet}", hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(0.)));
        n1bjets_vv[c][i][j] = new RooConstVar((std::string()+"n1bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{1 b-jet}", hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(1.)));
        n2bjets_vv[c][i][j] = new RooConstVar((std::string()+"n2bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{2 b-jets}",hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(2.)));
        n3bjets_vv[c][i][j] = new RooConstVar((std::string()+"n3bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{3 b-jets}",hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(3.)));
        n4bjets_vv[c][i][j] = new RooConstVar((std::string()+"n4bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{4 b-jets}",hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(4.)));
        //        n5bjets_vv[c][i][j] = new RooConstVar((std::string()+"n5bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_vvjets[c][i][j]->Integral(hNbtaggedJets_vvjets[c][i][j]->FindBin(5.),-1));
        n5bjets_vv[c][i][j] = new RooConstVar((std::string()+"n5bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(5.)/*,-1*/));
      }
    }
  }
  
  
  
  printf("Categories\n");
  
    // C o n s t r u c t   a   c a t e g o r y   w i t h   l a b e l s    a n d   i n d e c e s
    // -----------------------------------------------------------------------------------------
  
  RooCategory nbjets("nbjets","Number of b-jets");
  nbjets.defineType("N0bjet", 0);
  nbjets.defineType("N1bjet", 1);
  nbjets.defineType("N2bjets",2);
  nbjets.defineType("N3bjets",3);
  nbjets.defineType("N4bjets",4);
  nbjets.defineType("N5bjets",5);
  
  RooCategory jetMultCat("jetMult","Jet multiplicity bin");
  for (UInt_t i=0; i<NbOfJetBins; i++) {
    jetMultCat.defineType(jetSuffix[i], i);
  }
  
  RooCategory wpCat("WP","b-tagging Working Point");
  for (UInt_t i=0; i<NbOfWP; i++) {
    wpCat.defineType(wpSuffix[i], i);
  }
  
  RooCategory wpJetCat("WP_Jets","Product category for jet mult and WP");
  for (UInt_t i=0; i<NbOfJetBins; i++) {
    for (UInt_t j=0; j<NbOfWP; j++) {
      wpJetCat.defineType((std::string()+jetSuffix[i]+wpSuffix[j]).c_str(), i*NbOfWP+j);
    }
  }
  
  printf("Construct formula\n");
  
    // C o n s t r u c t   f o r m u l a s
    // -------------------------------------------------------------------------------
  RooFormulaVar ****p0bjets_tt = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****p1bjets_tt = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****p2bjets_tt = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****p3bjets_tt = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****p4bjets_tt = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****p5bjets_tt = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****p0bjets_v  = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****p1bjets_v  = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****p2bjets_v  = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****p3bjets_v  = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****p4bjets_v  = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****p5bjets_v  = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****geq1b_tt   = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****geq1b_tt_n = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****geq1b_v    = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****geq1b_v_n  = new RooFormulaVar***[NbOfChannels];
  RooFormulaVar ****geq1b_tot_n= new RooFormulaVar***[NbOfChannels];
  //  RooFormulaVar ****geq1b_tt_n_add = new RooFormulaVar***[NbOfChannels];
  //  RooFormulaVar ****geq1b_v_n_add  = new RooFormulaVar***[NbOfChannels];
  RooAddition ***geq1b_tt_jInc_n = new RooAddition**[NbOfChannels];
  RooAddition ***geq1b_v_jInc_n  = new RooAddition**[NbOfChannels];
  RooAddition ***geq1b_tot_jInc_n= new RooAddition**[NbOfChannels];
  RooArgList ***tt_geq1b_jInc = new RooArgList**[NbOfChannels];
  RooArgList ***v_geq1b_jInc  = new RooArgList**[NbOfChannels];
  RooArgList ***tot_geq1b_jInc= new RooArgList**[NbOfChannels];

  for (UInt_t c=0; c<NbOfChannels; c++) {
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    p0bjets_tt[c] = new RooFormulaVar**[NbOfJetBins];
    p1bjets_tt[c] = new RooFormulaVar**[NbOfJetBins];
    p2bjets_tt[c] = new RooFormulaVar**[NbOfJetBins];
    p3bjets_tt[c] = new RooFormulaVar**[NbOfJetBins];
    p4bjets_tt[c] = new RooFormulaVar**[NbOfJetBins];
    p5bjets_tt[c] = new RooFormulaVar**[NbOfJetBins];
    p0bjets_v[c]  = new RooFormulaVar**[NbOfJetBins];
    p1bjets_v[c]  = new RooFormulaVar**[NbOfJetBins];
    p2bjets_v[c]  = new RooFormulaVar**[NbOfJetBins];
    p3bjets_v[c]  = new RooFormulaVar**[NbOfJetBins];
    p4bjets_v[c]  = new RooFormulaVar**[NbOfJetBins];
    p5bjets_v[c]  = new RooFormulaVar**[NbOfJetBins];
    geq1b_tt[c]   = new RooFormulaVar**[NbOfJetBins];
    geq1b_tt_n[c] = new RooFormulaVar**[NbOfJetBins];
    geq1b_v[c]    = new RooFormulaVar**[NbOfJetBins];
    geq1b_v_n[c]  = new RooFormulaVar**[NbOfJetBins];
    geq1b_tot_n[c]= new RooFormulaVar**[NbOfJetBins];
    //    geq1b_tt_n_add[c] = new RooFormulaVar**[NbOfJetBins];
    //    geq1b_v_n_add[c]  = new RooFormulaVar**[NbOfJetBins];
    tt_geq1b_jInc[c] = new RooArgList*[NbOfWP];
    v_geq1b_jInc[c]  = new RooArgList*[NbOfWP];
    tot_geq1b_jInc[c]= new RooArgList*[NbOfWP];
    for (UInt_t j=0; j<NbOfWP; j++) {
      tt_geq1b_jInc[c][j] = new RooArgList();
      v_geq1b_jInc[c][j] = new RooArgList();
      tot_geq1b_jInc[c][j] = new RooArgList();
    }
    
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      p0bjets_tt[c][i] = new RooFormulaVar*[NbOfWP];
      p1bjets_tt[c][i] = new RooFormulaVar*[NbOfWP];
      p2bjets_tt[c][i] = new RooFormulaVar*[NbOfWP];
      p3bjets_tt[c][i] = new RooFormulaVar*[NbOfWP];
      p4bjets_tt[c][i] = new RooFormulaVar*[NbOfWP];
      p5bjets_tt[c][i] = new RooFormulaVar*[NbOfWP];
      p0bjets_v[c][i]  = new RooFormulaVar*[NbOfWP];
      p1bjets_v[c][i]  = new RooFormulaVar*[NbOfWP];
      p2bjets_v[c][i]  = new RooFormulaVar*[NbOfWP];
      p3bjets_v[c][i]  = new RooFormulaVar*[NbOfWP];
      p4bjets_v[c][i]  = new RooFormulaVar*[NbOfWP];
      p5bjets_v[c][i]  = new RooFormulaVar*[NbOfWP];
      geq1b_tt[c][i]   = new RooFormulaVar*[NbOfWP];
      geq1b_tt_n[c][i] = new RooFormulaVar*[NbOfWP];
      geq1b_v[c][i]    = new RooFormulaVar*[NbOfWP];
      geq1b_v_n[c][i]  = new RooFormulaVar*[NbOfWP];
      geq1b_tot_n[c][i]= new RooFormulaVar*[NbOfWP];
      //      geq1b_tt_n_add[c][i] = new RooFormulaVar*[NbOfWP];
      //      geq1b_v_n_add[c][i]  = new RooFormulaVar*[NbOfWP];
      if (i==NbOfJetBins-1) {
        geq1b_tt_jInc_n[c] = new RooAddition*[NbOfWP];
        geq1b_v_jInc_n[c]  = new RooAddition*[NbOfWP];
        geq1b_tot_jInc_n[c]= new RooAddition*[NbOfWP];
      }
      for (UInt_t j=0;j<NbOfWP; j++) {
        p0bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p0bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p0bjets_tt",
                                                (std::string()+"(1-@0)*(1-@0)*pow((1-@1),@2-2)*@5+(1-@0)*pow((1-@1),@2-1)*@4+pow((1-@1),@2)*@3"
                                                 + "+(pow((1-@0),3)*pow((1-@1),@2-3))*@6").c_str(),RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i]),*(e3bq[c][i])));
        p1bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p1bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p1bjets_tt",
                                                (std::string()+"(2*@0*(1-@0)*pow(1-@1,@2-2)+(1-@0)*(1-@0)*(@2-2)*@1*pow(1-@1,@2-3))*@5+(@0*pow(1-@1,@2-1)+(1-@0)*(@2-1)*@1*pow(1-@1,@2-2))*@4+(@2*@1*pow(1-@1,@2-1))*@3"
                                                 + "+(3*@0*pow(1-@0,2)*pow(1-@1, @2-3)"
                                                 + "+pow(1-@0,3)*(@2-3)*pow(@1,1)*pow(1-@1,@2-3-1))*@6").c_str(),RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i]),*(e3bq[c][i])));
        p2bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p2bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p2bjets_tt",
                                                (std::string()+"(@0*@0*pow(1-@1,@2-2)+2*@0*(1-@0)*(@2-2)*@1*pow(1-@1,@2-3)+(1-@0)*(1-@0)*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4))*@5+(@0*(@2-1)*@1*pow(1-@1,@2-2)+(1-@0)*((@2-1)*(@2-2)/2)*@1*@1*pow(1-@1,@2-3))*@4+((@2*(@2-1)/2)*@1*@1*pow(1-@1,@2-2))*@3"
                                                 + "+(3*pow(@0,2)*pow(1-@0,1)*pow(1-@1,@2-3)"
                                                 + "+3*pow(@0,1)*pow(1-@0,2)*(@2-3-1)*pow(@1,1)*pow(1-@1,@2-3-1)"
                                                 + "+pow(1-@0,3)*((@2-3)*(@2-3-1)/2)*pow(@1,2)*pow(1-@1,@2-3-2))*@6").c_str(),RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i]),*(e3bq[c][i])));
        //        p3bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p3bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p3bjets_tt","(@0*@0*(@2-2)*@1*pow(1-@1,@2-3)+2*@0*(1-@0)*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4)+(@2>4 ? pow((1-@0),2)*((@2-2)*(@2-3)*(@2-4)/6)*pow(@1,3)*pow((1-@1),@2-5) : 0 ))*@5+(@0*((@2-1)*(@2-2)/2)*@1*@1*pow(1-@1,@2-3)+(1-@0)*((@2-1)*(@2-2)*(@2-3)/6)*pow(@1,3)*pow(1-@1,@2-4))*@4+((@2*(@2-1)*(@2-2)/6)*pow(@1,3)*pow(1-@1,@2-3))*@3",RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i])));
        p3bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p3bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p3bjets_tt",
                                                (std::string()+"(@0*@0*(@2-2)*@1*pow(1-@1,@2-3)+2*@0*(1-@0)*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4)+pow((1-@0),2)*((@2-2)*(@2-3)*(@2-4)/6)*pow(@1,3)*pow((1-@1),@2-5))*@5+(@0*((@2-1)*(@2-2)/2)*@1*@1*pow(1-@1,@2-3)+(1-@0)*((@2-1)*(@2-2)*(@2-3)/6)*pow(@1,3)*pow(1-@1,@2-4))*@4+((@2*(@2-1)*(@2-2)/6)*pow(@1,3)*pow(1-@1,@2-3))*@3"
                                                 + "+(pow(@0,3)*pow(1-@1,@2-3)"
                                                 + "+3*pow(@0,2)*pow(1-@0,1)*(@2-3)*pow(@1,1)*pow(1-@1,@2-3-1)"
                                                 + "+3*pow(@0,1)*pow(1-@0,2)*((@2-3)*(@2-4)/2)*pow(@1,2)*pow(1-@1,@2-3-2)"
                                                 + "+pow(1-@0,3)*((@2-3)*(@2-4)*(@2-5)/6)*pow(@1,3)*pow(1-@1,@2-3-3))*@6").c_str(),RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i]),*(e3bq[c][i])));
        p4bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p4bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p4bjets_tt",
                                                (std::string()+"(@0*@0*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4)+2*@0*(1-@0)*((@2-2)*(@2-3)*(@2-4)/6)*pow(@1,3)*pow(1-@1,@2-5)+pow(1-@0,2)*((@2-2)*(@2-3)*(@2-4)*(@2-5)/24)*pow(@1,4)*pow((1-@1),@2-6))*@5+(@0*((@2-1)*(@2-2)*(@2-3)/6)*pow(@1,3)*pow(1-@1,@2-4)+(1-@0)*((@2-1)*(@2-2)*(@2-3)*(@2-4)/24)*pow(@1,4)*pow(1-@1,@2-5))*@4+((@2*(@2-1)*(@2-2)*(@2-3)/24)*pow(@1,4)*pow(1-@1,@2-4))*@3"
                                                 + "+(pow(@0,3)*(@2-3)*pow(@1,1)*pow(1-@1,@2-3-1)"
                                                 + "+3*pow(@0,2)*pow(1-@0,1)*((@2-3)*(@2-4)/2)*pow(@1,2)*pow(1-@1,@2-3-2)"
                                                 + "+3*pow(@0,1)*pow(1-@0,2)*((@2-3)*(@2-4)*(@2-5)/6)*pow(@1,3)*pow(1-@1,@2-3-3)"
                                                 + "+pow(1-@0,3)*((@2-3)*(@2-4)*(@2-5)*(@2-6)/24)*pow(@1,4)*pow(1-@1,@2-3-4))*@6").c_str(),RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i]),*(e3bq[c][i])));
        p5bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p5bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p5bjets_tt",
                                                (std::string()+"(@0*@0*((@2-2)*(@2-3)*(@2-4)/6)*pow(@1,3)*pow(1-@1,@2-5)+2*@0*(1-@0)*((@2-2)*(@2-3)*(@2-4)*(@2-5)/24)*pow(@1,4)*pow(1-@1,@2-6)+pow(1-@0,2)*((@2-2)*(@2-3)*(@2-4)*(@2-5)*(@2-6)/120)*pow(@1,5)*pow((1-@1),@2-7))*@5+(@0*((@2-1)*(@2-2)*(@2-3)*(@2-4)/24)*pow(@1,4)*pow(1-@1,@2-5)+(1-@0)*((@2-1)*(@2-2)*(@2-3)*(@2-4)*(@2-5)/120)*pow(@1,5)*pow(1-@1,@2-6))*@4+((@2*(@2-1)*(@2-2)*(@2-3)*(@2-4)/120)*pow(@1,5)*pow(1-@1,@2-5))*@3"
                                                 + "+(pow(@0,3)*((@2-3)*(@2-4)/2)*pow(@1,2)*pow(1-@1,@2-3-2)"
                                                 + "+3*pow(@0,2)*pow(1-@0,1)*((@2-3)*(@2-4)*(@2-5)/6)*pow(@1,3)*pow(1-@1,@2-3-3)"
                                                 + "+3*pow(@0,1)*pow(1-@0,2)*((@2-3)*(@2-4)*(@2-5)*(@2-6)/24)*pow(@1,4)*pow(1-@1,@2-3-4)"
                                                 + "+pow(1-@0,3)*((@2-3)*(@2-4)*(@2-5)*(@2-6)*(@2-7)/120)*pow(@1,5)*pow(1-@1,@2-3-5))*@6").c_str(),RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i]),*(e3bq[c][i])));
        
        p0bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p0bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p0bjets_v" ,"pow(1-@0,@1)",RooArgList(*(euds[c][j]),*(n[c][i])));
        p1bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p1bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p1bjets_v" ,"@1*@0*pow(1-@0,@1-1)",RooArgList(*(euds[c][j]),*(n[c][i])));
        p2bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p2bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p2bjets_v" ,"(@1*(@1-1)/2)*@0*@0*pow(1-@0,@1-2)",RooArgList(*(euds[c][j]),*(n[c][i])));
        p3bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p3bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p3bjets_v" ,"((@1)*(@1-1)*(@1-2)/6)*pow(@0,3)*pow(1-@0,@1-3)",RooArgList(*(euds[c][j]),*(n[c][i])));
        p4bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p4bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p4bjets_v" ,"((@1)*(@1-1)*(@1-2)*(@1-3)/24)*pow(@0,4)*pow(1-@0,@1-4)",RooArgList(*(euds[c][j]),*(n[c][i])));
        p5bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p5bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p5bjets_v" ,"((@1)*(@1-1)*(@1-2)*(@1-3)*(@1-4)/120)*pow(@0,5)*pow(1-@0,@1-5)",RooArgList(*(euds[c][j]),*(n[c][i])));

        geq1b_tt[c][i][j]   = new RooFormulaVar((std::string()+"geq1b_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "geq1b_tt" ,"@0+@1+@2+@3+@4",RooArgList(*(p1bjets_tt[c][i][j]),*(p2bjets_tt[c][i][j]),*(p3bjets_tt[c][i][j]),*(p4bjets_tt[c][i][j]),*(p5bjets_tt[c][i][j])));
        geq1b_tt_n[c][i][j] = new RooFormulaVar((std::string()+"geq1b_tt_n"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "geq1b_tt_n" ,"@0*(@2+@3+@4+@5+@6)/(@1+@2+@3+@4+@5+@6)",RooArgList(*(Ntt[c][i]),*(p0bjets_tt[c][i][j]),*(p1bjets_tt[c][i][j]),*(p2bjets_tt[c][i][j]),*(p3bjets_tt[c][i][j]),*(p4bjets_tt[c][i][j]),*(p5bjets_tt[c][i][j])));
        geq1b_v[c][i][j]    = new RooFormulaVar((std::string()+"geq1b_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "geq1b_v" ,"@0+@1+@2+@3+@4",RooArgList(*(p1bjets_v[c][i][j]),*(p2bjets_v[c][i][j]),*(p3bjets_v[c][i][j]),*(p4bjets_v[c][i][j]),*(p5bjets_v[c][i][j])));
        geq1b_v_n[c][i][j]  = new RooFormulaVar((std::string()+"geq1b_v_n"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "geq1b_v_n" ,"@0*(@2+@3+@4+@5+@6)/(@1+@2+@3+@4+@5+@6)",RooArgList(*(Nv[c][i]),*(p0bjets_v[c][i][j]),*(p1bjets_v[c][i][j]),*(p2bjets_v[c][i][j]),*(p3bjets_v[c][i][j]),*(p4bjets_v[c][i][j]),*(p5bjets_v[c][i][j])));

        RooRealVar p_st_k((std::string()+"p_st_k"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),(std::string()+"p_st_k"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),(hNbtaggedJets_stjets[c][i][j]->Integral()-hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(0.)))/hNbtaggedJets_stjets[c][i][j]->Integral());
        RooRealVar p_vv_k((std::string()+"p_vv_k"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),(std::string()+"p_vv_k"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),(hNbtaggedJets_vvjets[c][i][j]->Integral()-hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(0.)))/hNbtaggedJets_vvjets[c][i][j]->Integral());
        RooRealVar p_vb_k((std::string()+"p_vb_k"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),(std::string()+"p_vb_k"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),(hNbtaggedJets_vbjets[c][i][j]->Integral()-hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(0.)))/hNbtaggedJets_vbjets[c][i][j]->Integral());
        p_st_k.setConstant();
        p_vv_k.setConstant();
        p_vb_k.setConstant();
        geq1b_tot_n[c][i][j]= new RooFormulaVar((std::string()+"geq1b_tot_n"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), (std::string()+"geq1b_tot_n"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "@0+@1+@2*@3+@4*@5+@6*@7",
                                                RooArgList( *(geq1b_tt_n[c][i][j]), *(geq1b_v_n[c][i][j]),
                                                           *((RooRealVar*) Nst[c][i]), p_st_k,
                                                           *((RooRealVar*) Nvv[c][i]), p_vv_k,
                                                           *((RooRealVar*) Nvb[c][i]), p_vb_k));

        tt_geq1b_jInc[c][j]->add(*(geq1b_tt_n[c][i][j]));
        v_geq1b_jInc[c][j]->add(*(geq1b_v_n[c][i][j]));
        tot_geq1b_jInc[c][j]->add(*(geq1b_tot_n[c][i][j]));

        if (i==NbOfJetBins-1) {
          // !!!! FINIR !!!!
          geq1b_tt_jInc_n[c][j] = new RooAddition((std::string()+"geq1b_tt_jInc_n"+channelSuffix[c]+wpSuf[j]).c_str(), "geq1b_tt_jInc_n", *tt_geq1b_jInc[c][j]);
          geq1b_v_jInc_n[c][j]  = new RooAddition((std::string()+"geq1b_v_jInc_n"+channelSuffix[c]+wpSuf[j]).c_str(), "geq1b_v_jInc_n", *v_geq1b_jInc[c][j]);
          geq1b_tot_jInc_n[c][j] = new RooAddition((std::string()+"geq1b_tot_jInc_n"+channelSuffix[c]+wpSuf[j]).c_str(), "geq1b_tot_jInc_n", *tot_geq1b_jInc[c][j]);
        }
      }
    }
  }
  
  printf("Construct PDFs\n");
    // C o n s t r u c t   p . d . f 's
    // -------------------------------------------------------------------------------
  RooGenericPdf**** pbjets_tt = new RooGenericPdf***[NbOfChannels];
  RooGenericPdf**** pbjets_v = new RooGenericPdf***[NbOfChannels];
  RooExtendPdf****  pbjets_tt_ext = new RooExtendPdf***[NbOfChannels];
  RooExtendPdf****  pbjets_v_ext  = new RooExtendPdf***[NbOfChannels];
  RooGenericPdf**** pbjets_st = new RooGenericPdf***[NbOfChannels];
  RooGenericPdf**** pbjets_vb = new RooGenericPdf***[NbOfChannels];
  RooGenericPdf**** pbjets_qcd = new RooGenericPdf***[NbOfChannels];
  RooGenericPdf**** pbjets_vv = new RooGenericPdf***[NbOfChannels];
  RooAddPdf**** model = new RooAddPdf***[NbOfChannels];
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    pbjets_tt[c] = new RooGenericPdf**[NbOfJetBins];
    pbjets_v[c] = new RooGenericPdf**[NbOfJetBins];
    pbjets_tt_ext[c] = new RooExtendPdf**[NbOfJetBins];
    pbjets_v_ext[c]  = new RooExtendPdf**[NbOfJetBins];
    pbjets_st[c] = new RooGenericPdf**[NbOfJetBins];
    pbjets_vb[c] = new RooGenericPdf**[NbOfJetBins];
    pbjets_qcd[c] = new RooGenericPdf**[NbOfJetBins];
    pbjets_vv[c] = new RooGenericPdf**[NbOfJetBins];
    model[c] = new RooAddPdf**[NbOfJetBins];
    
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      pbjets_tt[c][i] = new RooGenericPdf*[NbOfWP];
      pbjets_v[c][i] = new RooGenericPdf*[NbOfWP];
      pbjets_tt_ext[c][i] = new RooExtendPdf*[NbOfWP];
      pbjets_v_ext[c][i]  = new RooExtendPdf*[NbOfWP];
      pbjets_st[c][i] = new RooGenericPdf*[NbOfWP];
      pbjets_vb[c][i] = new RooGenericPdf*[NbOfWP];
      pbjets_qcd[c][i] = new RooGenericPdf*[NbOfWP];
      pbjets_vv[c][i] = new RooGenericPdf*[NbOfWP];
      model[c][i] = new RooAddPdf*[NbOfWP];
      for (UInt_t j=0;j<NbOfWP; j++) {
        //        printf("\\geq 1 b :  V : $ %lf \\pm %lf $ \n");

        // new RooAddPdf("geq1b_tt","geq1b_tt", RooArgList(*(p0bjets_tt[c][i][j]),*(p1bjets_tt[c][i][j]),*(p2bjets_tt[c][i][j]),*(p3bjets_tt[c][i][j]),*(p4bjets_tt[c][i][j]),*(p5bjets_tt[c][i][j])),
        //                                        ,RooArgList(*(Ntt[c][i]),*(Ntt[c][i]),*(Ntt[c][i]),*(Ntt[c][i]),*(Ntt[c][i]),*(Ntt[c][i])));
        //        geq1b_tt[c][i][j] = new RooAddPdf("geq1b_tt","geq1b_tt", RooArgList(*(p0bjets_tt[c][i][j]),*(p1bjets_tt[c][i][j]),*(p2bjets_tt[c][i][j]),*(p3bjets_tt[c][i][j]),*(p4bjets_tt[c][i][j]),*(p5bjets_tt[c][i][j])),
        //                                          ,RooArgList(*(Ntt[c][i]),*(Ntt[c][i]),*(Ntt[c][i]),*(Ntt[c][i]),*(Ntt[c][i]),*(Ntt[c][i])));
        

        pbjets_tt[c][i][j] = new RooGenericPdf("pbjets_tt","pbjets_tt","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(p0bjets_tt[c][i][j]),*(p1bjets_tt[c][i][j]),*(p2bjets_tt[c][i][j]),*(p3bjets_tt[c][i][j]),*(p4bjets_tt[c][i][j]),*(p5bjets_tt[c][i][j])));
        pbjets_tt_ext[c][i][j] = new RooExtendPdf("pbjets_tt_ext","pbjets_tt_ext",*(pbjets_tt[c][i][j]),*(Ntt[c][i]));
        
        pbjets_v[c][i][j] = new RooGenericPdf("pbjets_v","pbjets_v","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(p0bjets_v[c][i][j]),*(p1bjets_v[c][i][j]),*(p2bjets_v[c][i][j]),*(p3bjets_v[c][i][j]),*(p4bjets_v[c][i][j]),*(p5bjets_v[c][i][j])));
        pbjets_v_ext[c][i][j] = new RooExtendPdf("pbjets_v_ext","pbjets_v_ext",*(pbjets_v[c][i][j]),*(Nv[c][i]));
        
        pbjets_st[c][i][j] = new RooGenericPdf("pbjets_st","pbjets_st","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(n0bjets_st[c][i][j]),*(n1bjets_st[c][i][j]),*(n2bjets_st[c][i][j]),*(n3bjets_st[c][i][j]),*(n4bjets_st[c][i][j]),*(n5bjets_st[c][i][j])));
        
        pbjets_vb[c][i][j] = new RooGenericPdf("pbjets_vb","pbjets_vb","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(n0bjets_vb[c][i][j]),*(n1bjets_vb[c][i][j]),*(n2bjets_vb[c][i][j]),*(n3bjets_vb[c][i][j]),*(n4bjets_vb[c][i][j]),*(n5bjets_vb[c][i][j])));
        pbjets_qcd[c][i][j] = new RooGenericPdf("pbjets_qcd","pbjets_qcd","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(n0bjets_qcd[c][i][j]),*(n1bjets_qcd[c][i][j]),*(n2bjets_qcd[c][i][j]),*(n3bjets_qcd[c][i][j]),*(n4bjets_qcd[c][i][j]),*(n5bjets_qcd[c][i][j])));
        
        pbjets_vv[c][i][j] = new RooGenericPdf("pbjets_vv","pbjets_vv","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(n0bjets_vv[c][i][j]),*(n1bjets_vv[c][i][j]),*(n2bjets_vv[c][i][j]),*(n3bjets_vv[c][i][j]),*(n4bjets_vv[c][i][j]),*(n5bjets_vv[c][i][j])));
        
          //RooAddPdf model("model","model",RooArgList(pbjets_tt,pbjets_v),RooArgList(Ntt,Nv));
        model[c][i][j] = new RooAddPdf("model","model",
                                       (estimateQCDApart==kTRUE ? RooArgList(*(pbjets_tt[c][i][j]),*(pbjets_v[c][i][j]),*(pbjets_st[c][i][j]),*(pbjets_vb[c][i][j]),*(pbjets_qcd[c][i][j]),*(pbjets_vv[c][i][j]))
                                        : RooArgList(*(pbjets_tt[c][i][j]),*(pbjets_v[c][i][j]),*(pbjets_st[c][i][j]),*(pbjets_vb[c][i][j])/*,*(pbjets_qcd[c][i][j])*/,*(pbjets_vv[c][i][j]))),
                                       (estimateQCDApart==kTRUE ? RooArgList(*(Ntt[c][i]),*(Nv[c][i]),*(Nst[c][i]),*(Nvb[c][i]),*(Nqcd[c][i]),*(Nvv[c][i]))
                                       : RooArgList(*(Ntt[c][i]),*(Nv[c][i]),*(Nst[c][i]),*(Nvb[c][i])/*,*(Nqcd[c][i])*/,*(Nvv[c][i]))));
      } 
    }
    
    printf("For channel #%d\n", c);
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      printf("eXbq_%d(0,1,2,3) = (%lf, %lf, %lf, %lf)\n", i, e0bq[c][i]->getVal(), e1bq[c][i]->getVal(), e2bq[c][i]->getVal(), e3bq[c][i]->getVal());
    }
    printf("\n");
    
    printf("\n\nConstraint values (values not always applied as constraint ....) : \n");
    for (UInt_t i=0; i<NbOfWP; i++) {
      printf("WP %d : \n*******\n e_b = %lf +/- %lf\ne_uds  = %lf +/- %lf\ne_udsc = %lf +/- %lf\n", i, Eb_const_mean[c][i], Eb_const_error[c][i], Euds_const_mean[c][i], Euds_const_error[c][i], Eudsc_const_mean[c][i], Eudsc_const_error[c][i]);
    }
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      printf("Nst[%d]    = %lf +/- %lf\nN_Vbb[%d]  = %lf +/- %lf\n", i, Nstjets[c][i], Nstjets[c][i]*Nstjets_uncert[c][i],i, Nvbjets[c][i], Nvbjets[c][i]*Nvbjets_uncert[c][i]);
    }
    printf("\n");
  }
  
    // Contrainsts on free parameters
  RooGaussian ***Eb_constraint = new RooGaussian**[NbOfChannels];
  RooGaussian ***Eudsc_constraint = new RooGaussian**[NbOfChannels];
  RooGaussian ***Euds_constraint = new RooGaussian**[NbOfChannels];
  
  RooGaussian ***Nst_constraint = new RooGaussian**[NbOfChannels];
  RooGaussian ***Nvb_constraint = new RooGaussian**[NbOfChannels];
  RooGaussian ***Nqcd_constraint = new RooGaussian**[NbOfChannels];
  RooGaussian ***Nvv_constraint = new RooGaussian**[NbOfChannels];
  RooProdPdf ****model_constraint = new RooProdPdf***[NbOfChannels];
  
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    Eb_constraint[c] = new RooGaussian*[NbOfWP];
    Eudsc_constraint[c] = new RooGaussian*[NbOfWP];
    Euds_constraint[c] = new RooGaussian*[NbOfWP];
    for (UInt_t wp=0; wp<NbOfWP; wp++) {
      Eb_constraint[c][wp] = new RooGaussian((std::string()+"Eb_constraint"+wpSuffix[wp]).c_str(),"Eb_constraint",*(eb[c][wp]),RooConst(Eb_const_mean[c][wp]),RooConst(Eb_const_error[c][wp]));
      Eudsc_constraint[c][wp] = new RooGaussian((std::string()+"Eudsc_constraint"+wpSuffix[wp]).c_str(),"Eudsc_constraint",*(eudsc[c][wp]),RooConst(Eudsc_const_mean[c][wp]),RooConst(Eudsc_const_error[c][wp]));
      Euds_constraint[c][wp] = new RooGaussian((std::string()+"Euds_constraint"+wpSuffix[wp]).c_str(),"Euds_constraint",*(euds[c][wp]),RooConst(Euds_const_mean[c][wp]),RooConst(Euds_const_error[c][wp]));
    }
    Nst_constraint[c] = new RooGaussian*[NbOfJetBins];
    Nvb_constraint[c] = new RooGaussian*[NbOfJetBins];
    Nqcd_constraint[c] = new RooGaussian*[NbOfJetBins];
    Nvv_constraint[c] = new RooGaussian*[NbOfJetBins];
    model_constraint[c] = new RooProdPdf**[NbOfJetBins];
    
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      Nst_constraint[c][i]   = new RooGaussian((std::string()+"Nst_constraint"+jetSuffix[i]).c_str(),"Nst_constraint",*(Nst[c][i]),RooConst(Nstjets[c][i]),RooConst(Nstjets[c][i]*Nstjets_uncert[c][i]));
      Nvb_constraint[c][i]   = new RooGaussian((std::string()+"Nvbb_constraint"+jetSuffix[i]).c_str(),"Nvbb_constraint",*(Nvb[c][i]),RooConst(Nvbjets[c][i]),RooConst(Nvbjets[c][i]*Nvbjets_uncert[c][i]));
      Nqcd_constraint[c][i]  = new RooGaussian((std::string()+"Nqcd_constraint"+jetSuffix[i]).c_str(),"Nqcd_constraint",*(Nqcd[c][i]),RooConst(Nqcd_mean[c][i]),RooConst(Nqcd_mean[c][i]*Nqcd_uncert[c][i]/*+100000.*/));
      Nvv_constraint[c][i]   = new RooGaussian((std::string()+"Nvv_constraint"+jetSuffix[i]).c_str(),"Nvv_constraint",*(Nvv[c][i]),RooConst(Nvvjets[c][i]),RooConst(Nvvjets[c][i]*Nvvjets_uncert[c][i]/*+100000.*/));
      model_constraint[c][i] = new RooProdPdf*[NbOfWP];
      for (UInt_t j=0 ; j<NbOfWP; j++) {
        model_constraint[c][i][j] = new RooProdPdf((std::string()+"model_constraint"+jetSuffix[i]+wpSuffix[j]).c_str(),"model with constraint",(estimateQCDApart==kTRUE ?
          RooArgSet(*(model[c][i][j]),*(Eb_constraint[c][j]),*(Euds_constraint[c][j]),*(Eudsc_constraint[c][j]),*(Nst_constraint[c][i]),*(Nvb_constraint[c][i]),*(Nqcd_constraint[c][i]),*(Nvv_constraint[c][i])) :
          RooArgSet(*(model[c][i][j]),*(Eb_constraint[c][j]),*(Euds_constraint[c][j]),*(Eudsc_constraint[c][j]),*(Nst_constraint[c][i]),*(Nvb_constraint[c][i])/*,*(Nqcd_constraint[c][i])*/,*(Nvv_constraint[c][i])))
) ;
      }
    }
  }
  
  printf("Before building RooSimultaneous\n");
  
    //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,Eb_constraint)) ;
    //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,Eudsc_constraint)) ;
    //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,Euds_constraint)) ;
  
  
  /*
   RooSimultaneous simPdf("SeveralJetMult","simultaneous pdf",jetMultCat) ;
   for(UInt_t i=0; i<NbOfJetBins; i++) {
   RooSimultaneous jetMultSim((std::string()+"SeveralJetMult"+jetSuffix[i]).c_str(),"simultaneous pdf",wpCat) ;
   for(UInt_t j=0; j<NbOfWP; j++) {
   jetMultSim.addPdf(*(model_constraint[i][j]),wpSuffix[j]) ;
   }
   simPdf.addPdf(jetMultSim,jetSuffix[i]);
   }  
   */
  RooSimultaneous** simPdf = new RooSimultaneous*[NbOfChannels];
  for (UInt_t c=0; c<NbOfChannels; c++) {
    simPdf[c] = NULL;
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    simPdf[c] = new RooSimultaneous((std::string()+"SeveralJetMult"+channelSuffix[c]).c_str(),"simultaneous pdf",wpJetCat) ;
    for(UInt_t i=0; i<NbOfJetBins; i++) {
      for(UInt_t j=0; j<NbOfWP; j++) {
        simPdf[c]->addPdf(*(model_constraint[c][i][j]),(std::string()+jetSuffix[i]+wpSuffix[j]).c_str()) ;
      }
    }
  }
  printf("After building RooSimultaneous\n");
    // C r e a t e   d a t a s e t 
    // -------------------------------------------------------------------------------
  
  RooDataSet** data = new RooDataSet*[NbOfChannels];
  for (UInt_t c=0; c<NbOfChannels; c++) {
    data[c] = NULL;
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    data[c] = new RooDataSet((std::string()+"data"+channelSuffix[c]).c_str(),(std::string()+"data"+channelSuffix[c]).c_str(),RooArgSet(nbjets,wpJetCat)) ;
    /*
     for (int j=0; j<NbOfJetBins; j++) {
     jetMultCat.setLabel(jetSuffix[j]);
     for (int wp=0; wp<NbOfWP; wp++) {
     wpCat.setLabel(wpSuffix[wp]);
     nbjets.setLabel("N0bjet"); 
     for (Int_t i=0 ; i<hNbtaggedJets[j][wp]->GetBinContent(hNbtaggedJets[j][wp]->FindBin(0.)) ; i++) { data.add(RooArgSet(nbjets,jetMultCat,wpCat));}
     nbjets.setLabel("N1bjet"); 
     for (Int_t i=0 ; i<hNbtaggedJets[j][wp]->GetBinContent(hNbtaggedJets[j][wp]->FindBin(1.)) ; i++) { data.add(RooArgSet(nbjets,jetMultCat,wpCat));}
     nbjets.setLabel("N2bjets"); 
     for (Int_t i=0 ; i<hNbtaggedJets[j][wp]->GetBinContent(hNbtaggedJets[j][wp]->FindBin(2.)) ; i++) { data.add(RooArgSet(nbjets,jetMultCat,wpCat));}
     nbjets.setLabel("N3bjets"); 
     for (Int_t i=0 ; i<hNbtaggedJets[j][wp]->Integral(hNbtaggedJets[j][wp]->FindBin(3.), hNbtaggedJets[j][wp]->FindBin(10.)) ; i++) { data.add(RooArgSet(nbjets,jetMultCat,wpCat));}
     }
     }*/
    for (UInt_t j=0; j<NbOfJetBins; j++) {
      for (UInt_t wp=0; wp<NbOfWP; wp++) {
        wpJetCat.setLabel((std::string()+jetSuffix[j]+wpSuffix[wp]).c_str());
        nbjets.setLabel("N0bjet"); 
        for (Int_t i=0 ; i<hNbtaggedJets[c][j][wp]->GetBinContent(hNbtaggedJets[c][j][wp]->FindBin(0.)) ; i++) { data[c]->add(RooArgSet(nbjets,wpJetCat));}
        nbjets.setLabel("N1bjet"); 
        for (Int_t i=0 ; i<hNbtaggedJets[c][j][wp]->GetBinContent(hNbtaggedJets[c][j][wp]->FindBin(1.)) ; i++) { data[c]->add(RooArgSet(nbjets,wpJetCat));}
        nbjets.setLabel("N2bjets"); 
        for (Int_t i=0 ; i<hNbtaggedJets[c][j][wp]->GetBinContent(hNbtaggedJets[c][j][wp]->FindBin(2.)) ; i++) { data[c]->add(RooArgSet(nbjets,wpJetCat));}
        nbjets.setLabel("N3bjets"); 
        for (Int_t i=0 ; i<hNbtaggedJets[c][j][wp]->GetBinContent(hNbtaggedJets[c][j][wp]->FindBin(3.)) ; i++) { data[c]->add(RooArgSet(nbjets,wpJetCat));}
        nbjets.setLabel("N4bjets"); 
        for (Int_t i=0 ; i<hNbtaggedJets[c][j][wp]->GetBinContent(hNbtaggedJets[c][j][wp]->FindBin(4.)) ; i++) { data[c]->add(RooArgSet(nbjets,wpJetCat));}
        nbjets.setLabel("N5bjets"); 
        for (Int_t i=0 ; i<hNbtaggedJets[c][j][wp]->Integral(hNbtaggedJets[c][j][wp]->FindBin(5.), hNbtaggedJets[c][j][wp]->FindBin(10.)) ; i++) { data[c]->add(RooArgSet(nbjets,wpJetCat));}
      }
    }
      //data.Print("v");
    
    Roo1DTable* table = data[c]->table(RooArgSet(nbjets,wpJetCat));
      //  Roo1DTable* table = data.table(RooArgSet(nbjets,jetMultCat,wpCat));
    table->Print("v");
  }
  
  
    // F i t   t h e   d a t a   a n d   c o n s t r u c t   t h e   l i k e l i h o o d   f u n c t i o n
    // ----------------------------------------------------------------------------------------------
  /*
   //RooAbsReal* nll = model.createNLL(data);//,Optimize(0));
   RooAbsReal* nll = model_constraint.createNLL(data,NumCPU(4));//,Optimize(0));
   
   RooMinimizer minimizer(*nll);
   
   //minimizer.optimizeConst(0) ; DO NOT SET IT TO TRUE, WILL NOT CONVERGE OTHERWISE
   minimizer.setPrintLevel(-1);
   //minimizer.setNoWarn();
   
   // Set algorithm
   minimizer.minimize("Minuit2", "Combined");
   //minimizer.minimize("GSLMultiMin", "ConjugateFR");
   //minimizer.minimize("GSLMultiMin", "BFGS2");
   minimizer.minos();
   */	
  
    //eb
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    for (UInt_t wp=0; wp<NbOfWP; wp++) {
      printf("Efficiencies before fit    :    e_b[%d][%d] = %lf   -   e_uds[%d][%d] = %lf   -   e_udsc[%d][%d] = %lf\n", c, wp, init_Eb[c][wp], c, wp, init_Euds[c][wp], c, wp, init_Eudsc[c][wp]);
    }
  }
  
  
  
  
  
  RooFitResult **fit_result = new RooFitResult*[NbOfChannels];
    //  if(!doPE)
    //  fit_result = model_constraint.fitTo(data,Constrain(RooArgSet(eb,euds,eudsc,Nst,Nvb)),Save(1),Extended(1),Minos(1),Strategy(2),Optimize(0));
  RooArgSet** constraintsList = new RooArgSet*[NbOfChannels];
  RooArgSet** consList = new RooArgSet*[NbOfChannels];
  for (UInt_t c=0; c<NbOfChannels; c++) {
    fit_result[c] = NULL;
    
      //    if (!((((c+1)&muonMask)!=0) && (((c+1)&electronMask)!=0))) continue;
    
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    
    constraintsList[c] = new RooArgSet();
    for (UInt_t wp=0; wp<NbOfWP; wp++){
      constraintsList[c]->add(*(eb[c][wp]));
      constraintsList[c]->add(*(euds[c][wp]));
      constraintsList[c]->add(*(eudsc[c][wp]));
    }
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      constraintsList[c]->add(*(Nst[c][i]));
      constraintsList[c]->add(*(Nvb[c][i]));
        //      constraintsList[c]->add(*(Nqcd[c][i]));
      constraintsList[c]->add(*(Nvv[c][i]));
    }
    
    consList[c] = new RooArgSet();
    for (UInt_t wp=0; wp<NbOfWP; wp++){
      consList[c]->add(*(Eb_constraint[c][wp]));
      consList[c]->add(*(Euds_constraint[c][wp]));
      consList[c]->add(*(Eudsc_constraint[c][wp]));
    }
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      consList[c]->add(*(Nst_constraint[c][i]));
      consList[c]->add(*(Nvb_constraint[c][i]));
        //      consList[c]->add(*(Nqcd_constraint[c][i]));
      consList[c]->add(*(Nvv_constraint[c][i]));
    }
    printf("\n\n\nTTTT Input values for Fit for channel #%d\n", c);
    for (UInt_t j=0; j<NbOfWP; j++) {
      printf("  %s     : %lf\n  %s  : %lf\n  %s : %lf\n", 
             eb[c][j]->getTitle().Data(), eb[c][j]->getVal(),
             euds[c][j]->getTitle().Data(), euds[c][j]->getVal(),
             eudsc[c][j]->getTitle().Data(), eudsc[c][j]->getVal());
    }
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      printf("  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s : %lf\n",
             Ntt[c][i]->getTitle().Data(), Ntt[c][i]->getVal(),
             Nv[c][i]->getTitle().Data(), Nv[c][i]->getVal(),
             Nst[c][i]->getTitle().Data(), Nst[c][i]->getVal(),
             Nvb[c][i]->getTitle().Data(), Nvb[c][i]->getVal(),
             Nvv[c][i]->getTitle().Data(), Nvv[c][i]->getVal(),
             Nqcd[c][i]->getTitle().Data(), Nqcd[c][i]->getVal());
    }
    printf("\n\n\n");
    
    fit_result[c] = simPdf[c]->fitTo(*(data[c]),Constrain(*(constraintsList[c])),Save(1),Extended(1),Minos(1),Strategy(2),Optimize(0),NumCPU(4));
    
    printf("\n\n\nTTTT Output values from Fit for channel #%d\n", c);
    for (UInt_t j=0; j<NbOfWP; j++) {
      printf("  %s     : %lf\n  %s  : %lf\n  %s : %lf\n", 
             eb[c][j]->getTitle().Data(), eb[c][j]->getVal(),
             euds[c][j]->getTitle().Data(), euds[c][j]->getVal(),
             eudsc[c][j]->getTitle().Data(), eudsc[c][j]->getVal());
    }
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      printf("  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s : %lf\n",
             Ntt[c][i]->getTitle().Data(), Ntt[c][i]->getVal(),
             Nv[c][i]->getTitle().Data(), Nv[c][i]->getVal(),
             Nst[c][i]->getTitle().Data(), Nst[c][i]->getVal(),
             Nvb[c][i]->getTitle().Data(), Nvb[c][i]->getVal(),
             Nvv[c][i]->getTitle().Data(), Nvv[c][i]->getVal(),
             Nqcd[c][i]->getTitle().Data(), Nqcd[c][i]->getVal());
    }
    printf("\n\n\n");
    
      // ,Verbose(kTRUE),PrintLevel(3),PrintEvalErrors(3));
      //  if(!doPE)   fit_result = model_constraint.fitTo(data,Constrain(RooArgSet(euds,Nst,Nvb)),Save(1),Extended(1),Minos(1),Strategy(2));
    
      //fit_result[c] = simPdf[c]->fitTo(*(data[c]),NumCPU(4),Save(1),Extended(1),Minos(1),Strategy(2),ExternalConstraints(* (consList[c])),Optimize(0));
    
      //RooFitResult* fit_result = minimizer.save();
    
    if(fit_result[c]) {
      fit_result[c]->Print("v");
      fit_result[c]->correlationMatrix().Print() ;
 
      for (UInt_t i=0; i<NbOfJetBins ; i++) {
        for (UInt_t j=0; j<NbOfWP ; j++) {
          printf("\\geq 1 b (channel %d ; jet %d ; WP %d) : TT : $ %lf \\pm %lf $ \n", c, i, j, geq1b_tt[c][i][j]->getVal(),    geq1b_tt[c][i][j]->getPropagatedError(*(fit_result[c])));
          printf("\\geq 1 b (channel %d ; jet %d ; WP %d) : TTn: $ %lf \\pm %lf $ \n", c, i, j, geq1b_tt_n[c][i][j]->getVal(),  geq1b_tt_n[c][i][j]->getPropagatedError(*(fit_result[c])));
          printf("\\geq 1 b (channel %d ; jet %d ; WP %d) :  V : $ %lf \\pm %lf $ \n", c, i, j, geq1b_v[c][i][j]->getVal(),     geq1b_v[c][i][j]->getPropagatedError(*(fit_result[c])));
          printf("\\geq 1 b (channel %d ; jet %d ; WP %d) :  Vn: $ %lf \\pm %lf $ \n", c, i, j, geq1b_v_n[c][i][j]->getVal(),   geq1b_v_n[c][i][j]->getPropagatedError(*(fit_result[c])));
        }
      }
      for (UInt_t j=0; j<NbOfWP ; j++) {
        printf("\\geq 1 b (channel %d ; WP %d) : TTnInc: $ %lf \\pm %lf $ \n", c, j, geq1b_tt_jInc_n[c][j]->getVal(),  geq1b_tt_jInc_n[c][j]->getPropagatedError(*(fit_result[c])));
        printf("\\geq 1 b (channel %d ; WP %d) :  VnInc: $ %lf \\pm %lf $ \n", c, j, geq1b_v_jInc_n[c][j]->getVal(),   geq1b_v_jInc_n[c][j]->getPropagatedError(*(fit_result[c])));
        printf("\\geq 1 b (channel %d ; WP %d) :totnInc: $ %lf \\pm %lf $ \n", c, j, geq1b_tot_jInc_n[c][j]->getVal(),   geq1b_tot_jInc_n[c][j]->getPropagatedError(*(fit_result[c])));
      }
      for (UInt_t j=0; j<NbOfWP ; j++) {
        Double_t tot = 0.;
        printf("\n\nTotal numbers forFor WP %d : %s\n\n", j, wpSuffix[j]);
        printf("\\begin{tabular}{r|ccc|l}\n");
        printf(" & $t\\bar{t}$-like & $V$-like & total est. & Data \\\\ \\hline \n");
        for (UInt_t i=0; i<NbOfJetBins ; i++) {
          if (i==NbOfJetBins-1 && !isLastBinExclusive) {
            printf(" $\\geq %d$ jets & $", MinNbOfJets+i);
          } else {
            printf(" $= %d$ jets & $", MinNbOfJets+i);
          }
          printfValErr(geq1b_tt_n[c][i][j]->getVal(),  geq1b_tt_n[c][i][j]->getPropagatedError(*(fit_result[c])));
          printf(" $ & $ ");
          printfValErr(geq1b_v_n[c][i][j]->getVal(),  geq1b_v_n[c][i][j]->getPropagatedError(*(fit_result[c])));
          printf(" $ & $ ");
          printfValErr(geq1b_tot_n[c][i][j]->getVal(),  geq1b_tot_n[c][i][j]->getPropagatedError(*(fit_result[c])));
          printf(" $ & $ ");
          printf(" %.0lf $ \\\\ \n", hNbtaggedJets[c][i][j]->Integral()-hNbtaggedJets[c][i][j]->GetBinContent(hNbtaggedJets[c][i][j]->FindBin(0.)));

          if (i==NbOfJetBins-1 && !isLastBinExclusive) {
            printf("%% $\\\\geq %d$ jets & $", MinNbOfJets+i);
          } else {
            printf(" $= %d$ jets & $", MinNbOfJets+i);
          }
          printf("%% $ %.0lf \\pm %.0lf $ & $ %.0lf \\pm %.0lf $ & $ %.0lf \\pm %.0lf $ & $ %.0lf $ \\\\ \n",
                 geq1b_tt_n[c][i][j]->getVal(),  geq1b_tt_n[c][i][j]->getPropagatedError(*(fit_result[c])),
                 geq1b_v_n[c][i][j]->getVal(),  geq1b_v_n[c][i][j]->getPropagatedError(*(fit_result[c])),
                 geq1b_tot_n[c][i][j]->getVal(),  geq1b_tot_n[c][i][j]->getPropagatedError(*(fit_result[c])),
                 hNbtaggedJets[c][i][j]->Integral()-hNbtaggedJets[c][i][j]->GetBinContent(hNbtaggedJets[c][i][j]->FindBin(0.)));
          tot += hNbtaggedJets[c][i][j]->Integral()-hNbtaggedJets[c][i][j]->GetBinContent(hNbtaggedJets[c][i][j]->FindBin(0.));
        }
        printf("\\hline \n");
        printf(" & $ ");
        printfValErr(geq1b_tt_jInc_n[c][j]->getVal(),  geq1b_tt_jInc_n[c][j]->getPropagatedError(*(fit_result[c])));
        printf(" $ & $ ");
        printfValErr(geq1b_v_jInc_n[c][j]->getVal(),  geq1b_v_jInc_n[c][j]->getPropagatedError(*(fit_result[c])));
        printf(" $ & $ ");
        printfValErr(geq1b_tot_jInc_n[c][j]->getVal(),  geq1b_tot_jInc_n[c][j]->getPropagatedError(*(fit_result[c])));
        printf(" $ & $ ");
        printf(" %.0lf $ \\\\ \n", tot);

        printf("%% & $ %.0lf \\pm %.0lf $ & %.0lf \\pm %.0lf $ & %.0lf \\pm %.0lf $ & $ %.0lf $ \\\\ \n",
               geq1b_tt_jInc_n[c][j]->getVal(),  geq1b_tt_jInc_n[c][j]->getPropagatedError(*(fit_result[c])),
               geq1b_v_jInc_n[c][j]->getVal(),  geq1b_v_jInc_n[c][j]->getPropagatedError(*(fit_result[c])),
               geq1b_tot_jInc_n[c][j]->getVal(),  geq1b_tot_jInc_n[c][j]->getPropagatedError(*(fit_result[c])),
               tot);
        
        printf("\\end{tabular}\n\n\n");
      }

    }
  }
  
  
  
  char label[100];
  
  
  TH1F ****histoTT_PDF_Fit = new TH1F***[NbOfChannels];
  TH1F ****histoV_PDF_Fit = new TH1F***[NbOfChannels];
  TH1F ****histoVbb_PDF_Fit = new TH1F***[NbOfChannels];
  TH1F ****histoST_PDF_Fit = new TH1F***[NbOfChannels];
  TH1F ****histoQCD_PDF_Fit = new TH1F***[NbOfChannels];
  TH1F ****histoVV_PDF_Fit = new TH1F***[NbOfChannels];
  TH1F ****histoTot_PDF_Fit = new TH1F***[NbOfChannels];
  TH1F ****histoTotC_PDF_Fit = new TH1F***[NbOfChannels];
  TH1F ****histoTotS_PDF_Fit = new TH1F***[NbOfChannels];
  for (UInt_t c=0; c<NbOfChannels; c++) {
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    histoTT_PDF_Fit[c] = new TH1F**[NbOfJetBins];
    histoV_PDF_Fit[c] = new TH1F**[NbOfJetBins];
    histoVbb_PDF_Fit[c] = new TH1F**[NbOfJetBins];
    histoST_PDF_Fit[c] = new TH1F**[NbOfJetBins];
    histoQCD_PDF_Fit[c] = new TH1F**[NbOfJetBins];
    histoVV_PDF_Fit[c] = new TH1F**[NbOfJetBins];
    histoTot_PDF_Fit[c] = new TH1F**[NbOfJetBins];
    histoTotC_PDF_Fit[c] = new TH1F**[NbOfJetBins];
    histoTotS_PDF_Fit[c] = new TH1F**[NbOfJetBins];
    for(UInt_t i=0; i<NbOfJetBins ; i++) {
      histoTT_PDF_Fit[c][i] = new TH1F*[NbOfWP];
      histoV_PDF_Fit[c][i] = new TH1F*[NbOfWP];
      histoVbb_PDF_Fit[c][i] = new TH1F*[NbOfWP];
      histoST_PDF_Fit[c][i] = new TH1F*[NbOfWP];
      histoQCD_PDF_Fit[c][i] = new TH1F*[NbOfWP];
      histoVV_PDF_Fit[c][i] = new TH1F*[NbOfWP];
      histoTot_PDF_Fit[c][i] = new TH1F*[NbOfWP];
      histoTotC_PDF_Fit[c][i] = new TH1F*[NbOfWP];
      histoTotS_PDF_Fit[c][i] = new TH1F*[NbOfWP];
      RooArgSet *ras = NULL;
        //ras = new RooArgSet(*(Ntt[i]),*(Nv[i]),*(Nst[i]),*(Nvb[i]));
      for (UInt_t j=0; j<NbOfWP ; j++) {
        histoTT_PDF_Fit[c][i][j] = new TH1F((std::string()+"histoTT_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoTT_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoV_PDF_Fit[c][i][j] = new TH1F((std::string()+"histoV_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoV_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoVbb_PDF_Fit[c][i][j] = new TH1F((std::string()+"histoVbb_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoVbb_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoST_PDF_Fit[c][i][j] = new TH1F((std::string()+"histoST_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoST_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoQCD_PDF_Fit[c][i][j] = new TH1F((std::string()+"histoQCD_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoQCD_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoVV_PDF_Fit[c][i][j] = new TH1F((std::string()+"histoVV_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoVV_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoTot_PDF_Fit[c][i][j] = new TH1F((std::string()+"histoTot_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoTot_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoTotC_PDF_Fit[c][i][j] = new TH1F((std::string()+"histoTotC_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoTotC_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoTotS_PDF_Fit[c][i][j] = new TH1F((std::string()+"histoTotS_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoTotS_PDF_Fit"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        wpJetCat.setLabel((std::string()+jetSuffix[i]+wpSuffix[j]).c_str());
        if ( fit_result[c]==NULL) {
          continue;
        }
        
        for (UInt_t bj=0; bj<6; bj++) {
          if (bj==0 || bj==1) {
            sprintf(label, "N%dbjet", bj);
          } else {
            sprintf(label, "N%dbjets", bj);
          }
          nbjets.setLabel(label); 
          histoTT_PDF_Fit[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoTT_PDF_Fit[c][i][j]->SetBinContent(bj+1, pbjets_tt_ext[c][i][j]->getVal());
          histoTT_PDF_Fit[c][i][j]->SetBinError(bj+1, pbjets_tt_ext[c][i][j]->getPropagatedError(*(fit_result[c])));
          histoV_PDF_Fit[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoV_PDF_Fit[c][i][j]->SetBinContent(bj+1, pbjets_v_ext[c][i][j]->getVal());
          histoV_PDF_Fit[c][i][j]->SetBinError(bj+1, pbjets_v_ext[c][i][j]->getPropagatedError(*(fit_result[c])));
          histoVbb_PDF_Fit[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoVbb_PDF_Fit[c][i][j]->SetBinContent(bj+1, pbjets_vb[c][i][j]->getVal());
            //histoVbb_PDF_Fit[c][i][j]->SetBinError(bj+1, pbjets_v_ext[i][j]->getPropagatedError(*fit_result[c]));
          histoST_PDF_Fit[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoST_PDF_Fit[c][i][j]->SetBinContent(bj+1, pbjets_st[c][i][j]->getVal());
            //histoST_PDF_Fit[c][i][j]->SetBinError(bj+1, pbjets_v_ext[i][j]->getPropagatedError(*fit_result[c]));
          histoQCD_PDF_Fit[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoQCD_PDF_Fit[c][i][j]->SetBinContent(bj+1, pbjets_qcd[c][i][j]->getVal());
            //histoVbb_PDF_Fit[c][i][j]->SetBinError(bj+1, pbjets_v_ext[i][j]->getPropagatedError(*fit_result[c]));
          histoVV_PDF_Fit[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoVV_PDF_Fit[c][i][j]->SetBinContent(bj+1, pbjets_vv[c][i][j]->getVal());
            //histoST_PDF_Fit[c][i][j]->SetBinError(bj+1, pbjets_v_ext[i][j]->getPropagatedError(*fit_result[c]));
          histoTot_PDF_Fit[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoTot_PDF_Fit[c][i][j]->SetBinContent(bj+1, model[c][i][j]->getVal(ras));
          histoTot_PDF_Fit[c][i][j]->SetBinError(bj+1, model[c][i][j]->getPropagatedError(*(fit_result[c])));
          histoTotC_PDF_Fit[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoTotC_PDF_Fit[c][i][j]->SetBinContent(bj+1, model_constraint[c][i][j]->getVal(ras));
          histoTotC_PDF_Fit[c][i][j]->SetBinError(bj+1, model_constraint[c][i][j]->getPropagatedError(*(fit_result[c])));
        }
        
        for (UInt_t bj=0; bj<6; bj++) {
          if (bj==0 || bj==1) {
            sprintf(label, "N%dbjet", bj);
          } else {
            sprintf(label, "N%dbjets", bj);
          }
          nbjets.setLabel(label); 
          histoTotS_PDF_Fit[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoTotS_PDF_Fit[c][i][j]->SetBinContent(bj+1, pbjets_tt_ext[c][i][j]->getVal()*pbjets_tt_ext[c][i][j]->expectedEvents(NULL)/histoTT_PDF_Fit[c][i][j]->Integral()
                                                    +pbjets_v_ext[c][i][j]->getVal()*pbjets_v_ext[c][i][j]->expectedEvents(NULL)/histoV_PDF_Fit[c][i][j]->Integral()
                                                    +Nvb[c][i]->getVal()*pbjets_vb[c][i][j]->getVal()/histoVbb_PDF_Fit[c][i][j]->Integral()
                                                    +Nst[c][i]->getVal()*pbjets_st[c][i][j]->getVal()/histoST_PDF_Fit[c][i][j]->Integral()
                                                    /*                                                    +Nqcd[c][i]->getVal()*pbjets_qcd[c][i][j]->getVal()/histoQCD_PDF_Fit[c][i][j]->Integral() */
                                                    +Nvv[c][i]->getVal()*pbjets_vv[c][i][j]->getVal()/histoVV_PDF_Fit[c][i][j]->Integral());
            //histoTotS_PDF_Fit[i][j]->SetBinError(2, model_constraint[i][j]->getPropagatedError(*fit_result));
        }
        
        
        for (int bj=0; bj<6; bj++) {
          if (bj==0 || bj==1) {
            sprintf(label, "N%dbjet", bj);
          } else {
            sprintf(label, "N%dbjets", bj);
          }
          nbjets.setLabel(label); 
          histoTotS_PDF_Fit[c][i][j]->SetBinError(bj+1, model_constraint[c][i][j]->getPropagatedError(*(fit_result[c])) * histoTotS_PDF_Fit[c][i][j]->Integral()/histoTot_PDF_Fit[c][i][j]->Integral());
        }
        
        
        histoTT_PDF_Fit[c][i][j]->Scale(pbjets_tt_ext[c][i][j]->expectedEvents(NULL)/histoTT_PDF_Fit[c][i][j]->Integral());
        histoV_PDF_Fit[c][i][j]->Scale(pbjets_v_ext[c][i][j]->expectedEvents(NULL)/histoV_PDF_Fit[c][i][j]->Integral());
        histoTot_PDF_Fit[c][i][j]->Scale(model[c][i][j]->expectedEvents(NULL)/histoTot_PDF_Fit[c][i][j]->Integral());
        histoTotC_PDF_Fit[c][i][j]->Scale(model_constraint[c][i][j]->expectedEvents(NULL)/histoTotC_PDF_Fit[c][i][j]->Integral());
          //      histoTotS_PDF_Fit[i][j]->Scale(model_constraint[i][j]->expectedEvents(NULL)/histoTotC_PDF_Fit[i][j]->Integral());
        
        /*   
         Double_t integral = 0.;
         integral = histoTT_PDF_Fit[i][j]->Integral();  histoTT_PDF_Fit[i][j]->Scale(Ntt[i]->getVal()/integral);
         integral = histoV_PDF_Fit[i][j]->Integral();   histoV_PDF_Fit[i][j]->Scale(Nv[i]->getVal()/integral);
         integral = histoTot_PDF_Fit[i][j]->Integral(); if (Nst[i]!=NULL && Nvb[i]!=NULL) {
         histoTot_PDF_Fit[i][j]->Scale((Ntt[i]->getVal()+Nv[i]->getVal()+Nst[i]->getVal()+Nvb[i]->getVal())/integral);
         } else {
         histoTot_PDF_Fit[i][j]->Scale((Ntt[i]->getVal()+Nv[i]->getVal())/integral);
         }
         */
          //      histoTT_PDF_Fit[i][j]->Print();
        printf(" PDF Fit TT (norm : %lf ; expectedEvents : %lf) : ", pbjets_tt_ext[c][i][j]->getNorm(), pbjets_tt_ext[c][i][j]->expectedEvents(NULL));
        for (int k=1; k<=4 ; k++) {
          printf(" %lf ", histoTT_PDF_Fit[c][i][j]->GetBinContent(k));
        }
        printf("\nNtt[i] : %lf\n",Ntt[c][i]->getVal());
          //      histoV_PDF_Fit[i][j]->Print();
        printf(" PDF Fit  V (norm : %lf ; expectedEvents : %lf): ", pbjets_v_ext[c][i][j]->getNorm(), pbjets_v_ext[c][i][j]->expectedEvents(NULL));
        for (int k=1; k<=4 ; k++) {
          printf(" %lf ", histoV_PDF_Fit[c][i][j]->GetBinContent(k));
        }
        printf("\nNv[i] : %lf\n",Nv[c][i]->getVal());
          //      histoTot_PDF_Fit[i][j]->Print();
        printf(" PDF Fit Tot (norm : %lf ; expectedEvents : %lf) : ", model[c][i][j]->getNorm(), model[c][i][j]->expectedEvents(NULL));
        for (int k=1; k<=4 ; k++) {
          printf(" %lf ", histoTot_PDF_Fit[c][i][j]->GetBinContent(k));
        }
        printf("\n");
        printf(" PDF Fit Tot constraint (norm : %lf ; expectedEvents : %lf) : ", model_constraint[c][i][j]->getNorm(), model_constraint[c][i][j]->expectedEvents(NULL));
        for (int k=1; k<=4 ; k++) {
          printf(" %lf ", histoTotC_PDF_Fit[c][i][j]->GetBinContent(k));
        }
        printf("\n");
        printf(" PDF Fit Tot summed (norm : %lf ; expectedEvents : %lf) : ", -1., -1.);
        for (int k=1; k<=4 ; k++) {
          printf(" %lf ", histoTotS_PDF_Fit[c][i][j]->GetBinContent(k));
        }
        printf("\n");
      }
      delete ras;
      
    }
  }
  
  
  double** save_tt = new double*[NbOfChannels];
  double** save_v  = new double*[NbOfChannels];
  double** save_st = new double*[NbOfChannels];
  double** save_vb = new double*[NbOfChannels];
  double** save_qcd = new double*[NbOfChannels];
  double** save_vv = new double*[NbOfChannels];
  
  double** save_eb = new double*[NbOfChannels];
  double** save_et = new double*[NbOfChannels];
  double** save_ev = new double*[NbOfChannels];
  
  printf("\n\nSummary of the fit results : \n\n\n");
 
  {
    /* Fit output to root file for DD-reweighting */
    TFile *outFile = new TFile (("FitOutput_"+rootFileName).c_str(), "RECREATE");
    outFile->cd();


    //       samples_names[0].push_back("TTJets"/*+CChannel[channel]*/);
    //    samples_labelBQ[0].push_back("TTJetsInc"/*+CChannel[channel]*/);
    //    files[0].push_back(TFile::Open((leg3Dir+"TTJets"/*+CChannel[channel]*/+".root").c_str(), "READ"));

    std::vector<std::vector<TH2D*> > vOutHistRoofit(NbOfChannels, std::vector<TH2D*>() ) ;
    //std::vector<std::vector<TH2D*> > vOutHistClass(NbOfChannels, std::vector<TH2D*>() ) ;

    for(UInt_t c=0; c<NbOfChannels; c++) {
      Int_t idx = 0;
      //      if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
      //      if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
      //std::vector<TH2D*> tmpVecHist_class;
      std::vector<TH2D*> tmpVecHist;
      for(std::vector<std::vector<std::string> >::iterator iter = samples_labelBQ.begin() ; iter!=samples_labelBQ.end() ; iter++) {
        std::string name;
        if (idx==0) {
          name = "TT-like" ;
        } else if (idx==1) {
          name = "V-like" ;
        } else if (idx==2) {
          name = "V-bb" ;
        } else if (idx==3) {
          name = "Single_top" ;
        } else if (idx==4) {
          name = "VV" ;
        } else {
          name = "!!!!  non supported catogory !!!!!" ;
          printf("!!!!  non supported catogory !!!!! : %d !!!!!\n", idx);
          continue ;
        }
        char xaxisName[50];
        char yaxisName[50];

        //TH2D* classHist = new TH2D( (name+channelSuffix[c]+"_fitOutput_classMethods").c_str(), (name+"_fitOutput_classMethods").c_str(), 6, 0., 6., NbOfJetBins, 0., NbOfJetBins);
        TH2D* roofitHist = new TH2D( (name+channelSuffix[c]+"_fitOutput").c_str(), (name+"_fitOutput").c_str(), 6, 0., 6., NbOfJetBins, 0., NbOfJetBins);
        for (UInt_t ix=0 ; ix<6 ; ix++) {
          sprintf(xaxisName, "%d b-tags", (Int_t) ix);
          //classHist->GetXaxis()->SetBinLabel(ix+1, xaxisName);
          roofitHist->GetXaxis()->SetBinLabel(ix+1, xaxisName);
        }
        for (UInt_t iy=0 ; iy<NbOfJetBins ; iy++) {
          if (iy==NbOfJetBins-1) { 
            sprintf(yaxisName, "geq %d jets", MinNbOfJets+iy);
          } else {
            sprintf(yaxisName, "%d jets", MinNbOfJets+iy);
          }
          //classHist->GetYaxis()->SetBinLabel(iy+1,yaxisName);
          roofitHist->GetYaxis()->SetBinLabel(iy+1,yaxisName);
        }
        //tmpVecHist_class.push_back(classHist);
        tmpVecHist.push_back(roofitHist);
        idx++ ;
      }
      //vOutHistClass[c] = tmpVecHist_class ;   
      vOutHistRoofit[c] = tmpVecHist ;   
    }
    //for (std::vector<std::vector<TH2D*> >::iterator iter=vOutHistClass.begin(); iter!=vOutHistClass.end(); iter++) {
    //for (std::vector<TH2D*>::iterator it=iter->begin(); it!=iter->end(); it++) {
    //printf(" --Class--> %s\n", (*it)->GetName());
    //}
    //}
    for (std::vector<std::vector<TH2D*> >::iterator iter=vOutHistRoofit.begin(); iter!=vOutHistRoofit.end(); iter++) {
      for (std::vector<TH2D*>::iterator it=iter->begin(); it!=iter->end(); it++) {
        printf(" --Roofit-> %s\n", (*it)->GetName());
      }
    }
    /* Printing the output of the data-fitted model (if isOnData) and MC input  */

    for(UInt_t c=0; c<NbOfChannels; c++) { 
      if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
      if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
      printf("\nOn MC (after the fit if you are really patient): \n");
      printf("For channel : %s\n", channelSuffix[c]);
      printf("\n\\begin{tabular}{r|ccccc|c|c}\n");
      for (UInt_t j=0; j<NbOfWP; j++) {
        if (NbOfWP!=1) {
          printf("\\multicolumn{8}{l}{For working point \\verb|%s| :}\\\\ \n", wpSuffix[j]);
          printf("\\hline");
        }
        for (UInt_t i=0; i<NbOfJetBins; i++) {
          std::string jts = jetSuffix[i];
          if (jts.find("_")!=std::string::npos) {
            jts=jts.substr(jts.find("_")+1,std::string::npos);
          }
          if (jts.find("Inc")!=std::string::npos) {
            jts = " \\geq "+jts.substr(0,jts.find("Inc"))+jts.substr(jts.find("Inc")+3,std::string::npos);
          } else if (jts.find("Exc")!=std::string::npos) {
            jts = " = "+jts.substr(0,jts.find("Exc"))+jts.substr(jts.find("Exc")+3,std::string::npos);
          }
          printf("\\multicolumn{8}{l}{For jet bin $ %s $ :}\\\\ \n", jts.c_str());
          printf("$\\#b-tags$ & $t\\bar{t}-like$ & $V-like$ & Single Top & $W-bb$ & $VV$ & total (MC) & total (Data) \\\\ \\hline \n");
          for (UInt_t k=0; k<6; k++) { //because up to 5 b-tags in the equations
            printf("$ %d $ & $ %.0lf $ & $ %.0lf $ & $ %.0lf $ & $ %.0lf $ & $ %.0lf $ & $ %.0lf $ & $ %.0lf $ \\\\ \n", k,
                   hNbtaggedJets_ttlike[c][i][j]->GetBinContent(hNbtaggedJets_ttlike[c][i][j]->FindBin(k)),
                   hNbtaggedJets_vlike[c][i][j]->GetBinContent(hNbtaggedJets_vlike[c][i][j]->FindBin(k)),
                   hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k)),
                   hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k)),
                   hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k)),
                   hNbtaggedJets_ttlike[c][i][j]->GetBinContent(hNbtaggedJets_ttlike[c][i][j]->FindBin(k))+
                   hNbtaggedJets_vlike[c][i][j]->GetBinContent(hNbtaggedJets_vlike[c][i][j]->FindBin(k))+
                   hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))+
                   hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))+
                   hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k)),
                   hNbtaggedJets[c][i][j]->GetBinContent(hNbtaggedJets[c][i][j]->FindBin(k))
                   );
          }
        }
      }
      printf("\\end{tabular}\n\n");
    }
    for(UInt_t c=0; c<NbOfChannels; c++) { 
      if (!semiMuon && (((c+1)&muonMask)!=0)) continue ;
      if (!semiElectron && (((c+1)&electronMask)!=0)) continue ;
      printf("\nOn the fit : \n");
      printf("For channel : %s\n", channelSuffix[c]);
      printf("\n\\begin{tabular}{r|ccccc|c|c}\n");
      RooArgList fVal = fit_result[c]->floatParsFinal() ;
      VJetEstimation vj;
      std::vector<Double_t> va, vb, vc, vd;
      for (UInt_t k=0; k<MinNbOfJets+NbOfJetBins; k++) { va.push_back((k<MinNbOfJets?0.:e0bq[c][k-MinNbOfJets]->getVal())); /*printf("%lf ",va[k]);*/ }; //printf("\n");
      for (UInt_t k=0; k<MinNbOfJets+NbOfJetBins; k++) { vb.push_back((k<MinNbOfJets?0.:e1bq[c][k-MinNbOfJets]->getVal())); /*printf("%lf ",vb[k]);*/ }; //printf("\n");
      for (UInt_t k=0; k<MinNbOfJets+NbOfJetBins; k++) { vc.push_back((k<MinNbOfJets?0.:e2bq[c][k-MinNbOfJets]->getVal())); /*printf("%lf ",vc[k]);*/ }; //printf("\n");
      for (UInt_t k=0; k<MinNbOfJets+NbOfJetBins; k++) { vd.push_back((k<MinNbOfJets?0.:e3bq[c][k-MinNbOfJets]->getVal())); /*printf("%lf ",vd[k]);*/ }; //printf("\n");
      vj.SetTTEffbq(va, vb, vc, vd);
      for (UInt_t j=0; j<NbOfWP; j++) {
        if (NbOfWP!=1) {
          printf("\\multicolumn{8}{l}{For working point \\verb|%s| :}\\\\ \n", wpSuffix[j]);
          printf("\\hline");
        }
        for (UInt_t i=0; i<NbOfJetBins; i++) {
                     std::string jts = jetSuffix[i];
          if (jts.find("_")!=std::string::npos) {
            jts=jts.substr(jts.find("_")+1,std::string::npos);
          }
          if (jts.find("Inc")!=std::string::npos) {
            jts = " \\geq "+jts.substr(0,jts.find("Inc"))+jts.substr(jts.find("Inc")+3,std::string::npos);
          } else if (jts.find("Exc")!=std::string::npos) {
            jts = " = "+jts.substr(0,jts.find("Exc"))+jts.substr(jts.find("Exc")+3,std::string::npos);
          }
          printf("\\multicolumn{8}{l}{For jet bin $ %s $ :}\\\\ \n", jts.c_str());
          printf("$\\#b-tags$ & $t\\bar{t}-like$ & $V-like$ & Single Top & $W-bb$ & $VV$ & total (est.) & total (Data) \\\\ \\hline \n");
          for (UInt_t k=0; k<6; k++) { //because up to 5 b-tags in the equations

                        Double_t __ntt = ((RooRealVar*) fVal.find((std::string()+"Ntt"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal() ;
                        Double_t __ntt_err = ((RooRealVar*) fVal.find((std::string()+"Ntt"+channelSuffix[c]+jtSuf[i]).c_str()))->getError() ;
                        Double_t __eb = ((RooRealVar*) fVal.find((std::string()+"eb"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal() ;
                        Double_t __eb_err = ((RooRealVar*) fVal.find((std::string()+"eb"+channelSuffix[c]+wpSuf[j]).c_str()))->getError() ;
                        Double_t __eudsc = ((RooRealVar*) fVal.find((std::string()+"eudsc"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal() ;
                        Double_t __eudsc_err = ((RooRealVar*) fVal.find((std::string()+"eudsc"+channelSuffix[c]+wpSuf[j]).c_str()))->getError() ;
printf("%%%%%%  Ntt :    indirect (fit_result) %lf \\pm %lf \t direct (variable itself) %lf \\pm %lf\n", __ntt, __ntt_err, Ntt[c][i]->getVal(), Ntt[c][i]->getError());
printf("%%%%%%  e_b :    indirect (fit_result) %lf \\pm %lf \t direct (variable itself) %lf \\pm %lf\n", __eb, __eb_err, eb[c][j]->getVal(), eb[c][j]->getError());
printf("%%%%%%  e_udsc : indirect (fit_result) %lf \\pm %lf \t direct (variable itself) %lf \\pm %lf\n", __eudsc, __eudsc_err, eudsc[c][j]->getVal(), eudsc[c][j]->getError());

            RooFormulaVar *tmptt = new RooFormulaVar((std::string()+"tmp").c_str(), (std::string()+"tmp").c_str(), "@0*@1/(@2+@3+@4+@5+@6+@7)", RooArgList(*(Ntt[c][i]), *((k==0 ? p0bjets_tt : (k==1 ? p1bjets_tt : (k==2 ? p2bjets_tt : (k==3 ? p3bjets_tt : (k==4 ? p4bjets_tt : p5bjets_tt)))))[c][i][j]), *(p0bjets_tt[c][i][j]), *(p1bjets_tt[c][i][j]), *(p2bjets_tt[c][i][j]), *(p3bjets_tt[c][i][j]), *(p4bjets_tt[c][i][j]), *(p5bjets_tt[c][i][j])));
            vOutHistRoofit[c][0]->SetBinContent(k+1, i+1, tmptt->getVal());
            vOutHistRoofit[c][0]->SetBinError(k+1, i+1, tmptt->getPropagatedError(*(fit_result[c])));

            RooFormulaVar *tmpv = new RooFormulaVar((std::string()+"vtmp").c_str(), (std::string()+"vtmp").c_str(), "@0*@1/(@2+@3+@4+@5+@6+@7)", RooArgList(*(Nv[c][i]), *((k==0 ? p0bjets_v : (k==1 ? p1bjets_v : (k==2 ? p2bjets_v : (k==3 ? p3bjets_v : (k==4 ? p4bjets_v : p5bjets_v)))))[c][i][j]), *(p0bjets_v[c][i][j]), *(p1bjets_v[c][i][j]), *(p2bjets_v[c][i][j]), *(p3bjets_v[c][i][j]), *(p4bjets_v[c][i][j]), *(p5bjets_v[c][i][j])));
            vOutHistRoofit[c][1]->SetBinContent(k+1, i+1, tmpv->getVal());
            vOutHistRoofit[c][1]->SetBinError(k+1, i+1, tmpv->getPropagatedError(*(fit_result[c])));

            //            RooFormulaVar * tmptot = new RooFormulaVar((std::string()+"tottmp").c_str(), (std::string()+"tottmp").c_str(), "@0+@1+@2*@3+@4*@5+@6*@7", RooArgList(*tmptt, *tmpv,
            //                                                                                                                                                                 *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str())), RooConstVar("p_st_k","p_st_k",(hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())),
            //                                                                                                                                                                 *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str())), RooConstVar("p_vv_k","p_vv_k",(hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())),
            //                                                                                                                                                                 *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str())), RooConstVar("p_vb_k","p_vb_k",(hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral()))
            //                                                                                                                                                                 ));
            RooRealVar p_st_k("p_st_k","p_st_k",(hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())) ;
            RooRealVar p_vv_k("p_vv_k","p_vv_k",(hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())) ;
            RooRealVar p_vb_k("p_vb_k","p_vb_k",(hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())) ;
            p_st_k.setConstant();
            p_vv_k.setConstant();
            p_vb_k.setConstant();
            RooFormulaVar * tmptot = new RooFormulaVar((std::string()+"tottmp").c_str(), (std::string()+"tottmp").c_str(), "@0+@1+@2*@3+@4*@5+@6*@7", RooArgList( *tmptt, *tmpv,
                                                                                                                                                                  *(Nst[c][i]), p_st_k,
                                                                                                                                                                  *(Nvv[c][i]), p_vv_k,
                                                                                                                                                                  *(Nvb[c][i]), p_vb_k
                                                                                                                                                                  ));
            printf("%d & $ ", k);
            printfValErr(tmptt->getVal(), tmptt->getPropagatedError(*(fit_result[c])));
            printf(" $ & $ ");
            printfValErr(tmpv->getVal(), tmpv->getPropagatedError(*(fit_result[c])));
            printf(" $ & $ ");
            printfValErr( (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                          *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                          (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                          *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getError());
            printf(" $ & $ ");
            printfValErr( (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                          *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                          (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                          *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getError());
            printf(" $ & $ ");
            printfValErr( (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                          *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                          (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                          *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError() );
            printf(" $ & $ ");
            printfValErr(tmptot->getVal(), tmptot->getPropagatedError(*(fit_result[c])));
            printf(" $ & $ ");
            printf(" %.0lf ", hNbtaggedJets[c][i][j]->GetBinContent(hNbtaggedJets[c][i][j]->FindBin(k)));
            printf(" $ \\\\ \n");

            printf("%% %d & $ %.0lf \\pm %.0lf $ & $ %.0lf \\pm %.0lf $ & $ %.0lf \\pm %.0lf $ & $ %.0lf \\pm %.0lf $ & $ %.0lf \\pm %.0lf $  \\\\ \n", k,
                   /*(k==0 ? vj.Ntt_0bjet(__ntt,__eb,__eudsc, MinNbOfJets+i) :
                     (k==1) ? vj.Ntt_1bjet(__ntt,__eb,__eudsc, MinNbOfJets+i) :
                     vj.Ntt_2bjets(__ntt,__eb,__eudsc, MinNbOfJets+i)),
                     (k==0 ? vj.Ntt_err_0bjet(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
                     (k==1 ? vj.Ntt_err_1bjet(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
                     vj.Ntt_err_2bjets(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i))),
                     (k==0 ? VJetEstimation::Nv_0bjet : (k==1 ? VJetEstimation::Nv_1bjet : VJetEstimation::Nv_2bjets))
                     (((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal(),
                     MinNbOfJets+i),
                     (k==0 ? VJetEstimation::Nv_err_0bjet : (k==1 ? VJetEstimation::Nv_err_1bjet : VJetEstimation::Nv_err_2bjets))
                     (((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getError(),
                     MinNbOfJets+i),*/
                   tmptt->getVal(), tmptt->getPropagatedError(*(fit_result[c])),
                   tmpv->getVal(), tmpv->getPropagatedError(*(fit_result[c])),
                   /* */
                   (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                   (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                   (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError()
                   );

            printf("%% %d & $ %lf \\pm %lf $ & $ %lf \\pm %lf $ & $ %lf \\pm %lf $ & $ %lf \\pm %lf $ & $ %lf \\pm %lf $  \\\\ \n", k,
                   /*(k==0 ? vj.Ntt_0bjet(__ntt,__eb,__eudsc, MinNbOfJets+i) :
                     (k==1) ? vj.Ntt_1bjet(__ntt,__eb,__eudsc, MinNbOfJets+i) :
                     vj.Ntt_2bjets(__ntt,__eb,__eudsc, MinNbOfJets+i)),
                     (k==0 ? vj.Ntt_err_0bjet(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
                     (k==1 ? vj.Ntt_err_1bjet(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
                     vj.Ntt_err_2bjets(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i))),
                     (k==0 ? VJetEstimation::Nv_0bjet : (k==1 ? VJetEstimation::Nv_1bjet : VJetEstimation::Nv_2bjets))
                     (((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal(),
                     MinNbOfJets+i),
                     (k==0 ? VJetEstimation::Nv_err_0bjet : (k==1 ? VJetEstimation::Nv_err_1bjet : VJetEstimation::Nv_err_2bjets))
                     (((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getError(),
                     MinNbOfJets+i),*/
                   tmptt->getVal(), tmptt->getPropagatedError(*(fit_result[c])),
                   tmpv->getVal(), tmpv->getPropagatedError(*(fit_result[c])),
                   /* */
                   (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                   (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                   (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError()
                   );

             printf("%%%% %d & $ %lf \\pm %lf $ & $ %lf \\pm %lf $ & $ %lf \\pm %lf $ & $ %lf \\pm %lf $ & $ %lf \\pm %lf $  \\\\ \n", k,
                    (k==0 ? vj.Ntt_0bjet(__ntt,__eb,__eudsc, MinNbOfJets+i) :
                     (k==1 ? vj.Ntt_1bjet(__ntt,__eb,__eudsc, MinNbOfJets+i) :
                      (k==2 ? vj.Ntt_2bjets(__ntt,__eb,__eudsc, MinNbOfJets+i) : 
                       vj.Ntt_3bjets(__ntt,__eb,__eudsc, MinNbOfJets+i)))),
                    (k==0 ? vj.Ntt_err_0bjet(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
                     (k==1 ? vj.Ntt_err_1bjet(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
                      (k==2 ? vj.Ntt_err_2bjets(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
                       vj.Ntt_err_3bjets(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i)))),
                     (k==0 ? VJetEstimation::Nv_0bjet : (k==1 ? VJetEstimation::Nv_1bjet : (k==2 ? VJetEstimation::Nv_2bjets : VJetEstimation::Nv_3bjets)))
                     (((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal(),
                     MinNbOfJets+i),
                     (k==0 ? VJetEstimation::Nv_err_0bjet : (k==1 ? VJetEstimation::Nv_err_1bjet : (k==2 ? VJetEstimation::Nv_err_2bjets : VJetEstimation::Nv_err_3bjets)))
                     (((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getError(),
                     MinNbOfJets+i),/*
                   tmptt->getVal(), tmptt->getPropagatedError(*(fit_result[c])),
                   tmpv->getVal(), tmpv->getPropagatedError(*(fit_result[c])),
                    */
                   (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                   (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                   (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError()
                   );

             printf("%%%%%% %d & $ %lf \\pm %lf $ & $ %lf \\pm %lf $ & $ %lf \\pm %lf $ & $ %lf \\pm %lf $ & $ %lf \\pm %lf $  \\\\ \n", k,
                    (k==0 ? vj.Ntt_0bjet_e3(__ntt,__eb,__eudsc, MinNbOfJets+i) :
                     (k==1 ? vj.Ntt_1bjet_e3(__ntt,__eb,__eudsc, MinNbOfJets+i) :
                      (k==2 ? vj.Ntt_2bjets_e3(__ntt,__eb,__eudsc, MinNbOfJets+i) : 
                       vj.Ntt_3bjets_e3(__ntt,__eb,__eudsc, MinNbOfJets+i)))),
                    (k==0 ? vj.Ntt_err_0bjet(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
                     (k==1 ? vj.Ntt_err_1bjet(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
                      (k==2 ? vj.Ntt_err_2bjets(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
                       vj.Ntt_err_3bjets(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i)))),
                     (k==0 ? VJetEstimation::Nv_0bjet : (k==1 ? VJetEstimation::Nv_1bjet : (k==2 ? VJetEstimation::Nv_2bjets : VJetEstimation::Nv_3bjets)))
                     (((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal(),
                     MinNbOfJets+i),
                     (k==0 ? VJetEstimation::Nv_err_0bjet : (k==1 ? VJetEstimation::Nv_err_1bjet : (k==2 ? VJetEstimation::Nv_err_2bjets : VJetEstimation::Nv_err_3bjets)))
                     (((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal(),
                     ((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getError(),
                     MinNbOfJets+i),/*
                   tmptt->getVal(), tmptt->getPropagatedError(*(fit_result[c])),
                   tmpv->getVal(), tmpv->getPropagatedError(*(fit_result[c])),
                    */
                   (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                   (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
                   (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
                   (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                   *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError()
                   );

           // vOutHistClass[c][0]->SetBinContent(k+1, i+1, (k==0 ? vj.Ntt_0bjet(__ntt,__eb,__eudsc, MinNbOfJets+i) :
            //(k==1) ? vj.Ntt_1bjet(__ntt,__eb,__eudsc, MinNbOfJets+i) :
            //vj.Ntt_2bjets(__ntt,__eb,__eudsc, MinNbOfJets+i)));
            //vOutHistClass[c][0]->SetBinError(k+1, i+1, (k==0 ? vj.Ntt_err_0bjet(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
            //(k==1 ? vj.Ntt_err_1bjet(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i) :
            //vj.Ntt_err_2bjets(__ntt, __ntt_err, __eb, __eb_err, __eudsc, __eudsc_err,  MinNbOfJets+i))));

            //vOutHistClass[c][1]->SetBinContent(k+1, i+1, (k==0 ? VJetEstimation::Nv_0bjet : (k==1 ? VJetEstimation::Nv_1bjet : VJetEstimation::Nv_2bjets))
            //(((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
            //((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal(),
            //MinNbOfJets+i));
            //vOutHistClass[c][1]->SetBinError(k+1, i+1, (k==0 ? VJetEstimation::Nv_err_0bjet : (k==1 ? VJetEstimation::Nv_err_1bjet : VJetEstimation::Nv_err_2bjets))
            //(((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal(),
            //((RooRealVar*) fVal.find((std::string()+"Nv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError(),
            //((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getVal(),
            //((RooRealVar*) fVal.find((std::string()+"euds"+channelSuffix[c]+wpSuf[j]).c_str()))->getError(),
            //MinNbOfJets+i));
 
 
            //vOutHistClass[c][2]->SetBinContent(k+1, i+1, (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
            //*((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal());
            //vOutHistClass[c][2]->SetBinError(k+1, i+1, (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
            //*((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getError());
            //vOutHistClass[c][3]->SetBinContent(k+1, i+1, (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
            //*((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal());
            //vOutHistClass[c][3]->SetBinError(k+1, i+1, (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
            //*((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getError());
            //vOutHistClass[c][4]->SetBinContent(k+1, i+1, (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
            //*((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal());
            //vOutHistClass[c][4]->SetBinError(k+1, i+1, (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
            //*((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError());

            vOutHistRoofit[c][2]->SetBinContent(k+1, i+1, (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                                                *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal());
            vOutHistRoofit[c][2]->SetBinError(k+1, i+1, (hNbtaggedJets_vbjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(k))/hNbtaggedJets_vbjets[c][i][j]->Integral())
                                              *((RooRealVar*) fVal.find((std::string()+"Nvb"+channelSuffix[c]+jtSuf[i]).c_str()))->getError());
            vOutHistRoofit[c][3]->SetBinContent(k+1, i+1, (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                                                *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal());
            vOutHistRoofit[c][3]->SetBinError(k+1, i+1, (hNbtaggedJets_stjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(k))/hNbtaggedJets_stjets[c][i][j]->Integral())
                                              *((RooRealVar*) fVal.find((std::string()+"Nst"+channelSuffix[c]+jtSuf[i]).c_str()))->getError());
            vOutHistRoofit[c][4]->SetBinContent(k+1, i+1, (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                                                *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getVal());
            vOutHistRoofit[c][4]->SetBinError(k+1, i+1, (hNbtaggedJets_vvjets[c][i][j]->Integral()==0.?0.:hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(k))/hNbtaggedJets_vvjets[c][i][j]->Integral())
                                              *((RooRealVar*) fVal.find((std::string()+"Nvv"+channelSuffix[c]+jtSuf[i]).c_str()))->getError());

            tmptot->Delete();
          }
        }
        printf("\\end{tabular}\n\n");
      }
    }
    outFile->Write();
    outFile->Close();

  }


  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  
  fout->cd();
  
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    
    save_tt[c] = new double[NbOfJetBins];
    save_v [c] = new double[NbOfJetBins];
    save_st[c] = new double[NbOfJetBins];
    save_vb[c] = new double[NbOfJetBins];
    save_qcd[c] = new double[NbOfJetBins];
    save_vv[c] = new double[NbOfJetBins];
    
    save_eb[c] = new double[NbOfWP];
    save_et[c] = new double[NbOfWP];
    save_ev[c] = new double[NbOfWP];
    
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    printf("\n\nFor channel #%d\n\n", c);
    if(fit_result[c]) {
      fit_result[c]->Print("v");
      fit_result[c]->correlationMatrix().Print() ;
      fit_result[c]->Write((std::string()+"FitResult"+channelSuffix[c]).c_str());
    }
      //Save the fit values (to reload them for PE)
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      save_tt[c][i] = Ntt[c][i]->getVal();
      save_v [c][i] = Nv[c][i]->getVal();
      save_st[c][i] = Nst[c][i]->getVal();
      save_vb[c][i] = Nvb[c][i]->getVal();
      save_qcd[c][i] = Nqcd[c][i]->getVal();
      save_vv[c][i] = Nvv[c][i]->getVal();
    }
    for (UInt_t wp=0; wp<NbOfWP; wp++) {
      save_eb[c][wp] = eb[c][wp]->getVal();
      save_et[c][wp] = eudsc[c][wp]->getVal();
      save_ev[c][wp] = euds[c][wp]->getVal();
    }
      //Reset to MC values (for the plot)    
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      hMC_Nominal_N[c][i]->SetBinContent(1, Nttlike[c][i]);   Ntt[c][i]->setVal(Nttlike[c][i]);
      hMC_Nominal_N[c][i]->SetBinContent(2, Nvlike[c][i]);    Nv[c][i]->setVal(Nvlike[c][i]);
      hMC_Nominal_N[c][i]->SetBinContent(3, Nstjets[c][i]);   Nst[c][i]->setVal(Nstjets[c][i]);
      hMC_Nominal_N[c][i]->SetBinContent(4, Nvbjets[c][i]);   Nvb[c][i]->setVal(Nvbjets[c][i]);
      hMC_Nominal_N[c][i]->SetBinContent(5, Nqcd_mean[c][i]); Nqcd[c][i]->setVal(Nqcd_mean[c][i]);
      hMC_Nominal_N[c][i]->SetBinContent(6, Nvvjets[c][i]);   Nvv[c][i]->setVal(Nvvjets[c][i]);
        //      hMC_Nominal_N[c][i]->SetBinContent(7, n[c][i]);   n[c][i]->setVal(Nvvjets[c][i]);
    }
    for (UInt_t wp=0; wp<NbOfWP; wp++) {
      hMC_Nominal_E[c][wp]->SetBinContent(1, init_Eb[c][wp]);    eb[c][wp]->setVal(init_Eb[c][wp]);
      hMC_Nominal_E[c][wp]->SetBinContent(2, init_Eudsc[c][wp]); eudsc[c][wp]->setVal(init_Eudsc[c][wp]);
      hMC_Nominal_E[c][wp]->SetBinContent(3, init_Euds[c][wp]);  euds[c][wp]->setVal(init_Euds[c][wp]);
    }
  }
  
  TH1F ****histoTT_PDF_MC = new TH1F***[NbOfChannels];
  TH1F ****histoV_PDF_MC = new TH1F***[NbOfChannels];
  TH1F ****histoVbb_PDF_MC = new TH1F***[NbOfChannels];
  TH1F ****histoST_PDF_MC = new TH1F***[NbOfChannels];
  TH1F ****histoQCD_PDF_MC = new TH1F***[NbOfChannels];
  TH1F ****histoVV_PDF_MC = new TH1F***[NbOfChannels];
  TH1F ****histoTot_PDF_MC = new TH1F***[NbOfChannels];
  TH1F ****histoTotC_PDF_MC = new TH1F***[NbOfChannels];
  TH1F ****histoTotS_PDF_MC = new TH1F***[NbOfChannels];
  for (UInt_t c=0; c<NbOfChannels; c++) {
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    histoTT_PDF_MC[c] = new TH1F**[NbOfJetBins];
    histoV_PDF_MC[c] = new TH1F**[NbOfJetBins];
    histoVbb_PDF_MC[c] = new TH1F**[NbOfJetBins];
    histoST_PDF_MC[c] = new TH1F**[NbOfJetBins];
    histoQCD_PDF_MC[c] = new TH1F**[NbOfJetBins];
    histoVV_PDF_MC[c] = new TH1F**[NbOfJetBins];
    histoTot_PDF_MC[c] = new TH1F**[NbOfJetBins];
    histoTotC_PDF_MC[c] = new TH1F**[NbOfJetBins];
    histoTotS_PDF_MC[c] = new TH1F**[NbOfJetBins];
    for(UInt_t i=0; i<NbOfJetBins ; i++) {
      histoTT_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoV_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoVbb_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoST_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoQCD_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoVV_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoTot_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoTotC_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoTotS_PDF_MC[c][i] = new TH1F*[NbOfWP];
      for (UInt_t j=0; j<NbOfWP ; j++) {
        histoTT_PDF_MC[c][i][j] = new TH1F((std::string()+"histoTT_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+channelSuffix[c]+"histoTT_PDF_MC"+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoV_PDF_MC[c][i][j] = new TH1F((std::string()+"histoV_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+channelSuffix[c]+"histoV_PDF_MC"+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoVbb_PDF_MC[c][i][j] = new TH1F((std::string()+"histoVbb_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoVbb_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoST_PDF_MC[c][i][j] = new TH1F((std::string()+"histoST_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoST_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoQCD_PDF_MC[c][i][j] = new TH1F((std::string()+"histoQCD_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoQCD_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoVV_PDF_MC[c][i][j] = new TH1F((std::string()+"histoVV_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+"histoVV_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoTot_PDF_MC[c][i][j] = new TH1F((std::string()+"histoTot_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+channelSuffix[c]+"histoTot_PDF_MC"+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoTotC_PDF_MC[c][i][j] = new TH1F((std::string()+"histoTotC_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+channelSuffix[c]+"histoTotC_PDF_MC"+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        histoTotS_PDF_MC[c][i][j] = new TH1F((std::string()+"histoTotS_PDF_MC"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), (std::string()+channelSuffix[c]+"histoTotS_PDF_MC"+jetSuffix[i]+wpSuffix[j]).c_str(), 10, 0, 10);
        wpJetCat.setLabel((std::string()+jetSuffix[i]+wpSuffix[j]).c_str());
        for (int bj=0; bj<6; bj++) {
          if (bj==0 || bj==1) {
            sprintf(label, "N%dbjet", bj);
          } else {
            sprintf(label, "N%dbjets", bj);
          }
          nbjets.setLabel(label); 
          histoTT_PDF_MC[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoTT_PDF_MC[c][i][j]->SetBinContent(bj+1, pbjets_tt_ext[c][i][j]->getVal());
          histoV_PDF_MC[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoV_PDF_MC[c][i][j]->SetBinContent(bj+1, pbjets_v_ext[c][i][j]->getVal());
          histoVbb_PDF_MC[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoVbb_PDF_MC[c][i][j]->SetBinContent(bj+1, pbjets_vb[c][i][j]->getVal());
            //histoVbb_PDF_MC[c][i][j]->SetBinError(bj+1, pbjets_v_ext[i][j]->getPropagatedError(*fit_result[c]));
          histoST_PDF_MC[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoST_PDF_MC[c][i][j]->SetBinContent(bj+1, pbjets_st[c][i][j]->getVal());
            //histoST_PDF_MC[c][i][j]->SetBinError(bj+1, pbjets_v_ext[i][j]->getPropagatedError(*fit_result[c]));
          histoQCD_PDF_MC[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoQCD_PDF_MC[c][i][j]->SetBinContent(bj+1, pbjets_qcd[c][i][j]->getVal());
            //histoVbb_PDF_MC[c][i][j]->SetBinError(bj+1, pbjets_v_ext[i][j]->getPropagatedError(*fit_result[c]));
          histoVV_PDF_MC[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoVV_PDF_MC[c][i][j]->SetBinContent(bj+1, pbjets_vv[c][i][j]->getVal());
            //histoST_PDF_MC[c][i][j]->SetBinError(bj+1, pbjets_v_ext[i][j]->getPropagatedError(*fit_result[c]));
          histoTot_PDF_MC[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoTot_PDF_MC[c][i][j]->SetBinContent(bj+1, model[c][i][j]->getVal());
          histoTotC_PDF_MC[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoTotC_PDF_MC[c][i][j]->SetBinContent(bj+1, model_constraint[c][i][j]->getVal());
        }
        
        
        
        for (int bj=0; bj<6; bj++) {
          if (bj==0 || bj==1) {
            sprintf(label, "N%dbjet", bj);
          } else {
            sprintf(label, "N%dbjets", bj);
          }
          nbjets.setLabel(label); 
          histoTotS_PDF_MC[c][i][j]->GetXaxis()->SetBinLabel(bj+1, nbjets.getLabel());
          histoTotS_PDF_MC[c][i][j]->SetBinContent(bj+1, pbjets_tt_ext[c][i][j]->getVal()*pbjets_tt_ext[c][i][j]->expectedEvents(NULL)/histoTT_PDF_MC[c][i][j]->Integral()
                                                   +pbjets_v_ext[c][i][j]->getVal()*pbjets_v_ext[c][i][j]->expectedEvents(NULL)/histoV_PDF_MC[c][i][j]->Integral()
                                                   +Nvb[c][i]->getVal()*pbjets_vb[c][i][j]->getVal()/histoVbb_PDF_MC[c][i][j]->Integral()
                                                   +Nst[c][i]->getVal()*pbjets_st[c][i][j]->getVal()/histoST_PDF_MC[c][i][j]->Integral()
                                                   /*                                                   +Nqcd[c][i]->getVal()*pbjets_qcd[c][i][j]->getVal()/histoQCD_PDF_MC[c][i][j]->Integral() */
                                                   +Nvv[c][i]->getVal()*pbjets_vv[c][i][j]->getVal()/histoVV_PDF_MC[c][i][j]->Integral());
            //histoTotS_PDF_MC[i][j]->SetBinError(1, sqrt(pow(pbjets_tt_ext[i][j]->getPropagatedError(*fit_result),2)));
        }
        
        
        /* // Should propagate the range for the variable and not the fit result.
         nbjets.setLabel("N0bjet"); 
         histoTotS_PDF_MC[c][i][j]->SetBinError(1, model_constraint[c][i][j]->getPropagatedError(*(fit_result[c])) * histoTotS_PDF_MC[c][i][j]->Integral()/histoTot_PDF_MC[c][i][j]->Integral());
         nbjets.setLabel("N1bjet"); 
         histoTotS_PDF_MC[c][i][j]->SetBinError(2, model_constraint[c][i][j]->getPropagatedError(*(fit_result[c])) * histoTotS_PDF_MC[c][i][j]->Integral()/histoTot_PDF_MC[c][i][j]->Integral());
         nbjets.setLabel("N2bjets"); 
         histoTotS_PDF_MC[c][i][j]->SetBinError(3, model_constraint[c][i][j]->getPropagatedError(*(fit_result[c])) * histoTotS_PDF_MC[c][i][j]->Integral()/histoTot_PDF_MC[c][i][j]->Integral());
         nbjets.setLabel("N3bjets"); 
         histoTotS_PDF_MC[c][i][j]->SetBinError(4, model_constraint[c][i][j]->getPropagatedError(*(fit_result[c])) * histoTotS_PDF_MC[c][i][j]->Integral()/histoTot_PDF_MC[c][i][j]->Integral());
         */
        
        histoTT_PDF_MC[c][i][j]->Scale(pbjets_tt_ext[c][i][j]->expectedEvents(NULL)/histoTT_PDF_MC[c][i][j]->Integral());
        histoV_PDF_MC[c][i][j]->Scale(pbjets_v_ext[c][i][j]->expectedEvents(NULL)/histoV_PDF_MC[c][i][j]->Integral());
        histoTot_PDF_MC[c][i][j]->Scale(model[c][i][j]->expectedEvents(NULL)/histoTot_PDF_MC[c][i][j]->Integral());
        histoTotC_PDF_MC[c][i][j]->Scale(model_constraint[c][i][j]->expectedEvents(NULL)/histoTotC_PDF_MC[c][i][j]->Integral());
        
        /*
         Double_t integral = 0.;
         integral = histoTT_PDF_MC[i][j]->Integral();  histoTT_PDF_MC[i][j]->Scale(Ntt[i]->getVal()/integral);
         integral = histoV_PDF_MC[i][j]->Integral();   histoV_PDF_MC[i][j]->Scale(Nv[i]->getVal()/integral);
         integral = histoTot_PDF_MC[i][j]->Integral(); histoTot_PDF_MC[i][j]->Scale((Ntt[i]->getVal()+Nv[i]->getVal()+Nst[i]->getVal()+Nvb[i]->getVal())/integral);
         */
          //      histoTT_PDF_MC[i][j]->Print();
        printf(" PDF MC TT (norm : %lf ; expectedEvents : %lf) : ", pbjets_tt_ext[c][i][j]->getNorm(), pbjets_tt_ext[c][i][j]->expectedEvents(NULL));
        for (int k=1; k<=6 ; k++) {
          printf(" %lf ", histoTT_PDF_MC[c][i][j]->GetBinContent(k));
        }
        printf("\nNtt[i] : %lf\n",Ntt[c][i]->getVal());
          //      histoV_PDF_MC[i][j]->Print();
        printf(" PDF MC V (norm : %lf ; expectedEvents : %lf) : ",  pbjets_v_ext[c][i][j]->getNorm(), pbjets_v_ext[c][i][j]->expectedEvents(NULL));
        for (int k=1; k<=6 ; k++) {
          printf(" %lf ", histoV_PDF_MC[c][i][j]->GetBinContent(k));
        }
        printf("\nNv[i] : %lf\n",Nv[c][i]->getVal());
          //      histoTot_PDF_MC[i][j]->Print();
        printf(" PDF MC Tot (norm : %lf ; expectedEvents : %lf) : ", model[c][i][j]->getNorm(), model[c][i][j]->expectedEvents(NULL));
        for (int k=1; k<=6 ; k++) {
          printf(" %lf ", histoTot_PDF_MC[c][i][j]->GetBinContent(k));
        }
        printf("\n");
        printf(" PDF MC Tot constraint (norm : %lf ; expectedEvents : %lf) : ", model_constraint[c][i][j]->getNorm(), model_constraint[c][i][j]->expectedEvents(NULL));
        for (int k=1; k<=6 ; k++) {
          printf(" %lf ", histoTotC_PDF_MC[c][i][j]->GetBinContent(k));
        }
        printf("\n");
        printf(" PDF MC Tot summed (norm : %lf ; expectedEvents : %lf) : ", -1., -1.);
        for (int k=1; k<=6 ; k++) {
          printf(" %lf ", histoTotS_PDF_MC[c][i][j]->GetBinContent(k));
        }
        printf("\n");
        
      }
    }
  }
  

  cout << endl << endl << endl << endl << endl << endl << "********************************************" << endl;
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to perform the fit" << endl;
  cout << "********************************************" << endl;
  cout << " Starting the pseudo-experiments now          " << endl;
  cout << "********************************************" << endl;
  
  /*
   //Reset to fit values for the PE
   for (int c=0; c<NbOfChannels; c++) {
   if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
   if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
   for (int i=0; i<NbOfJetBins; i++) {
   Ntt[c][i]->setVal(save_tt[c][i]);
   Nv[c][i]->setVal(save_v [c][i]);
   Nst[c][i]->setVal(save_st[c][i]);
   Nvb[c][i]->setVal(save_vb[c][i]);
   Nqcd[c][i]->setVal(save_qcd[c][i]);
   Nvv[c][i]->setVal(save_vv[c][i]);
   }
   for (int wp=0; wp<NbOfWP; wp++) {
   eb[c][wp]->setVal(save_eb[c][wp]);
   eudsc[c][wp]->setVal(save_et[c][wp]);
   euds[c][wp]->setVal(save_ev[c][wp]);
   }
   }
   */
  
  
    //    RooMCStudy* mcstudy = new RooMCStudy(simPdf,RooArgSet( nbjets ,jetMultCat, wpCat),Constrain(constraintsList),Binned(kTRUE),Silence(),Extended(kTRUE),FitOptions(Optimize(0),Save(0),Minos(1),Extended(1),Strategy(2),PrintEvalErrors(-1)));
  RooMCStudy** mcstudy = new RooMCStudy*[NbOfChannels];
  for (UInt_t c=0; c<NbOfChannels; c++) {
    mcstudy[c] = NULL;
    
      //    if (!((((c+1)&muonMask)!=0) && (((c+1)&electronMask)!=0))) continue;
    
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    printf("\n\n\nTTTT Input values for RooMCStudy for channel #%d ,\t\t  ready to generate %lf pseudo-experiments of %lf events each\n", c, (double) NbOfPE, (double) nTot[c]);
    for (UInt_t j=0; j<NbOfWP; j++) {
      printf("  %s     : %lf\n  %s  : %lf\n  %s : %lf\n", 
             eb[c][j]->getTitle().Data(), eb[c][j]->getVal(),
             euds[c][j]->getTitle().Data(), euds[c][j]->getVal(),
             eudsc[c][j]->getTitle().Data(), eudsc[c][j]->getVal());
    }
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      printf("  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s : %lf\n",
             Ntt[c][i]->getTitle().Data(), Ntt[c][i]->getVal(),
             Nv[c][i]->getTitle().Data(), Nv[c][i]->getVal(),
             Nst[c][i]->getTitle().Data(), Nst[c][i]->getVal(),
             Nvb[c][i]->getTitle().Data(), Nvb[c][i]->getVal(),
             Nvv[c][i]->getTitle().Data(), Nvv[c][i]->getVal(),
             Nqcd[c][i]->getTitle().Data(), Nqcd[c][i]->getVal());
    }
    printf("\n\n\n");
    
    
    mcstudy[c] = new RooMCStudy(*(simPdf[c]),RooArgSet( nbjets, wpJetCat),Constrain(* (constraintsList[c])),Binned(kTRUE),Silence(),Extended(kTRUE),FitOptions(Optimize(0),Save(0),Minos(1),Extended(1),Strategy(2),PrintEvalErrors(-1)),FitOptions(NumCPU(1)));
    
    
    Bool_t keepGenData=kFALSE;
    
      // Generate and fit NbOfPE samples of Poisson(nExpected) events
    if(doPE){
        //    mcstudy->generateAndFit(NbOfPE,nExpected,keepGenData) ;
      mcstudy[c]->generateAndFit(NbOfPE,nTot[c],keepGenData) ;
      
      printf("\n\n\nTTTT Output values from RooMCStudy for channel #%d\n", c);
      for (UInt_t j=0; j<NbOfWP; j++) {
        printf("  %s     : %lf\n  %s  : %lf\n  %s : %lf\n", 
               eb[c][j]->getTitle().Data(), eb[c][j]->getVal(),
               euds[c][j]->getTitle().Data(), euds[c][j]->getVal(),
               eudsc[c][j]->getTitle().Data(), eudsc[c][j]->getVal());
      }
      for (UInt_t i=0; i<NbOfJetBins; i++) {
        printf("  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s  : %lf\n  %s : %lf\n",
               Ntt[c][i]->getTitle().Data(), Ntt[c][i]->getVal(),
               Nv[c][i]->getTitle().Data(), Nv[c][i]->getVal(),
               Nst[c][i]->getTitle().Data(), Nst[c][i]->getVal(),
               Nvb[c][i]->getTitle().Data(), Nvb[c][i]->getVal(),
               Nvv[c][i]->getTitle().Data(), Nvv[c][i]->getVal(),
               Nqcd[c][i]->getTitle().Data(), Nqcd[c][i]->getVal());
      }
      printf("\n\n\n");
      
    }
  }
  
  
    // C o n s t r u c t   p r o f i l e   l i k e l i h o o d   i n   f r a c
    // -----------------------------------------------------------------------
  
    // The profile likelihood estimator on nll for Ntt will minimize nll w.r.t
    // all floating parameters except frac for each evaluation
  /*
   RooAbsReal* pll_Ntt = nll->createProfile(Ntt) ;
   
   // Plot the profile likelihood in frac
   RooPlot* frame1 = Ntt.frame(Title("ProfileLL in Ntt")) ;
   pll_Ntt->plotOn(frame1,LineColor(kRed)) ;
   */
    // P l o t   f i t   r e s u l t s 
    // ---------------------------------------------------------------------
  
    // Saving efficiencies
  /*
   ebTTMC_histo->Write();
   ebVMC_histo->Write();
   eudsMC_histo->Write();
   eudscMC_histo->Write();
   */
  
  for (UInt_t c=0; c<NbOfChannels; c++) {
    
      //    if (!((((c+1)&muonMask)!=0) && (((c+1)&electronMask)!=0))) continue;
    
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    for(UInt_t i=0; i<NbOfJetBins; i++) {
      for(UInt_t j=0; j<NbOfWP ; j++) {
        double min = FLT_MAX;
        double max = FLT_MIN;
        double mmin = 0.;
        double mmax = 0.;
        TCanvas *cTT = new TCanvas((std::string()+"distributions_ttlike"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),
                                   (std::string()+"distributions_ttlike"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
        cTT->cd();
        TLegend *ttLeg = new TLegend(0.8,0.8,0.99,0.99);
        mmax = hNbtaggedJets_ttlike[c][i][j]->GetBinContent(hNbtaggedJets_ttlike[c][i][j]->GetMaximumBin());
        mmin = hNbtaggedJets_ttlike[c][i][j]->GetBinContent(hNbtaggedJets_ttlike[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        hNbtaggedJets_ttlike[c][i][j]->SetMarkerStyle(20);
        hNbtaggedJets_ttlike[c][i][j]->SetMarkerSize(0.8);//1.5
        hNbtaggedJets_ttlike[c][i][j]->SetMarkerColor(kBlack);
        hNbtaggedJets_ttlike[c][i][j]->SetLineColor(kBlack);
        ttLeg->AddEntry(hNbtaggedJets_ttlike[c][i][j], "#b-tags","LP");
        mmax = histoTT_PDF_Fit[c][i][j]->GetBinContent(histoTT_PDF_Fit[c][i][j]->GetMaximumBin());
        mmin = histoTT_PDF_Fit[c][i][j]->GetBinContent(histoTT_PDF_Fit[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        histoTT_PDF_Fit     [c][i][j]->SetMarkerStyle(21);
        histoTT_PDF_Fit     [c][i][j]->SetMarkerSize(0.8);//1.
        histoTT_PDF_Fit     [c][i][j]->SetMarkerColor(kRed);
        histoTT_PDF_Fit     [c][i][j]->SetLineColor(kRed);
        ttLeg->AddEntry(histoTT_PDF_Fit     [c][i][j], "PDF after fit", "LP");
        mmax = histoTT_PDF_MC[c][i][j]->GetBinContent(histoTT_PDF_MC[c][i][j]->GetMaximumBin());
        mmin = histoTT_PDF_MC[c][i][j]->GetBinContent(histoTT_PDF_MC[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        histoTT_PDF_MC      [c][i][j]->SetMarkerStyle(22);
        histoTT_PDF_MC      [c][i][j]->SetMarkerSize(0.8);//1
        histoTT_PDF_MC      [c][i][j]->SetMarkerColor(kGreen);
        histoTT_PDF_MC      [c][i][j]->SetLineColor(kGreen);
        ttLeg->AddEntry(histoTT_PDF_MC      [c][i][j], "PDF on MC","LP");
        cTT->DrawFrame(0., 0.9*min, 10., 1.1*max);
        hNbtaggedJets_ttlike[c][i][j]->Draw("e0same");
        histoTT_PDF_Fit     [c][i][j]->Draw("e0same");
        histoTT_PDF_MC      [c][i][j]->Draw("e0same");
        ttLeg->Draw("same");
        cTT->Print((std::string()+postfix+"_"+cTT->GetName()+".pdf").c_str());
        cTT->Write();
        
        
        min = FLT_MAX;
        max = FLT_MIN;
        TCanvas *cV = new TCanvas((std::string()+"distributions_vlike"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),
                                  (std::string()+"distributions_vlike"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
        cV->cd();
        TLegend *vLeg = new TLegend(0.8,0.8,0.99,0.99);
        mmax = hNbtaggedJets_vlike[c][i][j]->GetBinContent(hNbtaggedJets_vlike[c][i][j]->GetMaximumBin());
        mmin = hNbtaggedJets_vlike[c][i][j]->GetBinContent(hNbtaggedJets_vlike[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        hNbtaggedJets_vlike [c][i][j]->SetMarkerStyle(20);
        hNbtaggedJets_vlike [c][i][j]->SetMarkerSize(0.8);//1.5
        hNbtaggedJets_vlike [c][i][j]->SetMarkerColor(kBlack);
        hNbtaggedJets_vlike [c][i][j]->SetLineColor(kBlack);
        vLeg->AddEntry(hNbtaggedJets_vlike[c][i][j], "#b-tags","LP");
        mmax = histoTT_PDF_Fit[c][i][j]->GetBinContent(histoTT_PDF_Fit[c][i][j]->GetMaximumBin());
        mmin = histoTT_PDF_Fit[c][i][j]->GetBinContent(histoTT_PDF_Fit[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        histoV_PDF_Fit      [c][i][j]->SetMarkerStyle(21);
        histoV_PDF_Fit      [c][i][j]->SetMarkerSize(0.8);//1
        histoV_PDF_Fit      [c][i][j]->SetMarkerColor(kRed);
        histoV_PDF_Fit      [c][i][j]->SetLineColor(kRed);
        vLeg->AddEntry(histoV_PDF_Fit     [c][i][j], "PDF after fit", "LP");
        mmax = histoV_PDF_MC[c][i][j]->GetBinContent(histoV_PDF_MC[c][i][j]->GetMaximumBin());
        mmin = histoV_PDF_MC[c][i][j]->GetBinContent(histoV_PDF_MC[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        histoV_PDF_MC       [c][i][j]->SetMarkerStyle(22);
        histoV_PDF_MC       [c][i][j]->SetMarkerSize(1);
        histoV_PDF_MC       [c][i][j]->SetMarkerColor(kGreen);
        histoV_PDF_MC       [c][i][j]->SetLineColor(kGreen);
        vLeg->AddEntry(histoV_PDF_MC      [c][i][j], "PDF on MC","LP");
        mmax = histoV_PDF_MC[c][i][j]->GetBinContent(histoV_PDF_MC[c][i][j]->GetMaximumBin());
        mmin = histoV_PDF_MC[c][i][j]->GetBinContent(histoV_PDF_MC[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        cV->DrawFrame(0., 0.9*min, 10., 1.1*max);
        hNbtaggedJets_vlike [c][i][j]->Draw("e0same");
        histoV_PDF_Fit      [c][i][j]->Draw("e0same");
        histoV_PDF_MC       [c][i][j]->Draw("e0same");
        vLeg->Draw("same");
        cV->Print((std::string()+postfix+"_"+cV->GetName()+".pdf").c_str());
        cV->Write();
        
        min = FLT_MAX;
        max = FLT_MIN;
        TCanvas *cTot = new TCanvas((std::string()+"distributions_totlike"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),
                                    (std::string()+"distributions_totlike"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
        cTot->cd();
        TLegend *totLeg = new TLegend(0.8,0.8,0.99,0.99);
        hNbtaggedJets    [c][i][j]->SetMarkerStyle(20);
        hNbtaggedJets    [c][i][j]->SetMarkerSize(1.5);
        hNbtaggedJets    [c][i][j]->SetMarkerColor(kBlack);
        hNbtaggedJets    [c][i][j]->SetLineColor(kBlack);
        totLeg->AddEntry(hNbtaggedJets[c][i][j], "#b-tags","LP");
        mmax = hNbtaggedJets[c][i][j]->GetBinContent(hNbtaggedJets[c][i][j]->GetMaximumBin());
        mmin = hNbtaggedJets[c][i][j]->GetBinContent(hNbtaggedJets[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        /*
         histoTot_PDF_Fit [c][i][j]->SetMarkerStyle(22);
         histoTot_PDF_Fit [c][i][j]->SetMarkerSize(1);
         histoTot_PDF_Fit [c][i][j]->SetMarkerColor(kRed);
         histoTot_PDF_Fit [c][i][j]->SetLineColor(kRed);
         totLeg->AddEntry(histoTot_PDF_Fit     [c][i][j], "PDF after fit (RooFit)", "LP");
         mmax = histoTot_PDF_Fit[c][i][j]->GetBinContent(histoTot_PDF_MC[c][i][j]->GetMaximumBin());
         mmin = histoTot_PDF_Fit[c][i][j]->GetBinContent(histoTot_PDF_MC[c][i][j]->GetMinimumBin());
         if (mmin<min) min = mmin;
         if (mmax>max) max = mmax;
         histoTot_PDF_MC  [c][i][j]->SetMarkerStyle(23);
         histoTot_PDF_MC  [c][i][j]->SetMarkerSize(1);
         histoTot_PDF_MC  [c][i][j]->SetMarkerColor(kGreen);
         histoTot_PDF_MC  [c][i][j]->SetLineColor(kGreen);
         totLeg->AddEntry(histoTot_PDF_MC      [c][i][j], "PDF on MC (RooFit)","LP");
         mmax = histoTot_PDF_MC[c][i][j]->GetBinContent(histoTot_PDF_MC[c][i][j]->GetMaximumBin());
         mmin = histoTot_PDF_MC[c][i][j]->GetBinContent(histoTot_PDF_MC[c][i][j]->GetMinimumBin());
         if (mmin<min) min = mmin;
         if (mmax>max) max = mmax;
         */
        /*
         histoTotC_PDF_Fit [i][j]->SetMarkerStyle(2);
         histoTotC_PDF_Fit [i][j]->SetMarkerColor(2);
         histoTotC_PDF_Fit [i][j]->SetLineColor(2);
         totLeg->AddEntry(histoTotC_PDF_Fit     [i][j], "PDF after fit (constraint)", "LP");
         mmax = histoTotC_PDF_Fit[i][j]->GetBinContent(histoTotC_PDF_Fit[i][j]->GetMaximumBin());
         mmin = histoTotC_PDF_Fit[i][j]->GetBinContent(histoTotC_PDF_Fit[i][j]->GetMinimumBin());
         if (mmin<min) min = mmin;
         if (mmax>max) max = mmax;
         */
        histoTotS_PDF_Fit [c][i][j]->SetMarkerStyle(21);
        histoTotS_PDF_Fit [c][i][j]->SetMarkerSize(1);
        histoTotS_PDF_Fit [c][i][j]->SetMarkerColor(kRed);
        histoTotS_PDF_Fit [c][i][j]->SetLineColor(kRed);
        totLeg->AddEntry(histoTotS_PDF_Fit     [c][i][j], "PDF after fit", "LP");
        mmax = histoTotS_PDF_Fit[c][i][j]->GetBinContent(histoTotS_PDF_Fit[c][i][j]->GetMaximumBin());
        mmin = histoTotS_PDF_Fit[c][i][j]->GetBinContent(histoTotS_PDF_Fit[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        
        histoTotS_PDF_MC [c][i][j]->SetMarkerStyle(22);
        histoTotS_PDF_MC [c][i][j]->SetMarkerSize(1);
        histoTotS_PDF_MC [c][i][j]->SetMarkerColor(kGreen);
        histoTotS_PDF_MC [c][i][j]->SetLineColor(kGreen);
        totLeg->AddEntry(histoTotS_PDF_MC     [c][i][j], "PDF on MC", "LP");
        mmax = histoTotS_PDF_MC[c][i][j]->GetBinContent(histoTotS_PDF_MC[c][i][j]->GetMaximumBin());
        mmin = histoTotS_PDF_MC[c][i][j]->GetBinContent(histoTotS_PDF_MC[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        
        cTot->DrawFrame(0., 0.9*min, 10., 1.1*max);
        hNbtaggedJets    [c][i][j]->Draw("e0same");
          //        histoTot_PDF_MC  [c][i][j]->Draw("same");
          //        histoTot_PDF_Fit [c][i][j]->Draw("same");
        histoTotS_PDF_MC  [c][i][j]->Draw("e0same");
        histoTotS_PDF_Fit  [c][i][j]->Draw("e0same");
        totLeg->Draw("same");
        cTot->Print((std::string()+postfix+"_"+cTot->GetName()+".pdf").c_str());
        cTot->Write();
      }
    }
    
    
    
    for (UInt_t i=0; i<NbOfJetBins; i++) {
      for (UInt_t wp=0; wp<NbOfWP; wp++) {
        TCanvas *c1 = new TCanvas((std::string()+"PDFsplitted_for_Fit"+jetSuffix[i]+wpSuffix[wp]).c_str(), (std::string()+"PDF_for_Fit"+jetSuffix[i]+wpSuffix[wp]).c_str()) ;
        c1->cd();
          //frame1->Write();
          //hNbtaggedJets_ttjets->Draw("COLZ");
        hNbtaggedJets_ttjets[c][i][wp]->Scale(1./IntLumi[c]);
        hNbtaggedJets_ttjets[c][i][wp]->Write();
        hNbtaggedJets_ttjets[c][i][wp]->DrawNormalized()->SetLineColor(kRed);
        
          //hNbtaggedJets_wjets->Draw("COLZ");
        hNbtaggedJets_wjets[c][i][wp]->Scale(1./IntLumi[c]);
        hNbtaggedJets_wjets[c][i][wp]->Write();
          //        hWLFg_In_WJets[c][i][wp]->DrawNormalized("same")->SetLineColor(kAzure);
          //        hWcx_In_WJets[c][i][wp]->DrawNormalized("same")->SetLineColor(kViolet);
          //hNbtaggedJets_zjets->Draw("COLZ");
        hNbtaggedJets_zjets[c][i][wp]->Scale(1./IntLumi[c]);
        hNbtaggedJets_zjets[c][i][wp]->Write();
        hNbtaggedJets_zjets[c][i][wp]->DrawNormalized("same")->SetLineColor(kTeal);
          //hNbtaggedJets_stjets->Draw("COLZ");
        hNbtaggedJets_stjets[c][i][wp]->Scale(1./IntLumi[c]);
        hNbtaggedJets_stjets[c][i][wp]->Write();
        hNbtaggedJets_stjets[c][i][wp]->DrawNormalized("same")->SetLineColor(kBlack);
          //hNbtaggedJets_vbjets->Draw("COLZ");
        hNbtaggedJets_vbjets[c][i][wp]->Scale(1./IntLumi[c]);
        hNbtaggedJets_vbjets[c][i][wp]->Write();
        hNbtaggedJets_vbjets[c][i][wp]->DrawNormalized("same")->SetLineColor(kCyan);
          //hNbtaggedJets_stjets->Draw("COLZ");
        hNbtaggedJets_vvjets[c][i][wp]->Scale(1./IntLumi[c]);
        hNbtaggedJets_vvjets[c][i][wp]->Write();
        hNbtaggedJets_vvjets[c][i][wp]->DrawNormalized("same")->SetLineColor(kYellow);
          //hNbtaggedJets_vbjets->Draw("COLZ");
        hNbtaggedJets_qcd[c][i][wp]->Scale(1./IntLumi[c]);
        hNbtaggedJets_qcd[c][i][wp]->Write();
        if (hNbtaggedJets_qcd[c][i][wp]->Integral() != 0.) {
          hNbtaggedJets_qcd[c][i][wp]->DrawNormalized("same")->SetLineColor(kPink);
        } else {
          hNbtaggedJets_qcd[c][i][wp]->DrawCopy("same")->SetLineColor(kPink);
        }
        if (addRareProcesses == kTRUE) {
          hNbtaggedJets_rare[c][i][wp]->Scale(1./IntLumi[c]);
          hNbtaggedJets_rare[c][i][wp]->Write();
          hNbtaggedJets_rare[c][i][wp]->DrawNormalized("same")->SetLineColor(kGreen);
        }
        if (signalFile.compare("")!=0) {
          hNbtaggedJets_signal[c][i][wp]->Scale(1./IntLumi[c]);
          hNbtaggedJets_signal[c][i][wp]->Write();
          hNbtaggedJets_signal[c][i][wp]->DrawNormalized("same")->SetLineColor(kMagenta);
        }
         
        c1->Update();
        c1->Write();
        
        if (!isOnData) {
          hNbtaggedJets[c][i][wp]->Scale(1./IntLumi[c]);
        }
        hNbtaggedJets[c][i][wp]->Write();
        hNbtaggedJets_ttlike[c][i][wp]->Scale(1./IntLumi[c]);
        hNbtaggedJets_ttlike[c][i][wp]->Write();
        hNbtaggedJets_vlike[c][i][wp]->Scale(1./IntLumi[c]);
        hNbtaggedJets_vlike[c][i][wp]->Write();
        TCanvas *c2 = new TCanvas((std::string()+"PDF_for_Fit"+jetSuffix[i]+wpSuffix[wp]).c_str(), (std::string()+"PDF_for_Fit"+jetSuffix[i]+wpSuffix[wp]).c_str()) ;
        c2->cd();
        hNbtaggedJets_ttlike[c][i][wp]->DrawNormalized()->SetLineColor(kRed);
        hNbtaggedJets_vlike[c][i][wp]->DrawNormalized("same")->SetLineColor(kBlue);
        hNbtaggedJets_stjets[c][i][wp]->DrawNormalized("same")->SetLineColor(kBlack);
        hNbtaggedJets_vbjets[c][i][wp]->DrawNormalized("same")->SetLineColor(kCyan);
        hNbtaggedJets_vvjets[c][i][wp]->DrawNormalized("same")->SetLineColor(kYellow);
        if (hNbtaggedJets_qcd[c][i][wp]->Integral() != 0.) {
          hNbtaggedJets_qcd[c][i][wp]->DrawNormalized("same")->SetLineColor(kPink);
        } else {
          hNbtaggedJets_qcd[c][i][wp]->DrawCopy("same")->SetLineColor(kPink);
        }
        if (addRareProcesses == kTRUE) {
          hNbtaggedJets_rare[c][i][wp]->DrawNormalized("same")->SetLineColor(kMagenta);
        }
        if (signalFile.compare("")!=0) {
          hNbtaggedJets_signal[c][i][wp]->DrawNormalized("same")->SetLineColor(kMagenta);
        }
        c2->Update();
        c2->Write();
        
        /*
         THStack hStackWJets((std::string()+"hStackWJets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[wp]).c_str(),"hStackWJets;b-jet multiplicity;");
         hWLFg_In_WJets[c][i][wp]->SetFillColor(5);
         hStackWJets.Add(hWLFg_In_WJets[c][i][wp]);
         hWcx_In_WJets[c][i][wp]->SetFillColor(4);
         hStackWJets.Add(hWcx_In_WJets[c][i][wp]);
         hWbb_In_WJets[c][i][wp]->SetFillColor(2);
         hStackWJets.Add(hWbb_In_WJets[c][i][wp]);
         hStackWJets.Write();
         
         //Save histograms save to normalisation of 1pb-1
         hWbb_In_WJets[c][i][wp]->Scale(1./IntLumi[c]);
         hWbb_In_WJets[c][i][wp]->Write();
         hWcx_In_WJets[c][i][wp]->Scale(1./IntLumi[c]);
         hWcx_In_WJets[c][i][wp]->Write();
         hWLFg_In_WJets[c][i][wp]->Scale(1./IntLumi[c]);
         hWLFg_In_WJets[c][i][wp]->Write();
         */
        
        
        histoTT_PDF_MC[c][i][wp]->Write();
        histoV_PDF_MC[c][i][wp]->Write();
        histoVbb_PDF_MC[c][i][wp]->Write();
        histoST_PDF_MC[c][i][wp]->Write();
        histoTot_PDF_MC[c][i][wp]->Write();
        histoTotC_PDF_MC[c][i][wp]->Write();
        histoTotS_PDF_MC[c][i][wp]->Write();
        
        histoTT_PDF_Fit[c][i][wp]->Write();
        histoV_PDF_Fit[c][i][wp]->Write();
        histoVbb_PDF_Fit[c][i][wp]->Write();
        histoST_PDF_Fit[c][i][wp]->Write();
        histoTot_PDF_Fit[c][i][wp]->Write();
        histoTotC_PDF_Fit[c][i][wp]->Write();
        histoTotS_PDF_Fit[c][i][wp]->Write();
        
        
      }
      eXbq__[c][i]->Write();
      
      /*  
       RooDataSet* dd = mcstudy->genParDataSet();
       dd->fillHistogram(hist, RooArgList());
       */
      TH1* Ntt_Nv_Gen = mcstudy[c]->fitParDataSet().createHistogram(*(Ntt[c][i]),*(Nv[c][i])) ;
      Ntt_Nv_Gen->Write();
    }  
    /*
     for (int sampleNum =0 ; sampleNum<nExpected ; sampleNum++) {
     RooAbsData* ab = mcstudy->genData(sampleNum);
     }
     */
    
    if(doPE){
      for (UInt_t i=0; i<NbOfJetBins; i++) {
        RooPlot* frame_Nttpull = mcstudy[c]->plotPull(* Ntt[c][i],-4,4,100) ;
        frame_Nttpull->SetName((std::string()+"myRooPlot_Nttpull"+channelSuffix[c]+jetSuffix[i]).c_str());
        frame_Nttpull->Write();
        RooPlot* frame_Ntt = mcstudy[c]->plotParam(* Ntt[c][i]);//,Binning(100,NbOfTTlike[NbOfJets-3]*(1-0.2),NbOfTTlike[NbOfJets-3]*(1+0.2))) ;
        frame_Ntt->SetName((std::string()+"myRooPlot_Ntt"+channelSuffix[c]+jetSuffix[i]).c_str());
        frame_Ntt->Write();
        RooPlot* frame_NttErr = mcstudy[c]->plotError(* Ntt[c][i]);
        frame_NttErr->SetName((std::string()+"myRooPlot_NttErr"+channelSuffix[c]+jetSuffix[i]).c_str());
        frame_NttErr->Write();
        
        RooPlot* frame_Nvpull  = mcstudy[c]->plotPull(* Nv[c][i],-4,4,100) ;
        frame_Nvpull->SetName((std::string()+"myRooPlot_Nvpull"+channelSuffix[c]+jetSuffix[i]).c_str());
        frame_Nvpull->Write();
        RooPlot* frame_Nv = mcstudy[c]->plotParam(* Nv[c][i]);//,Binning(100,NbOfVlike[NbOfJets-3]*(1-0.2),NbOfVlike[NbOfJets-3]*(1+0.2))) ;
        frame_Nv->SetName((std::string()+"myRooPlot_Nv"+channelSuffix[c]+jetSuffix[i]).c_str());
        frame_Nv->Write();
        RooPlot* frame_NvErr = mcstudy[c]->plotError(* Nv[c][i]);
        frame_NvErr->SetName((std::string()+"myRooPlot_NvErr"+channelSuffix[c]+jetSuffix[i]).c_str());
        frame_NvErr->Write();
        
        try {
          RooPlot* frame_Nstpull  = mcstudy[c]->plotPull(* Nst[c][i],-4,4,100) ;
          frame_Nstpull->SetName((std::string()+"myRooPlot_Nstpull"+channelSuffix[c]+jetSuffix[i]).c_str());
          frame_Nstpull->Write();
          RooPlot* frame_Nst = mcstudy[c]->plotParam(* Nst[c][i],FrameBins(100)) ;
          frame_Nst->SetName((std::string()+"myRooPlot_Nst"+channelSuffix[c]+jetSuffix[i]).c_str());
          frame_Nst->Write();
        } catch (...) { ; }
        try {
          RooPlot* frame_Nvbpull  = mcstudy[c]->plotPull(* Nvb[c][i],-4,4,100) ;
          frame_Nvbpull->SetName((std::string()+"myRooPlot_Nvbpull"+channelSuffix[c]+jetSuffix[i]).c_str());
          frame_Nvbpull->Write();
          RooPlot* frame_Nvb = mcstudy[c]->plotParam(* Nvb[c][i],FrameBins(100)) ;
          frame_Nvb->SetName((std::string()+"myRooPlot_Nvb"+channelSuffix[c]+jetSuffix[i]).c_str());
          frame_Nvb->Write();
        } catch (...) { ; }
        try {
          RooPlot* frame_Nvvpull  = mcstudy[c]->plotPull(* Nvv[c][i],-4,4,100) ;
          frame_Nvvpull->SetName((std::string()+"myRooPlot_Nvvpull"+channelSuffix[c]+jetSuffix[i]).c_str());
          frame_Nvvpull->Write();
          RooPlot* frame_Nvv = mcstudy[c]->plotParam(* Nvv[c][i],FrameBins(100)) ;
          frame_Nvv->SetName((std::string()+"myRooPlot_Nvv"+channelSuffix[c]+jetSuffix[i]).c_str());
          frame_Nvv->Write();
        } catch (...) { ; }
        /*
         try {
         RooPlot* frame_Nqcdpull  = mcstudy[c]->plotPull(* Nqcd[c][i],-4,4,100) ;
         frame_Nqcdpull->SetName((std::string()+"myRooPlot_Nqcdpull"+channelSuffix[c]+jetSuffix[i]).c_str());
         frame_Nqcdpull->Write();
         RooPlot* frame_Nqcd = mcstudy[c]->plotParam(* Nqcd[c][i],FrameBins(100)) ;
         frame_Nqcd->SetName((std::string()+"myRooPlot_Nqcd"+channelSuffix[c]+jetSuffix[i]).c_str());
         frame_Nqcd->Write();
         } catch (...) { ; }
         */
        hMC_Nominal_N[c][i]->Write();
      }
      
      RooPlot* frame_NLL = mcstudy[c]->plotNLL();
      frame_NLL->SetName("myRooPlot_NLL");
      frame_NLL->Write();
      
      for (UInt_t j=0; j<NbOfWP; j++) {
        if(!eb[c][j]->isConstant()){
          RooPlot* frame_ebpull = mcstudy[c]->plotPull(*(eb[c][j]),-4,4,100) ;
          frame_ebpull->SetName((std::string()+"myRooPlot_Ebpull"+channelSuffix[c]+wpSuffix[j]).c_str());
          frame_ebpull->Write();
          RooPlot* frame_eb = mcstudy[c]->plotParam(*(eb[c][j]),FrameBins(100)) ;
          frame_eb->SetName((std::string()+"myRooPlot_eb"+channelSuffix[c]+wpSuffix[j]).c_str());
          frame_eb->Write();
        }
        if(!eudsc[c][j]->isConstant()){
          RooPlot* frame_eudscpull = mcstudy[c]->plotPull(*(eudsc[c][j]),-4,4,100) ;
          frame_eudscpull->SetName((std::string()+"myRooPlot_Eudscpull"+channelSuffix[c]+wpSuffix[j]).c_str());
          frame_eudscpull->Write();
          RooPlot* frame_eudsc = mcstudy[c]->plotParam(*(eudsc[c][j]),FrameBins(100)) ;
          frame_eudsc->SetName((std::string()+"myRooPlot_eudsc"+channelSuffix[c]+wpSuffix[j]).c_str());
          frame_eudsc->Write();
        }
        if(!euds[c][j]->isConstant()){
          RooPlot* frame_eudspull  = mcstudy[c]->plotPull(*(euds[c][j]),-4,4,100) ;
          frame_eudspull->SetName((std::string()+"myRooPlot_Eudspull"+channelSuffix[c]+wpSuffix[j]).c_str());
          frame_eudspull->Write();
          RooPlot* frame_euds = mcstudy[c]->plotParam(*(euds[c][j]),FrameBins(100)) ;
          frame_euds->SetName((std::string()+"myRooPlot_euds"+channelSuffix[c]+wpSuffix[j]).c_str());
          frame_euds->Write();
        }
        hMC_Nominal_E[c][j]->Write();
      }
    }
    
    histoNJ_ttjets[c]->Write();
    histoNJ_wjets[c]->Write();
    histoNJ_zjets[c]->Write();
    histoNJ_stjets[c]->Write();
  }
  
  
  
  /*
   TCanvas *c = new TCanvas("RooFit_FitPlot", "Fit Plot");
   c->cd();
   RooPlot *fitPlot=new RooPlot(0.,4.);
   fitPlot->updateNormVars(RooArgSet(nbjets));
   data.plotOn(fitPlot);
   model.plotOn(fitPlot);
   fitPlot->Draw();
   
   c->Write();
   */
  
    // Closing files
  fout->Close();
  
    // Free memory (if needed)
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;
  
  return 0;
}

/*
 gcc -O $(root-config --cflags) $(root-config --libs)  -lMathMore -lFoam -lMinuit -lRooFitCore -lRooFit -L/sandbox/cmss/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_2_8/external/slc5_amd64_gcc434/lib -ldcap -o StopSearches_VJetsBckgdEst_NoFiles_SevJets.exe StopSearches_VJetsBckgdEst_NoFiles_SevJets.cc
 g++ -O -g $(root-config --cflags) $(root-config --libs)  -lMathMore -lFoam -lMinuit -lRooFitCore -lRooFit -o StopSearches_VJetsBckgdEst.exe UserCode/caebergs/Stop_Systematics/VJetEstimation.cc StopSearches_VJetsBckgdEst_8TeV.cc
 */

int main (int argc, char *argv[])
{
  if (argc>=3) {
    Double_t signal_xs = 0.;
    signal_xs = atof(argv[2]);
    return StopSearches_VJetsBckgdEst_NoFiles_SevJets(std::string(argv[1]), signal_xs);
  } else {
    return StopSearches_VJetsBckgdEst_NoFiles_SevJets();
  }
}

