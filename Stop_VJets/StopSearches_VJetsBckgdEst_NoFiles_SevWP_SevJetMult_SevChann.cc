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

  //TITI
#include "/user/caebergs/VJetEstimation/StopBckg/MistagFuncs.C"
#include "/user/caebergs/VJetEstimation/LumiReWeighting.h"
#include "/user/caebergs/VJetEstimation/StopBckg/SFlightFuncs.C"

#include "/user/caebergs/VJetEstimation/Muon/AnaTree.h"
#include "/user/caebergs/VJetEstimation/Electrons/AnaTree.h"

using namespace reweight;
/**
 Results from BTV-11-003, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG link to https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/eff_b_c-ttbar_payload.txt 
 In Range : 0<=|eta|<=2.4 , 30<=p_T<=200  (approximately)
 */
/** TCHE */
Float_t eff_b_ttbar_FtCM(Float_t /*b-disc*/ x) {
  return -3.67153247396e-07*x*x*x*x +  -2.81599797034e-05*x*x*x +  0.00293190163243*x*x +  -0.0849600849778*x +  0.928524440715 ;
};
Float_t eff_b_ttbar_FtCM_pluserr(Float_t /*b-disc*/ x) {
  return 3.03337430722e-06*x*x*x*x + -0.000171604835897*x*x*x + 0.00474711667943*x*x + -0.0929933040514*x + 0.978347619293 ;
};
Float_t eff_b_ttbar_MC(Float_t /*b-disc*/ x) {
  return 3.90732786802e-06*x*x*x*x +  -0.000239934437355*x*x*x +  0.00664986827287*x*x +  -0.112578996016*x +  1.00775721404 ;
};
Float_t eff_c_ttbar_MC(Float_t /*b-disc*/ x) {
  return 0.343760640168*exp(-0.00315525164823*x*x*x + 0.0805427315196*x*x + -0.867625139194*x + 1.44815935164 ) ;
};



  // BTV-11-004
const Int_t mistag_eta_binning_loose_NUM = 5;
Double_t mistag_eta_binning_loose[mistag_eta_binning_loose_NUM]  = { 0., 0.5, 1., 1.5, 2.4 };
const Int_t mistag_eta_binning_medium_NUM = 4;
Double_t mistag_eta_binning_medium[mistag_eta_binning_medium_NUM] = { 0., 0.8, 1.6, 2.4 };
const Int_t mistag_eta_binning_tight_NUM = 0;
Double_t mistag_eta_binning_tight[mistag_eta_binning_tight_NUM]   = {};
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
       */
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



int StopSearches_VJetsBckgdEst_NoFiles_SevJets()
{
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
  
  
  TF1** sf_light_mean = new TF1*[mistag_eta_binning_medium_NUM-1];
  for (Int_t i=1; i<mistag_eta_binning_medium_NUM; i++) {
    sf_light_mean[i-1] = NULL;
    sf_light_mean[i-1] = GetSFlmean("TCHE","M",mistag_eta_binning_medium[i-1], mistag_eta_binning_medium[i]);
  }
  
  TF1** sf_light_mean_WPloose = new TF1*[mistag_eta_binning_loose_NUM-1];
  for (Int_t i=1; i<mistag_eta_binning_loose_NUM; i++) {
    sf_light_mean_WPloose[i-1] = NULL;
    sf_light_mean_WPloose[i-1] = GetSFlmean("TCHE","L",mistag_eta_binning_loose[i-1], mistag_eta_binning_loose[i]);
  }
  
  
  clock_t start = clock();
  
  cout << "***************************************************************" << endl;
  cout << " Beginning of the program for the V+Jets background estimation " << endl;
  cout << "***************************************************************" << endl;
  
  bool verbose = true;
  bool areDataChannelTreatedExclusive = false;
  bool runOnPNFS = true;
  bool Wbb_from_WJets = true;
  bool doPE = true;
  bool isOnData = true;
  bool reloadHisto = false;
  bool isWbbFromWbbMeas = false; //else, take the measure cross-section ; not using the Flavour History Path
  bool semiMuon = true;
  bool semiElectron = false;
  
  
  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////Configuration ////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  std::string WbbShapeRootfile = "../StopBckg_Outputtttttt.root" ; //For Flavour history path ???
  
    //Output ROOT file
  string postfix = "_VJetsBckgdEst";
  string channelpostfix = "_";
  if (semiMuon) {
    channelpostfix = channelpostfix + "SemiMuon"; 
  }
  if (semiElectron) {
    channelpostfix = channelpostfix + "SemiElectron"; 
  }
  string comment = "_Constants_euds";
  if (isOnData) {
    comment = comment + "_Data" ;
  }
  string rootFileName ("StopSearches"+postfix+channelpostfix+comment+".root");
  std::string histoFileName = rootFileName ;
  if (reloadHisto) {
    rootFileName = "_" + rootFileName;
  }
  
  //  std::string albertoFileName(".root");
  

  TFile *pileUpRewightingFile = new TFile("../PileUpWeights2011.root","read");
  TH1F * numInteractions = (TH1F*) pileUpRewightingFile->Get("data_pileup");
  
  /*
   LumiReWeighting LumiWeights;
   LumiWeights = reweight::LumiReWeighting("../pileup_MC_Fall11.root", "../pileup_2011Data_UpToRun180252.root", "pileup", "pileup");
   */
  
  /** TO BE ADAPTED **/
  
  
  
  TChain *dataTree_Muon = new TChain("AnaTree");
  TChain *dataTree_Electron = new TChain("AnaTree");
  if (!reloadHisto) {
    if (runOnPNFS && isOnData) {
      dataTree_Muon->Add("dcap:///dataTree.root");
      dataTree_Electron->Add("dcap:///dataTree.root");
      
        //    dataTree->Add("dcap:///pnfs/iihe/cms/store/user/caebergs/AnaTree/Data2011/2011B_ReReco_19Nov2011_V1_AOD_official_total/Treefiles_from_step2_lumi400/treeFile_Data_SingleMu_2011B_19Nov_rereco_lumi400_99_1_ZpG.root");
      
        //    dataTree->Add("dcap:///pnfs/iihe/cms/store/user/caebergs/AnaTree/Data2011/2011A_ReReco_08Nov2011_V1_AOD_official_part2/Treefiles_from_step2_lumi400/treeFile_Data_SingleMu_2011A_08Nov_rereco_part2_lumi400_100_1_kDg.root");
      
        //    dataTree->Add("dcap:///pnfs/iihe/cms/store/user/caebergs/AnaTree/Data2011/2011A_ReReco_08Nov2011_V1_AOD_official_part1/Treefiles_from_step2_lumi400/treeFile_Data_SingleMu_2011A_08Nov_rereco_part1_lumi400_100_1_UOV.root");
      
    } else {
      ;  //      data->Add("/user/nstrobbe/DATA/stop/skims_443/TTbar_skim.root");
    }
  }
  TEntryList *entrylist_data= NULL;//new TEntryList("entrylist_data", "entrylist_data");
  if(verbose && isOnData && !reloadHisto){
    cout<<"Number of events for :"<<endl;
    if (semiMuon) {
      cout<<" - data : Muon = "<<dataTree_Muon->GetEntries()<<endl;
    }
    if (semiElectron) {
      cout<<" - data : Electron = "<<dataTree_Electron->GetEntries()<<endl;
    }
    
  }
  
  TChain *ttjets_Muon = new TChain("AnaTree");
  TChain *ttjets_Electron = new TChain("AnaTree");
  if (!reloadHisto) {
    if (runOnPNFS) {
      ttjets->Add("dcap:///ttjetsTree.root");
        //    ttjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/big_TTbar_skim.root");
    } else {
      ttjets_Muon->Add("/user/nstrobbe/DATA/stop/skims_443/TTbar_skim.root");
      ttjets_Electron->Add("/user/nstrobbe/DATA/stop/skims_443/TTbar_skim.root");
    }
  }
  
  TEntryList *entrylist_ttjets = NULL;//new TEntryList("entrylist_ttjets", "entrylist_ttjets");
  vector<double>   *vjtDiscri_pf_ttjets = 0;
  
  TChain  *wjets_Muon = new TChain("AnaTree");
  TChain  *wjets_Electron = new TChain("AnaTree");
  if (!reloadHisto) {
    if (runOnPNFS) {
      wjets->Add("dcap:///wjetsTree.root");
        //    wjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/WJets_skim.root");
    } else {
      wjets_Muon->Add("/user/nstrobbe/DATA/stop/skims_443/WJets_skim.root");
      wjets_Electron->Add("/user/nstrobbe/DATA/stop/skims_443/WJets_skim.root");
    }
  }
  
  TEntryList *entrylist_wjets  = NULL;//new TEntryList("entrylist_wjets", "entrylist_wjets");
  vector<double>   *vjtDiscri_pf_wjets  = 0;
  vector<double>   *vjtFlav_pf_wjets    = 0;
  vector<double>   *vjtFlav_pf_zjets    = 0;
  vector<double>   *vjtFlav_pf_ttjets    = 0;
  
  TChain  *zjets_Muon = new TChain("AnaTree");
  TChain  *zjets_Electron = new TChain("AnaTree");
  if (!reloadHisto) {
    if (runOnPNFS) {
      zjets->Add("dcap:///zjetsTree.root");
        //    zjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/ZJets_skim.root");
    } else {
      zjets_Muon->Add("/user/nstrobbe/DATA/stop/skims_443/ZJets_skim.root");
      zjets_Electron->Add("/user/nstrobbe/DATA/stop/skims_443/ZJets_skim.root");
    }
  }
  
  TEntryList *entrylist_zjets  = NULL;//new TEntryList("entrylist_zjets", "entrylist_zjets");
  vector<double>   *vjtDiscri_pf_zjets  = 0;
  
  TChain *stjets_Muon = new TChain("AnaTree");
  TChain *stjets_Electron = new TChain("AnaTree");
  if (!reloadHisto) {
    if (runOnPNFS) {
      stjets->Add("dcap:///stjetsTree.root");
        //    stjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/T_tW_skim.root");
        //    stjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/Tbar_tW_skim.root");
    } else {
      stjets_Muon->Add("/user/nstrobbe/DATA/stop/skims_443/T_tW_skim.root");
      stjets_Electron->Add("/user/nstrobbe/DATA/stop/skims_443/T_tW_skim.root");
      stjets_Muon->Add("/user/nstrobbe/DATA/stop/skims_443/Tbar_tW_skim.root");
      stjets_Electron->Add("/user/nstrobbe/DATA/stop/skims_443/Tbar_tW_skim.root");
    }
  }
  
  TEntryList *entrylist_stjets = NULL;//new TEntryList("entrylist_stjets", "entrylist_stjets");
  vector<double>   *vjtDiscri_pf_stjets = 0;
  
    //  TChain *vbjets = new TChain("AnaTree");
    //  vbjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/XXX_skim.root");
    //  TEntryList *entrylist_vbjets = 0;
    //  vector<double>   *vjtDiscri_pf_vbjets = 0;
  
    //  TBranch        *b_vjtDiscri_pf = 0;
  
  
  
  
  
  if(verbose && !reloadHisto){
    cout<<"Number of events (entries in chain) for :"<<endl;
    cout<<" - ttjets = "<<ttjets_Muon->GetEntries()<< "   ";
    cout<<" - ttjets = "<<ttjets_Electron->GetEntries()<<endl;
    cout<<" -  wjets = "<< wjets_Muon->GetEntries()<< "   ";
    cout<<" -  wjets = "<< wjets_Electron->GetEntries()<<endl;
    cout<<" -  zjets = "<< zjets_Muon->GetEntries()<< "   ";
    cout<<" -  zjets = "<< zjets_Electron->GetEntries()<<endl;
    cout<<" - stjets = "<<stjets_Muon->GetEntries()<< "   ";
    cout<<" - stjets = "<<stjets_Electron->GetEntries()<<endl;
      //    cout<<" - vbjets ="<<vbjets->GetEntries()<<endl;
  }
  
    // Processes XS (pb) taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSections
  /*** TO BE UPDATED TO THE MEASURED XS WHENEVER POSSIBLE ***/
  float XS_ttjets = 165.; //NNLL resummations
  float XS_wjets  = 31314.;
  float XS_zjets  = 3048.;
  //  float XS_stjets = 10.6 + 4.6/*NNLLresummations*/ + 64.6;
  /* //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSections
  float XS_stjets_s = 4.6; float XS_stjets_s_err = 0.06; //NNLL resummations
  float XS_stjets_t = 64.6; float XS_stjets_t_err = 3.4;
  float XS_stjets_twDR = 10.6;    float XS_stjets_twDR_err = 0.8;
*/
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopSigma
  float XS_stjets = 15.74 + 4.63 + 64.57;
  float XS_stjets_s = 3.19 ; float XS_stjets_s_bar = 1.44 ; float XS_stjets_s_err = 0.07+0.19;
  float XS_stjets_t = 41.92 ;  float XS_stjets_t_bar = 22.65 ; float XS_stjets_t_err = 2.09+1.74;
  float XS_stjets_twDR = 15.74;    float XS_stjets_twDR_err = 0.4+1.14;
  /*
    //Kidonakis : 1109.3231 ; take only the biggest error when asymmetric
      float XS_stjets_s = 3.17+1.42; float XS_stjets_s_err = 0.06+0.13 + 0.01+0.07;
float XS_stjets_t = 41.7+22.5; float XS_stjets_t_err = 1.6+0.8 + 0.5+0.9;
float XS_stjets_twDR = 15.6;    float XS_stjets_twDR_err = 0.2+0.6;
*/

  //  float XS_vbjets = 35.5; //OLD
  float XS_vbjets = 35.3; //OLD
  
  
  printf("AnaTree objects instanciated\n");
  
  AnaTree_Muon *anaMu_data   = NULL;
  AnaTree_Muon *anaMu_ttjets = NULL;
  AnaTree_Muon *anaMu_wjets  = NULL;
  AnaTree_Muon *anaMu_zjets  = NULL;
  AnaTree_Muon *anaMu_stjets = NULL;
  if (semiMuon && !reloadHisto) {
    printf("a\n");
    if (isOnData) {
      anaMu_data   = new AnaTree_Muon((TTree*) dataTree_Muon);
        //      anaMu_data->Init((TTree*) dataTree_Muon);
        //      printf("Init forced for anaMu_data\n");
    }
    printf("b\n");
    anaMu_ttjets = new AnaTree_Muon((TTree*) ttjets_Muon);
      //    anaMu_ttjets->Init((TTree*) dataTree);
    printf("c\n");
    
    anaMu_wjets  = new AnaTree_Muon((TTree*) wjets_Muon);
    printf("d\n");
    anaMu_zjets  = new AnaTree_Muon((TTree*) zjets_Muon);
    printf("e\n");
    anaMu_stjets = new AnaTree_Muon((TTree*) stjets_Muon);
    printf("f\n");
  }
  
  printf("AnaTree objects instanciated\n");
  AnaTree_Electrons *anaEl_data   = NULL;
  AnaTree_Electrons *anaEl_ttjets = NULL;
  AnaTree_Electrons *anaEl_wjets  = NULL;
  AnaTree_Electrons *anaEl_zjets  = NULL;
  AnaTree_Electrons *anaEl_stjets = NULL;
  if (semiElectron && !reloadHisto) {
    if (isOnData) {
      anaEl_data   = new AnaTree_Electrons((TTree*) dataTree_Electron);
    }
    anaEl_ttjets = new AnaTree_Electrons((TTree*) ttjets_Electron);
    anaEl_wjets  = new AnaTree_Electrons((TTree*) wjets_Electron);
    anaEl_zjets  = new AnaTree_Electrons((TTree*) zjets_Electron);
    anaEl_stjets = new AnaTree_Electrons((TTree*) stjets_Electron);
  }
  
  
  printf("AnaTree objects instanciated\n");
  
  const int NbOfPE = 10000;
  const int NbOfJetBins = 1;
  const int MinNbOfJets = 4;
  const bool isLastBinExclusive = true;
  const int NbOfWP = 2;
  const int muonMask = (1<<0);
  const int electronMask = (1<<1);
  const int NbOfChannels = muonMask + electronMask;
  
    // Normalization factors for a given int. lumi.
  /** TO BE ADAPTED **/
  bool skimmedTree = false;
  double* NormFact_ttjets = new double[NbOfChannels];
  double* NormFact_wjets = new double[NbOfChannels];
  double* NormFact_zjets = new double[NbOfChannels];
  double* NormFact_stjets = new double[NbOfChannels];
  double* NormFact_qcd = new double[NbOfChannels];
  double* NormFact_vvjets = new double[NbOfChannels];
  for (int c=0; c<NbOfChannels; c++) {
    NormFact_ttjets[c] = 0.;
    NormFact_wjets[c] = 0.;
    NormFact_zjets[c] = 0.;
    NormFact_stjets[c] = 0.;
    NormFact_qcd[c] = 0.;
    NormFact_vvjets[c] = 0.;

  }
    //  double NormFact_vbjets = 0;
  
  double* preselEff_ttjets_NUM = new double[NbOfChannels] ;
  double* preselEff_ttjets_DEN = new double[NbOfChannels] ;
  double* preselEff_wjets_NUM = new double[NbOfChannels] ;
  double* preselEff_wjets_DEN = new double[NbOfChannels] ;
  double* preselEff_zjets_NUM = new double[NbOfChannels] ;
  double* preselEff_zjets_DEN = new double[NbOfChannels] ;
  double* preselEff_stjets_s_NUM = new double[NbOfChannels] ;  double* preselEff_stjets_s_NUM_bar = new double[NbOfChannels] ;
  double* preselEff_stjets_s_DEN = new double[NbOfChannels] ;  double* preselEff_stjets_s_DEN_bar = new double[NbOfChannels] ;
  double* preselEff_stjets_t_NUM = new double[NbOfChannels] ;  double* preselEff_stjets_t_NUM_bar = new double[NbOfChannels] ;
  double* preselEff_stjets_t_DEN = new double[NbOfChannels] ;  double* preselEff_stjets_t_DEN_bar = new double[NbOfChannels] ;
  double* preselEff_stjets_twDR_NUM = new double[NbOfChannels] ;
  double* preselEff_stjets_twDR_DEN = new double[NbOfChannels] ;
  for (int c=0; c<NbOfChannels; c++) {
    preselEff_ttjets_NUM[c] = 0. ;
    preselEff_ttjets_DEN[c] = 0. ;
    preselEff_wjets_NUM[c] = 0. ;
    preselEff_wjets_DEN[c] = 0. ;
    preselEff_zjets_NUM[c] = 0. ;
    preselEff_zjets_DEN[c] = 0. ;
    preselEff_stjets_s_NUM[c] = 0. ;    preselEff_stjets_s_NUM_bar[c] = 0. ;
    preselEff_stjets_s_DEN[c] = 0. ;    preselEff_stjets_s_DEN_bar[c] = 0. ;
    preselEff_stjets_t_NUM[c] = 0. ;    preselEff_stjets_t_NUM_bar[c] = 0. ;
    preselEff_stjets_t_DEN[c] = 0. ;    preselEff_stjets_t_DEN_bar[c] = 0. ;
    preselEff_stjets_twDR_NUM[c] = 0. ;
    preselEff_stjets_twDR_DEN[c] = 0. ;
  }
  
  
  /* //OLD skim trees
   double preselEff_ttjets = (1.*(664244+676929+665934+675885+677845+678195+677301+677914+677282+677672+675988+540827))/(4906594+5000000+4917487+5000000+5000000+5000000+5000000+5000000+5000000+5000000+5000000+4000000);
   double preselEff_wjets = (1.*(5028+4996+5019+5022+4947+5011+4970+4707+4930+5063+5078+4624+4980+4997+5168+4975+1575))/(5000000+5000000+5000000+5000000+5000000+5000000+5000000+4727845+5000000+5000000+5000000+4599579+5000000+5000000+5000000+5000000+1500000);
   double preselEff_zjets = (1.*(8791+8048+8776+8890+8729+8846+8807+2643))/(5000000+4599900+5000000+5000000+5000000+5000000+5000000+1459126);
   double preselEff_stjets_s_NUM = (1.*6784);   double preselEff_stjets_s_DEN =259971;
   double preselEff_stjets_t_NUM = (1.*42385);   double preselEff_stjets_t_DEN =2700161;
   double preselEff_stjets_twDR_NUM =  (1.*18970);   double preselEff_stjets_twDR_DEN =1200000;
   */
    //Strong trees
  /*  
   preselEff_ttjets_NUM[muonMask-1] = (1.*(94717+97788+98252+98149+94935+97817+98246+97704+98517+97593+98215+98204+98037+98538+98348+98196+97980+98014+97852+97819+97504+97999+98023+59281));
   preselEff_ttjets_DEN[muonMask-1] = 1. * (2406594+2500000+2500000+2500000+2417487+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+1500000);
   preselEff_wjets_NUM[muonMask-1] = (1.*(383+365+383+375+383+358+380+380+374+372+363+381+367+349+357+321+392+365+380+366+365+364+307+344+382+360+383+403+371+403+417+384+242));
   preselEff_wjets_DEN[muonMask-1] = 1.*(2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2227845+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2099579+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+1500000);
   preselEff_zjets_NUM[muonMask-1] = (1.*(282+253+229+287+283+300+301+313+284+255+284+262+256+300+140));
   preselEff_zjets_DEN[muonMask-1] = 1.*(2500000+2500000+2099900+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+1459126);
   preselEff_stjets_s_NUM[muonMask-1] = (1.*(545));   preselEff_stjets_s_DEN[muonMask-1] =259971;
   preselEff_stjets_t_NUM[muonMask-1] = (1.*(2028+1657+1674));   preselEff_stjets_t_DEN[muonMask-1] =(1500000+1200161+1200161);
   preselEff_stjets_twDR_NUM[muonMask-1] = 13988 ;   preselEff_stjets_twDR_DEN[muonMask-1] =814390;
   */
  
  /* elFixed, relIso cut at skim level */
  /* find . -name "dat_*.txt" -exec grep "Selection" -H {} \; | grep "Zjets" | grep "muon" | sed "s/.* \([0-9]*\)\/\([0-9]*\)/\1+/"    */
  preselEff_ttjets_NUM[muonMask-1] = 97124+97120+97481+97580+97560+78229+94301+94053+97154+97397+97595+97674+97161+97498+97346+97372+96925+97610+97902+97603+97883+96850+97343+97201;
  preselEff_ttjets_DEN[muonMask-1] = 2500000+2500000+2500000+2500000+2500000+2000000+2417487+2406594+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000;
  preselEff_wjets_NUM[muonMask-1] = 367+372+365+241+383+358+355+383+371+305+361+402+364+380+382+321+344+358+360+377+384+365+378+365+392+382+374+417+374+383+364+402+371;
  preselEff_wjets_DEN[muonMask-1] = 2500000+2500000+2500000+1500000+2500000+2500000+2227845+2500000+2500000+2099579+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000;
  /**
Old Zjets with bug of electrons not present in endcaps
  preselEff_zjets_NUM[muonMask-1] = 297+299+252+165+239+284+253+264+277+289+146+284+315+249+304;
  preselEff_zjets_DEN[muonMask-1] = 2500000+2500000+2500000+1500000+2099900+2500000+2500000+2500000+2500000+2500000+1459126+2500000+2500000+2500000+2500000;
  */
  preselEff_zjets_NUM[muonMask-1] = 150+143+98+130+159+149+140+116+143+132+118+142+129+120+116+134+158+148+156+128+139+132+138+143+156+145+156+143+146;
  preselEff_zjets_DEN[muonMask-1] = 1250000+1250000+1099900+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1209126+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000;

  preselEff_stjets_s_NUM[muonMask-1] = 544 ;
  preselEff_stjets_s_NUM_bar[muonMask-1] = 544 ;
  preselEff_stjets_s_DEN[muonMask-1] = 259971 ;
  preselEff_stjets_s_DEN_bar[muonMask-1] = 137980 ;
  preselEff_stjets_t_NUM[muonMask-1] = 1671+2024+1656 ;
  preselEff_stjets_t_NUM_bar[muonMask-1] = 1671+2024+1656 ;
  preselEff_stjets_t_DEN[muonMask-1] = 1200000+1500000+1200161 ;
  preselEff_stjets_t_DEN_bar[muonMask-1] = 1944822 ;
  preselEff_stjets_twDR_NUM[muonMask-1] = 13922 ;
  preselEff_stjets_twDR_DEN[muonMask-1] = 814390 ;
  
  
  
    ///  OLD ELECTRON TREES
  /*  cat selection_efficiency_skim_electron.txt | grep "Zjets" | grep "electron" | sed -E "s/^.* : ([0-9]+)\/([0-9]+)$/\2 + /"*/
  /*
   preselEff_ttjets_NUM[electronMask-1] = 335787+335680+335449+334903+335879+267695+328718+334078+328928+335000+334647+334572 ;
   preselEff_ttjets_DEN[electronMask-1] = 5000000+5000000+5000000+5000000+5000000+4000000+4906594+5000000+4917487+5000000+5000000+5000000 ;
   preselEff_wjets_NUM[electronMask-1] = 2621+2646+2782+2694+2694+2588+815+2820+2652+2770+2499+2685+2679+2656+2728+2629+2514 ;
   preselEff_wjets_DEN[electronMask-1] = 5000000+5000000+5000000+5000000+5000000+5000000+1500000+5000000+5000000+5000000+4599579+5000000+5000000+5000000+5000000+5000000+4727845 ;
   preselEff_zjets_NUM[electronMask-1] = 6088+1939+6553+6435+6533+6546+6566+6645 ;
   preselEff_zjets_DEN[electronMask-1] = 4599900+1459126+5000000+5000000+5000000+5000000+5000000+5000000 ;
   preselEff_stjets_s_NUM[electronMask-1] = 2264 ;   preselEff_stjets_s_DEN[electronMask-1] = 259971 ;
   preselEff_stjets_t_NUM[electronMask-1] = 17691+7935 ;   preselEff_stjets_t_DEN[electronMask-1] = 2700161+1200000 ;
   preselEff_stjets_twDR_NUM[electronMask-1] = 34229 ;   preselEff_stjets_twDR_DEN[electronMask-1] = 814390 ;
   */
  
  /* elFixed, relIso cut at skim level */
  /* find . -name "dat_*.txt" -exec grep "Selection" -H {} \; | grep "Zjets" | grep "muon" | sed "s/.* \([0-9]*\)\/\([0-9]*\)/\1+/"    */
  preselEff_ttjets_NUM[electronMask-1] = 100933+101122+101132+101105+100861+80629+97432+97756+101130+100335+101323+101156+101073+101123+100853+101027+101620+101143+101489+101163+101259+101485+101113+101111;
  preselEff_ttjets_DEN[electronMask-1] = 2500000+2500000+2500000+2500000+2500000+2000000+2417487+2406594+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000;
  preselEff_wjets_NUM[electronMask-1] = 353+376+387+214+381+359+295+332+371+295+364+381+334+400+381+352+384+373+361+375+363+398+359+335+372+377+352+389+379+385+372+379+362;
  preselEff_wjets_DEN[electronMask-1] = 2500000+2500000+2500000+1500000+2500000+2500000+2227845+2500000+2500000+2099579+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000+2500000;
  /**
Old Zjets with bug of electrons not present in endcaps
  preselEff_zjets_NUM[electronMask-1] = 423+449+490+244+378+437+471+487+461+421+282+459+423+438+444;
  preselEff_zjets_DEN[electronMask-1] = 2500000+2500000+2500000+1500000+2099900+2500000+2500000+2500000+2500000+2500000+1459126+2500000+2500000+2500000+2500000;
  */
  preselEff_zjets_NUM[electronMask-1] = 305+294+243+286+329+275+291+287+298+277+275+284+290+292+295+295+267+305+276+284+315+292+287+266+291+295+299+285+283;
  preselEff_zjets_DEN[electronMask-1] = 1250000+1250000+1099900+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1209126+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000+1250000;
  
  
  preselEff_stjets_s_NUM[electronMask-1] = 666 ;
  preselEff_stjets_s_NUM_bar[muonMask-1] = 666 ;
  preselEff_stjets_s_DEN[electronMask-1] = 259971 ;
  preselEff_stjets_s_DEN_bar[electronMask-1] = 137980 ;
  preselEff_stjets_t_NUM[electronMask-1] = 1963+2325+2001 ;
  preselEff_stjets_t_NUM_bar[electronMask-1] = 1963+2325+2001 ;
  preselEff_stjets_t_DEN[electronMask-1] = 1200000+1500000+1200161 ;
  preselEff_stjets_t_DEN_bar[electronMask-1] = 1944822 ;
  preselEff_stjets_twDR_NUM[electronMask-1] = 14038 ;
  preselEff_stjets_twDR_DEN[electronMask-1] = 814390 ;
  
  for (int c=0; c<NbOfChannels; c++) {
    if (((c+1) & muonMask)!=0 && (c+1)!=muonMask) {
      preselEff_ttjets_NUM[c] += preselEff_ttjets_NUM[muonMask-1] ;
      preselEff_ttjets_DEN[c] += preselEff_ttjets_DEN[muonMask-1] ;
      preselEff_wjets_NUM[c] += preselEff_wjets_NUM[muonMask-1] ;
      preselEff_wjets_DEN[c] += preselEff_wjets_DEN[muonMask-1] ;
      preselEff_zjets_NUM[c] += preselEff_zjets_NUM[muonMask-1] ;
      preselEff_zjets_DEN[c] += preselEff_zjets_DEN[muonMask-1] ;
      preselEff_stjets_s_NUM[c] += preselEff_stjets_s_NUM[muonMask-1] ;
      preselEff_stjets_s_DEN[c] += preselEff_stjets_s_DEN[muonMask-1] ;       preselEff_stjets_s_DEN_bar[c] += preselEff_stjets_s_DEN_bar[muonMask-1] ;
      preselEff_stjets_t_NUM[c] += preselEff_stjets_t_NUM[muonMask-1] ;
      preselEff_stjets_t_DEN[c] += preselEff_stjets_t_DEN[muonMask-1] ;       preselEff_stjets_t_DEN_bar[c] += preselEff_stjets_t_DEN_bar[muonMask-1] ;
      preselEff_stjets_twDR_NUM[c] += preselEff_stjets_twDR_NUM[muonMask-1] ;
      preselEff_stjets_twDR_DEN[c] += preselEff_stjets_twDR_DEN[muonMask-1] ;
    }
    if (((c+1) & electronMask)!=0 && (c+1)!=electronMask) {
      preselEff_ttjets_NUM[c] += preselEff_ttjets_NUM[electronMask-1] ;
      preselEff_ttjets_DEN[c] += preselEff_ttjets_DEN[electronMask-1] ;
      preselEff_wjets_NUM[c] += preselEff_wjets_NUM[electronMask-1] ;
      preselEff_wjets_DEN[c] += preselEff_wjets_DEN[electronMask-1] ;
      preselEff_zjets_NUM[c] += preselEff_zjets_NUM[electronMask-1] ;
      preselEff_zjets_DEN[c] += preselEff_zjets_DEN[electronMask-1] ;
      preselEff_stjets_s_NUM[c] += preselEff_stjets_s_NUM[electronMask-1] ;
      preselEff_stjets_s_DEN[c] += preselEff_stjets_s_DEN[electronMask-1] ;      preselEff_stjets_s_DEN_bar[c] += preselEff_stjets_s_DEN_bar[electronMask-1] ;
      preselEff_stjets_t_NUM[c] += preselEff_stjets_t_NUM[electronMask-1] ;
      preselEff_stjets_t_DEN[c] += preselEff_stjets_t_DEN[electronMask-1] ;      preselEff_stjets_t_DEN_bar[c] += preselEff_stjets_t_DEN_bar[electronMask-1] ;
      preselEff_stjets_twDR_NUM[c] += preselEff_stjets_twDR_NUM[electronMask-1] ;
      preselEff_stjets_twDR_DEN[c] += preselEff_stjets_twDR_DEN[electronMask-1] ;
    }
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
    //WW : (1.*(9030+9181+9083+8362+2624))/(1000000+1000000+1000000+925857+300000)
    //WZ : (1.*(8176+7778+8221+7401+3405))/(1000000+949999+1000000+915172+400000)
    //ZZ : (1.*(6549+6488+6585+6598+1295))/(1000000+1000000+1000000+990973+200000)
  
  
  
  
  
  printf("Event weights for skimmed trees\n");
  
  if (!reloadHisto) {
    if(skimmedTree) {
      ttjets_Muon->SetBranchAddress("eventWeight",&NormFact_ttjets);
        //      ttjets_Electron->SetBranchAddress("eventWeight",&NormFact_ttjets_Electron);
      wjets_Muon->SetBranchAddress("eventWeight",&NormFact_wjets);
      zjets_Muon->SetBranchAddress("eventWeight",&NormFact_zjets);
      stjets_Muon->SetBranchAddress("eventWeight",&NormFact_stjets);
        //vbjets->SetBranchAddress("eventWeight",&NormFact_vbjets);
    }
    else{
      for (int c=0; c<NbOfChannels; c++) {
        /* Preselection (skimming) efficiency has to be accounted for at filling of the histograms. */
        NormFact_ttjets[c] = XS_ttjets/((((c+1)&muonMask!=0?1:0)*ttjets_Muon->GetEntries())
                                        +(((c+1)&electronMask!=0?1:0)*ttjets_Electron->GetEntries()));
          //      NormFact_ttjets_Electron = preselEff_ttjets * XS_ttjets/(ttjets_Electron->GetEntries());
        NormFact_wjets[c]  = XS_wjets/((((c+1)&muonMask!=0?1:0)*wjets_Muon->GetEntries())
                                       +(((c+1)&electronMask!=0?1:0)*wjets_Electron->GetEntries()));
        NormFact_zjets[c]  = XS_zjets/((((c+1)&muonMask!=0?1:0)*zjets_Muon->GetEntries())
                                       +(((c+1)&electronMask!=0?1:0)*zjets_Electron->GetEntries()));
        NormFact_stjets[c] = XS_stjets/((((c+1)&muonMask!=0?1:0)*stjets_Muon->GetEntries())
                                        +(((c+1)&electronMask!=0?1:0)*stjets_Electron->GetEntries()));
          //NormFact_vbjets = XS_vbjets/vbjets->GetEntries();
      }
    }
  }
  
  
  TFile* histoFile = new TFile(histoFileName.c_str(), "read");
  
  
  printf("Parameters jobs definition\n");
  
    // E v e n t   s e l e c t i o n   c u t s
    // ---------------------------------------
  
  
  
    // Retrieve int. lumi from InputTree (?)
  double *IntLumi = new double[NbOfChannels];
  for (int c=0; c<NbOfChannels; c++) {
    IntLumi[c] = 0.;
  }
  for (int c=0; c<NbOfChannels; c++) {
    if (((c+1) & muonMask)!=0) {
      IntLumi[c] += 5049.92;
    }
    if (((c+1) & electronMask)!=0) {
      IntLumi[c] += 5050.821;
    }
  }
  
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
  
  char** channelSuffix = new char*[NbOfChannels];
  for (int i=0; i<NbOfChannels; i++) {
    channelSuffix[i] = new char[50];
    sprintf(channelSuffix[i], "_%s", (std::string()
                                      +(((i+1)&(1<<0))?"SemiMuon":"")
                                      +(((i+1)&(1<<1))?"SemiElectron":"")).c_str());
  }    
  
  
  char** jetSuffix = new char*[NbOfJetBins];
  char** jtSuf = new char*[NbOfJetBins];
  for (int i=0; i<NbOfJetBins; i++) {
    jetSuffix[i] = new char[50];
    jtSuf[i] = new char[5];
    if (i==NbOfJetBins-1 && !isLastBinExclusive) {
      sprintf(jetSuffix[i], "_geq%djets", MinNbOfJets+i);
      sprintf(jtSuf[i], "_%dj", MinNbOfJets+i);
    } else {
      sprintf(jetSuffix[i], "_%djets", MinNbOfJets+i);
      sprintf(jtSuf[i], "_%dj", MinNbOfJets+i);
    }
  }    
  
  char** wpSuffix = new char*[NbOfWP];
  char** wpSuf = new char*[NbOfWP];
  for (int i=0; i<NbOfWP; i++) {
    wpSuffix[i] = new char[50];
    wpSuf[i] = new char[5];
    sprintf(wpSuffix[i], "_WP%d", i);
    sprintf(wpSuf[i], "_wp%d", i);
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
    //    TCut cuts = (skimmedTree ? pfnjet[JetIdx] : singlep+pfnjet[JetIdx]);
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
  
  for (int c=0; c<NbOfChannels; c++) {
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
    
      // Applying the scale factor if from MC
      // Applying the value if from Data
    /**
     Results from BTV-11-003, https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG link to https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/eff_b_c-ttbar_payload.txt 
     In Range : 0<=|eta|<=2.4 , 30<=p_T<=200  (approximately)
     */
    /** TCHE */
    for (int j=0; j<NbOfWP; j++) {
      if (isOnData) {
        Eb_const_mean[c][j] = eff_b_ttbar_FtCM(btagCuts[j]);
        Eb_const_error[c][j] = eff_b_ttbar_FtCM_pluserr(btagCuts[j])-eff_b_ttbar_FtCM(btagCuts[j]);
      }
      else {
        Eb_const_mean[c][j] = eff_b_ttbar_MC(btagCuts[j]);
        Eb_const_error[c][j] = Eb_const_mean[c][j] * 0.05;
      }
    }
    for (int j=0; j<NbOfWP; j++) {
        //Eb_const_mean[j]     = 0.75;//, 0.795};     /** TO BE CALCULATED BY FUNCT FROM BTV-11-003 **/
        //    Eb_const_error[j]    = Eb_const_mean*0.1;//, 0.795*0.1}; /** TO BE CALCULATED BY FUNCT FROM BTV-11-003 **/
      Euds_const_mean[c][j]   = Eb_const_mean[c][j] / 6.; //, 0.000};     /** TO BE DERIVED FROM MC **/
      Euds_const_error[c][j]  = Euds_const_mean[c][j] * .15;//, 0.000};     /** TO BE DERIVED FROM MC **/
      Eudsc_const_mean[c][j]  = Euds_const_mean[c][j] * 1.3;//, 0.179};     /** TO BE CALCULATED BY FUNCT FROM BTV-11-003 **/
      Eudsc_const_error[c][j] = Eudsc_const_mean[c][j] * .15;//, 0.179*0.1}; /** TO BE CALCULATED BY FUNCT FROM BTV-11-003 **/
      
      init_Eb[c][j] = Eb_const_mean[c][j];
      init_Euds[c][j] = Euds_const_mean[c][j];
      init_Eudsc[c][j] = Eudsc_const_mean[c][j];
    }
  }
  
  
  printf("Defining (or reloading) histograms and useful MC information\n");
  
  char name[50];
  TH1F*** eXbq__ = new TH1F**[NbOfChannels];
  for (int c=0; c<NbOfChannels ; c++) {
    eXbq__[c] = new TH1F*[NbOfJetBins];
    if (reloadHisto) {
      for (int i=0; i<NbOfJetBins; i++) {
        eXbq__[c][i] = (TH1F*) histoFile->Get((std::string()+"eXbq"+channelSuffix[c]+jetSuffix[i]).c_str());
      }
    } else {
      for (int i=0; i<NbOfJetBins; i++) {
        eXbq__[c][i] = new TH1F((std::string()+"eXbq"+channelSuffix[c]+jetSuffix[i]).c_str(), (std::string()+"eXbq"+channelSuffix[c]+jetSuffix[i]).c_str(), 4, 0, 4);
      }
    }
  }
  
  
  TH1I* histoTrig_ttjets = new TH1I("TrigPattern_ttjets","Pattern trigger for ttjets;Trigger pathes combination",512,0,512);
  TH1I* histoTrig_wjets  = new TH1I("TrigPattern_wjets_lightAndWbb","Pattern trigger for wjets (light + Wbb);Trigger pathes combination",512,0,512);
  TH1I* histoTrig_zjets  = new TH1I("TrigPattern_zjets","Pattern trigger for zjets;Trigger pathes combination",512,0,512);
  TH1I* histoTrig_stjets = new TH1I("TrigPattern_stjets","Pattern trigger for stjets;Trigger pathes combination",512,0,512);
  TH1I** histoNJ_data   = new TH1I*[NbOfChannels]; for (int c=0; c<NbOfChannels; c++) { histoNJ_data[c]   = new TH1I((std::string()+"JetMult_data"+channelSuffix[c]).c_str(),"JetMult_datas",50,0,50); }
  TH1I** histoNJ_ttjets = new TH1I*[NbOfChannels]; for (int c=0; c<NbOfChannels; c++) { histoNJ_ttjets[c] = new TH1I((std::string()+"JetMult_ttjets"+channelSuffix[c]).c_str(),"JetMult_ttjets",50,0,50); }
  TH1I** histoNJ_wjets  = new TH1I*[NbOfChannels]; for (int c=0; c<NbOfChannels; c++) { histoNJ_wjets[c]  = new TH1I((std::string()+"JetMult_wjets"+channelSuffix[c]).c_str(),"JetMult_wjets",50,0,50); }
  TH1I** histoNJ_zjets  = new TH1I*[NbOfChannels]; for (int c=0; c<NbOfChannels; c++) { histoNJ_zjets[c]  = new TH1I((std::string()+"JetMult_zjets"+channelSuffix[c]).c_str(),"JetMult_zjets",50,0,50); }
  TH1I** histoNJ_stjets = new TH1I*[NbOfChannels]; for (int c=0; c<NbOfChannels; c++) { histoNJ_stjets[c] = new TH1I((std::string()+"JetMult_stjets"+channelSuffix[c]).c_str(),"JetMult_stjets",50,0,50); }
  
    // Counters for the b-tagged jet multiplicity
  int NbtaggedJets = 0;//new int[NbOfChannels];
  int nB_Flavoured_Jets = 0;//new int[NbOfChannels];
  int nC_Flavoured_Jets = 0;//new int[NbOfChannels];
  int nLightAndGluon_Flavoured_Jets = 0;//new int[NbOfChannels];
                                        // Histograms for the b-tagged jet multiplicity
  TH1F*** hMC_Nominal_N = new TH1F**[NbOfChannels]; for (int c=0; c<NbOfChannels; c++) { hMC_Nominal_N[c] = new TH1F*[NbOfJetBins]; for (int i=0; i<NbOfJetBins; i++) {
                                                                                           hMC_Nominal_N[c][i] = new TH1F((std::string()+"MC_Nominal_N"+channelSuffix[c]+jetSuffix[i]).c_str() ,"MC_Nominal_N;", 20, 0, 20);
                                                                                         }
  }
  
 TH1F*** hMC_Nominal_E = new TH1F**[NbOfChannels]; for (int c=0; c<NbOfChannels; c++) { hMC_Nominal_E[c] = new TH1F*[NbOfWP]; for (int j=0; j<NbOfWP; j++) {
                                                                                          hMC_Nominal_E[c][j] = new TH1F((std::string()+"MC_Nominal_eff"+channelSuffix[c]+wpSuffix[j]).c_str() ,"MC_Nominal_eff", 20, 0, 20);
                                                                                        }
 }
 
  TH1F**** hNbtaggedJets_ttjets = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hNbtaggedJets_ttjets[c] = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_ttjets[c][i] = new TH1F*[NbOfWP]; } }
  TH1F**** hNbtaggedJets_wjets  = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hNbtaggedJets_wjets[c]  = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_wjets[c][i]  = new TH1F*[NbOfWP]; } }
  TH1F**** hNbtaggedJets_zjets  = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hNbtaggedJets_zjets[c]  = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_zjets[c][i]  = new TH1F*[NbOfWP]; } }
  TH1F**** hNbtaggedJets_stjets = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hNbtaggedJets_stjets[c] = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_stjets[c][i] = new TH1F*[NbOfWP]; } }
  TH1F**** hNbtaggedJets_vbjets = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hNbtaggedJets_vbjets[c] = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_vbjets[c][i] = new TH1F*[NbOfWP]; } }
  TH1F**** hNbtaggedJets_qcd    = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hNbtaggedJets_qcd[c]    = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_qcd[c][i] = new TH1F*[NbOfWP]; } }
  TH1F**** hNbtaggedJets_vvjets = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hNbtaggedJets_vvjets[c] = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) { hNbtaggedJets_vvjets[c][i] = new TH1F*[NbOfWP]; } }
  TH1F**** hWbb_In_WJets        = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hWbb_In_WJets[c]        = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) { hWbb_In_WJets[c][i]  = new TH1F*[NbOfWP]; } }
  TH1F**** hWcx_In_WJets        = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hWcx_In_WJets[c]        = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) { hWcx_In_WJets[c][i]  = new TH1F*[NbOfWP]; } }
  TH1F**** hWLFg_In_WJets       = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hWLFg_In_WJets[c]       = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) { hWLFg_In_WJets[c][i] = new TH1F*[NbOfWP]; } }
  TH1F**** hNbtaggedJets        = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hNbtaggedJets[c] = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) hNbtaggedJets[c][i] = new TH1F*[NbOfWP]; }
  for(int c=0; c<NbOfChannels; c++) { 
      //    NbtaggedJets[c] = 0;
      //    nB_Flavoured_Jets[c] = 0;
      //    nC_Flavoured_Jets[c] = 0;
      //    nLightAndGluon_Flavoured_Jets[c] = 0;
    for (int i=0; i<NbOfJetBins; i++) {
      for (int j=0; j<NbOfWP; j++) {
        hNbtaggedJets[c][i][j] = NULL;
        hNbtaggedJets_ttjets[c][i][j] = NULL;
        hNbtaggedJets_wjets[c][i][j] = NULL;
        hNbtaggedJets_zjets[c][i][j] = NULL;
        hNbtaggedJets_stjets[c][i][j] = NULL;
        hNbtaggedJets_vbjets[c][i][j] = NULL;
        hNbtaggedJets_qcd[c][i][j] = NULL;
        hNbtaggedJets_vvjets[c][i][j] = NULL;
        hWbb_In_WJets[c][i][j] = NULL;
        hWcx_In_WJets[c][i][j] = NULL;
        hWLFg_In_WJets[c][i][j] = NULL;
      }
    }
  }
  
  TH1F**** hNbtaggedJets_ttlike = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hNbtaggedJets_ttlike[c] = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) hNbtaggedJets_ttlike[c][i] = new TH1F*[NbOfWP]; }
  TH1F**** hNbtaggedJets_vlike  = new TH1F***[NbOfChannels]; for(int c=0; c<NbOfChannels; c++) { hNbtaggedJets_vlike[c] = new TH1F**[NbOfJetBins]; for(int i=0; i<NbOfJetBins ; i++) hNbtaggedJets_vlike[c][i] = new TH1F*[NbOfWP]; }
  for(int c=0; c<NbOfChannels; c++) { 
    for (int i=0; i<NbOfJetBins; i++) {
      for (int j=0; j<NbOfWP; j++) {
        hNbtaggedJets_ttlike[c][i][j] = new TH1F((std::string()+"hNbtaggedJets_ttlike"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),";B-tagged jet multiplicity",10,0,10);
        hNbtaggedJets_vlike[c][i][j]  = new TH1F((std::string()+"hNbtaggedJets_vlike"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str() ,";B-tagged jet multiplicity",10,0,10);
      }
    }
  }
  
  TH2F *ebVMC_histo = NULL; 
  TH2F *ebTTMC_histo = NULL;
  TH2F *eudsMC_histo = NULL; 
  TH2F *eudscMC_histo = NULL;
  {
    
    if (reloadHisto) {
      if (histoFile->IsZombie()) {
        printf("The Histo file is Zombie ---> unable to reload the histograms\n");
          //exit(1);
      } else {
        ebVMC_histo = (TH2F*) histoFile->Get("ebV");
        ebTTMC_histo = (TH2F*) histoFile->Get("ebTT");
        eudsMC_histo = (TH2F*) histoFile->Get("eudsMC");
        eudscMC_histo = (TH2F*) histoFile->Get("eudscMC");
        for(int c=0; c<NbOfChannels; c++) { 
          for (int i=0; i<NbOfJetBins; i++) {
            for (int j=0; j<NbOfWP; j++) {
              if (reloadHisto) {
                hNbtaggedJets[c][i][j] = (TH1F*) histoFile->Get((std::string()+"hNbtaggedJets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
              } else {
                hNbtaggedJets[c][i][j]        = new TH1F((std::string()+"hNbtaggedJets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str()       ,";B-tagged jet multiplicity",10,0,10);
              }
              if (!isOnData) {
                hNbtaggedJets[c][i][j]->Scale(IntLumi[c]);
              }
              hNbtaggedJets_ttjets[c][i][j] = (TH1F*) histoFile->Get((std::string()+"hNbtaggedJets_ttjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
              hNbtaggedJets_wjets[c][i][j]  = (TH1F*) histoFile->Get((std::string()+"hNbtaggedJets_wjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
              hNbtaggedJets_zjets[c][i][j]  = (TH1F*) histoFile->Get((std::string()+"hNbtaggedJets_zjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
              hNbtaggedJets_stjets[c][i][j] = (TH1F*) histoFile->Get((std::string()+"hNbtaggedJets_stjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
              hNbtaggedJets_vbjets[c][i][j] = (TH1F*) histoFile->Get((std::string()+"hNbtaggedJets_vbjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
              
              hWbb_In_WJets[c][i][j] = (TH1F*) histoFile->Get((std::string()+"hNbtaggedJets_wjets_Wbb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
              hWcx_In_WJets[c][i][j] = (TH1F*) histoFile->Get((std::string()+"hNbtaggedJets_wjets_Wcx"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
              hWLFg_In_WJets[c][i][j] = (TH1F*) histoFile->Get((std::string()+"hNbtaggedJets_wjets_WLFg"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
              
              hNbtaggedJets_ttjets[c][i][j]->Scale(IntLumi[c]);
              hNbtaggedJets_wjets[c][i][j]->Scale(IntLumi[c]);
              hNbtaggedJets_zjets[c][i][j]->Scale(IntLumi[c]);
              hNbtaggedJets_stjets[c][i][j]->Scale(IntLumi[c]);
              hNbtaggedJets_vbjets[c][i][j]->Scale(IntLumi[c]);
              hWbb_In_WJets[c][i][j]->Scale(IntLumi[c]);
              hWcx_In_WJets[c][i][j]->Scale(IntLumi[c]);
              hWLFg_In_WJets[c][i][j]->Scale(IntLumi[c]);
            }
          }
        }
      }
    } else {
      ebVMC_histo = new TH2F("ebV", "ebV", NbOfChannels, 0, NbOfChannels, NbOfWP, 0, NbOfWP);
      ebTTMC_histo = new TH2F("ebTT", "ebTT", NbOfChannels, 0, NbOfChannels, NbOfWP, 0, NbOfWP);
      eudsMC_histo = new TH2F("eudsMC", "eudsMC", NbOfChannels, 0, NbOfChannels, NbOfWP, 0, NbOfWP);
      eudscMC_histo = new TH2F("eudscMC", "eudscMC", NbOfChannels, 0, NbOfChannels, NbOfWP, 0, NbOfWP);
      for(int c=0; c<NbOfChannels; c++) { 
        for (int i=0; i<NbOfJetBins; i++) {
          for (int j=0; j<NbOfWP; j++) {
            hNbtaggedJets[c][i][j]        = new TH1F((std::string()+"hNbtaggedJets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str()       ,";B-tagged jet multiplicity",10,0,10);
            hNbtaggedJets_ttjets[c][i][j] = new TH1F((std::string()+"hNbtaggedJets_ttjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),";B-tagged jet multiplicity",10,0,10);
            hNbtaggedJets_wjets[c][i][j]  = new TH1F((std::string()+"hNbtaggedJets_wjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), ";B-tagged jet multiplicity",10,0,10);
            hNbtaggedJets_zjets[c][i][j]  = new TH1F((std::string()+"hNbtaggedJets_zjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(), ";B-tagged jet multiplicity",10,0,10);
            hNbtaggedJets_stjets[c][i][j] = new TH1F((std::string()+"hNbtaggedJets_stjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),";B-tagged jet multiplicity",10,0,10);
            hNbtaggedJets_vbjets[c][i][j] = new TH1F((std::string()+"hNbtaggedJets_vbjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),";B-tagged jet multiplicity",10,0,10);
            
            hWbb_In_WJets[c][i][j]  = new TH1F((std::string()+"hNbtaggedJets_wjets_Wbb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),";B-tagged jet multiplicity",10,0,10);
            hWcx_In_WJets[c][i][j]  = new TH1F((std::string()+"hNbtaggedJets_wjets_Wcx"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),";B-tagged jet multiplicity",10,0,10);
            hWLFg_In_WJets[c][i][j] = new TH1F((std::string()+"hNbtaggedJets_wjets_WLFg"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),";B-tagged jet multiplicity",10,0,10);
          }
        }
      }
    }
  }
    //  histoFile->Close();
  
  {
    // Loading Alberto's histograms
    const int NbOfFiles = 12;
    const char* albertoFileNames[NbOfFiles] = { "../Leg3ForVJets_20121024/QCDmu20_mu.root", "../Leg3ForVJets_20121024/QCDel30to80_el.root", "../Leg3ForVJets_20121024/QCDel80to170_el.root", "../Leg3ForVJets_20121024/QCDel170to250_el.root", "../Leg3ForVJets_20121024/QCDel250to350_el.root", "../Leg3ForVJets_20121024/QCDel350_el.root", "../Leg3ForVJets_20121024/WW_el.root", "../Leg3ForVJets_20121024/WW_mu.root", "../btagMiguel/WZ_el.root", "../Leg3ForVJets_20121024/WZ_mu.root", "../Leg3ForVJets_20121024/ZZ_el.root", "../Leg3ForVJets_20121024/ZZ_mu.root" };
    const float xs[NbOfFiles] = { 84679.3, 3625840., 142813.8, 3263.436, 368.01, 55.016, 43., 43., 18.2, 18.2, 5.9, 5.9 };
    
    TFile **alberto = new TFile*[NbOfFiles];
    char histoName[200] = "";
    for (int f=0; f<NbOfFiles; f++) {
      alberto[f] = new TFile(albertoFileNames[f]);
    }
    for(int c=0; c<NbOfChannels; c++) { 
      for (int i=0; i<NbOfJetBins; i++) {
        for (int j=0; j<NbOfWP; j++) {
          for (int f=0; f<NbOfFiles; f++) {
            if ((std::string(albertoFileNames[f]).find("_mu.root")!=std::string::npos) || std::string(albertoFileNames[f]).find("Mu.root")!=std::string::npos) {
              if ((((c+1) & muonMask) != 0) && semiMuon) {
                if (i+1==NbOfJetBins && !isLastBinExclusive) {
                  sprintf(histoName, "B_Jet_Multiplicity_Mu_%djInc", MinNbOfJets+i);
                } else {
                  sprintf(histoName, "B_Jet_Multiplicity_Mu_%djExc", MinNbOfJets+i);
                }
                printf("histoName[%d][%d][%d][%d] : %s    (%s)\n", c, i, j, f, histoName, albertoFileNames[f]);
                if (std::string(albertoFileNames[f]).find("QCD")!=std::string::npos) {
                  if (hNbtaggedJets_qcd[c][i][j]==NULL) {
                    hNbtaggedJets_qcd[c][i][j] = (TH1F*) alberto[f]->Get(histoName);
                    hNbtaggedJets_qcd[c][i][j]->SetDirectory(NULL);
                    hNbtaggedJets_qcd[c][i][j]->SetName((std::string()+"hNbtaggedJets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
                    hNbtaggedJets_qcd[c][i][j]->Reset();
                  }
                  hNbtaggedJets_qcd[c][i][j]->Add(((TH1F*)alberto[f]->Get(histoName)), IntLumi[muonMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2));
//                  for (int k=1; k<=hNbtaggedJets_qcd[c][i][j]->GetNbinsX() ; k++)
//                    hNbtaggedJets_qcd   [c][i][j]->SetBinContent(k, 1.);
                  printf("factor QCD [%d][%d][%d][%d] %f * %lf(%lf)  = %lf       ////    %lf\n", c, i, j, f, IntLumi[muonMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2), ((TH1F*)alberto[f]->Get(histoName))->Integral(), ((TH1F*)alberto[f]->Get(histoName))->GetEntries(), IntLumi[muonMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2)*((TH1F*)alberto[f]->Get(histoName))->Integral(), hNbtaggedJets_qcd[c][i][j]->Integral());
}
                if ((std::string(albertoFileNames[f]).find("WW")!=std::string::npos) || (std::string(albertoFileNames[f]).find("WZ")!=std::string::npos) || (std::string(albertoFileNames[f]).find("ZZ")!=std::string::npos)) {
                  if (hNbtaggedJets_vvjets[c][i][j]==NULL) {
                    hNbtaggedJets_vvjets[c][i][j] = (TH1F*) alberto[f]->Get(histoName);
                    hNbtaggedJets_vvjets[c][i][j]->SetDirectory(NULL);
                    hNbtaggedJets_vvjets[c][i][j]->SetName((std::string()+"hNbtaggedJets_vvjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
                    hNbtaggedJets_vvjets[c][i][j]->Reset();
                  }
                  hNbtaggedJets_vvjets[c][i][j]->Add(((TH1F*)alberto[f]->Get(histoName)), IntLumi[muonMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2));
//                  for (int k=1; k<=hNbtaggedJets_vvjets[c][i][j]->GetNbinsX() ; k++)
//                    hNbtaggedJets_vvjets[c][i][j]->SetBinContent(k, 1.);
                  printf("factor VV [%d][%d][%d][%d] %f * %lf(%lf)  = %lf       ////    %lf\n", c, i, j, f, IntLumi[muonMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2), ((TH1F*)alberto[f]->Get(histoName))->Integral(), ((TH1F*)alberto[f]->Get(histoName))->GetEntries(), IntLumi[muonMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2)*((TH1F*)alberto[f]->Get(histoName))->Integral(), hNbtaggedJets_vvjets[c][i][j]->Integral());
                }
              }
            } else { //if not an ONLY muon file, it is a ONLY electron file
              if ((((c+1) & electronMask) != 0) && semiElectron) {
                if (i+1==NbOfJetBins && !isLastBinExclusive) {
                  sprintf(histoName, "B_Jet_Multiplicity_El_%djInc", MinNbOfJets+i);
                } else {
                  sprintf(histoName, "B_Jet_Multiplicity_El_%djExc", MinNbOfJets+i);
                }
                printf("genericEl histoName[%d][%d][%d][%d] : %s    (%s)\n", c, i, j, f, histoName, albertoFileNames[f]);
                if (std::string(albertoFileNames[f]).find("QCD")!=std::string::npos) {
                  /*
  if (i+1==NbOfJetBins && !isLastBinExclusive) {
                    sprintf(histoName, "B_Jet_Multiplicity_Mu_%djInc", MinNbOfJets+i);
                  } else {
                    sprintf(histoName, "B_Jet_Multiplicity_Mu_%djExc", MinNbOfJets+i);
                  }
                printf("QCD histoName[%d][%d][%d][%d] : %s\n", c, i, j, f, histoName);
                printf("histoName : %s\n", histoName);
                  */
                  if (hNbtaggedJets_qcd[c][i][j]==NULL) {
                    hNbtaggedJets_qcd[c][i][j] = (TH1F*) alberto[f]->Get(histoName);
                    hNbtaggedJets_qcd[c][i][j]->SetDirectory(NULL);
                    hNbtaggedJets_qcd[c][i][j]->SetName((std::string()+"hNbtaggedJets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
                    hNbtaggedJets_qcd[c][i][j]->Reset();
                  }
                  hNbtaggedJets_qcd   [c][i][j]->Add(((TH1F*)alberto[f]->Get(histoName)), IntLumi[electronMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2));
//                  for (int k=1; k<=hNbtaggedJets_qcd[c][i][j]->GetNbinsX() ; k++)
//                    hNbtaggedJets_qcd   [c][i][j]->SetBinContent(k, 1.);
                  printf("factor QCD [%d][%d][%d][%d] %f * %lf(%lf)  = %lf       ////    %lf\n", c, i, j, f, IntLumi[electronMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2), ((TH1F*)alberto[f]->Get(histoName))->Integral(), ((TH1F*)alberto[f]->Get(histoName))->GetEntries(), IntLumi[electronMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2)*((TH1F*)alberto[f]->Get(histoName))->Integral(), hNbtaggedJets_qcd[c][i][j]->Integral());
                }
                if ((std::string(albertoFileNames[f]).find("WW")!=std::string::npos) || (std::string(albertoFileNames[f]).find("WZ")!=std::string::npos) || (std::string(albertoFileNames[f]).find("ZZ")!=std::string::npos)) {
                  if (hNbtaggedJets_vvjets[c][i][j]==NULL) {
                    hNbtaggedJets_vvjets[c][i][j] = (TH1F*) alberto[f]->Get(histoName);
                    hNbtaggedJets_vvjets[c][i][j]->SetDirectory(NULL);
                    hNbtaggedJets_vvjets[c][i][j]->SetName((std::string()+"hNbtaggedJets_vvjets"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str());
                    hNbtaggedJets_vvjets[c][i][j]->Reset();
                  }
                  hNbtaggedJets_vvjets[c][i][j]->Add(((TH1F*)alberto[f]->Get(histoName)), IntLumi[electronMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2));
//                  for (int k=1; k<=hNbtaggedJets_vvjets[c][i][j]->GetNbinsX() ; k++)
//                    hNbtaggedJets_vvjets[c][i][j]->SetBinContent(k, 1.);
                  printf("factor VV [%d][%d][%d][%d] %f * %lf(%lf)  = %lf\t       ///\t    %lf\n", c, i, j, f, IntLumi[electronMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2), ((TH1F*)alberto[f]->Get(histoName))->Integral(), ((TH1F*)alberto[f]->Get(histoName))->GetEntries(), IntLumi[electronMask-1]*xs[f]/((TH1F*)alberto[f]->Get("Entries"))->GetBinContent(2)*((TH1F*)alberto[f]->Get(histoName))->Integral(), hNbtaggedJets_vvjets[c][i][j]->Integral());
                }
              }
            }
          }
        }
      }
    }
    /*
    for(int c=0; c<NbOfChannels; c++) { 
      for (int i=0; i<NbOfJetBins; i++) {
        if (hNbtaggedJets_vvjets[c][i][0]!=NULL)
          printf("Before closing : hNbtaggedJets_vvjets[%d][%d] = %lf\n", c, i, hNbtaggedJets_vvjets[c][i][0]->Integral());
        if (hNbtaggedJets_qcd[c][i][0]!=NULL)
          printf("Before closing : hNbtaggedJets_qcd   [%d][%d] = %lf\n", c, i, hNbtaggedJets_qcd[c][i][0]->Integral());
      }
    }
    */
    for (int f=0; f<NbOfFiles; f++) {
      alberto[f]->Close();
    }
  }


  /*
  for(int c=0; c<NbOfChannels; c++) { 
    for (int i=0; i<NbOfJetBins; i++) {
      if (hNbtaggedJets_vvjets[c][i][0]!=NULL)
        printf("After closing : hNbtaggedJets_vvjets[%d][%d] = %lf\n", c, i, hNbtaggedJets_vvjets[c][i][0]->Integral());
      if (hNbtaggedJets_qcd[c][i][0]!=NULL)
        printf("After closing : hNbtaggedJets_qcd   [%d][%d] = %lf\n", c, i, hNbtaggedJets_qcd[c][i][0]->Integral());
    }
  }
*/
    


  vector<int>::iterator Idx;
  vector<int> FixedVarIdx;
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
  if (! reloadHisto) {
    if(verbose){
      cout<<"Number of events (entries in chain) for :"<<endl;
      cout<<" - ttjets = "<<ttjets_Muon->GetEntries()<<"   ";
      cout<<" - ttjets = "<<ttjets_Electron->GetEntries()<<endl;
      cout<<" -  wjets = "<< wjets_Muon->GetEntries()<<"   ";
      cout<<" -  wjets = "<< wjets_Electron->GetEntries()<<endl;
      cout<<" -  zjets = "<< zjets_Muon->GetEntries()<<"    ";
      cout<<" -  zjets = "<< zjets_Electron->GetEntries()<<endl;
      cout<<" - stjets = "<<stjets_Muon->GetEntries()<<"    ";
      cout<<" - stjets = "<<stjets_Electron->GetEntries()<<endl;
        //    cout<<" - vbjets ="<<vbjets->GetEntries()<<endl;
    }
  }
  
  Long64_t NbOfEvtsToProcess = 0;
  
    // VERIFIER : ICI JE PRENDS COMME SI C_ETAIT MESURE SUR DONNEES
  Double_t** number_events_mistagLight = new Double_t*[NbOfChannels];
  Double_t** number_jets_mistagLight = new Double_t*[NbOfChannels];
  Double_t** sum_events_mistagLight = new Double_t*[NbOfChannels];
  Double_t** sum_jets_mistagLight = new Double_t*[NbOfChannels];
  Double_t** number_events_SFLight = new Double_t*[NbOfChannels];
  Double_t** number_jets_SFLight = new Double_t*[NbOfChannels];
  Double_t** sum_events_SFLight = new Double_t*[NbOfChannels];
  Double_t** sum_jets_SFLight = new Double_t*[NbOfChannels];
  
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
  
  for(int c=0; c<NbOfChannels; c++) {
    number_events_mistagLight[c] = new Double_t[NbOfWP];
    number_jets_mistagLight[c] = new Double_t[NbOfWP];
    sum_events_mistagLight[c] = new Double_t[NbOfWP];
    sum_jets_mistagLight[c] = new Double_t[NbOfWP];
    number_events_SFLight[c] = new Double_t[NbOfWP];
    number_jets_SFLight[c] = new Double_t[NbOfWP];
    sum_events_SFLight[c] = new Double_t[NbOfWP];
    sum_jets_SFLight[c] = new Double_t[NbOfWP];
    
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
    for (int i=0; i<NbOfWP; i++) {
      number_events_mistagLight[c][i] = 0.;
      number_jets_mistagLight[c][i] = 0.;
      sum_events_mistagLight[c][i] = 0.;
      sum_jets_mistagLight[c][i] = 0.;
      number_events_SFLight[c][i] = 0.;
      number_jets_SFLight[c][i] = 0.;
      sum_events_SFLight[c][i] = 0.;
      sum_jets_SFLight[c][i] = 0.;
      
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
  
  float **Nttlike = new float*[NbOfChannels];
  double *NttlikeTot = new double[NbOfChannels];
  float **Nvlike = new float*[NbOfChannels];
  double *NvlikeTot = new double[NbOfChannels];
  for (int c=0; c<NbOfChannels; c++) {
    Nttlike[c] = new float[NbOfJetBins];
    NttlikeTot[c] = 0.;
    Nvlike[c] = new float[NbOfJetBins];
    NvlikeTot[c] = 0.;
  }
  
  printf("Ready to go through all the events (if not reloading histograms)\n");
  
    //cout<<" - data : Muon = "<<
    //        anaMu_data->Show(0);
    //<<endl;
  
  if (isOnData && !reloadHisto) {
    //    NbOfEvtsToProcess = entrylist_data->GetN();
    NbOfEvtsToProcess = 0;
    Long64_t NbOfEvtsToProcess_Muon = dataTree_Muon->GetEntries();
    Long64_t NbOfEvtsToProcess_Electron = dataTree_Electron->GetEntries();
//    if (semiMuon) {
    NbOfEvtsToProcess += NbOfEvtsToProcess_Muon;
//    }
//    if (semiElectron){
    NbOfEvtsToProcess += NbOfEvtsToProcess_Electron;
//    }
    
    
    if(verbose) cout<<endl<<"Processing data Tree (Nb of events = "<<NbOfEvtsToProcess<<")"<<endl;
    std::string currentFileName = "";
    for (Long64_t el = 0; el<NbOfEvtsToProcess; el++) {
      if(el%2000 == 0) cout<<"-- Processing the "<<el<<"th event ("<<100*(el)/(NbOfEvtsToProcess)<<"%)"<<flush<<"\r";
        //      dataTree->GetEntry(entrylist_data->GetEntry(el));
      
        // Begin Selection part
      bool triggerMu=false;
      bool triggerEl=false;
      int NbOfFiredChannels = 0;
      if (((semiMuon && ((el>=0) && (el<NbOfEvtsToProcess_Muon)))/* || (!areDataChannelTreatedExclusive)*/)) {
        size_t k=0;
        /*
         while ( (!triggerMu) && (k < anaMu_data->trigNames->size()) ) {
         if ( 160431<=anaMu_data->run && anaMu_data->run<=165633 && (anaMu_data->trigNames->at(k)).find("HLT_Mu17_TriCentralJet30_v")==0 )
         triggerMu = true; 
         else if ( 165970<=anaMu_data->run && anaMu_data->run<=173198 && (anaMu_data->trigNames->at(k)).find("HLT_IsoMu17_TriCentralJet30_v")==0 )
         triggerMu = true;
         else if ( 173236<=anaMu_data->run && anaMu_data->run<=178380 && (anaMu_data->trigNames->at(k)).find("HLT_IsoMu17_eta2p1_TriCentralJet30_v")==0 )
         triggerMu = true;
         else if ( 178420<=anaMu_data->run && anaMu_data->run<=180252 && (anaMu_data->trigNames->at(k)).find("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v")==0 )
         triggerMu = true;
         k++;
         }
         */
        string token1("HLT_IsoMu17_eta2p1_TriCentralJet30");
        string token2("HLT_Mu17_TriCentralJet30");
        string token3("HLT_IsoMu17_TriCentralJet30");
        string token4("HLT_IsoMu17_eta2p1_TriCentralPFJet30");
        
        size_t found1,found2,found3,found4;
        
          //found = trigger.find(token);
//        try {
          anaMu_data->GetEntry(el);
//        }
//        catch (...) {
//          printf("\nLast event : %l\n\n", el);
//        }
/*
          if (std::string(dataTree_Muon->GetCurrentFile()->GetName()) !=currentFileName) {
            TObjArray *arr = dataTree_Muon->GetListOfFiles();
            printf("\n");
            for (int kk=0; kk<arr->GetSize(); kk++) {
              if (arr->At(kk) != NULL)
                printf("%s , ", ((TChainElement*) arr->At(kk))->GetTitle());
            }
            printf("\n\nChanging file from   %s   to   %s \n\n", currentFileName.c_str(), dataTree_Muon->GetCurrentFile()->GetName());
            currentFileName = dataTree_Muon->GetCurrentFile()->GetName();
          }
*/        
         
          //        dataTree_Muon->GetEntry(el);
          //                string trign = (*( anaMu_data->trigNames))[0];
        
          //        for (std::vector<std::string>::iterator it = anaMu_data->trigNames->begin(); it!= anaMu_data->trigNames->end(); it++)
          //          {
          //            printf(" %s ,", it->c_str());
          //          }
        std::vector<std::string> trigNames = *(anaMu_data->trigNames);
        Int_t run = anaMu_data->run;

        while ( (!triggerMu) && (k < trigNames.size()) ) {
          found1 =  trigNames.at(k).find(token1);
          found2 =  trigNames.at(k).find(token2);
          found3 =  trigNames.at(k).find(token3);
          found4 =  trigNames.at(k).find(token4);
          
          if( run >= 160431 && run <= 165633 && found2 != string::npos )
            triggerMu = true;
          else if( run >= 165970 && run <= 173198 && found3 != string::npos )
            triggerMu = true;
          else if( run >= 173236 && run <= 178380 && found1 != string::npos )
            triggerMu = true;
          else if( run >= 178420 && run <= 180252 && found4 != string::npos )
            triggerMu = true;
          k++;
        }
        
        k=0;
      }
      if (((semiElectron && ((el>=NbOfEvtsToProcess_Muon) && (el<NbOfEvtsToProcess_Muon+NbOfEvtsToProcess_Electron)))/* || (!areDataChannelTreatedExclusive)*/)) {
        size_t k=0;
          //        printf("\"Dont forget to implement me\", Mr semiElectron DATA\n");
        string token1("HLT_Ele25_CaloIdVT_TrkIdT_TriCentralJet30");
        string token2("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30");
        string token3("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30");
        string token4("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30");
          //        string tok("HLT_Ele25");
        
        size_t  found1,found2,found3,found4;
          //found = trigger.find(token);
//        try {
          anaEl_data->GetEntry(el-NbOfEvtsToProcess_Muon);
          /*      }
        catch (...) {
          printf("\nLast event : %l\n\n", el);
        }
        
        if (std::string(dataTree_Electron->GetCurrentFile()->GetName()) !=currentFileName) {
          printf("Changing file from   %s   to   %s \n", currentFileName.c_str(), dataTree_Electron->GetCurrentFile()->GetName());
          currentFileName = dataTree_Electron->GetCurrentFile()->GetName();
        }
          */
        std::vector<std::string> trigNames = *(anaEl_data->trigNames);
        Int_t run = anaEl_data->run;

        while ( (!triggerEl) && (k < trigNames.size()) ) {
          found1 = trigNames.at(k).find(token1);
          found2 = trigNames.at(k).find(token2);
          found3 = trigNames.at(k).find(token3);
          found4 = trigNames.at(k).find(token4);
          if( run >= 160403 /*160431*/ && run <= 163869 && found2 != string::npos )
            triggerEl = true;
          else if( run >= 165088 && run <= 165633 && found1 != string::npos )
            triggerEl = true;
          else if( run >= 165970 && run <= 178380 && found3 != string::npos )
            triggerEl = true;
          else if( run >= 178420 && run <= 180252 && found4 != string::npos )
            triggerEl = true;
          k++;
        }
        
        k=0;
      }
        //      if (triggerMu || triggerEl)
        //        printf("trigger Mu:%d El:%d\n", (int) triggerMu, (int) triggerEl);
      
        // End (trigger) Selection part
      if (semiMuon && triggerMu && anaMu_data->CutMu_iso0p1_3jets(-1)==0) {
        NbOfFiredChannels++;
        Int_t njt = anaMu_data->Njt_pf;
        Int_t njt_idx = njt-MinNbOfJets;
        if (njt_idx>=NbOfJetBins && isLastBinExclusive) {
          continue;
        } else if (njt_idx>=NbOfJetBins) {
          njt_idx=NbOfJetBins-1;
        }
        if (njt<MinNbOfJets) {
          continue;
        }
        std::vector<double> discri = *(anaMu_data->vjtDiscri_pf);
        std::vector<double> pt = *(anaMu_data->vjtPt_pf);
        std::vector<double> eta = *(anaMu_data->vjtEta_pf);
        for (int c=0; c<NbOfChannels; c++){
          if ( ((c+1) & muonMask) == 0)
            continue;
          for (int wp=0; wp<NbOfWP; wp++) {
            NbtaggedJets = 0;
            for(UInt_t i = 0; i<discri.size(); ++i)
            {
              if(pt.at(i)<30.) { printf("Jet looser for tt_jets\n"); }
              
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
              
              if(discri[i]>btagCuts[wp]) {
                NbtaggedJets++;
              }
            }
            if (njt_idx>=0) {
                //          printf("fill [%d,%d,%d] %d \t", c, njt_idx,wp, NbtaggedJets);
              hNbtaggedJets[c][njt_idx][wp]->Fill(NbtaggedJets);
            }
          }
          histoNJ_data[c]->Fill(njt);
        }
      }
      if (semiElectron && triggerEl && anaEl_data->Cut_iso0p1_3jets(-1)==0) {
        NbOfFiredChannels++;
        Int_t njt = anaEl_data->Njt_pf;
        Int_t njt_idx = njt-MinNbOfJets;
        if (njt_idx>=NbOfJetBins && isLastBinExclusive) {
          continue;
        } else if (njt_idx>=NbOfJetBins) {
          njt_idx=NbOfJetBins-1;
        }
        if (njt<MinNbOfJets) {
          continue;
        }
        std::vector<double> discri = *(anaEl_data->vjtDiscri_pf);
        std::vector<double> pt = *(anaEl_data->vjtPt_pf);
        std::vector<double> eta = *(anaEl_data->vjtEta_pf);
        for (int c=0; c<NbOfChannels; c++){
          if ( ((c+1) & electronMask) == 0)
            continue;
          for (int wp=0; wp<NbOfWP; wp++) {
            NbtaggedJets = 0;
            for(UInt_t i = 0; i<discri.size(); ++i)
            {
              if(pt.at(i)<30.) { printf("Jet looser for tt_jets\n"); }
              
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
              
              if(discri[i]>btagCuts[wp]) {
                NbtaggedJets++;
              }
            }
            if (njt_idx>=0) {
                //          printf("fill [%d,%d,%d] %d \t", c, njt_idx,wp, NbtaggedJets);
              hNbtaggedJets[c][njt_idx][wp]->Fill(NbtaggedJets);
            }
          }
          histoNJ_data[c]->Fill(njt);
        }
      }
      if (NbOfFiredChannels > 1) {
        printf("Several channels fired for the same entry !!!!\n");
      }
    }
    cout << endl;
  }
  
  
    // For tt+jets events 
  if (!reloadHisto) {
      //    NbOfEvtsToProcess = entrylist_ttjets->GetN();
    Long64_t n_Mu = ttjets_Muon->GetEntriesFast();
    Long64_t n_El = ttjets_Electron->GetEntriesFast();
    NbOfEvtsToProcess = n_Mu + n_El; //CHECK
                                     //             printf("bl0\n");
    if(verbose) cout<<endl<<"Processing tt+jets Tree (Nb of events = "<<NbOfEvtsToProcess<<")"<<endl;
      //             printf("bl1\n");
    for (Long64_t el = 0; el<NbOfEvtsToProcess; el++) {
        //double lumiWeight = LumiWeights.ITweight( (int) truenuminter_ttjets );
      double lumiWeight = 1.;
        //            printf("bl2\n");
      if (el%10000 == 0) cout<<"-- Processing the "<<el<<"th event ("<<100*(el)/(NbOfEvtsToProcess)<<"%)"<<flush<<"\r";
        //      ttjets->LoadTree(entrylist_ttjets->GetEntry(el));
        //        ttjets->GetEntry(entrylist_ttjets->GetEntry(el));
        //   printf("bl3\n");
      /* 
       if (semiMuon) {
       anaMu_ttjets->GetEntry(el);
       }
       if (semiElectron) {
       anaEl_ttjets->GetEntry(el);
       }
       */      
        // Begin Selection part
      bool triggerMu=false;
      bool triggerEl=false;
      std::string trigPat = "" ;
        //           printf("bli\n");
      if (semiMuon && el>=0 && el<n_Mu) {
          //          anaMu_ttjets->Init(ttjets);
        anaMu_ttjets->GetEntry(el);
        size_t k=0;
          //          printf("bliaaa %lu    %lu      %s\n", (unsigned long) anaMu_ttjets, (unsigned long) anaMu_ttjets->trigNames, ttjets->GetCurrentFile()->GetName());
        std::vector<std::string> trigNames = *(anaMu_ttjets->trigNames);
        while ( (!triggerMu) && (k < trigNames.size()) ) {  
            //  printf("blili\n");
          if ( ( trigNames.at(k).find("HLT_IsoMu17_eta2p1_TriCentralJet30_v")==0)||(trigNames[k].find("HLT_Mu17_TriCentralJet30_v")==0) ) {
            triggerMu = true;
          }
          k++;
        }
        k=0;
        while ((k < trigNames.size()) ) {  
          if ( ((trigNames.at(k)).find("HLT_IsoMu17_eta2p1_TriCentralJet30_v")==0)||((trigNames.at(k)).find("HLT_Mu17_TriCentralJet30_v")==0) ) {
            trigPat = trigPat + trigNames.at(k) + " && ";
          }
          k++;
        }
      }
        //    printf("blo\n");
      
        //        FillTriggerPattern(histoTrig_ttjets, trigPat, 1.);
        //        printf("bly    mu : %d      ele : %d\t\t", (int) triggerMu, (int) triggerEl);
        // End (trigger) Selection part
      if (semiMuon && triggerMu && anaMu_ttjets->CutMu_iso0p1_3jets(-1)==0) {
          //          printf("\n");
        
        lumiWeight = numInteractions->GetBinContent(numInteractions->FindBin(anaMu_ttjets->truenuminter));
        Int_t njt = anaMu_ttjets->Njt_pf;
        Int_t njt_idx = njt - MinNbOfJets;
        if (njt_idx>=NbOfJetBins && isLastBinExclusive) {
          continue;
        } else if (njt_idx >= NbOfJetBins) {
          njt_idx=NbOfJetBins-1;
        }
        
        if (njt < MinNbOfJets) {
          continue;
        }
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
      
      if (semiElectron && el>=n_Mu && el<n_Mu+n_El) {
          //          anaEl_ttjets->Init(ttjets);
        anaEl_ttjets->GetEntry(el-n_Mu);
          //          printf("\"Dont forget to implement me\", Mr semiElectron TTJETS\n");
        size_t k=0;
        std::vector<std::string> trigNames = *(anaEl_ttjets->trigNames);
        while ( (!triggerEl) && (k < trigNames.size()) ) {  
            //            printf("blu\n");
          if ( ((trigNames.at(k)).find("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30")==0)||((trigNames.at(k)).find("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30")==0) ) {
            triggerEl = true;
          }
          k++;
        }
        k=0;
        while ((k < trigNames.size()) ) {  
          if ( ((trigNames.at(k)).find("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30")==0)||((trigNames.at(k)).find("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30")==0) ) {
            trigPat = trigPat + trigNames.at(k) + " && ";
          }
          k++;
        }
      }
      if (semiElectron && triggerEl && anaEl_ttjets->Cut_iso0p1_3jets(-1)==0) {
          //           printf("\n");
        lumiWeight = numInteractions->GetBinContent(numInteractions->FindBin(anaEl_ttjets->truenuminter));
        Int_t njt = anaEl_ttjets->Njt_pf;
        Int_t njt_idx = njt - MinNbOfJets;
        if (njt_idx>=NbOfJetBins && isLastBinExclusive) {
          continue;
        } else if (njt_idx >= NbOfJetBins) {
          njt_idx=NbOfJetBins-1;
        }
        
        if (njt < MinNbOfJets) {
          continue;
        }
        Double_t MCweight = anaEl_ttjets->MCweight;
        std::vector<double> discri = *(anaEl_ttjets->vjtDiscri_pf);
        std::vector<double> flav = *(anaEl_ttjets->vjtFlav_pf);
        std::vector<double> pt = *(anaEl_ttjets->vjtPt_pf);
        std::vector<double> eta = *(anaEl_ttjets->vjtEta_pf);
        double selWeight = (IntLumi[electronMask-1] * preselEff_ttjets_NUM[electronMask-1]) / (preselEff_ttjets_DEN[electronMask-1]*n_El);
        for (int c=0; c<NbOfChannels; c++){
          if ( ((c+1) &electronMask) == 0)
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
      
      
      /*
       } else {
       printf("Problem with selection on number of jets at the TCut/EntryList level\n");
       }*/
    }
    
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
        /*
         printf("----->>>>>> Mistag rate in tt+jets sample (with NormFactor) WP %d : %lf mean %lf up %lf down %lf = %lf / %lf\n", j,
         nEffUdsc_pass[c][j] / nEffUdsc_tot[c][j],
         WilsonScoreIntervalMean(nEffUdsc_pass[c][j]*NormFact_ttjets*IntLumi[c],nEffUdsc_tot[c][j]*NormFact_ttjets*IntLumi[c]),
         WilsonScoreIntervalHigh(nEffUdsc_pass[c][j]*NormFact_ttjets*IntLumi[c],nEffUdsc_tot[c][j]*NormFact_ttjets*IntLumi[c]),
         WilsonScoreIntervalLow(nEffUdsc_pass[c][j]*NormFact_ttjets*IntLumi[c],nEffUdsc_tot[c][j]*NormFact_ttjets*IntLumi[c]),
         nEffUdsc_pass[c][j]*NormFact_ttjets*IntLumi[c], nEffUdsc_tot[c][j]*NormFact_ttjets*IntLumi[c]);
         */
        
        
        for (int i=0; i<NbOfJetBins; i++) {
          hNbtaggedJets_ttjets[c][i][j]->Scale(XS_ttjets    ); // Selection, on skimmed dataset, to correct the number of entries of the tree
        }
      }
    }
    cout<< endl;
  }
  
  
    // For w+jets events 
  if (!reloadHisto) {
      //      NbOfEvtsToProcess = entrylist_wjets->GetN();
    Long64_t n_Mu = wjets_Muon->GetEntriesFast();
    Long64_t n_El = wjets_Electron->GetEntriesFast();
    NbOfEvtsToProcess = n_Mu + n_El;
    if(verbose) cout<<endl<<"Processing W+jets Tree (Nb of events = "<<NbOfEvtsToProcess<<")"<<endl;
    for (int c=0; c<NbOfChannels; c++) {
      for (int j=0; j<NbOfWP; j++) {
        nEffUds_pass_tmp[c][j] = 0.;
        nEffUds_tot_tmp[c][j] = 0.;
      }
    }
    for (Long64_t el = 0; el<NbOfEvtsToProcess;el++) {
      double lumiWeight = 1.;
      if (el%10000 == 0) cout<<"-- Processing the "<<el<<"th event ("<<100*(el)/(NbOfEvtsToProcess)<<"%)"<<flush<<"\r";
        //      ttjets->LoadTree(entrylist_wjets->GetEntry(el));
        //        ttjets->GetEntry(entrylist_wjets->GetEntry(el));
      
        // Begin Selection part
      bool triggerMu=false;
      bool triggerEl=false;
      std::string trigPat = "" ;
      
        //FillTriggerPattern(histoTrig_wjets, trigPat, 1.);
      
        // End (trigger) Selection part
      if (semiMuon && el>=0 && el<n_Mu) {
          //          anaMu_wjets->Init(wjets);
        anaMu_wjets->GetEntry(el);
        size_t k=0;
        while ( (!triggerMu) && (k < anaMu_wjets->trigNames->size()) ) {  
          if ( ((anaMu_wjets->trigNames->at(k)).find("HLT_IsoMu17_eta2p1_TriCentralJet30_v")==0)||((anaMu_wjets->trigNames->at(k)).find("HLT_Mu17_TriCentralJet30_v")==0) ) {
            triggerMu = true;
          }
          k++;
        }
        k=0;
        while ((k < anaMu_wjets->trigNames->size()) ) {  
          if ( ((anaMu_wjets->trigNames->at(k)).find("HLT_IsoMu17_eta2p1_TriCentralJet30_v")==0)||((anaMu_wjets->trigNames->at(k)).find("HLT_Mu17_TriCentralJet30_v")==0) ) {
            trigPat = trigPat + anaMu_wjets->trigNames->at(k) + " && ";
          }
          k++;
        }
      }
      if (semiMuon && triggerMu && anaMu_wjets->CutMu_iso0p1_3jets(-1)==0) {
        lumiWeight = numInteractions->GetBinContent(numInteractions->FindBin(anaMu_wjets->truenuminter));
        Int_t njt_idx = anaMu_wjets->Njt_pf-MinNbOfJets;
        if (njt_idx>=NbOfJetBins && isLastBinExclusive) {
          continue;
        } else if (njt_idx>=NbOfJetBins) {
          njt_idx=NbOfJetBins-1;
        }
        if (anaMu_wjets->Njt_pf<MinNbOfJets) {
          continue;
        }
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
              /*
               if ((*vjtDiscri_pf_wjets)[i]>btagCut && std::fabs((*vjtFlav_pf_wjets)[i])!=5.) {
               nEffUds_pass+=MCweight_wjets*lumiWeight*NormFact_wjets;
               }
               if (std::fabs((*vjtFlav_pf_wjets)[i])!=5.) {
               nEffUds_tot+=MCweight_wjets*lumiWeight*NormFact_wjets;
               }
               */  
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
          histoNJ_wjets[c]->Fill(anaMu_wjets->Njt_pf, anaMu_wjets->MCweight*lumiWeight * selWeight);
        }
      }
      
      if (semiElectron && el>=n_Mu && el<n_Mu+n_El) {
          //          anaEl_wjets->Init(wjets);
        anaEl_wjets->GetEntry(el-n_Mu);
          //          printf("\"Dont forget to implement me\", Mr semiElectron WJETS\n");
        size_t k=0;
        while ( (!triggerEl) && (k < anaEl_wjets->trigNames->size()) ) {  
          if ( ((anaEl_wjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30")==0)||((anaEl_wjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30")==0) ) {
            triggerEl = true;
          }
          k++;
        }
        k=0;
        while ((k < anaEl_wjets->trigNames->size()) ) {  
          if ( ((anaEl_wjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30")==0)||((anaEl_wjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30")==0) ) {
            trigPat = trigPat + anaEl_wjets->trigNames->at(k) + " && ";
          }
          k++;
        }
      }
      if (semiElectron && triggerEl && anaEl_wjets->Cut_iso0p1_3jets(-1)==0) {
        lumiWeight = numInteractions->GetBinContent(numInteractions->FindBin(anaEl_wjets->truenuminter));
        Int_t njt_idx = anaEl_wjets->Njt_pf-MinNbOfJets;
        if (njt_idx>=NbOfJetBins && isLastBinExclusive) {
          continue;
        } else if (njt_idx>=NbOfJetBins) {
          njt_idx=NbOfJetBins-1;
        }
        if (anaEl_wjets->Njt_pf<MinNbOfJets) {
          continue;
        }
        double selWeight= (IntLumi[electronMask-1] * preselEff_wjets_NUM[electronMask-1]) / (preselEff_wjets_DEN[electronMask-1]*n_El);
        for (int c=0; c<NbOfChannels; c++){
          if ( ((c+1) &electronMask) == 0)
            continue;
          for (int wp=0; wp<NbOfWP; wp++) {
            nB_Flavoured_Jets = 0;
            nC_Flavoured_Jets = 0;
            nLightAndGluon_Flavoured_Jets = 0;
            
            NbtaggedJets = 0;
            for(UInt_t i = 0; i<anaEl_wjets->vjtDiscri_pf->size(); ++i)
            {
              if(anaEl_wjets->vjtPt_pf->at(i)<30.) { printf("Jet looser for w_jets\n"); }
              if((*(anaEl_wjets->vjtDiscri_pf))[i]>btagCuts[wp]) NbtaggedJets++;
              if (std::fabs((*(anaEl_wjets->vjtFlav_pf))[i])==5.) {
                nB_Flavoured_Jets++;
              } else if (std::fabs((*(anaEl_wjets->vjtFlav_pf))[i])==4.) {
                nC_Flavoured_Jets++;
              } else {
                nLightAndGluon_Flavoured_Jets++;
              }
              /*
               if ((*vjtDiscri_pf_wjets)[i]>btagCut && std::fabs((*vjtFlav_pf_wjets)[i])!=5.) {
               nEffUds_pass+=MCweight_wjets*lumiWeight*NormFact_wjets;
               }
               if (std::fabs((*vjtFlav_pf_wjets)[i])!=5.) {
               nEffUds_tot+=MCweight_wjets*lumiWeight*NormFact_wjets;
               }
               */  
            }
            if (nB_Flavoured_Jets == 0) {
              for(UInt_t i = 0; i<anaEl_wjets->vjtDiscri_pf->size(); ++i)
              {
                if ((*(anaEl_wjets->vjtDiscri_pf))[i]>btagCuts[wp]) {
                  nEffUds_pass_tmp[c][wp]+=anaEl_wjets->MCweight*lumiWeight;
                }
                nEffUds_tot_tmp[c][wp]+=anaEl_wjets->MCweight*lumiWeight;
              }
            } else {
              for(UInt_t i = 0; i<anaEl_wjets->vjtDiscri_pf->size(); ++i)
              {
                if (std::fabs((*(anaEl_wjets->vjtFlav_pf))[i])==5.) {
                  if((*(anaEl_wjets->vjtDiscri_pf))[i]>btagCuts[wp]) {
                    nEffb_w_pass[c][wp] += anaEl_wjets->MCweight*lumiWeight;
                  }
                  nEffb_w_tot[c][wp] += anaEl_wjets->MCweight*lumiWeight;
                }
              }
            }
            hNbtaggedJets_wjets[c][njt_idx][wp]->Fill(NbtaggedJets, anaEl_wjets->MCweight*lumiWeight * selWeight);
            if (nB_Flavoured_Jets>0) {
              hWbb_In_WJets[c][njt_idx][wp]->Fill(NbtaggedJets, anaEl_wjets->MCweight*lumiWeight * selWeight);
            } else if (nC_Flavoured_Jets>0) {
              hWcx_In_WJets[c][njt_idx][wp]->Fill(NbtaggedJets, anaEl_wjets->MCweight*lumiWeight * selWeight);
            } else {
              hWLFg_In_WJets[c][njt_idx][wp]->Fill(NbtaggedJets, anaEl_wjets->MCweight*lumiWeight * selWeight);
            }
          }
          histoNJ_wjets[c]->Fill(anaEl_wjets->Njt_pf, anaEl_wjets->MCweight*lumiWeight * selWeight);
        }
      }
    }
    
    
    
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
        /*
         printf("----->>>>>> Mistag rate in W+jets sample (with NormFact) WP %d : %lf mean %lf up %lf down %lf = %lf / %lf\n",
         j, nEffUds_pass_tmp[c][j] / nEffUds_tot_tmp[c][j],
         WilsonScoreIntervalMean(nEffUds_pass_tmp[c][j]*NormFact_wjets*IntLumi[c],nEffUds_tot_tmp[c][j]*NormFact_wjets*IntLumi[c]),
         WilsonScoreIntervalHigh(nEffUds_pass_tmp[c][j]*NormFact_wjets*IntLumi[c],nEffUds_tot_tmp[c][j]*NormFact_wjets*IntLumi[c]),
         WilsonScoreIntervalLow(nEffUds_pass_tmp[c][j]*NormFact_wjets*IntLumi[c],nEffUds_tot_tmp[c][j]*NormFact_wjets*IntLumi[c]),
         nEffUds_pass_tmp[c][j]*NormFact_wjets*IntLumi[c], nEffUds_tot_tmp[c][j]*NormFact_wjets*IntLumi[c]);
         */
        /*
         nEffUds_pass_tmp[c][j] *= NormFact_wjets;
         nEffUds_tot_tmp[c][j]  *= NormFact_wjets;
         */
        
        nEffUds_pass[c][j] += nEffUds_pass_tmp[c][j];
        nEffUds_tot[c][j]  += nEffUds_tot_tmp[c][j];
        
        for (int i=0; i<NbOfJetBins; i++) {
          hNbtaggedJets_wjets[c][i][j]->Scale(XS_wjets /*NormFact_wjets*IntLumi[c]*/ ); // Selection, on skimmed dataset, to correct the number of entries of the tree
        }
      }
    }
    cout<< endl;
  }
  
  
  
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
          hWbb_In_WJets[c][i][j]->Scale(XS_wjets /* NormFact_wjets*IntLumi[c]*/*normChangeWbb);
          hWcx_In_WJets[c][i][j]->Scale(XS_wjets /* NormFact_wjets*IntLumi[c]*/);
          hWLFg_In_WJets[c][i][j]->Scale(XS_wjets /* NormFact_wjets*IntLumi[c]*/);
        }
      }
    }
  }
  

  
    // For z+jets events 
  if (!reloadHisto) {
      //      NbOfEvtsToProcess = entrylist_zjets->GetN();
    Long64_t n_Mu = zjets_Muon->GetEntriesFast();
    Long64_t n_El = zjets_Electron->GetEntriesFast();
    NbOfEvtsToProcess = n_Mu + n_El;
    if(verbose) cout<<endl<<"Processing Z+jets Tree (Nb of events = "<<NbOfEvtsToProcess<<")"<<endl;
    for (int c=0; c<NbOfChannels; c++) {
      for (int j=0; j<NbOfWP; j++) {
        nEffUds_pass_tmp[c][j] = 0.;
        nEffUds_tot_tmp[c][j] = 0.;
      }
    }
    for (Long64_t el = 0; el<NbOfEvtsToProcess;el++) {
      double lumiWeight = 1.;
      if (el%10000 == 0) cout<<"-- Processing the "<<el<<"th event ("<<100*(el)/(NbOfEvtsToProcess)<<"%)"<<flush<<"\r";
        //      ttjets->LoadTree(entrylist_zjets->GetEntry(el));
        //        ttjets->GetEntry(entrylist_zjets->GetEntry(el));
      
        // Begin Selection part
      bool triggerMu=false;
      bool triggerEl=false;
      std::string trigPat = "" ;
      if (semiMuon && el>=0 && el<n_Mu) {
          //          anaMu_zjets->Init(zjets);
        anaMu_zjets->GetEntry(el);
        size_t k=0;
        while ( (!triggerMu) && (k < anaMu_zjets->trigNames->size()) ) {  
          if ( ((anaMu_zjets->trigNames->at(k)).find("HLT_IsoMu17_eta2p1_TriCentralJet30_v")==0)||((anaMu_zjets->trigNames->at(k)).find("HLT_Mu17_TriCentralJet30_v")==0) ) {
            triggerMu = true;
          }
          k++;
        }
        k=0;
        while ((k < anaMu_zjets->trigNames->size()) ) {  
          if ( ((anaMu_zjets->trigNames->at(k)).find("HLT_IsoMu17_eta2p1_TriCentralJet30_v")==0)||((anaMu_zjets->trigNames->at(k)).find("HLT_Mu17_TriCentralJet30_v")==0) ) {
            trigPat = trigPat + anaMu_zjets->trigNames->at(k) + " && ";
          }
          k++;
        }
      }
      
        //    FillTriggerPattern(histoTrig_zjets, trigPat, 1.);
        // End (trigger) Selection part
      
      
        // End Selection part
      if (semiMuon && triggerMu && anaMu_zjets->CutMu_iso0p1_3jets(-1)==0) {
        lumiWeight = numInteractions->GetBinContent(numInteractions->FindBin(anaMu_zjets->truenuminter));
        Int_t njt_idx = anaMu_zjets->Njt_pf-MinNbOfJets;
        if (njt_idx>=NbOfJetBins && isLastBinExclusive) {
          continue;
        } else if (njt_idx>=NbOfJetBins) {
          njt_idx=NbOfJetBins-1;
        }
        if (anaMu_zjets->Njt_pf<MinNbOfJets) {
          continue;
        }
        double selWeight= (IntLumi[muonMask-1] * preselEff_zjets_NUM[muonMask-1]) / (preselEff_zjets_DEN[muonMask-1]*n_Mu);
        for (int c=0; c<NbOfChannels; c++) {
          if ( ((c+1) &muonMask) == 0)
            continue;
          for (int wp=0; wp<NbOfWP; wp++) {
            NbtaggedJets = 0;
            nB_Flavoured_Jets = 0;
            for(UInt_t i = 0; i<anaMu_zjets->vjtDiscri_pf->size(); ++i) {
              if(anaMu_zjets->vjtPt_pf->at(i)<30.) { printf("Jet looser for z_jets\n"); }
              if((*(anaMu_zjets->vjtDiscri_pf))[i]>btagCuts[wp]) NbtaggedJets++;
              if (std::fabs((*(anaMu_zjets->vjtFlav_pf))[i])==5.) {
                nB_Flavoured_Jets++;
              }
              /*  
               if ((*vjtDiscri_pf_zjets)[i]>btagCut && std::fabs((*vjtFlav_pf_zjets)[i])!=5.) {
               nEffUds_pass+=MCweight_zjets*lumiWeight*NormFact_zjets;
               }
               if (std::fabs((*vjtFlav_pf_zjets)[i])!=5.) {
               nEffUds_tot+=MCweight_zjets*lumiWeight*NormFact_zjets;
               }
               */
            }
            if (nB_Flavoured_Jets == 0) {
              for(UInt_t i = 0; i<anaMu_zjets->vjtDiscri_pf->size(); ++i) {
                if ((*(anaMu_zjets->vjtDiscri_pf))[i]>btagCuts[wp]) {
                  nEffUds_pass_tmp[c][wp]+=anaMu_zjets->MCweight*lumiWeight;
                }
                nEffUds_tot_tmp[c][wp]+=anaMu_zjets->MCweight*lumiWeight;
              }
            }
            hNbtaggedJets_zjets[c][njt_idx][wp]->Fill(NbtaggedJets,anaMu_zjets->MCweight * lumiWeight * selWeight);
          }
          histoNJ_zjets[c]->Fill(anaMu_zjets->Njt_pf, anaMu_zjets->MCweight*lumiWeight * selWeight);
        }
      }
      
      if (semiElectron && el>=n_Mu && el<n_Mu+n_El) {
          //           anaEl_zjets->Init(zjets);
        anaEl_zjets->GetEntry(el-n_Mu);
          //          printf("\"Dont forget to implement me\", Mr semiElectron ZJETS\n");
        size_t k=0;
        while ( (!triggerEl) && (k < anaEl_zjets->trigNames->size()) ) {  
          if ( ((anaEl_zjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30")==0)||((anaEl_zjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30")==0) ) {
            triggerEl = true;
          }
          k++;
        }
        k=0;
        while ((k < anaEl_zjets->trigNames->size()) ) {  
          if ( ((anaEl_zjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30")==0)||((anaEl_zjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30")==0) ) {
            trigPat = trigPat + anaEl_zjets->trigNames->at(k) + " && ";
          }
          k++;
        }
      }
      if (semiElectron && triggerEl && anaEl_zjets->Cut_iso0p1_3jets(-1)==0) {
        lumiWeight = numInteractions->GetBinContent(numInteractions->FindBin(anaEl_zjets->truenuminter));
        Int_t njt_idx = anaEl_zjets->Njt_pf-MinNbOfJets;
        if (njt_idx>=NbOfJetBins && isLastBinExclusive) {
          continue;
        } else if (njt_idx>=NbOfJetBins) {
          njt_idx=NbOfJetBins-1;
        }
        if (anaEl_zjets->Njt_pf<MinNbOfJets) {
          continue;
        }
        double selWeight= (IntLumi[electronMask-1] * preselEff_zjets_NUM[electronMask-1]) / (preselEff_zjets_DEN[electronMask-1]*n_El);
        for (int c=0; c<NbOfChannels; c++) {
          if ( ((c+1) &electronMask) == 0)
            continue;
          for (int wp=0; wp<NbOfWP; wp++) {
            NbtaggedJets = 0;
            nB_Flavoured_Jets = 0;
            for(UInt_t i = 0; i<anaEl_zjets->vjtDiscri_pf->size(); ++i) {
              if(anaEl_zjets->vjtPt_pf->at(i)<30.) { printf("Jet looser for z_jets\n"); }
              if((*(anaEl_zjets->vjtDiscri_pf))[i]>btagCuts[wp]) NbtaggedJets++;
              if (std::fabs((*(anaEl_zjets->vjtFlav_pf))[i])==5.) {
                nB_Flavoured_Jets++;
              }
              /*  
               if ((*vjtDiscri_pf_zjets)[i]>btagCut && std::fabs((*vjtFlav_pf_zjets)[i])!=5.) {
               nEffUds_pass+=MCweight_zjets*lumiWeight*NormFact_zjets;
               }
               if (std::fabs((*vjtFlav_pf_zjets)[i])!=5.) {
               nEffUds_tot+=MCweight_zjets*lumiWeight*NormFact_zjets;
               }
               */
            }
            if (nB_Flavoured_Jets == 0) {
              for(UInt_t i = 0; i<anaEl_zjets->vjtDiscri_pf->size(); ++i) {
                if ((*(anaEl_zjets->vjtDiscri_pf))[i]>btagCuts[wp]) {
                  nEffUds_pass_tmp[c][wp]+=anaEl_zjets->MCweight*lumiWeight;
                }
                nEffUds_tot_tmp[c][wp]+=anaEl_zjets->MCweight*lumiWeight;
              }
            }
            hNbtaggedJets_zjets[c][njt_idx][wp]->Fill(NbtaggedJets,anaEl_zjets->MCweight * lumiWeight * selWeight);
          }
          histoNJ_zjets[c]->Fill(anaEl_zjets->Njt_pf, anaEl_zjets->MCweight*lumiWeight * selWeight);
        }
      }
      
      
    }
    
    
    for (int c=0; c<NbOfChannels; c++) {
      printf("\nFor channel #%d\n", c);
      for (int j=0; j<NbOfWP; j++) {
        printf("\n\n----->>>>>> Mistag rate in Z+jets sample WP %d : %lf mean %lf up %lf down %lf = %lf / %lf\n", j,
               nEffUds_pass_tmp[c][j] / nEffUds_tot_tmp[c][j],
               WilsonScoreIntervalMean(nEffUds_pass_tmp[c][j],nEffUds_tot_tmp[c][j]),
               WilsonScoreIntervalHigh(nEffUds_pass_tmp[c][j],nEffUds_tot_tmp[c][j]),
               WilsonScoreIntervalLow(nEffUds_pass_tmp[c][j],nEffUds_tot_tmp[c][j]),
               nEffUds_pass_tmp[j], nEffUds_tot_tmp[c][j]);
        /*
         printf("\n\n----->>>>>> Mistag rate in Z+jets sample (with NormFact) WP %d : %lf mean %lf up %lf down %lf = %lf / %lf\n", j,
         nEffUds_pass_tmp[c][j] / nEffUds_tot_tmp[c][j],
         WilsonScoreIntervalMean(nEffUds_pass_tmp[c][j]*NormFact_zjets*IntLumi[c],nEffUds_tot_tmp[c][j]*NormFact_zjets*IntLumi[c]),
         WilsonScoreIntervalHigh(nEffUds_pass_tmp[c][j]*NormFact_zjets*IntLumi[c],nEffUds_tot_tmp[c][j]*NormFact_zjets*IntLumi[c]),
         WilsonScoreIntervalLow(nEffUds_pass_tmp[c][j]*NormFact_zjets*IntLumi[c],nEffUds_tot_tmp[c][j]*NormFact_zjets*IntLumi[c]),
         nEffUds_pass_tmp[c][j], nEffUds_tot_tmp[c][j]);
         */
        /*       nEffUds_pass_tmp[c][j] *= NormFact_zjets;
         nEffUds_tot_tmp[c][j]  *= NormFact_zjets;
         */
        nEffUds_pass[c][j] += nEffUds_pass_tmp[c][j];
        nEffUds_tot[c][j]  += nEffUds_tot_tmp[c][j];
      }
      for (int i=0; i<NbOfJetBins; i++) {
        for (int wp=0; wp<NbOfWP; wp++) {
          hNbtaggedJets_zjets[c][i][wp]->Scale(XS_zjets/*NormFact_zjets*IntLumi[c]*/    ); // Selection, on skimmed dataset, to correct the number of entries of the tree);
        }
      }
    }
    cout <<endl;
  }
  
    // For st+jets events 
  if (!reloadHisto) {
    Long64_t n_Mu = wjets_Muon->GetEntriesFast();
    Long64_t n_El = wjets_Electron->GetEntriesFast();
    NbOfEvtsToProcess = n_Mu + n_El;
    int* n = new int[NbOfChannels];
    double* Norm = new double[NbOfChannels];
    for (int c=0; c<NbOfChannels; c++) {
      n[c] = 0;
      Norm[c] = 0.;
    }
    if(verbose) cout<<endl<<"Processing st+jets Tree (Nb of events = "<<NbOfEvtsToProcess<<")"<<endl;
    for (Long64_t el = 0; el<NbOfEvtsToProcess;el++) {
      double lumiWeight = 1.;
      if (el%10000 == 0) cout<<"-- Processing the "<<el<<"th event ("<<100*(el)/(NbOfEvtsToProcess)<<"%)"<<flush<<"\r";
        //      ttjets->LoadTree(entrylist_stjets->GetEntry(el));
        //        ttjets->GetEntry(entrylist_stjets->GetEntry(el));
      
        // Begin Selection part
      bool triggerMu=false;
      bool triggerEl=false;
      std::string trigPat = "" ;
      if (semiMuon && el>=0 && el<n_Mu) {
          //           anaMu_stjets->Init(stjets);
        anaMu_stjets->GetEntry(el);
        size_t k=0;
        while ( (!triggerMu) && (k < anaMu_stjets->trigNames->size()) ) {  
          if ( ((anaMu_stjets->trigNames->at(k)).find("HLT_IsoMu17_eta2p1_TriCentralJet30_v")==0)||((anaMu_stjets->trigNames->at(k)).find("HLT_Mu17_TriCentralJet30_v")==0) ) {
            triggerMu = true;
          }
          k++;
        }
        k=0;
        while ((k < anaMu_stjets->trigNames->size()) ) {  
          if ( ((anaMu_stjets->trigNames->at(k)).find("HLT_IsoMu17_eta2p1_TriCentralJet30_v")==0)||((anaMu_stjets->trigNames->at(k)).find("HLT_Mu17_TriCentralJet30_v")==0) ) {
            trigPat = trigPat + anaMu_stjets->trigNames->at(k) + " && ";
          }
          k++;
        }
      }
      
        //FillTriggerPattern(histoTrig_stjets, trigPat, 1.);
      
        // End (trigger) Selection part
      if (semiMuon && triggerMu && anaMu_stjets->CutMu_iso0p1_3jets(-1)==0) {
        lumiWeight = numInteractions->GetBinContent(numInteractions->FindBin(anaMu_stjets->truenuminter));
        Int_t njt_idx = anaMu_stjets->Njt_pf-MinNbOfJets;
        if (njt_idx>=NbOfJetBins && isLastBinExclusive) {
          continue;
        } else if (njt_idx>=NbOfJetBins) {
          njt_idx=NbOfJetBins-1;
        }
        if (anaMu_stjets->Njt_pf<MinNbOfJets) {
          continue;
        }
          //          double selWeight= (IntLumi[muonMask-1] * preselEff_stjets_NUM[muonMask-1]) / (preselEff_stjets_DEN[muonMask-1]*n_Mu);
        for (int c=0; c<NbOfChannels; c++){
          if ( ((c+1) & muonMask) == 0)
            continue;
          for (int wp=0; wp<NbOfWP; wp++) {
            NbtaggedJets = 0;
            for(UInt_t i = 0; i<anaMu_stjets->vjtDiscri_pf->size(); ++i) {
              if(anaMu_stjets->vjtPt_pf->at(i)<30.) { printf("Jet looser for st_jets\n"); }
              if((*(anaMu_stjets->vjtDiscri_pf))[i]>btagCuts[wp]) NbtaggedJets++;
            }
            
            std::string filename = stjets_Muon->GetCurrentFile()->GetName();
              //      printf("--> %s\n", filename.c_str());
            if (skimmedTree) {
              hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaMu_stjets->MCweight * NormFact_stjets[c]);
            } else {
              if (filename.find("_s_")!=std::string::npos && filename.find("Tbar_")==std::string::npos) {
                  //                  printf("s : %d, %d, %lf, %lf, %lf, %lf\n", njt_idx, NbtaggedJets,  lumiWeight,  anaMu_stjets->MCweight, XS_stjets_s, preselEff_stjets_s_DEN);
                Norm[c]+=IntLumi[muonMask-1]*XS_stjets_s/ preselEff_stjets_s_DEN[muonMask-1] ;
                n[c]++;
                hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaMu_stjets->MCweight * (IntLumi[muonMask-1] * XS_stjets_s)/(preselEff_stjets_s_DEN[muonMask-1]));
              } else if (filename.find("_s_")!=std::string::npos && filename.find("Tbar_")!=std::string::npos) {
                Norm[c]+=IntLumi[muonMask-1]*XS_stjets_s_bar/ preselEff_stjets_s_DEN_bar[muonMask-1] ;
                n[c]++;
                hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaMu_stjets->MCweight * (IntLumi[muonMask-1] * XS_stjets_s_bar)/(preselEff_stjets_s_DEN_bar[muonMask-1]));
              } else if (filename.find("_t_")!=std::string::npos && filename.find("Tbar_")==std::string::npos) {
                  //                  printf("t\n");
                Norm[c]+=IntLumi[muonMask-1]*XS_stjets_t/preselEff_stjets_t_DEN[muonMask-1] ;
                n[c]++;
                hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaMu_stjets->MCweight * (IntLumi[muonMask-1] * XS_stjets_t)/(preselEff_stjets_t_DEN[muonMask-1]));
              } else if (filename.find("_t_")!=std::string::npos && filename.find("Tbar_")!=std::string::npos) {
                  //                  printf("t\n");
                Norm[c]+=IntLumi[muonMask-1]*XS_stjets_t_bar/preselEff_stjets_t_DEN_bar[muonMask-1] ;
                n[c]++;
                hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaMu_stjets->MCweight * (IntLumi[muonMask-1] * XS_stjets_t_bar)/(preselEff_stjets_t_DEN_bar[muonMask-1]));
              } else  if (filename.find("_twDR_")!=std::string::npos) {
                  //                  printf("twDR\n");
                Norm[c]+=IntLumi[muonMask-1]*XS_stjets_twDR/preselEff_stjets_twDR_DEN[muonMask-1] ;
                n[c]++;
                hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaMu_stjets->MCweight * (IntLumi[muonMask-1] * XS_stjets_twDR)/(preselEff_stjets_twDR_DEN[muonMask-1]));
              } else {
                printf("You have a problem with Single-top\n");
              }
              if (wp==0) {
                NormFact_stjets[c] = Norm[c] / n[c];
              }
            }
          }
          histoNJ_stjets[c]->Fill(anaMu_stjets->Njt_pf, anaMu_stjets->MCweight*lumiWeight);
        }
      }
      
      if (semiElectron && el>=n_Mu && el<n_Mu+n_El) {
          //          anaEl_stjets->Init(stjets);
        anaEl_stjets->GetEntry(el-n_Mu);
          //          printf("\"Dont forget to implement me\", Mr semiElectron STJETS\n");
        size_t k=0;
        while ( (!triggerEl) && (k < anaEl_stjets->trigNames->size()) ) {  
          if ( ((anaEl_stjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30")==0)||((anaEl_stjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30")==0) ) {
            triggerEl = true;
          }
          k++;
        }
        k=0;
        while ((k < anaEl_stjets->trigNames->size()) ) {  
          if ( ((anaEl_stjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_TrkIdT_CentralTriJet30")==0)||((anaEl_stjets->trigNames->at(k)).find("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralJet30")==0) ) {
            trigPat = trigPat + anaEl_stjets->trigNames->at(k) + " && ";
          }
          k++;
        }
      }
      if (semiElectron && triggerEl && anaEl_stjets->Cut_iso0p1_3jets(-1)==0) {
        lumiWeight = numInteractions->GetBinContent(numInteractions->FindBin(anaEl_stjets->truenuminter));
        Int_t njt_idx = anaEl_stjets->Njt_pf-MinNbOfJets;
        if (njt_idx>=NbOfJetBins && isLastBinExclusive) {
          continue;
        } else if (njt_idx>=NbOfJetBins) {
          njt_idx=NbOfJetBins-1;
        }
        if (anaEl_stjets->Njt_pf<MinNbOfJets) {
          continue;
        }
          //          double selWeight= (IntLumi[electronMask-1] * preselEff_stjets_NUM[electronMask-1]) / (preselEff_stjets_DEN[electronMask-1]*n_El);
        for (int c=0; c<NbOfChannels; c++){
          if ( ((c+1) & electronMask) == 0)
            continue;
          for (int wp=0; wp<NbOfWP; wp++) {
            NbtaggedJets = 0;
            for(UInt_t i = 0; i<anaEl_stjets->vjtDiscri_pf->size(); ++i) {
              if(anaEl_stjets->vjtPt_pf->at(i)<30.) { printf("Jet looser for st_jets\n"); }
              if((*(anaEl_stjets->vjtDiscri_pf))[i]>btagCuts[wp]) NbtaggedJets++;
            }
            
            std::string filename = stjets_Electron->GetCurrentFile()->GetName();
              //      printf("--> %s\n", filename.c_str());
            if (skimmedTree) {
              hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaEl_stjets->MCweight * NormFact_stjets[c]);
            } else {
              if (filename.find("_s_")!=std::string::npos && filename.find("Tbar_")==std::string::npos) {
                Norm[c]+=IntLumi[electronMask-1]*XS_stjets_s/preselEff_stjets_s_DEN[electronMask-1] ;
                n[c]++;
                hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaEl_stjets->MCweight * (IntLumi[electronMask-1] * XS_stjets_s)/(preselEff_stjets_s_DEN[electronMask-1]));
              } else if (filename.find("_s_")!=std::string::npos && filename.find("Tbar_")!=std::string::npos) {
                Norm[c]+=IntLumi[electronMask-1]*XS_stjets_s_bar/preselEff_stjets_s_DEN_bar[electronMask-1] ;
                n[c]++;
                hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaEl_stjets->MCweight * (IntLumi[electronMask-1] * XS_stjets_s_bar)/(preselEff_stjets_s_DEN_bar[electronMask-1]));
              } else if (filename.find("_t_")!=std::string::npos && filename.find("Tbar_")==std::string::npos) {
                Norm[c]+=IntLumi[electronMask-1]*XS_stjets_t/preselEff_stjets_t_DEN[electronMask-1] ;
                n[c]++;
                hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaEl_stjets->MCweight * (IntLumi[electronMask-1] * XS_stjets_t)/(preselEff_stjets_t_DEN[electronMask-1]));
              } else if (filename.find("_t_")!=std::string::npos && filename.find("Tbar_")!=std::string::npos) {
                Norm[c]+=IntLumi[electronMask-1]*XS_stjets_t_bar/preselEff_stjets_t_DEN_bar[electronMask-1] ;
                n[c]++;
                hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaEl_stjets->MCweight * (IntLumi[electronMask-1] * XS_stjets_t_bar)/(preselEff_stjets_t_DEN_bar[electronMask-1]));
              } else  if (filename.find("_twDR_")!=std::string::npos) {
                Norm[c]+=IntLumi[electronMask-1]*XS_stjets_twDR/preselEff_stjets_twDR_DEN[electronMask-1] ;
                n[c]++;
                hNbtaggedJets_stjets[c][njt_idx][wp]->Fill(NbtaggedJets, lumiWeight * anaEl_stjets->MCweight * (IntLumi[electronMask-1] * XS_stjets_twDR)/(preselEff_stjets_twDR_DEN[electronMask-1]));
              } else {
                printf("You have a problem with Single-top\n");
              }
              if (wp==0) {
                NormFact_stjets[c] = Norm[c] / n[c];
              }
            }
          }
          histoNJ_stjets[c]->Fill(anaEl_stjets->Njt_pf, anaEl_stjets->MCweight*lumiWeight);
        }
      }
    }
    
    printf("\n");
    for (int c=0; c<NbOfChannels; c++){
      for (int wp=0; wp<NbOfWP; wp++) {
        for (int i=0; i<NbOfJetBins; i++) {
          printf("Int Lumi : %lf  : integral wp%d jetBin%d = %lf \n", IntLumi[c], wp, i, hNbtaggedJets_stjets[c][i][wp]->Integral());
          hNbtaggedJets_stjets[c][i][wp]->Scale(/*IntLumi[c]*/ 1.); // Selection, on skimmed dataset, to correct the number of entries of the tree
        }
      }
    }
    cout<< endl;
  }
  
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
  
  for (int c=0; c<NbOfChannels; c++) {
    for (int i=0; i<NbOfJetBins; i++) {
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
            hNbtaggedJets_vlike[c][i][wp]->Fill(hNbtaggedJets_qcd[c][i][wp]->GetBinLowEdge(k), hNbtaggedJets_qcd[c][i][wp]->GetBinContent(k)/*0.*/);
          }
        }
      }
        // Set true numbers of tt-like and v-like events
      Nttlike[c][i] = hNbtaggedJets_ttlike[c][i][0]->Integral();
      NttlikeTot[c] += Nttlike[c][i];
      Nvlike[c][i]  = hNbtaggedJets_vlike[c][i][0]->Integral();
      NvlikeTot[c] += Nvlike[c][i];
    }
  }
  
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
  for (int c=0; c<NbOfChannels; c++) {
    for (int wp=0; wp<NbOfWP; wp++) {
      if (reloadHisto) {
        init_Euds[c][wp] =  eudsMC_histo->GetBinContent(c+1, wp+1);
        init_Eudsc[c][wp] = eudscMC_histo->GetBinContent(c+1, wp+1);
      } else {
        init_Euds[c][wp] = nEffUds_pass[c][wp]/nEffUds_tot[c][wp] ;
        init_Eudsc[c][wp]= nEffUdsc_pass[c][wp]/nEffUdsc_tot[c][wp] ;
        if (isOnData) {
// Applying scale factors for light from BTV-11-004
//init_Euds[c][wp] *= sum_jets_SFLight[c][wp] / number_jets_SFLight[c][wp];
          printf("Applying Scale Factor [%d][%d] :   sum : %lf \t\t n_evts : %lf\t\t SF : %lf \n", c, wp, sum_events_SFLight[c][wp], number_events_SFLight[c][wp], sum_events_SFLight[c][wp] / number_events_SFLight[c][wp]);
          init_Euds[c][wp] *= sum_events_SFLight[c][wp] / number_events_SFLight[c][wp];
          //          init_Eudsc[c][wp] *= sum_events_SFLight[c][wp] / number_events_SFLight[c][wp];
          
/* // Mistag rate from BTV-11-004
//           init_Euds[c][wp] = sum_jets_mistagLight[c][wp] / number_jets_mistagLight[c][wp];
init_Euds[c][wp] = sum_events_mistagLight[c][wp] / number_events_mistagLight[c][wp];
*/
        }
      }

      Euds_const_mean[c][wp] = init_Euds[c][wp];
      Euds_const_error[c][wp] = WilsonScoreIntervalHigh(nEffUds_pass[c][wp]*IntLumi[c],nEffUds_tot[c][wp]*IntLumi[c]) - init_Euds[c][wp];
      if ( init_Euds[c][wp] - WilsonScoreIntervalLow(nEffUds_pass[c][wp]*IntLumi[c],nEffUds_tot[c][wp]*IntLumi[c]) > Euds_const_error[c][wp] )
        Euds_const_error[c][wp] = init_Euds[c][wp] - WilsonScoreIntervalLow(nEffUds_pass[c][wp]*IntLumi[c],nEffUds_tot[c][wp]*IntLumi[c]);
      //        Euds_const_error[wp] *= 20.;
      //        if (Euds_const_error[wp]/Euds_const_mean[wp] >= 0.15)  
      Euds_const_error[c][wp] = Euds_const_mean[c][wp] * 0.15;
      
      
      //Euds_const_error = WilsonScoreIntervalMean(nEffUds_pass,nEffUds_tot);
      Eudsc_const_mean[c][wp] = init_Eudsc[c][wp] ;
      Eudsc_const_error[c][wp] = WilsonScoreIntervalHigh(nEffUdsc_pass[c][wp]*IntLumi[c],nEffUdsc_tot[c][wp]*IntLumi[c]) - init_Eudsc[c][wp];
      if ( init_Eudsc[c][wp] - WilsonScoreIntervalLow(nEffUdsc_pass[c][wp],nEffUdsc_tot[c][wp]) > Eudsc_const_error[c][wp] )
        Eudsc_const_error[c][wp] = init_Eudsc[c][wp] - WilsonScoreIntervalLow(nEffUdsc_pass[c][wp],nEffUdsc_tot[c][wp]);
      //Eudsc_const_error = WilsonScoreIntervalMean(nEffUdsc_pass,nEffUdsc_tot);
      
      Eudsc_const_error[c][wp] = Eudsc_const_mean[c][wp] * .15 ;
      init_Eb[c][wp] = Eb_const_mean[c][wp] ;
      
        /*
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
         */
        
      ebTTMC_histo->SetBinContent(c+1, wp+1, nEffb_tt_pass[c][wp]/nEffb_tt_tot[c][wp]);
      ebVMC_histo->SetBinContent(c+1, wp+1, nEffb_w_pass[c][wp]/nEffb_w_tot[c][wp]);
      eudsMC_histo->SetBinContent(c+1, wp+1, init_Euds[c][wp]);
      eudscMC_histo->SetBinContent(c+1, wp+1, init_Eudsc[c][wp]);
      
    }    
  }
  
  
  for (int c=0; c<NbOfChannels; c++) {
    if (Wbb_from_WJets) {
      for (int i=0; i<NbOfJetBins; i++) {
        for (int wp=0; wp<NbOfWP; wp++) {
          for (int k =0; k<10; k++) {
            hNbtaggedJets_vbjets[c][i][wp]->SetBinContent(k+1, hWbb_In_WJets[c][i][wp]->GetBinContent(hWbb_In_WJets[c][i][wp]->FindBin(k)));
          }
        }
      }
    } else {
      int Nbins = 10 ;
      Int_t Min = 0 ;
      Int_t Max = 10 ;
      TFile *f = new TFile(WbbShapeRootfile.c_str(), "READ");
      if (!f->IsZombie()) {
        char name[100];
        for (UInt_t j=0; j<12; j++) {
          if ((j!=1)&&(j!=2)&&(j!=5)) {
            continue;
          }
          for (int i=0; i<NbOfJetBins; i++) {
            for (int wp=0; wp<NbOfWP; wp++) {
              sprintf(name, (std::string()+"HistoryFlavor--FlavorHistoryPath_%d"+channelSuffix[c]+jetSuffix[i]+wpSuffix[wp]).c_str(), j);
              TH1I *h = (TH1I*) f->Get(name);
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
    for (int c=0; c<NbOfChannels; c++) {
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
      for (int i=0; i<NbOfJetBins; i++) {
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
          Nqcd_mean[c][i] = hNbtaggedJets_qcd[c][i][0]->Integral()/*0.*/;
        }
        Nqcd_uncert[c][i] = 0.30 ;
        NmultijetsTot[c] += Nqcd_mean[c][i];
      }
    }
    
  
    // Sum over all datasets
  if (!isOnData && !reloadHisto) {
    for (int c=0; c<NbOfChannels; c++) {
      for (int wp=0; wp<NbOfWP; wp++) {
        for (int i=0; i<NbOfJetBins; i++) {
          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_ttlike[c][i][wp]);
          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_vlike[c][i][wp]);
          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_stjets[c][i][wp]);
          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_vbjets[c][i][wp]);


//          if (hNbtaggedJets_qcd[c][i][wp]!=NULL) {
//            for(int k=0; k<hNbtaggedJets_qcd[c][i][wp]->GetNbinsX(); k++) {
////                          hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_qcd[c][i][wp]);
//              hNbtaggedJets[c][i][wp]->Fill(hNbtaggedJets_qcd[c][i][wp]->GetBinLowEdge(k), hNbtaggedJets_qcd[c][i][wp]->GetBinContent(k)/*0.*/);
//            }
//          }

          if (hNbtaggedJets_vvjets[c][i][wp]!=NULL) {
            for(int k=0; k<hNbtaggedJets_vvjets[c][i][wp]->GetNbinsX(); k++) {
//                            hNbtaggedJets[c][i][wp]->Add(hNbtaggedJets_vvjets[c][i][wp]);
              hNbtaggedJets[c][i][wp]->Fill(hNbtaggedJets_vvjets[c][i][wp]->GetBinLowEdge(k), hNbtaggedJets_vvjets[c][i][wp]->GetBinContent(k)/*0.*/);
            }
          }
        }
      }
    }
  }
  
    // Total nb of events (rescaled to the required int. lumi.)
  Int_t** nExpected = new Int_t*[NbOfChannels];
  Int_t* nTot = new Int_t[NbOfChannels];
  for (int c=0; c<NbOfChannels; c++) {
    nExpected[c] = new Int_t[NbOfJetBins];
    nTot[c] = 0;
    for (int i=0; i<NbOfJetBins; i++) {
      nExpected[c][i] = 0;
      for(int j=0; j<NbOfWP; j++) {
        nExpected[c][i] += (int) hNbtaggedJets[c][i][0]->Integral();
      }
      nTot[c] += nExpected[c][i];
    }
    
    for(int i=0; i<NbOfJetBins; i++) {
      eXbq__[c][i]->Scale(1./eXbq__[c][i]->Integral());
    }
  }
  
  for (int c=0; c<NbOfChannels; c++) {
    for (int wp=0; wp<NbOfWP; wp++) {
      for (int i=0; i<NbOfJetBins; i++) {
        for (int k=1; k<=hNbtaggedJets_ttlike[c][i][wp]->GetNbinsX(); k++)
          hNbtaggedJets_ttlike[c][i][wp]->SetBinError(k, TMath::Sqrt(hNbtaggedJets_ttlike[c][i][wp]->GetBinContent(k)));
        for (int k=1; k<=hNbtaggedJets_vlike[c][i][wp]->GetNbinsX(); k++)
          hNbtaggedJets_vlike[c][i][wp]->SetBinError(k, TMath::Sqrt(hNbtaggedJets_vlike[c][i][wp]->GetBinContent(k)));
        for (int k=1; k<=hNbtaggedJets_stjets[c][i][wp]->GetNbinsX(); k++)
          hNbtaggedJets_stjets[c][i][wp]->SetBinError(k, TMath::Sqrt(hNbtaggedJets_stjets[c][i][wp]->GetBinContent(k)));
        for (int k=1; k<=hNbtaggedJets_vbjets[c][i][wp]->GetNbinsX(); k++)
          hNbtaggedJets_vbjets[c][i][wp]->SetBinError(k, TMath::Sqrt(hNbtaggedJets_vbjets[c][i][wp]->GetBinContent(k)));
        for (int k=1; k<=hNbtaggedJets[c][i][wp]->GetNbinsX(); k++)
          hNbtaggedJets[c][i][wp]->SetBinError(k, TMath::Sqrt(hNbtaggedJets[c][i][wp]->GetBinContent(k)));
      }
    }
  }
  
  
  if(verbose){
    cout<<"***************************************"<<endl;
    cout<<"***************************************"<<endl;
    for (int c=0; c<NbOfChannels; c++) {
      cout<<"- Number of entries / NormFactors "<<endl;
      cout<<"-- for tt-like processes : "<<endl;
      for (int i=0; i<NbOfJetBins; i++) {
        cout<<"--- tt+jets (jet bin "<< i << ") = "<<hNbtaggedJets_ttjets[c][i][0]->Integral()<<"/ "<<NormFact_ttjets[c]<<endl;
      }
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NttlikeTot[c]<<endl;
      cout<<"---------------------------------"<<endl;
      cout<<"-- for v-like processes : "<<endl;
      for (int i=0; i<NbOfJetBins; i++) {
        cout<<"--- W_LFg+jets (jet bin "<< i << ") = "<<hWLFg_In_WJets[c][i][0]->Integral()<<"/ "<<NormFact_wjets[c]<<endl;
        cout<<"--- W_cx+jets (jet bin "<< i << ") = "<<hWcx_In_WJets[c][i][0]->Integral()<<"/ "<<NormFact_wjets[c]<<endl;
      }
      for (int i=0; i<NbOfJetBins; i++) {
        cout<<"--- Z+jets (jet bin "<< i << ") = "<<hNbtaggedJets_zjets[c][i][0]->Integral()<<"/ "<<NormFact_zjets[c]<<endl;
      }
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NvlikeTot[c]<<endl;
      cout<<"---------------------------------"<<endl;
      cout<<"-- for background processes : "<<endl;
      for (int i=0; i<NbOfJetBins; i++) {
        cout<<"--- st+jets (jet bin "<< i << ") = "<<hNbtaggedJets_stjets[c][i][0]->Integral()<<"/ "<<NormFact_stjets[c]<<endl;
      }
        //    cout<<"--- st+jets = "<<hNbtaggedJets_stjets->GetEntries()<<"/ "<<NormFact_stjets<<endl;
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NstjetsTot[c]<<endl;
      for (int i=0; i<NbOfJetBins; i++) {
        cout<<"--- vb+jets (jet bin "<< i << ") = "<<hNbtaggedJets_vbjets[c][i][0]->Integral()<<endl;//"/ "<<NormFact_vbjets<<endl;
      }
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NvbjetsTot[c]<<endl;
      for (int i=0; i<NbOfJetBins; i++) {
        if (hNbtaggedJets_qcd[c][i][0]!=NULL) {
          cout<<"--- multijets (jet bin "<< i << ") = "<<hNbtaggedJets_qcd[c][i][0]->Integral()<<"/ "<<NormFact_qcd[c]<<endl;
        }
      }
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NmultijetsTot[c]<<endl;
      for (int i=0; i<NbOfJetBins; i++) {
        if (hNbtaggedJets_vvjets[c][i][0]!=NULL) {
          cout<<"--- VV+jets (jet bin "<< i << ") = "<<hNbtaggedJets_vvjets[c][i][0]->Integral()<<"/ "<<NormFact_vvjets[c]<<endl;
        }
      }
      cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<NvvjetsTot[c]<<endl;
      cout<<"---------------------------------"<<endl;
      if (isOnData) {
        cout<<"-- for Data : "<<endl;
        double sum = 0.;
        double sumMC = 0.;
        for (int i=0; i<NbOfJetBins; i++) {
          sum += hNbtaggedJets[c][i][0]->Integral();
          double tmpMCsum = hNbtaggedJets_ttjets[c][i][0]->Integral() + hWLFg_In_WJets[c][i][0]->Integral()+hWcx_In_WJets[c][i][0]->Integral() + hNbtaggedJets_zjets[c][i][0]->Integral() + hNbtaggedJets_stjets[c][i][0]->Integral() + hNbtaggedJets_vbjets[c][i][0]->Integral();
          sumMC += tmpMCsum;
          cout<<"--- total (jet bin "<< i << ") = "<<hNbtaggedJets[c][i][0]->Integral() << "         and for MC : " << tmpMCsum <<endl;
        }
        cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<< sum /*nTot[c]/NbOfWP*/ << "          and for MC : " << sumMC << endl;
      } else {
        cout<<"-- for all processes : "<<endl;
        for (int i=0; i<NbOfJetBins; i++) {
          cout<<"--- total (jet bin "<< i << ") = "<<hNbtaggedJets[c][i][0]->Integral()<<endl;
        }
        cout<<"--- total nb of events (for "<<IntLumi[c]<<" /pb) = "<<nTot[c]/NbOfWP<<endl;
      }
      cout<<"***************************************"<<endl;
    }
    cout<<"***************************************"<<endl;
  }
  
  if (isOnData) {
    printf("Saving total b-tagged jets distribution plot\n");
    TFile *f = new TFile((std::string()+postfix+"_TotalHisto_hNbtaggedJets.root").c_str(), "recreate");
    f->cd();
    for (int c=0; c<NbOfChannels; c++) {
      for (int i=0; i<NbOfJetBins; i++) {
        for (int j=0; j<NbOfWP; j++) {
          hNbtaggedJets[c][i][j]->Write();
        }
      }
      histoNJ_data[c]->Write();
    }
    f->Close();
  }
  
  
    // C r e a t e   o b s e r v a b l e
    // ---------------------------------
  printf("Efficiencies\n");
  
  
  RooRealVar ***eb = new RooRealVar**[NbOfChannels];
  RooRealVar ***eudsc = new RooRealVar**[NbOfChannels];
  RooRealVar ***euds = new RooRealVar**[NbOfChannels];
  
  for (int c=0; c<NbOfChannels; c++) {
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    eb[c] = new RooRealVar*[NbOfWP];
    eudsc[c] = new RooRealVar*[NbOfWP];
    euds[c] = new RooRealVar*[NbOfWP];
    
    for (int j=0; j<NbOfWP; j++) {
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

  for (int c=0; c<NbOfChannels; c++) {
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
    n0bjets_st[c] = new RooConstVar**[NbOfJetBins]; n1bjets_st[c] = new RooConstVar**[NbOfJetBins]; n2bjets_st[c] = new RooConstVar**[NbOfJetBins];
    n3bjets_st[c] = new RooConstVar**[NbOfJetBins]; n4bjets_st[c] = new RooConstVar**[NbOfJetBins]; n5bjets_st[c] = new RooConstVar**[NbOfJetBins];
    n0bjets_vb[c] = new RooConstVar**[NbOfJetBins]; n1bjets_vb[c] = new RooConstVar**[NbOfJetBins]; n2bjets_vb[c] = new RooConstVar**[NbOfJetBins];
    n3bjets_vb[c] = new RooConstVar**[NbOfJetBins]; n4bjets_vb[c] = new RooConstVar**[NbOfJetBins]; n5bjets_vb[c] = new RooConstVar**[NbOfJetBins];
    n0bjets_qcd[c] = new RooConstVar**[NbOfJetBins]; n1bjets_qcd[c] = new RooConstVar**[NbOfJetBins]; n2bjets_qcd[c] = new RooConstVar**[NbOfJetBins];
    n3bjets_qcd[c] = new RooConstVar**[NbOfJetBins]; n4bjets_qcd[c] = new RooConstVar**[NbOfJetBins]; n5bjets_qcd[c] = new RooConstVar**[NbOfJetBins];
    n0bjets_vv[c] = new RooConstVar**[NbOfJetBins]; n1bjets_vv[c] = new RooConstVar**[NbOfJetBins]; n2bjets_vv[c] = new RooConstVar**[NbOfJetBins];
    n3bjets_vv[c] = new RooConstVar**[NbOfJetBins]; n4bjets_vv[c] = new RooConstVar**[NbOfJetBins]; n5bjets_vv[c] = new RooConstVar**[NbOfJetBins];

    printf("Other variables\n");
    
    for (int i=0; i<NbOfJetBins; i++) {
      
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
      for (int j=0; j<NbOfWP; j++) {
        n0bjets_st[c][i][j] = new RooConstVar((std::string()+"n0bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{0 b-jet}", hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(0.)));
        n1bjets_st[c][i][j] = new RooConstVar((std::string()+"n1bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{1 b-jet}", hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(1.)));
        n2bjets_st[c][i][j] = new RooConstVar((std::string()+"n2bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{2 b-jets}",hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(2.)));
        n3bjets_st[c][i][j] = new RooConstVar((std::string()+"n3bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{3 b-jets}",hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(3.)));
        n4bjets_st[c][i][j] = new RooConstVar((std::string()+"n4bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{4 b-jets}",hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(4.)));
        n5bjets_st[c][i][j] = new RooConstVar((std::string()+"n5bjets_st"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_stjets[c][i][j]->GetBinContent(hNbtaggedJets_stjets[c][i][j]->FindBin(5.)));

        n0bjets_vb[c][i][j] = new RooConstVar((std::string()+"n0bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{0 b-jet}", hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(0.)));
        n1bjets_vb[c][i][j] = new RooConstVar((std::string()+"n1bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{1 b-jet}", hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(1.)));
        n2bjets_vb[c][i][j] = new RooConstVar((std::string()+"n2bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{2 b-jets}",hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(2.)));
        n3bjets_vb[c][i][j] = new RooConstVar((std::string()+"n3bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{3 b-jets}",hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(3.)));
        n4bjets_vb[c][i][j] = new RooConstVar((std::string()+"n4bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{4 b-jets}",hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(4.)));
        n5bjets_vb[c][i][j] = new RooConstVar((std::string()+"n5bjets_vb"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_vbjets[c][i][j]->GetBinContent(hNbtaggedJets_vbjets[c][i][j]->FindBin(5.)));

        n0bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n0bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{0 b-jet}", hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(0.)));
        n1bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n1bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{1 b-jet}", hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(1.)));
        n2bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n2bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{2 b-jets}",hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(2.)));
        n3bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n3bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{3 b-jets}",hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(3.)));
        n4bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n4bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{4 b-jets}",hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(4.)));
        n5bjets_qcd[c][i][j] = new RooConstVar((std::string()+"n5bjets_qcd"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_qcd[c][i][j]->GetBinContent(hNbtaggedJets_qcd[c][i][j]->FindBin(5.)));

        n0bjets_vv[c][i][j] = new RooConstVar((std::string()+"n0bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{0 b-jet}", hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(0.)));
        n1bjets_vv[c][i][j] = new RooConstVar((std::string()+"n1bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{1 b-jet}", hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(1.)));
        n2bjets_vv[c][i][j] = new RooConstVar((std::string()+"n2bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{2 b-jets}",hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(2.)));
        n3bjets_vv[c][i][j] = new RooConstVar((std::string()+"n3bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{3 b-jets}",hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(3.)));
        n4bjets_vv[c][i][j] = new RooConstVar((std::string()+"n4bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{4 b-jets}",hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(4.)));
        n5bjets_vv[c][i][j] = new RooConstVar((std::string()+"n5bjets_vv"+channelSuffix[c]+jetSuffix[i]+wpSuffix[j]).c_str(),"n^{5 b-jets}",hNbtaggedJets_vvjets[c][i][j]->GetBinContent(hNbtaggedJets_vvjets[c][i][j]->FindBin(5.)));
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
  for (int i=0; i<NbOfJetBins; i++) {
    jetMultCat.defineType(jetSuffix[i], i);
  }
  
  RooCategory wpCat("WP","b-tagging Working Point");
  for (int i=0; i<NbOfWP; i++) {
    wpCat.defineType(wpSuffix[i], i);
  }
  
  RooCategory wpJetCat("WP_Jets","Product category for jet mult and WP");
  for (int i=0; i<NbOfJetBins; i++) {
    for (int j=0; j<NbOfWP; j++) {
      wpJetCat.defineType((std::string()+jetSuffix[i]+wpSuffix[j]).c_str(), i*NbOfWP+j);
    }
  }
  
  printf("Construct formula\n");
  
    // C o n s t r u c t   f o r m u l a s
    // -------------------------------------------------------------------------------
  RooFormulaVar ****p0bjets_tt = new RooFormulaVar***[NbOfJetBins];
  RooFormulaVar ****p1bjets_tt = new RooFormulaVar***[NbOfJetBins];
  RooFormulaVar ****p2bjets_tt = new RooFormulaVar***[NbOfJetBins];
  RooFormulaVar ****p3bjets_tt = new RooFormulaVar***[NbOfJetBins];
  RooFormulaVar ****p4bjets_tt = new RooFormulaVar***[NbOfJetBins];
  RooFormulaVar ****p5bjets_tt = new RooFormulaVar***[NbOfJetBins];
  RooFormulaVar ****p0bjets_v  = new RooFormulaVar***[NbOfJetBins];
  RooFormulaVar ****p1bjets_v  = new RooFormulaVar***[NbOfJetBins];
  RooFormulaVar ****p2bjets_v  = new RooFormulaVar***[NbOfJetBins];
  RooFormulaVar ****p3bjets_v  = new RooFormulaVar***[NbOfJetBins];
  RooFormulaVar ****p4bjets_v  = new RooFormulaVar***[NbOfJetBins];
  RooFormulaVar ****p5bjets_v  = new RooFormulaVar***[NbOfJetBins];

  for (int c=0; c<NbOfChannels; c++) {
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

    for (int i=0; i<NbOfJetBins; i++) {
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
      for (int j=0;j<NbOfWP; j++) {
        p0bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p0bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p0bjets_tt","(1-@0)*(1-@0)*pow((1-@1),@2-2)*@5+(1-@0)*pow((1-@1),@2-1)*@4+pow((1-@1),@2)*@3",RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i])));
        p1bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p1bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p1bjets_tt","(2*@0*(1-@0)*pow(1-@1,@2-2)+(1-@0)*(1-@0)*(@2-2)*@1*pow(1-@1,@2-3))*@5+(@0*pow(1-@1,@2-1)+(1-@0)*(@2-1)*@1*pow(1-@1,@2-2))*@4+(@2*@1*pow(1-@1,@2-1))*@3",RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i])));
        p2bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p2bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p2bjets_tt","(@0*@0*pow(1-@1,@2-2)+2*@0*(1-@0)*(@2-2)*@1*pow(1-@1,@2-3)+(1-@0)*(1-@0)*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4))*@5+(@0*(@2-1)*@1*pow(1-@1,@2-2)+(1-@0)*((@2-1)*(@2-2)/2)*@1*@1*pow(1-@1,@2-3))*@4+((@2*(@2-1)/2)*@1*@1*pow(1-@1,@2-2))*@3",RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i])));
        //        p3bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p3bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p3bjets_tt","(@0*@0*(@2-2)*@1*pow(1-@1,@2-3)+2*@0*(1-@0)*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4)+(@2>4 ? pow((1-@0),2)*((@2-2)*(@2-3)*(@2-4)/6)*pow(@1,3)*pow((1-@1),@2-5) : 0 ))*@5+(@0*((@2-1)*(@2-2)/2)*@1*@1*pow(1-@1,@2-3)+(1-@0)*((@2-1)*(@2-2)*(@2-3)/6)*pow(@1,3)*pow(1-@1,@2-4))*@4+((@2*(@2-1)*(@2-2)/6)*pow(@1,3)*pow(1-@1,@2-3))*@3",RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i])));
        p3bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p3bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p3bjets_tt","(@0*@0*(@2-2)*@1*pow(1-@1,@2-3)+2*@0*(1-@0)*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4)+pow((1-@0),2)*((@2-2)*(@2-3)*(@2-4)/6)*pow(@1,3)*pow((1-@1),@2-5))*@5+(@0*((@2-1)*(@2-2)/2)*@1*@1*pow(1-@1,@2-3)+(1-@0)*((@2-1)*(@2-2)*(@2-3)/6)*pow(@1,3)*pow(1-@1,@2-4))*@4+((@2*(@2-1)*(@2-2)/6)*pow(@1,3)*pow(1-@1,@2-3))*@3",RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i])));
        p4bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p4bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p4bjets_tt","(@0*@0*((@2-2)*(@2-3)/2)*@1*@1*pow(1-@1,@2-4)+2*@0*(1-@0)*((@2-2)*(@2-3)*(@2-4)/6)*pow(@1,3)*pow(1-@1,@2-5)+pow(1-@0,2)*((@2-2)*(@2-3)*(@2-4)*(@2-5)/24)*pow(@1,4)*pow((1-@1),@2-6))*@5+(@0*((@2-1)*(@2-2)*(@2-3)/6)*pow(@1,3)*pow(1-@1,@2-4)+(1-@0)*((@2-1)*(@2-2)*(@2-3)*(@2-4)/24)*pow(@1,4)*pow(1-@1,@2-5))*@4+((@2*(@2-1)*(@2-2)*(@2-3)/24)*pow(@1,4)*pow(1-@1,@2-4))*@3",RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i])));
        p5bjets_tt[c][i][j] = new RooFormulaVar((std::string()+"p5bjets_tt"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(),"p5bjets_tt","(@0*@0*((@2-2)*(@2-3)*(@2-4)/6)*pow(@1,3)*pow(1-@1,@2-5)+2*@0*(1-@0)*((@2-2)*(@2-3)*(@2-4)*(@2-5)/24)*pow(@1,4)*pow(1-@1,@2-6)+pow(1-@0,2)*((@2-2)*(@2-3)*(@2-4)*(@2-5)*(@2-6)/120)*pow(@1,5)*pow((1-@1),@2-7))*@5+(@0*((@2-1)*(@2-2)*(@2-3)*(@2-4)/24)*pow(@1,4)*pow(1-@1,@2-5)+(1-@0)*((@2-1)*(@2-2)*(@2-3)*(@2-4)*(@2-5)/120)*pow(@1,5)*pow(1-@1,@2-6))*@4+((@2*(@2-1)*(@2-2)*(@2-3)*(@2-4)/120)*pow(@1,5)*pow(1-@1,@2-5))*@3",RooArgList(*(eb[c][j]),*(eudsc[c][j]),*(n[c][i]),*(e0bq[c][i]),*(e1bq[c][i]),*(e2bq[c][i])));
        
        p0bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p0bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p0bjets_v" ,"pow(1-@0,@1)",RooArgList(*(euds[c][j]),*(n[c][i])));
        p1bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p1bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p1bjets_v" ,"@1*@0*pow(1-@0,@1-1)",RooArgList(*(euds[c][j]),*(n[c][i])));
        p2bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p2bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p2bjets_v" ,"(@1*(@1-1)/2)*@0*@0*pow(1-@0,@1-2)",RooArgList(*(euds[c][j]),*(n[c][i])));
        p3bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p3bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p3bjets_v" ,"((@1)*(@1-1)*(@1-2)/6)*pow(@0,3)*pow(1-@0,@1-3)",RooArgList(*(euds[c][j]),*(n[c][i])));
        p4bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p4bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p4bjets_v" ,"((@1)*(@1-1)*(@1-2)*(@1-3)/24)*pow(@0,4)*pow(1-@0,@1-4)",RooArgList(*(euds[c][j]),*(n[c][i])));
        p5bjets_v[c][i][j] = new RooFormulaVar((std::string()+"p5bjets_v"+channelSuffix[c]+jtSuf[i]+wpSuf[j]).c_str(), "p5bjets_v" ,"((@1)*(@1-1)*(@1-2)*(@1-3)*(@1-4)/120)*pow(@0,5)*pow(1-@0,@1-5)",RooArgList(*(euds[c][j]),*(n[c][i])));
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
  
  for (int c=0; c<NbOfChannels; c++) {
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
    
    for (int i=0; i<NbOfJetBins; i++) {
      pbjets_tt[c][i] = new RooGenericPdf*[NbOfWP];
      pbjets_v[c][i] = new RooGenericPdf*[NbOfWP];
      pbjets_tt_ext[c][i] = new RooExtendPdf*[NbOfWP];
      pbjets_v_ext[c][i]  = new RooExtendPdf*[NbOfWP];
      pbjets_st[c][i] = new RooGenericPdf*[NbOfWP];
      pbjets_vb[c][i] = new RooGenericPdf*[NbOfWP];
      pbjets_qcd[c][i] = new RooGenericPdf*[NbOfWP];
      pbjets_vv[c][i] = new RooGenericPdf*[NbOfWP];
      model[c][i] = new RooAddPdf*[NbOfWP];
      for (int j=0;j<NbOfWP; j++) {
        pbjets_tt[c][i][j] = new RooGenericPdf("pbjets_tt","pbjets_tt","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(p0bjets_tt[c][i][j]),*(p1bjets_tt[c][i][j]),*(p2bjets_tt[c][i][j]),*(p3bjets_tt[c][i][j]),*(p4bjets_tt[c][i][j]),*(p5bjets_tt[c][i][j])));
        pbjets_tt_ext[c][i][j] = new RooExtendPdf("pbjets_tt_ext","pbjets_tt_ext",*(pbjets_tt[c][i][j]),*(Ntt[c][i]));
        
        pbjets_v[c][i][j] = new RooGenericPdf("pbjets_v","pbjets_v","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(p0bjets_v[c][i][j]),*(p1bjets_v[c][i][j]),*(p2bjets_v[c][i][j]),*(p3bjets_v[c][i][j]),*(p4bjets_v[c][i][j]),*(p5bjets_v[c][i][j])));
        pbjets_v_ext[c][i][j] = new RooExtendPdf("pbjets_v_ext","pbjets_v_ext",*(pbjets_v[c][i][j]),*(Nv[c][i]));
        
        pbjets_st[c][i][j] = new RooGenericPdf("pbjets_st","pbjets_st","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(n0bjets_st[c][i][j]),*(n1bjets_st[c][i][j]),*(n2bjets_st[c][i][j]),*(n3bjets_st[c][i][j]),*(n4bjets_st[c][i][j]),*(n5bjets_st[c][i][j])));
        
        pbjets_vb[c][i][j] = new RooGenericPdf("pbjets_vb","pbjets_vb","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(n0bjets_vb[c][i][j]),*(n1bjets_vb[c][i][j]),*(n2bjets_vb[c][i][j]),*(n3bjets_vb[c][i][j]),*(n4bjets_vb[c][i][j]),*(n5bjets_vb[c][i][j])));
        pbjets_qcd[c][i][j] = new RooGenericPdf("pbjets_qcd","pbjets_qcd","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(n0bjets_qcd[c][i][j]),*(n1bjets_qcd[c][i][j]),*(n2bjets_qcd[c][i][j]),*(n3bjets_qcd[c][i][j]),*(n4bjets_qcd[c][i][j]),*(n5bjets_qcd[c][i][j])));
        
        pbjets_vv[c][i][j] = new RooGenericPdf("pbjets_vv","pbjets_vv","(nbjets==0)*@1+(nbjets==1)*@2+(nbjets==2)*@3+(nbjets==3)*@4+(nbjets==4)*@5+(nbjets==5)*@6",RooArgList(nbjets,*(n0bjets_vv[c][i][j]),*(n1bjets_vv[c][i][j]),*(n2bjets_vv[c][i][j]),*(n3bjets_vv[c][i][j]),*(n4bjets_vv[c][i][j]),*(n5bjets_vv[c][i][j])));
        
          //RooAddPdf model("model","model",RooArgList(pbjets_tt,pbjets_v),RooArgList(Ntt,Nv));
        model[c][i][j] = new RooAddPdf("model","model",RooArgList(*(pbjets_tt[c][i][j]),*(pbjets_v[c][i][j]),*(pbjets_st[c][i][j]),*(pbjets_vb[c][i][j])/*,*(pbjets_qcd[c][i][j])*/,*(pbjets_vv[c][i][j])),RooArgList(*(Ntt[c][i]),*(Nv[c][i]),*(Nst[c][i]),*(Nvb[c][i])/*,*(Nqcd[c][i])*/,*(Nvv[c][i])));
      } 
    }
    
    printf("For channel #%d\n", c);
    for (int i=0; i<NbOfJetBins; i++) {
      printf("eXbq_%d(0,1,2) = (%lf, %lf, %lf)\n", i, e0bq[c][i]->getVal(), e1bq[c][i]->getVal(), e2bq[c][i]->getVal());
    }
    printf("\n");
    
    printf("\n\nConstraint values (values not always applied as constraint ....) : \n");
    for (int i=0; i<NbOfWP; i++) {
      printf("WP %d : \n*******\n e_b = %lf +/- %lf\ne_uds  = %lf +/- %lf\ne_udsc = %lf +/- %lf\n", i, Eb_const_mean[c][i], Eb_const_error[c][i], Euds_const_mean[c][i], Euds_const_error[c][i], Eudsc_const_mean[c][i], Eudsc_const_error[c][i]);
    }
    for (int i=0; i<NbOfJetBins; i++) {
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

  
  for (int c=0; c<NbOfChannels; c++) {
    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    Eb_constraint[c] = new RooGaussian*[NbOfWP];
    Eudsc_constraint[c] = new RooGaussian*[NbOfWP];
    Euds_constraint[c] = new RooGaussian*[NbOfWP];
    for (int wp=0; wp<NbOfWP; wp++) {
      Eb_constraint[c][wp] = new RooGaussian((std::string()+"Eb_constraint"+wpSuffix[wp]).c_str(),"Eb_constraint",*(eb[c][wp]),RooConst(Eb_const_mean[c][wp]),RooConst(Eb_const_error[c][wp]));
      Eudsc_constraint[c][wp] = new RooGaussian((std::string()+"Eudsc_constraint"+wpSuffix[wp]).c_str(),"Eudsc_constraint",*(eudsc[c][wp]),RooConst(Eudsc_const_mean[c][wp]),RooConst(Eudsc_const_error[c][wp]));
      Euds_constraint[c][wp] = new RooGaussian((std::string()+"Euds_constraint"+wpSuffix[wp]).c_str(),"Euds_constraint",*(euds[c][wp]),RooConst(Euds_const_mean[c][wp]),RooConst(Euds_const_error[c][wp]));
    }
    Nst_constraint[c] = new RooGaussian*[NbOfJetBins];
    Nvb_constraint[c] = new RooGaussian*[NbOfJetBins];
    Nqcd_constraint[c] = new RooGaussian*[NbOfJetBins];
    Nvv_constraint[c] = new RooGaussian*[NbOfJetBins];
    model_constraint[c] = new RooProdPdf**[NbOfJetBins];
    
    for (int i=0; i<NbOfJetBins; i++) {
      Nst_constraint[c][i]   = new RooGaussian((std::string()+"Nst_constraint"+jetSuffix[i]).c_str(),"Nst_constraint",*(Nst[c][i]),RooConst(Nstjets[c][i]),RooConst(Nstjets[c][i]*Nstjets_uncert[c][i]));
      Nvb_constraint[c][i]   = new RooGaussian((std::string()+"Nvbb_constraint"+jetSuffix[i]).c_str(),"Nvbb_constraint",*(Nvb[c][i]),RooConst(Nvbjets[c][i]),RooConst(Nvbjets[c][i]*Nvbjets_uncert[c][i]));
      Nqcd_constraint[c][i]  = new RooGaussian((std::string()+"Nqcd_constraint"+jetSuffix[i]).c_str(),"Nqcd_constraint",*(Nqcd[c][i]),RooConst(Nqcd_mean[c][i]),RooConst(Nqcd_mean[c][i]*Nqcd_uncert[c][i]/*+100000.*/));
      Nvv_constraint[c][i]   = new RooGaussian((std::string()+"Nvv_constraint"+jetSuffix[i]).c_str(),"Nvv_constraint",*(Nvv[c][i]),RooConst(Nvvjets[c][i]),RooConst(Nvvjets[c][i]*Nvvjets_uncert[c][i]/*+100000.*/));
      model_constraint[c][i] = new RooProdPdf*[NbOfWP];
      for (int j=0 ; j<NbOfWP; j++) {
        model_constraint[c][i][j] = new RooProdPdf((std::string()+"model_constraint"+jetSuffix[i]+wpSuffix[j]).c_str(),"model with constraint",RooArgSet(*(model[c][i][j]),*(Eb_constraint[c][j]),*(Euds_constraint[c][j]),*(Eudsc_constraint[c][j]),*(Nst_constraint[c][i]),*(Nvb_constraint[c][i])/*,*(Nqcd_constraint[c][i])*/,*(Nvv_constraint[c][i]))) ;
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
  for (int c=0; c<NbOfChannels; c++) {
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
  for (int c=0; c<NbOfChannels; c++) {
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
    for (int j=0; j<NbOfJetBins; j++) {
      for (int wp=0; wp<NbOfWP; wp++) {
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
  
  for (int c=0; c<NbOfChannels; c++) {
    for (int wp=0; wp<NbOfWP; wp++) {
      printf("Efficiencies before fit    :    e_b[%d][%d] = %lf   -   e_uds[%d][%d] = %lf   -   e_udsc[%d][%d] = %lf\n", c, wp, init_Eb[c][wp], c, wp, init_Euds[c][wp], c, wp, init_Eudsc[c][wp]);
    }
  }
  
  
  
  
  
  RooFitResult **fit_result = new RooFitResult*[NbOfChannels];
    //  if(!doPE)
    //  fit_result = model_constraint.fitTo(data,Constrain(RooArgSet(eb,euds,eudsc,Nst,Nvb)),Save(1),Extended(1),Minos(1),Strategy(2),Optimize(0));
  RooArgSet** constraintsList = new RooArgSet*[NbOfChannels];
  RooArgSet** consList = new RooArgSet*[NbOfChannels];
  for (int c=0; c<NbOfChannels; c++) {
    fit_result[c] = NULL;
    
//    if (!((((c+1)&muonMask)!=0) && (((c+1)&electronMask)!=0))) continue;

    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
        
    constraintsList[c] = new RooArgSet();
    for (int wp=0; wp<NbOfWP; wp++){
      constraintsList[c]->add(*(eb[c][wp]));
      constraintsList[c]->add(*(euds[c][wp]));
      constraintsList[c]->add(*(eudsc[c][wp]));
    }
    for (int i=0; i<NbOfJetBins; i++) {
      constraintsList[c]->add(*(Nst[c][i]));
      constraintsList[c]->add(*(Nvb[c][i]));
//      constraintsList[c]->add(*(Nqcd[c][i]));
      constraintsList[c]->add(*(Nvv[c][i]));
}
    
    consList[c] = new RooArgSet();
    for (int wp=0; wp<NbOfWP; wp++){
      consList[c]->add(*(Eb_constraint[c][wp]));
      consList[c]->add(*(Euds_constraint[c][wp]));
      consList[c]->add(*(Eudsc_constraint[c][wp]));
    }
    for (int i=0; i<NbOfJetBins; i++) {
      consList[c]->add(*(Nst_constraint[c][i]));
      consList[c]->add(*(Nvb_constraint[c][i]));
//      consList[c]->add(*(Nqcd_constraint[c][i]));
      consList[c]->add(*(Nvv_constraint[c][i]));
    }
    printf("\n\n\nTTTT Input values for Fit for channel #%d\n", c);
    for (int j=0; j<NbOfWP; j++) {
      printf("  %s     : %lf\n  %s  : %lf\n  %s : %lf\n", 
             eb[c][j]->getTitle().Data(), eb[c][j]->getVal(),
             euds[c][j]->getTitle().Data(), euds[c][j]->getVal(),
             eudsc[c][j]->getTitle().Data(), eudsc[c][j]->getVal());
    }
    for (int i=0; i<NbOfJetBins; i++) {
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
    for (int j=0; j<NbOfWP; j++) {
      printf("  %s     : %lf\n  %s  : %lf\n  %s : %lf\n", 
             eb[c][j]->getTitle().Data(), eb[c][j]->getVal(),
             euds[c][j]->getTitle().Data(), euds[c][j]->getVal(),
             eudsc[c][j]->getTitle().Data(), eudsc[c][j]->getVal());
    }
    for (int i=0; i<NbOfJetBins; i++) {
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
  for (int c=0; c<NbOfChannels; c++) {
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
    for(int i=0; i<NbOfJetBins ; i++) {
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
      for (int j=0; j<NbOfWP ; j++) {
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
        
        for (int bj=0; bj<6; bj++) {
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
        
        for (int bj=0; bj<6; bj++) {
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
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  TDirectory *myDir = 0;
  
  fout->cd();
 

  for (int c=0; c<NbOfChannels; c++) {

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
    for (int i=0; i<NbOfJetBins; i++) {
      save_tt[c][i] = Ntt[c][i]->getVal();
      save_v [c][i] = Nv[c][i]->getVal();
      save_st[c][i] = Nst[c][i]->getVal();
      save_vb[c][i] = Nvb[c][i]->getVal();
      save_qcd[c][i] = Nqcd[c][i]->getVal();
      save_vv[c][i] = Nvv[c][i]->getVal();
    }
    for (int wp=0; wp<NbOfWP; wp++) {
      save_eb[c][wp] = eb[c][wp]->getVal();
      save_et[c][wp] = eudsc[c][wp]->getVal();
      save_ev[c][wp] = euds[c][wp]->getVal();
    }
    //Reset to MC values (for the plot)    
    for (int i=0; i<NbOfJetBins; i++) {
      hMC_Nominal_N[c][i]->SetBinContent(1, Nttlike[c][i]);   Ntt[c][i]->setVal(Nttlike[c][i]);
      hMC_Nominal_N[c][i]->SetBinContent(2, Nvlike[c][i]);    Nv[c][i]->setVal(Nvlike[c][i]);
      hMC_Nominal_N[c][i]->SetBinContent(3, Nstjets[c][i]);   Nst[c][i]->setVal(Nstjets[c][i]);
      hMC_Nominal_N[c][i]->SetBinContent(4, Nvbjets[c][i]);   Nvb[c][i]->setVal(Nvbjets[c][i]);
      hMC_Nominal_N[c][i]->SetBinContent(5, Nqcd_mean[c][i]); Nqcd[c][i]->setVal(Nqcd_mean[c][i]);
      hMC_Nominal_N[c][i]->SetBinContent(6, Nvvjets[c][i]);   Nvv[c][i]->setVal(Nvvjets[c][i]);
//      hMC_Nominal_N[c][i]->SetBinContent(7, n[c][i]);   n[c][i]->setVal(Nvvjets[c][i]);
    }
    for (int wp=0; wp<NbOfWP; wp++) {
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
  for (int c=0; c<NbOfChannels; c++) {
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
    for(int i=0; i<NbOfJetBins ; i++) {
      histoTT_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoV_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoVbb_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoST_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoQCD_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoVV_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoTot_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoTotC_PDF_MC[c][i] = new TH1F*[NbOfWP];
      histoTotS_PDF_MC[c][i] = new TH1F*[NbOfWP];
      for (int j=0; j<NbOfWP ; j++) {
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
  for (int c=0; c<NbOfChannels; c++) {
    mcstudy[c] = NULL;

//    if (!((((c+1)&muonMask)!=0) && (((c+1)&electronMask)!=0))) continue;

    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    printf("\n\n\nTTTT Input values for RooMCStudy for channel #%d ,\t\t  ready to generate %lf pseudo-experiments of %lf events each\n", c, (double) NbOfPE, (double) nTot[c]);
    for (int j=0; j<NbOfWP; j++) {
      printf("  %s     : %lf\n  %s  : %lf\n  %s : %lf\n", 
             eb[c][j]->getTitle().Data(), eb[c][j]->getVal(),
             euds[c][j]->getTitle().Data(), euds[c][j]->getVal(),
             eudsc[c][j]->getTitle().Data(), eudsc[c][j]->getVal());
    }
    for (int i=0; i<NbOfJetBins; i++) {
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
    for (int j=0; j<NbOfWP; j++) {
      printf("  %s     : %lf\n  %s  : %lf\n  %s : %lf\n", 
             eb[c][j]->getTitle().Data(), eb[c][j]->getVal(),
             euds[c][j]->getTitle().Data(), euds[c][j]->getVal(),
             eudsc[c][j]->getTitle().Data(), eudsc[c][j]->getVal());
    }
    for (int i=0; i<NbOfJetBins; i++) {
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
  ebTTMC_histo->Write();
  ebVMC_histo->Write();
  eudsMC_histo->Write();
  eudscMC_histo->Write();
  
  
  for (int c=0; c<NbOfChannels; c++) {

//    if (!((((c+1)&muonMask)!=0) && (((c+1)&electronMask)!=0))) continue;

    if (!semiMuon && (((c+1)&muonMask)!=0)) continue;
    if (!semiElectron && (((c+1)&electronMask)!=0)) continue;
    for(int i=0; i<NbOfJetBins; i++) {
      for(int j=0; j<NbOfWP ; j++) {
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
        hNbtaggedJets_ttlike[c][i][j]->SetMarkerSize(1.5);
        hNbtaggedJets_ttlike[c][i][j]->SetMarkerColor(kBlack);
        hNbtaggedJets_ttlike[c][i][j]->SetLineColor(kBlack);
        ttLeg->AddEntry(hNbtaggedJets_ttlike[c][i][j], "#b-tags","LP");
        mmax = histoTT_PDF_Fit[c][i][j]->GetBinContent(histoTT_PDF_Fit[c][i][j]->GetMaximumBin());
        mmin = histoTT_PDF_Fit[c][i][j]->GetBinContent(histoTT_PDF_Fit[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        histoTT_PDF_Fit     [c][i][j]->SetMarkerStyle(21);
        histoTT_PDF_Fit     [c][i][j]->SetMarkerSize(1);
        histoTT_PDF_Fit     [c][i][j]->SetMarkerColor(kRed);
        histoTT_PDF_Fit     [c][i][j]->SetLineColor(kRed);
        ttLeg->AddEntry(histoTT_PDF_Fit     [c][i][j], "PDF after fit", "LP");
        mmax = histoTT_PDF_MC[c][i][j]->GetBinContent(histoTT_PDF_MC[c][i][j]->GetMaximumBin());
        mmin = histoTT_PDF_MC[c][i][j]->GetBinContent(histoTT_PDF_MC[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        histoTT_PDF_MC      [c][i][j]->SetMarkerStyle(22);
        histoTT_PDF_MC      [c][i][j]->SetMarkerSize(1);
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
        hNbtaggedJets_vlike [c][i][j]->SetMarkerSize(1.5);
        hNbtaggedJets_vlike [c][i][j]->SetMarkerColor(kBlack);
        hNbtaggedJets_vlike [c][i][j]->SetLineColor(kBlack);
        vLeg->AddEntry(hNbtaggedJets_vlike[c][i][j], "#b-tags","LP");
        mmax = histoTT_PDF_Fit[c][i][j]->GetBinContent(histoTT_PDF_Fit[c][i][j]->GetMaximumBin());
        mmin = histoTT_PDF_Fit[c][i][j]->GetBinContent(histoTT_PDF_Fit[c][i][j]->GetMinimumBin());
        if (mmin<min) min = mmin;
        if (mmax>max) max = mmax;
        histoV_PDF_Fit      [c][i][j]->SetMarkerStyle(21);
        histoV_PDF_Fit      [c][i][j]->SetMarkerSize(1);
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
    
    
    
    for (int i=0; i<NbOfJetBins; i++) {
      for (int wp=0; wp<NbOfWP; wp++) {
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
        hWLFg_In_WJets[c][i][wp]->DrawNormalized("same")->SetLineColor(kAzure);
        hWcx_In_WJets[c][i][wp]->DrawNormalized("same")->SetLineColor(kViolet);
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
        hNbtaggedJets_qcd[c][i][wp]->DrawNormalized("same")->SetLineColor(kPink);
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
        hNbtaggedJets_qcd[c][i][wp]->DrawNormalized("same")->SetLineColor(kPink);
        c2->Update();
        c2->Write();
                
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
      for (int i=0; i<NbOfJetBins; i++) {
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
      
      for (int j=0; j<NbOfWP; j++) {
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
  
  printf("Triggers for TTjets\n");
  for (int i=1; i<=histoTrig_ttjets->GetNbinsX(); i++) {
    if (histoTrig_ttjets->GetBinContent(i)!=0) {
      printf(" -> %s : %lf\n", histoTrig_ttjets->GetXaxis()->GetBinLabel(i), histoTrig_ttjets->GetBinContent(i));
    }
  }
  histoTrig_ttjets->Write();
  printf("\nTriggers for Wjets\n");
  for (int i=1; i<=histoTrig_wjets->GetNbinsX(); i++) {
    if (histoTrig_wjets->GetBinContent(i)!=0) {
      printf(" -> %s : %lf\n", histoTrig_wjets->GetXaxis()->GetBinLabel(i), histoTrig_wjets->GetBinContent(i));
    }
  }
  histoTrig_wjets->Write();
  printf("\nTriggers for Zjets\n");
  for (int i=1; i<=histoTrig_zjets->GetNbinsX(); i++) {
    if (histoTrig_zjets->GetBinContent(i)!=0) {
      printf(" -> %s : %lf\n", histoTrig_zjets->GetXaxis()->GetBinLabel(i), histoTrig_zjets->GetBinContent(i));
    }
  }
  histoTrig_zjets->Write();
  printf("\nTriggers for STjets\n");
  for (int i=1; i<=histoTrig_stjets->GetNbinsX(); i++) {
    if (histoTrig_stjets->GetBinContent(i)!=0) {
      printf(" -> %s : %lf\n", histoTrig_stjets->GetXaxis()->GetBinLabel(i), histoTrig_stjets->GetBinContent(i));
    }
  }
  histoTrig_stjets->Write();
  
  
  
  
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
 gcc -O $(root-config --cflags) $(root-config --libs)  -lMathMore -lFoam -lMinuit -lRooFitCore -lRooFit -o StopSearches_VJetsBckgdEst_NoFiles_SevWP.exe StopSearches_VJetsBckgdEst_NoFiles_SevWP.cc
 */

int main (int argc, char *argv[])
{
  return StopSearches_VJetsBckgdEst_NoFiles_SevJets();
}

