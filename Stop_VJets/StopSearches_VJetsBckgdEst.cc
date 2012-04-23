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

#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooNLLVar.h"

#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooTable.h"
#include "Roo1DTable.h"

using namespace RooFit;
using namespace std;

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

int StopSearches_VJetsBckgdEst()
{

  clock_t start = clock();

  cout << "***************************************************************" << endl;
  cout << " Beginning of the program for the V+Jets background estimation " << endl;
  cout << "***************************************************************" << endl;

  bool verbose = true;
  bool runOnPNFS = false;
  bool Wbb_from_WJets = true;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::string WbbShapeRootfile = "StopBckg_Outputtttttt.root";
  
  /** TO BE ADAPTED **/
  TChain *ttjets = new TChain("AnaTree");
  if (runOnPNFS) {
    ttjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/big_TTbar_skim.root");
  } else {
    ttjets->Add("/user/nstrobbe/DATA/stop/skims_443/TTbar_skim.root");
  }
  TEntryList *entrylist_ttjets;
  vector<double>   *vjtDiscri_pf_ttjets = 0;

  TChain  *wjets = new TChain("AnaTree");
  if (runOnPNFS) {
    wjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/WJets_skim.root");
  } else {
    wjets->Add("/user/nstrobbe/DATA/stop/skims_443/WJets_skim.root");
  }
  TEntryList *entrylist_wjets  = 0;
  vector<double>   *vjtDiscri_pf_wjets  = 0;
  vector<double>   *vjtFlav_pf_wjets    = 0;

  TChain  *zjets = new TChain("AnaTree");
  if (runOnPNFS) {
    zjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/ZJets_skim.root");
  } else {
    zjets->Add("/user/nstrobbe/DATA/stop/skims_443/ZJets_skim.root");
  }
  TEntryList *entrylist_zjets  = 0;
  vector<double>   *vjtDiscri_pf_zjets  = 0;

  TChain *stjets = new TChain("AnaTree");
  if (runOnPNFS) {
    stjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/T_tW_skim.root");
    stjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/Tbar_tW_skim.root");
  } else {
    stjets->Add("/user/nstrobbe/DATA/stop/skims_443/T_tW_skim.root");
    stjets->Add("/user/nstrobbe/DATA/stop/skims_443/Tbar_tW_skim.root");
  }
  TEntryList *entrylist_stjets = 0;
  vector<double>   *vjtDiscri_pf_stjets = 0;

//  TChain *vbjets = new TChain("AnaTree");
//  vbjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/XXX_skim.root");
//  TEntryList *entrylist_vbjets = 0;
//  vector<double>   *vjtDiscri_pf_vbjets = 0;

//  TBranch        *b_vjtDiscri_pf = 0;


  if(verbose){
    cout<<"Number of events for :"<<endl;
    cout<<" - ttjets = "<<ttjets->GetEntries()<<endl;
    cout<<" -  wjets = "<< wjets->GetEntries()<<endl;
    cout<<" -  zjets = "<< zjets->GetEntries()<<endl;
    cout<<" - stjets = "<<stjets->GetEntries()<<endl;
//    cout<<" - vbjets ="<<vbjets->GetEntries()<<endl;
  }

  // Processes XS (pb) taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSections
  /*** TO BE UPDATED TO THE MEASURED XS WHENEVER POSSIBLE ***/
  float XS_ttjets = 165.;
  float XS_wjets  = 31314.;
  float XS_zjets  = 3048.;
  float XS_stjets = 10.6;
  float XS_vbjets = 35.5;
  
  // Retrieve int. lumi from InputTree (?)
  const double IntLumi = 5000;
  
  // Normalization factors for a given int. lumi.
  /** TO BE ADAPTED **/
  bool skimmedTree = true;
  double NormFact_ttjets = 0;
  double NormFact_wjets  = 0;
  double NormFact_zjets  = 0;
  double NormFact_stjets = 0;
//  double NormFact_vbjets = 0;
  
  if(skimmedTree) {
    ttjets->SetBranchAddress("eventWeight",&NormFact_ttjets);
     wjets->SetBranchAddress("eventWeight",&NormFact_wjets);
     zjets->SetBranchAddress("eventWeight",&NormFact_zjets);
    stjets->SetBranchAddress("eventWeight",&NormFact_stjets);
    //vbjets->SetBranchAddress("eventWeight",&NormFact_vbjets);
  }
  else{
    NormFact_ttjets = XS_ttjets/ttjets->GetEntries();
    NormFact_wjets  = XS_wjets/wjets->GetEntries();
    NormFact_zjets  = XS_zjets/zjets->GetEntries();
    NormFact_stjets = XS_stjets/stjets->GetEntries();
//    NormFact_vbjets = XS_vbjets/vbjets->GetEntries();
  }

  //Output ROOT file
  string postfix = "_VJetsBckgdEst";
  string channelpostfix = "_SemiMuon";
  string comment = "_Constants_euds";
  string rootFileName ("StopSearches"+postfix+channelpostfix+comment+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  TDirectory *myDir = 0;

  // E v e n t   s e l e c t i o n   c u t s
  // ---------------------------------------

  const int NbOfPE = 1000;
  const int NbOfJetBins = 4;
  const int MinNbOfJets = 3;
  int JetIdx = 1;
  /** B-tagging working points ; https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP#B_tagging_Operating_Points_for_3
  TrackCountingHighEff  	  TCHEL 1.7
  TrackCountingHighEff 	    TCHEM 3.3
  TrackCountingHighEff 	    TCHET 10.2 (not supported)
  CombinedSecondaryVertex  	CSVL  0.244
  CombinedSecondaryVertex 	CSVM 	0.679
  CombinedSecondaryVertex 	CSVT 	0.898 
  **/
  float btagCut = 1.7;
  // Muon channel
  TCut muon = "(muPt[0]>20.)&&(TMath::Abs(muEta[0])<=2.1)";
  // Electron channel
  //TCut electron = "(elPt[0]>0.)&&(elSCenergy[0]>=17.0)";
  // Lepton selection
  TCut lveto = "NmuIdIso==1";

  // Channel selection
  TCut singlep = muon+lveto;
    
  // Jet selection
  TCut pfnjet[NbOfJetBins];
  pfnjet[0] = "Njt_pf==3";
  pfnjet[1] = "Njt_pf==4";
  pfnjet[2] = "Njt_pf==5";
  pfnjet[3] = "Njt_pf>=6";
  //TCut pf4jet = "Njt_pf>3";
  //TCut btag   = "Njtbtag_pf>0";
  
  // Selection cuts
  TCut cuts = (skimmedTree ? pfnjet[JetIdx] : singlep+pfnjet[JetIdx]);

  // S e t t i n g s    f o r   m a x i m u m   l i k e l i h o o d   e s t i m a t i o n
  // ------------------------------------------------------------------------------------

  // Probabilities to have X b-quarks in the tt-like final state **DERIVED FROM MC SIMULATIONS**
  const float e0bq_[NbOfJetBins] = {0.0247448,0.00873498,0.00638495,0.00734619};
  const float e1bq_[NbOfJetBins] = {0.379931,0.179928,0.145577,0.134986};
  const float e2bq_[NbOfJetBins] = {0.594488,0.808598,0.83759,0.836547};

  // Initial values for minimization
  float init_Nttlike[NbOfJetBins] = {0.,24500.,9850.,4015.};
  float init_Nvlike[NbOfJetBins]  = {0.,23400.,3900., 750.};
  float init_Eb    = 0.8;  /** TO BE ADAPTED IF USED AS FIXED VARIABLE **/
  float init_Eudsc = 0.2;  /** TO BE ADAPTED IF USED AS FIXED VARIABLE **/
  float init_Euds  = 0.179; /** TO BE ADAPTED IF USED AS FIXED VARIABLE **/

  // Mean and error values for constraints
  float Eb_const_mean     = 0.795;     /** TO BE CALCULATED BY FUNCT FROM BTV-11-003 **/
  float Eb_const_error    = 0.795*0.1; /** TO BE CALCULATED BY FUNCT FROM BTV-11-003 **/
  float Eudsc_const_mean  = 0.000;     /** TO BE DERIVED FROM MC **/
  float Eudsc_const_error = 0.000;     /** TO BE DERIVED FROM MC **/
  float Euds_const_mean   = 0.179;     /** TO BE CALCULATED BY FUNCT FROM BTV-11-003 **/
  float Euds_const_error  = 0.179*0.1; /** TO BE CALCULATED BY FUNCT FROM BTV-11-003 **/


  
    // Counters for the b-tagged jet multiplicity
  int NbtaggedJets = 0;
  int nB_Flavoured_Jets = 0;
  int nC_Flavoured_Jets = 0;
  int nLightAndGluon_Flavoured_Jets = 0;
  // Histograms for the b-tagged jet multiplicity
  TH1F* hNbtaggedJets_ttjets = new TH1F("hNbtaggedJets_ttjets",";B-tagged jet multiplicity",10,0,10);
  TH1F* hNbtaggedJets_wjets  = new TH1F("hNbtaggedJets_wjets", ";B-tagged jet multiplicity",10,0,10);
  TH1F* hNbtaggedJets_zjets  = new TH1F("hNbtaggedJets_zjets", ";B-tagged jet multiplicity",10,0,10);
  TH1F* hNbtaggedJets_stjets = new TH1F("hNbtaggedJets_stjets",";B-tagged jet multiplicity",10,0,10);
  TH1F* hNbtaggedJets_vbjets = new TH1F("hNbtaggedJets_vbjets",";B-tagged jet multiplicity",10,0,10);

  TH1F* hNbtaggedJets_ttlike = new TH1F("hNbtaggedJets_ttlike",";B-tagged jet multiplicity",10,0,10);
  TH1F* hNbtaggedJets_vlike  = new TH1F("hNbtaggedJets_vlike" ,";B-tagged jet multiplicity",10,0,10);
  TH1F* hNbtaggedJets        = new TH1F("hNbtaggedJets"       ,";B-tagged jet multiplicity",10,0,10);

  TH1F* hWbb_In_WJets        = new TH1F("hNbtaggedJets_wjets_Wbb",";B-tagged jet multiplicity",10,0,10);
  TH1F* hWcx_In_WJets        = new TH1F("hNbtaggedJets_wjets_Wcx",";B-tagged jet multiplicity",10,0,10);
  TH1F* hWLFg_In_WJets       = new TH1F("hNbtaggedJets_wjets_WLFg",";B-tagged jet multiplicity",10,0,10);
  
  vector<int>::iterator Idx;
  vector<int> FixedVarIdx;
  //FixedVarIdx.push_back(0); // Fix eb, b-tagging efficiency
  //FixedVarIdx.push_back(1); // Fix eudsc, mis-tagging efficiency for tt-like events
  //FixedVarIdx.push_back(2); // Fix euds, mis-tagging efficiency for v-like events

  string wp[3]  = {"loose","medium","tight"};
  bool doPE = false;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Analysis ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // E v e n t   s e l e c t i o n
  // -----------------------------
  ttjets->SetBranchAddress("vjtDiscri_pf", &vjtDiscri_pf_ttjets);//, &b_vjtDiscri_pf);
   wjets->SetBranchAddress("vjtDiscri_pf", &vjtDiscri_pf_wjets);//, &b_vjtDiscri_pf);
   zjets->SetBranchAddress("vjtDiscri_pf", &vjtDiscri_pf_zjets);//, &b_vjtDiscri_pf);
  stjets->SetBranchAddress("vjtDiscri_pf", &vjtDiscri_pf_stjets);//, &b_vjtDiscri_pf);
//  vbjets->SetBranchAddress("vjtDiscri_pf", &vjtDiscri_pf_vbjets);//, &b_vjtDiscri_pf);

  wjets->SetBranchAddress("vjtFlav_pf", &vjtFlav_pf_wjets);//, &vjtFlav_pf);

  ttjets->Draw(">>entrylist_ttjets", cuts, "entrylist");
   wjets->Draw(">>entrylist_wjets",  cuts, "entrylist");
   zjets->Draw(">>entrylist_zjets",  cuts, "entrylist");
  stjets->Draw(">>entrylist_stjets", cuts, "entrylist");
//  vbjets->Draw(">>entrylist_vbjets", cuts, "entrylist");
   
  entrylist_ttjets = (TEntryList*)gDirectory->Get("entrylist_ttjets");
  entrylist_wjets  = (TEntryList*)gDirectory->Get("entrylist_wjets");
  entrylist_zjets  = (TEntryList*)gDirectory->Get("entrylist_zjets");
  entrylist_stjets = (TEntryList*)gDirectory->Get("entrylist_stjets");
//  entrylist_vbjets = (TEntryList*)gDirectory->Get("entrylist_vbjets");
//BBB
  ttjets->SetEntryList(entrylist_ttjets);
   wjets->SetEntryList(entrylist_wjets);
   zjets->SetEntryList(entrylist_zjets);
  stjets->SetEntryList(entrylist_stjets);
//  vbjets->SetEntryList(entrylist_vbjets);

  Long64_t NbOfEvtsToProcess = 0;
  // For tt+jets events 
  NbOfEvtsToProcess = entrylist_ttjets->GetN();
  if(verbose) cout<<"Processing tt+jets Tree (Nb of events = "<<NbOfEvtsToProcess<<")"<<endl;
  for (Long64_t el = 0; el<NbOfEvtsToProcess; el++) {
    if(el%100 == 0) cout<<"-- Processing the "<<el<<"th event ("<<100*(el)/(NbOfEvtsToProcess)<<"%)"<<flush<<"\r";
    ttjets->GetEntry(entrylist_ttjets->GetEntry(el));
    NbtaggedJets = 0;
    for(UInt_t i = 0; i<vjtDiscri_pf_ttjets->size(); ++i)
    {
      if((*vjtDiscri_pf_ttjets)[i]>btagCut) NbtaggedJets++;
    }
    hNbtaggedJets_ttjets->Fill(NbtaggedJets);
  }
  hNbtaggedJets_ttjets->Scale(NormFact_ttjets*IntLumi);
  // For w+jets events 
  NbOfEvtsToProcess = entrylist_wjets->GetN();
  if(verbose) cout<<"Processing W+jets Tree (Nb of events = "<<NbOfEvtsToProcess<<")"<<endl;
  for (Long64_t el = 0; el<entrylist_wjets->GetN();el++) {
    if(el%100 == 0) cout<<"-- Processing the "<<el<<"th event ("<<100*(el)/(NbOfEvtsToProcess)<<"%)"<<flush<<"\r";
    wjets->GetEntry(entrylist_wjets->GetEntry(el));
    NbtaggedJets = 0;
    nB_Flavoured_Jets = 0;
    nC_Flavoured_Jets = 0;
    nLightAndGluon_Flavoured_Jets = 0;
    for(UInt_t i = 0; i<vjtDiscri_pf_wjets->size(); ++i)
    {
      if((*vjtDiscri_pf_wjets)[i]>btagCut) NbtaggedJets++;
      if (std::fabs((*vjtFlav_pf_wjets)[i])==5.) {
        nB_Flavoured_Jets++;
      } else if (std::fabs((*vjtFlav_pf_wjets)[i])==4.) {
        nC_Flavoured_Jets++;
      } else {
        nLightAndGluon_Flavoured_Jets++;
      }
    }
    hNbtaggedJets_wjets->Fill(NbtaggedJets);
    if (nB_Flavoured_Jets>0) {
      hWbb_In_WJets->Fill(NbtaggedJets);
    } else if (nC_Flavoured_Jets>0) {
      hWcx_In_WJets->Fill(NbtaggedJets);
    } else {
      hWLFg_In_WJets->Fill(NbtaggedJets);
    }
  }
  hNbtaggedJets_wjets->Scale(NormFact_wjets*IntLumi);
  // For z+jets events 
  NbOfEvtsToProcess = entrylist_zjets->GetN();
  if(verbose) cout<<"Processing Z+jets Tree (Nb of events = "<<NbOfEvtsToProcess<<")"<<endl;
  for (Long64_t el = 0; el<entrylist_zjets->GetN();el++) {
    if(el%100 == 0) cout<<"-- Processing the "<<el<<"th event ("<<100*(el)/(NbOfEvtsToProcess)<<"%)"<<flush<<"\r";
    zjets->GetEntry(entrylist_zjets->GetEntry(el));
    NbtaggedJets = 0;
    for(UInt_t i = 0; i<vjtDiscri_pf_zjets->size(); ++i)
    {
      if((*vjtDiscri_pf_zjets)[i]>btagCut) NbtaggedJets++;
    }
    hNbtaggedJets_zjets->Fill(NbtaggedJets);


    
  }
  hNbtaggedJets_zjets->Scale(NormFact_zjets*IntLumi);
  // For st+jets events 
  NbOfEvtsToProcess = entrylist_stjets->GetN();
  if(verbose) cout<<"Processing st+jets Tree (Nb of events = "<<NbOfEvtsToProcess<<")"<<endl;
  for (Long64_t el = 0; el<entrylist_stjets->GetN();el++) {
    if(el%100 == 0) cout<<"-- Processing the "<<el<<"th event ("<<100*(el)/(NbOfEvtsToProcess)<<"%)"<<flush<<"\r";
    stjets->GetEntry(entrylist_stjets->GetEntry(el));
    NbtaggedJets = 0;
    for(UInt_t i = 0; i<vjtDiscri_pf_stjets->size(); ++i)
    {
      if((*vjtDiscri_pf_stjets)[i]>btagCut) NbtaggedJets++;
    }
    hNbtaggedJets_stjets->Fill(NbtaggedJets);
  }
  hNbtaggedJets_stjets->Scale(NormFact_stjets*IntLumi);

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
  // - tt-like processes
  hNbtaggedJets_ttlike->Add(hNbtaggedJets_ttjets);
  // - v-like processes
  hNbtaggedJets_vlike->Add(hNbtaggedJets_wjets);
  hNbtaggedJets_vlike->Add(hNbtaggedJets_zjets);

  // Set true numbers of tt-like and v-like events
  float Nttlike = hNbtaggedJets_ttlike->Integral();
  float Nvlike  = hNbtaggedJets_vlike->Integral();

  if (Wbb_from_WJets) {
    for (int k =0; k<10; k++) {
      hNbtaggedJets_vbjets->SetBinContent(k+1, hWbb_In_WJets->GetBinContent(hWbb_In_WJets->FindBin(k)));
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
        if (JetIdx==MinNbOfJets+NbOfJetBins-1) {
          sprintf(name, "HistoryFlavor--AtLeast_%d_jets--FlavorHistoryPath_%d", JetIdx+MinNbOfJets, j);
        } else {
          sprintf(name, "HistoryFlavor--%d_jets--FlavorHistoryPath_%d", JetIdx+MinNbOfJets, j);
        }
        TH1I *h = (TH1I*) f->Get(name);
        for (UInt_t k=1; k<=Nbins; k++) {
          hNbtaggedJets_vbjets->Fill(k-1,h->GetBinContent(k));
        }
      }
      f->Close();
    }  
  }

    // Set numbers of background events and associated uncertainties
  float Nstjets = hNbtaggedJets_stjets->Integral();
  float Nstjets_uncert = 0.3; // 30% uncertainty on the st+jets XS ** TO BE ADAPTED **
  float Nvbjets = hNbtaggedJets_vbjets->Integral();
  float Nvbjets_uncert = 0.3; // 30% uncertainty on the vb+jets XS ** TO BE ADAPTED **
  
  
  
  
  // Sum over all datasets  
  hNbtaggedJets->Add(hNbtaggedJets_ttlike);
  hNbtaggedJets->Add(hNbtaggedJets_vlike);
  hNbtaggedJets->Add(hNbtaggedJets_stjets);
//  hNbtaggedJets->Add(hNbtaggedJets_vbjets);

  // Total nb of events (rescaled to the required int. lumi.)
  int nExpected = (int)hNbtaggedJets->Integral();

 
  if(verbose){
    cout<<"***************************************"<<endl;
    cout<<"- Number of entries / NormFactors "<<endl;
    cout<<"-- for tt-like processes : "<<endl;
    cout<<"--- tt+jets = "<<hNbtaggedJets_ttjets->GetEntries()<<"/ "<<NormFact_ttjets<<endl;
    cout<<"--- total nb of events (for "<<IntLumi<<" /pb) = "<<Nttlike<<endl;
    cout<<"---------------------------------"<<endl;
    cout<<"-- for v-like processes : "<<endl;
    cout<<"--- w+jets = "<<hNbtaggedJets_wjets->GetEntries()<<"/ "<<NormFact_wjets<<endl;
    cout<<"--- z+jets = "<<hNbtaggedJets_zjets->GetEntries()<<"/ "<<NormFact_zjets<<endl;
    cout<<"--- total nb of events (for "<<IntLumi<<" /pb) = "<<Nvlike<<endl;
    cout<<"---------------------------------"<<endl;
    cout<<"-- for background processes : "<<endl;
    cout<<"--- st+jets = "<<hNbtaggedJets_stjets->GetEntries()<<"/ "<<NormFact_stjets<<endl;
    cout<<"--- total nb of events (for "<<IntLumi<<" /pb) = "<<Nstjets<<endl;
//    cout<<"--- vb+jets = "<<hNbtaggedJets_vbjets->GetEntries()<<"/ "<<NormFact_vbjets<<endl;
//    cout<<"--- total nb of events (for "<<IntLumi<<" /pb) = "<<Nvbjets<<endl;
    cout<<"---------------------------------"<<endl;
    cout<<"-- for all processes : "<<endl;
    cout<<"--- total = "<<hNbtaggedJets->GetEntries()<<endl;
    cout<<"--- total nb of events (for "<<IntLumi<<" /pb) = "<<nExpected<<endl;
    cout<<"***************************************"<<endl;
  }

  // C r e a t e   o b s e r v a b l e
  // ---------------------------------

  // Estimation parameters
  RooRealVar Ntt("Ntt","N_{t#bar{t}-like}",init_Nttlike[JetIdx],0.0,init_Nttlike[JetIdx]*3);
  RooRealVar Nv("Nv","N_{V-like}",init_Nvlike[JetIdx],0.0,init_Nvlike[JetIdx]*4);

  // Background normalization factors
  RooRealVar Nst("Nst","N_{single top}",Nstjets,0.0,Nstjets*4);
  RooRealVar Nvb("Nvb","N_{V+b-jets}",Nvbjets,0.0,Nvbjets*4);

	RooRealVar  eb(    "eb",  "#epsilon_{b-tag}",      init_Eb,0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 0 );
	if(Idx != FixedVarIdx.end()) eb.setConstant(kTRUE) ;
	RooRealVar  eudsc("eudsc","#epsilon_{mis-tag}",    init_Eudsc,0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 1 );
	if(Idx != FixedVarIdx.end()) eudsc.setConstant(kTRUE) ;
	RooRealVar euds( "euds",  "#epsilon^{,}_{mis-tag}",init_Euds,0.0,1.0);
	Idx = find( FixedVarIdx.begin(), FixedVarIdx.end(), 2 );
	if(Idx != FixedVarIdx.end()) euds.setConstant(kTRUE) ;
  
  // Declare constants
  
	RooConstVar n("n","number of selected jets",MinNbOfJets+JetIdx) ;
	RooConstVar e0bq("e0bq","e0bq",e0bq_[JetIdx]);
	RooConstVar e1bq("e1bq","e1bq",e1bq_[JetIdx]);
	RooConstVar e2bq("e2bq","e2bq",e2bq_[JetIdx]);
  //RooConstVar e3bq("e3bq","e3bq",e3bq_[njets]);
  
	RooConstVar n0bjets_st("n0bjets_st","n^{0 b-jet}", hNbtaggedJets_stjets->GetBinContent(hNbtaggedJets_stjets->FindBin(0.)));
	RooConstVar n1bjets_st("n1bjets_st","n^{1 b-jet}", hNbtaggedJets_stjets->GetBinContent(hNbtaggedJets_stjets->FindBin(1.)));
	RooConstVar n2bjets_st("n2bjets_st","n^{2 b-jets}",hNbtaggedJets_stjets->GetBinContent(hNbtaggedJets_stjets->FindBin(2.)));
	RooConstVar n3bjets_st("n3bjets_st","n^{3 b-jets}",hNbtaggedJets_stjets->GetBinContent(hNbtaggedJets_stjets->FindBin(3.)));

	RooConstVar n0bjets_vb("n0bjets_vb","n^{0 b-jet}", hNbtaggedJets_vbjets->GetBinContent(hNbtaggedJets_vbjets->FindBin(0.)));
	RooConstVar n1bjets_vb("n1bjets_vb","n^{1 b-jet}", hNbtaggedJets_vbjets->GetBinContent(hNbtaggedJets_vbjets->FindBin(1.)));
	RooConstVar n2bjets_vb("n2bjets_vb","n^{2 b-jets}",hNbtaggedJets_vbjets->GetBinContent(hNbtaggedJets_vbjets->FindBin(2.)));
	RooConstVar n3bjets_vb("n3bjets_vb","n^{3 b-jets}",hNbtaggedJets_vbjets->GetBinContent(hNbtaggedJets_vbjets->FindBin(3.)));

  // C o n s t r u c t   a   c a t e g o r y   w i t h   l a b e l s    a n d   i n d e c e s
  // -----------------------------------------------------------------------------------------

	RooCategory nbjets("nbjets","Number of b-jets");
	nbjets.defineType("N0bjet", 0);
	nbjets.defineType("N1bjet", 1);
	nbjets.defineType("N2bjets",2);
	nbjets.defineType("N3bjets",3);

  // C o n s t r u c t   f o r m u l a s
  // -------------------------------------------------------------------------------
  
  RooFormulaVar p0bjets_tt("p0bjets_tt","p0bjets_tt","(1-eb)*(1-eb)*pow((1-eudsc),n-2)*e2bq+(1-eb)*pow((1-eudsc),n-1)*e1bq+pow((1-eudsc),n)*e0bq",RooArgList(eb,eudsc,n,e0bq,e1bq,e2bq));
  RooFormulaVar p1bjets_tt("p1bjets_tt","p1bjets_tt","(2*eb*(1-eb)*pow(1-eudsc,n-2)+(1-eb)*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3))*e2bq+(eb*pow(1-eudsc,n-1)+(1-eb)*(n-1)*eudsc*pow(1-eudsc,n-2))*e1bq+(n*eudsc*pow(1-eudsc,n-1))*e0bq",RooArgList(eb,eudsc,n,e0bq,e1bq,e2bq));
  RooFormulaVar p2bjets_tt("p2bjets_tt","p2bjets_tt","(eb*eb*pow(1-eudsc,n-2)+2*eb*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3)+(1-eb)*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4))*e2bq+(eb*(n-1)*eudsc*pow(1-eudsc,n-2)+(1-eb)*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*e1bq+((n*(n-1)/2)*eudsc*eudsc*pow(1-eudsc,n-2))*e0bq",RooArgList(eb,eudsc,n,e0bq,e1bq,e2bq));
  RooFormulaVar p3bjets_tt("p3bjets_tt","p3bjets_tt","(eb*eb*(n-2)*eudsc*pow(1-eudsc,n-3)+2*eb*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4)+(n>4 ? pow((1-eb),2)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow((1-eudsc),n-5) : 0 ))*e2bq+(eb*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3)+(1-eb)*((n-1)*(n-2)*(n-3)/6)*pow(eudsc,3)*pow(1-eudsc,n-4))*e1bq+((n*(n-1)*(n-2)/6)*pow(eudsc,3)*pow(1-eudsc,n-3))*e0bq",RooArgList(eb,eudsc,n,e0bq,e1bq,e2bq));
  
  RooFormulaVar p0bjets_v ("p0bjets_v", "p0bjets_v" ,"pow(1-euds,n)",RooArgList(euds,n));
  RooFormulaVar p1bjets_v ("p1bjets_v", "p1bjets_v" ,"n*euds*pow(1-euds,n-1)",RooArgList(euds,n));
  RooFormulaVar p2bjets_v ("p2bjets_v", "p2bjets_v" ,"(n*(n-1)/2)*euds*euds*pow(1-euds,n-2)",RooArgList(euds,n));
  RooFormulaVar p3bjets_v ("p3bjets_v", "p3bjets_v" ,"((n)*(n-1)*(n-2)/6)*pow(euds,3)*pow(1-euds,n-3)",RooArgList(euds,n));
  
  // C o n s t r u c t   p . d . f 's
  // -------------------------------------------------------------------------------
  
	RooGenericPdf pbjets_tt("pbjets_tt","pbjets_tt","(nbjets==0)*p0bjets_tt+(nbjets==1)*p1bjets_tt+(nbjets==2)*p2bjets_tt+(nbjets==3)*p3bjets_tt",RooArgList(nbjets,p0bjets_tt,p1bjets_tt,p2bjets_tt,p3bjets_tt));
	RooExtendPdf  pbjets_tt_ext("pbjets_tt_ext","pbjets_tt_ext",pbjets_tt,Ntt);
  
	RooGenericPdf pbjets_v("pbjets_v","pbjets_v","(nbjets==0)*p0bjets_v+(nbjets==1)*p1bjets_v+(nbjets==2)*p2bjets_v+(nbjets==3)*p3bjets_v",RooArgList(nbjets,p0bjets_v,p1bjets_v,p2bjets_v,p3bjets_v));
	RooExtendPdf  pbjets_v_ext("pbjets_v_ext","pbjets_v_ext",pbjets_v,Nv);
  
  RooGenericPdf pbjets_st("pbjets_st","pbjets_st","(nbjets==0)*n0bjets_st+(nbjets==1)*n1bjets_st+(nbjets==2)*n2bjets_st+(nbjets==3)*n3bjets_st",RooArgList(nbjets,n0bjets_st,n1bjets_st,n2bjets_st,n3bjets_st));

  RooGenericPdf pbjets_vb("pbjets_vb","pbjets_vb","(nbjets==0)*n0bjets_vb+(nbjets==1)*n1bjets_vb+(nbjets==2)*n2bjets_vb+(nbjets==3)*n3bjets_vb",RooArgList(nbjets,n0bjets_vb,n1bjets_vb,n2bjets_vb,n3bjets_vb));

    //RooAddPdf model("model","model",RooArgList(pbjets_tt,pbjets_v),RooArgList(Ntt,Nv));
  RooAddPdf model("model","model",RooArgList(pbjets_tt,pbjets_v,pbjets_st,pbjets_vb),RooArgList(Ntt,Nv,Nst,Nvb));
	
  // Contrainsts on free parameters
  //RooGaussian Eb_constraint("Eb_constraint","Eb_constraint",eb,RooConst(Eb_const_mean),RooConst(Eb_const_error));
  //RooGaussian Eudsc_constraint("Eudsc_constraint","Eudsc_constraint",eudsc,RooConst(XXX),RooConst(XXX));
  RooGaussian Euds_constraint("Euds_constraint","Euds_constraint",euds,RooConst(Euds_const_mean),RooConst(Euds_const_error));

  RooGaussian Nst_constraint("Nst_constraint","Nst_constraint",Nst,RooConst(Nstjets),RooConst(Nstjets*Nstjets_uncert));
  
  RooGaussian Nvb_constraint("Nvbb_constraint","Nvbb_constraint",Nvb,RooConst(Nvbjets),RooConst(Nvbjets*Nvbjets_uncert));
  
  //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,Eb_constraint)) ;
  //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,Eudsc_constraint)) ;
  //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,Euds_constraint)) ;
  RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,Euds_constraint,Nst_constraint,Nvb_constraint)) ;

  // C r e a t e   d a t a s e t 
  // -------------------------------------------------------------------------------

	RooDataSet data("data","data",RooArgSet(nbjets)) ;
  for (Int_t i=0 ; i<hNbtaggedJets->GetBinContent(hNbtaggedJets->FindBin(0.)) ; i++) { nbjets.setLabel("N0bjet")  ; data.add(RooArgSet(nbjets));}
  for (Int_t i=0 ; i<hNbtaggedJets->GetBinContent(hNbtaggedJets->FindBin(1.)) ; i++) { nbjets.setLabel("N1bjet")  ; data.add(RooArgSet(nbjets));}
  for (Int_t i=0 ; i<hNbtaggedJets->GetBinContent(hNbtaggedJets->FindBin(2.)) ; i++) { nbjets.setLabel("N2bjets") ; data.add(RooArgSet(nbjets));}
  for (Int_t i=0 ; i<hNbtaggedJets->GetBinContent(hNbtaggedJets->FindBin(3.)) ; i++) { nbjets.setLabel("N3bjets") ; data.add(RooArgSet(nbjets));}
  
  //data.Print("v");
  Roo1DTable* table = data.table(nbjets);
  table->Print("v");

  // F i t   t h e   d a t a   a n d   c o n s t r u c t   t h e   l i k e l i h o o d   f u n c t i o n
  // ----------------------------------------------------------------------------------------------

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
	
  if(doPE) RooFitResult* fit_result = model_constraint.fitTo(data,Constrain(RooArgSet(Euds_constraint,Nst_constraint)),Save(1),Extended(1),Minos(1),Strategy(2));//,Constrain(RooArgSet(Euds_constraint,Nst_constraint)));
  //,ExternalConstraints(RooArgSet(Eb_constraint,Eudsc_constraint,Euds_constraint)));

	RooFitResult* fit_result = minimizer.save();
	
	fit_result->Print("v");


  RooMCStudy* mcstudy = new RooMCStudy(model,nbjets,Constrain(RooArgSet(Euds_constraint,Nst_constraint)),Binned(kTRUE),Silence(),Extended(kTRUE),FitOptions(Optimize(0),Save(0),Minos(1),Extended(1),Strategy(2),PrintEvalErrors(-1)));

  // Generate and fit NbOfPE samples of Poisson(nExpected) events
  if(doPE){
    mcstudy->generateAndFit(NbOfPE,nExpected) ;
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
  fout->cd();

  //frame1->Write();
  //hNbtaggedJets_ttjets->Draw("COLZ");
  hNbtaggedJets_ttjets->Write();
  //hNbtaggedJets_wjets->Draw("COLZ");
  hNbtaggedJets_wjets->Write();
  //hNbtaggedJets_zjets->Draw("COLZ");
  hNbtaggedJets_zjets->Write();
  //hNbtaggedJets_stjets->Draw("COLZ");
  hNbtaggedJets_stjets->Write();
  //hNbtaggedJets_vbjets->Draw("COLZ");
  //hNbtaggedJets_vbjets->Write();

  hNbtaggedJets->Write();
  
  hWbb_In_WJets->Write();
  hWcx_In_WJets->Write();
  hWLFg_In_WJets->Write();

if(doPE){
  RooPlot* frame_Nttpull = mcstudy->plotPull(Ntt,-4,4,100) ;
  frame_Nttpull->SetName("myRooPlot_Nttpull");
  frame_Nttpull->Write();
  RooPlot* frame_Ntt = mcstudy->plotParam(Ntt);//,Binning(100,NbOfTTlike[NbOfJets-3]*(1-0.2),NbOfTTlike[NbOfJets-3]*(1+0.2))) ;
  frame_Ntt->SetName("myRooPlot_Ntt");
  frame_Ntt->Write();

  RooPlot* frame_Nvpull  = mcstudy->plotPull(Nv,-4,4,100) ;
  frame_Nvpull->SetName("myRooPlot_Nvpull");
  frame_Nvpull->Write();
  RooPlot* frame_Nv = mcstudy->plotParam(Nv);//,Binning(100,NbOfVlike[NbOfJets-3]*(1-0.2),NbOfVlike[NbOfJets-3]*(1+0.2))) ;
  frame_Nv->SetName("myRooPlot_Nv");
  frame_Nv->Write();

  if(!eb.isConstant()){
	RooPlot* frame_ebpull = mcstudy->plotPull(eb,-4,4,100) ;
	frame_ebpull->SetName("myRooPlot_Ebpull");
	frame_ebpull->Write();
  	RooPlot* frame_eb = mcstudy->plotParam(eb,Binning(100)) ;
  	frame_eb->SetName("myRooPlot_eb");
  	frame_eb->Write();
  }
  if(!eudsc.isConstant()){
	RooPlot* frame_eudscpull = mcstudy->plotPull(eudsc,-4,4,100) ;
	frame_eudscpull->SetName("myRooPlot_Eudscpull");
	frame_eudscpull->Write();
  	RooPlot* frame_eudsc = mcstudy->plotParam(eudsc,Binning(100)) ;
  	frame_eudsc->SetName("myRooPlot_eudsc");
  	frame_eudsc->Write();
  }
  if(!euds.isConstant()){
	RooPlot* frame_eudspull  = mcstudy->plotPull(euds,-4,4,100) ;
	frame_eudspull->SetName("myRooPlot_Eudspull");
	frame_eudspull->Write();
  	RooPlot* frame_euds = mcstudy->plotParam(euds,Binning(100)) ;
  	frame_euds->SetName("myRooPlot_euds");
  	frame_euds->Write();
  }
}
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
 gcc -O $(root-config --cflags) $(root-config --libs)  -lMathMore -lFoam -lMinuit -lRooFitCore -lRooFit -L/sandbox/cmss/slc5_amd64_gcc434/cms/cmssw/CMSSW_4_2_8/external/slc5_amd64_gcc434/lib -ldcap -o StopSearches_VJetsBckgdEst.exe StopSearches_VJetsBckgdEst.cc
 */

int main (int argc, char *argv[])
{
  return StopSearches_VJetsBckgdEst();
}
