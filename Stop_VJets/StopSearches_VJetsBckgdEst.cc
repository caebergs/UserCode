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
#include "TEventList.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TMarker.h"
#include "TPaveStats.h"

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

//int main (int argc, char *argv[])
int StopSearches_VJetsBckgdEst()
{

  clock_t start = clock();

  cout << "***************************************************************" << endl;
  cout << " Beginning of the program for the V+Jets background estimation " << endl;
  cout << "***************************************************************" << endl;

  bool verbose = true;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////Configuration ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /** TO BE ADAPTED **/
  TChain *ttjets = new TChain("AnaTree");
  ttjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/big_TTbar_skim.root");
  TChain *wjets = new TChain("AnaTree");
  wjets->Add("dcap:///pnfs/iihe/cms/store/user/nstrobbe/Stop_443/muon/4J_0B/skims/WJets_skim.root");
  
  if(verbose){
    cout<<"Number of events for :"<<endl;
    cout<<" - ttjets ="<<ttjets->GetEntries()<<endl;
    cout<<" -  wjets ="<< wjets->GetEntries()<<endl;
  }

  // Processes XS (pb)
  float XS_ttjets = 165.;
  float XS_wjets = 31314.;
  
  // Retrieve int. lumi from InputTree (?)
  const double IntLumi = 5000;
  
  // Normalization factors for a given int. lumi.
  /** TO BE ADAPTED **/
  bool skimmedTree = true;
  double NormFact_ttjets = 0;
  double NormFact_wjets = 0;
  
  if(skimmedTree) {
    ttjets->SetBranchAddress("eventWeight",&NormFact_ttjets);
    ttjets->GetEntry(0);
    NormFact_ttjets *= IntLumi;
    wjets->SetBranchAddress("eventWeight",&NormFact_wjets);
    wjets->GetEntry(0);
    NormFact_wjets *= IntLumi;
  }
  else{
    NormFact_ttjets = IntLumi*XS_ttjets/ttjets->GetEntries();
    NormFact_wjets  = IntLumi*XS_wjets/wjets->GetEntries();
  }

  //Output ROOT file
  string postfix = "_VJetsEstimation";
  string channelpostfix = "_SemiMuon";
  string comment = "_Constants_euds";
  string rootFileName ("StopSearches"+postfix+channelpostfix+comment+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");
  TDirectory *myDir = 0;

  // E v e n t   s e l e c t i o n   c u t s
  // ---------------------------------------

  const int NbOfJetBins = 4;
  const int MinNbOfJets = 3;
  int JetIdx = 1;

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

/*
  TEntryList *entrylist_ttjets = 0;
  TEntryList *entrylist_wjets  = 0;

  vector<double>  *vjtPt_pf;
  vector<double>  *vjtPx_pf;
  vector<double>  *vjtPy_pf;
  vector<double>  *vjtPz_pf;
  vector<double>  *vjtE_pf;
  vector<double>  *vjtJECUn_pf;
  vector<double>  *vjtEta_pf;
  vector<double>  *vjtPhi_pf;
  vector<double>  *vjtDiscri_pf;
*/
/*
  TBranch        *b_vjtPt_pf;   //!
  TBranch        *b_vjtPx_pf;   //!
  TBranch        *b_vjtPy_pf;   //!
  TBranch        *b_vjtPz_pf;   //!
  TBranch        *b_vjtE_pf;   //!
  TBranch        *b_vjtJECUn_pf;   //!
  TBranch        *b_vjtEta_pf;   //!
  TBranch        *b_vjtPhi_pf;   //!
  TBranch        *b_vjtDiscri_pf;   //!
*/
/*
  fChain->SetBranchAddress("vjtPt_pf", &vjtPt_pf);//, &b_vjtPt_pf);
  fChain->SetBranchAddress("vjtPx_pf", &vjtPx_pf);//, &b_vjtPx_pf);
  fChain->SetBranchAddress("vjtPy_pf", &vjtPy_pf);//, &b_vjtPy_pf);
  fChain->SetBranchAddress("vjtPz_pf", &vjtPz_pf);//, &b_vjtPz_pf);
  fChain->SetBranchAddress("vjtE_pf", &vjtE_pf);//, &b_vjtE_pf);
  fChain->SetBranchAddress("vjtJECUn_pf", &vjtJECUn_pf);//, &b_vjtJECUn_pf);
  fChain->SetBranchAddress("vjtEta_pf", &vjtEta_pf);//, &b_vjtEta_pf);
  fChain->SetBranchAddress("vjtPhi_pf", &vjtPhi_pf);//, &b_vjtPhi_pf);
  fChain->SetBranchAddress("vjtDiscri_pf", &vjtDiscri_pf);//, &b_vjtDiscri_pf);
*/
  // S e t t i n g s    f o r   m a x i m u m   l i k e l i h o o d   e s t i m a t i o n
  // ------------------------------------------------------------------------------------

  // Probabilities to have X b-quarks in the tt-like final state **DERIVED FROM MC SIMULATIONS**
  const float e0bq_[NbOfJetBins] = {0.0247448,0.00873498,0.00638495,0.00734619};
  const float e1bq_[NbOfJetBins] = {0.379931,0.179928,0.145577,0.134986};
  const float e2bq_[NbOfJetBins] = {0.594488,0.808598,0.83759,0.836547};

  // Initial values for minimization
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


  // Histograms for the b-tagged jet multiplicity
  TH1F* hNbtaggedJets_ttjets = 0;
  TH1F* hNbtaggedJets_wjets  = 0;
  TH1F* hNbtaggedJets_ttlike = new TH1F("hNbtaggedJets_ttlike",";B-tagged jet multiplicity",10,0,10);
  TH1F* hNbtaggedJets_vlike  = new TH1F("hNbtaggedJets_vlike" ,";B-tagged jet multiplicity",10,0,10);
  TH1F* hNbtaggedJets        = new TH1F("hNbtaggedJets" ,";B-tagged jet multiplicity",10,0,10);

  vector<int>::iterator Idx;
  vector<int> FixedVarIdx;
  //FixedVarIdx.push_back(0); // Fix eb, b-tagging efficiency
  //FixedVarIdx.push_back(1); // Fix eudsc, mis-tagging efficiency for tt-like events
  FixedVarIdx.push_back(2); // Fix euds, mis-tagging efficiency for v-like events

  string wp[3]  = {"loose","medium","tight"};

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Analysis ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // E v e n t   s e l e c t i o n
  // -----------------------------
/*  
  ttjets->Draw(">>entrylist_ttjets", cuts, "entrylist");
   wjets->Draw(">>entrylist_wjets",  cuts, "entrylist");
   
  entrylist_ttjets = (TEntryList*)gDirectory->Get("entrylist_ttjets");
  entrylist_wjets  = (TEntryList*)gDirectory->Get("entrylist_wjets");
*/
  // Fill 2D histograms Nb b-tagged jets vs Nb jets
  ttjets->Draw("Njtbtag_pf>>hNbtaggedJets_ttjets(10,0,10)",cuts);
   wjets->Draw("Njtbtag_pf>>hNbtaggedJets_wjets(10,0,10)",cuts);

  // Retrieve the b-tagged jet multiplicity histograms
  hNbtaggedJets_ttjets = (TH1F*)gDirectory->Get("hNbtaggedJets_ttjets");
  hNbtaggedJets_wjets  = (TH1F*)gDirectory->Get("hNbtaggedJets_wjets");

  // Set title axis for histograms
  hNbtaggedJets_ttjets->SetTitle(";B-tagged jet multiplicity");
  hNbtaggedJets_wjets->SetTitle(";B-tagged jet multiplicity");

  // Select processes to be considered in the estimation
  // - tt-like processes
  hNbtaggedJets_ttlike->Add(hNbtaggedJets_ttjets,NormFact_ttjets);
  // - v-like processes
  hNbtaggedJets_vlike->Add(hNbtaggedJets_wjets,NormFact_wjets);

  hNbtaggedJets->Add(hNbtaggedJets_ttlike);
  hNbtaggedJets->Add(hNbtaggedJets_vlike);

  // Set true numbers of tt-like and v-like events
  float Nttlike = hNbtaggedJets_ttlike->GetEntries();
  float Nvlike = hNbtaggedJets_vlike->GetEntries();

  if(verbose){
    cout<<"***************************************"<<endl;
    cout<<"- Number of entries / NormFactors "<<endl;
    cout<<"-- for tt-like processes : "<<endl;
    cout<<"--- tt+jets = "<<hNbtaggedJets_ttjets->GetEntries()<<"/ "<<NormFact_ttjets<<endl;
    cout<<"--- total nb of events (for "<<IntLumi<<" /pb) = "<<Nttlike<<endl;
    cout<<"---------------------------------"<<endl;
    cout<<"-- for v-like processes : "<<endl;
    cout<<"--- w+jets = "<<hNbtaggedJets_wjets->GetEntries()<<"/ "<<NormFact_wjets<<endl;
    cout<<"--- total nb of events (for "<<IntLumi<<" /pb) = "<<Nvlike<<endl;
    cout<<"---------------------------------"<<endl;
    cout<<"-- for all processes : "<<endl;
    cout<<"--- total = "<<hNbtaggedJets->GetEntries()<<endl;
    cout<<"***************************************"<<endl;
  }

  // C r e a t e   o b s e r v a b l e
  // ---------------------------------

  RooRealVar Ntt("Ntt","N_{t#bar{t}-like}",Nttlike,0.0,Nttlike*3);
  RooRealVar Nv("Nv","N_{V-like}",Nvlike,0.0,Nvlike*4);

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
  
  RooAddPdf model("model","model",RooArgList(pbjets_tt,pbjets_v),RooArgList(Ntt,Nv));
	
  // Contrainsts on free parameters
  //RooGaussian Eb_constraint("Eb_constraint","Eb_constraint",eb,RooConst(btag_eff[NbOfJets-3][BtagIdx]),RooConst(btag_eff[NbOfJets-3][BtagIdx]*0.15));
  //RooGaussian Eudsc_constraint("Eudsc_constraint","Eudsc_constraint",eudsc,RooConst(mistag_eff[NbOfJets-3][BtagIdx]),RooConst(mistag_eff[NbOfJets-3][BtagIdx]*0.10));
  RooGaussian Euds_constraint("Euds_constraint","Euds_constraint",euds,RooConst(Euds_const_mean),RooConst(Euds_const_error));
  
  //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,Eb_constraint)) ;
  //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,Eudsc_constraint)) ;
  //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,Euds_constraint)) ;
  //RooProdPdf model_constraint("model_constraint","model with constraint",RooArgSet(model,Eb_constraint,Eudsc_constraint,Euds_constraint)) ;

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
/*  
  RooAbsReal* nll = model.createNLL(data);//,Optimize(0));
  
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
  RooFitResult* fit_result = model.fitTo(data,Save(1),Extended(1),Minos(1),Strategy(2));//,ExternalConstraints(RooArgSet(Euds_constraint)));
  //,ExternalConstraints(RooArgSet(Eb_constraint,Eudsc_constraint,Euds_constraint)));

	//RooFitResult* fit_result = minimizer.save();
	
	fit_result->Print("v");


  //RooMCStudy* mcstudy = new RooMCStudy(model,nbjets,Binned(kTRUE),Silence(),Extended(kTRUE),FitOptions(Optimize(0),Save(0),Minos(1),Extended(1),Strategy(2),PrintEvalErrors(-1)));

  // Generate and fit 1000 samples of Poisson(nExpected) events
  //mcstudy->generateAndFit(NbOfPE,(int)(Nttlike+Nvlike)) ;

  // P l o t   f i t   r e s u l t s 
  // ---------------------------------------------------------------------
  fout->cd();

  //hNbtaggedJets_ttjets->Draw("COLZ");
  hNbtaggedJets_ttjets->Write();
  //hNbtaggedJets_wjets->Draw("COLZ");
  hNbtaggedJets_wjets->Write();
  //hNbtaggedJets_wjets->Draw("COLZ");
  hNbtaggedJets->Write();
  
  // Closing files
  fout->Close();

  // Free memory (if needed)
  
  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}
