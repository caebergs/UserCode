#ifndef VJetEstimation_h
#define VJetEstimation_h

#include <iostream>
#include <iomanip>
#include <vector>

#include "Math/VectorUtil.h"
//#include "Math/Polynomial.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TList.h"
#include "TObject.h"

#include "TROOT.h"
/*
#include "TopTreeAnalysis/Content/interface/MCObsExpectation.h"
#include "TopTreeProducer/interface/TRootJet.h"
#include "TopTreeAnalysis/Content/interface/Dataset.h"
*/
  // RooFit librairies
#include "RooArgSet.h"
#include "RooAddition.h"
#include "RooCategory.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"

#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
  //#include "RooSimultaneous.h"

#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
  //#include "RooNLLVar.h"
  //#include "RooProfileLL.h"
#include "RooPlot.h"

#include "RooFitResult.h"

using namespace std;
  //using namespace TopTree;
using namespace RooFit ;

/**
 
 Usage :
 Constructor
 fill
 fillsummaryhistos
 ComputeEffFromMC
 ...
 SetInitialValues_...
 Unbinned...
 */
class VJetEstimation : public TObject
{
  
public:
  VJetEstimation();
	/** Method to book histograms */
      //////////////////////////////////////
    // Access to the class members
    //////////////////////////////////////
  
	
  
    // private functions
public:
	
    //Important:
    //in the following, njets denotes the index of the Njets_ array :
  
    ///////////////////////////////////////////////////////////////////	
    // Methods to calculates the estimators ...  (from equations)
    ///////////////////////////////////////////////////////////////////
  Double_t Nbjets   (const Double_t Ntt, const Double_t Nv, const Double_t Nvb, const Double_t eb, const Double_t eudsc,  Double_t euds, Int_t njets, Int_t nbjets);
  Double_t Ntt_bjets(const Double_t Ntt, const Double_t eb, const Double_t eudsc, const Int_t njets, const Int_t nbjets);
  static Double_t Nvb_bjets(const Double_t Nvb, const Double_t eb, const Double_t euds,  Int_t njets, Int_t nbjets);
  static Double_t Nv_bjets( Double_t Nv,  Double_t euds, Int_t njets, Int_t nbjets);
  Double_t Ntt_err_bjets(const Double_t Ntt, const Double_t Ntt_err, const Double_t eb, const Double_t eb_err, const Double_t eudsc, const Double_t eudsc_err, const Int_t njets, const Int_t nbjets);
  static Double_t Nv_err_bjets( Double_t Nv,  Double_t Nv_err,  Double_t euds, const Double_t euds_err, Int_t njets, Int_t nbjets);
  
  Double_t N0bjet   (const Double_t Ntt, const Double_t Nv, const Double_t Nvb, const Double_t eb, const Double_t eudsc,  Double_t euds, Int_t n);
  Double_t Ntt_0bjet(const Double_t Ntt, const Double_t eb, const Double_t eudsc, const Int_t n);
  Double_t Ntt_0bjet_e3(const Double_t Ntt, const Double_t eb, const Double_t eudsc, const Int_t n);
  static Double_t Nvb_0bjet(const Double_t Nvb, const Double_t eb, const Double_t euds,  Int_t n);
  static Double_t Nv_0bjet (const Double_t Nv,  Double_t euds,  Int_t n);
  Double_t Ntt_err_0bjet(const Double_t Ntt, const Double_t Ntt_err, const Double_t eb,  const  Double_t eb_err,   const Double_t eudsc, const Double_t eudsc_err, const Int_t n);
  static Double_t Nv_err_0bjet (const Double_t Nv,  Double_t Nv_err,  Double_t euds, const Double_t euds_err, Int_t n);
	
  Double_t N1bjet   (const Double_t Ntt, const Double_t Nv, const Double_t Nvb, const Double_t eb, const Double_t eudsc,  Double_t euds, Int_t n);
  Double_t Ntt_1bjet(const Double_t Ntt, const Double_t eb, const Double_t eudsc, const Int_t n);
  Double_t Ntt_1bjet_e3(const Double_t Ntt, const Double_t eb, const Double_t eudsc, const Int_t n);
  static Double_t Nvb_1bjet(const Double_t Nvb, const Double_t eb, const Double_t euds,  Int_t n);
  static Double_t Nv_1bjet (const Double_t Nv,  Double_t euds,  Int_t n);
  Double_t Ntt_err_1bjet(const Double_t Ntt, const Double_t Ntt_err, const Double_t eb,   const Double_t eb_err, const Double_t eudsc, const Double_t eudsc_err, const Int_t n);
  static Double_t Nv_err_1bjet (const Double_t Nv,  Double_t Nv_err,  Double_t euds, const Double_t euds_err, Int_t n);
	
  Double_t N2bjets   (const Double_t Ntt, const Double_t Nv, const Double_t Nvb, const Double_t eb, const Double_t eudsc,  Double_t euds, Int_t n);
  Double_t Ntt_2bjets(const Double_t Ntt, const Double_t eb, const Double_t eudsc, const Int_t n);
  Double_t Ntt_2bjets_e3(const Double_t Ntt, const Double_t eb, const Double_t eudsc, const Int_t n);
  static Double_t Nvb_2bjets(const Double_t Nvb, const Double_t eb, const Double_t euds,  Int_t n);
  static Double_t Nv_2bjets (const Double_t Nv,  Double_t euds,  Int_t n);
  Double_t Ntt_err_2bjets(const Double_t Ntt, const Double_t Ntt_err, const Double_t eb,   Double_t eb_err, const Double_t eudsc, const Double_t eudsc_err, Int_t n);
  static Double_t Nv_err_2bjets (const Double_t Nv,  Double_t Nv_err,  Double_t euds, const Double_t euds_err, Int_t n);
	
  Double_t N3bjets   (const Double_t Ntt, const Double_t Nv, const Double_t Nvb, const Double_t eb, const Double_t eudsc,  Double_t euds, Int_t n);
  Double_t Ntt_3bjets(const Double_t Ntt, const Double_t eb, const Double_t eudsc, Int_t n);
  Double_t Ntt_3bjets_e3(const Double_t Ntt, const Double_t eb, const Double_t eudsc, Int_t n);
  static Double_t Nvb_3bjets(const Double_t Nvb, const Double_t eb, const Double_t euds,  Int_t n);
  static Double_t Nv_3bjets (const Double_t Nv,  Double_t euds,  Int_t n);
  Double_t Ntt_err_3bjets(const Double_t Ntt, const Double_t Ntt_err, const Double_t eb,   Double_t eb_err, const Double_t eudsc, const Double_t eudsc_err, Int_t n);
  static Double_t Nv_err_3bjets (const Double_t Nv,  Double_t Nv_err,  Double_t euds, const Double_t euds_err, Int_t n);

  
  static Double_t probElemWbb(Double_t e_b, Double_t e_uds, Int_t Nbtags, Int_t n);
  static Double_t probElemWbb(Double_t e_b, Double_t e_uds, Int_t Nb, Int_t Nbtags, Int_t n);
  static Double_t biComb(UInt_t in1, UInt_t outof1, UInt_t in2, UInt_t outof2);
  static Double_t biCombD(UInt_t in1, UInt_t outof1, UInt_t in2, UInt_t outof2);

  
  void SetTTEffbq(std::vector<Double_t> Eff0bq, std::vector<Double_t> Eff1bq, std::vector<Double_t> Eff2bq, std::vector<Double_t> Eff3bq=std::vector<Double_t>(15, 0.));

  Double_t GetTTEff0bq(Int_t n) { return Eff0bq_[n]; };
  Double_t GetTTEff1bq(Int_t n) { return Eff1bq_[n]; };
  Double_t GetTTEff2bq(Int_t n) { return Eff2bq_[n]; };
  Double_t GetTTEff3bq(Int_t n) { return Eff3bq_[n]; };

private:

std::vector<Double_t> Eff0bq_;
std::vector<Double_t> Eff1bq_;
std::vector<Double_t> Eff2bq_;
std::vector<Double_t> Eff3bq_;
/*
   //SemiMuon
  static Double_t GetTTEff0bq(Int_t n) { Double_t a[7]={0,0,0,0, 0.022348, 0.016215, 0.013229 }; return a[n]; };
  static Double_t GetTTEff1bq(Int_t n) { Double_t a[7]={0,0,0,0, 0.283191, 0.238223, 0.211845 }; return a[n]; };
  static Double_t GetTTEff2bq(Int_t n) { Double_t a[7]={0,0,0,0, 0.686278, 0.725439, 0.732775 }; return a[n]; };
*/
/*
   //SemiElectron
  static Double_t GetTTEff0bq(Int_t n) { Double_t a[7]={0,0,0,0, 0.023797, 0.017614, 0.013601 }; return a[n]; };
  static Double_t GetTTEff1bq(Int_t n) { Double_t a[7]={0,0,0,0, 0.288904, 0.242580, 0.217566 }; return a[n]; };
  static Double_t GetTTEff2bq(Int_t n) { Double_t a[7]={0,0,0,0, 0.679199, 0.719735, 0.727453 }; return a[n]; };
*/
/*
   //SemiMuonSemiElectron
  static Double_t GetTTEff0bq(Int_t n) { Double_t a[7]={0,0,0,0, 0.022928, 0.016823, 0.013386 }; return a[n]; };
  static Double_t GetTTEff1bq(Int_t n) { Double_t a[7]={0,0,0,0, 0.284861, 0.240341, 0.214254 }; return a[n]; };
  static Double_t GetTTEff2bq(Int_t n) { Double_t a[7]={0,0,0,0, 0.684076, 0.722721, 0.730549 }; return a[n]; };
*/
//	ClassDef (VJetEstimation,1);
};

#endif
