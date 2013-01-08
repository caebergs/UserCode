#include "VJetEstimation.h"

ClassImp(VJetEstimation)

VJetEstimation::VJetEstimation(){
}

void VJetEstimation::SetTTEffbq(std::vector<Double_t> Eff0bq, std::vector<Double_t> Eff1bq, std::vector<Double_t> Eff2bq) {
 this->Eff0bq_ = Eff0bq ;
 this->Eff1bq_ = Eff1bq ;
 this->Eff2bq_ = Eff2bq ;
}

  //////////////////////////////////////
  // Access to the parameter estimated values
  //////////////////////////////////////


  ///////////////////////////////////////////////////////////////////	
  // Methods to calculates the estimators ...  (from equations)
  ///////////////////////////////////////////////////////////////////
Double_t VJetEstimation::Nbjets(Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t njets, Int_t nbjets) {
	Double_t Ntt_bjets_ = this->Ntt_bjets(Ntt,eb,eudsc,njets, nbjets);
	Double_t Nv_bjets_  = VJetEstimation::Nv_bjets( Nv, euds, njets, nbjets);
	Double_t Nvb_bjets_ = VJetEstimation::Nvb_bjets(Nvb,eb,euds, njets, nbjets);
	return (Ntt_bjets_+Nv_bjets_+Nvb_bjets_);
}

Double_t VJetEstimation::Ntt_bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t njets, Int_t nbjets) {
	if(nbjets == 0)      return this->Ntt_0bjet( Ntt, eb, eudsc, njets);
	else if(nbjets == 1) return this->Ntt_1bjet( Ntt, eb, eudsc, njets);
	else if(nbjets == 2) return this->Ntt_2bjets(Ntt, eb, eudsc, njets);
	else if(nbjets == 3) return this->Ntt_3bjets(Ntt, eb, eudsc, njets);
	else                 return -9999;
}

Double_t VJetEstimation::Ntt_err_bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb, Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t njets, Int_t nbjets) {
	if(nbjets == 0)      return this->Ntt_err_0bjet( Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 1) return this->Ntt_err_1bjet( Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 2) return this->Ntt_err_2bjets(Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 3) return this->Ntt_err_3bjets(Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else                 return -9999;
}

Double_t VJetEstimation::Nvb_bjets(Double_t Nvb, Double_t eb, Double_t euds, Int_t njets, Int_t nbjets) {
	if(nbjets == 0)      return VJetEstimation::Nvb_0bjet( Nvb, eb, euds, njets);
	else if(nbjets == 1) return VJetEstimation::Nvb_1bjet( Nvb, eb, euds, njets);
	else if(nbjets == 2) return VJetEstimation::Nvb_2bjets(Nvb, eb, euds, njets);
	else if(nbjets == 3) return VJetEstimation::Nvb_3bjets(Nvb, eb, euds, njets);
	else                 return -9999;
}

Double_t VJetEstimation::Nv_bjets(Double_t Nv, Double_t euds, Int_t njets, Int_t nbjets) {
	if(nbjets == 0)      return VJetEstimation::Nv_0bjet( Nv, euds, njets);
	else if(nbjets == 1) return VJetEstimation::Nv_1bjet( Nv, euds, njets);
	else if(nbjets == 2) return VJetEstimation::Nv_2bjets(Nv, euds, njets);
	else if(nbjets == 3) return VJetEstimation::Nv_3bjets(Nv, euds, njets);
	else                 return -9999;
}

Double_t VJetEstimation::Nv_err_bjets(Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t njets, Int_t nbjets) {
	if(nbjets == 0)      return VJetEstimation::Nv_err_0bjet( Nv, Nv_err, euds, euds_err, njets);
	else if(nbjets == 1) return VJetEstimation::Nv_err_1bjet( Nv, Nv_err, euds, euds_err, njets);
	else if(nbjets == 2) return VJetEstimation::Nv_err_2bjets(Nv, Nv_err, euds, euds_err, njets);
	else if(nbjets == 3) return VJetEstimation::Nv_err_3bjets(Nv, Nv_err, euds, euds_err, njets);
	else                 return -9999;
}

Double_t VJetEstimation::N0bjet(Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) {
	Double_t Ntt_0bjet_ = this->Ntt_0bjet(Ntt,eb,eudsc,n);
	Double_t Nv_0bjet_  = VJetEstimation::Nv_0bjet( Nv, euds, n);
	Double_t Nvb_0bjet_ = VJetEstimation::Nvb_0bjet(Nvb,eb,euds, n);
	return (Ntt_0bjet_+Nv_0bjet_+Nvb_0bjet_);
}

Double_t VJetEstimation::Ntt_0bjet(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) {
	Double_t  Ntt_0bjet_ = ((1-eb)*(1-eb)*pow((1-eudsc),n-2)*GetTTEff2bq(n)
                         +         (1-eb)*pow((1-eudsc),n-1)*GetTTEff1bq(n)
                         +                pow((1-eudsc),n)  *GetTTEff0bq(n))*Ntt;
	return (Ntt_0bjet_ < 0 ? 0 : Ntt_0bjet_);
}
Double_t VJetEstimation::Ntt_err_0bjet(Double_t Ntt, Double_t Ntt_err, Double_t eb, Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) {
	Double_t  Ntt_err_0bjet_ = pow((2*(eb-1)*pow((1-eudsc),n-2)*GetTTEff2bq(n)
                                 +    (-1)*pow((1-eudsc),n-1)*GetTTEff1bq(n))*Ntt*eb_err,2)
  + pow((pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)*(1-eb)*(1-eb)*GetTTEff2bq(n)
         +pow(-1.,n-1)*(n-1)*pow(eudsc-1,n-2)*       (1-eb)*GetTTEff1bq(n)
         +pow(-1.,n  )*(n  )*pow(eudsc-1,n-1)*       (1   )*GetTTEff0bq(n))*Ntt*eudsc_err,2)
  + pow(Ntt_0bjet(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_0bjet_ < 0 ? 0 : sqrt(Ntt_err_0bjet_));
}
Double_t VJetEstimation::Nvb_0bjet(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) {
	Double_t  Nvb_0bjet_ = (1-eb)*pow((1-euds),n-1)*Nvb;
	return (Nvb_0bjet_ < 0 ? 0 : Nvb_0bjet_);
}
Double_t VJetEstimation::Nv_0bjet (Double_t Nv, Double_t euds, Int_t n) {
	Double_t  Nv_0bjet_ = pow((1-euds),n)*Nv;
	return (Nv_0bjet_ < 0 ? 0 : Nv_0bjet_);
}
Double_t VJetEstimation::Nv_err_0bjet (Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) {
	Double_t  Nv_err_0bjet_ = pow(pow(-1.,n)*n*pow(euds-1,n-1)*Nv*euds_err,2)
  + pow(pow((1-euds),n)*Nv_err,2);
	return (Nv_err_0bjet_ < 0 ? 0 : sqrt(Nv_err_0bjet_));
}

Double_t VJetEstimation::N1bjet   (Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc, Double_t euds, Int_t n) {
	Double_t Ntt_1bjet_ = this->Ntt_1bjet(Ntt,eb,eudsc,n);
	Double_t Nv_1bjet_  = VJetEstimation::Nv_1bjet( Nv, euds, n);
	Double_t Nvb_1bjet_ = VJetEstimation::Nvb_1bjet(Nvb,eb,euds, n);
	return (Ntt_1bjet_+Nv_1bjet_+Nvb_1bjet_);
}
Double_t VJetEstimation::Ntt_1bjet(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) {
	Double_t  Ntt_1bjet_ =((2*eb*(1-eb)*pow(1-eudsc,n-2)+(1-eb)*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3))*GetTTEff2bq(n)
                        +          (eb*pow(1-eudsc,n-1)+       (1-eb)*(n-1)*eudsc*pow(1-eudsc,n-2))*GetTTEff1bq(n)
                        +                                                (n*eudsc*pow(1-eudsc,n-1))*GetTTEff0bq(n))*Ntt;
	return (Ntt_1bjet_ < 0 ? 0 : Ntt_1bjet_);
}
Double_t VJetEstimation::Ntt_err_1bjet(Double_t Ntt, Double_t Ntt_err, Double_t eb,  Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) {
	Double_t  Ntt_err_1bjet_ = pow(((2*(1-2*eb)*pow(1-eudsc,n-2)+2*(eb-1)*(n-2)*eudsc*pow(1-eudsc,n-3))*GetTTEff2bq(n)
                                 +           (pow(1-eudsc,n-1)+    (-1)*(n-1)*eudsc*pow(1-eudsc,n-2))*GetTTEff1bq(n))*Ntt*eb_err,2)
  + pow(((2*eb*(1-eb)*(n-2)*pow(-1.,n-2)*pow(eudsc-1,n-3)+2*(eb-1)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(1-eudsc,n-4)))*GetTTEff2bq(n)
         +(         eb*(n-1)*pow(-1.,n-1)*pow(eudsc-1,n-2)+  (1-eb)*(n-1)*(pow(1-eudsc,n-2)+eudsc*pow(-1.,n-2)*(n-2)*pow(1-eudsc,n-3)))*GetTTEff1bq(n)
         +(                                                        (n  )*(pow(1-eudsc,n-1)+eudsc*pow(-1.,n-1)*(n-1)*pow(1-eudsc,n-2)))*GetTTEff0bq(n))*Ntt*eudsc_err,2)
  + pow(Ntt_1bjet(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_1bjet_ < 0 ? 0 : sqrt(Ntt_err_1bjet_));
}
Double_t VJetEstimation::Nvb_1bjet(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) {
	Double_t  Nvb_1bjet_ = (eb*pow(1-euds,n-1)+(1-eb)*(n-1)*euds*pow(1-euds,n-2))*Nvb;
	return (Nvb_1bjet_ < 0 ? 0 : Nvb_1bjet_);
}
Double_t VJetEstimation::Nv_1bjet(Double_t Nv, Double_t euds, Int_t n) {
	Double_t  Nv_1bjet_ = n*euds*pow(1-euds,n-1)*Nv;
	return (Nv_1bjet_ < 0 ? 0 : Nv_1bjet_);
}
Double_t VJetEstimation::Nv_err_1bjet(Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) {
	Double_t  Nv_err_1bjet_ = pow(n*(pow(1-euds,n-1)+euds*pow(-1.,n-1)*(n-1)*pow(euds-1,n-2))*Nv*euds_err,2)
  + pow(Nv_1bjet(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_1bjet_ < 0 ? 0 : sqrt(Nv_err_1bjet_));
}

Double_t VJetEstimation::N2bjets(Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) {
	Double_t  Ntt_2bjets_ = this->Ntt_2bjets(Ntt,eb,eudsc,n);
	Double_t  Nv_2bjets_  = VJetEstimation::Nv_2bjets( Nv ,euds, n);
	Double_t  Nvb_2bjets_ = VJetEstimation::Nvb_2bjets(Nvb,eb,euds, n);
	return (Ntt_2bjets_+Nv_2bjets_+Nvb_2bjets_);
}
Double_t VJetEstimation::Ntt_2bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) {
	Double_t  Ntt_2bjets_ =((eb*eb*pow(1-eudsc,n-2)+2*eb*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3)+(1-eb)*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4))*GetTTEff2bq(n)
                         +                                 (eb*(n-1)*eudsc*pow(1-eudsc,n-2)+       (1-eb)*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*GetTTEff1bq(n)
                         +                                                                                   ((n*(n-1)/2)*eudsc*eudsc*pow(1-eudsc,n-2))*GetTTEff0bq(n))*Ntt;
	return (Ntt_2bjets_<0 ? 0 : Ntt_2bjets_);
}
Double_t VJetEstimation::Ntt_err_2bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb,  Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) {
	Double_t  Ntt_err_2bjets_  = pow(((2*eb*pow(1-eudsc,n-2)+2*(1-2*eb)*(n-2)*eudsc*pow(1-eudsc,n-3)+2*(eb-1)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4))*GetTTEff2bq(n)
                                   +                                 ((n-1)*eudsc*pow(1-eudsc,n-2)+    (-1)*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*GetTTEff1bq(n))*Ntt*eb_err,2)
  + pow(((pow(eb,2)*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)+2*eb*(1-eb)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+pow(1-eb,2)*((n-2)*(n-3)/2)*(2*eudsc*pow(1-eudsc,n-4)+pow(eudsc,2)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5)))*GetTTEff2bq(n)
         +                                                      (eb*(n-1)*(pow(1-eudsc,n-2)+eudsc*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3))+     (1-eb)*((n-1)*(n-2)/2)*(2*eudsc*pow(1-eudsc,n-3)+pow(eudsc,2)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))))*GetTTEff1bq(n)
        +                                                                                                                 			(((n)*(n-1)/2)*(2*eudsc*pow(1-eudsc,n-2)+pow(eudsc,2)*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)))*GetTTEff0bq(n)*Ntt*eudsc_err,2)
  + pow(Ntt_2bjets(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_2bjets_ < 0 ? 0 : sqrt(Ntt_err_2bjets_));
}

Double_t VJetEstimation::Nvb_2bjets(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) {
	Double_t  Nvb_2bjets_ = (eb*(n-1)*euds*pow(1-euds,n-2)+(1-eb)*((n-1)*(n-2)/2)*euds*euds*pow(1-euds,n-3))*Nvb;
	return (Nvb_2bjets_<0 ? 0 : Nvb_2bjets_);
}
Double_t VJetEstimation::Nv_2bjets (Double_t Nv, Double_t euds, Int_t n) {
	Double_t  Nv_2bjets_  = ((n*(n-1)/2)*euds*euds*pow((1-euds),n-2))*Nv;
	return (Nv_2bjets_<0 ? 0 : Nv_2bjets_);
}
Double_t VJetEstimation::Nv_err_2bjets(Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) {
	Double_t  Nv_err_2bjets_ = pow((n*(n-1)/2)*(2*euds*pow(1-euds,n-2)+pow(euds,2)*pow(-1.,n-2)*(n-2)*pow(euds-1,n-3))*Nv*euds_err,2)
  + pow(Nv_2bjets(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_2bjets_ < 0 ? 0 : sqrt(Nv_err_2bjets_));
}

Double_t VJetEstimation::N3bjets   (Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) {
  Double_t Ntt_3bjets_ = this->Ntt_3bjets(Ntt, eb, eudsc,n);
	Double_t Nv_3bjets_  = VJetEstimation::Nv_3bjets( Nv , euds, n);
	Double_t Nvb_3bjets_ = VJetEstimation::Nvb_3bjets(Nvb, eb, euds, n);
	return (Ntt_3bjets_+Nv_3bjets_+Nvb_3bjets_);
}
Double_t VJetEstimation::Ntt_3bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) {
	Double_t  Ntt_3bjets_ =((eb*eb*(n-2)*eudsc*pow(1-eudsc,n-3)+2*eb*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4)+(n>4 ? pow((1-eb),2)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow((1-eudsc),n-5) : 0 ))*GetTTEff2bq(n)
                         +                                              (eb*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3)+ (1-eb)*((n-1)*(n-2)*(n-3)/6)*pow(eudsc,3)*pow(1-eudsc,n-4))*GetTTEff1bq(n)
                         +                                                                                                          ((n*(n-1)*(n-2)/6)*pow(eudsc,3)*pow(1-eudsc,n-3))*GetTTEff0bq(n))*Ntt;
	return (Ntt_3bjets_<0 ? 0 : Ntt_3bjets_);
}
Double_t VJetEstimation::Ntt_err_3bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb,  Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) {
	Double_t  Ntt_err_3bjets_  = pow(((2*eb*(n-2)*eudsc*pow(1-eudsc,n-3)+2*(1-2*eb)*(n-2)*(n-3)/2*pow(eudsc,2)*pow(1-eudsc,n-4)+2*(eb-1)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow(1-eudsc,n-5))*GetTTEff2bq(n)
                                   +                                             ((n-1)*(n-2)/2*pow(eudsc,2)*pow(1-eudsc,n-3)+    (-1)*((n-1)*(n-2)*(n-3)/6)*pow(eudsc,3)*pow(1-eudsc,n-4))*GetTTEff1bq(n))*Ntt*eb_err,2)
  + pow(((pow(eb,2)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+2*eb*(1-eb)*((n-2)*(n-3)/2)*(2*eudsc*pow(1-eudsc,n-4)+pow(eudsc,2)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5))+pow(1-eb,2)*((n-2)*(n-3)*(n-4)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-5)+pow(eudsc,3)*pow(-1.,n-5)*(n-5)*pow(eudsc-1,n-6)))*GetTTEff2bq(n)
         +                                                     			            (eb*((n-1)*(n-2)/2)*(2*eudsc*pow(1-eudsc,n-3)+pow(eudsc,2)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+     (1-eb)*((n-1)*(n-2)*(n-3)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-4)+pow(eudsc,3)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5))))*GetTTEff1bq(n)
        +                                                                                                                 				      			      		        (((n)*(n-1)*(n-2)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-3)+pow(eudsc,3)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4)))*GetTTEff0bq(n)*Ntt*eudsc_err,2)
  + pow(Ntt_3bjets(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_3bjets_ < 0 ? 0 : sqrt(Ntt_err_3bjets_));
}
Double_t VJetEstimation::Nvb_3bjets(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) {
	Double_t  Nvb_3bjets_ = (eb*((n-1)*(n-2)/2)*euds*euds*pow(1-euds,n-3) + (1-eb)*((n-1)*(n-2)*(n-3)/6)*pow(euds,3)*pow(1-euds,n-4))*Nvb;
	return (Nvb_3bjets_<0 ? 0 : Nvb_3bjets_);
}
Double_t VJetEstimation::Nv_3bjets (Double_t Nv, Double_t euds, Int_t n) {
	Double_t  Nv_3bjets_  = ((n*(n-1)*(n-2)/6)*pow((euds),3)*pow((1-euds),n-3))*Nv;
	return (Nv_3bjets_<0 ? 0 : Nv_3bjets_);
}
Double_t VJetEstimation::Nv_err_3bjets(Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) {
	Double_t  Nv_err_3bjets_ = pow((n*(n-1)*(n-2)/6)*(3*pow(euds,2)*pow(1-euds,n-3)+pow(euds,3)*pow(-1.,n-3)*(n-3)*pow(euds-1,n-4))*Nv*euds_err,2)
  + pow(Nv_3bjets(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_3bjets_ < 0 ? 0 : sqrt(Nv_err_3bjets_));
}












Double_t VJetEstimation::probElemWbb(Double_t e_b, Double_t e_uds, Int_t Nbtags, Int_t n) {
    Double_t sum = 0.;
    for (Int_t k=0 ; k<=Nbtags ; k++) {
        for (Int_t L=k ; L<=n ; L++) {
            Double_t biComb = VJetEstimation::biComb(k, L, Nbtags-k, n-L);
            if (biComb == -1.) {
                biComb = VJetEstimation::biCombD(k, L, Nbtags-k, n-L);
            }
            sum += pow(e_b,k) * pow(1.-e_b,L-k) * biComb * pow(e_uds,Nbtags-k) * pow(1.-e_uds,n-L+Nbtags-k) ;
        }
    }
    return sum;
}

Double_t VJetEstimation::probElemWbb(Double_t e_b, Double_t e_uds, Int_t Nb, Int_t Nbtags, Int_t n) {
    Double_t sum = 0.;
    for (Int_t k=0 ; k<=Nbtags ; k++) {
        Int_t L = Nb;
        if (k<=L && L<=n) {
            Double_t biComb = VJetEstimation::biComb(k, L, Nbtags-k, n-L);
            if (biComb == -1.) {
                biComb = VJetEstimation::biCombD(k, L, Nbtags-k, n-L);
            }
            sum += pow(e_b,k) * pow(1.-e_b,L-k) * biComb * pow(e_uds,Nbtags-k) * pow(1.-e_uds,n-L+Nbtags-k) ;
        }
    }
    return sum;
}

Double_t VJetEstimation::biComb(UInt_t in1, UInt_t outof1, UInt_t in2, UInt_t outof2) {
    if (in1>outof1 || in2>outof2) {
        return -1.;
    }
    ULong_t iMax1 = in1;
    if (iMax1>(outof1/2)) {
        iMax1 = outof1-in1;
    }
    ULong_t iMax2 = in2;
    if (iMax2>(outof2/2)) {
        iMax2 = outof2-in2;
    }
    ULong_t tempNum = 1;
    ULong_t tempDen = 1;
    ULong_t numerator = 1;
    ULong_t denominator = 1;
    for (UInt_t i=1 ; i<=iMax1 ; i++) {
        tempNum = numerator ;
        tempDen = denominator ;
        numerator *= (outof1-i)+1 ;
        denominator *=i ;
        if ((numerator/((outof1-i)+1) != tempNum) || (denominator/i != tempDen)) {
            return -1.;
        }
    }
    for (UInt_t i=1 ; i<=iMax2 ; i++) {
        tempNum = numerator ;
        tempDen = denominator ;
        numerator *= (outof2-i)+1 ;
        denominator *= i ;
        if ((numerator/((outof2-i)+1) != tempNum) || (denominator/i != tempDen)) {
            return -1.;
        }
    }
    return ((Double_t)numerator)/((Double_t)denominator);
}

Double_t VJetEstimation::biCombD(UInt_t in1, UInt_t outof1, UInt_t in2, UInt_t outof2) {
    if (in1>outof1 || in2>outof2) {
        return -1.;
    }
    Double_t iMax1 = in1;
    if (iMax1>(outof1/2)) {
        iMax1 = outof1-in1;
    }
    Double_t iMax2 = in2;
    if (iMax2>(outof2/2)) {
        iMax2 = outof2-in2;
    }
    Double_t tempNum = 1;
    Double_t tempDen = 1;
    Double_t numerator = 1;
    Double_t denominator = 1;
    for (UInt_t i=1 ; i<=iMax1 ; i++) {
        tempNum = numerator ;
        tempDen = denominator ;
        numerator *= (outof1-i)+1 ;
        denominator *=i ;
    }
    for (UInt_t i=1 ; i<=iMax2 ; i++) {
        tempNum = numerator ;
        tempDen = denominator ;
        numerator *= (outof2-i)+1 ;
        denominator *= i ;
    }
    return numerator/denominator;
}
