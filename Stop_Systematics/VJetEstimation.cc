#include "../interface/VJetEstimation.h"

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
	Double_t Ntt_bjets = VJetEstimation::Ntt_bjets(Ntt,eb,eudsc,njets, nbjets);
	Double_t Nv_bjets  = VJetEstimation::Nv_bjets( Nv, euds, njets, nbjets);
	Double_t Nvb_bjets = VJetEstimation::Nvb_bjets(Nvb,eb,euds, njets, nbjets);
	return (Ntt_bjets+Nv_bjets+Nvb_bjets);
}

Double_t VJetEstimation::Ntt_bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t njets, Int_t nbjets) {
	if(nbjets == 0)      return VJetEstimation::Ntt_0bjet( Ntt, eb, eudsc, njets);
	else if(nbjets == 1) return VJetEstimation::Ntt_1bjet( Ntt, eb, eudsc, njets);
	else if(nbjets == 2) return VJetEstimation::Ntt_2bjets(Ntt, eb, eudsc, njets);
	else if(nbjets == 3) return VJetEstimation::Ntt_3bjets(Ntt, eb, eudsc, njets);
	else                 return -9999;
}

Double_t VJetEstimation::Ntt_err_bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb, Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t njets, Int_t nbjets) {
	if(nbjets == 0)      return VJetEstimation::Ntt_err_0bjet( Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 1) return VJetEstimation::Ntt_err_1bjet( Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 2) return VJetEstimation::Ntt_err_2bjets(Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
	else if(nbjets == 3) return VJetEstimation::Ntt_err_3bjets(Ntt, Ntt_err, eb, eb_err, eudsc, eudsc_err, njets);
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
	Double_t Ntt_0bjet = VJetEstimation::Ntt_0bjet(Ntt,eb,eudsc,n);
	Double_t Nv_0bjet  = VJetEstimation::Nv_0bjet( Nv, euds, n);
	Double_t Nvb_0bjet = VJetEstimation::Nvb_0bjet(Nvb,eb,euds, n);
	return (Ntt_0bjet+Nv_0bjet+Nvb_0bjet);
}

Double_t VJetEstimation::Ntt_0bjet(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) {
	Double_t  Ntt_0bjet = ((1-eb)*(1-eb)*pow((1-eudsc),n-2)*GetTTEff2bq(n)
                         +         (1-eb)*pow((1-eudsc),n-1)*GetTTEff1bq(n)
                         +                pow((1-eudsc),n)  *GetTTEff0bq(n))*Ntt;
	return (Ntt_0bjet < 0 ? 0 : Ntt_0bjet);
}
Double_t VJetEstimation::Ntt_err_0bjet(Double_t Ntt, Double_t Ntt_err, Double_t eb, Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) {
	Double_t  Ntt_err_0bjet = pow((2*(eb-1)*pow((1-eudsc),n-2)*GetTTEff2bq(n)
                                 +    (-1)*pow((1-eudsc),n-1)*GetTTEff1bq(n))*Ntt*eb_err,2)
  + pow((pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)*(1-eb)*(1-eb)*GetTTEff2bq(n)
         +pow(-1.,n-1)*(n-1)*pow(eudsc-1,n-2)*       (1-eb)*GetTTEff1bq(n)
         +pow(-1.,n  )*(n  )*pow(eudsc-1,n-1)*       (1   )*GetTTEff0bq(n))*Ntt*eudsc_err,2)
  + pow(Ntt_0bjet(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_0bjet < 0 ? 0 : sqrt(Ntt_err_0bjet));
}
Double_t VJetEstimation::Nvb_0bjet(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) {
	Double_t  Nvb_0bjet = (1-eb)*pow((1-euds),n-1)*Nvb;
	return (Nvb_0bjet < 0 ? 0 : Nvb_0bjet);
}
Double_t VJetEstimation::Nv_0bjet (Double_t Nv, Double_t euds, Int_t n) {
	Double_t  Nv_0bjet = pow((1-euds),n)*Nv;
	return (Nv_0bjet < 0 ? 0 : Nv_0bjet);
}
Double_t VJetEstimation::Nv_err_0bjet (Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) {
	Double_t  Nv_err_0bjet = pow(pow(-1.,n)*n*pow(euds-1,n-1)*Nv*euds_err,2)
  + pow(pow((1-euds),n)*Nv_err,2);
	return (Nv_err_0bjet < 0 ? 0 : sqrt(Nv_err_0bjet));
}

Double_t VJetEstimation::N1bjet   (Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc, Double_t euds, Int_t n) {
	Double_t Ntt_1bjet = VJetEstimation::Ntt_1bjet(Ntt,eb,eudsc,n);
	Double_t Nv_1bjet  = VJetEstimation::Nv_1bjet( Nv, euds, n);
	Double_t Nvb_1bjet = VJetEstimation::Nvb_1bjet(Nvb,eb,euds, n);
	return (Ntt_1bjet+Nv_1bjet+Nvb_1bjet);
}
Double_t VJetEstimation::Ntt_1bjet(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) {
	Double_t  Ntt_1bjet =((2*eb*(1-eb)*pow(1-eudsc,n-2)+(1-eb)*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3))*GetTTEff2bq(n)
                        +          (eb*pow(1-eudsc,n-1)+       (1-eb)*(n-1)*eudsc*pow(1-eudsc,n-2))*GetTTEff1bq(n)
                        +                                                (n*eudsc*pow(1-eudsc,n-1))*GetTTEff0bq(n))*Ntt;
	return (Ntt_1bjet < 0 ? 0 : Ntt_1bjet);
}
Double_t VJetEstimation::Ntt_err_1bjet(Double_t Ntt, Double_t Ntt_err, Double_t eb,  Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) {
	Double_t  Ntt_err_1bjet = pow(((2*(1-2*eb)*pow(1-eudsc,n-2)+2*(eb-1)*(n-2)*eudsc*pow(1-eudsc,n-3))*GetTTEff2bq(n)
                                 +           (pow(1-eudsc,n-1)+    (-1)*(n-1)*eudsc*pow(1-eudsc,n-2))*GetTTEff1bq(n))*Ntt*eb_err,2)
  + pow(((2*eb*(1-eb)*(n-2)*pow(-1.,n-2)*pow(eudsc-1,n-3)+2*(eb-1)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(1-eudsc,n-4)))*GetTTEff2bq(n)
         +(         eb*(n-1)*pow(-1.,n-1)*pow(eudsc-1,n-2)+  (1-eb)*(n-1)*(pow(1-eudsc,n-2)+eudsc*pow(-1.,n-2)*(n-2)*pow(1-eudsc,n-3)))*GetTTEff1bq(n)
         +(                                                        (n  )*(pow(1-eudsc,n-1)+eudsc*pow(-1.,n-1)*(n-1)*pow(1-eudsc,n-2)))*GetTTEff0bq(n))*Ntt*eudsc_err,2)
  + pow(Ntt_1bjet(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_1bjet < 0 ? 0 : sqrt(Ntt_err_1bjet));
}
Double_t VJetEstimation::Nvb_1bjet(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) {
	Double_t  Nvb_1bjet = (eb*pow(1-euds,n-1)+(1-eb)*(n-1)*euds*pow(1-euds,n-2))*Nvb;
	return (Nvb_1bjet < 0 ? 0 : Nvb_1bjet);
}
Double_t VJetEstimation::Nv_1bjet(Double_t Nv, Double_t euds, Int_t n) {
	Double_t  Nv_1bjet = n*euds*pow(1-euds,n-1)*Nv;
	return (Nv_1bjet < 0 ? 0 : Nv_1bjet);
}
Double_t VJetEstimation::Nv_err_1bjet(Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) {
	Double_t  Nv_err_1bjet = pow(n*(pow(1-euds,n-1)+euds*pow(-1.,n-1)*(n-1)*pow(euds-1,n-2))*Nv*euds_err,2)
  + pow(Nv_1bjet(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_1bjet < 0 ? 0 : sqrt(Nv_err_1bjet));
}

Double_t VJetEstimation::N2bjets(Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) {
	Double_t  Ntt_2bjets = VJetEstimation::Ntt_2bjets(Ntt,eb,eudsc,n);
	Double_t  Nv_2bjets  = VJetEstimation::Nv_2bjets( Nv ,euds, n);
	Double_t  Nvb_2bjets = VJetEstimation::Nvb_2bjets(Nvb,eb,euds, n);
	return (Ntt_2bjets+Nv_2bjets+Nvb_2bjets);
}
Double_t VJetEstimation::Ntt_2bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) {
	Double_t  Ntt_2bjets =((eb*eb*pow(1-eudsc,n-2)+2*eb*(1-eb)*(n-2)*eudsc*pow(1-eudsc,n-3)+(1-eb)*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4))*GetTTEff2bq(n)
                         +                                 (eb*(n-1)*eudsc*pow(1-eudsc,n-2)+       (1-eb)*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*GetTTEff1bq(n)
                         +                                                                                   ((n*(n-1)/2)*eudsc*eudsc*pow(1-eudsc,n-2))*GetTTEff0bq(n))*Ntt;
	return (Ntt_2bjets<0 ? 0 : Ntt_2bjets);
}
Double_t VJetEstimation::Ntt_err_2bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb,  Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) {
	Double_t  Ntt_err_2bjets  = pow(((2*eb*pow(1-eudsc,n-2)+2*(1-2*eb)*(n-2)*eudsc*pow(1-eudsc,n-3)+2*(eb-1)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4))*GetTTEff2bq(n)
                                   +                                 ((n-1)*eudsc*pow(1-eudsc,n-2)+    (-1)*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3))*GetTTEff1bq(n))*Ntt*eb_err,2)
  + pow(((pow(eb,2)*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)+2*eb*(1-eb)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+pow(1-eb,2)*((n-2)*(n-3)/2)*(2*eudsc*pow(1-eudsc,n-4)+pow(eudsc,2)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5)))*GetTTEff2bq(n)
         +                                                      (eb*(n-1)*(pow(1-eudsc,n-2)+eudsc*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3))+     (1-eb)*((n-1)*(n-2)/2)*(2*eudsc*pow(1-eudsc,n-3)+pow(eudsc,2)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))))*GetTTEff1bq(n)
        +                                                                                                                 			(((n)*(n-1)/2)*(2*eudsc*pow(1-eudsc,n-2)+pow(eudsc,2)*pow(-1.,n-2)*(n-2)*pow(eudsc-1,n-3)))*GetTTEff0bq(n)*Ntt*eudsc_err,2)
  + pow(Ntt_2bjets(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_2bjets < 0 ? 0 : sqrt(Ntt_err_2bjets));
}

Double_t VJetEstimation::Nvb_2bjets(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) {
	Double_t  Nvb_2bjets = (eb*(n-1)*euds*pow(1-euds,n-2)+(1-eb)*((n-1)*(n-2)/2)*euds*euds*pow(1-euds,n-3))*Nvb;
	return (Nvb_2bjets<0 ? 0 : Nvb_2bjets);
}
Double_t VJetEstimation::Nv_2bjets (Double_t Nv, Double_t euds, Int_t n) {
	Double_t  Nv_2bjets  = ((n*(n-1)/2)*euds*euds*pow((1-euds),n-2))*Nv;
	return (Nv_2bjets<0 ? 0 : Nv_2bjets);
}
Double_t VJetEstimation::Nv_err_2bjets(Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) {
	Double_t  Nv_err_2bjets = pow((n*(n-1)/2)*(2*euds*pow(1-euds,n-2)+pow(euds,2)*pow(-1.,n-2)*(n-2)*pow(euds-1,n-3))*Nv*euds_err,2)
  + pow(Nv_2bjets(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_2bjets < 0 ? 0 : sqrt(Nv_err_2bjets));
}

Double_t VJetEstimation::N3bjets   (Double_t Ntt, Double_t Nv, Double_t Nvb, Double_t eb, Double_t eudsc,  Double_t euds, Int_t n) {
  Double_t Ntt_3bjets = VJetEstimation::Ntt_3bjets(Ntt, eb, eudsc,n);
	Double_t Nv_3bjets  = VJetEstimation::Nv_3bjets( Nv , euds, n);
	Double_t Nvb_3bjets = VJetEstimation::Nvb_3bjets(Nvb, eb, euds, n);
	return (Ntt_3bjets+Nv_3bjets+Nvb_3bjets);
}
Double_t VJetEstimation::Ntt_3bjets(Double_t Ntt, Double_t eb, Double_t eudsc, Int_t n) {
	Double_t  Ntt_3bjets =((eb*eb*(n-2)*eudsc*pow(1-eudsc,n-3)+2*eb*(1-eb)*((n-2)*(n-3)/2)*eudsc*eudsc*pow(1-eudsc,n-4)+(n>4 ? pow((1-eb),2)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow((1-eudsc),n-5) : 0 ))*GetTTEff2bq(n)
                         +                                              (eb*((n-1)*(n-2)/2)*eudsc*eudsc*pow(1-eudsc,n-3)+ (1-eb)*((n-1)*(n-2)*(n-3)/6)*pow(eudsc,3)*pow(1-eudsc,n-4))*GetTTEff1bq(n)
                         +                                                                                                          ((n*(n-1)*(n-2)/6)*pow(eudsc,3)*pow(1-eudsc,n-3))*GetTTEff0bq(n))*Ntt;
	return (Ntt_3bjets<0 ? 0 : Ntt_3bjets);
}
Double_t VJetEstimation::Ntt_err_3bjets(Double_t Ntt, Double_t Ntt_err, Double_t eb,  Double_t eb_err, Double_t eudsc, Double_t eudsc_err, Int_t n) {
	Double_t  Ntt_err_3bjets  = pow(((2*eb*(n-2)*eudsc*pow(1-eudsc,n-3)+2*(1-2*eb)*(n-2)*(n-3)/2*pow(eudsc,2)*pow(1-eudsc,n-4)+2*(eb-1)*((n-2)*(n-3)*(n-4)/6)*pow(eudsc,3)*pow(1-eudsc,n-5))*GetTTEff2bq(n)
                                   +                                             ((n-1)*(n-2)/2*pow(eudsc,2)*pow(1-eudsc,n-3)+    (-1)*((n-1)*(n-2)*(n-3)/6)*pow(eudsc,3)*pow(1-eudsc,n-4))*GetTTEff1bq(n))*Ntt*eb_err,2)
  + pow(((pow(eb,2)*(n-2)*(pow(1-eudsc,n-3)+eudsc*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+2*eb*(1-eb)*((n-2)*(n-3)/2)*(2*eudsc*pow(1-eudsc,n-4)+pow(eudsc,2)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5))+pow(1-eb,2)*((n-2)*(n-3)*(n-4)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-5)+pow(eudsc,3)*pow(-1.,n-5)*(n-5)*pow(eudsc-1,n-6)))*GetTTEff2bq(n)
         +                                                     			            (eb*((n-1)*(n-2)/2)*(2*eudsc*pow(1-eudsc,n-3)+pow(eudsc,2)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4))+     (1-eb)*((n-1)*(n-2)*(n-3)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-4)+pow(eudsc,3)*pow(-1.,n-4)*(n-4)*pow(eudsc-1,n-5))))*GetTTEff1bq(n)
        +                                                                                                                 				      			      		        (((n)*(n-1)*(n-2)/6)*(3*pow(eudsc,2)*pow(1-eudsc,n-3)+pow(eudsc,3)*pow(-1.,n-3)*(n-3)*pow(eudsc-1,n-4)))*GetTTEff0bq(n)*Ntt*eudsc_err,2)
  + pow(Ntt_3bjets(Ntt,eb,eudsc,n)*Ntt_err/Ntt,2);
	return (Ntt_err_3bjets < 0 ? 0 : sqrt(Ntt_err_3bjets));
}
Double_t VJetEstimation::Nvb_3bjets(Double_t Nvb, Double_t eb, Double_t euds, Int_t n) {
	Double_t  Nvb_3bjets = (eb*((n-1)*(n-2)/2)*euds*euds*pow(1-euds,n-3) + (1-eb)*((n-1)*(n-2)*(n-3)/6)*pow(euds,3)*pow(1-euds,n-4))*Nvb;
	return (Nvb_3bjets<0 ? 0 : Nvb_3bjets);
}
Double_t VJetEstimation::Nv_3bjets (Double_t Nv, Double_t euds, Int_t n) {
	Double_t  Nv_3bjets  = ((n*(n-1)*(n-2)/6)*pow((euds),3)*pow((1-euds),n-3))*Nv;
	return (Nv_3bjets<0 ? 0 : Nv_3bjets);
}
Double_t VJetEstimation::Nv_err_3bjets(Double_t Nv, Double_t Nv_err, Double_t euds, Double_t euds_err, Int_t n) {
	Double_t  Nv_err_3bjets = pow((n*(n-1)*(n-2)/6)*(3*pow(euds,2)*pow(1-euds,n-3)+pow(euds,3)*pow(-1.,n-3)*(n-3)*pow(euds-1,n-4))*Nv*euds_err,2)
  + pow(Nv_3bjets(Nv,euds,n)*Nv_err/Nv,2);
	return (Nv_err_3bjets < 0 ? 0 : sqrt(Nv_err_3bjets));
}
