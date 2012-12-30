#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>


double ErrorOnRatio(double A, double sigmaA, double B, double sigmaB);
double ErrorOnProd( double N, double sigmaN, double R, double sigmaR);
double ErrorOnProd( double N, double sigmaN, double R, double sigmaR, double Nest, double sigmaNest);
double CorrCoef(double R, double sigmaR, double N, double sigmaN);

void TotalEstimatedNumbers_Errors_IndivChannel(string filename, int UseCase, int bin, bool UseWNJets, int channel){

  std::ifstream file;

  const int nrpoints = 15;

  double **corrMatrix = new double*[nrpoints];
  double *statUncert = new double[nrpoints];
  double *systUncert = new double[nrpoints];
  double *nbOfEvents = new double[nrpoints];

  for(int i=0; i<nrpoints; i++){
    corrMatrix[i] = new double[nrpoints];
    statUncert[i] = 0;
    for(int j=0; j<nrpoints; j++)
      corrMatrix[i][j] = 0;
  }

  string names[15] = {"$N_{st+jets}^{4jets}$","$N_{st+jets}^{5jets}$","$N_{st+jets}^{6jets}$","$\\Nttlike^{4jets}$","$\\Nttlike^{5jets}$","$\\Nttlike^{6jets}$","$\\Nvlike^{4jets}$","$\\Nvlike^{5jets}$","$\\Nvlike^{6jets}$","$N_{Vbb}^{4jets}$","$N_{Vbb}^{5jets}$","$N_{Vbb}^{6jets}$","$N_{VV}^{4jets}$","$N_{VV}^{5jets}$","$N_{VV}^{6jets}$"};
  
  {
    int i=0;
    file.open(filename.c_str());
    while (i<15){//(file.good()) {
      file >> corrMatrix[i][0] >> corrMatrix[i][1] >> corrMatrix[i][2] >> corrMatrix[i][3] >> corrMatrix[i][4] >> corrMatrix[i][5] >> corrMatrix[i][6] >> corrMatrix[i][7] >> corrMatrix[i][8] >> corrMatrix[i][9] >> corrMatrix[i][10] >> corrMatrix[i][11] >> corrMatrix[i][12] >> corrMatrix[i][13] >> corrMatrix[i][14];
      //cout<<std::fixed<<setprecision(1);
      cout<<names[i];
      for(int j=0;j<15;j++) cout<<"& $"<<corrMatrix[i][j]*100<<"$ ";
      cout<<"\\\\"<<endl;
      ++i;
    }
    file.close();
  }

  nbOfEvents[ 0] = (channel==0 ?   604 :   418);//  1030; // Nstjets 4 jets
  nbOfEvents[ 1] = (channel==0 ?   137 :   106);//   245; // Nstjets 5 jets
  nbOfEvents[ 2] = (channel==0 ?    31 :    27);//    60; // Nstjets 6 jets
  nbOfEvents[ 3] = (channel==0 ? 18409 : 12093);// 30249; // Nttlike 4 jets
  nbOfEvents[ 4] = (channel==0 ?  6854 :  4536);// 11321; // Nttlike 5 jets
  nbOfEvents[ 5] = (channel==0 ?  2437 :  1612);//  4030; // Nttlike 6 jets
  nbOfEvents[ 6] = (channel==0 ? 18555 : 14572);// 33343; // Nvlike  4 jets
  nbOfEvents[ 7] = (channel==0 ?  3016 :  2673);//  5745; // Nvlike  5 jets
  nbOfEvents[ 8] = (channel==0 ?   562 :   477);//  1052; // Nvlike  6 jets
  nbOfEvents[ 9] = (channel==0 ?  1029 :   797);//  1855; // Nvblike 4 jets
  nbOfEvents[10] = (channel==0 ?   213 :   168);//   390; // Nvblike 5 jets
  nbOfEvents[11] = (channel==0 ?    52 :    35);//    90; // Nvblike 6 jets
  nbOfEvents[12] = (channel==0 ?   253 :   216);//   470; // Nvvjets 4 jets
  nbOfEvents[13] = (channel==0 ?    40 :    34);//    74; // Nvvjets 5 jets
  nbOfEvents[14] = (channel==0 ?     7 :     5);//    12; // Nvvjets 6 jets
  
  statUncert[ 0] = (channel==0 ? 180 : 124);//  303; // Nstjets 4 jets
  statUncert[ 1] = (channel==0 ?  40 :  32);//   73; // Nstjets 5 jets
  statUncert[ 2] = (channel==0 ?  10 :   8);//   18; // Nstjets 6 jets
  statUncert[ 3] = (channel==0 ? 650 : 499);// 1010; // Nttlike 4 jets
  statUncert[ 4] = (channel==0 ? 193 : 160);//  313; // Nttlike 5 jets
  statUncert[ 5] = (channel==0 ?  72 :  59);//  103; // Nttlike 6 jets
  statUncert[ 6] = (channel==0 ? 613 : 467);//  982; // Nvlike  4 jets
  statUncert[ 7] = (channel==0 ? 168 : 137);// 266; // Nvlike  5 jets
  statUncert[ 8] = (channel==0 ?  57 :  47);//  85; // Nvlike  6 jets
  statUncert[ 9] = (channel==0 ? 290 : 226);// 502; // Nvblike 4 jets
  statUncert[10] = (channel==0 ?  62 :  51);// 111; // Nvblike 5 jets
  statUncert[11] = (channel==0 ?  16 :  10);//  25; // Nvblike 6 jets
  statUncert[12] = (channel==0 ?  77 :  65);// 141; // Nvvjets 4 jets
  statUncert[13] = (channel==0 ?  12 :  10);//  22; // Nvvjets 5 jets
  statUncert[14] = (channel==0 ?   2 :   2);//   6; // Nvvjets 6 jets

  systUncert[ 0] = (channel==0 ? 0 : 0);//  0;
  systUncert[ 1] = (channel==0 ? 0 : 0);//  0;
  systUncert[ 2] = (channel==0 ? 0 : 0);//  0;
  systUncert[ 3] = (channel==0 ? 0 : 0);//255;
  systUncert[ 4] = (channel==0 ? 0 : 0);// 65;
  systUncert[ 5] = (channel==0 ? 0 : 0);// 17;
  systUncert[ 6] = (channel==0 ? 0 : 0);//235;
  systUncert[ 7] = (channel==0 ? 0 : 0);// 66;
  systUncert[ 8] = (channel==0 ? 0 : 0);// 22;
  systUncert[ 9] = (channel==0 ? 0 : 0);//  0;
  systUncert[10] = (channel==0 ? 0 : 0);//  0;
  systUncert[11] = (channel==0 ? 0 : 0);//  0;
  systUncert[12] = (channel==0 ? 0 : 0);//  0;
  systUncert[13] = (channel==0 ? 0 : 0);//  0;
  systUncert[14] = (channel==0 ? 0 : 0);//  0;
  
  double Nttlike=0, Nttlike_err=0, Nttlike_stat_err=0, Nttlike_syst_err=0;
  double  Nvlike=0,  Nvlike_err=0,  Nvlike_stat_err=0,  Nvlike_syst_err=0;
  double Nvbjets=0, Nvbjets_err=0, Nvbjets_stat_err=0, Nvbjets_syst_err=0;
  double Nstjets=0, Nstjets_err=0, Nstjets_stat_err=0, Nstjets_syst_err=0;
  double Nvvjets=0, Nvvjets_err=0, Nvvjets_stat_err=0, Nvvjets_syst_err=0;
  double    Ntot=0,    Ntot_err=0,    Ntot_stat_err=0,    Ntot_syst_err=0;
  
  for(int i=0;i<3;i++){
    Nstjets           += nbOfEvents[i];
    Nstjets_stat_err  += pow(statUncert[i],2);
    Nstjets_syst_err  += pow(systUncert[i],2);
    Nstjets_err       += pow(statUncert[i],2) + pow(systUncert[i],2);

    Nttlike          += nbOfEvents[i+3];
    Nttlike_stat_err += pow(statUncert[i+3],2);
    Nttlike_syst_err += pow(systUncert[i+3],2);
    Nttlike_err      += pow(statUncert[i+3],2) + pow(systUncert[i+3],2);

    Nvlike           += nbOfEvents[i+6];
    Nvlike_stat_err  += pow(statUncert[i+6],2);
    Nvlike_syst_err  += pow(systUncert[i+6],2);
    Nvlike_err       += pow(statUncert[i+6],2) + pow(systUncert[i+6],2);

    Nvbjets           += nbOfEvents[i+9];
    Nvbjets_stat_err  += pow(statUncert[i+9],2);
    Nvbjets_syst_err  += pow(systUncert[i+9],2);
    Nvbjets_err       += pow(statUncert[i+9],2) + pow(systUncert[i+9],2);

    Nvvjets           += nbOfEvents[i+12];
    Nvvjets_stat_err  += pow(statUncert[i+12],2);
    Nvvjets_syst_err  += pow(systUncert[i+12],2);
    Nvvjets_err       += pow(statUncert[i+12],2) + pow(systUncert[i+12],2);

    for(int j=0;j<3;j++){
      if(i==j) continue;
        Nstjets_stat_err += corrMatrix[i][j]*statUncert[i]*statUncert[j];
        Nstjets_syst_err += corrMatrix[i][j]*systUncert[i]*systUncert[j];
        Nstjets_err      += corrMatrix[i][j]*sqrt(pow(statUncert[i],2)+pow(systUncert[i],2))*sqrt(pow(statUncert[j],2)+pow(systUncert[j],2));

        Nttlike_stat_err += corrMatrix[i+3][j+3]*statUncert[i+3]*statUncert[j+3];
        Nttlike_syst_err += corrMatrix[i+3][j+3]*systUncert[i+3]*systUncert[j+3];
        Nttlike_err      += corrMatrix[i+3][j+3]*sqrt(pow(statUncert[i+3],2)+pow(systUncert[i+3],2))*sqrt(pow(statUncert[j+3],2)+pow(systUncert[j+3],2));

        Nvlike_stat_err += corrMatrix[i+6][j+6]*statUncert[i+6]*statUncert[j+6];
        Nvlike_syst_err += corrMatrix[i+6][j+6]*systUncert[i+6]*systUncert[j+6];
        Nvlike_err      += corrMatrix[i+6][j+6]*sqrt(pow(statUncert[i+6],2)+pow(systUncert[i+6],2))*sqrt(pow(statUncert[j+6],2)+pow(systUncert[j+6],2));

        Nvbjets_stat_err += corrMatrix[i+9][j+9]*statUncert[i+9]*statUncert[j+9];
        Nvbjets_syst_err += corrMatrix[i+9][j+9]*systUncert[i+9]*systUncert[j+9];
        Nvbjets_err      += corrMatrix[i+9][j+9]*sqrt(pow(statUncert[i+9],2)+pow(systUncert[i+9],2))*sqrt(pow(statUncert[j+9],2)+pow(systUncert[j+9],2));

        Nvvjets_stat_err += corrMatrix[i+12][j+12]*statUncert[i+12]*statUncert[j+12];
        Nvvjets_syst_err += corrMatrix[i+12][j+12]*systUncert[i+12]*systUncert[j+12];
        Nvvjets_err      += corrMatrix[i+12][j+12]*sqrt(pow(statUncert[i+12],2)+pow(systUncert[i+12],2))*sqrt(pow(statUncert[j+12],2)+pow(systUncert[j+12],2));
    }
  }
  
  Nttlike_err = sqrt(Nttlike_err);
  Nttlike_stat_err = sqrt(Nttlike_stat_err);
  Nttlike_syst_err = sqrt(Nttlike_syst_err);
  cout<<"Nttlike = "<<Nttlike<<"\\pm"<<Nttlike_stat_err<<"\\pm"<<Nttlike_syst_err<<" = "<<Nttlike<<"\\pm"<<Nttlike_err<<endl;
  
  Nvlike_err  = sqrt(Nvlike_err);
  Nvlike_stat_err = sqrt(Nvlike_stat_err);
  Nvlike_syst_err = sqrt(Nvlike_syst_err);
  cout<<"Nvlike = "<<Nvlike<<"\\pm"<<Nvlike_stat_err<<"\\pm"<<Nvlike_syst_err<<" = "<<Nvlike<<"\\pm"<<Nvlike_err<<endl;

  Nvbjets_err  = sqrt(Nvbjets_err);
  Nvbjets_stat_err = sqrt(Nvbjets_stat_err);
  Nvbjets_syst_err = sqrt(Nvbjets_syst_err);
  cout<<"Nvbjets = "<<Nvbjets<<"\\pm"<<Nvbjets_stat_err<<"\\pm"<<Nvbjets_syst_err<<" = "<<Nvbjets<<"\\pm"<<Nvbjets_err<<endl;

  Nstjets_err  = sqrt(Nstjets_err);
  Nstjets_stat_err = sqrt(Nstjets_stat_err);
  Nstjets_syst_err = sqrt(Nstjets_syst_err);
  cout<<"Nstjets = "<<Nstjets<<"\\pm"<<Nstjets_stat_err<<"\\pm"<<Nstjets_syst_err<<" = "<<Nstjets<<"\\pm"<<Nstjets_err<<endl;

  Nvvjets_err  = sqrt(Nvvjets_err);
  Nvvjets_stat_err = sqrt(Nvvjets_stat_err);
  Nvvjets_syst_err = sqrt(Nvvjets_syst_err);
  cout<<"Nvvjets = "<<Nvvjets<<"\\pm"<<Nvvjets_stat_err<<"\\pm"<<Nvvjets_syst_err<<" = "<<Nvvjets<<"\\pm"<<Nvvjets_err<<endl;


  for(int i=0;i<9;i++){
    Ntot          += nbOfEvents[i];
    Ntot_stat_err += pow(statUncert[i],2);
    Ntot_syst_err += pow(systUncert[i],2);
    Ntot_err      += pow(statUncert[i],2)+pow(systUncert[i],2);

    for(int j=0;j<9;j++){
      if(i==j) continue;
      Ntot_stat_err += corrMatrix[i][j]*statUncert[i]*statUncert[j];
      Ntot_syst_err += corrMatrix[i][j]*systUncert[i]*systUncert[j];
      Ntot_err      += corrMatrix[i][j]*sqrt(pow(statUncert[i],2)+pow(systUncert[i],2))*sqrt(pow(statUncert[j],2)+pow(systUncert[j],2));
    }
  }
  
  Ntot_err      = sqrt(Ntot_err);
  Ntot_stat_err = sqrt(Ntot_stat_err);
  Ntot_syst_err = sqrt(Ntot_syst_err);
  cout<<"Nttlike+Nvlike+Nvbjets = "<<Ntot<<"\\pm"<<Ntot_stat_err<<"\\pm"<<Ntot_syst_err<<" = "<<Ntot<<"\\pm"<<Ntot_err<<endl;

  const int NbOfFiles = 15;
  string BckgdNames[NbOfFiles];
  string CChannel[2] = {"_mu","_el"};
  
  cout<<"Lepton channel : "<<(channel==0 ? CChannel[0] : CChannel[1])<<endl;
  
  std::string path = "$HOME/AnalysisCode/GentStopAnalysis/UGentCode/Leg3/v1/Samples_21102012/";
  TFile *inputfiles[NbOfFiles];
  // TT+jets
  BckgdNames[0] = "tt+jets";
  inputfiles[0] = TFile::Open((path+"TTbar"+CChannel[channel]+".root").c_str());

  // W+jets
  BckgdNames[1] = "w+jets";
  inputfiles[1] = TFile::Open((path+"Wjets"+CChannel[channel]+".root").c_str());

  // Z+jets
  BckgdNames[2] = "z+jets";
  inputfiles[2] = TFile::Open((path+"Zjets"+CChannel[channel]+".root").c_str());

  // st+jets
  BckgdNames[3] = "st+jets (s-ch)";
  inputfiles[3] = TFile::Open((path+"TopS"+CChannel[channel]+".root").c_str());
  
  BckgdNames[4] = "stbar+jets (s-ch)";
  inputfiles[4] = TFile::Open((path+"TopbarS"+CChannel[channel]+".root").c_str());
  
  BckgdNames[5] = "st+jets (t-ch)";
  inputfiles[5] = TFile::Open((path+"TopT"+CChannel[channel]+".root").c_str());
  
  BckgdNames[6] = "stbar+jets (t-ch)";
  inputfiles[6] = TFile::Open((path+"TopbarT"+CChannel[channel]+".root").c_str());
  
  BckgdNames[7] = "st+jets (tW-ch)";
  inputfiles[7] = TFile::Open((path+"TopTW"+CChannel[channel]+".root").c_str());
  
  BckgdNames[8] = "stbar+jets (tW-ch)";
  inputfiles[8] = TFile::Open((path+"TopTW"+CChannel[channel]+".root").c_str());

  // di-boson
  BckgdNames[9]  = "WW+jets";
  inputfiles[9]  = TFile::Open((path+"WW"+CChannel[channel]+".root").c_str());
  BckgdNames[10] = "WZ+jets";
  inputfiles[10] = TFile::Open((path+"WZ"+CChannel[channel]+".root").c_str());
  BckgdNames[11] = "ZZ+jets";
  inputfiles[11] = TFile::Open((path+"ZZ"+CChannel[channel]+".root").c_str());

  // W+jets exclusive
  BckgdNames[12] = "W2+jets";
  inputfiles[12] = TFile::Open((path+"W2Jets"+CChannel[channel]+".root").c_str());
  BckgdNames[13] = "W3+jets";
  inputfiles[13] = TFile::Open((path+"W3Jets"+CChannel[channel]+".root").c_str());
  BckgdNames[14] = "W4+jets";
  inputfiles[14] = TFile::Open((path+"W4Jets"+CChannel[channel]+".root").c_str());
  
  // multi-jets
  TFile *inputfiles_multijets[6];
  string BckgdNames_multijets[6];
  
  BckgdNames_multijets[0] = "QCDmu20";
  inputfiles_multijets[0] = TFile::Open((path+"QCDmu20_mu.root").c_str());
  BckgdNames_multijets[1] = "QCDel30to80";
  inputfiles_multijets[1] = TFile::Open((path+"QCDel30to80_el.root").c_str());
  BckgdNames_multijets[2] = "QCDel80to170";
  inputfiles_multijets[2] = TFile::Open((path+"QCDel80to170_el.root").c_str());
  BckgdNames_multijets[3] = "QCDel170to250";
  inputfiles_multijets[3] = TFile::Open((path+"QCDel170to250_el.root").c_str());
  BckgdNames_multijets[4] = "QCDel250to35";
  inputfiles_multijets[4] = TFile::Open((path+"QCDel250to350_el.root").c_str());
  BckgdNames_multijets[5] = "QCDel350";
  inputfiles_multijets[5] = TFile::Open((path+"QCDel350_el.root").c_str());


  const int NbOfUseCase = 4;
  // Use case : 
  // 0 : low DM
  // 1 : Int DM 1
  // 2 : Int DM 2
  // 3 : High DM
  //const int UseCase = 0;

  //string teffname[NbOfUseCase] = {"Eff_LowDM_bin", "Eff_IntDM1_bin", "Eff_IntDM2_bin", "Eff_HighDM_bin"};
  string teffname[NbOfUseCase] = {"Eff_LowDM", "Eff_IntDM1", "Eff_IntDM2", "Eff_HighDM"};
  string effname = teffname[UseCase];
  cout << "Use case : "<<teffname[UseCase]<<endl;

  //TEfficiency *tEff_ttbar, *tEff_wjets, *tEff_zjets, *tEff4_diboson;
  //TEfficiency *tEff;//[NbOfFiles];
  //int bin[NbOfUseCase] = {18,19,19,19};
  //int bin[NbOfUseCase] = {16,16,16,16};
  TEfficiency *tEff[NbOfFiles];
  TEfficiency *tDummy=0;
  TGraphAsymmErrors *tg=0;
  
  //bool UseWNJets = false;

  double x =0;
  double y[NbOfFiles]={0}, yerr_low[NbOfFiles]={0}, yerr_high[NbOfFiles]={0};
  double y_ttlike = 0, yerr_low_ttlike = 0, yerr_high_ttlike = 0;
  double y_vjets  = 0, yerr_low_vjets  = 0, yerr_high_vjets  = 0;
  double y_stjets = 0, yerr_low_stjets = 0, yerr_high_stjets = 0;
  double y_vvjets = 0, yerr_low_vvjets = 0, yerr_high_vvjets = 0;
  double y_multijets = 0, yerr_low_multijets = 0, yerr_high_multijets = 0;

  TList *pList = new TList();

  for(int i=0;i<NbOfFiles;i++){
    tEff[i] = (TEfficiency*)inputfiles[i]->Get(effname.c_str());
    tEff[i]->SetStatisticOption(TEfficiency::kBBayesian);

    y[i]         = tEff[i]->GetEfficiency(bin);
    yerr_low[i]  = tEff[i]->GetEfficiencyErrorLow(bin);
    yerr_high[i] = tEff[i]->GetEfficiencyErrorUp(bin);
    
    cout<<"R_{"<<BckgdNames[i]<<"} = "<<y[i]<<"^{+"<<yerr_high[i]<<"}_{-"<<yerr_low[i]<<"}"<<endl;
  }
  for(int i=0;i<6;i++){
    tDummy = (TEfficiency*)inputfiles_multijets[i]->Get(effname.c_str());
    tDummy->SetStatisticOption(TEfficiency::kBBayesian);
    
    cout<<"R_{"<<BckgdNames_multijets[i]<<"} = "<<tDummy->GetEfficiency(bin)<<"^{+"<<tDummy->GetEfficiencyErrorUp(bin)<<"}_{-"<<tDummy->GetEfficiencyErrorLow(bin)<<"}"<<endl;
  }

  cout<<"********************************************************"<<endl;
  cout<<"********************************************************"<<endl;

  //*******************************************************************************************	
  // tt+jets
  y_ttlike                   = y[0];
  yerr_low_ttlike            = yerr_low[0];
  yerr_high_ttlike           = yerr_high[0];

  double RTT_syst_err_up[4];
  RTT_syst_err_up[0]  = (channel==0 ? 0.000652 : 0.000936);
  RTT_syst_err_up[1]  = (channel==0 ? 0.000615 : 0.000648);
  RTT_syst_err_up[2]  = (channel==0 ? 0.000803 : 0.000923);
  RTT_syst_err_up[3]  = (channel==0 ? 0.000797 : 0.001488);
  double RTT_syst_err_low[4];
  RTT_syst_err_low[0] = (channel==0 ? 0.000825 : 0.001018);
  RTT_syst_err_low[1] = (channel==0 ? 0.000543 : 0.000657);
  RTT_syst_err_low[2] = (channel==0 ? 0.001207 : 0.001459);
  RTT_syst_err_low[3] = (channel==0 ? 0.001030 : 0.001045);

  cout<<std::fixed<<setprecision(2);
  cout<<"R_{TT} = ("<<y_ttlike*1000<<"^{+"<<yerr_high_ttlike*1000<<"}_{-"<<yerr_low_ttlike*1000<<"}^{+"<<RTT_syst_err_up[UseCase]*1000<<"}_{-"<<RTT_syst_err_low[UseCase]*1000<<"})\\cdot 10^{-3}"<<endl;

  //*******************************************************************************************
  // V+jets
  TGraphAsymmErrors *tg_vjets;
  double *w_vjets = 0;
  
  if(!UseWNJets){
    pList->Add((TEfficiency*)tEff[1]); // W+jets
    pList->Add((TEfficiency*)tEff[2]); // Z+jets
    // weights
    w_vjets = new double[2];
    w_vjets[0] = (channel==0 ? 18934 : 14970)/tEff[1]->GetTotalHistogram()->GetEntries(); // W+jets
    w_vjets[1] = (channel==0 ? 1400  :  2223)/tEff[2]->GetTotalHistogram()->GetEntries(); // Z+jets
  }
  else{
    pList->Add((TEfficiency*)tEff[12]); // W2+jets
    pList->Add((TEfficiency*)tEff[13]); // W3+jets
    pList->Add((TEfficiency*)tEff[14]); // W4+jets
    pList->Add((TEfficiency*)tEff[2]);    // Z+jets
    // weights
    w_vjets = new double[4];
    w_vjets[0] = (channel==0 ? 168.1  : 124.5  )/tEff[12]->GetTotalHistogram()->GetEntries();//1618.1;
    w_vjets[1] = (channel==0 ? 327.1  : 233.5  )/tEff[13]->GetTotalHistogram()->GetEntries();//343.2;
    w_vjets[2] = (channel==0 ? 16870. : 13643. )/tEff[14]->GetTotalHistogram()->GetEntries();//194.6;
    w_vjets[3] = (channel==0 ? 1400.  : 2223.  )/tEff[2]->GetTotalHistogram()->GetEntries();
  }
  tg_vjets = (TGraphAsymmErrors*)tDummy->Combine(pList,"mode",(UseWNJets ? 4:2),w_vjets);

  tg_vjets->GetPoint(bin-1,x,y_vjets);
  yerr_low_vjets  = tg_vjets->GetErrorYlow(bin-1);
  yerr_high_vjets = tg_vjets->GetErrorYhigh(bin-1);
  pList->Clear();
  /*******************************************************************************************************************************/
  /*******************************************************************************************************************************/
//  double RV_syst_err_up[4]  = {0.000275,0.000155,0.000262,0.000738}; // Up
//  double RV_syst_err_low[4] = {0.000117,0.000084,0.000589,0.000227}; // Down
//  double RV_syst_err_up[4]  = {0.000087,0.000083,0.000150,0.000307}; // Up
//  double RV_syst_err_low[4] = {0.000171,0.000140,0.000054,0.000000}; // Down

  double RV_syst_err_up[4]  = {0.0482035*y_vjets,0.0226508*y_vjets,0.0233062*y_vjets,0.0149049*y_vjets}; // Up
  double RV_syst_err_low[4] = {0.2622430*y_vjets,0.3797950*y_vjets,0.3956500*y_vjets,0.4722110*y_vjets}; // Down

  cout<<"R_{V} = ("<<y_vjets*1000<<"^{+"<<yerr_high_vjets*1000<<"}_{-"<<yerr_low_vjets*1000<<"}^{+"<<RV_syst_err_up[UseCase]*1000<<"}_{-"<<RV_syst_err_low[UseCase]*1000<<"})\\cdot 10^{-3}"<<endl;

  //*******************************************************************************************
  // st+jets
  TGraphAsymmErrors *tg_stjets;
  pList->Add((TEfficiency*)tEff[3]); // st+jets s-ch
  pList->Add((TEfficiency*)tEff[4]); // st+jets s-ch
  pList->Add((TEfficiency*)tEff[5]); // st+jets t-ch
  pList->Add((TEfficiency*)tEff[6]); // st+jets t-ch
  pList->Add((TEfficiency*)tEff[7]); // st+jets tW-ch
  pList->Add((TEfficiency*)tEff[8]); // st+jets tW-ch
  // weights
  double w_stjets[6];
  w_stjets[0] = (channel==0 ?   27.02 :  20.33)/tEff[3]->GetTotalHistogram()->GetEntries();
  w_stjets[1] = (channel==0 ?   12.40 :   8.83)/tEff[4]->GetTotalHistogram()->GetEntries();
  w_stjets[2] = (channel==0 ?  243.65 : 175.97)/tEff[5]->GetTotalHistogram()->GetEntries();
  w_stjets[3] = (channel==0 ?  133.98 :  94.45)/tEff[6]->GetTotalHistogram()->GetEntries();
  w_stjets[4] = (channel==0 ? 1115.05 : 870.01)/tEff[7]->GetTotalHistogram()->GetEntries();
  w_stjets[5] = (channel==0 ? 1115.05 : 870.01)/tEff[8]->GetTotalHistogram()->GetEntries();

  tg_stjets = (TGraphAsymmErrors*)tDummy->Combine(pList,"mode",6,w_stjets);
//  tg_stjets->GetPoint(bin[UseCase],x,y_stjets);
//  yerr_low_stjets  = tg_stjets->GetErrorYlow(bin[UseCase]);
//  yerr_high_stjets = tg_stjets->GetErrorYhigh(bin[UseCase]);
  tg_stjets->GetPoint(bin-1,x,y_stjets);
  yerr_low_stjets  = tg_stjets->GetErrorYlow(bin-1);
  yerr_high_stjets = tg_stjets->GetErrorYhigh(bin-1);
  pList->Clear();
  cout<<"R_{st} = ("<<y_stjets*1000<<"^{+"<<yerr_high_stjets*1000<<"}_{-"<<yerr_low_stjets*1000<<"})\\cdot 10^{-3}"<<endl;

  //*******************************************************************************************
  // di-bosons
  TGraphAsymmErrors *tg_vvjets;
  pList->Add((TEfficiency*)tEff[9]);  // WW+jets
  pList->Add((TEfficiency*)tEff[10]); // WZ+jets
  pList->Add((TEfficiency*)tEff[11]); // ZZ+jets
  // weights
  double w_vvjets[3];
  w_vvjets[0] = (channel==0 ? 225.10 : 175.97)/tEff[8]->GetTotalHistogram()->GetEntries();
  w_vvjets[1] = (channel==0 ?  68.48 :  67.74)/tEff[9]->GetTotalHistogram()->GetEntries();
  w_vvjets[2] = (channel==0 ?   7.12 :  11.31)/tEff[10]->GetTotalHistogram()->GetEntries();

  tg_vvjets = (TGraphAsymmErrors*)tDummy->Combine(pList,"mode",3,w_vvjets);
//  tg_vvjets->GetPoint(bin[UseCase],x,y_vvjets);
//  yerr_low_vvjets  = tg_vvjets->GetErrorYlow(bin[UseCase]);
//  yerr_high_vvjets = tg_vvjets->GetErrorYhigh(bin[UseCase]);
  tg_vvjets->GetPoint(bin-1,x,y_vvjets);
  yerr_low_vvjets  = tg_vvjets->GetErrorYlow(bin-1);
  yerr_high_vvjets = tg_vvjets->GetErrorYhigh(bin-1);
  pList->Clear();
  cout<<"R_{VV} = ("<<y_vvjets*1000<<"^{+"<<yerr_high_vvjets*1000<<"}_{-"<<yerr_low_vvjets*1000<<"})\\cdot 10^{-3}"<<endl;

  //*******************************************************************************************
  // Multi-jets
  TGraphAsymmErrors *tg_multijets;
  double qcd[6];
  qcd[0] = 1281.98;
  qcd[1] =    0.00;
  qcd[2] = 2205.11;
  qcd[3] =  517.63;
  qcd[4] =  131.28;
  qcd[5] =   36.42;
  vector<double> w_multijet;

  if(channel==0){ // muon channel (1281.98 events at QC+4jets level)
    tDummy = (TEfficiency*)inputfiles_multijets[0]->Get(effname.c_str());
    tDummy->SetStatisticOption(TEfficiency::kBBayesian);
    y_multijets         = tDummy->GetEfficiency(bin);
    yerr_low_multijets  = tDummy->GetEfficiencyErrorLow(bin);
    yerr_high_multijets = tDummy->GetEfficiencyErrorUp(bin);
    cout<<"R_{multijet} = ("<<y_multijets*1000<<"^{+"<<yerr_high_multijets*1000<<"}_{-"<<yerr_low_multijets*1000<<"})\\cdot 10^{-3}"<<endl;
  }
  else{ // electron channel
    for(int i=1;i<6;i++){
      tDummy = (TEfficiency*)inputfiles_multijets[i]->Get(effname.c_str());
      tDummy->SetStatisticOption(TEfficiency::kBBayesian);
      if(qcd[i]>0){
        cout<<"Multi-jet "<<i<<" : eff = "<<tDummy->GetEfficiency(bin)<<endl;
        pList->Add(tDummy);
        w_multijet.push_back(qcd[i]/tDummy->GetTotalHistogram()->GetEntries());
      }
    }
    tDummy = 0;
    if(w_multijet.size()>0){
      tg_multijets = (TGraphAsymmErrors*)tDummy->Combine(pList,"mode",w_multijet.size(),&w_multijet[0]);
      tg_multijets->GetPoint(bin-1,x,y_multijets);
      yerr_low_multijets  = tg_multijets->GetErrorYlow(bin-1);
      yerr_high_multijets = tg_multijets->GetErrorYhigh(bin-1);
      cout<<"R_{multijet} = ("<<y_multijets*1000<<"^{+"<<yerr_high_multijets*1000<<"}_{-"<<yerr_low_multijets*1000<<"})\\cdot 10^{-3}"<<endl;
    }
    else cout<<"No QCD contribution"<<endl;
  }
  w_multijet.clear();
  pList->Clear();

  double QCDEst = 4172;
  double QCDEst_err = 4172*0.5;

  double TotUncert_Up[15];
  double TotUncert_Low[15];
    
  // adding syst. uncert. on RTT ;
  yerr_high_ttlike = sqrt(pow(yerr_high_ttlike,2)+pow(RTT_syst_err_up[UseCase],2));
  yerr_low_ttlike = sqrt(pow(yerr_low_ttlike,2)+pow(RTT_syst_err_low[UseCase],2));
  // adding syst. uncert. on RV ;
  yerr_high_vjets = sqrt(pow(yerr_high_vjets,2)+pow(RV_syst_err_up[UseCase],2));
  yerr_low_vjets  = sqrt(pow(yerr_low_vjets,2)+pow(RV_syst_err_low[UseCase],2));

  TotUncert_Up[0] = ErrorOnProd(y_stjets, yerr_high_stjets, nbOfEvents[0],statUncert[0]);
  TotUncert_Up[1] = ErrorOnProd(y_stjets, yerr_high_stjets, nbOfEvents[1],statUncert[1]);
  TotUncert_Up[2] = ErrorOnProd(y_stjets, yerr_high_stjets, nbOfEvents[2],statUncert[2]);

  TotUncert_Up[3] = ErrorOnProd(y_ttlike,yerr_high_ttlike,nbOfEvents[3],sqrt(pow(statUncert[3],2)+pow(systUncert[3],2)));
  TotUncert_Up[4] = ErrorOnProd(y_ttlike,yerr_high_ttlike,nbOfEvents[4],sqrt(pow(statUncert[4],2)+pow(systUncert[4],2)));
  TotUncert_Up[5] = ErrorOnProd(y_ttlike,yerr_high_ttlike,nbOfEvents[5],sqrt(pow(statUncert[5],2)+pow(systUncert[5],2)));

  TotUncert_Up[6] = ErrorOnProd(y_vjets, yerr_high_vjets, nbOfEvents[6]+nbOfEvents[ 9],sqrt(pow(statUncert[6],2)+pow(systUncert[6],2)+pow(statUncert[ 9],2)+pow(systUncert[ 9],2)),QCDEst*0.70,QCDEst_err*0.70);
  TotUncert_Up[7] = ErrorOnProd(y_vjets, yerr_high_vjets, nbOfEvents[7]+nbOfEvents[10],sqrt(pow(statUncert[7],2)+pow(systUncert[7],2)+pow(statUncert[10],2)+pow(systUncert[10],2)),QCDEst*0.25,QCDEst_err*0.25);
  TotUncert_Up[8] = ErrorOnProd(y_vjets, yerr_high_vjets, nbOfEvents[8]+nbOfEvents[11],sqrt(pow(statUncert[8],2)+pow(systUncert[8],2)+pow(statUncert[11],2)+pow(systUncert[11],2)),QCDEst*0.05,QCDEst_err*0.05);

  TotUncert_Up[9]  = 0; // Vbb component already added to V+jets
  TotUncert_Up[10] = 0; // Vbb component already added to V+jets
  TotUncert_Up[11] = 0; // Vbb component already added to V+jets

  TotUncert_Up[12] = ErrorOnProd(y_vvjets, yerr_high_vvjets, nbOfEvents[12],statUncert[12]);
  TotUncert_Up[13] = ErrorOnProd(y_vvjets, yerr_high_vvjets, nbOfEvents[13],statUncert[13]);
  TotUncert_Up[14] = ErrorOnProd(y_vvjets, yerr_high_vvjets, nbOfEvents[14],statUncert[14]);

  TotUncert_Low[0] = ErrorOnProd(y_stjets, yerr_low_stjets, nbOfEvents[0],statUncert[0]);
  TotUncert_Low[1] = ErrorOnProd(y_stjets, yerr_low_stjets, nbOfEvents[1],statUncert[1]);
  TotUncert_Low[2] = ErrorOnProd(y_stjets, yerr_low_stjets, nbOfEvents[2],statUncert[2]);

  TotUncert_Low[3] = ErrorOnProd(y_ttlike,yerr_low_ttlike,nbOfEvents[3],sqrt(pow(statUncert[3],2)+pow(systUncert[3],2)));
  TotUncert_Low[4] = ErrorOnProd(y_ttlike,yerr_low_ttlike,nbOfEvents[4],sqrt(pow(statUncert[4],2)+pow(systUncert[4],2)));
  TotUncert_Low[5] = ErrorOnProd(y_ttlike,yerr_low_ttlike,nbOfEvents[5],sqrt(pow(statUncert[5],2)+pow(systUncert[5],2)));

  TotUncert_Low[6] = ErrorOnProd(y_vjets, yerr_low_vjets, nbOfEvents[6]+nbOfEvents[ 9],sqrt(pow(statUncert[6],2)+pow(systUncert[6],2)+pow(statUncert[ 9],2)+pow(systUncert[ 9],2)),QCDEst*0.70,QCDEst_err*0.70);
  TotUncert_Low[7] = ErrorOnProd(y_vjets, yerr_low_vjets, nbOfEvents[7]+nbOfEvents[10],sqrt(pow(statUncert[7],2)+pow(systUncert[7],2)+pow(statUncert[10],2)+pow(systUncert[10],2)),QCDEst*0.25,QCDEst_err*0.25);
  TotUncert_Low[8] = ErrorOnProd(y_vjets, yerr_low_vjets, nbOfEvents[8]+nbOfEvents[11],sqrt(pow(statUncert[8],2)+pow(systUncert[8],2)+pow(statUncert[11],2)+pow(systUncert[11],2)),QCDEst*0.05,QCDEst_err*0.05);

  TotUncert_Low[9]  = 0; // Vbb component already added to V+jets
  TotUncert_Low[10] = 0; // Vbb component already added to V+jets
  TotUncert_Low[11] = 0; // Vbb component already added to V+jets

  TotUncert_Low[12] = ErrorOnProd(y_vvjets, yerr_low_vvjets, nbOfEvents[12],statUncert[12]);
  TotUncert_Low[13] = ErrorOnProd(y_vvjets, yerr_low_vvjets, nbOfEvents[13],statUncert[13]);
  TotUncert_Low[14] = ErrorOnProd(y_vvjets, yerr_low_vvjets, nbOfEvents[14],statUncert[14]);

  double Ntot_err_up      = 0;
  double Ntot_err_low     = 0;
/*
  for(int i=0;i<15;i++){
    Ntot_err_up  += pow(TotUncert_Up[i],2);
    Ntot_err_low += pow(TotUncert_Low[i],2);
    for(int j=0;j<15;j++){
      if(i==j) continue;
      Ntot_err_up += corrMatrix[i][j]*TotUncert_Up[i]*TotUncert_Up[j];
      Ntot_err_low+= corrMatrix[i][j]*TotUncert_Low[i]*TotUncert_Low[j];
    }
  }
  Ntot_err_up = sqrt(Ntot_err_up);
  Ntot_err_low= sqrt(Ntot_err_low);
  cout<<"Ntotal_err_up = "<<Ntot_err_up<<endl;
  cout<<"Ntotal_err_low = "<<Ntot_err_low<<endl;
*/
  double Rx[5]          = {y_stjets,y_ttlike,y_vjets,0,y_vvjets};
  double Rx_err_up[5]   = {yerr_high_stjets,yerr_high_ttlike,yerr_high_vjets,0,yerr_high_vvjets};
  double Rx_err_down[5] = {yerr_low_stjets,yerr_low_ttlike,yerr_low_vjets,0,yerr_low_vvjets};
  double N[5]           = {Nstjets,Nttlike,Nvlike+Nvbjets-QCDEst,0,Nvvjets};
  double Nerr[5]        = {0};
  double FinalError[5][5];
  double MeanCorr[5][5] = {{0}};
  
  for(int i=0;i<15;i=i+3){
    for(int j=0;j<15;j=j+3){
      for(int k=0;k<3;k++)
        for(int l=0;l<3;l++) MeanCorr[i/3][j/3] += (nbOfEvents[i+k]/(nbOfEvents[i]+nbOfEvents[i+1]+nbOfEvents[i+2]))*(nbOfEvents[j+l]/(nbOfEvents[j]+nbOfEvents[j+1]+nbOfEvents[j+2]))*corrMatrix[i+k][j+l];
    }
  }
  for(int i=0;i<15;i=i+3){
    Nerr[i/3] = sqrt(pow(statUncert[i],2)+pow(systUncert[i],2)+pow(statUncert[i+1],2)+pow(systUncert[i+1],2)+pow(statUncert[i+2],2)+pow(systUncert[i+2],2)
                     +2*sqrt(pow(statUncert[i],2)  +pow(systUncert[i],2))*sqrt(pow(statUncert[i+1],2)  +pow(systUncert[i+1],2))*corrMatrix[i][i+1]
                     +2*sqrt(pow(statUncert[i],2)  +pow(systUncert[i],2))*sqrt(pow(statUncert[i+2],2)  +pow(systUncert[i+2],2))*corrMatrix[i][i+2]
                     +2*sqrt(pow(statUncert[i+1],2)+pow(systUncert[i+1],2))*sqrt(pow(statUncert[i+2],2)+pow(systUncert[i+2],2))*corrMatrix[i+1][i+2]);
  }
  double SigmaN_up[5];
  for(int i=0;i<5;i++) SigmaN_up[i] = sqrt(pow(Rx_err_up[i]*N[i],2) + pow(Rx[i]*Nerr[i],2));
  double SigmaN_down[5];
  for(int i=0;i<5;i++) SigmaN_down[i] = sqrt(pow(Rx_err_down[i]*N[i],2) + pow(Rx[i]*Nerr[i],2));

  Ntot_err_up  = 0;
  Ntot_err_low = 0;
  
  float Ntotal = 0;
  
  for(int i=0;i<5;i++){
    if(i==3) continue;
    Ntotal += Rx[i]*N[i];
    Ntot_err_up  += pow(SigmaN_up[i],2);
    Ntot_err_low += pow(SigmaN_down[i],2);
    for(int j=0;j<5;j++){
      if(j==i | j==3) continue;
      Ntot_err_up += MeanCorr[i][j]*CorrCoef(Rx[i],Rx_err_up[i],N[i],Nerr[i])*CorrCoef(Rx[j],Rx_err_up[j],N[j],Nerr[j])*SigmaN_up[i]*SigmaN_up[j];
      Ntot_err_low+= MeanCorr[i][j]*CorrCoef(Rx[i],Rx_err_down[i],N[i],Nerr[i])*CorrCoef(Rx[j],Rx_err_down[j],N[j],Nerr[j])*SigmaN_down[i]*SigmaN_down[j];
    }
  }
  Ntot_err_up = sqrt(Ntot_err_up);
  Ntot_err_low= sqrt(Ntot_err_low);
  string comments[5] = {"Rst * Nst = ","Rtt * Ntt = ","Rv * Nv = ","Rvbb * Nvbb = ","Rvv * Nvv = "};
  cout<<"********************************************************"<<endl;
  cout<<"********************************************************"<<endl;
  for(int i=0;i<5;i++){
    cout<<comments[i]<<Rx[i]*N[i]<<"^{+"<<SigmaN_up[i]<<"}_{-"<<SigmaN_down[i]<<"}"<<endl;
  }
/*
  cout<<"Nttlike*RTT = "<<Nttlike*y_ttlike<<endl;
  cout<<"(Nvlike-Nqcd+Nvbjets)*RV = "<<(Nvlike-QCDEst+Nvbjets)*y_vjets<<endl;
  cout<<"Nstjets*Rst = "<<Nstjets*y_stjets<<endl;
  cout<<"Nvvjets*RVV = "<<Nvvjets*y_vvjets<<endl;
  cout<<"Nmultijets*RQCD = "<<QCDEst*y_multijets<<endl;
  
  cout<<"Nttlike*RTT+(Nvlike-Nqcd+Nvbjets)*RV = "<<Nttlike*y_ttlike+(Nvlike-QCDEst+Nvbjets)*y_vjets<<endl;
  cout<<"Nttlike*RTT+(Nvlike-Nqcd+Nvbjets)*RV+Nstjets*Rst+Nvvjets*RVV = "<<Nttlike*y_ttlike+(Nvlike-QCDEst+Nvbjets)*y_vjets+Nstjets*y_stjets+Nvvjets*y_vvjets<<endl;
  cout<<"Total = "<<Nttlike*y_ttlike+(Nvlike-QCDEst+Nvbjets)*y_vjets+Nstjets*y_stjets+Nvvjets*y_vvjets+QCDEst*y_multijets<<endl;
*/  
  cout<<"Ntotal = "<<Ntotal<<"^{+"<<Ntot_err_up<<"}_{-"<<Ntot_err_low<<"} / (+"<<100*(Ntot_err_up/Ntotal)<<"% / -"<<100*(Ntot_err_low/Ntotal)<<"%)"<<endl;
  cout<<"********************************************************"<<endl;
  cout<<"********************************************************"<<endl;
  for(int i=0;i<NbOfFiles;i++){
    inputfiles[i]->Close();
  }

}

double ErrorOnRatio(double A, double sigmaA, double B, double sigmaB){
  if(A<=0 || B<=0) return 0;
  double sigmaC = (A/B)*sqrt(pow(sigmaA/A,2)+pow(sigmaB/B,2));
  return sigmaC;
}
double ErrorOnProd(double R, double sigmaR, double N, double sigmaN){
  if(R<=0 || N<=0) return 0;
  double sigmaC = (R*N)*sqrt(pow(sigmaR/R,2)+pow(sigmaN/N,2));
  return sigmaC;
}
double ErrorOnProd(double R, double sigmaR, double N, double sigmaN, double Nest, double sigmaNest){
  if(R<=0 || N<=0) return 0;
  double sigmaD = (R*(N-Nest))*sqrt(pow(sigmaR/R,2)+(sigmaN*sigmaN+sigmaNest*sigmaNest)*pow(1/(N-Nest),2));
  return sigmaD;
}
double CorrCoef(double R, double sigmaR, double N, double sigmaN){
  double coef = 1 / sqrt(1 + pow((sigmaN/sigmaR)*(N/R),2));
  return coef;
}

