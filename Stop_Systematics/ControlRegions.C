
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>

#include "VJetEstimation.h"

void ControlRegions(std::string filename, int UseCase, int bin, bool UseWNJets, int channelConf) {
  
  std::ifstream file;
  
  const UInt_t nrpoints = 15;
  const UInt_t nreffs = 3;
  const UInt_t NbOfUseCase = 4;
  const UInt_t NbOfControlRegions = 2;
  const UInt_t NbOfMVARegions = 2; 
  Double_t ilel = 5050.821;
  Double_t ilmu = 5049.92;
  Double_t Wn_SF = 1.34;
  Int_t njets = 0 ; // 0 for 4jExc, 1 for 5jExc , 2 for 6jInc
  double **corrMatrix = new double*[nrpoints+nreffs];
  double *statUncert = new double[nrpoints+nreffs];
  //    double *systUncert = new double[nrpoints+nreffs];
  double *nbOfEvents = new double[nrpoints+nreffs];
  
  for(UInt_t i=0; i<nrpoints+nreffs; i++){
    corrMatrix[i] = new double[nrpoints+nreffs];
    statUncert[i] = 0;
    for(UInt_t j=0; j<nrpoints+nreffs; j++)
      corrMatrix[i][j] = 0;
  }
  
  std::string names[nrpoints+nreffs] = {"$N_{st+jets}^{4jets}$","$N_{st+jets}^{5jets}$","$N_{st+jets}^{6jets}$","$\\Nttlike^{4jets}$","$\\Nttlike^{5jets}$","$\\Nttlike^{6jets}$","$\\Nvlike^{4jets}$","$\\Nvlike^{5jets}$","$\\Nvlike^{6jets}$","$N_{Vbb}^{4jets}$","$N_{Vbb}^{5jets}$","$N_{Vbb}^{6jets}$","$N_{VV}^{4jets}$","$N_{VV}^{5jets}$","$N_{VV}^{6jets}$"    ,"$\\epsilon_{b}$","$\\epsilon_{uds}$","$\\epsilon_{udsc}$"};
  
  {
    UInt_t i=0;
    file.open(filename.c_str());
    while (i<nrpoints+nreffs){//(file.good()) {
      file >> corrMatrix[i][0] >> corrMatrix[i][1] >> corrMatrix[i][2] >> corrMatrix[i][3] >> corrMatrix[i][4] >> corrMatrix[i][5] >> corrMatrix[i][6] >> corrMatrix[i][7] >> corrMatrix[i][8] >> corrMatrix[i][9] >> corrMatrix[i][10] >> corrMatrix[i][11] >> corrMatrix[i][12] >> corrMatrix[i][13] >> corrMatrix[i][14]          >> corrMatrix[i][15]  >> corrMatrix[i][16]  >> corrMatrix[i][17] ;
      //cout<<std::fixed<<setprecision(1);
      cout<<names[i];
      for(UInt_t j=0;j<nrpoints+nreffs;j++) cout<<"& $"<<corrMatrix[i][j]*100<<"$ ";
      cout<<"\\\\"<<endl;
      ++i;
    }
    file.close();
  }
  VJetEstimation vj ;
  
  if (channelConf==0) {
    //SemiMuon
    std::vector<Double_t> va, vb, vc;
    Double_t a[7]={0,0,0,0, 0.022348, 0.016215, 0.013229 }; for (int k=0; k<7; k++) { va.push_back(a[k]); }
    Double_t b[7]={0,0,0,0, 0.283191, 0.238223, 0.211845 }; for (int k=0; k<7; k++) { vb.push_back(b[k]); }
    Double_t c[7]={0,0,0,0, 0.686278, 0.725439, 0.732775 }; for (int k=0; k<7; k++) { vc.push_back(c[k]); }
    vj.SetTTEffbq(va, vb, vc);
  } else if (channelConf==1) {
    //SemiElectron
    std::vector<Double_t> va, vb, vc;
    Double_t a[7]={0,0,0,0, 0.023797, 0.017614, 0.013601 }; for (int k=0; k<7; k++) { va.push_back(a[k]); }
    Double_t b[7]={0,0,0,0, 0.288904, 0.242580, 0.217566 }; for (int k=0; k<7; k++) { vb.push_back(b[k]); }
    Double_t c[7]={0,0,0,0, 0.679199, 0.719735, 0.727453 }; for (int k=0; k<7; k++) { vc.push_back(c[k]); }
    vj.SetTTEffbq(va, vb, vc);
  } else if (channelConf==2) {
    //SemiMuonSemiElectron
    std::vector<Double_t> va, vb, vc;
    Double_t a[7]={0,0,0,0, 0.022928, 0.016823, 0.013386 }; for (int k=0; k<7; k++) { va.push_back(a[k]); }
    Double_t b[7]={0,0,0,0, 0.284861, 0.240341, 0.214254 }; for (int k=0; k<7; k++) { vb.push_back(b[k]); }
    Double_t c[7]={0,0,0,0, 0.684076, 0.722721, 0.730549 }; for (int k=0; k<7; k++) { vc.push_back(c[k]); }
    vj.SetTTEffbq(va, vb, vc);
  }
  
  
  
  nbOfEvents[ 0] =  1030; // Nstjets 4 jets
  nbOfEvents[ 1] =   245; // Nstjets 5 jets
  nbOfEvents[ 2] =    60; // Nstjets 6 jets
  nbOfEvents[ 3] = 30249; // Nttlike 4 jets
  nbOfEvents[ 4] = 11321; // Nttlike 5 jets
  nbOfEvents[ 5] =  4030; // Nttlike 6 jets
  nbOfEvents[ 6] = 33343; // Nvlike  4 jets
  nbOfEvents[ 7] =  5745; // Nvlike  5 jets
  nbOfEvents[ 8] =  1052; // Nvlike  6 jets
  nbOfEvents[ 9] =  1855; // Nvblike 4 jets
  nbOfEvents[10] =   390; // Nvblike 5 jets
  nbOfEvents[11] =    90; // Nvblike 6 jets
  nbOfEvents[12] =   470; // Nvvjets 4 jets
  nbOfEvents[13] =    74; // Nvvjets 5 jets
  nbOfEvents[14] =    12; // Nvvjets 6 jets
  
  nbOfEvents[15] =   0.; // eb
  nbOfEvents[16] =   0.; // euds
  nbOfEvents[17] =   0.; // eudsc
  
  
  statUncert[ 0] =  303; // Nstjets 4 jets
  statUncert[ 1] =   73; // Nstjets 5 jets
  statUncert[ 2] =   18; // Nstjets 6 jets
  statUncert[ 3] = 1010; // Nttlike 4 jets
  statUncert[ 4] =  313; // Nttlike 5 jets
  statUncert[ 5] =  103; // Nttlike 6 jets
  statUncert[ 6] =  982; // Nvlike  4 jets
  statUncert[ 7] =  266; // Nvlike  5 jets
  statUncert[ 8] =   85; // Nvlike  6 jets
  statUncert[ 9] =  502; // Nvblike 4 jets
  statUncert[10] =  111; // Nvblike 5 jets
  statUncert[11] =   25; // Nvblike 6 jets
  statUncert[12] =  141; // Nvvjets 4 jets
  statUncert[13] =   22; // Nvvjets 5 jets
  statUncert[14] =    6; // Nvvjets 6 jets
  
  statUncert[15] =  0.; // eb
  statUncert[16] =  0.; // euds
  statUncert[17] =  0.; // eudsc
  
  /*
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
   */
  
  std::string MVA[NbOfUseCase] = { "LowDM" , "IntDM1", "IntDM2" , "HighDM" };
  std::string controlRegions[NbOfControlRegions] = { "4jExc_0b", "4jExc_2b" } ;
  std::string MVAcut[NbOfMVARegions+1] = { "_00MVA05" , "_05MVA08" , "" } ;
  
  
  
  
  
  const int NbOfFiles = 15;
  std::string BckgdNames[NbOfFiles];
  std::string CChannel[2] = {"_mu","_el"};
  
  printf("UseWnJets : %d\n", UseWNJets);
  printf("NN UseCase : %d (%s)\n", UseCase, MVA[UseCase].c_str());
  cout<<"Lepton channel : "<<(channelConf==0 ? CChannel[0] : (channelConf==1?CChannel[1]:"_combined"))<<endl;
  cerr<<"aerr"<< endl;
  printf("a\n"); 
  std::string path = "$HOME/AnalysisCode/GentStopAnalysis/UGentCode/Leg3/v1/Samples_21102012/";
  printf("b\n");
  std::vector<std::vector<TFile *> > inputfiles(5, std::vector<TFile*>(0, NULL));
  printf("c\n");
  std::vector<std::vector<std::string> > listNames(5, std::vector<std::string>(0, ""));
  printf("d\n");
  std::vector<std::vector<Double_t> > weight_onMC(5, std::vector<Double_t>(0, -1.));
  printf("e\n");
  std::vector<std::vector<Double_t> > weight_VJet(5, std::vector<Double_t>(0, -1.));
  printf("f\n");
  
  std::vector<TFile *> datafiles(0, NULL);
  printf("g\n");
  std::vector<std::string> dataNames(0, "");
  printf("h\n");
  std::vector<Double_t> weight_data(0, -1.);
  printf("i\n");
  std::vector<Double_t> weight_data_est(0, -1.);
  printf("j\n");
  
  for (int i=0; (1<<i)<=channelConf+1; i++) {
    Double_t entries = -1.;
    cerr << "1err" <<endl; 
    printf("1\n");
    Int_t channel=(1<<i)-1 ;
    printf("2\n");
    if (((channel+1)&(channelConf+1))==0) {
      continue;
    }
    printf("channel : %d\n", channel) ;
    // TT+jets
    BckgdNames[1] = "TT-like";
    listNames[1].push_back("ttbar"+CChannel[channel]);
    inputfiles[1].push_back(TFile::Open((path+"TTbar"+CChannel[channel]+".root").c_str()));
    entries = ((TH1D*) (* inputfiles[1].rbegin())->Get("Entries"))->GetBinContent(2);
    weight_onMC[1].push_back(165. *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //TTbar
    weight_VJet[1].push_back( nbOfEvents[3*1+njets] ); //TTbar
    
    BckgdNames[2] = "V-like";
    
    
    if(UseWNJets==kTRUE) {
      // W+jets exclusive
      //  BckgdNames[1] = "W2+jets";
      listNames[2].push_back("W2Jets"+CChannel[channel]);
      inputfiles[2].push_back(TFile::Open((path+"W2Jets"+CChannel[channel]+".root").c_str()));
      entries = ((TH1D*) (* inputfiles[2].rbegin())->Get("Entries"))->GetBinContent(2);
      //  BckgdNames[1] = "W3+jets";
      weight_onMC[2].push_back(Wn_SF*1435. *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //W2jets
      weight_VJet[2].push_back( nbOfEvents[3*2+njets] ); //W2jets
      
      listNames[2].push_back("W3Jets"+CChannel[channel]);
      inputfiles[2].push_back(TFile::Open((path+"W3Jets"+CChannel[channel]+".root").c_str()));
      entries = ((TH1D*) (* inputfiles[2].rbegin())->Get("Entries"))->GetBinContent(2);
      //  BckgdNames[1] = "W4+jets";
      weight_onMC[2].push_back(Wn_SF*304.2 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //W3jets
      weight_VJet[2].push_back( nbOfEvents[3*2+njets] ); //W3jets
      
      listNames[2].push_back("W4Jets"+CChannel[channel]);
      inputfiles[2].push_back(TFile::Open((path+"W4Jets"+CChannel[channel]+".root").c_str()));
      entries = ((TH1D*) (* inputfiles[2].rbegin())->Get("Entries"))->GetBinContent(2);
      weight_onMC[2].push_back(Wn_SF*172.6 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //W4jets
      weight_VJet[2].push_back( nbOfEvents[3*2+njets] ); //W4jets
    } else {
      // W+jets
      listNames[2].push_back("Wjets"+CChannel[channel]);
      inputfiles[2].push_back(TFile::Open((path+"Wjets"+CChannel[channel]+".root").c_str()));
      entries = ((TH1D*) (* inputfiles[2].rbegin())->Get("Entries"))->GetBinContent(2);
      weight_onMC[2].push_back(31314. *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //Wjets
      weight_VJet[2].push_back( nbOfEvents[3*2+njets] ); //Wjets
    }
    
    // Z+jets
    //  BckgdNames[2] = "z+jets";
    listNames[2].push_back("Zjets"+CChannel[channel]);
    inputfiles[2].push_back(TFile::Open((path+"Zjets"+CChannel[channel]+".root").c_str()));
    entries = ((TH1D*) (* inputfiles[2].rbegin())->Get("Entries"))->GetBinContent(2);
    weight_onMC[2].push_back(3048. *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //Zjets
    weight_VJet[2].push_back( nbOfEvents[3*2+njets] ); //Zjets
                                                       // multi-jets
    if (((channel+1)&1)!=0) {
      //  BckgdNames[0] = "QCDmu20";
      listNames[2].push_back("QCDmu20");
      inputfiles[2].push_back(TFile::Open((path+"QCDmu20_mu.root").c_str()));
      entries = ((TH1D*) (* inputfiles[2].rbegin())->Get("Entries"))->GetBinContent(2);
      weight_onMC[2].push_back(2.966E8*2.855E-4 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //QCDmu_20
      weight_VJet[2].push_back( nbOfEvents[3*2+njets] ); //QCDmu_20
    }
    if (((channel+1)&2)!=0) {
      //  BckgdNames[1] = "QCDel30to80";
      listNames[2].push_back("QCDel30to80");
      inputfiles[2].push_back(TFile::Open((path+"QCDel30to80_el.root").c_str()));
      entries = ((TH1D*) (* inputfiles[2].rbegin())->Get("Entries"))->GetBinContent(2);
      weight_onMC[2].push_back(5.944E7*0.061 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //QCDel_30_80
      weight_VJet[2].push_back( nbOfEvents[3*2+njets] ); //QCDel_30_80
                                                         //  BckgdNames[2] = "QCDel80to170";
      listNames[2].push_back("QCDel80to170");
      inputfiles[2].push_back(TFile::Open((path+"QCDel80to170_el.root").c_str()));
      entries = ((TH1D*) (* inputfiles[2].rbegin())->Get("Entries"))->GetBinContent(2);
      weight_onMC[2].push_back(898200*0.159 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //QCDel_80_170
      weight_VJet[2].push_back( nbOfEvents[3*2+njets] ); //QCDel_80_170
                                                         //  BckgdNames[3] = "QCDel170to250";
      listNames[2].push_back("QCDel170to250");
      inputfiles[2].push_back(TFile::Open((path+"QCDel170to250_el.root").c_str()));
      entries = ((TH1D*) (* inputfiles[2].rbegin())->Get("Entries"))->GetBinContent(2);
      weight_onMC[2].push_back(22140*0.1474 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //QCDel_170_250
      weight_VJet[2].push_back( nbOfEvents[3*2+njets] ); //QCDel_170_250
                                                         //  BckgdNames[4] = "QCDel250to35";
      listNames[2].push_back("QCDel250to350");
      inputfiles[2].push_back(TFile::Open((path+"QCDel250to350_el.root").c_str()));
      entries = ((TH1D*) (* inputfiles[2].rbegin())->Get("Entries"))->GetBinContent(2);
      weight_onMC[2].push_back(2900*0.1269 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //QCDel_250_350
      weight_VJet[2].push_back( nbOfEvents[3*2+njets] ); //QCDel_250_350
                                                         //  BckgdNames[5] = "QCDel350";
      listNames[2].push_back("QCDel350");
      inputfiles[2].push_back(TFile::Open((path+"QCDel350_el.root").c_str()));
      entries = ((TH1D*) (* inputfiles[2].rbegin())->Get("Entries"))->GetBinContent(2);
      weight_onMC[2].push_back(520*0.1058 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //QCDel_350
      weight_VJet[2].push_back( nbOfEvents[3*2+njets] ); //QCDel_350
    }
    
    
    
    // st+jets
    BckgdNames[0] = "Single Top";
    listNames[0].push_back("TopS"+CChannel[channel]);
    inputfiles[0].push_back(TFile::Open((path+"TopS"+CChannel[channel]+".root").c_str()));
    entries = ((TH1D*) (* inputfiles[0].rbegin())->Get("Entries"))->GetBinContent(2);
    weight_onMC[0].push_back(3.19 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //ST_s
    weight_VJet[0].push_back( nbOfEvents[3*0+njets] ); //ST_s
    
    //  BckgdNames[4] = "stbar+jets (s-ch)";
    listNames[0].push_back("TopbarS"+CChannel[channel]);
    inputfiles[0].push_back(TFile::Open((path+"TopbarS"+CChannel[channel]+".root").c_str()));
    entries = ((TH1D*) (* inputfiles[0].rbegin())->Get("Entries"))->GetBinContent(2);
    weight_onMC[0].push_back(1.44 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //STbar_s
    weight_VJet[0].push_back( nbOfEvents[3*0+njets] ); //STbar_s
    
    //  BckgdNames[5] = "st+jets (t-ch)";
    listNames[0].push_back("TopT"+CChannel[channel]);
    inputfiles[0].push_back(TFile::Open((path+"TopT"+CChannel[channel]+".root").c_str()));
    entries = ((TH1D*) (* inputfiles[0].rbegin())->Get("Entries"))->GetBinContent(2);
    weight_onMC[0].push_back(41.92 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //ST_t
    weight_VJet[0].push_back( nbOfEvents[3*0+njets] ); //ST_t
    
    //  BckgdNames[6] = "stbar+jets (t-ch)";
    listNames[0].push_back("TopbarT"+CChannel[channel]);
    inputfiles[0].push_back(TFile::Open((path+"TopbarT"+CChannel[channel]+".root").c_str()));
    entries = ((TH1D*) (* inputfiles[0].rbegin())->Get("Entries"))->GetBinContent(2);
    weight_onMC[0].push_back(22.65 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //STbar_t
    weight_VJet[0].push_back( nbOfEvents[3*0+njets] ); //STbar_t
    
    //  BckgdNames[7] = "st+jets (tW-ch)";
    listNames[0].push_back("ToptW"+CChannel[channel]);
    inputfiles[0].push_back(TFile::Open((path+"TopTW"+CChannel[channel]+".root").c_str()));
    entries = ((TH1D*) (* inputfiles[0].rbegin())->Get("Entries"))->GetBinContent(2);
    weight_onMC[0].push_back(7.87 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //ST_tW
    weight_VJet[0].push_back( nbOfEvents[3*0+njets] ); //ST_tW
    
    //  BckgdNames[8] = "stbar+jets (tW-ch)";
    listNames[0].push_back("TopbarTW"+CChannel[channel]);
    inputfiles[0].push_back(TFile::Open((path+"TopTW"+CChannel[channel]+".root").c_str()));
    entries = ((TH1D*) (* inputfiles[0].rbegin())->Get("Entries"))->GetBinContent(2);
    weight_onMC[0].push_back(7.87 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //ST_tW
    weight_VJet[0].push_back( nbOfEvents[3*0+njets] ); //ST_tW
    
    // di-boson
    BckgdNames[4]  = "VV";
    listNames[4].push_back("WW"+CChannel[channel]);
    inputfiles[4].push_back(TFile::Open((path+"WW"+CChannel[channel]+".root").c_str()));
    entries = ((TH1D*) (* inputfiles[4].rbegin())->Get("Entries"))->GetBinContent(2);
    weight_onMC[4].push_back(43. *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //WW
    weight_VJet[4].push_back( nbOfEvents[3*4+njets] ); //WW
                                                       //  BckgdNames[3] = "WZ+jets";
    listNames[4].push_back("WZ"+CChannel[channel]);
    inputfiles[4].push_back(TFile::Open((path+"WZ"+CChannel[channel]+".root").c_str()));
    entries = ((TH1D*) (* inputfiles[4].rbegin())->Get("Entries"))->GetBinContent(2);
    weight_onMC[4].push_back(18.2 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //WZ
    weight_VJet[4].push_back( nbOfEvents[3*4+njets] ); //WZ
                                                       //  BckgdNames[11] = "ZZ+jets";
    listNames[4].push_back("ZZ"+CChannel[channel]);
    inputfiles[4].push_back(TFile::Open((path+"ZZ"+CChannel[channel]+".root").c_str()));
    entries = ((TH1D*) (* inputfiles[4].rbegin())->Get("Entries"))->GetBinContent(2);
    weight_onMC[4].push_back(5.9 *(((channel+1)&1)!=0 ? (entries==0. ? 0. : ilmu/entries) : 0.) + (((channel+1)&2)!=0 ? (entries==0. ? 0. : ilel/entries) : 0.) ); //ZZ
    weight_VJet[4].push_back( nbOfEvents[3*4+njets] ); //ZZ
    
    
    
    
    if (((channel+1)&1)!=0) {
      //  BckgdNames[0] = "MuHad_2011A";
      dataNames.push_back("MuHad_2011A");
      datafiles.push_back(TFile::Open((path+"MuHad_2011A_mu.root").c_str()));
      weight_data.push_back( 1 ); //MuHad_2011A
      weight_data_est.push_back( 1 ); //MuHad_2011A
      
      //  BckgdNames[0] = "MuHad_2011B";
      dataNames.push_back("MuHad_2011B");
      datafiles.push_back(TFile::Open((path+"MuHad_2011B_mu.root").c_str()));
      weight_data.push_back( 1 ); //MuHad_2011B
      weight_data_est.push_back( 1 ); //MuHad_2011B
    }
    if (((channel+1)&2)!=0) {
      //  BckgdNames[0] = "ElHad2011A";
      dataNames.push_back("ElHad2011A");
      datafiles.push_back(TFile::Open((path+"ElHad_2011A_el.root").c_str()));
      weight_data.push_back( 1 ); //ElHad2011A
      weight_data_est.push_back( 1 ); //ElHad2011A
      
      //  BckgdNames[0] = "ElHad2011B";
      dataNames.push_back("ElHad2011B");
      datafiles.push_back(TFile::Open((path+"ElHad_2011B_el.root").c_str()));
      weight_data.push_back( 1 ); //ElHad2011B
      weight_data_est.push_back( 1 ); //ElHad2011B
    }
    
  }
  
  // Printing input info (independent of the Control Regions looked at)
  {
    printf("Printing input info (independent of the Control Regions looked at)\n");
    Double_t sum;
    printf(" weight_onMC : \n");
    for (UInt_t i=0; i<5; i++) {
      sum = 0.;
      printf("  Process %d : %s : ", i, BckgdNames[i].c_str());
      std::vector<Double_t>::const_iterator beg = weight_onMC[i].begin();
      for (std::vector<Double_t>::const_iterator it = beg; it!=weight_onMC[i].end(); it++) {
        if (it==beg) {
          printf("%lf", *it);
        } else {
          printf(" ; %lf", *it);
        }
        //                sum += *it ;
      }
      //            printf(" = %lf\n", sum);
      printf("\n");
    }
    printf(" weight_VJet : \n");
    for (UInt_t i=0; i<5; i++) {
      sum = 0.;
      printf("  Process %d : %s : ", i, BckgdNames[i].c_str());
      std::vector<Double_t>::const_iterator beg = weight_VJet[i].begin();
      for (std::vector<Double_t>::const_iterator it = beg; it!=weight_VJet[i].end(); it++) {
        if (it==beg) {
          printf("%lf", *it);
        } else {
          printf(" ; %lf", *it);
        }
        //                sum += *it ;
      }
      //            printf(" = %lf\n", sum);
      printf("\n");
    }
  }
  
  
  //    for (UInt_t l=0; l<NbOfUseCase; l++) {
  for (UInt_t j=0; j<NbOfControlRegions; j++) {
    std::string wbb_bJetMult_filename = "" ;
    TFile *wbb_bJetMult_file = TFile::Open(wbb_bJetMult_filename.c_str()) ;
    std::string wbbSuff[3] = { "4jets" , "5jets" , "geq6jets" } ;
    std::string channelS[3] = { "SemiMuon" , "SemiElectron" , "SemiMuonSemiElectron" };
    TH1F* wbb_bJetMult_hist = (TH1F*) wbb_bJetMult_file->Get((std::string()+"hNbtaggedJets_wjets_Wbb_"+channelS[channelConf]+"_"+wbbSuff[0]+"_WP0").c_str()) ;
    printf(" %s : { ", wbb_bJetMult_hist->GetName());
    for (int kkk=1; kkk<=wbb_bJetMult_hist->GetNbinsX() ; kkk++) {
      if (kkk==1) {
        printf("%lf", wbb_bJetMult_hist->GetBinContent(kkk));
      } else {
        printf(" ; %lf", wbb_bJetMult_hist->GetBinContent(kkk));
      }
    }
    printf(" } --> Integral : %lf \n",  wbb_bJetMult_hist->Integral(0,-1)) ;

    for (UInt_t k=0; k<=NbOfMVARegions; k++) {
      std::string teffname = "Eff_" + MVA[UseCase] + "_" + controlRegions[j] + MVAcut[(k==NbOfMVARegions?0:k)] ;
      if (k==NbOfMVARegions) {
        printf("\n\n\nNo TEfficiency applied (info from %s)     ... gives your info on the Control Region ...\n", teffname.c_str());
      } else {
        printf("\n\n\nInvestigated TEfficiency name : %s     ... gives your info on the Control Region ...\n", teffname.c_str());
      }
      /**
       Pure MC and extracting numbers from the files for the efficiencies
       */
      TGraphAsymmErrors **tg_categ = new TGraphAsymmErrors*[5];
      TGraphAsymmErrors *tg_tot = NULL;
      std::vector<std::vector<Double_t> > weights(5, std::vector<Double_t>()) ;
      std::vector<std::vector<Double_t> > passed(5, std::vector<Double_t>()) ;
      std::vector<Double_t> all_weights, all_weights_forCombineV ;
      std::vector<Double_t> all_passed ;
      std::vector<std::vector<TH1D*> > bJetMult(5, std::vector<TH1D*>() );
      std::vector<TH1D*> bJetMult_Avg(5, NULL); //Average b-jet mult on the different sub-processes
      TList *tlist_tot = new TList();
      for (UInt_t i=0; i<5; i++) { //Number of categories for V+jets method
        tg_categ[i] = NULL ;
        TList *tlist = new TList();
        Double_t sumppp = 0.;
        Double_t sumw = 0.;
        Double_t sumwMC = 0.;
        for (UInt_t m=0; m<inputfiles[i].size(); m++) {
          TEfficiency *teff = NULL;
          teff = (TEfficiency*) inputfiles[i][m]->Get(teffname.c_str());
          if (k==NbOfMVARegions) {
            TEfficiency *teff2 = teff;
            teff = new TEfficiency("", "", bin, 0, 1);
            teff->SetTotalEvents(bin, teff2->GetTotalHistogram()->GetBinContent(1+bin));
            teff->SetPassedEvents(bin, teff2->GetTotalHistogram()->GetBinContent(1+bin));
            /*
            teff->SetTotalEvents(bin, (((UInt_t) 0)-((UInt_t)1))/2);
            teff->SetPassedEvents(bin, (((UInt_t) 0)-((UInt_t)1))/2);
             */
          } else {
          }
          //if (teff==NULL) {
          //continue ;
          //}
          std::string fichName = inputfiles[i][m]->GetName() ;
          std::string chan = ( fichName.find("_mu.root")==fichName.size()-8 ? "Mu" : ( fichName.find("_el.root")==fichName.size()-8 ? "El" : "" ) );
          printf("%s\n", (std::string()+"B_Jet_Multiplicity_"+chan+"_4jExc").c_str());
          TH1D* bJetDistr = (TH1D*) inputfiles[i][m]->Get((std::string()+"B_Jet_Multiplicity_"+chan+"_4jExc").c_str());
          bJetMult[i].push_back(bJetDistr);
          if (bJetDistr!=NULL) {
            if (bJetMult_Avg[i]==NULL) {
              bJetMult_Avg[i] = (TH1D*) bJetDistr->Clone();
              bJetMult_Avg[i]->Scale( weight_onMC[i][m] /*/((TH1D*) inputfiles[i][m]->Get("Entries"))->GetBinContent(2)*/ );
            } else {
              bJetMult_Avg[i]->Add(bJetDistr, weight_onMC[i][m] /*/((TH1D*) inputfiles[i][m]->Get("Entries"))->GetBinContent(2)*/ );
            }
          }
          Double_t w = weight_onMC[i][m]*teff->GetTotalHistogram()->GetBinContent(1+bin) ;
          Double_t ppp = weight_onMC[i][m]*teff->GetPassedHistogram()->GetBinContent(1+bin) ;
          sumw += w;
          sumppp += ppp;
          sumwMC += weight_onMC[i][m];
          weights[i].push_back(w); //Number of events (after the cut, before the TEff cut), pure MC
          all_weights.push_back(w * weight_VJet[i][m]); //number of events after cut (before the TEff cut) from MC and VJetEstimation (both taken into account)
          all_weights_forCombineV.push_back(weight_onMC[i][m] * weight_VJet[i][m]); //weight for combination (number of events after cut) from MC and VJetEstimation (both taken into account)
          passed[i].push_back((weight_VJet[i][m]==0. ? 0. : ppp));
          all_passed.push_back(ppp * weight_VJet[i][m]); //number of passed element
          if (teff==NULL) { printf("NULL TEfficiency\n"); }
          tlist->Add(teff);
          tlist_tot->Add(teff);                        
        }
        for (UInt_t kk=all_passed.size()-passed.size(); kk<all_passed.size(); kk++) {
          //                        bJetMult_Avg[i]->Scale(1./);
          all_passed[kk] /= sumppp;
          all_weights[kk] /= sumw;
          all_weights_forCombineV[kk] /= sumwMC;
        }
        if (weight_onMC[i].size()!=0) { 
          cerr << "Category combination"<< endl;
          tg_categ[i] = TEfficiency::Combine(tlist, "mode", weight_onMC[i].size(), & weight_onMC[i][0]);
        }
        tlist->Clear();
      }
      if (all_weights_forCombineV.size()!=0) {
        cerr << "Whole combination" << endl;
        tg_tot = TEfficiency::Combine(tlist_tot, "mode", all_weights_forCombineV.size(), & all_weights_forCombineV[0]);
        cerr << "Whole combined" << endl;
      }
      
      // Printing input values
      {
        printf("\n\nPrinting input values (for checks) : \n");
        Double_t sum;
        printf(" passed : \n");
        for (UInt_t i=0; i<5; i++) {
          sum = 0.;
          printf("  Process %d : %s : ", i, BckgdNames[i].c_str());
          std::vector<Double_t>::const_iterator beg = passed[i].begin();
          for (std::vector<Double_t>::const_iterator it = beg; it!=passed[i].end(); it++) {
            if (it==beg) {
              printf("%lf", *it);
            } else {
              printf(" + %lf", *it);
            }
            sum += *it ;
          }
          printf(" = %lf\n", sum);
        }
        printf(" weights : \n");
        for (UInt_t i=0; i<5; i++) {
          sum = 0.;
          printf("  Process %d : %s : ", i, BckgdNames[i].c_str());
          std::vector<Double_t>::const_iterator beg = weights[i].begin();
          for (std::vector<Double_t>::const_iterator it = beg; it!=weights[i].end(); it++) {
            if (it==beg) {
              printf("%lf", *it);
            } else {
              printf(" + %lf", *it);
            }
            sum += *it ;
          }
          printf(" = %lf\n", sum);
        }
        sum = 0.;
        printf(" all_passed : \n  ");
        std::vector<Double_t>::const_iterator beg = all_passed.begin();
        for (std::vector<Double_t>::const_iterator it = beg; it!=all_passed.end(); it++) {
          if (it==beg) {
            printf("%lf", *it);
          } else {
            printf(" + %lf", *it);
          }
          sum += *it ;
        }
        printf(" = %lf\n", sum);
        
        sum = 0.;
        printf(" all_weights : \n  ");
        beg = all_weights.begin();
        for (std::vector<Double_t>::const_iterator it = beg; it!=all_weights.end(); it++) {
          if (it==beg) {
            printf("%lf", *it);
          } else {
            printf(" + %lf", *it);
          }
          sum += *it ;
        }
        printf(" = %lf\n", sum);
        
        sum = 0.;
        printf(" all_weights_forCombineV : \n  ");
        beg = all_weights_forCombineV.begin();
        for (std::vector<Double_t>::const_iterator it = beg; it!=all_weights_forCombineV.end(); it++) {
          if (it==beg) {
            printf("%lf", *it);
          } else {
            printf(" + %lf", *it);
          }
          sum += *it ;
        }
        printf(" = %lf\n", sum);
        
        for (UInt_t i=0; i<5; i++) {
          if (bJetMult_Avg[i]==NULL) {
            continue;
          }
          printf(" bJetMult_Avg[%d] : { ", i);
          for (int kkk=1; kkk<=bJetMult_Avg[i]->GetNbinsX() ; kkk++) {
            if (kkk==1) {
              printf("%lf", bJetMult_Avg[i]->GetBinContent(kkk));
            } else {
              printf(" ; %lf", bJetMult_Avg[i]->GetBinContent(kkk));
            }
          }
          printf(" } --> Integral : %lf \n",  bJetMult_Avg[i]->Integral(0,-1)) ;
        }
      }
      
      // Printing results for pure MC (non reweighted)
      {
        Double_t tot = 0.;
        Double_t vlike_plus_bb = 0.;
        printf("\n\nPrinting results for pure MC (non reweighted) : \n");
        for (UInt_t i=0; i<5; i++) {
          printf(" Process %d : %s : ", i, BckgdNames[i].c_str());
          Double_t sum = 0.;
          for (UInt_t m=0; m<inputfiles[i].size(); m++) {
            if (bJetMult[i][m] == NULL) {
              continue ;
            }
            Double_t bBinFrac = bJetMult[i][m]->GetBinContent((j==0 ? 1 : (j==1 ? 3 : 0))) / bJetMult[i][m]->Integral(0,-1);
            if (m==0) {
              printf(" %lf", bBinFrac * passed[i][m]);
            } else {
              printf(" + %lf", bBinFrac * passed[i][m]);
            }
            sum += bBinFrac * passed[i][m] ;
          }
          printf(" = %lf\n", sum);
          if (i==2/* || i==3*/) {
            vlike_plus_bb += sum ;
          }
          tot += sum ;
        }
        printf("  V-like + Wbb category (sum) : %lf\n", vlike_plus_bb);
        printf("  Total : %lf\n", tot);
      }
      /**
       Estimated with V+jets on data (and R_X from MC)
       */
      {
        Double_t tot=0., bTot=0.;
        Double_t vlike_plus_bb = 0., bVlike_plus_bb=0.;
        Double_t ntt=0., ntt_err=0., nv=0., nv_err=0. ;
        if (j==0) {
          ntt     = vj.Ntt_0bjet(nbOfEvents[3*1+njets], nbOfEvents[nrpoints+0], nbOfEvents[nrpoints+2], 4+njets);
          ntt_err = vj.Ntt_err_0bjet(nbOfEvents[3*1+njets], statUncert[3*1+njets], nbOfEvents[nrpoints+0], statUncert[nrpoints+0], nbOfEvents[nrpoints+2], statUncert[nrpoints+2], 4+njets);
          nv      = vj.Nv_0bjet(nbOfEvents[3*2+njets], nbOfEvents[nrpoints+1], 4+njets);
          nv_err  = vj.Nv_err_0bjet(nbOfEvents[3*2+njets], statUncert[3*2+njets], nbOfEvents[nrpoints+1], statUncert[nrpoints+1], 4+njets);
        } else if (j==1) {
          ntt     = vj.Ntt_2bjets(nbOfEvents[3*1+njets], nbOfEvents[nrpoints+0], nbOfEvents[nrpoints+2], 4+njets);
          ntt_err = vj.Ntt_err_2bjets(nbOfEvents[3*1+njets], statUncert[3*1+njets], nbOfEvents[nrpoints+0],  statUncert[nrpoints+0], nbOfEvents[nrpoints+2], statUncert[nrpoints+2], 4+njets);
          nv      = vj.Nv_2bjets(nbOfEvents[3*2+njets], nbOfEvents[nrpoints+1], 4+njets);
          nv_err  = vj.Nv_err_2bjets(nbOfEvents[3*2+njets], statUncert[3*2+njets], nbOfEvents[nrpoints+1], statUncert[nrpoints+1], 4+njets);
        }
        printf("\n\nPrinting results for the estimated channels : \n");
        Double_t ytt = 0., ytemp=0., ytterr=0., dummy;
        if (tg_categ[1] != NULL) {
          tg_categ[1]->GetPoint(bin, dummy,ytt);
          ytterr = tg_categ[1]->GetErrorYlow(bin);
          ytemp = tg_categ[1]->GetErrorYhigh(bin);
          if (ytemp>ytterr) {
            ytterr = ytemp;
          }
        }
        if (k==NbOfMVARegions) {
          ytt = 1. ; ytterr = 0. ;
        }
        tot += ntt*ytt ;
        bTot += nbOfEvents[3*1+njets]*ytt ;
        printf("  Ntt = ( %lf \\pm %lf ) * ( %lf \\pm %lf ) = %lf \\pm %lf \t\t\t\t //// Ntt tot. for jet mult. : %lf \n", ntt, ntt_err, ytt, ytterr, ntt*ytt, (ntt_err/ntt + ytterr/ytt) *ntt*ytt   , nbOfEvents[3*1+njets]);
        printf("b Ntt = ( %lf \\pm %lf ) * ( %lf \\pm %lf ) = %lf \\pm %lf \t\t\t\t //// Ntt tot. for jet mult. : %lf \n", nbOfEvents[3*1+njets], statUncert[3*1+njets], ytt, ytterr, nbOfEvents[3*1+njets]*ytt, (statUncert[3*1+njets]/nbOfEvents[3*1+njets] + ytterr/ytt) *nbOfEvents[3*1+njets]*ytt   , nbOfEvents[3*1+njets]);
        ytemp=0.;
        Double_t yv = 0., yverr=0.;
        if (tg_categ[2] != NULL) {
          tg_categ[2]->GetPoint(bin, dummy,yv);
          yverr = tg_categ[2]->GetErrorYlow(bin);
          ytemp = tg_categ[2]->GetErrorYhigh(bin);
          if (ytemp>yverr) {
            yverr = ytemp;
          }
        }
        if (k==NbOfMVARegions) {
          yv = 1. ; yverr = 0. ;
        }
        tot += nv*yv ;
        bTot += nbOfEvents[3*2+njets]*yv ;
        vlike_plus_bb += nv*yv;
        bVlike_plus_bb += nbOfEvents[3*2+njets]*yv;
        printf("  Nv = ( %lf \\pm %lf ) * ( %lf \\pm %lf ) = %lf \\pm %lf \t\t\t\t //// Nv tot. for jet mult. : %lf \n", nv, nv_err, yv, yverr, nv*yv, (nv_err/nv + yverr/yv) *nv*yv   , nbOfEvents[3*2+njets]);
        printf("b Nv = ( %lf \\pm %lf ) * ( %lf \\pm %lf ) = %lf \\pm %lf \t\t\t\t //// Nv tot. for jet mult. : %lf \n", nbOfEvents[3*2+njets], statUncert[3*2+njets], yv, yverr, nbOfEvents[3*2+njets]*yv, (statUncert[3*2+njets]/nbOfEvents[3*2+njets] + yverr/yv) *nbOfEvents[3*2+njets]*yv   , nbOfEvents[3*2+njets]);
        for (UInt_t i=0; i<5; i++) {
          printf("  Process %d : %s : ", i, BckgdNames[i].c_str());
          Double_t sum = 0.;
          Double_t y = 0., yerr=0.;
          ytemp=0.;
          if (tg_categ[i] != NULL) {
            tg_categ[i]->GetPoint(bin, dummy,y);   
            yerr = tg_categ[i]->GetErrorYlow(bin);
            ytemp = tg_categ[i]->GetErrorYhigh(bin);
            if (ytemp>yerr) {
              yerr = ytemp;
            }
          }
          if (k==NbOfMVARegions) {
            y = 1. ; yerr = 0. ;
          }
          if (bJetMult_Avg[i] == NULL && i!=3) {
            printf("\n");
            continue ;
          }
          Double_t bBinFrac = -1.;
          if (i!=3) {
            bBinFrac = bJetMult_Avg[i]->GetBinContent((j==0 ? 1 : (j==1 ? 3 : 0))) / bJetMult_Avg[i]->Integral(0,-1);
            sum = y * bBinFrac * nbOfEvents[3*i+njets] ;
            if (i!=1 && i!=2) {
              tot += sum;
              bTot += y * nbOfEvents[3*i+njets] ;
            }
          } else if (i==3) {
            // Same factors as V-like
            y = yv ;
            yerr = yverr ;
            if (k==NbOfMVARegions) {
              y = 1. ; yerr = 0. ;
            }
            bBinFrac = bJetMult_Avg[2]->GetBinContent((j==0 ? 1 : (j==1 ? 3 : 0))) / bJetMult_Avg[2]->Integral(0,-1);
            sum = yv * bBinFrac * nbOfEvents[3*3+njets];
            
            Double_t factWbbEff = VJetEstimation::probElemWbb(nbOfEvents[nrpoints], nbOfEvents[nrpoints+1], (j==1 ? 2 : 0), 4+njets)
            - VJetEstimation::probElemWbb(nbOfEvents[nrpoints], nbOfEvents[nrpoints+1], 0, (j==1 ? 2 : 0), 4+njets) ;
            /*
             Double_t factWbbEff = 0.; // Taken from model
             if (j==0) {
             factWbbEff = pow(1.-nbOfEvents[nrpoints+1]+1.-nbOfEvents[nrpoints+0], 4.+njets) - pow(1.-nbOfEvents[nrpoints+1], 4.+njets) ;
             } else if (j==1) {
             factWbbEff = pow(1.-nbOfEvents[nrpoints+1]+1.-nbOfEvents[nrpoints+0], 4.+njets) - pow(1.-nbOfEvents[nrpoints+1], 4.+njets) ;
             }
             */
            vlike_plus_bb += sum ;
            printf("[factWbb:est(%lf) Vest(%lf) VjetsMethods(%lf)]\t\t\t", factWbbEff, bBinFrac, wbb_bJetMult_hist->GetBinContent((j==0 ? 1 : (j==1 ? 3 : 0))) / wbb_bJetMult_hist->Integral(0,-1));
            tot += sum ;
            bTot += yv * nbOfEvents[3*3+njets] ;
            bVlike_plus_bb += yv * nbOfEvents[3*3+njets] ;
          }
          printf("%lf * %lf * %lf = %lf(bF) * ( %lf \\pm %lf ) = %lf \\pm %lf \n", y, bBinFrac, nbOfEvents[3*i+njets], bBinFrac, y*nbOfEvents[3*i+njets], y*nbOfEvents[3*i+njets]*(statUncert[3*i+njets]/nbOfEvents[3*i+njets] + yerr/y), sum, sum * (statUncert[3*i+njets]/nbOfEvents[3*i+njets] + yerr/y) );
          
        }
        printf("  V-like + Wbb category (sum) : %lf\n", vlike_plus_bb);
        printf("b V-like + Wbb category (sum) : %lf\n", bVlike_plus_bb);
        printf("  Total : %lf\n", tot);
        printf("b Total : %lf\n", bTot);
      }
      
      /**
       Observed in data (only totals ...)
       */
      {
        printf("\n\nPrinting results for pure data (total numbers only ...) : \n");
        Double_t sum = 0.;
        printf("  From B_Jet_Multiplicities \n");
        for (UInt_t i=0; i<datafiles.size(); i++) {
          std::string fichName = datafiles[i]->GetName() ;
          std::string chan = ( fichName.find("_mu.root")==fichName.size()-8 ? "Mu" : ( fichName.find("_el.root")==fichName.size()-8 ? "El" : "" ) );
          TH1D* histo = (TH1D*) datafiles[i]->Get((std::string()+"B_Jet_Multiplicity_"+chan+"_4jExc").c_str());
          if (histo==NULL) {
            continue;
          }
          Double_t bBinFrac = histo->GetBinContent((j==0 ? 1 : (j==1 ? 3 : 0)));
          if (i==0) {
            printf(" %lf", bBinFrac );
          } else {
            printf(" + %lf", bBinFrac);
          }
          sum += bBinFrac ;
        }
        printf(" = %lf\n", sum);
        printf("  From TEfficiencies \n");
        printf("   total : \n");
        sum=0.;
        for (UInt_t i=0; i<datafiles.size(); i++) {
          TEfficiency *teff = (TEfficiency*) datafiles[i]->Get(teffname.c_str());
          Double_t bBinFrac = teff->GetTotalHistogram()->GetBinContent(1+bin);
          if (i==0) {
            printf(" %lf", bBinFrac );
          } else {
            printf(" + %lf", bBinFrac);
          }
          sum += bBinFrac ;
        }
        printf(" = %lf\n", sum);
        printf("   passed : \n");
        sum=0.;
        for (UInt_t i=0; i<datafiles.size(); i++) {
          TEfficiency *teff = (TEfficiency*) datafiles[i]->Get(teffname.c_str());
          Double_t bBinFrac = teff->GetPassedHistogram()->GetBinContent(1+bin);
          if (i==0) {
            printf(" %lf", bBinFrac );
          } else {
            printf(" + %lf", bBinFrac);
          }
          sum += bBinFrac ;
        }
        printf(" = %lf\n", sum);
      }
    }
    wbb_bJetMult_file->Close();
  }
  //   }
  
  
  for (std::vector<std::vector<TFile*> >::iterator it = inputfiles.begin() ; it != inputfiles.end() ; it++ ) {
    for (std::vector<TFile*>::iterator it2 = it->begin() ; it2 != it->end() ; it2++ ) {
      (*it2)->Close();
    }
  }
  
  printf(" --> End of the program.\n") ;    
}
