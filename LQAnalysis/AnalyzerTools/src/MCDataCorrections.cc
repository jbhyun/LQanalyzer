#include "MCDataCorrections.h"

// std includes
#include <iostream>
#include <vector>
#include <string>

// ROOT includes
#include <TFile.h>

/// local includes
#include "Reweight.h"


using namespace snu;
using namespace std;

MCDataCorrections::MCDataCorrections() {
  
  corr_isdata=false;
  k_period=-1;

  string lqdir = getenv("LQANALYZER_DIR");
  rc =  new RoccoR(lqdir + "/data/rochester/80X/rcdata.2016.v3");  
  
  deg_etaptmap_leg1.clear();
  deg_etaptmap_leg2.clear();


  CorrectionMap.clear();
  CorrectionMapGraph.clear();
  
  string pileupdir = getenv("PILEUPFILEDIR");
  
  SetupDoubleEGTriggerSF(1,  "/data/Trigger/80X/HWW_HLT_DoubleEleLegHighPt.txt");                                                                       
  SetupDoubleEGTriggerSF(2,  "/data/Trigger/80X/HWW_HLT_DoubleEleLegLowPt.txt");      

  FillCorrectionHists();
  reweightPU = new Reweight((pileupdir + "/DataPileUp_2016_Moriond.root").c_str());       

}


MCDataCorrections::MCDataCorrections(bool isdata) {
  corr_isdata=isdata;
  MCDataCorrections();
}
  

MCDataCorrections::~MCDataCorrections(){
  delete rc;
  delete reweightPU;
  for( std::map<float, std::vector<float>* >::iterator mit = deg_etaptmap_leg1.begin(); mit != deg_etaptmap_leg1.end(); mit++){
    delete mit->second;
  }
  for( std::map<float, std::vector<float>* >::iterator mit = deg_etaptmap_leg2.begin();mit != deg_etaptmap_leg2.end();mit++){
    delete mit->second;
  }
  CorrectionMap.clear();
  CorrectionMapGraph.clear();
  deg_etaptmap_leg1.clear();
  deg_etaptmap_leg2.clear();


}


void MCDataCorrections::SetPeriod(int period){
  k_period=period;
}

void MCDataCorrections::SetIsData(bool isdata){
  corr_isdata=isdata;
}

void MCDataCorrections::PrintSummary(){
  //// summarize results
}



void  MCDataCorrections::SetupDoubleEGTriggerSF(int ileg, string sleg){

  string analysisdir = getenv("LQANALYZER_DIR");

  string file = analysisdir + sleg;
  ifstream trigsf_file(file.c_str());

  string sline;
  float etaboundary(-999.);
  while(getline(trigsf_file,sline) ){
    std::istringstream is( sline );
    if(TString(sline).Contains("##")) continue;

    float eta_min;
    is >> eta_min;
    if(etaboundary != eta_min) {
      std::vector<float>* vpt = new std::vector<float>();
      if(ileg==1)deg_etaptmap_leg1[eta_min] = vpt;
      if(ileg==2)deg_etaptmap_leg2[eta_min] = vpt;
    }
    etaboundary = eta_min;

    if(sline=="END") break;
  }

  std::map<float, std::vector<float>* >::iterator it;
  etaboundary = -999.;
  ifstream trigsf_file2(file.c_str());

  while(getline(trigsf_file2,sline) ){

    std::istringstream is( sline );
    if(TString(sline).Contains("##")) continue;
    float eta_min, eta_max, pt_min, pt_max, _sf, _sferr;
    is >> eta_min;
    
    if(ileg==1)it = deg_etaptmap_leg1.find(eta_min);
    if(ileg==2)it = deg_etaptmap_leg2.find(eta_min);
    if(ileg==1){
      if(it == deg_etaptmap_leg1.end()){
	cout << "Error in map for double eg trigger SF" << endl;
	exit(EXIT_FAILURE);
      }
    }
    if(ileg==2){
      if(it == deg_etaptmap_leg2.end()){
        cout << "Error in map for double eg trigger SF" << endl;
        exit(EXIT_FAILURE);
      }
    }
    is >> eta_max;
    is >> pt_min;
    is >> pt_max;
    is >> _sf;
    is >> _sferr;


    vector<float>* v_sf = it->second;
    v_sf->push_back(_sf);
    if(ileg==1)deg_etaptmap_leg1[eta_min] = v_sf;
    if(ileg==2)deg_etaptmap_leg2[eta_min] = v_sf;

    if(sline=="END") break;
  }

  return;
}





void MCDataCorrections::FillCorrectionHists(){

  string file = string(getenv("LQANALYZER_DIR")) + "/CATConfig/CattupleConfig/corrections_"+string(getenv("CATVERSION"))+".txt" ;
  cout << "Correction file = " << file << endl;
  ifstream corr_file(file.c_str());
  string sline;
  while(getline(corr_file,sline) ){
    std::istringstream is( sline );
    if(sline.find("###")!=string::npos) {  
      string tmp;
      is >> tmp;
      string label;
      is >> label;
      string dirname;
      is >> dirname;
      string filename;
      is >> filename;
      string histsname;
      is >> histsname;
      string histtype;
      is >> histtype;
      FillCorrectionHist(label,dirname,filename,histsname,histtype);
    }
    if(sline=="#END") break;
  }
}

void MCDataCorrections::FillCorrectionHist(string label, string dirname, string filename, string histsname, string histtype){

  TDirectory* origDir = gDirectory;

  if(TString(histtype).Contains("TGraphAsymmErrors")) {
    TFile *infile_sf = TFile::Open((string(getenv(dirname.c_str()))+ "/" + filename).c_str());
    CheckFile(infile_sf);
    TDirectory* tempDir = getTemporaryDirectory();
    tempDir->cd();
    TGraphAsymmErrors* tmp = (TGraphAsymmErrors*)infile_sf->Get(histsname.c_str())->Clone();
    CorrectionMapGraph[label] = tmp;
    infile_sf->Close();
    delete infile_sf;
    origDir->cd();
  }

  if(!TString(histtype).Contains("TH2")) return;
  TFile *infile_sf = TFile::Open((string(getenv(dirname.c_str()))+ "/" + filename).c_str());
  CheckFile(infile_sf);
  TDirectory* tempDir = getTemporaryDirectory();
  tempDir->cd();

  //==== for leg-by-leg strategy, we should get efficiencies, not scale factors
  TString tmp_histsname = histsname;
  TString tmp_filename = filename;
  if(tmp_filename.Contains("MuonTriggerEfficiency")){
    cout << "[MuonTriggerEfficiency]" << endl;
    TString tmp_label = label;
    cout << "label = "  << label << endl;
    //==== Data
    TH2F* tmp_Data =  dynamic_cast<TH2F*> (( infile_sf->Get((histsname+"_Data").c_str()))->Clone());
    CorrectionMap[label+"_Data"] = tmp_Data;
    cout << "CorrectionMap["<<label+"_Data" <<"] = " << histsname+"_Data" << " (from " << getenv(dirname.c_str())<< "/" << filename<<")" << endl;
    //==== MC
    TH2F* tmp_MC =  dynamic_cast<TH2F*> (( infile_sf->Get((histsname+"_MC").c_str()))->Clone());
    CorrectionMap[label+"_MC"] = tmp_MC;
    cout << "CorrectionMap["<<label+"_MC" <<"] = " << histsname+"_MC" << " (from " << getenv(dirname.c_str())<< "/" << filename<<")" << endl;
  }
  else{
    TH2F* tmp =  dynamic_cast<TH2F*> (( infile_sf->Get(histsname.c_str()))->Clone());
    CorrectionMap[label] = tmp;
    cout << "CorrectionMap["<<label <<"] = " << histsname << " (from " << getenv(dirname.c_str())<< "/" << filename<<")" << endl;
  }
  infile_sf->Close();
  delete infile_sf;
  origDir->cd();

}


double MCDataCorrections::MuonTrackingEffScaleFactor(vector<snu::KMuon> mu){
  float sf= 1.;
  if(corr_isdata) return 1.;
  for(vector<KMuon>::iterator itmu=mu.begin(); itmu!=mu.end(); ++itmu) {
    float mueta=itmu->Eta();
    if(CheckCorrectionGraph("TRACKING_EFF")){
      sf*= GetCorrectionGraph("TRACKING_EFF")->Eval(mueta);
    }
  }
  return sf;
}

double MCDataCorrections::MuonISOScaleFactor(TString muid, vector<snu::KMuon> mu,int sys){

  if(k_period < 0) {
    /// If k_period < 0 then using ALL data periods and use weighted SF
    
    double lumi_periodB = 5.929001722;
    double lumi_periodC = 2.645968083;
    double lumi_periodD = 4.35344881;
    double lumi_periodE = 4.049732039;
    double lumi_periodF = 3.157020934;
    double lumi_periodG = 7.549615806;
    double lumi_periodH = 8.545039549 + 0.216782873;
    double total_lumi = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG+lumi_periodH);

    double WeightBtoF = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF)/total_lumi;
    double WeightGtoH = (lumi_periodG+lumi_periodH)/total_lumi;
    double SF_bf= MuonISOScaleFactorPeriodDependant(muid,mu, 1, sys);
    double SF_gh= MuonISOScaleFactorPeriodDependant(muid,mu, 7, sys);
    
    double SF_weight = WeightBtoF*SF_bf + WeightGtoH*SF_gh;
    return SF_weight;
  }
  return MuonISOScaleFactorPeriodDependant(muid, mu, k_period, sys);

}

double MCDataCorrections::MuonISOScaleFactorPeriodDependant(TString muid, vector<snu::KMuon> mu , int cat_period, int sys){
  float sf= 1.;
  float sferr=1.;
  if(corr_isdata) return 1.;

  // ref https://indico.cern.ch/event/595070/contributions/2405098/attachments/1388788/2114715/TnPIso12Dic2016.pdf 

  TString tag = "";
  if(cat_period < 6) tag = "_BCDEF";
  else tag = "_GH";


  for(vector<KMuon>::iterator itmu=mu.begin(); itmu!=mu.end(); ++itmu) {
    float mupt=itmu->MiniAODPt();
    if(itmu->MiniAODPt() >120. )mupt=119.;
    if(itmu->MiniAODPt() < 20.) mupt=21.;


    if(CheckCorrectionHist("ISO"+tag +"_"+ muid)){

      sferr = double(sys)*GetCorrectionHist("ISO"+tag +"_"+ muid)->GetBinError( GetCorrectionHist("ISO"+tag +"_"+ muid)->FindBin( fabs(itmu->Eta()), mupt) );
      
      sf*= (1. + sferr)*GetCorrectionHist("ISO"+tag +"_"+ muid)->GetBinContent( GetCorrectionHist("ISO"+tag +"_"+ muid)->FindBin( fabs(itmu->Eta()), mupt) );
    }
  }

  return sf;
}




double MCDataCorrections::MuonScaleFactor(TString muid, vector<snu::KMuon> mu,int sys){

  if(k_period < 0) {
    /// If k_period < 0 then using ALL data periods and use weighted SF                                                                                                                                                                                                        

    double lumi_periodB = 5.929001722;
    double lumi_periodC = 2.645968083;
    double lumi_periodD = 4.35344881;
    double lumi_periodE = 4.049732039;
    double lumi_periodF = 3.157020934;
    double lumi_periodG = 7.549615806;
    double lumi_periodH = 8.545039549 + 0.216782873;
    double total_lumi = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG+lumi_periodH);

    double WeightBtoF = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF)/total_lumi;
    double WeightGtoH = (lumi_periodG+lumi_periodH)/total_lumi;
    double SF_bf= MuonScaleFactorPeriodDependant(muid,mu, 1, sys);
    double SF_gh= MuonScaleFactorPeriodDependant(muid,mu, 7, sys);

    double SF_weight = WeightBtoF*SF_bf + WeightGtoH*SF_gh;
    return SF_weight;
  }
  return MuonScaleFactorPeriodDependant(muid, mu, k_period, sys);

}


double MCDataCorrections::MuonScaleFactorPeriodDependant(TString muid, vector<snu::KMuon> mu,int cat_period, int sys){
  float sf= 1.;
  float sferr=1.;
  if(corr_isdata) return 1.;
  
  // ref https://indico.cern.ch/event/595070/contributions/2405095/attachments/1388822/2114847/MC_12_12_2016.pdf

  TString tag = "";
  if(cat_period < 6) tag = "_BCDEF";
  else tag = "_GH";

  double min_pt = 20., max_pt = 120.;
  
  if(mu.size() == 0) return 1.;
  for(vector<KMuon>::iterator itmu=mu.begin(); itmu!=mu.end(); ++itmu) {
    float mupt=itmu->MiniAODPt();

    if(itmu->MiniAODPt() < min_pt) mupt = min_pt+1.;
    if(itmu->MiniAODPt() >= max_pt) mupt = max_pt-1.;

    if(CheckCorrectionHist("ID" +tag+ "_"+ muid)){
      sferr = double(sys)*GetCorrectionHist("ID" +tag+ "_"+ muid)->GetBinError( GetCorrectionHist("ID" +tag+ "_"+ muid)->FindBin( fabs(itmu->Eta()), mupt) );
      
      sf*=  (1. + sferr)* GetCorrectionHist("ID" +tag+ "_"+ muid)->GetBinContent( GetCorrectionHist("ID" +tag+ "_"+ muid)->FindBin( fabs(itmu->Eta()), mupt) );
    }
  }

  return sf;
}

double MCDataCorrections::TriggerScaleFactor( vector<snu::KElectron> el, vector<snu::KMuon> mu,  TString trigname, int direction){

  if(k_period < 0) {
    /// If k_period < 0 then using ALL data periods and use weighted SF                                                                                                        
    
    double lumi_periodB = 5.929001722;
    double lumi_periodC = 2.645968083;
    double lumi_periodD = 4.35344881;
    double lumi_periodE = 4.049732039;
    double lumi_periodF = 3.157020934;
    double lumi_periodG = 7.549615806;
    double lumi_periodH = 8.545039549 + 0.216782873;
    double total_lumi = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG+lumi_periodH);

    double WeightBtoF = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF)/total_lumi;
    double WeightGtoH = (lumi_periodG+lumi_periodH)/total_lumi;
    double SF_bf= TriggerScaleFactorPeriodDependant(el, mu, trigname, direction, 1);
    double SF_gh= TriggerScaleFactorPeriodDependant(el, mu, trigname, direction, 7);


    double SF_weight = WeightBtoF*SF_bf + WeightGtoH*SF_gh;
    return SF_weight;
  }
  return TriggerScaleFactorPeriodDependant(el, mu, trigname, direction,k_period);

}


double MCDataCorrections::TriggerSFFullLumi_LegByLeg( vector<snu::KElectron> EleColl, vector<snu::KMuon> MuColl, TString TrigName, TString EleID, TString MuID, TString Option){
  //For people who doesn't consider using partial dataset and partial period SF.
  //And for whom prefers lumi weighted way to random period approach.
  //Option : 1. SystUp  : Systematic uncertainty + / SystDown : Systematic Uncertainty down (Not implemented yet)
  //         2. MCEff   : return MC efficiency of the trigger for the lepton selection
  //            DataEff : return Data efficiency of the trigger for the lepton selection
  //            (Defalt is to return scale factor, not the trigger efficiency)
  
  int SystDirection=0;
  if     (Option.Contains("SystUp"))   SystDirection=1;
  else if(Option.Contains("SystDown")) SystDirection=-1; 

  std::vector<TString> KeyVector=GetTrigSFKeys(MuID, EleID, TrigName);

  if(TrigName.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")){

    if(KeyVector.size()!=4){ return 0.;} //ID+TrigName Did not match 

    if(EleColl.size()==1 && MuColl.size()==1){
      TString Leg1DataEffSrcName=KeyVector.at(0), Leg1MCEffSrcName=KeyVector.at(1);
      TString Leg2DataEffSrcName=KeyVector.at(2), Leg2MCEffSrcName=KeyVector.at(3);

      if( !(CheckCorrectionHist(Leg1DataEffSrcName) && CheckCorrectionHist(Leg1MCEffSrcName)
            && CheckCorrectionHist(Leg2DataEffSrcName) && CheckCorrectionHist(Leg2MCEffSrcName)) ) return 0.;

      float EffMuBd_down=8., EffMuBd_up=200., EffEleBd_down=23., EffEleBd_up=200.;
      float mupt=MuColl.at(0).MiniAODPt(), elpt=EleColl.at(0).Pt();
      float mueta=EleColl.at(0).Eta(),     eleta=EleColl.at(0).Eta();

      //efficiency below threshold is too uncertain. Do not apply SF below threshold
      if(mupt<EffMuBd_down) { return 1.;} else if(mupt>EffMuBd_up){ mupt=EffMuBd_up-1.; }
      if(elpt<EffEleBd_down){ return 1.;} else if(elpt>EffEleBd_up){ elpt=EffEleBd_up-1.; }


      float MuLegEff_Data=1., MuLegEff_MC=1., EleLegEff_Data=1., EleLegEff_MC=1.;
      MuLegEff_Data  = GetCorrectionHist(Leg1DataEffSrcName)->GetBinContent( GetCorrectionHist(Leg1DataEffSrcName)->FindBin(mupt, mueta));
      MuLegEff_MC    = GetCorrectionHist(Leg1MCEffSrcName)  ->GetBinContent( GetCorrectionHist(Leg1MCEffSrcName)  ->FindBin(mupt, mueta));
      EleLegEff_Data = GetCorrectionHist(Leg2DataEffSrcName)->GetBinContent( GetCorrectionHist(Leg2DataEffSrcName)->FindBin(elpt, eleta));
      EleLegEff_MC   = GetCorrectionHist(Leg2MCEffSrcName)  ->GetBinContent( GetCorrectionHist(Leg2MCEffSrcName)  ->FindBin(elpt, eleta));

      float DataEff=1., MCEff=1.;
      DataEff=MuLegEff_Data*EleLegEff_Data;
      MCEff  =MuLegEff_MC*EleLegEff_MC;
      
      if     (Option.Contains("MCEff"))   return MCEff;
      else if(Option.Contains("DataEff")) return DataEff;

      if(MCEff==0.) return 0.;
      
      return DataEff/MCEff;

    }//Mu8Ele23 1e1mu ends
    else if(EleColl.size()==1 && MuColl.size()==2){
      TString Leg1DataEffSrcName=KeyVector.at(0), Leg1MCEffSrcName=KeyVector.at(1);
      TString Leg2DataEffSrcName=KeyVector.at(2), Leg2MCEffSrcName=KeyVector.at(3);

      if( !(  CheckCorrectionHist(Leg1DataEffSrcName) && CheckCorrectionHist(Leg1MCEffSrcName)
           && CheckCorrectionHist(Leg2DataEffSrcName) && CheckCorrectionHist(Leg2MCEffSrcName) ) ) return 0.;

      float EffMuBd_down=8., EffMuBd_up=200., EffEleBd_down=23., EffEleBd_up=200.;
      float mu1pt=MuColl.at(0).MiniAODPt(), mu2pt=MuColl.at(1).MiniAODPt(), elpt=EleColl.at(0).Pt();
      float mu1eta=EleColl.at(0).Eta(),     mu2eta=MuColl.at(1).Eta(),      eleta=EleColl.at(0).Eta();

      //efficiency below threshold is too uncertain. Do not apply SF below threshold
      if(mu1pt<EffMuBd_down){ return 1.;} else if(mu1pt>EffMuBd_up){ mu1pt=EffMuBd_up-1.; }
      if(mu2pt<EffMuBd_down){ return 1.;} else if(mu2pt>EffMuBd_up){ mu2pt=EffMuBd_up-1.; }
      if(elpt<EffEleBd_down){ return 1.;} else if(elpt>EffEleBd_up){ elpt=EffEleBd_up-1.; }


      float Mu1LegEff_Data=1., Mu1LegEff_MC=1., Mu2LegEff_Data=1., Mu2LegEff_MC=1., EleLegEff_Data=1., EleLegEff_MC=1.;
      Mu1LegEff_Data = GetCorrectionHist(Leg1DataEffSrcName)->GetBinContent( GetCorrectionHist(Leg1DataEffSrcName)->FindBin(mu1pt, mu1eta));
      Mu1LegEff_MC   = GetCorrectionHist(Leg1MCEffSrcName)  ->GetBinContent( GetCorrectionHist(Leg1MCEffSrcName)  ->FindBin(mu1pt, mu1eta));
      Mu2LegEff_Data = GetCorrectionHist(Leg1DataEffSrcName)->GetBinContent( GetCorrectionHist(Leg1DataEffSrcName)->FindBin(mu2pt, mu2eta));
      Mu2LegEff_MC   = GetCorrectionHist(Leg1MCEffSrcName)  ->GetBinContent( GetCorrectionHist(Leg1MCEffSrcName)  ->FindBin(mu2pt, mu2eta));

      EleLegEff_Data = GetCorrectionHist(Leg2DataEffSrcName)->GetBinContent( GetCorrectionHist(Leg1DataEffSrcName)->FindBin(elpt, eleta));
      EleLegEff_MC   = GetCorrectionHist(Leg2MCEffSrcName)  ->GetBinContent( GetCorrectionHist(Leg1MCEffSrcName)  ->FindBin(elpt, eleta));

      float DataEff=1., MCEff=1.;
      DataEff= EleLegEff_Data*(1.-(1.-Mu1LegEff_Data)*(1.-Mu2LegEff_Data));
      MCEff  = EleLegEff_MC*(1.-(1.-Mu1LegEff_MC)*(1.-Mu2LegEff_MC));
      
      if(MCEff==0.) return 0.;
      
      if     (Option.Contains("MCEff"))   return MCEff;
      else if(Option.Contains("DataEff")) return DataEff;

      return DataEff/MCEff;

    }//Mu8Ele23 1e2mu ends
  }//Mu8Ele23 Trigger ends
  else if(TrigName.Contains("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")){

    if(KeyVector.size()!=6){ return 0.;} //ID+TrigName Did not match for some list. Since KeyVectorSize=2*Nfilter;

    if(EleColl.size()==1 && MuColl.size()==1){
      TString Leg1DataEffSrcName=KeyVector.at(0), Leg1MCEffSrcName=KeyVector.at(1);
      TString Leg2DataEffSrcName=KeyVector.at(2), Leg2MCEffSrcName=KeyVector.at(3);
      TString DzDataEffSrcName  =KeyVector.at(4), DzMCEffSrcName  =KeyVector.at(5);

      if( !(  CheckCorrectionHist(Leg1DataEffSrcName) && CheckCorrectionHist(Leg1MCEffSrcName)
           && CheckCorrectionHist(Leg2DataEffSrcName) && CheckCorrectionHist(Leg2MCEffSrcName)
           && CheckCorrectionHist(DzDataEffSrcName)   && CheckCorrectionHist(DzMCEffSrcName)  ) ) return 0.;

      float EffMuBd_down=12., EffMuBd_up=200., EffEleBd_down=23., EffEleBd_up=200.;
      float mupt=MuColl.at(0).MiniAODPt(), elpt=EleColl.at(0).Pt();
      float mueta=EleColl.at(0).Eta(),     eleta=EleColl.at(0).Eta();

      //efficiency below threshold is too uncertain. Do not apply SF below threshold
      if(mupt<EffMuBd_down) { return 1.;} else if(mupt>EffMuBd_up){ mupt=EffMuBd_up-1.; }
      if(elpt<EffEleBd_down){ return 1.;} else if(elpt>EffEleBd_up){ elpt=EffEleBd_up-1.; }


      float MuLegEff_Data=1., MuLegEff_MC=1., EleLegEff_Data=1., EleLegEff_MC=1., DzEff_Data=1., DzEff_MC=1.;
      MuLegEff_Data  = GetCorrectionHist(Leg1DataEffSrcName)->GetBinContent( GetCorrectionHist(Leg1DataEffSrcName)->FindBin(mupt, mueta));
      MuLegEff_MC    = GetCorrectionHist(Leg1MCEffSrcName)  ->GetBinContent( GetCorrectionHist(Leg1MCEffSrcName)  ->FindBin(mupt, mueta));
      EleLegEff_Data = GetCorrectionHist(Leg2DataEffSrcName)->GetBinContent( GetCorrectionHist(Leg2DataEffSrcName)->FindBin(elpt, eleta));
      EleLegEff_MC   = GetCorrectionHist(Leg2MCEffSrcName)  ->GetBinContent( GetCorrectionHist(Leg2MCEffSrcName)  ->FindBin(elpt, eleta));
      DzEff_Data     = GetCorrectionHist(DzDataEffSrcName)  ->GetBinContent( GetCorrectionHist(DzDataEffSrcName)  ->FindBin(mueta, eleta));
      DzEff_MC       = GetCorrectionHist(DzMCEffSrcName)    ->GetBinContent( GetCorrectionHist(DzMCEffSrcName)    ->FindBin(mueta, eleta));


      float DataEff=1., MCEff=1.;
      DataEff = MuLegEff_Data*EleLegEff_Data*DzEff_Data;
      MCEff   = MuLegEff_MC*EleLegEff_MC*DzEff_MC;
      
      if(MCEff==0.) return 0.;

      if     (Option.Contains("MCEff"))   return MCEff;
      else if(Option.Contains("DataEff")) return DataEff;
      
      return DataEff/MCEff;

    }//Mu12Ele23Dz 1e1mu ends
    else if(EleColl.size()==1 && MuColl.size()==2){
      TString Leg1DataEffSrcName=KeyVector.at(0), Leg1MCEffSrcName=KeyVector.at(1);
      TString Leg2DataEffSrcName=KeyVector.at(2), Leg2MCEffSrcName=KeyVector.at(3);
      TString DzDataEffSrcName  =KeyVector.at(4), DzMCEffSrcName  =KeyVector.at(5);

      if( !(  CheckCorrectionHist(Leg1DataEffSrcName) && CheckCorrectionHist(Leg1MCEffSrcName)
           && CheckCorrectionHist(Leg2DataEffSrcName) && CheckCorrectionHist(Leg2MCEffSrcName)
           && CheckCorrectionHist(DzDataEffSrcName)   && CheckCorrectionHist(DzMCEffSrcName)   ) ) return 0.;

      float EffMuBd_down=12., EffMuBd_up=200., EffEleBd_down=23., EffEleBd_up=200.;
      float mu1pt=MuColl.at(0).MiniAODPt(), mu2pt=MuColl.at(1).MiniAODPt(), elpt=EleColl.at(0).Pt();
      float mu1eta=EleColl.at(0).Eta(),     mu2eta=MuColl.at(1).Eta(),      eleta=EleColl.at(0).Eta();

      //efficiency below threshold is too uncertain. Do not apply SF below threshold
      if(mu1pt<EffMuBd_down){ return 1.;} else if(mu1pt>EffMuBd_up){ mu1pt=EffMuBd_up-1.; }
      if(mu2pt<EffMuBd_down){ return 1.;} else if(mu2pt>EffMuBd_up){ mu2pt=EffMuBd_up-1.; }
      if(elpt<EffEleBd_down){ return 1.;} else if(elpt>EffEleBd_up){ elpt=EffEleBd_up-1.; }


      float Mu1LegEff_Data=1., Mu1LegEff_MC=1., Mu2LegEff_Data=1., Mu2LegEff_MC=1., EleLegEff_Data=1., EleLegEff_MC=1.;
      float DzEffMu1Ele_Data=1., DzEffMu1Ele_MC=1., DzEffMu2Ele_Data=1., DzEffMu2Ele_MC=1.;
      Mu1LegEff_Data = GetCorrectionHist(Leg1DataEffSrcName)->GetBinContent( GetCorrectionHist(Leg1DataEffSrcName)->FindBin(mu1pt, mu1eta));
      Mu1LegEff_MC   = GetCorrectionHist(Leg1MCEffSrcName)  ->GetBinContent( GetCorrectionHist(Leg1MCEffSrcName)  ->FindBin(mu1pt, mu1eta));
      Mu2LegEff_Data = GetCorrectionHist(Leg1DataEffSrcName)->GetBinContent( GetCorrectionHist(Leg1DataEffSrcName)->FindBin(mu2pt, mu2eta));
      Mu2LegEff_MC   = GetCorrectionHist(Leg1MCEffSrcName)  ->GetBinContent( GetCorrectionHist(Leg1MCEffSrcName)  ->FindBin(mu2pt, mu2eta));

      EleLegEff_Data = GetCorrectionHist(Leg2DataEffSrcName)->GetBinContent( GetCorrectionHist(Leg1DataEffSrcName)->FindBin(elpt, eleta));
      EleLegEff_MC   = GetCorrectionHist(Leg2MCEffSrcName)  ->GetBinContent( GetCorrectionHist(Leg1MCEffSrcName)  ->FindBin(elpt, eleta));

      DzEffMu1Ele_Data= GetCorrectionHist(DzDataEffSrcName) ->GetBinContent( GetCorrectionHist(DzDataEffSrcName) ->FindBin(mu1eta, eleta));
      DzEffMu1Ele_MC  = GetCorrectionHist(DzMCEffSrcName)   ->GetBinContent( GetCorrectionHist(DzMCEffSrcName)   ->FindBin(mu1eta, eleta));
      DzEffMu2Ele_Data= GetCorrectionHist(DzDataEffSrcName) ->GetBinContent( GetCorrectionHist(DzDataEffSrcName) ->FindBin(mu2eta, eleta));
      DzEffMu2Ele_MC  = GetCorrectionHist(DzMCEffSrcName)   ->GetBinContent( GetCorrectionHist(DzMCEffSrcName)   ->FindBin(mu2eta, eleta));



      float DataEff=1., MCEff=1.;
      DataEff= EleLegEff_Data*(1.-(1.-Mu1LegEff_Data*DzEffMu1Ele_Data)*(1.-Mu2LegEff_Data*DzEffMu2Ele_Data));
      MCEff  = EleLegEff_MC*(1.-(1.-Mu1LegEff_MC*DzEffMu1Ele_MC)*(1.-Mu2LegEff_MC*DzEffMu2Ele_MC));

      if     (Option.Contains("MCEff"))   return MCEff;
      else if(Option.Contains("DataEff")) return DataEff;
      
      if(MCEff==0.) return 0.;
      
      return DataEff/MCEff;

    }//Mu12Ele23Dz 1e2mu ends
  }//Mu12Ele23Dz Trigger ends
  else if(TrigName.Contains("HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v")){
  
    if(MuID.Contains("HN_TRI_TIGHT")){
      //Lumi uncertainty is 2~3%. There's no point calculating with more numbers
      float Lumi_BtoF=20.13, Lumi_GtoH=16.31, TotalLumi=36.44;
      float SF_bf=TriggerScaleFactorPeriodDependant(EleColl, MuColl, TrigName, SystDirection, 1);
      float SF_gh=TriggerScaleFactorPeriodDependant(EleColl, MuColl, TrigName, SystDirection, 7);
  
      float SF_weight = SF_bf*Lumi_BtoF/TotalLumi + SF_gh*Lumi_GtoH/TotalLumi;


      return SF_weight;
    }

  }
    
  return 1.;

}




std::vector<TString> MCDataCorrections::GetTrigSFKeys(TString MuID, TString EleID, TString TrigName){
  
  //Types of SF related data: 1. DataEff, 2. MCEff, 3. SF (Sigl Lep Trigger)
  TString LabelEffData="", LabelEffMC="", LabelSF="";

  //Types of SF related data: 1. Leg1Data, 2. Leg1MC, 3. Leg2Data 4. Leg2MC, 5. DzData 6. DzMC (DiLep Trigger)
  TString LabelEffLeg1Data="", LabelEffLeg2Data="", LabelEffLeg1MC="", LabelEffLeg2MC="", LabelEffDzData="", LabelEffDzMC="";
  std::vector<TString> KeyVectors;

  bool SiglLepTrig=false, DiLepTrig=false;

  if(TrigName.Contains("HLT_IsoMu24") || TrigName.Contains("HLT_IsoTkMu24")){
    if(MuID=="MUON_POG_TIGHT"){
      LabelEffData="MUON_MU24_TRIGGER_DATA_EFF", LabelEffMC="MUON_MU24_TRIGGER_MC_EFF", LabelSF="MUON_MU24_TRIGGER";
      SiglLepTrig=true;
    }
  }
  else if(TrigName.Contains("HLT_Ele25_eta2p1_WPTight_Gsf")){
    if(EleID=="ELECTRON_POG_TIGHT"){
      LabelSF="ELECTRON_ELE25_TRIGGER_SF";
      SiglLepTrig=true;
    }
  }
  else if(TrigName.Contains("HLT_Ele27_WPTight_Gsf")){
    if(EleID=="ELECTRON_POG_TIGHT"){
      LabelSF="ELECTRON_ELE27_TRIGGER_SF";
      SiglLepTrig=true;
    }
  }
  else if(TrigName.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")){
    if(EleID=="ELECTRON_HWA_TEST1" && MuID=="HN_TRI_TEST1"){
      LabelEffLeg2Data="MU8ELE23_TRIGGEREFF_ELE23_DATA", LabelEffLeg2MC="MU8ELE23_TRIGGEREFF_ELE23_MC";
      LabelEffLeg1Data="MU8ELE23_TRIGGEREFF_MU8_DATA",   LabelEffLeg1MC="MU8ELE23_TRIGGEREFF_MU8_MC";
      DiLepTrig=true;
    }
  }
  else if(TrigName.Contains("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v")){
    if(EleID=="ELECTRON_HWA_TEST1" && MuID=="HN_TRI_TEST1"){
      LabelEffLeg2Data="MU8ELE23_TRIGGEREFF_ELE23_DATA", LabelEffLeg2MC="MU8ELE23_TRIGGEREFF_ELE23_MC";
      LabelEffLeg1Data="MU8ELE23_TRIGGEREFF_MU8_DATA",   LabelEffLeg1MC="MU8ELE23_TRIGGEREFF_MU8_MC";
      LabelEffDzData="MU8ELE23_TRIGGEREFF_DZ_DATA",      LabelEffDzMC="MU8ELE23_TRIGGEREFF_DZ_MC";
      DiLepTrig=true;
    }
  }


  if(SiglLepTrig){
    KeyVectors.push_back(LabelEffData);
    KeyVectors.push_back(LabelEffMC);
    KeyVectors.push_back(LabelSF);
  }
  else if(DiLepTrig){
    KeyVectors.push_back(LabelEffLeg1Data);
    KeyVectors.push_back(LabelEffLeg1MC);
    KeyVectors.push_back(LabelEffLeg2Data);
    KeyVectors.push_back(LabelEffLeg2MC);
    KeyVectors.push_back(LabelEffDzData);
    KeyVectors.push_back(LabelEffDzMC);
  }

  return KeyVectors;

}



double MCDataCorrections::TriggerScaleFactorPeriodDependant( vector<snu::KElectron> el, vector<snu::KMuon> mu,  TString trigname, int cat_period, int direction){

  /// Currently only muon scale factors (single) are added by pog
  
  if(corr_isdata) return 1.;

  TString s_ptthreshold="";
  float f1_ptthreshold=-999.;
  float f2_ptthreshold=-999.;
  TString DataEffSrcName="", MCEffSrcName="";

  
  TString tag = "";
  if(cat_period < 6) tag = "_BCDEF";
  else tag = "_GH";
  

  if(trigname == "HLT_Ele25_eta2p1_WPTight_Gsf" || trigname == "HLT_Ele27_WPTight_Gsf" || trigname == "HLT_Ele32_eta2p1_WPTight_Gsf"){
    //https://indico.cern.ch/event/604912/contributions/2490011/attachments/1418869/2173471/2017.02.27_EGM_Ele25-and-Ele27-trigger-SF_v1.pdf
    // HLT_Ele27_WPTight_Gsf and HLT_Ele25_eta2p1_WPTight_Gsf are unprescaled
    //https://indico.cern.ch/event/604911/contributions/2474230/attachments/1411693/2159538/2017.02.13_EGM_Electron-trigger-SF-2016_v1.pdf
    if(trigname == "HLT_Ele32_eta2p1_WPTight_Gsf") {
      s_ptthreshold="32";
      f1_ptthreshold = 35.;
    }
    else if(trigname == "HLT_Ele25_eta2p1_WPTight_Gsf"){
      s_ptthreshold="25";
      f1_ptthreshold = 25.;
    }
    else {
      s_ptthreshold="27";
      f1_ptthreshold = 28.;
    }
    f2_ptthreshold = 200.;
    
    for(unsigned int iel = 0; iel < el.size(); iel++){
      if(!el.at(iel).TriggerMatched(trigname)) continue;
      float elpt = el.at(iel).Pt();
      if(elpt >  f2_ptthreshold) elpt = (f2_ptthreshold-1.);
      if(elpt < f1_ptthreshold) return 1.;
      if(CheckCorrectionHist(("ELECTRON_ELE"+s_ptthreshold+"_TRIGGER_SF" + tag).Data())){
	float sferr = double(direction)*GetCorrectionHist(("ELECTRON_ELE"+s_ptthreshold+"_TRIGGER_SF" + tag).Data())->GetBinError(GetCorrectionHist(("ELECTRON_ELE"+s_ptthreshold+"_TRIGGER_SF" + tag).Data())->FindBin( elpt,  fabs(el.at(iel).Eta())));
														 
	return  (1. + sferr)*GetCorrectionHist(("ELECTRON_ELE"+s_ptthreshold+"_TRIGGER_SF" + tag).Data())->GetBinContent(GetCorrectionHist(("ELECTRON_ELE"+s_ptthreshold+"_TRIGGER_SF" + tag).Data())->FindBin( elpt,  fabs(el.at(0).Eta())));
      }
    }// el loop
    return 1.;
  }
  if(trigname == "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"){
    // AN-16-172 HWW analysis (uses POGTight ID)

    

  }
  if(trigname.Contains("HLT_IsoMu") || trigname.Contains("HLT_IsoTkMu")){
    if(mu.size()==1){
      if (trigname.Contains("24")){
	  s_ptthreshold = "24";
	  f1_ptthreshold = 24.;
	  f2_ptthreshold = 500.;
      }      
      else if (trigname.Contains("50")){
	s_ptthreshold= "50";
	f1_ptthreshold= 50.;
	f2_ptthreshold= 800.;
      }

      // G+H    https://twiki.cern.ch/twiki/pub/CMS/MuonWorkInProgressAndPagResults/2016.12.14_MuonPOGTriggerSF_KPLee_v1.pdf 
      // BtoF https://indico.cern.ch/event/608200/contributions/2452382/attachments/1401679/2156728/2017.02.05_MuonPOGTriggerSF_KPLee_v1.pdf 
      
      float mupt=mu.at(0).MiniAODPt();
      if(mupt >  f2_ptthreshold) mupt = (f2_ptthreshold-1.);
      if(mupt < f1_ptthreshold) return 1.;
      if(CheckCorrectionHist(("MUON_MU"+s_ptthreshold+"_TRIGGER"+tag).Data())){
	float sferr = double(direction)*GetCorrectionHist(("MUON_MU"+s_ptthreshold+"_TRIGGER"+tag).Data())->GetBinError(GetCorrectionHist(("MUON_MU"+s_ptthreshold+"_TRIGGER"+tag).Data())->FindBin( fabs(mu.at(0).Eta()), mupt) );
	return  (1. + sferr)*GetCorrectionHist(("MUON_MU"+s_ptthreshold+"_TRIGGER"+tag).Data())->GetBinContent( GetCorrectionHist(("MUON_MU"+s_ptthreshold+"_TRIGGER"+tag).Data())->FindBin(  fabs(mu.at(0).Eta()), mupt) );
      }
    }//One Muon case to SinglMuTrig Ends
    else if(mu.size()>=2){
      if(trigname.Contains("24")){
        s_ptthreshold  ="24";
        f1_ptthreshold = 24.;
        f2_ptthreshold = 500.;
      }
      DataEffSrcName = "MUON_MU"+s_ptthreshold+"_TRIGGER_DATA_EFF"+tag;
      MCEffSrcName   = "MUON_MU"+s_ptthreshold+"_TRIGGER_MC_EFF"+tag;
      if( !( CheckCorrectionHist(DataEffSrcName) && CheckCorrectionHist(MCEffSrcName) ) ) return 0.;

      std::vector<float> muptColl;
        for(unsigned int i=0; i<mu.size(); i++){
          float mupt=mu.at(i).MiniAODPt();
          if     (mupt>f2_ptthreshold){ muptColl.push_back(f2_ptthreshold - 1.); }
          else if(mupt<f1_ptthreshold){ muptColl.push_back(f1_ptthreshold - 1.); }
          else                          muptColl.push_back(mupt);
        }
      std::vector<float> muetaColl;
      for(unsigned int i=0; i<mu.size(); i++){ muetaColl.push_back(mu.at(i).Eta()); };

      std::vector<float> dataeffColl;
      for(unsigned int i=0; i<mu.size(); i++){ dataeffColl.push_back( GetCorrectionHist(DataEffSrcName)->GetBinContent( GetCorrectionHist(DataEffSrcName)->FindBin(muptColl.at(i), fabs(muetaColl.at(i))) ) ); };

      std::vector<float> mceffColl;
      for(unsigned int i=0; i<mu.size(); i++){ mceffColl.push_back( GetCorrectionHist(MCEffSrcName)->GetBinContent( GetCorrectionHist(MCEffSrcName)->FindBin(muptColl.at(i), fabs(muetaColl.at(i))) ) ); };
        
      float DataFailProb =1.; for(int i=0; i<mu.size(); i++){ DataFailProb *= (1.-dataeffColl.at(i)); };
      float MCFailProb   =1.; for(int i=0; i<mu.size(); i++){ MCFailProb   *= (1.-mceffColl.at(i));   };
      
      if((1.-MCFailProb)==0.) return 0.;
      float SF = (1.-DataFailProb)/(1.-MCFailProb);

      return SF;
      

    }//2 Muon case to SinglMuTrig ends.
  }//SingleMuTrig case ends.

  return 1.;


}



double  MCDataCorrections::GetDoubleEGTriggerEff(vector<snu::KElectron> el){
  
  if(corr_isdata) return 1.;
  

  // https://twiki.cern.ch/twiki/pub/CMS/HWW2016TriggerAndIdIsoScaleFactorsResults/AN-16-172_temp.pdf
  if(el.size() <  2.) return 0.;
  if(el.size() == 2.){
    double eff_tl = GetEffDEG1(el[0]) * GetEffDEG2(el[1]);
    double eff_lt = GetEffDEG2(el[0]) * GetEffDEG1(el[1]);
    
    double evt_eff = eff_tl + (1. - eff_lt) * eff_tl;
    return evt_eff;
  }
  /// This is approx correct to set nel > 2 as same form for nel=2 (Fix for correct form. when time permitss)
  double eff_tl = GetEffDEG1(el[0]) * GetEffDEG2(el[1]);
  double eff_lt = GetEffDEG2(el[0]) * GetEffDEG1(el[1]);

  double evt_eff = eff_tl + (1. - eff_lt) * eff_tl;
  return evt_eff;


  return 1.;
}


double MCDataCorrections::GetEffDEG1(snu::KElectron el){

  float el_eta = el.Eta();
  float el_pt = el.Pt();
  if(el_pt < 10.) return 0.;
  if(el_pt >= 100.) el_pt = 99.;

  if(fabs(el_eta ) > 2.5) return 0.;

  vector<float> etabins;
  etabins.push_back(2.1);
  etabins.push_back(1.6);
  etabins.push_back(1.4);
  etabins.push_back(0.8);
  etabins.push_back(0.);
  etabins.push_back(-0.8);
  etabins.push_back(-1.4);
  etabins.push_back(-1.6);
  etabins.push_back(-2.1);
  etabins.push_back(-2.5);

  vector<float> ptbins;
  ptbins.push_back(23.);
  ptbins.push_back(24.);
  ptbins.push_back(25.);
  ptbins.push_back(26.);
  ptbins.push_back(30.);
  ptbins.push_back(35.);
  ptbins.push_back(40.);
  ptbins.push_back(45.);
  ptbins.push_back(50.);
  ptbins.push_back(100.);

  for(unsigned int i = 0 ; i <  etabins.size(); i++){
    if(el_eta  > etabins[i]) {
      std::map<float, std::vector<float>* >::iterator it = deg_etaptmap_leg1.find(etabins[i]);

      if(it == deg_etaptmap_leg1.end()){
        cout << "Error in map for double 1 eg trigger SF" << endl;
        exit(EXIT_FAILURE);
      }
      for(unsigned int ipt = 0; ipt < ptbins.size(); ipt++){
        if(el_pt < ptbins[ipt]) return it->second->at(ipt);
      }
    }
  }
  return 1.;
}



double MCDataCorrections::GetEffDEG2(snu::KElectron el){
  float el_eta = el.Eta();
  float el_pt = el.Pt();
  if(el_pt < 10.) return 0.;
  if(el_pt >= 100.) el_pt = 99.;

  if(fabs(el_eta ) > 2.5) return 0.;

  vector<float> etabins;
  etabins.push_back(2.1);
  etabins.push_back(1.6);
  etabins.push_back(1.4);
  etabins.push_back(0.8);
  etabins.push_back(0.);
  etabins.push_back(-0.8);
  etabins.push_back(-1.4);
  etabins.push_back(-1.6);
  etabins.push_back(-2.1);
  etabins.push_back(-2.5);
  
  // low pt leg
  vector<float> ptbins;
  ptbins.push_back(12.);
  ptbins.push_back(13.);
  ptbins.push_back(15.);
  ptbins.push_back(18.);
  ptbins.push_back(22.);
  ptbins.push_back(30.);
  ptbins.push_back(35.);
  ptbins.push_back(40.);
  ptbins.push_back(50.);
  ptbins.push_back(100.);

  for(unsigned int i = 0 ; i <  etabins.size(); i++){
    if(el_eta > etabins[i]) {
      std::map<float, std::vector<float>* >::iterator it = deg_etaptmap_leg2.find(etabins[i]);
      if(it == deg_etaptmap_leg2.end()){
        cout << "Error in map for double eg 2 trigger SF" << endl;
        exit(EXIT_FAILURE);
      }
      for(unsigned int ipt = 0; ipt < ptbins.size(); ipt++){
        if(el_pt < ptbins[ipt]) return it->second->at(ipt);
      }
    }
  }
  return 1.;
}

//==== Trigger Scale Factor, by "leg-by-leg" strategy

double MCDataCorrections::TriggerEfficiencyLegByLeg(std::vector<snu::KElectron> el, std::vector<snu::KMuon> mu, int TriggerCategory, int DataOrMC, int direction){

  if(k_period < 0) {
    /// If k_period < 0 then using ALL data periods and use weighted SF                                                                                                        

    double lumi_periodB = 5.929001722;
    double lumi_periodC = 2.645968083;
    double lumi_periodD = 4.35344881;
    double lumi_periodE = 4.049732039;
    double lumi_periodF = 3.157020934;
    double lumi_periodG = 7.549615806;
    double lumi_periodH = 8.545039549 + 0.216782873;
    double total_lumi = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG+lumi_periodH);

    double WeightBtoF = (lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF)/total_lumi;
    double WeightGtoH = (lumi_periodG+lumi_periodH)/total_lumi;
    double Eff_bf = TriggerEfficiencyLegByLegPeriodDependant(el, mu, TriggerCategory, 1, DataOrMC, direction);
    double Eff_gh = TriggerEfficiencyLegByLegPeriodDependant(el, mu, TriggerCategory, 7, DataOrMC, direction);

    double Eff_weight = WeightBtoF*Eff_bf + WeightGtoH*Eff_gh;
    //cout << "[MCDataCorrections::TriggerEfficiencyLegByLeg] WeightBtoF = " << WeightBtoF << endl;
    //cout << "[MCDataCorrections::TriggerEfficiencyLegByLeg] WeightGtoH = " << WeightGtoH << endl;
    //cout << "[MCDataCorrections::TriggerEfficiencyLegByLeg] Eff_bf = " << Eff_bf << endl;
    //cout << "[MCDataCorrections::TriggerEfficiencyLegByLeg] Eff_gh = " << Eff_gh << endl;

    return Eff_weight;
  }
  return TriggerEfficiencyLegByLegPeriodDependant(el, mu, TriggerCategory, k_period, DataOrMC, direction);

}

double MCDataCorrections::TriggerEfficiencyLegByLegPeriodDependant(std::vector<snu::KElectron> el,std::vector<snu::KMuon> mu, int TriggerCategory, int catperiod, int DataOrMC, int direction){

  //==== Now, only for tri MUON case..
  //==== 1) TriggerCategory = 0
  //==== - HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v
  //==== - HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v
  //==== ==> Leg1 : Mu17_TrkIsoVVL
  //====     Leg2 : Mu8_TrkIsoVVL or TkMu8_TrkIsoVVL

  if(TriggerCategory==0){

    if(mu.size()<2) return 1.;

    TString leg1 = "MU17";
    TString leg2 = "MU8_OR_TKMU8";

    double faileff(1.);

    for(unsigned int i=0; i<mu.size()-1; i++){
      snu::KMuon mu1 = mu.at(i);
      for(unsigned j=i+1; j<mu.size(); j++){
        snu::KMuon mu2 = mu.at(j);
        double dimueff = TriggerEfficiency_DiMuon_passing_DoubleMuonTrigger(mu1, mu2, "MU17", "MU8_OR_TKMU8", DataOrMC, catperiod);
        //cout << "[MCDataCorrections::TriggerEfficiencyLegByLegPeriodDependant] dimueff = " << dimueff << endl;
        
        faileff *= (1.-dimueff);
      }
    }
    bool debug(false);
    if(debug){cout << "Direction = " << direction << "n_el " << el.size() << endl;
}
    return 1.-faileff;

  }
  else{
    return 1.;
  }

}

double MCDataCorrections::TriggerEfficiency_DiMuon_passing_DoubleMuonTrigger(snu::KMuon mu1, snu::KMuon mu2, TString leg1, TString leg2, int DataOrMC, int catperiod){

  TString tag = "";
  if(catperiod < 6) tag = "_BCDEF";
  else tag = "_GH";

  TString sample="";
  if(DataOrMC==0) sample = "_Data";
  else sample = "_MC";

  double eta1 = abs(mu1.Eta());
  double pt1 = mu1.Pt();
  if(pt1>=120.) pt1 = 119.;
  if(pt1<10.) pt1 = 10.1;
  double eta2 = abs(mu2.Eta());
  double pt2 = mu2.Pt();
  if(pt2>120.) pt2 = 119.;
  if(pt2<10.) pt2 = 10.1;

  //cout << "[MCDataCorrections::TriggerEfficiency_DiMuon_passing_DoubleMuonTrigger] tag = " << tag << ", sample = " << sample << endl;

  if( (leg1=="MU17" && leg2=="MU8_OR_TKMU8") || (leg2=="MU17" && leg1=="MU8_OR_TKMU8") ){

    TH2F *hist_leg1 = GetCorrectionHist("MUON_"+leg1+"_TRIGGER"+tag+sample);
    TH2F *hist_leg2 = GetCorrectionHist("MUON_"+leg2+"_TRIGGER"+tag+sample);

    double eff_mu1leg1 = hist_leg1->GetBinContent( hist_leg1->FindBin(eta1,pt1) );
    double eff_mu2leg2 = hist_leg2->GetBinContent( hist_leg2->FindBin(eta2,pt2) );

    double eff_mu1leg2 = hist_leg2->GetBinContent( hist_leg2->FindBin(eta1,pt1) );
    double eff_mu2leg1 = hist_leg1->GetBinContent( hist_leg1->FindBin(eta2,pt2) );

    //cout << "[MCDataCorrections::TriggerEfficiency_DiMuon_passing_DoubleMuonTrigger] pt1 = " << pt1 << ", eta1 = " << eta1 << endl;
    //cout << "[MCDataCorrections::TriggerEfficiency_DiMuon_passing_DoubleMuonTrigger] => " << leg1 << " : " << eff_mu1leg1 << endl;
    //cout << "[MCDataCorrections::TriggerEfficiency_DiMuon_passing_DoubleMuonTrigger] => " << leg2 << " : " << eff_mu1leg2 << endl;
    //cout << "[MCDataCorrections::TriggerEfficiency_DiMuon_passing_DoubleMuonTrigger] pt2 = " << pt2 << ", eta2 = " << eta2 << endl;
    //cout << "[MCDataCorrections::TriggerEfficiency_DiMuon_passing_DoubleMuonTrigger] => " << leg1 << " : " << eff_mu2leg1 << endl;
    //cout << "[MCDataCorrections::TriggerEfficiency_DiMuon_passing_DoubleMuonTrigger] => " << leg2 << " : " << eff_mu2leg2 << endl;

    double eff = 1.-(1.-eff_mu1leg1*eff_mu2leg2)*(1.-eff_mu1leg2*eff_mu2leg1); 

    return eff;
        

  } 
  else{
    return 1.;
  }
}




double MCDataCorrections::ElectronScaleFactor( TString elid, vector<snu::KElectron> el, int sys){
  float sf= 1.;
  if(corr_isdata) return 1.;

  //http://fcouderc.web.cern.ch/fcouderc/EGamma/scaleFactors/Moriond17/approval/EleID/passingVeto80X/egammaEffi.txt_egammaPlots.pdf
  std::string sid= "";
  
  if(elid.Contains("HN")) elid = "ELECTRON_MVA_80";
  if(elid.Contains("ELECTRON_MVA_90")) elid = "ELECTRON_MVA_80";

  for(vector<KElectron>::iterator itel=el.begin(); itel!=el.end(); ++itel) {
    float elpt=itel->Pt();
    if(elpt > 500.) elpt= 499.;
    if(elpt < 10.) elpt= 11;
    float unc = 0.02; //// Check this
    
    if(CheckCorrectionHist("ID_" + elid)){
      int bin =  GetCorrectionHist("ID_" + elid)->FindBin(fabs(itel->SCEta()), elpt);
      sf *= GetCorrectionHist("ID_" + elid)->GetBinContent(bin);
      float err =  GetCorrectionHist("EL_RECO")->GetBinError(bin);
      err = sqrt (pow(err, 2.) + pow(unc, 2.));
      if(sys == 1)sf *= (1. + err);
      if(sys == -1)sf *= (1. - err);
    }
  }

 
  return sf;
}

double MCDataCorrections::ElectronRecoScaleFactor(vector<snu::KElectron> el, int sys){
  
  // https://indico.cern.ch/event/604907/contributions/2452907/attachments/1401460/2139067/RecoSF_ApprovalMoriond17_25Jan2017.pdf

  if(corr_isdata) return 1.;
  float sf= 1.;
  for(vector<KElectron>::iterator itel=el.begin(); itel!=el.end(); ++itel) {
    float elpt= itel->Pt() ;
    float unc = 0.;
    if(elpt < 20.) unc=0.01;
    if(elpt >  80.) unc=0.01;
    if(itel->Pt() > 500.) elpt=499.;
    if(itel->Pt() < 25.) elpt=25.;
    if(CheckCorrectionHist("EL_RECO")){
      int bin =  GetCorrectionHist("EL_RECO")->FindBin(fabs(itel->SCEta()), elpt);
      sf *= GetCorrectionHist("EL_RECO")->GetBinContent(bin);
      if(fabs(sys) == 1){
	float err =  GetCorrectionHist("EL_RECO")->GetBinError(bin);
	err = sqrt (pow(err, 2.) + pow(unc, 2.));
	if(sys == 1)sf *= (1. + err);
	if(sys == -1)sf *= (1. - err);
      }
    }
  }
  
  return sf;
}


float MCDataCorrections::UserPileupWeight(snu::KEvent ev, int nj){
  
  if(corr_isdata) return 1.;
  return reweightPU->GetWeight(ev.nVertices(),TString(getenv("CATVERSION")), nj);
}


float MCDataCorrections::PileupWeightByPeriod(snu::KEvent ev){
  
  if(k_period < 0) return CatPileupWeight(ev,0);
  return ev.PeriodPileUpWeight(k_period);
  
}

float MCDataCorrections::CatPileupWeight(snu::KEvent ev, int syst){
  snu::KEvent::syst_dir sdir= snu::KEvent::central;
  if(syst < 0) sdir = snu::KEvent::down;
  if(syst > 0) sdir = snu::KEvent::up;
  
  return ev.PileUpWeight(sdir);

}



void MCDataCorrections::CheckFile(TFile* file){

  if(!file) exit(EXIT_FAILURE);
  return;
}


bool MCDataCorrections::CheckCorrectionHist(TString label){
  map<TString, TH2F*>::iterator  mapit = CorrectionMap.find(label);
  if (mapit!= CorrectionMap.end()){
    return true;
  }
  else return false;
}
bool MCDataCorrections::CheckCorrectionGraph(TString label){
  map<TString, TGraphAsymmErrors*>::iterator  mapit = CorrectionMapGraph.find(label);
  if (mapit!= CorrectionMapGraph.end()){
    return true;
  }
  else return false;
}

TGraphAsymmErrors* MCDataCorrections::GetCorrectionGraph(TString label){
  map<TString, TGraphAsymmErrors*>::iterator mapit = CorrectionMapGraph.find(label);
  if (mapit!= CorrectionMapGraph.end()){
    return mapit->second;
  }
  else{
    exit(0);
  }
}

TH2F* MCDataCorrections::GetCorrectionHist(TString label){
  map<TString, TH2F*>::iterator mapit = CorrectionMap.find(label);
  if (mapit!= CorrectionMap.end()){
    return mapit->second;
  }
  else{
    exit(0);
  }
}



void MCDataCorrections::CorrectMuonMomentum(vector<snu::KMuon>& k_muons, vector<snu::KTruth> truth){
  

  /// This function only works to correct muons if rochester correction was not applied in  selection code

  for(std::vector<snu::KMuon>::iterator it = k_muons.begin(); it != k_muons.end(); it++){
    double scalefactor = 1.;
    if(it->IsRochesterCorrected()) return;
    if (corr_isdata) scalefactor = rc->kScaleDT(float(it->Charge()), it->Pt(), it->Eta(), it->Phi(),0,0);
    else {
      //gRandom->SetSeed(1111.);
	double u1 = gRandom->Rndm();
	double u2 = gRandom->Rndm();

	unsigned int mu_index = it->MCTruthIndex();
	float genpt(-999.);
	if(mu_index > 0 && mu_index < truth.size()) {
	  if(fabs(truth.at(mu_index).PdgId() ) == 13) genpt = truth.at(mu_index).Pt();
	}
	

	if ( genpt> 0.)  scalefactor = rc->kScaleFromGenMC(float(it->Charge()), it->Pt(), it->Eta(), it->Phi(), it->ActiveLayer(), genpt, u1,0, 0);
	else scalefactor = rc->kScaleAndSmearMC(float(it->Charge()), it->Pt(), it->Eta(), it->Phi(), it->ActiveLayer(), u1, u2, 0,0);
    }
    it->SetRelIso(0.3,it->RelMiniAODIso03()/scalefactor);
    it->SetRelIso(0.4,it->RelMiniAODIso04()/scalefactor);
    it->SetPtEtaPhiM( (scalefactor*it->Pt() ), it->Eta(), it->Phi(), it->M());
  }  
  
  
  /*
    
    -------------------------------------------------------------------------------------
    Following variations are provided currently
    -------------------------------------------------------------------------------------
    set        members     comment
    Default  0          1           default, reference based on madgraph sample, with pt distribution reweighted to data in |Y| bins. 
    Stat     1          100         pre-generated stat. variations; can include more (e.g. 400). RMS over these would give stat. uncertainty
    Zpt      2          1           derived without reweighting reference to data. Difference wrt central can be taken as systematics
    Ewk      3          1           Ad-hoc weighting applied to reference to change sw2 and Zw. Use for x-check for now. Need to compare with Powheg (may become default). 
    CorDm    4          5           varied profile mass window; one can take maximum deviation as systematics
    FitDm    5          5           varied fitting mass window; one can take maximum deviation as systematics
    LHEw     6          0           (please ignore for now)
    Run      7          7           derived from B,C,D,E,F,G,H; for cross-check study. No dramatic run denendence observed. MC's should be split by lumi if applied this way. 
    AMCNLO   8          1           derived from limited stat. amcnlo sample; for cross-check study
    -------------------------------------------------------------------------------------
  */

}

float MCDataCorrections::GetCorrectedMuonMomentum(snu::KMuon muon, std::vector<snu::KTruth> truth){
  double scalefactor = 1.;
  if (corr_isdata) scalefactor = rc->kScaleDT(float(muon.Charge()), muon.Pt(), muon.Eta(), muon.Phi(),0,0);
  else {
    //gRandom->SetSeed(1111);
    double u1 = gRandom->Rndm();
    double u2 = gRandom->Rndm();

    int mu_index = muon.MCTruthIndex();
    float genpt(-999.);

    if(mu_index > 0&& mu_index < truth.size()) {
      if(fabs(truth.at(mu_index).PdgId() ) == 13) genpt = truth.at(mu_index).Pt();
    }
    
    if ( genpt> 0.)  scalefactor = rc->kScaleFromGenMC(float(muon.Charge()), muon.Pt(), muon.Eta(), muon.Phi(), muon.ActiveLayer(), genpt, u1,0, 0);
    else scalefactor = rc->kScaleAndSmearMC(float(muon.Charge()), muon.Pt(), muon.Eta(), muon.Phi(), muon.ActiveLayer(), u1, u2, 0,0);
  }
  
  return (scalefactor*muon.Pt());
}

vector<TLorentzVector> MCDataCorrections::MakeTLorentz(vector<snu::KElectron> el){

  vector<TLorentzVector> tl_el;
  for(vector<KElectron>::iterator itel=el.begin(); itel!=el.end(); ++itel) {
    TLorentzVector tmp_em;
    tmp_em.SetPtEtaPhiM((*itel).Pt(),(*itel).Eta(),(*itel).Phi(),(*itel).M());
    tl_el.push_back(tmp_em);
  }
  return tl_el;
}

vector<TLorentzVector> MCDataCorrections::MakeTLorentz(vector<snu::KMuon> mu){
  
  vector<TLorentzVector> tl_mu;
  for(vector<KMuon>::iterator itmu=mu.begin(); itmu!=mu.end(); ++itmu) {
    TLorentzVector tmp_mu;
    tmp_mu.SetPtEtaPhiM((*itmu).Pt(),(*itmu).Eta(),(*itmu).Phi(),(*itmu).M());
    tl_mu.push_back(tmp_mu);
  }
  return tl_mu;
}


vector<TLorentzVector> MCDataCorrections::MakeTLorentz(vector<snu::KJet> j){

  vector<TLorentzVector> tl_jet;
  for(vector<KJet>::iterator itj=j.begin(); itj!=j.end(); ++itj) {
    TLorentzVector tmp_j;
    tmp_j.SetPtEtaPhiM((*itj).Pt(),(*itj).Eta(),(*itj).Phi(),(*itj).M());
    tl_jet.push_back(tmp_j);
  }
  return tl_jet;
}


vector<TLorentzVector> MCDataCorrections::MakeTLorentz(vector<snu::KFatJet> j){

  vector<TLorentzVector> tl_jet;
  for(vector<KFatJet>::iterator itj=j.begin(); itj!=j.end(); ++itj) {
    TLorentzVector tmp_j;
    tmp_j.SetPtEtaPhiM((*itj).Pt(),(*itj).Eta(),(*itj).Phi(),(*itj).M());
    tl_jet.push_back(tmp_j);
  }
  return tl_jet;
}

TDirectory* MCDataCorrections::getTemporaryDirectory(void) const
{

  
  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:               
    std::stringstream dirname;
    dirname << "MCDataCorrections_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:                                        
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }

  return tempDir;
}



