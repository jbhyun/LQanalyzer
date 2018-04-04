/***************************************************************************
 * @Project: LQAnalyzer Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes 
#include "AnalyzerCore.h"
#include "EventBase.h"

//Plotting                                                      
#include "MuonPlots.h"
#include "ElectronPlots.h"
#include "JetPlots.h"
#include "SignalPlotsEE.h"
#include "SignalPlotsMM.h"
#include "SignalPlotsEM.h"
#include "TriLeptonPlots.h"
#include "HNpairPlotsMM.h"
#include "HNTriLeptonPlots.h"

//ROOT includes
#include <TFile.h>
#include "TStyle.h"



AnalyzerCore::AnalyzerCore() : LQCycleBase(), n_cutflowcuts(0), MCweight(-999.),reset_lumi_mask(false),changed_target_lumi(false), k_reset_period(false), a_mcperiod(-1),comp_file_firstev(true) {

  k_debugmode=false;
  IDSetup=false;  
  setupDDBkg=false;

  //// set rare collections to false. Can be set truth in analysis code
  k_usegenjet=true;
  k_usetruth=true;
  k_usephotons=true;
  k_usefatjet=true;

  fake_path="";
  fake_configured=true;
  self_configured=false;

  TH1::SetDefaultSumw2(true);  
  /// clear list of triggers stored in KTrigger
  triggerlist.clear();
  cutflow_list.clear();
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("AnalyzerCore");

  Message("In AnalyzerCore constructor", INFO);
  
  /////////////////////////////////////////////////////////////////////// 
  //////// For HN analysis  /////////////////////////////////////////////  
  //////////////////////////////////////////////////////////////////////  
  //// MC Data corrections
  //////////////////////////////////////////////////////////////////////                                                                                                   
  
  //mcdata_correction = new MCDataCorrections();
  
  string lqdir =  getenv("LQANALYZER_DIR");

  string username = getenv("USER");

  if( TString(getenv("CATDEBUG")) == "True") k_debugmode=true;

  
  cout << "                                                  " << endl;
  cout << "   ########    ###       ###   ###  ###           " << endl;
  cout << "   ########    ####      ###   ###  ###           " << endl;
  cout << "   ###         #####     ###   ###  ###           " << endl;
  cout << "   ###         ### ##    ###   ###  ###           " << endl;
  cout << "   ########    ###  ##   ###   ###  ###           " << endl;
  cout << "   ########    ###   ##  ###   ###  ###           " << endl;
  cout << "        ###    ###    ## ###   ###  ###           " << endl;
  cout << "        ###    ###     #####   ###  ###           " << endl;
  cout << "   ########    ###      ####   ########           " << endl;
  cout << "   ########    ###       ###   ########           " << endl;
  cout << "                                                  " << endl;

  cout << "##################################################" << endl;
  compmap.clear();
  compmap2.clear();
  
  if((TString(getenv("USER")) == "jskim" || TString(getenv("USER")) =="shjeon")){
    //==== HN Gen Matching Class
    m_HNgenmatch = new HNGenMatching();
  }

  cout <<  "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;


}

bool AnalyzerCore::CheckEventComparison(TString user, TString label, TString user2, TString label2, bool switchorder){
  
  if(!switchorder){
    if(compmap.size() ==0) compmap=CheckEventComparisonList(user, label, user2, label2);
  }
  else  if(compmap2.size() ==0) compmap2=CheckEventComparisonList(user2, label2, user, label);

  if(!switchorder){
    for(map<int,int>::iterator mit = compmap.begin(); mit != compmap.end(); mit++){
      if(mit->second == eventbase->GetEvent().RunNumber() && mit->first == eventbase->GetEvent().EventNumber() ) return true;
    }
  }
  else{
    for(map<int,int>::iterator mit = compmap2.begin(); mit != compmap2.end(); mit++){
      if(mit->second == eventbase->GetEvent().RunNumber() && mit->first == eventbase->GetEvent().EventNumber() ) return true;
    }

  }
  
  return false;
}

map<int, int> AnalyzerCore::CheckEventComparisonList(TString user, TString label, TString user2, TString label2){

    map<int, int> diffmap;
    map<int, int> list1;
    map<int, int> list2;

    cout << "CheckEventComparisonList " <<  "/data1/LQAnalyzer_rootfiles_for_analysis/EventComparisons/"+  user+"/" + label + ".txt" << endl;
    if(1){
      ifstream comp(( "/data1/LQAnalyzer_rootfiles_for_analysis/EventComparisons/"+  user+"/" + label + ".txt"));
      if(!comp) {
	cout << "file " << "/data1/LQAnalyzer_rootfiles_for_analysis/EventComparisons/"+  user+"/" + label + ".txt not found" << endl;
	exit(EXIT_FAILURE);
      }
      
      string lline;
      while(getline(comp,lline) ){
	std::istringstream is( lline );
	TString blank1;
	TString blank2;
	int run;
	TString tmp;
	int ev;
	float met;
	is >>blank1;
	is >>blank2;
	is >> run;
	is >> tmp;
	is >> ev;
	is >> met;
	is >> tmp;
	if(blank2!=user) break;
	list1[ev] =run;
	continue;
      }
    }
    if(1){
      ifstream comp(( "/data1/LQAnalyzer_rootfiles_for_analysis/EventComparisons/"+  user2+"/" + label2 + ".txt"));
      if(!comp) {
	cout << "file " << "/data1/LQAnalyzer_rootfiles_for_analysis/EventComparisons/"+  user2+"/" + label2 + ".txt not found" << endl;
	exit(EXIT_FAILURE);
      }
      
      string lline;
      while(getline(comp,lline) ){
	std::istringstream is( lline );
	TString blank1;
	TString blank2;
	int run;
	TString tmp;
	int ev;
	float met;
	is >>blank1;
	is >>blank2;
	is >> run;
	is >> tmp;
	is >> ev;
	is >> met;
	is >> tmp;
	if(blank2!=user) break;
	list2[ev] =run;
	continue;
      }
    }
    
    for(map<int,int>::iterator mit = list1.begin(); mit != list1.end(); mit++){
      bool found=false;
      for(map<int,int>::iterator mit2 = list2.begin(); mit2 != list2.end(); mit2++){
	if(mit2->first == mit->first && mit2->second==mit2->second) found=true;
      }
      if(!found) {

	diffmap[mit->first] = mit->second;
      }
    }


    return diffmap;
}

void AnalyzerCore::FillEventComparisonFile(TString label){
  
  //// Make TEX file                                                                                                                                                        
  ofstream ofile_tex;
  string lqdir = getenv("LQANALYZER_DIR");

  label = label + k_tag_name +"_" +k_sample_name;
  //label = label + k_sample_name;
  string compfile = "/data1/LQAnalyzer_rootfiles_for_analysis/EventComparisons/"+  string(getenv("USER")) + "/"+string(label)+ ".txt";     

  //if(comp_file_firstev)    ofile_tex.open(compfile.c_str());
  // else
  ofile_tex.open(compfile.c_str(),ios::out | ios::app);

  ofile_tex.setf(ios::fixed,ios::floatfield);
  ofile_tex << "[ "<<getenv("USER") << " " << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << " " <<  eventbase->GetEvent().MET()<<" ]"<< endl;
  ofile_tex.close();
  comp_file_firstev=false;
  
}

vector<TString >  AnalyzerCore::GetHNDiLepElTriggers(){

  vector<TString> triglist;
  if(isData){
    if(k_channel.Contains("SingleElectron"))triglist.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");
    if(k_channel.Contains("DoubleEG")) triglist.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  }
  else{
    triglist.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");
    triglist.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
  }
  return triglist;

}
bool AnalyzerCore::FailHNDataSetCheck(){

  bool _singleEG =(k_channel.Contains("SingleElectron"));
  bool _singleMuon =(k_channel.Contains("SingleMuon"));
  TString analysis_trigger_eg="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";
  TString analysis_trigger_muon="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v";
  if(isData && _singleEG && PassTrigger(analysis_trigger_eg)) return false;
  if(isData && _singleMuon && PassTrigger(analysis_trigger_muon)) return false;

  return true;
}

void AnalyzerCore::setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);

  tdrStyle->SetEndErrorSize(2);
  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.1);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

  // For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.07, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);

  // For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:
  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  
  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

  

}

float AnalyzerCore::GetConvWeight(snu::KMuon mu){

  if(mu.Pt() < 10) return  0.460224;
  else if(mu.Pt() < 15) return 0.476428;
  else if(mu.Pt() < 20) return 0.531144;
  else if(mu.Pt() < 25) return 0.57826;
  else if(mu.Pt() < 30) return 0.591419;
  else if(mu.Pt() < 35) return 0.64385;
  else if(mu.Pt() < 45) return 0.641256;
  else if(mu.Pt() < 60) return 0.724696;
  else if(mu.Pt() < 100) return 0.727273;
  else return 0.730769;

}

void AnalyzerCore::SetupLuminosityMap(bool initialsetup, TString forceperiod){
  if(isData) return ;

  TString lumitriggerpath="";
  TString singleperiod = getenv("CATAnalyzerPeriod");
  
  if(!initialsetup) {
    singleperiod=forceperiod;
    trigger_lumi_map_cat2016.clear();
  }
  string lqdir = getenv("LQANALYZER_DIR");
  
  if(singleperiod.Contains("None")){
    lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+".txt";
  }
  else{
    if(singleperiod== "B")   lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_272007_275376.txt";
    else if(singleperiod=="C")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_275657_276283.txt";
    else if(singleperiod=="D")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_276315_276811.txt";
    else if(singleperiod=="E")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_276831_277420.txt";
    else if(singleperiod=="F")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_277772_278808.txt";
    else if(singleperiod=="G")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_280919_284044.txt";
    else if(singleperiod.Contains("H"))lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_278820_280385.txt";
    else if(singleperiod=="GH")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_280919_280385.txt";
    else {  cerr << "Wrong period setting in SetupLuminosityMap"<< endl;  exit(EXIT_FAILURE);}
  }
  ifstream triglumi2016(lumitriggerpath.Data());
  if(!triglumi2016) {
    cerr << "Did not find period  "+ lumitriggerpath + " exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  
  string sline2016;
  

  if(k_debugmode)cout << "Trigname : Lumi pb-1" << endl;
  while(getline(triglumi2016,sline2016) ){
    std::istringstream is( sline2016 );
    
    string trigname;
    float trig_lumi;
    is >> trigname;
    if(trigname=="###" ) continue;
    
    
    if(trigname=="END") break;
    if(!TString(trigname).Contains("Lumi")){
      is >> trig_lumi;
      if(k_debugmode)cout << trigname << " " << trig_lumi << endl;
    }
    else{
      string tmp;
      is >> tmp;
      is >> trig_lumi;
      if(k_debugmode)cout << trigname << " " << trig_lumi << endl;
    }
    trigger_lumi_map_cat2016[TString(trigname)] = trig_lumi;
    continue;
  }
  triglumi2016.close();
}


bool  AnalyzerCore::Check(float val){

  if(abs(val) == 999.) return false;

  return true;
}

float AnalyzerCore::MC_CR_Correction(int syst){
  
  float fsyst = 0.;
  if(syst==1) fsyst=1.;
  if(syst==-1) fsyst=-1.;

  ///  updated 2 Oct

  if(k_sample_name.Contains("WZTo3LNu_powheg")) return 0.988021 +  fsyst*0.0652167;
  if(k_sample_name.Contains("ZGto2LG")) return  0.83069 + fsyst*0.174719;
  if(k_sample_name.Contains("WGtoLNuG")) return 1.;
  if(k_sample_name.Contains("ZZTo4L_powheg")) return 0.94117 + fsyst*0.105665;
  if(k_sample_name.Contains("ggZZto2e2mu")) return 0.94117 + fsyst*0.105665;
  if(k_sample_name.Contains("ggZZto2e2nu")) return 0.94117 + fsyst*0.105665;
  if(k_sample_name.Contains("ggZZto2e2tau")) return 0.94117 + fsyst*0.105665;
  if(k_sample_name.Contains("ggZZto2mu2nu")) return 0.94117 + fsyst*0.105665;
  if(k_sample_name.Contains("ggZZto2mu2tau")) return 0.94117 + fsyst*0.105665;
  if(k_sample_name.Contains("ggZZto4e")) return 0.94117 + fsyst*0.105665;
  if(k_sample_name.Contains("ggZZto4mu")) return 0.94117 + fsyst*0.105665;
  if(k_sample_name.Contains("ggZZto4tau")) return 0.94117 + fsyst*0.105665;
  
  return 1.;
}

float AnalyzerCore::GetTriggerPrescaleCorrection(TString triggername){
  float corr_trig=1.;

  if(triggername.Contains( "HLT_Mu3_PFJet40_v")) corr_trig = 0.728;
  if(triggername.Contains("HLT_Mu8_TrkIsoVVL_v")) corr_trig = 1.399;
  return corr_trig;

}


float AnalyzerCore::GetKFactor(){
  
  if(k_sample_name.Contains("ZZTo4L_powheg")) {
    // Physics Letters B 750 (2015) 407–410 --New
    return 1.15;
  }
  if(k_sample_name.Contains("WZTo3LNu_powheg")){
    //Physics Letters B 761 (2016) 179–183
    //http://dx.doi.org/10.1016/j.physletb.2016.08.017 
    return 1.109;
  }
  if(k_sample_name.Contains("ttZ") && !k_sample_name.Contains("To")){
    return 839.3/780.;
  }
  if(k_sample_name.Contains("ttW") && !k_sample_name.Contains("To")){
    return 600.8/610.;
  }

  if(k_sample_name.Contains("ZZTo2L2Nu_Powheg")) {
    // http://arxiv.org/abs/1405.2219                                                                                     
    //  1.16[3] brings pp->ZZ from NLO to NNLO                                                                            
    return 1.16;
  }
  if(k_sample_name.Contains("ZZTo2L2Q_Powheg")) {
    // http://arxiv.org/abs/1405.2219                                                                                        //  1.16[3] brings pp->ZZ from NLO to NNLO 
    
    return 1.16;
  }

  if(k_sample_name.Contains("ggZZto")){
    //  1.67 brings gg->ZZ from LO to NLO (http://arxiv.org/abs/1509.06734)
    return 1.67;
  }
  if(k_sample_name.Contains("ggHtoZZ")){
    return 1.67;
    //AN2016_359
  }

  return 1.;
    
}

void  AnalyzerCore::CorrectedMETRochester( std::vector<snu::KMuon> muall){

  /// function returns corrected met + can be used to set event met to corrected met

  float met_x =eventbase->GetEvent().PFMETx(); 
  float met_y =eventbase->GetEvent().PFMETy();
  
  float px_orig(0.), py_orig(0.),px_corrected(0.), py_corrected(0.);
  for(unsigned int im=0; im < muall.size() ; im++){
      
      px_orig+= muall.at(im).MiniAODPt()*TMath::Cos(muall.at(im).Phi());
      py_orig+= muall.at(im).MiniAODPt()*TMath::Sin(muall.at(im).Phi());
      px_corrected += muall.at(im).Px();
      py_corrected += muall.at(im).Py();
      
  }
  
  if(!eventbase->GetEvent().PropagatedRochesterToMET()){
    met_x = met_x + px_orig - px_corrected;	
    met_y = met_y + py_orig - py_corrected;	
  }
  
  if(!eventbase->GetEvent().PropagatedRochesterToMET()){
    snu::KEvent tempev = eventbase->GetEvent();
    tempev.SetMET(snu::KEvent::pfmet,  sqrt(met_x*met_x + met_y*met_y), TMath::ATan2(met_y,met_x), eventbase->GetEvent().SumET());
    tempev.SetPFMETx(met_x);
    tempev.SetPFMETy(met_y);
    tempev.SetPropagatedRochesterToMET(true);
    eventbase->SetEventBase(tempev);
  }

  return;
}   
void  AnalyzerCore::CorrectedMETJMR( std::vector<snu::KFatJet>  fjetall, std::vector<snu::KJet>  jetall){

  /// function returns corrected met + can be used to set event met to corrected met                                                                                                                                          

  float met_x =eventbase->GetEvent().PFMETx();
  float met_y =eventbase->GetEvent().PFMETy();

  float px_orig(0.), py_orig(0.),px_corrected(0.), py_corrected(0.);
  float px_orig_ak4(0.), py_orig_ak4(0.),px_corrected_ak4(0.), py_corrected_ak4(0.);
  for(unsigned int ij=0; ij < fjetall.size() ; ij++){

    px_orig+=  fjetall.at(ij).MiniAODPt()*TMath::Cos( fjetall.at(ij).Phi());
    py_orig+=  fjetall.at(ij).MiniAODPt()*TMath::Sin( fjetall.at(ij).Phi());
    px_corrected += fjetall.at(ij).Px();
    py_corrected += fjetall.at(ij).Py();
    
  }
  
  for(unsigned int ij=0; ij < jetall.size() ; ij++){
    for(unsigned int fij=0; fij < fjetall.size() ; fij++){
      if (jetall[ij].DeltaR(fjetall[fij]) < 0.8){
	px_orig_ak4+=  (jetall.at(ij).Pt()/jetall.at(ij).SmearedRes())*TMath::Cos( jetall.at(ij).Phi());
	py_orig_ak4+=  (jetall.at(ij).Pt()/jetall.at(ij).SmearedRes())*TMath::Sin( jetall.at(ij).Phi());
	px_corrected_ak4 += jetall.at(ij).Px();
	py_corrected_ak4 += jetall.at(ij).Py();
      }
    }
  }

  if(!eventbase->GetEvent().PropagatedJMRToMET()){
    met_x = met_x + px_orig - px_corrected - px_orig_ak4 + px_corrected_ak4;
    met_y = met_y + py_orig - py_corrected - py_orig_ak4 + py_corrected_ak4;
  }
  
  
  if(!eventbase->GetEvent().PropagatedRochesterToMET()){
    snu::KEvent tempev = eventbase->GetEvent();
    tempev.SetMET(snu::KEvent::pfmet,  sqrt(met_x*met_x + met_y*met_y), TMath::ATan2(met_y,met_x), eventbase->GetEvent().SumET());
    tempev.SetPFMETx(met_x);
    tempev.SetPFMETy(met_y);
    tempev.SetPropagatedJMRToMET(true);
    eventbase->SetEventBase(tempev);
  }

  return;
}





void  AnalyzerCore::CorrectedMETElectron(int sys, std::vector<snu::KElectron> elall,  double& OrignialMET, double& OriginalMETPhi){

  if(sys==0) return;

  float met_x =eventbase->GetEvent().PFMETx();
  float met_y =eventbase->GetEvent().PFMETy();

  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int iel=0; iel < elall.size() ; iel++){


    px_orig+= elall.at(iel).Px();
    py_orig+= elall.at(iel).Py();
    if(sys==1){
      px_shifted += elall.at(iel).Px()*elall.at(iel).PtShiftedUp();
      py_shifted += elall.at(iel).Py()*elall.at(iel).PtShiftedUp();
    }
    if(sys==-1){
      px_shifted += elall.at(iel).Px()*elall.at(iel).PtShiftedDown();
      py_shifted += elall.at(iel).Py()*elall.at(iel).PtShiftedDown();
    }


  }
  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;
  OrignialMET =  sqrt(met_x*met_x + met_y*met_y);
  OriginalMETPhi = TMath::ATan2(met_y,met_x);
  

}

void  AnalyzerCore::CorrectedMETMuon( int sys, std::vector<snu::KMuon> muall,   double& OrignialMET, double& OriginalMETPhi){
  
  if(sys==0) return;
  
  float met_x1 = OrignialMET*TMath::Cos(OriginalMETPhi);
  float met_y1 = OrignialMET*TMath::Sin(OriginalMETPhi);

  cout << "MET " << OrignialMET << " " << eventbase->GetEvent().PFMET() << endl;

  float met_x =eventbase->GetEvent().PFMETx();
  float met_y =eventbase->GetEvent().PFMETy();
  
  cout << met_x1 << " " << met_x << endl;

  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int imu=0; imu < muall.size() ; imu++){
    
    px_orig+= muall.at(imu).Px();
    py_orig+= muall.at(imu).Py();
    if(sys==1){
      px_shifted += muall.at(imu).Px()*muall.at(imu).PtShiftedUp();
      py_shifted += muall.at(imu).Py()*muall.at(imu).PtShiftedUp();
    }
    if(sys==-1){
      px_shifted += muall.at(imu).Px()*muall.at(imu).PtShiftedDown();
      py_shifted += muall.at(imu).Py()*muall.at(imu).PtShiftedDown();
    }  
  }
  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;
  
  OrignialMET =  sqrt(met_x*met_x + met_y*met_y);
  OriginalMETPhi = TMath::ATan2(met_y,met_x);
}



void  AnalyzerCore::CorrectedMETJES(int sys, vector<snu::KJet> jetall, vector<snu::KFatJet> fjetall,  double& OrignialMET, double& OriginalMETPhi){

  if(sys==0) return;

  float met_x =eventbase->GetEvent().PFMETx();
  float met_y =eventbase->GetEvent().PFMETy();

  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int ij=0; ij < jetall.size() ; ij++){


    px_orig+= jetall.at(ij).Px();
    py_orig+= jetall.at(ij).Py();
    if(sys==1){

      px_shifted += jetall.at(ij).Px()*jetall.at(ij).ScaledUpEnergy();
      py_shifted += jetall.at(ij).Py()*jetall.at(ij).ScaledUpEnergy();

    }
    if(sys==-1){
      px_shifted += jetall.at(ij).Px()*jetall.at(ij).ScaledDownEnergy();
      py_shifted += jetall.at(ij).Py()*jetall.at(ij).ScaledDownEnergy();

    }

  }
  for(unsigned int ij=0; ij < fjetall.size() ; ij++){


    px_orig+= fjetall.at(ij).Px();
    py_orig+= fjetall.at(ij).Py();
    if(sys==1){

      px_shifted += fjetall.at(ij).Px()*fjetall.at(ij).ScaledUpEnergy();
      py_shifted += fjetall.at(ij).Py()*fjetall.at(ij).ScaledUpEnergy();

    }
    if(sys==-1){
      px_shifted += jetall.at(ij).Px()*fjetall.at(ij).ScaledDownEnergy();
      py_shifted += jetall.at(ij).Py()*fjetall.at(ij).ScaledDownEnergy();

    }

  }

  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;

  OrignialMET =  sqrt(met_x*met_x + met_y*met_y);
  OriginalMETPhi = TMath::ATan2(met_y,met_x);


}



void  AnalyzerCore::CorrectedMETJMS(int sys, vector<snu::KFatJet> fjetall,  double& OrignialMET, double& OriginalMETPhi){

  if(sys==0) return;

  float met_x =eventbase->GetEvent().PFMETx();
  float met_y =eventbase->GetEvent().PFMETy();

  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int ij=0; ij < fjetall.size() ; ij++){


    px_orig+= fjetall.at(ij).Px();
    py_orig+= fjetall.at(ij).Py();
    if(sys==1){

      px_shifted += fjetall.at(ij).Px()*fjetall.at(ij).ScaledMassUp();
      py_shifted += fjetall.at(ij).Py()*fjetall.at(ij).ScaledMassUp();

    }
    if(sys==-1){
      px_shifted += fjetall.at(ij).Px()*fjetall.at(ij).ScaledMassDown();
      py_shifted += fjetall.at(ij).Py()*fjetall.at(ij).ScaledMassDown();

    }

  }

  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;

  OrignialMET =  sqrt(met_x*met_x + met_y*met_y);
  OriginalMETPhi = TMath::ATan2(met_y,met_x);


}



void AnalyzerCore::CorrectedMETJER(int sys, vector<snu::KJet> jetall, vector<snu::KFatJet> fjetall,   double& OrignialMET, double& OriginalMETPhi){


  float met_x =eventbase->GetEvent().PFMETx();
  float met_y =eventbase->GetEvent().PFMETy();

  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int ij=0; ij < jetall.size() ; ij++){


    px_orig+= jetall.at(ij).Px();
    py_orig+= jetall.at(ij).Py();
    if(sys==1){
      px_shifted += jetall.at(ij).Px()*(jetall.at(ij).SmearedResUp() / jetall.at(ij).SmearedRes()  );
      py_shifted += jetall.at(ij).Py()*(jetall.at(ij).SmearedResUp() / jetall.at(ij).SmearedRes() );
      
    }
    if(sys==-1){
      px_shifted += jetall.at(ij).Px()*(jetall.at(ij).SmearedResDown()/jetall.at(ij).SmearedRes() );
      py_shifted += jetall.at(ij).Py()*(jetall.at(ij).SmearedResDown()/jetall.at(ij).SmearedRes()) ;
      
    }
    
  }
  for(unsigned int ij=0; ij < fjetall.size() ; ij++){


    px_orig+= fjetall.at(ij).Px();
    py_orig+= fjetall.at(ij).Py();
    if(sys==1){
      px_shifted += fjetall.at(ij).Px()*(fjetall.at(ij).SmearedResUp()/ fjetall.at(ij).SmearedRes());
      py_shifted += fjetall.at(ij).Py()*(fjetall.at(ij).SmearedResUp()/fjetall.at(ij).SmearedRes());
    }
    if(sys==-1){
      px_shifted += fjetall.at(ij).Px()*(fjetall.at(ij).SmearedResDown()/fjetall.at(ij).SmearedRes());
      py_shifted += fjetall.at(ij).Py()*(fjetall.at(ij).SmearedResDown()/fjetall.at(ij).SmearedRes());
    }
  }


  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;


  OrignialMET =  sqrt(met_x*met_x + met_y*met_y);
  OriginalMETPhi = TMath::ATan2(met_y,met_x);


}


void AnalyzerCore::CorrectedMETJMR(int sys, vector<snu::KFatJet> fjetall,   double& OrignialMET, double& OriginalMETPhi){


  float met_x =eventbase->GetEvent().PFMETx();
  float met_y =eventbase->GetEvent().PFMETy();

  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int ij=0; ij < fjetall.size() ; ij++){


    px_orig+= fjetall.at(ij).Px();
    py_orig+= fjetall.at(ij).Py();
    if(sys==1){
      px_shifted += fjetall.at(ij).Px()*(fjetall.at(ij).SmearedMassResUp());
      py_shifted += fjetall.at(ij).Py()*(fjetall.at(ij).SmearedMassResUp());
    }
    if(sys==-1){
      px_shifted += fjetall.at(ij).Px()*(fjetall.at(ij).SmearedMassResDown());
      py_shifted += fjetall.at(ij).Py()*(fjetall.at(ij).SmearedMassResDown());
    }
  }
  
  
  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;


  OrignialMET =  sqrt(met_x*met_x + met_y*met_y);
  OriginalMETPhi = TMath::ATan2(met_y,met_x);


}




float AnalyzerCore::GetFatJetSF(snu::KFatJet fjet, float tau21cut, int sys){
  
  float fsys = -1;
  if(sys > 0) fsys =1;
  if(sys==0) fsys=0.;
  if(tau21cut == 0.45){
    if((fjet.Tau2()/fjet.Tau1())  < 0.45)     return 0.88 + fsys*0.1;
    else return 1.;
  }
  if(tau21cut == 0.6){
    if((fjet.Tau2()/fjet.Tau1()) < 0.6)     return 1.11 + fsys*0.08;
    else return 1.;
  }
  else return 1.;
  
}


vector<snu::KFatJet>  AnalyzerCore::GetCorrectedFatJet(vector<snu::KFatJet>   fjets){

  vector<snu::KFatJet>  corr_fatjets;
  
  for(unsigned int ifj=0; ifj < fjets.size(); ifj++){
    snu::KFatJet fjet = fjets[ifj];
    float L1corr = fjet.L1JetCorr();
    
    TLorentzVector v;
    v.SetPtEtaPhiM(fjet.Pt(), fjet.Eta(), fjet.Phi(), fjet.M());
    
    if(L1corr==0.) L1corr = 0.95;
    /// remove L1 correction (only L2L3 used)
    v=v* (1./L1corr);
    
    /// smear mass with JMR central
    v=v*fjet.SmearedRes();
    fjet.SetPrunedMass(fjet.PrunedMass() * fjet.SmearedMassRes()/L1corr);
    snu::KFatJet fjet_corr(fjet);
    if(fjet_corr.MiniAODPt() <0)fjet_corr.SetMiniAODPt(fjet_corr.Pt());
    fjet_corr.SetPtEtaPhiM(v.Pt(), v.Eta(), v.Phi(), v.M());
    
    corr_fatjets.push_back(fjet_corr);
  }

  return corr_fatjets;
}

snu::KJet AnalyzerCore::GetCorrectedJetCloseToLepton(snu::KElectron el, snu::KJet jet, bool usem){
  //jet_LepAwareJECv2 = (raw_jet * L1 - lepton) * L2L3Res + lepton
  
  float rawpt= jet.RawPt();  
  float rawe= jet.RawE();
  float L1corr = jet.L1JetCorr();
  float l2l3res = jet.L2L3ResJetCorr();
  float leppt = el.Pt();
  float lepe = el.E();
  float corr_pt = (rawpt*L1corr - leppt)*l2l3res + leppt;
  float corr_e = (rawe*L1corr - lepe)*l2l3res + lepe;
  
  snu::KJet jet_corr(jet);
  
  if(usem){
    TLorentzVector v;
    v.SetPtEtaPhiM(jet.Pt(), jet.Eta(), jet.Phi(), jet.M());
    v=v*(corr_pt/jet.Pt());
    jet_corr.SetPtEtaPhiM(v.Pt(), v.Eta(), v.Phi(), v.M());
    return jet_corr;
  }
  else jet_corr.SetPtEtaPhiE(corr_pt, jet.Eta(), jet.Phi(),corr_e);

  return jet_corr;
}

snu::KJet AnalyzerCore::GetCorrectedJetCloseToLepton(snu::KMuon mu, snu::KJet jet){

  //jet_LepAwareJECv2 = (raw_jet * L1 - lepton) * L2L3Res + lepton                                                                                                          
  
  float rawpt= jet.RawPt();
  float rawe= jet.RawE();
  float L1corr = jet.L1JetCorr();
  float l2l3res = jet.L2L3ResJetCorr();
  float leppt = mu.Pt();
  float lepe = mu.E();
  float corr_pt = (rawpt*L1corr - leppt)*l2l3res + leppt;
  float corr_e = (rawe*L1corr - lepe)*l2l3res + lepe;

  snu::KJet jet_corr(jet);
  TLorentzVector v;
  v.SetPtEtaPhiM(jet.Pt(), jet.Eta(), jet.Phi(), jet.M());
  v=v*(corr_pt/jet.Pt());
  jet_corr.SetPtEtaPhiM(v.Pt(), v.Eta(), v.Phi(), v.M());
  return jet_corr;
}


float AnalyzerCore::GetPtRelLepTJet(snu::KElectron electron, std::vector<snu::KJet> jets, bool usecorrectedpt){

  if(jets.size() == 0) return -999.;
  if(electron.Pt() < 10.) return -999.;
  snu::KParticle closejet;
  float mindR=0.7;

  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( electron.DeltaR(jets.at(ijet)) < mindR){
      if(usecorrectedpt)closejet= GetCorrectedJetCloseToLepton(electron,jets.at(ijet));
      else closejet=jets.at(ijet);
      mindR=electron.DeltaR(jets.at(ijet));
    }
  }

  if(mindR==0.7) return 0.;

  FillHist(("ptrel_dr"),mindR, weight, 0., 4., 100);

  TVector3 el3=  electron.Vect();
  TVector3 jet3= closejet.Vect();
  TVector3 lepjetrel = jet3-el3;
  FillHist(("ptrel_lepjetmag"),lepjetrel.Mag(), weight, 0., 100., 100);
  FillHist(("ptrel_crosslepjetmag"), (lepjetrel.Cross(el3)).Mag(), weight, 0., 100., 100);
  float ptrel = (lepjetrel.Cross(el3)).Mag()/ lepjetrel.Mag();

  return ptrel;
}


float AnalyzerCore::GetPtRelLepTJet(snu::KMuon muon, std::vector<snu::KJet> jets,bool usecorrectedpt){

  if(jets.size() == 0) return -999.;
  if(muon.Pt() < 10.) return -999.;
  snu::KParticle closejet;
  float mindR=0.7;

  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( muon.DeltaR(jets.at(ijet)) < mindR){
      if(usecorrectedpt)closejet= GetCorrectedJetCloseToLepton(muon,jets.at(ijet));
      else closejet = jets.at(ijet);
      mindR=muon.DeltaR(jets.at(ijet));
    }
  }
  
  if(mindR==0.7) return 0.;
  
  FillHist(("ptrel_dr"),mindR, weight, 0., 4., 100);
  
  TVector3 el3=  muon.Vect();
  TVector3 jet3= closejet.Vect();
  TVector3 lepjetrel = jet3-el3;
  FillHist(("ptrel_lepjetmag"),lepjetrel.Mag(), weight, 0., 100., 100);
  FillHist(("ptrel_crosslepjetmag"), (lepjetrel.Cross(el3)).Mag(), weight, 0., 100., 100);
  float ptrel = (lepjetrel.Cross(el3)).Mag()/ lepjetrel.Mag();
  
  return ptrel;
}


float AnalyzerCore::GetJetsCloseToLeptonPt(snu::KElectron electron, std::vector<snu::KJet> jets,bool usecorrectedpt){

  float mindR=.4;
  float jetpT=-999.;

  if(electron.Pt() < 10.) return 0.;

  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( electron.DeltaR(jets.at(ijet)) < mindR){
      mindR=electron.DeltaR(jets.at(ijet));
      if(usecorrectedpt)      jetpT=GetCorrectedJetCloseToLepton(electron,jets.at(ijet)).Pt();
      else jetpT=jets.at(ijet).Pt();
    }
  }

  return jetpT;
}



float AnalyzerCore::GetJetsCloseToLeptonPt(snu::KMuon muon, std::vector<snu::KJet> jets,bool usecorrectedpt){
  float mindR=.4;
  float jetpT=-999.;

  if(muon.Pt() < 10.) return 0.;

  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( muon.DeltaR(jets.at(ijet)) < mindR){
      mindR=muon.DeltaR(jets.at(ijet));
      if(usecorrectedpt)jetpT=GetCorrectedJetCloseToLepton(muon,jets.at(ijet)).Pt();
      else jetpT=jets.at(ijet).Pt();
    }
  }
  return jetpT;
}





float AnalyzerCore::MassDrop(snu::KElectron electron, std::vector<snu::KJet> jets,bool usecorrectedpt){
  if(jets.size() == 0) return -999.;
  snu::KParticle closejet;
  float mindR=0.7;
  if(electron.Pt() < 10.) return -999.;

  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( electron.DeltaR(jets.at(ijet)) < mindR){
      if(usecorrectedpt)closejet= GetCorrectedJetCloseToLepton(electron,jets.at(ijet));
      else  closejet= jets.at(ijet);
      
      mindR=mindR;
    }
  }

  if(mindR >= 0.7)  return -999.;

  snu::KParticle lj = closejet+electron;

  return (lj.M() - closejet.M());


}

float AnalyzerCore::MassDrop(snu::KMuon muon, std::vector<snu::KJet> jets,bool usecorrectedpt){
  if(jets.size() == 0) return -999.;
  snu::KParticle closejet;
  float mindR=.7;

  if(muon.Pt() < 10.) return -999.;
  for(unsigned int ijet=0; ijet < jets.size(); ijet++){
    if( muon.DeltaR(jets.at(ijet)) < mindR){
      if(usecorrectedpt)closejet= GetCorrectedJetCloseToLepton(muon,jets.at(ijet));
      else  closejet= jets.at(ijet);
      mindR=mindR;
    }
  }
  if(mindR >= 0.7)  return -999.;

  snu::KParticle lj = closejet+muon;

  return (lj.M() - closejet.M());


}



// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
///// FUNCTION USED TO CREATE BTAG EFFICIENCIES USED BY BTAGSF.cxx CLass
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void AnalyzerCore::MakeBTagEfficiencyPlots(){
  GetJetTaggerEfficiences("CSVv2M",snu::KJet::CSVv2, snu::KJet::Medium); 
  GetJetTaggerEfficiences("CSVv2T",snu::KJet::CSVv2, snu::KJet::Tight); 
  GetJetTaggerEfficiences("cMVAv2M",snu::KJet::cMVAv2, snu::KJet::Medium); 
  GetJetTaggerEfficiences("cMVAv2L",snu::KJet::cMVAv2, snu::KJet::Loose); 
  GetJetTaggerEfficiences("cMVAv2T",snu::KJet::cMVAv2, snu::KJet::Tight);       
}
void AnalyzerCore::GetJetTaggerEfficiences(TString taggerWP, KJet::Tagger tag,  KJet::WORKING_POINT wp){
  // taken frmo https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#Example_code_AN5  (USE HADRON FLAVOUR)

  // Make tagger efficiency file for btagging 
  Float_t ptbins[7] = { 20., 40., 60., 80., 100., 120., 3000.};
  Float_t etabins[5] = { 0., 0.6, 1.2, 1.8, 2.4};


  for(unsigned int ij =0; ij < GetJets("JET_NOCUT").size(); ij++){
    
    snu::KJet j = GetJets("JET_NOCUT").at(ij);
    int hadronFlavor = j.HadronFlavour();
    if( abs(hadronFlavor)==5 ){
      FillHist("h2_BTaggingEff_"+taggerWP+"_Denom_b" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
      if( j.IsBTagged( tag, wp)) FillHist("h2_BTaggingEff_"+taggerWP+"_Num_b" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
    }
    else if( abs(hadronFlavor)==4 ){
      FillHist("h2_BTaggingEff_"+taggerWP+"_Denom_c" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
      if( j.IsBTagged( tag, wp))          FillHist("h2_BTaggingEff_"+taggerWP+"_Num_c" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
    }
    else {
      FillHist("h2_BTaggingEff_"+taggerWP+"_Denom_udsg" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
      if( j.IsBTagged( tag, wp) )         FillHist("h2_BTaggingEff_"+taggerWP+"_Num_udsg" , j.Pt(),  j.Eta(), weight, ptbins, 6 , etabins, 4);
    }
  }
}


double AnalyzerCore::GetIsoCut(snu::KElectron el, TString cutlable){

  if(fabs(el.SCEta())<1.479){
    if(cutlable.Contains("b050")) return 0.05;
    if(cutlable.Contains("b0525")) return 0.0525;
    if(cutlable.Contains("b055")) return 0.055;
    if(cutlable.Contains("b060")) return 0.06;
    if(cutlable.Contains("b065")) return 0.065;
    if(cutlable.Contains("b075")) return 0.075;
    if(cutlable.Contains("b100")) return 0.1;
    if(cutlable.Contains("b125")) return 0.125;
  }
  else { 
    if(cutlable.Contains("e050")) return 0.05;
    if(cutlable.Contains("e0525")) return 0.0525;
    if(cutlable.Contains("e055")) return 0.055;
    if(cutlable.Contains("e060")) return 0.06;
    if(cutlable.Contains("e065")) return 0.065;
    if(cutlable.Contains("e075")) return 0.075;
    if(cutlable.Contains("e100")) return 0.1;
    if(cutlable.Contains("e125")) return 0.125;
  }   
  
  return 0.;
}
  
double AnalyzerCore::GetDXYCut(snu::KElectron el, TString cutlable){

  if(fabs(el.SCEta())<1.479){
    if(cutlable.Contains("b050")) return 0.05;
    if(cutlable.Contains("b040")) return 0.04;
    if(cutlable.Contains("b030")) return 0.03;
    if(cutlable.Contains("b025")) return 0.025;
    if(cutlable.Contains("b020")) return 0.02;
    if(cutlable.Contains("b019")) return 0.019;
    if(cutlable.Contains("b018")) return 0.018;      
    if(cutlable.Contains("b017")) return 0.017;
    if(cutlable.Contains("b016")) return 0.016;
    if(cutlable.Contains("b015")) return 0.015;
    if(cutlable.Contains("b014")) return 0.014;
    if(cutlable.Contains("b013")) return 0.013;
    if(cutlable.Contains("b012")) return 0.012;
    if(cutlable.Contains("b011")) return 0.011;
    if(cutlable.Contains("b010")) return 0.01;
  }
  else{
    if(cutlable.Contains("e100")) return 0.05;
    if(cutlable.Contains("e050")) return 0.05;
    if(cutlable.Contains("e040")) return 0.04;
    if(cutlable.Contains("e035")) return 0.035;
    if(cutlable.Contains("e025")) return 0.025;
    if(cutlable.Contains("e020")) return 0.02;
    if(cutlable.Contains("e019")) return 0.019;
    if(cutlable.Contains("e018")) return 0.018;
    if(cutlable.Contains("e017")) return 0.017;
    if(cutlable.Contains("e016")) return 0.016;
    if(cutlable.Contains("e015")) return 0.015;
    if(cutlable.Contains("e014")) return 0.014;
    if(cutlable.Contains("e013")) return 0.013;
    if(cutlable.Contains("e012")) return 0.012;
    if(cutlable.Contains("e011")) return 0.011;
    if(cutlable.Contains("e010")) return 0.01;
    
  }
  return 0.;
}

int AnalyzerCore::GetDataPeriod(){
  
  /// returns 1 for peiord B.... 7 for period H
  if(eventbase->GetEvent().RunNumber() < 272007) return -1;
  else  if(eventbase->GetEvent().RunNumber() <= 275376) return 1; 
  else  if(eventbase->GetEvent().RunNumber() <= 276283) return 2; 
  else  if(eventbase->GetEvent().RunNumber() <= 276811) return 3; 
  else  if(eventbase->GetEvent().RunNumber() <= 277420) return 4; 
  else  if(eventbase->GetEvent().RunNumber() <= 278808) return 5; 
  else  if(eventbase->GetEvent().RunNumber() <= 280385) return 6; 
  else  if(eventbase->GetEvent().RunNumber() <= 284044) return 7; 
  else return -1;
}

int AnalyzerCore::GetPeriod(){

  if(isData) return GetDataPeriod();
  else return GetMCPeriod();

}

int AnalyzerCore::GetMCPeriod(){
  /// This function returns a period B-H for MC events. 
  /// It uses a random number and retrunds a period based on the luminosity of each period
  /// It assumes the trigger used is unprescaled
  if(isData) return -1;
  if(!k_reset_period) return a_mcperiod;
  k_reset_period=false;
  


  a_mcperiod = k_mcperiod;
  return a_mcperiod;
  
}

int AnalyzerCore::GetMCPeriodRandom(){
  if(isData) return -1;
  double r =gRandom->Rndm(); /// random number between 0 and 1
  
  // https://docs.google.com/spreadsheets/d/1rWM3AlFKO8IJVaeoQkWZYWwSvicQ1QCXYSzH74QyZqE/edit?alt=json#gid=1689385956
  // using single muon luminosities (luminosities differ slightky for each dataset but difference is neglibable)
  double lumi_periodB = 5.929001722;
  double lumi_periodC = 2.645968083;
  double lumi_periodD = 4.35344881;
  double lumi_periodE = 4.049732039;
  double lumi_periodF = 3.157020934;
  double lumi_periodG = 7.549615806;
  double lumi_periodH = 8.545039549 + 0.216782873;
  double total_lumi = (lumi_periodB+lumi_periodC + lumi_periodD + lumi_periodE + lumi_periodF + lumi_periodG + lumi_periodH) ;
  
  vector<double> cum_lumi;
  cum_lumi.push_back(lumi_periodB/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC)/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC+lumi_periodD)/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE)/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF)/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG)/total_lumi); 
  cum_lumi.push_back((lumi_periodB+lumi_periodC+lumi_periodD+lumi_periodE+lumi_periodF+lumi_periodG+lumi_periodH)/total_lumi); 
  
  /// returns an int

  /// r = 1       |  period B
  /// r = 2       |  period C
  /// r = 3       |  period D
  /// r = 4       |  period E
  /// r = 5       |  period F
  /// r = 6       |  period G
  /// r = 7       |  period H

  for(unsigned int i=0; i < cum_lumi.size(); i++){
    if ( r < cum_lumi.at(i)) {
      a_mcperiod =  (i+1);
      return a_mcperiod;
    }
  }

  /// return period H is for some reason r > cum_lumi.at(max) 'should not happen'
  return  cum_lumi.size();

}


void AnalyzerCore::SetupSelectionJet(std::string path_sel){
  Message("SetupSelectionJet", DEBUG);
  ifstream jetselconfig(path_sel.c_str());
  if(!jetselconfig) {
    cerr << "Did not find "+ path_sel+", exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  else{
    cout << "Found " + path_sel << endl;
   
  }
  string jetline;
  int ncuts=0;
  vector<TString> cutnames;
  while(getline(jetselconfig,jetline) ){
    vector<pair<TString,TString> > string_jetsel;
    vector<pair<TString,float> > float_jetsel;
    TString idlabel;
    if (TString(jetline).Contains("webpage")) continue;
    if (TString(jetline).Contains("###")) continue;
    if (TString(jetline).Contains("ncut")) {
      std::istringstream is( jetline );
      string tmp;
      int itmp;
      is >> tmp;
      is >> itmp;
      ncuts = 2*(itmp +1);
      continue;
    }
    if (TString(jetline).Contains("ptmin")) {
      std::istringstream is( jetline );
      string tmp;
      for (int x =0; x < ncuts; x++){
        is >> tmp;
        cutnames.push_back(TString(tmp));
      }
    }
    else{
      std::istringstream is( jetline );
      string tmp;
      float tmpf;
      for (int x =0; x < ncuts; x++){
        if ( x%2 ==0) {
          is >> tmp;
          continue;
        }

        if (x > 6 && x < 18){
          is >> tmp;
          string_jetsel.push_back(make_pair(cutnames.at(x),tmp) );
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;

        }
        else if ( x ==1) {is >> idlabel; cout << "idlabel=" << idlabel << endl;}
        else {
          is >> tmpf;
          float_jetsel.push_back(make_pair(cutnames.at(x),tmpf));
          cout << "Setup: float " << cutnames.at(x) << " = " <<tmpf << endl;

        }
      }
    }
    std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfJet.find(idlabel);

    if(idlabel!=""){
      
      if(fit != selectionIDMapfJet.end()){
	cerr << "Repeated ID " <<idlabel<< endl;
	exit(EXIT_FAILURE);
	
      }
    }

    selectionIDMapsJet[idlabel] = string_jetsel;
    selectionIDMapfJet[idlabel] = float_jetsel;
  }
}


void AnalyzerCore::SetupSelectionFatJet(std::string path_sel){
  Message("SetupSelectionJet", DEBUG);
  ifstream jetselconfig(path_sel.c_str());
  if(!jetselconfig) {
    cerr << "Did not find "+ path_sel+", exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  string jetline;
  int ncuts=0;
  vector<TString> cutnames;
  while(getline(jetselconfig,jetline) ){
    vector<pair<TString,TString> > string_jetsel;
    vector<pair<TString,float> > float_jetsel;
    TString idlabel;
    if (TString(jetline).Contains("webpage")) continue;
    if (TString(jetline).Contains("###")) continue;
    if (TString(jetline).Contains("ncut")) {
      std::istringstream is( jetline );
      string tmp;
      int itmp;
      is >> tmp;
      is >> itmp;
      ncuts = 2*(itmp +1);
      continue;
    }
    if (TString(jetline).Contains("ptmin")) {
      std::istringstream is( jetline );
      string tmp;
      for (int x =0; x < ncuts; x++){
        is >> tmp;
        cutnames.push_back(TString(tmp));
      }
    }
    else{
      std::istringstream is( jetline );
      string tmp;
      float tmpf;
      for (int x =0; x < ncuts; x++){
        if ( x%2 ==0) {
          is >> tmp;
          continue;
        }

        if (x > 8 && x < 20){
          is >> tmp;
          string_jetsel.push_back(make_pair(cutnames.at(x),tmp) );
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;

        }
        else if ( x ==1) {is >> idlabel;cout << "idlabel=" << idlabel << endl;}
        else {
          is >> tmpf;
          float_jetsel.push_back(make_pair(cutnames.at(x),tmpf));
          cout << "Setup: float " << cutnames.at(x) << " = " <<tmpf << endl;
	  
        }
      }
    }

    std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfFatJet.find(idlabel);

    if(idlabel!=""){
      if(fit != selectionIDMapfFatJet.end()){
	cerr << "Repeated ID " <<idlabel<< endl;
	exit(EXIT_FAILURE);
	
      }
    }

    selectionIDMapsFatJet[idlabel] = string_jetsel;
    selectionIDMapfFatJet[idlabel] = float_jetsel;
  }
}

void AnalyzerCore::SetupSelectionMuon(std::string path_sel){
  Message("SetupSelectionMuon", DEBUG);
  ifstream muonselconfig(path_sel.c_str());
  if(!muonselconfig) {
    cerr << "Did not find "+ path_sel+", exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  string muonline;
  int ncuts=0;
  vector<TString> cutnames;
  while(getline(muonselconfig,muonline) ){
    vector<pair<TString,TString> > string_muonsel;
    vector<pair<TString,float> > float_muonsel;
    TString idlabel;
    if (TString(muonline).Contains("webpage")) continue;
    if (TString(muonline).Contains("###")) continue;
    if (TString(muonline).Contains("ncut")) {
      std::istringstream is( muonline );
      string tmp;
      int itmp;
      is >> tmp;
      is >> itmp;
      ncuts = 2*(itmp +1);
      continue;
    }
    if (TString(muonline).Contains("ptmin")) {
      std::istringstream is( muonline );
      string tmp;
      for (int x =0; x < ncuts; x++){
	is >> tmp;
	cutnames.push_back(TString(tmp));
      }
    }
    else{
      std::istringstream is( muonline );
      string tmp;
      float tmpf;
      for (int x =0; x < ncuts; x++){
	if ( x%2 ==0) {
	  is >> tmp;
	  continue;
	}

	if (x > 12 && x < 18){
	  is >> tmp;
	  string_muonsel.push_back(make_pair(cutnames.at(x),tmp) );
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;

	}
	else if ( x ==1) {is >> idlabel;cout << "idlabel=" << idlabel << endl;}
	else if (x > 30 && x < 32){
          is >> tmp;
          string_muonsel.push_back(make_pair(cutnames.at(x),tmp) );
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;

        }
	else {
	  is >> tmpf;
	  float_muonsel.push_back(make_pair(cutnames.at(x),tmpf));
          cout << "Setup: float " << cutnames.at(x) << " = " <<tmpf << endl;

	}
      }
    }
    
    std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfMuon.find(idlabel);

    if(idlabel!=""){
      if(fit != selectionIDMapfMuon.end()){
	cerr << "Repeated ID " << idlabel << endl;
	exit(EXIT_FAILURE);
	
      }
    }
    selectionIDMapsMuon[idlabel] = string_muonsel;
    selectionIDMapfMuon[idlabel] = float_muonsel;
  }
}


void AnalyzerCore::SetupSelectionElectron(std::string path_sel){

  Message("SetupSelectionElectron", DEBUG);

  //// currently hard coded set of cuts

  ifstream elselconfig(path_sel.c_str());
  if(!elselconfig) {
    cerr << "Did not find "+ path_sel+", exiting ..." << endl;
    exit(EXIT_FAILURE);
  }
  string elline;
  vector<TString> cutnames;
  int ncuts=0;
  while(getline(elselconfig,elline) ){
    vector<pair<TString,TString> > string_elsel;
    vector<pair<TString,float> > float_elsel;
    TString idlabel;
    if (TString(elline).Contains("webpage")) continue;
    if (TString(elline).Contains("###")) continue;
    if (TString(elline).Contains("ncut")) {
      std::istringstream is( elline );
      string tmp;
      int itmp;
      is >> tmp;
      is >> itmp;
      ncuts = 2*(itmp +1);
      continue;
    }

    if (TString(elline).Contains("ptmin")) {
      std::istringstream is( elline );
      string tmp;
      for (int x =0; x < ncuts; x++){
        is >> tmp;
        cutnames.push_back(TString(tmp));
      }
    }
    else{
      std::istringstream is( elline );
      string tmp;
      float tmpf;

      for (int x =0; x < ncuts; x++){
        if ( x%2 ==0) {
          is >> tmp;
          continue;
        }

        if (x > 10 && x < 18){
          is >> tmp;
	  cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;

          string_elsel.push_back(make_pair(cutnames.at(x),tmp) );
        }
	else  if (x > 26 && x < 30){
          is >> tmp;
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;
          string_elsel.push_back(make_pair(cutnames.at(x),tmp) );
        }
	else  if (x > 34 && x < 46){
          is >> tmp;
          cout << "Setup: string " << cutnames.at(x) << " = " <<tmp << endl;
          string_elsel.push_back(make_pair(cutnames.at(x),tmp) );
        }
        else if ( x ==1) {is >> idlabel;cout << "idlabel=" << idlabel << endl;}
        else {
          is >> tmpf;
          cout << "Setup: float " << cutnames.at(x) << " = " <<tmpf << endl;

          float_elsel.push_back(make_pair(cutnames.at(x),tmpf));
        }
      }
    }

    std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfElectron.find(idlabel);

    if(idlabel!=""){
      if(fit != selectionIDMapfElectron.end()){
	cerr << "Repeated ID " <<idlabel<< endl;
	exit(EXIT_FAILURE);
	
      }
    }
    selectionIDMapsElectron[idlabel] = string_elsel;
    selectionIDMapfElectron[idlabel] = float_elsel;
  }
}




std::map<TString,BTagSFUtil*> AnalyzerCore::SetupBTagger(std::vector<TString> taggers, std::vector<TString> wps){

  //// Btagging code for Moriond17 samples
  
  //// Current use is for HN analyses only.
  //// HN analysis uses method 2 a) from twiki:
  //// https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
  

  /// If users want to use another method please add additional function and modify BTag/ dir in LQLAnalysis/src/

  std::map<TString,BTagSFUtil*>  tmpmap;
  for(std::vector<TString>::const_iterator it = taggers.begin(); it != taggers.end(); it++){
    for(std::vector<TString>::const_iterator it2 = wps.begin(); it2 != wps.end(); it2++){
      if (it->Contains("CSVv2")){
        tmpmap[*it + "_" + *it2 + "_lf"]          = new BTagSFUtil("incl"  , (*it + "_BtoF").Data(), (*it + "_GtoH").Data(), it2->Data());
        tmpmap[*it + "_" + *it2 + "_hf"]          = new BTagSFUtil("mujets", (*it + "_BtoF").Data(), (*it + "_GtoH").Data(), it2->Data());
        tmpmap[*it + "_" + *it2 + "_lf_systup"]   = new BTagSFUtil("incl"  , (*it + "_BtoF").Data(), (*it + "_GtoH").Data(), it2->Data(),  3);  
        tmpmap[*it + "_" + *it2 + "_hf_systup"]   = new BTagSFUtil("mujets", (*it + "_BtoF").Data(), (*it + "_GtoH").Data(), it2->Data(),  1);  
        tmpmap[*it + "_" + *it2 + "_lf_systdown"] = new BTagSFUtil("incl"  , (*it + "_BtoF").Data(), (*it + "_GtoH").Data(), it2->Data(), -3); 
        tmpmap[*it + "_" + *it2 + "_hf_systdown"] = new BTagSFUtil("mujets", (*it + "_BtoF").Data(), (*it + "_GtoH").Data(), it2->Data(), -1);
	// tmpmap[*it +  "_" + *it2 + "_hfcomb"]= new BTagSFUtil("comb", it->Data(), it2->Data());                /// SWITCH ON IF USER NEEDS THIS METHOD  
	// tmpmap[*it +  "_" + *it2 + "iterativefit"]= new BTagSFUtil("iterativefit", it->Data(), it2->Data());   /// SWITCH ON IF USER NEEDS THIS METHOD
      }
      if (it->Contains("cMVA")){
	tmpmap[*it  + "_" + *it2 + "_lf"]= new BTagSFUtil("incl", (*it + "_BtoF").Data(),  (*it + "_GtoH").Data(), it2->Data());
	tmpmap[*it  +  "_" + *it2 + "_hf"]= new BTagSFUtil("ttbar", (*it + "_BtoF").Data(),  (*it + "_GtoH").Data(), it2->Data());

	// tmpmap[*it +  "_" + *it2 + "iterativefit"]= new BTagSFUtil("iterativefit", it->Data(), it2->Data());   /// SWITCH ON IF USER NEEDS THIS METHOD       
      }
      if (it->Contains("DeepCVS")){
	/// TO ADD If requested
      }
    }
  }

  return tmpmap;
}


//################################################################################################  
//@@@@@@@@@@@@@@@@@@@  ANALYSIS FUNCTIONALITY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    

float AnalyzerCore::GetDiLepMass(std::vector<snu::KElectron> electrons){

  if(electrons.size() != 2) return 0.;
  snu::KParticle p = electrons.at(0) + electrons.at(1);
  return p.M();
}

float AnalyzerCore::GetDiLepMass(std::vector<snu::KMuon> muons){

  if(muons.size() != 2) return 0.;
  snu::KParticle p = muons.at(0) + muons.at(1);
  return p.M();
}



float AnalyzerCore::GetMasses(TString svariable, std::vector<snu::KElectron> electrons, std::vector<snu::KJet> jets,  std::vector<snu::KFatJet> fatjets, vector<int> ijets, bool lowmass){
  if(electrons.size() != 2) return 0.;


  // variable 1 = lljj                                                                                                                                                               
  // variable 2 = l1jj                                                                                                                                                               
  // variable 3 = l2jj                                                                                                                                                               
  // variable 4 = llj                                                                                                                                                                
  // variable 5 = jj                                                                                                                                                                 
  // variable 6 = contra JJ mass                                                                                                                                                     

  int variable (-1);
  if(svariable == "lljj") variable = 1;
  else if(svariable == "l1jj") variable = 2;
  else if(svariable == "l2jj") variable = 3;
  else if(svariable == "llj") variable = 4;
  else if(svariable == "l1j") variable = 7;
  else if(svariable == "l2j") variable = 8;
  else if(svariable == "jj") variable = 5;
  else if(svariable == "contMT") variable = 6;
  else if(svariable == "llfj") variable = -1;
  else if(svariable == "l1fj") variable = -2;
  else if(svariable == "l2fj") variable = -3;
  else if(svariable == "fj") variable = -4;
  else return -999.;

  snu::KFatJet fatjet;
  float dMFatJet=9999.;
  for(UInt_t emme=0; emme<fatjets.size(); emme++){
    if(fabs(fatjets[emme].PrunedMass() -  80.4) < dMFatJet){
      dMFatJet=fatjets[emme].PrunedMass();
      fatjet=fatjets[emme];
    }
  }
  fatjet.SetPtEtaPhiM(fatjet.Pt(), fatjet.E(), fatjet.Phi(), fatjet.PrunedMass());
  if(variable==-1) return (electrons[0] + electrons[1] + fatjet).M();
  if(variable==-2) return (electrons[0] + fatjet).M();
  if(variable==-3) return (electrons[1] + fatjet).M();
  if(variable==-4) return fatjet.PrunedMass();

  if(jets.size() == 1){
    if(variable==4) return (electrons[0] + electrons[1] + jets[0]).M();
    if(variable==7) return (electrons[0]  + jets[0]).M();
    if(variable==8) return (electrons[1] + jets[0]).M();

  }
  if(jets.size() < 2) return -999.;


  float dijetmass_tmp=999.;
  float dijetmass=9990000.;
  int m=-999;
  int n=-999;
  for(UInt_t emme=0; emme<jets.size(); emme++){
    for(UInt_t enne=1; enne<jets.size(); enne++) {
      if(emme == enne) continue;
      if(lowmass)   dijetmass_tmp = (jets[emme]+jets[enne]+electrons[0] + electrons[1]).M();
      else dijetmass_tmp = (jets[emme]+jets[enne]).M();
      if ( fabs(dijetmass_tmp-80.4) < fabs(dijetmass-80.4) ) {
        dijetmass = dijetmass_tmp;
        m = emme;
        n = enne;
      }
    }
  }

  if(ijets.size() ==2){
    if(ijets[0] != 0){
      ijets.push_back(m);
      ijets.push_back(n);
    }
  }

  if(variable==1) return (electrons[0] + electrons[1] + jets[m]+jets[n]).M();
  if(variable==2) return (electrons[0]  + jets[m]+jets[n]).M();
  if(variable==3) return (electrons[1] + jets[m]+jets[n]).M();
  if(variable==5) return (jets[m]+jets[n]).M();
  if(variable==6) {
    float dPhi = fabs(TVector2::Phi_mpi_pi(jets[m].Phi() - jets[n].Phi()));
    float contramass=2*jets[m].Pt()*jets[n].Pt()*(1+cos(dPhi));
    contramass=sqrt(contramass);
    return contramass;
  }

  return 0.;



}

float AnalyzerCore::GetMasses(TString svariable, std::vector<snu::KMuon> muons, std::vector<snu::KJet> jets,  std::vector<snu::KFatJet> fatjets, vector<int> ijets, bool lowmass){
  
  if(muons.size() != 2) return 0.;


  // variable 1 = lljj
  // variable 2 = l1jj
  // variable 3 = l2jj
  // variable 4 = llj
  // variable 5 = jj
  // variable 6 = contra JJ mass

  int variable (-1);
  if(svariable == "lljj") variable = 1;
  else if(svariable == "l1jj") variable = 2;
  else if(svariable == "l2jj") variable = 3;
  else if(svariable == "llj") variable = 4;
  else if(svariable == "l1j") variable = 7;
  else if(svariable == "l2j") variable = 8;
  else if(svariable == "jj") variable = 5;
  else if(svariable == "contMT") variable = 6;
  else if(svariable == "llfj") variable = -1;
  else if(svariable == "l1fj") variable = -2;
  else if(svariable == "l2fj") variable = -3;
  else if(svariable == "fj") variable = -4;
  else return -999.;

  snu::KFatJet fatjet;
  float dMFatJet=9999.;
  for(UInt_t emme=0; emme<fatjets.size(); emme++){
    if(fabs(fatjets[emme].PrunedMass() -  80.4) < dMFatJet){
      dMFatJet=fatjets[emme].PrunedMass();
      fatjet=fatjets[emme];
    }
  }

  if(variable==-1) return (muons[0] + muons[1] + fatjet).M();
  if(variable==-2) return (muons[0] + fatjet).M();
  if(variable==-3) return (muons[1] + fatjet).M();
  if(variable==-4) return fatjet.PrunedMass();

  if(jets.size() == 1){
    if(variable==4) return (muons[0] + muons[1] + jets[0]).M();
    if(variable==7) return (muons[0]  + jets[0]).M();
    if(variable==8) return (muons[1] + jets[0]).M();
    
  }
  if(jets.size() < 2) return -999.;


  float dijetmass_tmp=999.;
  float dijetmass=9990000.;
  int m=-999;
  int n=-999;
  for(UInt_t emme=0; emme<jets.size(); emme++){
    for(UInt_t enne=1; enne<jets.size(); enne++) {
      if(emme == enne) continue;
      if(lowmass)   dijetmass_tmp = (jets[emme]+jets[enne]+muons[0] + muons[1]).M();
      else dijetmass_tmp = (jets[emme]+jets[enne]).M();
      if ( fabs(dijetmass_tmp-80.4) < fabs(dijetmass-80.4) ) {
	dijetmass = dijetmass_tmp;
	m = emme;
	n = enne;
      }
    }
  }
  
  if(ijets.size() ==2){
    if(ijets[0] != 0){
      ijets.push_back(m);
      ijets.push_back(n);
    }
  }

  if(variable==1) return (muons[0] + muons[1] + jets[m]+jets[n]).M();
  if(variable==2) return (muons[0]  + jets[m]+jets[n]).M();
  if(variable==3) return (muons[1] + jets[m]+jets[n]).M();
  if(variable==5) return (jets[m]+jets[n]).M();
  if(variable==6) {
    float dPhi = fabs(TVector2::Phi_mpi_pi(jets[m].Phi() - jets[n].Phi()));
    float contramass=2*jets[m].Pt()*jets[n].Pt()*(1+cos(dPhi));
    contramass=sqrt(contramass);
    return contramass;
  }
  
  return 0.;
}



bool AnalyzerCore::EtaRegion(TString reg,  std::vector<snu::KElectron> electrons){
  if(electrons.size() != 2) return false;

  if(reg.Contains("EE")){
    if(fabs(electrons.at(0).Eta()) < 1.5) return false;
    if(fabs(electrons.at(1).Eta()) < 1.5) return false;
    return true;
  }
  if(reg.Contains("BB")){
    if(fabs(electrons.at(0).Eta()) > 1.5) return false;
    if(fabs(electrons.at(1).Eta()) > 1.5) return false;
    return true;
  }
  if(reg.Contains("EB")){
    if(fabs(electrons.at(0).Eta()) > 1.5){
      if(fabs(electrons.at(1).Eta()) > 1.5) return false;
    }
    if(fabs(electrons.at(0).Eta()) < 1.5){
      if(fabs(electrons.at(1).Eta()) < 1.5) return false;
    }
    return true;
  }
  return false;
  
}


bool AnalyzerCore::EtaRegion(TString reg,  std::vector<snu::KMuon> muons){
  
  if(muons.size() != 2) return false;
  if(reg.Contains("EE")){
    if(muons.at(0).Eta() < 1.5) return false;
    if(muons.at(1).Eta() < 1.5) return false;
    return true;
  }
  if(reg.Contains("BB")){
    if(muons.at(0).Eta() > 1.5) return false;
    if(muons.at(1).Eta() > 1.5) return false;
    return true;
  }
  if(reg.Contains("EB")){
    if(muons.at(0).Eta() > 1.5){
      if(muons.at(1).Eta() > 1.5) return false;
    }
    if(muons.at(0).Eta() < 1.5){
      if(muons.at(1).Eta() < 1.5) return false;
    }
    return true;
  }
  return false;
}



//################################################################################################ 
//@@@@@@@@@@@@@@@@@@@  GET SKTREE  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    


std::vector<snu::KFatJet> AnalyzerCore::GetFatJets(BaseSelection::ID jetid, float ptcut, float etacut){
  return GetFatJets(GetStringID(jetid), ptcut, etacut);
}


std::vector<snu::KJet> AnalyzerCore::GetJets(BaseSelection::ID jetid, float ptcut, float etacut){
  return GetJets(GetStringID(jetid), ptcut, etacut);
}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(BaseSelection::ID muonid, float ptcut, float etacut){
  return GetMuons(GetStringID(muonid), ptcut, etacut);
}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(BaseSelection::ID muonid, bool keepfake, float ptcut, float etacut){
  return GetMuons(GetStringID(muonid), keepfake, ptcut, etacut);
}

std::vector<snu::KElectron> AnalyzerCore::GetElectrons(BaseSelection::ID electronid, float ptcut, float etacut){
  return GetElectrons(true, true, GetStringID(electronid), ptcut, etacut);
}




std::vector<snu::KElectron> AnalyzerCore::GetElectrons(bool keepcf, bool keepfake,   BaseSelection::ID electronid, float ptcut, float etacut){
  return GetElectrons(keepcf, keepfake, GetStringID(electronid), ptcut, etacut);
}



TString AnalyzerCore::GetStringID(BaseSelection::ID id){
  if(id == BaseSelection::ELECTRON_POG_VETO ) return "ELECTRON_POG_VETO";
  if(id == BaseSelection::ELECTRON_POG_LOOSE) return "ELECTRON_POG_LOOSE";
  if(id == BaseSelection::ELECTRON_POG_MEDIUM) return "ELECTRON_POG_MEDIUM";
  if(id == BaseSelection::ELECTRON_POG_TIGHT           ) return "ELECTRON_POG_TIGHT";
  if(id == BaseSelection::ELECTRON_POG_MVATrig) return "ELECTRON_POG_MVATrig";
  if(id == BaseSelection::ELECTRON_POG_MVANonTrig    ) return "ELECTRON_POG_MVANonTrig";
  if(id == BaseSelection::ELECTRON_ECAL_FIDUCIAL     ) return "ELECTRON_ECAL_FIDUCIAL";
  if(id == BaseSelection::ELECTRON_HN_VETO           ) return "ELECTRON_HN_VETO";
  if(id == BaseSelection::ELECTRON_HN_TIGHT          ) return "ELECTRON_HN_TIGHT";
  if(id == BaseSelection::ELECTRON_HN_FAKELOOSE      ) return "ELECTRON_HN_FAKELOOSE";
  if(id == BaseSelection::ELECTRON_HN_FAKELOOSE_NOD0 ) return "ELECTRON_HN_FAKELOOSE_NOD0";
  if(id == BaseSelection::ELECTRON_TOP_VETO          ) return "ELECTRON_TOP_VETO";
  if(id == BaseSelection::ELECTRON_TOP_LOOSE         ) return "ELECTRON_TOP_LOOSE";
  if(id == BaseSelection::ELECTRON_TOP_TIGHT         ) return "ELECTRON_TOP_TIGHT";
  if(id == BaseSelection::ELECTRON_PTETA             ) return "ELECTRON_PTETA";
  if(id == BaseSelection::ELECTRON_NOCUT             ) return "ELECTRON_NOCUT";
  if(id == BaseSelection::MUON_POG_LOOSE             ) return "MUON_POG_LOOSE";
  if(id == BaseSelection::MUON_POG_MEDIUM            ) return "MUON_POG_MEDIUM";
  if(id == BaseSelection::MUON_POG_TIGHT             ) return "MUON_POG_TIGHT";
  if(id == BaseSelection::MUON_HN_VETO               ) return "MUON_HN_VETO";
 if(id == BaseSelection::MUON_HN_FAKELOOSE          ) return "MUON_HN_FAKELOOSE";
 if(id == BaseSelection::MUON_HN_TIGHT              ) return "MUON_HN_TIGHT";
 if(id == BaseSelection::MUON_FAKELOOSE             ) return "MUON_FAKELOOSE";
 if(id == BaseSelection::MUON_TOP_VETO              ) return "MUON_TOP_VETO";
 if(id == BaseSelection::MUON_TOP_LOOSE             ) return "MUON_TOP_LOOSE";
 if(id == BaseSelection::MUON_TOP_TIGHT             ) return "MUON_TOP_TIGHT";
 if(id == BaseSelection::MUON_PTETA                 ) return "MUON_PTETA";
 if(id == BaseSelection::MUON_NOCUT                 ) return "MUON_NOCUT";
 if(id == BaseSelection::PFJET_LOOSE                ) return "PFJET_LOOSE";
 if(id == BaseSelection::PFJET_MEDIUM               ) return "PFJET_MEDIUM";
 if(id == BaseSelection::PFJET_TIGHT                ) return "PFJET_TIGHT";
 if(id == BaseSelection::JET_HN                     ) return "JET_HN";
 if(id == BaseSelection::JET_HN_TChannel            ) return "JET_HN_TChannel";
 if(id == BaseSelection::JET_NOLEPTONVETO           ) return "JET_NOLEPTONVETO";
 if(id == BaseSelection::JET_LOOSE                  ) return "JET_LOOSE";
 if(id == BaseSelection::JET_TIGHT                  ) return "JET_TIGHT";
 if(id == BaseSelection::PHOTON_POG_LOOSE           ) return "PHOTON_POG_LOOSE";
 if(id == BaseSelection::PHOTON_POG_MEDIUM          ) return "PHOTON_POG_MEDIUM";
 if(id == BaseSelection::PHOTON_POG_TIGHT           ) return "PHOTON_POG_TIGHT";
 cerr << " ID [--] not found" << endl; exit(EXIT_FAILURE);

}

std::vector<snu::KJet> AnalyzerCore::GetJets(TString jetid,float ptcut, float etacut){
  
  std::vector<snu::KJet> jetColl;
  
  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsJet.find(jetid);
  std::map<TString, vector<pair<TString,float > > >::iterator fit = selectionIDMapfJet.find(jetid);
  if(it== selectionIDMapsJet.end()){
    cerr << "Jet ID ["+jetid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {
    TString muontag="";
    TString eltag="";
    for (unsigned int i=0; i  < it->second.size(); i++){
      if ( it->second.at(i).first == "remove_near_muonID") muontag =  it->second.at(i).second;
      if ( it->second.at(i).first == "remove_near_electronID") eltag =  it->second.at(i).second;
    }
    
    if (muontag.Contains("NONE") && eltag.Contains("NONE"))  eventbase->GetJetSel()->SelectJets(jetColl, it->second, fit->second, ptcut,etacut);
    else if (muontag.Contains("NONE") || eltag.Contains("NONE")) {    cerr << "cannot choose to remove jets near only one lepton" << endl; exit(EXIT_FAILURE);}
    else eventbase->GetJetSel()->SelectJets(jetColl, GetMuons(muontag), GetElectrons(eltag) , it->second, fit->second, ptcut,etacut);
  }
  return jetColl;
  
}

std::vector<snu::KJet> AnalyzerCore::GetJetsWFT(TString jetid,TString fatjetid,float ptcut, float etacut){

  std::vector<snu::KJet> jetColl;

  std::vector<snu::KFatJet> fatjetColl = GetFatJets(fatjetid);
  
  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsJet.find(jetid);
  std::map<TString, vector<pair<TString,float > > >::iterator fit = selectionIDMapfJet.find(jetid);
  if(it== selectionIDMapsJet.end()){
    cerr << "Jet ID ["+jetid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {
    TString muontag="";
    TString eltag="";
    for (unsigned int i=0; i  < it->second.size(); i++){
      if ( it->second.at(i).first == "remove_near_muonID") muontag =  it->second.at(i).second;
      if ( it->second.at(i).first == "remove_near_electronID") eltag =  it->second.at(i).second;
    }

    if (muontag.Contains("NONE") && eltag.Contains("NONE"))  eventbase->GetJetSel()->SelectJets(jetColl, it->second, fit->second, ptcut,etacut);
    else if (muontag.Contains("NONE") || eltag.Contains("NONE")) {    cerr << "cannot choose to remove jets near only one lepton" << endl; exit(EXIT_FAILURE);}
    else eventbase->GetJetSel()->SelectJets(jetColl, fatjetColl,GetMuons(muontag), GetElectrons(eltag) , it->second, fit->second, ptcut,etacut);
  }
  return jetColl;

}




std::vector<snu::KFatJet> AnalyzerCore::GetFatJets(TString fatjetid, float ptcut, float etacut){

  std::vector<snu::KFatJet> fatjetColl;

  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsFatJet.find(fatjetid);
  std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfFatJet.find(fatjetid);
  if(it== selectionIDMapsFatJet.end()){
    cerr << "FatJet ID ["+fatjetid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {
    TString muontag="";
    TString eltag="";
    for (unsigned int i=0; i  < it->second.size(); i++){
      if ( it->second.at(i).first == "remove_near_muonID") muontag =  it->second.at(i).second;
      if ( it->second.at(i).first == "remove_near_electronID") eltag =  it->second.at(i).second;
    }
    if (muontag.Contains("NONE") && eltag.Contains("NONE"))  eventbase->GetFatJetSel()->SelectFatJets(fatjetColl, it->second, fit->second, ptcut,etacut);
    else if (muontag.Contains("NONE") || eltag.Contains("NONE")) {    cerr << "cannot choose to remove jets near only one lepton" << endl; exit(EXIT_FAILURE);}
    else eventbase->GetFatJetSel()->SelectFatJets(fatjetColl, GetMuons(muontag), GetElectrons(eltag) ,  it->second, fit->second, ptcut,etacut);
    
  }

  return fatjetColl;

}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(TString muid, float ptcut, float etacut){
  return GetMuons(muid, true, ptcut, etacut);
}

std::vector<snu::KMuon> AnalyzerCore::GetMuons(TString muid, bool keepfakes, float ptcut, float etacut){

  std::vector<snu::KMuon> muonColl;
  
  if(muid.Contains("NONE")) return muonColl;

  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsMuon.find(muid);
  std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfMuon.find(muid);
  
  if(1){
    //// This method was 10% faster in processing 40K MC events accessing 1000 vectors                                                                                          
    if(it== selectionIDMapsMuon.end()){
      cerr << "Muon ID ["+muid+"] not found" << endl; exit(EXIT_FAILURE);
    }
    if(fit== selectionIDMapfMuon.end()){
      cerr << "Muon ID ["+muid+"] not found" << endl; exit(EXIT_FAILURE);
    }

    if (ptcut == -999.)  eventbase->GetMuonSel()->SelectMuons(muonColl ,muid, it->second, fit->second);
    else eventbase->GetMuonSel()->SelectMuons(muonColl ,muid, it->second, fit->second, ptcut, etacut);
    return  GetTruePrompt(muonColl, keepfakes);
    
  }
  
  return  GetTruePrompt(muonColl, keepfakes);
  
}


std::vector<snu::KElectron> AnalyzerCore::GetElectrons(TString elid,float ptcut, float etacut){

  if(k_classname.Contains("HNDiElectronOpt")){
    if(k_running_chargeflip)  return GetElectrons( true, false,   elid, ptcut, etacut);
    return GetElectrons(false, false,elid, ptcut, etacut);
  }
  /// if cf flag set and MC keep CF and conversion electrons
  if(k_running_chargeflip)  return GetElectrons( true, false, elid, ptcut, etacut);
  
  
  //// by default keep all electrons
  return GetElectrons(true,  true ,elid, ptcut, etacut);

}



std::vector<snu::KElectron> AnalyzerCore::GetElectrons(bool keepcf, bool keepfake,  TString elid,float ptcut, float etacut){

  std::vector<snu::KElectron> electronColl;
  
  if(elid.Contains("NONE")) return electronColl;

  
  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsElectron.find(elid);
  std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfElectron.find(elid);
  
  if(1){
    //// This method was 10% faster in processing 40K MC events accessing 1000 vectors
    if(it== selectionIDMapsElectron.end()){
      cerr << "Electron ID ["+elid+"] not found" << endl; exit(EXIT_FAILURE);
    }
    if(fit== selectionIDMapfElectron.end()){
      cerr << "Electron ID ["+elid+"] not found" << endl; exit(EXIT_FAILURE);
    }
    
    if(ptcut == -999.)eventbase->GetElectronSel()->SelectElectrons(electronColl,elid, it->second, fit->second);
    else eventbase->GetElectronSel()->SelectElectrons(electronColl,elid, it->second, fit->second,ptcut, etacut);
    
    return  GetTruePrompt(electronColl, keepcf, keepfake);
  }

  else{
    if(it== selectionIDMapsElectron.end()){
      cerr << "Electron ID ["+elid+"] not found" << endl; exit(EXIT_FAILURE);
    }
    else {
      
      bool check_cc(false);
      bool check_cv(false);
      TString el_id="";
      vector<pair<TString,TString> > v_string =  it->second;
      for(unsigned int iv=0; iv < v_string.size(); iv++){
	if(v_string[iv].second == "false") continue;
	if(v_string[iv].first.Contains("(POG)")) el_id=v_string[iv].first;
	if(v_string[iv].first.Contains("(MVA)") ) el_id=v_string[iv].first;
	if(v_string[iv].first == "GsfCtfScPix") check_cc=true;
	if(v_string[iv].first == "convveto") check_cv=true;
      }
      
      
      vector<pair<TString,float> > v_float =  fit->second;
      float isomax_b(-999.); 
      float isomax_e (-999.);
      float dxymax_b (-999.);
      float dxymax_e (-999.);
      float dzmax_b (-999.);
      float dzmax_e (-999.);
      
      float dxysigmax(-999.);
      float dxysigmin(-999.);
      
      for(unsigned int iv=0; iv < v_float.size(); iv++){
	if(!Check(v_float[iv].second)) continue;
	if(v_float[iv].first == "ptmin") {
	  if (ptcut == -999.) ptcut =v_float[iv].second;
	}
	if(v_float[iv].first == "|etamax|"){
	  if (etacut == -999.) etacut = v_float[iv].second;
	}
	if(v_float[iv].first == "isomax03_b")isomax_b =v_float[iv].second;
	if(v_float[iv].first == "isomax03_e") isomax_e=v_float[iv].second;
	if(v_float[iv].first == "|dxymax_b|") dxymax_b=v_float[iv].second;
	if(v_float[iv].first == "|dxymax_e|") dxymax_e=v_float[iv].second;
	if(v_float[iv].first == "|dzmax_b|") dzmax_b=v_float[iv].second;
	if(v_float[iv].first == "|dzmax_e|") dzmax_e=v_float[iv].second;
	if(v_float[iv].first == "|dxysigmax|") dxysigmax=v_float[iv].second;
	if(v_float[iv].first == "|dxysigmin|") dxysigmin=v_float[iv].second;
      }
      
      
      eventbase->GetElectronSel()->SelectElectrons(electronColl,elid, el_id,check_cc,check_cv, isomax_b,isomax_e,dxymax_b,dxymax_e,dzmax_b,dzmax_e, dxysigmax,dxysigmin, ptcut, etacut);
    }
  }
  
  //if(elid == "ELECTRON_HN_TIGHT"){
    /// This is the vector of electrons with optimie cuts
    //std::vector<snu::KElectron> _electronColl;
    //if(k_running_nonprompt) eventbase->GetElectronSel()->SelectElectrons(_electronColl, BaseSelection::ELECTRON_HN_FAKELOOSE_NOD0, 15., 2.5);
    //else eventbase->GetElectronSel()->SelectElectrons(_electronColl,BaseSelection::ELECTRON_HN_TIGHT, 15., 2.5);
    //electronColl =ShiftElectronEnergy(_electronColl, k_running_chargeflip);
  // }
   

  return  GetTruePrompt(electronColl, keepcf, keepfake); 

}




bool AnalyzerCore::HasCloseBJet(snu::KElectron el, KJet::Tagger tag, KJet::WORKING_POINT wp, int period){

  std::vector<snu::KJet> alljets = GetJets("JET_NOLEPTONVETO");

  if(period == 0) {
    Message("period not set in AnalyzerCore::HasCloseBJet. Will assign mcperiod for you but this may not give correct behaviour", WARNING);
    period=GetPeriod();
  }

  bool cl = false;
  for(unsigned int ij =0; ij < alljets.size(); ij++){
    if(el.DeltaR(alljets.at(ij)) < 0.5){

      if(IsBTagged(alljets.at(ij), tag, wp, period))cl = true;

    }
  }

  return cl;
}



bool AnalyzerCore::TriggerMatch(TString trigname, vector<snu::KMuon> mu){
  
  if(mu.size() == 2){
    if(!mu.at(0).TriggerMatched(trigname)) return false;
    if(!mu.at(1).TriggerMatched(trigname)) return false;
  }
  return true;
}

bool AnalyzerCore::Is2015Analysis(){
  if(TString(getenv("running2015")).Contains("True")) return true;
  else return false;
}



float AnalyzerCore::WeightByTrigger(vector<TString> triggernames, float tlumi){
  
  if(isData){
    for(unsigned int i=0; i < triggernames.size() ; i++){
      //// code here sets weight to -99999. if user tries to use incorrect dataset in data 
      /// datasets for each trigger in 2015 can be seen in following googledoc:
      //// https://docs.google.com/spreadsheets/d/1BkgAHCC4UtP5sddTZ5G5iWY16BxleuK7rqT-Iz2LHiM/pubhtml?gid=0&single=true
       if(!k_channel.Contains("DoubleMuon")){
	if(triggernames.at(i).Contains("HLT_Mu17_Mu8_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_Mu8_SameSign_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu20_Mu10_SameSign_DZ_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu8_TrkIsoVVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_TripleMu_12_10_5")) return -99999.;
      }
      if(!k_channel.Contains("SingleMuon")){
	if(triggernames.at(i).Contains("HLT_Mu8_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu20_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu50_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoMu22_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoTkMu22_")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoMu20_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoTkMu20_")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoMu24_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_IsoTkMu24_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu24_eta2p1_v" )) return -99999.;
      }
      if(!k_channel.Contains("DoubleEG")){
	if(triggernames.at(i).Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;
        if(triggernames.at(i).Contains("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return -99999.;

      }
      if(!k_channel.Contains("SingleElectron")){
	// Single Electon                                                                                                                                                       
	if(triggernames.at(i).Contains("HLT_Ele23_WPLoose_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele27_WPTight_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele27_eta2p1_WPLoose_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele25_eta2p1_WPTight_Gsf_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Ele45_WPLoose_Gsf_v")) return -99999.;
      }
      if(!k_channel.Contains("MuonEG")){
	if(triggernames.at(i).Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"))  return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v")) return -99999.;
	if(triggernames.at(i).Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v")) return -99999.;
      }
    }
    return 1.;
  }


  float trigps= -1.;
  for(unsigned int i=0; i < triggernames.size() ; i++){
    if(WeightByTrigger(triggernames.at(i), tlumi) > trigps) trigps = WeightByTrigger(triggernames.at(i),tlumi) ;
  }

  //if(trigps  > 1.) {m_logger << ERROR << "Error in getting weight for trigger prescale. It cannot be > 1, this means trigger lumi >> total lumi"  << LQLogger::endmsg; exit(0);}
  if(trigps  < 0.) {m_logger << ERROR << "Error in getting weight for trigger " << triggernames.at(0) << " prescale. It cannot be < 0, this means trigger lumi >> total lumi"  << LQLogger::endmsg; exit(0);}
  
  return trigps;
 
}

float AnalyzerCore::WeightByTrigger(TString triggername, float tlumi){

  /// Function applies weight to MC 
  /// Depends on trigger 
  /// 

  if(isData) return 1.;
  //  brilcalc lumi -u /pb 
  // --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json 
  // -i jsonfiles/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt 
  //--hltpath ""
  
  /// NUMBERS FROM GetLumi_Triggers.py SCRIPT
  
  /// In v766 path lumi is corrected for removal of bad beamspot LS
  // https://github.com/vallot/CATTools/commit/aae3e60b194b1bacf2595a33c8fa27f411dac16b
  for(map<TString, float>::iterator mit = trigger_lumi_map_cat2016.begin(); mit != trigger_lumi_map_cat2016.end(); mit++){
    float corr_trig = GetTriggerPrescaleCorrection(triggername);
    if(triggername.Contains(mit->first)) return (corr_trig*mit->second / tlumi);
  }
  m_logger << ERROR << "Error in getting weight for trigger  " << triggername << "  prescale. Trigname is not correct or not in map"  << LQLogger::endmsg; exit(0);

  return 1.;
}



void AnalyzerCore::AddTriggerToList(TString triggername){
  
  triggerlist.push_back(triggername);
}

AnalyzerCore::~AnalyzerCore(){
  

  Message("In AnalyzerCore Destructor" , INFO);

  trigger_lumi_map_cat2016.clear();
  compmap.clear();
  for(map<TString, TH1*>::iterator it = maphist.begin(); it!= maphist.end(); it++){
    delete it->second;
  }
  maphist.clear();

  for(map<TString, TH2*>::iterator it = maphist2D.begin(); it!= maphist2D.end(); it++){
    delete it->second;
  }
  maphist2D.clear();
  
  for(map<TString, TH3*>::iterator it = maphist3D.begin(); it!= maphist3D.end(); it++){
    delete it->second;
  }
  maphist3D.clear();


  for(map<TString, MuonPlots*>::iterator it = mapCLhistMu.begin(); it != mapCLhistMu.end(); it++){
    delete it->second;
  }
  mapCLhistMu.clear();
  

  for(map<TString, JetPlots*>::iterator it = mapCLhistJet.begin(); it != mapCLhistJet.end(); it++){
    delete it->second;
  }
  mapCLhistJet.clear();

  for(map<TString, ElectronPlots*>::iterator it = mapCLhistEl.begin(); it != mapCLhistEl.end(); it++){
    delete it->second;
  }
  mapCLhistEl.clear();

  for(map<TString, SignalPlotsEE*>::iterator it = mapCLhistSigEE.begin(); it != mapCLhistSigEE.end(); it++){
    delete it->second;
  }
  mapCLhistSigEE.clear();

  for(map<TString, SignalPlotsMM*>::iterator it = mapCLhistSigMM.begin(); it != mapCLhistSigMM.end(); it++){
    delete it->second;
  }
  mapCLhistSigMM.clear();


  for(map<TString, SignalPlotsEM*>::iterator it = mapCLhistSigEM.begin(); it != mapCLhistSigEM.end(); it++){
    delete it->second;
  }
  mapCLhistSigEM.clear();

  
  for(map<TString, TriLeptonPlots*>::iterator it = mapCLhistTriLep.begin(); it != mapCLhistTriLep.end(); it++){
    delete it->second;
  }
  mapCLhistTriLep.clear();
  
  for(map<TString, HNpairPlotsMM*>::iterator it = mapCLhistHNpairMM.begin(); it != mapCLhistHNpairMM.end(); it++){
    delete it->second;
  }
  mapCLhistHNpairMM.clear();

  for(map<TString,TNtupleD*>::iterator it = mapntp.begin(); it!= mapntp.end(); it++){ 
    delete it->second;
  }
  mapntp.clear();


  if(!k_classname.Contains("SKTreeMaker")){
    for(std::map<TString,BTagSFUtil*>::iterator it = MapBTagSF.begin(); it!= MapBTagSF.end(); it++){
      delete it->second;
    }
    MapBTagSF.clear();
  }

  for(map<TString, HNTriLeptonPlots*>::iterator it = mapCLhistHNTriLep.begin(); it != mapCLhistHNTriLep.end(); it++){
    delete it->second;
  }
  mapCLhistHNTriLep.clear();

  //// New class functions for databkg+corrections
  if(k_classname == "SKTreeMaker")   delete mcdata_correction;
  if(!k_classname.Contains("SKTreeMaker")){
    delete mcdata_correction;
  }

  if(k_running_nonprompt || k_running_chargeflip)  delete m_datadriven_bkg;

  //==== HN Gen Matching
  if((TString(getenv("USER")) == "jskim" || TString(getenv("USER")) =="shjeon")){
    delete m_HNgenmatch;
  }
}

//###
//###   IMPORTANT BASE FUNCTION: SETS UP EVENT FOR ALL CYCLES
//###

void AnalyzerCore::SetupID(){

  if(IDSetup) return;

  // IF STANDARD SKTREE CODE NO NEED FOR ID
  if(!k_classname.Contains("HN")){
    if(k_classname.Contains("SKTreeMaker")){
      IDSetup=true;
    return;
    }
  }

  string lqdir =  getenv("LQANALYZER_DIR");

  string username = getenv("USER");

  SetupSelectionMuon(lqdir + "/CATConfig/SelectionConfig/muons.sel");
  SetupSelectionMuon(lqdir + "/CATConfig/SelectionConfig/user_muons.sel");

  SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/electrons.sel");
  SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/user_electrons.sel");
  //if(k_classname.Contains("HNDiElectron"))SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/"+username+"_electrons.sel");
  //if(k_classname.Contains("FakeRateCalculator_El")) SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/"+username+"_electrons.sel");
  //if(k_classname.Contains("ElectronTypes")) SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/"+username+"_electrons.sel");
  SetupSelectionJet(lqdir + "/CATConfig/SelectionConfig/jets.sel");
  SetupSelectionJet(lqdir + "/CATConfig/SelectionConfig/user_jets.sel");

  SetupSelectionFatJet(lqdir + "/CATConfig/SelectionConfig/fatjets.sel");
  SetupSelectionFatJet(lqdir + "/CATConfig/SelectionConfig/user_fatjets.sel");

  IDSetup=true;
  if(k_debugmode){
    for( map<TString,vector<pair<TString,float> > >::iterator it=  selectionIDMapfMuon.begin() ; it !=  selectionIDMapfMuon.end(); it++){
      cout << it->first << endl;
      for (unsigned int i=0 ; i < it->second.size(); i++){
        cout << it->second.at(i).first << " " << it->second.at(i).second << endl;
      }
    }
  }
  Message("SetupSelection DONE", DEBUG);

}

void AnalyzerCore::ConfigureFake(){
  
  /// switch
  if(!k_running_nonprompt) return;
  self_configured=true;
  fake_configured = false;
}


bool AnalyzerCore::ConfigureFakeHists(TString path, std::map<TString, std::pair<std::pair<TString,TString>  ,std::pair<float,TString> > > fake_hists){

  if(!k_running_nonprompt) return false;
  if(fake_configured) return false;

  fake_configured=true;
  fake_path=path;

  bool setupok= m_datadriven_bkg->SetupFake(path, fake_hists);
  if(setupok) return true;
  else{
    cerr << "Trying to configure fakes using ConfigureFakeHists, but fake code is already setup" << endl;
    exit(EXIT_FAILURE);

  }
  return true;
}


void AnalyzerCore::SetupDDBkg(){

  if(k_running_nonprompt || k_running_chargeflip)m_datadriven_bkg = new DataDrivenBackgrounds(self_configured);

  setupDDBkg=true;
  
  /// setup correction class at teh same time
  /// not needed for sktreemaker
  // save time as code does not need to setup and save files
  
  string lqdir =  getenv("LQANALYZER_DIR");

  // List of working points                                                                                                                                                                                                                                                  

  if(k_classname == "SKTreeMaker")  mcdata_correction = new MCDataCorrections();

  if(!k_classname.Contains("SKTreeMaker")){

    mcdata_correction = new MCDataCorrections();
    
    std::vector<TString> vtaggers;
    vtaggers.push_back("CSVv2Moriond17_2017_1_26");
    vtaggers.push_back("cMVAv2Moriond17_2017_1_26");

    std::vector<TString> v_wps;
    v_wps.push_back("Loose");
    v_wps.push_back("Medium");
    v_wps.push_back("Tight");
    MapBTagSF = SetupBTagger(vtaggers,v_wps);
    
    if(1){
      ifstream runlumi((lqdir + "/data/Luminosity/"+getenv("yeartag")+"/lumi_catversion_" + getenv("CATVERSION")+".txt").c_str());
      if(!runlumi) {
	cerr << "Did not find "+lqdir + "/data/Luminosity/"+getenv("yeartag")+"/lumi_catversion_" + getenv("CATVERSION")+".txt'), exiting ..." << endl;
	exit(EXIT_FAILURE);
      }
      string lline;
      int x=1;
      while(getline(runlumi,lline) ){
	std::istringstream is( lline );

	string trigname;
	float trig_lumi;
	int run;
	is >> trigname;
	if(trigname=="###" ) continue;
	if(trigname=="END") break;
	if(trigname=="run" ){
	  is >> run;
	  is >> trig_lumi;
	  if(k_debugmode)
	    cout << "Run number "<< run <<" ; Muon trigger (unprescaled) luminosity  " << trig_lumi << ";" << endl;

	  mapLumi2016[run] = trig_lumi;
	  continue;
	}
	if(trigname=="block" ){
	  is >> run;
	  is >> trig_lumi;
	  if(k_debugmode)cout << "mapLumi[" << run <<" ] = " << trig_lumi << ";" << endl;

	  mapLumiPerBlock2016[run] = trig_lumi;
	  ostringstream ss;
	  ss << x;
	  mapLumiNamePerBlock2016[run]="Lumi"+ss.str();
	  x++;
	  continue;
	}

      }
      runlumi.close();
    }

    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    
    if(k_debugmode)cout << "Reading Luminosity File" << endl;

    SetupLuminosityMap(true);


  }

}

void AnalyzerCore::SetUpEvent(Long64_t entry, float ev_weight) throw( LQError ) {
  
  if(!IDSetup)   SetupID();
  if(!setupDDBkg)SetupDDBkg();
  
//  if(k_running_nonprompt&&fake_configured &&!self_configured){
//    cout << "Setting up fakes(Def)" << endl;
//    m_datadriven_bkg->SetupFake();self_configured=true; }



  Message("In SetUpEvent(Long64_t entry) " , DEBUG);
  m_logger << DEBUG << "This is entry " << entry << LQLogger::endmsg;
  if (!fChain) throw LQError( "Chain is not initialized",  LQError::SkipCycle );     
  
  if(LQinput){

    m_logger << DEBUG << "k_isdata = " << k_isdata << " and isData = " << isData << LQLogger::endmsg;
    if(k_isdata != isData) throw LQError( "!!! Event is confused. It does not know if it is data or MC", LQError::SkipCycle );
  }
  else isData = k_isdata;
  
  if (!(entry % output_interval)) {
    m_logger << INFO <<  "Processing entry " << entry <<  "/" << nentries << LQLogger::endmsg;

  }
  
  snu::KEvent eventinfo = GetEventInfo();
  
  if(k_isdata){
    if(ev_weight!=1.) Message("ERROR in setting weights. This is Data...", INFO);
    MCweight=1.;
    weight = 1.;
  }
  else {
    MCweight = eventinfo.MCWeight(); 
    weight= ev_weight; 
  }
 //
  // creates object that stores all SKTree classes	
  //                                                                                                        
  
  snu::KTrigger triggerinfo = GetTriggerInfo(triggerlist);
    
  std::vector<snu::KJet> skjets= GetAllJets();
  std::vector<snu::KFatJet> skfatjets;
  std::vector<snu::KGenJet> skgenjets;
  if(k_usegenjet) skgenjets=GetAllGenJets();
  if(k_usefatjet) skfatjets= GetAllFatJets();

  /// np == numberof particles you want to store at truth info. 30 is default unless running nocut sktree OR signal
  int np =  AssignnNumberOfTruth();
  
  //// Tmp vect to speed jobs when not used
  
  vector<snu::KPhoton> photons ;
  if(k_usephotons) photons =  GetAllPhotons();
  vector<snu::KTruth> gen;
  if(k_usetruth) gen=GetTruthParticles(np);
  
  vector<snu::KMuon>  muons = GetAllMuons();

  if(!k_classname.Contains("SKTreeMaker")){
    if(k_skim=="FLATCAT"){
      SetCorrectedMomentum(muons, gen);
    }
  }
  LQEvent lqevent(muons, GetAllElectrons(), photons, skjets,skfatjets, skgenjets, gen, triggerinfo,eventinfo);
  

  //  eventbase is master class to use in analysis 
  //
  
  eventbase = new EventBase(lqevent);
  
  /*eventbase->GetElectronSel()->SetIDSMap(selectionIDMapsElectron);
    eventbase->GetElectronSel()->SetIDFMap(selectionIDMapfElectron);*/
  eventbase->GetJetSel()->SetIDSMap(selectionIDMapsJet);
  eventbase->GetJetSel()->SetIDFMap(selectionIDMapfJet);
  eventbase->GetFatJetSel()->SetIDSMap(selectionIDMapsFatJet);
  eventbase->GetFatJetSel()->SetIDFMap(selectionIDMapfFatJet);

  bool setupFullIDinSelection=false;
  /// setting this to true makes the code run slower

  if(k_running_nonprompt || k_running_chargeflip || setupFullIDinSelection){
    //m_datadriven_bkg needs ID maps to get isTight for IDs 
    eventbase->GetElectronSel()->SetIDSMap(selectionIDMapsElectron); 
    eventbase->GetElectronSel()->SetIDFMap(selectionIDMapfElectron);
    eventbase->GetMuonSel()->SetIDSMap(selectionIDMapsMuon);
    eventbase->GetMuonSel()->SetIDFMap(selectionIDMapfMuon);
    m_datadriven_bkg->SetEventBase(eventbase);
  }
  
  if(!k_isdata){
    if(!changed_target_lumi){
      changed_target_lumi=true;
    }
  }
  

  /// Setup correction class
  k_reset_period=true;
  

  if(!k_classname.Contains("SKTreeMaker")){
    mcdata_correction->SetPeriod(GetPeriod());
    mcdata_correction->SetIsData(isData);
  }
  if (k_classname == "SKTreeMaker"){
    mcdata_correction->SetPeriod(GetPeriod());
    mcdata_correction->SetIsData(isData);
    
  }
  
  

}


int AnalyzerCore::VersionStamp(TString cversion){
  
  if(cversion.Contains("v7-4-4")) return 1;
  else if(cversion.Contains("v7-4-5")) return 2;
  else if(cversion.Contains("v7-6-2") || cversion.Contains("v7-6-3") || cversion.Contains("v7-6-4")   ) return 3;
  else if((cversion.Contains("v7-6-5") || cversion.Contains("v7-6-6"))) return 4;
  else if((cversion.Contains("v8-0-1"))) return 5;
  else if((cversion.Contains("v8-0-2"))) return 6;
  else if((cversion.Contains("v8-0-3"))) return 7;
  else if((cversion.Contains("v8-0-4") || cversion.Contains("v8-0-5"))) return 8;
  
  return 5;
 
}



snu::KTruth AnalyzerCore::GetTruthMatchedParticle(snu::KElectron el){
  
  if(el.MCTruthIndex() >  eventbase->GetTruth().size()){
    snu::KTruth tr;
    return tr;
  }
  return   eventbase->GetTruth()[el.MCTruthIndex()];
}


int AnalyzerCore::AssignnNumberOfTruth(){
  int np = 1000;
  if(k_classname.Contains("SKTreeMaker")) np = 1000;
  if(k_classname.Contains("SKTreeMakerDiLep")) np = 0;
  if(k_classname.Contains("SKTreeMakerTriLep")) np = 1000;

  if(k_classname.Contains("SKTreeMaker")){
    if(k_sample_name.Contains("QCD") && !k_sample_name.Contains("mad")) np = 0;
  }

  /// List of signal samples
  /// G.Yu needs to add signal here
  
  if(IsSignal()) np = 1000;
  
  return np;
}



bool AnalyzerCore::IsSignal(){

  if(isData) return false;
  if(k_sample_name.Contains("Majornana")) return true;
  if(k_sample_name.Contains("Tchannel")) return true;
  if(k_sample_name.Contains("HNE")) return true;
  if(k_sample_name.Contains("HNM")) return true;
  if(k_sample_name.Contains("HNDilepton"))  return false;
  if(k_sample_name.Contains("MM")) return true;
  
  return false;
}


void AnalyzerCore::ClassInfo(){
  
  /*if(eventinfo.CatVersion().empty()){ 
    m_logger << INFO << "Catuple version is v7-4-X. Only basic infomation is available." << LQLogger::endmsg;
    
  }
    
  else if(TString(eventinfo.CatVersion()).Contains("v7-6-2")){
    m_logger << INFO <<  "Running on catuples version " << eventinfo.CatVersion() << LQLogger::endmsg;
    
    
    }*/

}

float AnalyzerCore::SumPt( std::vector<snu::KJet> particles){

  float sumpt=0.;
  for(std::vector<snu::KJet>::iterator it = particles.begin(); it != particles.end(); it++){
    sumpt += it->Pt();
  }
  return sumpt;
}


float AnalyzerCore::SumPt( std::vector<snu::KFatJet> particles){

  float sumpt=0.;
  for(std::vector<snu::KFatJet>::iterator it = particles.begin(); it != particles.end(); it++){
    sumpt += it->Pt();
  }
  return sumpt;
}

  

float AnalyzerCore::GetLT(std::vector<snu::KMuon> muons){
  float lt=0.;
  for(unsigned int i = 0; i < muons.size(); i++){
    lt+= muons[i].Pt();
  }

  return lt;
}

float AnalyzerCore::GetLT(std::vector<snu::KElectron> electrons){
  float lt=0.;
  for(unsigned int i = 0; i < electrons.size(); i++){
    lt+= electrons[i].Pt();
  }
  
  return lt;
}


bool AnalyzerCore::IsDiEl(){
  if(isData) return false;
  int iel(0);
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    if(fabs(eventbase->GetTruth().at(ig).PdgId()) == 11) iel++;
  }
  if(iel >1) return true;
  else return false;
}

bool AnalyzerCore::ISCF(snu::KElectron el){
  
  if(el.GetType() == 4) return true;
  if(el.GetType() == 5)return true;
  if(el.GetType() == 6)return true;
  if(el.GetType() == 13&&el.MCMatched())return true;
  if(el.GetType() == 19)return true;
  if(el.GetType() == 20)return true;
  if(el.GetType() == 21)return true;
  return false;
}

bool AnalyzerCore::IsInternalConversion(snu::KMuon mu){

  if(isData) return false;

  bool conv=false;
  std::vector<snu::KTruth> truthColl= eventbase->GetTruth();

  if(GetLeptonType(mu,truthColl )== 4 ||  GetLeptonType(mu,truthColl )==  5) {
    int tr_index= mu.MCTruthIndex();
    while(fabs(eventbase->GetTruth().at(tr_index).PdgId()) == 13 || fabs(eventbase->GetTruth().at(tr_index).PdgId()) == 22){
      tr_index = eventbase->GetTruth().at(tr_index).IndexMother();
    }
    if(fabs(eventbase->GetTruth().at(tr_index).PdgId()) == 23 || fabs(eventbase->GetTruth().at(tr_index).PdgId()) == 24 || fabs(eventbase->GetTruth().at(tr_index).PdgId()) == 15) conv=true;
  }
  
  return conv;

}

bool AnalyzerCore::IsInternalConversion(snu::KElectron el){

  if(isData) return false;
  std::vector<snu::KTruth> truthColl= eventbase->GetTruth();

  bool conv=false;
  if(GetLeptonType(el,truthColl )== 4 ||  GetLeptonType(el,truthColl )==  5) {
    int tr_index= el.MCTruthIndex();
    while(fabs(eventbase->GetTruth().at(tr_index).PdgId()) == 11 || fabs(eventbase->GetTruth().at(tr_index).PdgId()) == 22){
      tr_index = eventbase->GetTruth().at(tr_index).IndexMother();
    }
    if(fabs(eventbase->GetTruth().at(tr_index).PdgId()) == 23 || fabs(eventbase->GetTruth().at(tr_index).PdgId()) == 24 || fabs(eventbase->GetTruth().at(tr_index).PdgId()) == 15) conv=true;
  }

  return conv;

}

bool AnalyzerCore::IsExternalConversion(snu::KElectron el){

  if(isData) return false;

  std::vector<snu::KTruth> truthColl= eventbase->GetTruth();

  bool conv=false;
  if(GetLeptonType(el,truthColl )== -5 ||  GetLeptonType(el,truthColl )==  -6) {
    conv=true;
  }
  if(el.GetType()==40) conv=true;
  return conv;

}



bool AnalyzerCore::TruthMatched(snu::KElectron el, bool keepCF){
  bool pass=false;
  if(!keepCF && ISCF(el)) return false;
  if(keepCF && ISCF(el)) return true;
  if((keepCF && !ISCF(el)) || !keepCF) {
    
    if(el.GetType() ==1)   pass=true; /// Z/W
    if(el.GetType() ==2)   pass=true; /// Z/W
    if(el.GetType() ==3)   pass=true; /// Z/W 
    if(el.GetType() == 11) pass=true; /// Tau
    if(el.GetType() == 14) pass=true; /// Z*
    if(el.GetType() == 15) pass=true; /// W*
    if(el.GetType() == 17) pass=true; /// * CF
    if(el.GetType() == 18) pass=true; /// * CF  
    //if(el.GetType() == 23) pass=true;
    if(el.GetType() == 35) pass=true;
    if(el.GetType() == 40) pass=true;
    if(el.GetType() == 16) pass=true;
  }

  return pass;
}

int AnalyzerCore::IsFakeEvent(vector<snu::KElectron> els ){
  
  if( ((k_sample_name.Contains("TT"))&&SameCharge(els))) return 2;
  if( ((k_sample_name.Contains("DY"))&&SameCharge(els))) return 2;
  if( ((k_sample_name.Contains("WJets"))&&els.size()==2)) return 2;
  if( ((k_sample_name.Contains("qcd"))&&els.size()>0)) return els.size();
  if( ((k_sample_name.Contains("QCD"))&&els.size()>1)) return els.size();
  
  return -1;

}

int  AnalyzerCore::IsFakeEvent(vector<snu::KMuon> mus ){

  if( ((k_sample_name.Contains("TT"))&&SameCharge(mus))) return 2;
  if( ((k_sample_name.Contains("DY"))&&SameCharge(mus))) return 2;
  if( ((k_sample_name.Contains("WJets"))&&mus.size()==2)) return 2;
  if( ((k_sample_name.Contains("qcd"))&&mus.size()>0)) return mus.size();
  if( ((k_sample_name.Contains("QCD"))&&mus.size()>1)) return mus.size();

  return -1;



}


bool AnalyzerCore::NonPrompt(snu::KElectron el){
  
  if(el.GetType() == 7) return true;
  return false;

}
bool AnalyzerCore::NonPrompt(snu::KMuon mu){

  if(mu.GetType() == 2) return true;
  return false;

}


bool AnalyzerCore::AllPrompt(std::vector<snu::KMuon> muons, int method){
  
  if(isData) return true;
  bool allprompt=true;
  std::vector<snu::KTruth> truthColl= eventbase->GetTruth();


  for(unsigned int im = 0; im < muons.size(); im++){
    if(method==0 && !(TruthMatched(muons[im]) && muons[im].MCMatched())) allprompt=false;
    int LepType=GetLeptonType(muons[im],truthColl);
    if(method==1 && (LepType<=0)) allprompt=false;
  }
  return allprompt;
}
  

bool AnalyzerCore::TruthMatched(std::vector<snu::KElectron> el, bool tightdxy, bool allowCF){
  
  bool pass=false;
  
  return pass;
}

bool AnalyzerCore::TruthMatched(snu::KMuon mu){

  bool pass=false;
  
  if(mu.GetType() ==1) pass=true;
  if(mu.GetType() ==6) pass=true;
  if(mu.GetType() ==8) pass=true;
  if(mu.GetType() ==9) pass=true;
  if(mu.GetType() ==10) pass=true; // Virtual photon(m<5) to dimuon. FR > 0.4
  if(mu.GetType() ==12) pass=true;
  if(mu.GetType() ==25) pass=true;  //// HAS CHANGED AFTER SKTREE were remade
  //if(mu.GetType() ==28) pass=true;
  // 23>?? 28?? sshould these be included?
  return pass;
}

vector<int> AnalyzerCore::GetVirtualMassIndex(int mode, int pdgid){
  
  vector<int> indexZ;
  vector<int> indexG;
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    if(indexZ.size()==2&&mode==1) return indexZ;
    if(indexZ.size()>1) continue;
    if(fabs(eventbase->GetTruth().at(ig).PdgId()) == pdgid){
      int index_m=eventbase->GetTruth().at(ig).IndexMother() ;
      if(eventbase->GetTruth().at(ig).GenStatus() ==1  || eventbase->GetTruth().at(ig).GenStatus() ==23 ){
	int daughter=ig;
	while(fabs(eventbase->GetTruth().at(index_m).PdgId()) == pdgid){
	  daughter=index_m;
	  index_m=eventbase->GetTruth().at(index_m).IndexMother();
	}
	if(eventbase->GetTruth().at(index_m).PdgId() == 23 || fabs(eventbase->GetTruth().at(index_m).PdgId()) < 6 ){
	  cout << "daughter = " << daughter << endl;
	  indexZ.push_back(daughter);
	  for(unsigned int ig2=0; ig2 < eventbase->GetTruth().size(); ig2++){
	    if(eventbase->GetTruth().at(ig2).IndexMother() <= 0)continue;
	    if(ig2 == daughter) continue;
	    if(fabs(eventbase->GetTruth().at(ig2).PdgId()) == pdgid){
	      cout << eventbase->GetTruth().at(ig2).IndexMother() << " ind " << index_m << endl;
	      if(eventbase->GetTruth().at(ig2).IndexMother()==index_m)           indexZ.push_back(ig2);
	    }
	  }
	}
      }
    }
  }
  if(mode==1) return indexZ;
  if(indexZ.size()!=2) return indexG;

  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
  
  if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;

    if(fabs(eventbase->GetTruth().at(ig).PdgId()) == pdgid){
      if(indexZ[0] == ig || indexZ[1] == ig ) continue;
      
      if(eventbase->GetTruth().at(ig).GenStatus() ==1){

	int index_m=eventbase->GetTruth().at(ig).IndexMother() ;
	
	while(fabs(eventbase->GetTruth().at(index_m).PdgId()) == pdgid){
	  index_m=eventbase->GetTruth().at(index_m).IndexMother();
	}
	
	for(unsigned int ig2=0; ig2 < eventbase->GetTruth().size(); ig2++){
	  
	  if(eventbase->GetTruth().at(ig2).IndexMother() <= 0)continue;
	  if(eventbase->GetTruth().at(ig2).IndexMother() >= int(eventbase->GetTruth().size()))continue;
	  if(fabs(eventbase->GetTruth().at(ig2).PdgId()) == pdgid){
	    if(indexZ[0] == ig2 || indexZ[1] == ig2 ) continue;
	    if(ig==ig2) continue;
	    if(eventbase->GetTruth().at(ig2).GenStatus() ==1){
	      int index_m2=eventbase->GetTruth().at(ig2).IndexMother() ;
	      while(fabs(eventbase->GetTruth().at(index_m2).PdgId()) == pdgid){
		index_m2=eventbase->GetTruth().at(index_m2).IndexMother();
	      }
	      if(index_m2 == index_m){
		indexG.push_back(ig);
		indexG.push_back(ig2);
		return indexG;
	      }
	    }
	  }
	}
      }
    }
  }
   

}
float AnalyzerCore::GetVirtualMass(int pdg, bool includenu, bool includeph){
  if(isData) return -999.;
  vector<KTruth> es1;
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){

    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    
    if(fabs(eventbase->GetTruth().at(ig).PdgId()) == pdg){
      if(eventbase->GetTruth().at(ig).GenStatus() ==1){
	int index_m=eventbase->GetTruth().at(ig).IndexMother() ;
	while(fabs(eventbase->GetTruth().at(index_m).PdgId()) == pdg){
	  index_m=eventbase->GetTruth().at(index_m).IndexMother();

	}
	if(eventbase->GetTruth().at(index_m).PdgId() == 23 || eventbase->GetTruth().at(index_m).PdgId() ==22){
	  es1.push_back(eventbase->GetTruth().at(ig));
	}
      }
    }
    else   if(includenu){
      if(fabs(eventbase->GetTruth().at(ig).PdgId()) == pdg+1){
	if(eventbase->GetTruth().at(ig).GenStatus() ==1){
	  es1.push_back(eventbase->GetTruth().at(ig));
	}
      }
    }

    else if(includeph){
      if(eventbase->GetTruth().at(ig).GenStatus() ==1){
	es1.push_back(eventbase->GetTruth().at(ig));
      }
    }
  }

  if(!includeph){
    if(!includenu){
      if(es1.size()==2){
	cout << "Mother = " << eventbase->GetTruth().at(es1[0].IndexMother()).PdgId() << endl;

	snu::KParticle ll = es1[0]  + es1[1];
	return ll.M();
      }
    }
    else  if(es1.size()==2){
      snu::KParticle ll = es1[0]  + es1[1];
      return ll.M();
    }
  }
  else{
    snu::KParticle ll;
    for(unsigned int iel=0; iel < es1.size(); iel++){
      ll+= es1[iel];
    }
    return ll.M();
  }
  
  return -999.;

}




float AnalyzerCore::GetVirtualMassConv(int cmindex,int nconvindx){

  if(isData) return -999.;
  cout << "cmindex = " << cmindex << endl;
  vector<KTruth> es1;
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){

    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    
    if(eventbase->GetTruth().at(ig).IndexMother() == eventbase->GetTruth().at(cmindex).IndexMother()){
      if(fabs(eventbase->GetTruth().at(ig).PdgId()) == 11){
	if(eventbase->GetTruth().at(ig).GenStatus() ==1){
	  es1.push_back(eventbase->GetTruth().at(ig));
	}
      }
    }
  }
  
  if(es1.size()==3){
    if(nconvindx==0){
      snu::KParticle ll = es1[0]  + es1[1];
      return ll.M();
    }
    if(nconvindx==1){
      snu::KParticle ll = es1[0]  + es1[2];
      return ll.M();
    }
    if(nconvindx==2){
      snu::KParticle ll = es1[1]  + es1[2];
      return ll.M();
    }
  }
  return -999.;
}



void AnalyzerCore::TruthPrintOut(){
  if(isData) return;
  m_logger << INFO<< "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  cout << "Particle Index |  PdgId  | GenStatus   | Mother PdgId |  Part_Eta | Part_Pt | Part_Phi | Mother Index |   " << endl;



  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){

    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    //if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;                                                                                                                                                                              
    if (eventbase->GetTruth().at(ig).PdgId() == 2212)  cout << ig << " | " << eventbase->GetTruth().at(ig).PdgId() << "  |               |         |        |         |       |         |" << endl;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size())){
      cout << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  ---  |   " << eventbase->GetTruth().at(ig).Eta() << " | " << eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi()	   << " |  --- |" << eventbase->GetTruth().at(ig).ReadStatusFlag(7) <<  endl;

    }
    else{
      cout << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  " << eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()<< " |   " << eventbase->GetTruth().at(ig).Eta() << " | " <<	eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi()<< " |   " << eventbase->GetTruth().at(ig).IndexMother()  << " " << eventbase->GetTruth().at(ig).ReadStatusFlag(7) <<  " " <<  eventbase->GetTruth().at(ig).M() <<endl;
    }
  }

}

void AnalyzerCore::TruthPrintOut(snu::KMuon muon){
  if(isData) return;
  m_logger << INFO<< "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  cout << "Particle Index |  PdgId  | GenStatus   | Mother PdgId |  Part_Eta | Part_Pt | Part_Phi | Mother Index |   " << endl;



  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){

    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    //if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;                                                                                                                                                                              
    if (eventbase->GetTruth().at(ig).PdgId() == 2212)  cout << ig << " | " << eventbase->GetTruth().at(ig).PdgId() << "  |               |         |        |         |       |         |" << endl;
    double dr = sqrt( pow(fabs(  muon.Eta() - eventbase->GetTruth().at(ig).Eta()),2.0) +  pow( fabs(TVector2::Phi_mpi_pi( muon.Phi() - eventbase->GetTruth().at(ig).Phi())),2.0));											
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size())){

      
      if(dr < 0.4)  cout << "MATCHE TO MUON: " << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  ---  |   " << eventbase->GetTruth().at(ig).Eta() << " | " << eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi()	   << " |  --- |" << eventbase->GetTruth().at(ig).ReadStatusFlag(7) <<  endl;
      else cout << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  ---  |   " << eventbase->GetTruth().at(ig).Eta() << " | " << eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi()	   << " |  --- |" << eventbase->GetTruth().at(ig).ReadStatusFlag(7) <<  endl; 
      
    }
    else{
      if(dr < 0.4)  cout << "MATCHE TO MUON: " << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  " << eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()<< " |   " << eventbase->GetTruth().at(ig).Eta() << " | " <<	eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi()<< " |   " << eventbase->GetTruth().at(ig).IndexMother()  << " " << eventbase->GetTruth().at(ig).ReadStatusFlag(7) <<  endl;
      else cout << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  " << eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()<< " |   " << eventbase->GetTruth().at(ig).Eta() << " | " <<	eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi()<< " |   " << eventbase->GetTruth().at(ig).IndexMother()  << " " << eventbase->GetTruth().at(ig).ReadStatusFlag(7) <<  endl;
    }
  }
}

void AnalyzerCore::TruthPrintOut(snu::KElectron electron){
  if(isData) return;
  m_logger << INFO<< "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  cout << "Particle Index |  PdgId  | GenStatus   | Mother PdgId |  Part_Eta | Part_Pt | Part_Phi | Mother Index |   " << endl;



  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){

    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    //if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    if (eventbase->GetTruth().at(ig).PdgId() == 2212)  cout << ig << " | " << eventbase->GetTruth().at(ig).PdgId() << "  |               |         |        |         |       |         |" << endl;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size())){
      cout << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  ---  |   " << eventbase->GetTruth().at(ig).Eta() << " | " << eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi()<< " |  --- |" << eventbase->GetTruth().at(ig).ReadStatusFlag(7) <<  endl;

    }
    else{
      cout << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  " << eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()<< " |   " << eventbase->GetTruth().at(ig).Eta() << " | " << eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi()<< " |   " << eventbase->GetTruth().at(ig).IndexMother()  << " " << eventbase->GetTruth().at(ig).ReadStatusFlag(7) <<  endl; 
    }
  }

}


bool AnalyzerCore::isPrompt(long pdgid) {
  /// mother pdgid
  pdgid = abs(pdgid);
  if (pdgid == 24) return true; // Z
  else if (pdgid == 23) return true; // W
  else if (pdgid == 15) return true; // taus
  else if (pdgid == 90) return true; // N
  else return false;
}

void AnalyzerCore::EndEvent()throw( LQError ){
  
  if(self_configured && !fake_configured) {
    cerr << "Setting up own fake files but no hists given" << endl;
    exit(EXIT_FAILURE);
  }

  delete eventbase;                                                                                                            

}
  

void AnalyzerCore::ListTriggersAvailable(){
  cout << "Set of triggers you can use are: " << endl;
  for(unsigned int i=0; i < eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames().size(); i++){
    cout << eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames().at(i)<< " has prescale " << eventbase->GetTrigger().GetHLTInsideDatasetTriggerPrescales().at(i)<< endl;
  }
  return;
}


bool AnalyzerCore::PassJets(std::vector<snu::KJet> jets, std::vector<snu::KFatJet> fatjets ){

  if(jets.size() < 2 && fatjets.size() == 0) return false;
  
  return true;
}

//################################################################################################
//@@@@@@@@@@@@@@@@@@@  TRIGGER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     

bool AnalyzerCore::PassTrigger(vector<pair<TString,TString> > list){
  
  vector<TString> triglist;
  for(vector<pair<TString,TString>  >::iterator it = list.begin(); it != list.end(); it++){
    if(k_channel.Contains(it->second)){
      for(unsigned int i=0; i <  triglist.size(); i++){
	if(PassTrigger(it->first)) return false;
      }
      if(PassTrigger(it->first)) return true;
      triglist.push_back(it->first);
    }
  }
  return false;
}



bool AnalyzerCore::PassTrigger(TString trig){
  vector<TString> list;
  list.push_back(trig);
  int  prescaler=1.;
  return TriggerSelector(list, eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames(), eventbase->GetTrigger().GetHLTInsideDatasetTriggerDecisions(), eventbase->GetTrigger().GetHLTInsideDatasetTriggerPrescales(), prescaler);
 
 
}


bool AnalyzerCore::PassTriggerOR(vector<TString> list){

  int  prescaler=1.;
  return TriggerSelector(list, eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames(), eventbase->GetTrigger().GetHLTInsideDatasetTriggerDecisions(), eventbase->GetTrigger().GetHLTInsideDatasetTriggerPrescales(), prescaler);
  
}


////###############################################################################################
/// @@@@@@@@@@@@@@@@@@@@@@@@@ MISC   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

TDirectory* AnalyzerCore::GetTemporaryDirectory(void) const
{
  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:                                              
    std::stringstream dirname;
    dirname << "CATAnalzer_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:                                                                       
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }

  return tempDir;

}


std::vector<TLorentzVector> AnalyzerCore::MakeTLorentz( std::vector<snu::KElectron> el){return mcdata_correction->MakeTLorentz(el);}
std::vector<TLorentzVector> AnalyzerCore::MakeTLorentz( std::vector<snu::KMuon> mu){return mcdata_correction->MakeTLorentz(mu);}
std::vector<TLorentzVector> AnalyzerCore::MakeTLorentz( std::vector<snu::KJet> jet){return mcdata_correction->MakeTLorentz(jet);}
std::vector<TLorentzVector> AnalyzerCore::MakeTLorentz( std::vector<snu::KFatJet> jet){return mcdata_correction->MakeTLorentz(jet);}

/// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2  FUNCTIONS USED BY TRILEPTON ANALYSIS: jskim

//==== HN Trilepton stuffs 

void AnalyzerCore::PutNuPz(TLorentzVector *nu, double Pz){
  double Px, Py;
  Px = nu->Px();
  Py = nu->Py();
  nu->SetPxPyPzE(Px, Py, Pz, TMath::Sqrt(Px*Px+Py*Py+Pz*Pz));
}

void AnalyzerCore::PutNuPz(snu::KParticle *nu, double Pz){
  double Px, Py;
  Px = nu->Px();
  Py = nu->Py();
  nu->SetPxPyPzE(Px, Py, Pz, TMath::Sqrt(Px*Px+Py*Py+Pz*Pz));
}

double AnalyzerCore::solveqdeq(double W_mass, TLorentzVector l1l2l3, double MET, double METphi, TString pm){
  TLorentzVector met;
  met.SetPxPyPzE(MET*cos(METphi),
		 MET*sin(METphi),
		 0,
		 MET);

  Double_t d = (W_mass*W_mass)-(l1l2l3.M())*(l1l2l3.M())+2.0*l1l2l3.Px()*met.Px()+2.0*l1l2l3.Py()*met.Py();
  Double_t a = l1l2l3.E()*l1l2l3.E() - l1l2l3.Pz()*l1l2l3.Pz();
  Double_t b = d*l1l2l3.Pz();
  Double_t c = l1l2l3.E()*l1l2l3.E()*met.E()*met.E()-d*d/4.0;
  if(b*b-4*a*c<0){
    return b/(2*a);
  }
  else{
    if(pm=="p") return (b+TMath::Sqrt(b*b-4*a*c))/(2*a);
    else if(pm=="m")  return (b-TMath::Sqrt(b*b-4*a*c))/(2*a);
    else return 0;
  }
}

int AnalyzerCore::find_mlmet_closest_to_W(snu::KParticle lep[], snu::KParticle MET, int n_lep){
  double m_diff[n_lep];
  double m_diff_min = 999999999.;
  int outindex = 0;
  for(int i=0; i<n_lep; i++){
    double dphi = lep[i].DeltaPhi(MET);
    double mt2 = 2.*lep[i].Pt()*MET.Pt()*(1.-TMath::Cos(dphi));
    m_diff[i] = fabs( sqrt(mt2) - 80.385 );
    if( m_diff[i] < m_diff_min ){
      m_diff_min = m_diff[i];
      outindex = i;
    }
  }
  return outindex;
  
}

double AnalyzerCore::MT(TLorentzVector a, TLorentzVector b){
  double dphi = a.DeltaPhi(b);
  return TMath::Sqrt( 2.*a.Pt()*b.Pt()*(1.- TMath::Cos(dphi) ) );
  
}

bool AnalyzerCore::GenMatching(snu::KParticle gen, snu::KParticle reco, double maxDeltaR, double maxPtDiff){
  bool matched = true;
  if(reco.Pt() > 200.) maxPtDiff = 0.1;
  
  if( gen.DeltaR(reco) >= maxDeltaR ) matched = false;
  if( fabs(gen.Pt() - reco.Pt()) / reco.Pt() >= maxPtDiff ) matched = false;
  
  return matched;
  
}

std::vector<snu::KMuon> AnalyzerCore::GetHNTriMuonsByLooseRelIso(double LooseRelIsoMax, bool keepfake){
 
  std::vector<snu::KMuon> muontriLooseColl_raw = GetMuons("MUON_HN_TRI_VLOOSE", keepfake);
  std::vector<snu::KMuon> muontriLooseColl;
  muontriLooseColl.clear();
  for(unsigned int j=0; j<muontriLooseColl_raw.size(); j++){
    snu::KMuon this_muon = muontriLooseColl_raw.at(j);
    if( this_muon.RelIso04() < LooseRelIsoMax ) muontriLooseColl.push_back( this_muon );
  }
  return muontriLooseColl;
  
}

void AnalyzerCore::PrintTruth(){
  std::vector<snu::KTruth> truthColl;
  eventbase->GetTruthSel()->Selection(truthColl);
  
  cout << "=========================================================" << endl;
  cout << "truth size = " << truthColl.size() << endl;
  //cout << "index" << '\t' << "pdgid" << '\t' << "mother" << '\t' << "mother pid" << endl;
  cout<<"Idx"<<'\t'<<"PID"<<'\t'<<"MIdx"<<'\t'<<"MPID"<<'\t'<<"GenSt"<<'\t'<<"pt"<<'\t'<<"eta"<<'\t'<<"phi"<<endl;
  for(int i=2; i<truthColl.size(); i++){
    //cout << i << '\t' << truthColl.at(i).PdgId() << '\t' << truthColl.at(i).IndexMother() << '\t' << truthColl.at( truthColl.at(i).IndexMother() ).PdgId() << endl;
    cout<<i<<'\t'<<truthColl.at(i).PdgId()<<'\t'<<truthColl.at(i).IndexMother()<<'\t'<<truthColl.at(truthColl.at(i).IndexMother()).PdgId()<<'\t'<<truthColl.at(i).GenStatus()<<'\t'<<truthColl.at(i).Pt()<<'\t'<<truthColl.at(i).Eta()<<'\t'<<truthColl.at(i).Phi()<<endl;
  }

}


void AnalyzerCore::Message(TString message, LQMsgType type){
  m_logger <<  type << message << LQLogger::endmsg;
}



////###############################################################################################                                                                           
 /// @@@@@@@@@@@@@@@@@@@@@@@@@ HIST   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                                                             

void AnalyzerCore::MakeCleverHistograms(histtype type, TString clhistname ){
  
  if(type==sighist_e|| type==sighist_ee || type==sighist_eee) {
    map<TString, SignalPlotsEE*>::iterator fit = mapCLhistSigEE.find(clhistname);
    if (fit != mapCLhistSigEE.end()) return;
  }

  if(type==sighist_m|| type==sighist_mm || type==sighist_mmm) {
    map<TString, SignalPlotsMM*>::iterator fit = mapCLhistSigMM.find(clhistname);
    if (fit != mapCLhistSigMM.end()) return;
  }
  
  //// ELECTRON PLOTs                                                                                          
  if(type==elhist) mapCLhistEl[clhistname] = new ElectronPlots(clhistname);
  //// BaseSelection::MUON PLOTs                                                                                              
  if(type==muhist) mapCLhistMu[clhistname] = new MuonPlots(clhistname);
  /// JET PLOTs                                                                                                
  if(type==jethist) mapCLhistJet[clhistname] = new JetPlots(clhistname);
  /// Signal plots                                                                                             
  if(type==sighist_e)  mapCLhistSigEE[clhistname] = new SignalPlotsEE(clhistname, 1);
  if(type==sighist_ee)  mapCLhistSigEE[clhistname] = new SignalPlotsEE(clhistname, 2);
  if(type==sighist_eee)  mapCLhistSigEE[clhistname] = new SignalPlotsEE(clhistname, 3);
  if(type==sighist_eeee)  mapCLhistSigEE[clhistname] = new SignalPlotsEE(clhistname,-1);
  if(type==sighist_m)  mapCLhistSigMM[clhistname] = new SignalPlotsMM(clhistname,1);
  if(type==sighist_mm)  mapCLhistSigMM[clhistname] = new SignalPlotsMM(clhistname,2);
  if(type==sighist_mmm)  mapCLhistSigMM[clhistname] = new SignalPlotsMM(clhistname,3);
  if(type==sighist_mmmm)  mapCLhistSigMM[clhistname] = new SignalPlotsMM(clhistname,-1);
  if(type==sighist_em)  mapCLhistSigEM[clhistname] = new SignalPlotsEM(clhistname);

  if(type==trilephist)  mapCLhistTriLep[clhistname] = new TriLeptonPlots(clhistname);
  if(type==hnpairmm) mapCLhistHNpairMM[clhistname] = new HNpairPlotsMM(clhistname);
  if(type==hntrilephist)  mapCLhistHNTriLep[clhistname] = new HNTriLeptonPlots(clhistname);
    
  return;
}

void AnalyzerCore::MakeHistograms(){
  //// Additional plots to make                                                                                
  maphist.clear();
  maphist2D.clear();
  maphist3D.clear();

    
}

void AnalyzerCore::MakeHistograms(TString hname, int nbins, float xbins[], TString label){
  maphist[hname] =  new TH1D(hname.Data(),hname.Data(),nbins,xbins);
  maphist[hname]->GetXaxis()->SetTitle(label);
}

void AnalyzerCore::MakeHistograms(TString hname, int nbins, float xmin, float xmax, TString label){

  maphist[hname] =  new TH1D(hname.Data(),hname.Data(),nbins,xmin,xmax);
  //maphist[hname]->GetXaxis()->SetTitle("TEST");
  maphist[hname]->GetXaxis()->SetTitle(label);
}


void AnalyzerCore::MakeHistograms2D(TString hname, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, TString label, TString labely) {

  maphist2D[hname] =  new TH2D(hname.Data(),hname.Data(),nbinsx,xmin,xmax, nbinsy,ymin,ymax);
  maphist2D[hname]->GetXaxis()->SetTitle(label);
  maphist2D[hname]->GetYaxis()->SetTitle(labely);
}


void AnalyzerCore::MakeHistograms3D(TString hname, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, int nbinsz, float zmin, float zmax, TString label) {
  
  maphist3D[hname] =  new TH3D(hname.Data(),hname.Data(),nbinsx,xmin,xmax, nbinsy,ymin,ymax, nbinsz, zmin, zmax);
  maphist3D[hname]->GetXaxis()->SetTitle(label);
}

void AnalyzerCore::MakeHistograms2D(TString hname, int nbinsx,  float xbins[], int nbinsy,  float ybins[], TString label, TString labely) {

  maphist2D[hname] =  new TH2D(hname.Data(),hname.Data(),nbinsx , xbins, nbinsy,ybins);
  maphist2D[hname]->GetXaxis()->SetTitle(label);
  maphist2D[hname]->GetYaxis()->SetTitle(labely);
}


void AnalyzerCore::FillHist(TString histname, float value, float w, float xbins[], int nbins , TString label){
  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;

  if(GetHist(histname)) GetHist(histname)->Fill(value, w);
  
  else{
    if (nbins < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms(histname, nbins, xbins, label);
    if(GetHist(histname)) GetHist(histname)->GetXaxis()->SetTitle(label);
    if(GetHist(histname)) GetHist(histname)->Fill(value, w);
  }

}


void AnalyzerCore::FillHistPerLumi(TString histname, float value, float w, float xmin, float xmax,int nbins, int nlumibins){
  
  if(!Is2015Analysis()){
    if(nlumibins==10){
      
      if(!GetHist(histname+"_perlumi")) {
	MakeHistograms(histname+"_perlumi", 9, 0., 9.);
	int nbin=0;

	for(std::map<int,TString>::iterator it = mapLumiNamePerBlock.begin(); it != mapLumiNamePerBlock.end(); it++){
	  nbin++;
	  GetHist(histname+"_perlumi")->GetXaxis()->SetBinLabel(nbin,it->second);
	}
      }
      
      for(map<int,TString>::iterator it = mapLumiNamePerBlock.begin(); it != mapLumiNamePerBlock.end(); it++){
	
        if(!GetHist(histname+"_"+it->second)) {
          MakeHistograms(histname+"_"+it->second, nbins, xmin, xmax);
	}
      } 


      for(map<int,TString>::iterator it = mapLumiNamePerBlock.begin(); it != mapLumiNamePerBlock.end(); it++){
	if(eventbase->GetEvent().RunNumber()  < it->first) {
	  map<int,float>::iterator it2 = mapLumiPerBlock.find(it->first);
	  
	  if(isData){
	    float neww= w /it2->second;
	    if(GetHist(histname+"_perlumi")) GetHist(histname+"_perlumi")->Fill(it->second, neww);
	    if(GetHist(histname+"_"+it->second)) GetHist(histname+"_"+it->second)->Fill(value,neww);
	  }
	  else{
	    //	    float neww = w * (it2->second/TargetLumi);
	    float neww = w /TargetLumi;
	    if(GetHist(histname+"_perlumi")) GetHist(histname+"_perlumi")->Fill(it->second, neww);
	    if(GetHist(histname+"_"+it->second)) GetHist(histname+"_"+it->second)->Fill(value,neww);
	  }
	} 
      }
    }// nbins
  }

}

void AnalyzerCore::FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins , TString label){
  
  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist(histname)) GetHist(histname)->Fill(value, w);  
  else{
    if (nbins < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms(histname, nbins, xmin, xmax, label);
    if(GetHist(histname)) GetHist(histname)->GetXaxis()->SetTitle(label);
    if(GetHist(histname)) GetHist(histname)->Fill(value, w);
  }
  
}

void AnalyzerCore::FillHist(TString histname, float value1, float value2, float w, float xmin, float xmax, int nbinsx, float ymin, float ymax, int nbinsy , TString label, TString labely){

  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist2D(histname)) GetHist2D(histname)->Fill(value1,value2, w);
  else{
    if (nbinsx < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms2D(histname, nbinsx, xmin, xmax,nbinsy, ymin, ymax , label, labely);
    if(GetHist2D(histname)) GetHist2D(histname)->GetXaxis()->SetTitle(label);
    if(GetHist2D(histname)) GetHist2D(histname)->GetYaxis()->SetTitle(labely);
    if(GetHist2D(histname)) GetHist2D(histname)->Fill(value1,value2, w);
  }

}

void AnalyzerCore::FillHist(TString histname, float valuex, float valuey, float w, float xbins[], int nxbins, float ybins[], int nybins , TString label, TString labely){
  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist2D(histname)) GetHist2D(histname)->Fill(valuex,valuey, w);

  else{
    if (nxbins < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms2D(histname, nxbins, xbins, nybins, ybins , label,labely);
    if(GetHist2D(histname)) GetHist2D(histname)->GetXaxis()->SetTitle(label);
    
    if(GetHist2D(histname)) GetHist2D(histname)->Fill(valuex, valuey, w);
  }

}


void AnalyzerCore::FillHist(TString histname, float value1, float value2,float value3, float w, float xmin, float xmax, int nbinsx, float ymin, float ymax, int nbinsy , float zmin, float zmax, int nbinxz, TString label){

  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist3D(histname)) GetHist3D(histname)->Fill(value1,value2,value3, w);
  else{
    if (nbinsx < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms3D(histname, nbinsx, xmin, xmax,nbinsy, ymin, ymax , nbinxz, zmin, zmax,label);
    if(GetHist3D(histname)) GetHist3D(histname)->GetXaxis()->SetTitle(label);
    if(GetHist3D(histname)) GetHist3D(histname)->Fill(value1,value2,value3, w);
  }

}


void AnalyzerCore::FillHist(TString histname, float value, float w , TString label){

  if(GetHist(histname)){
    GetHist(histname)->Fill(value, w);  /// Plots Z peak        
    GetHist(histname)->GetXaxis()->SetTitle(label);
  }
  
  else m_logger << INFO << histname << " was NOT found. Will add the histogram to the hist map on first event." << LQLogger::endmsg;
  
  
  return;
}

void AnalyzerCore::FillCLHist(histtype type, TString hist, vector<snu::KMuon> muons, double w){

  if(type==muhist){
    map<TString, MuonPlots*>::iterator mupit = mapCLhistMu.find(hist);
    if(mupit != mapCLhistMu.end()) mupit->second->Fill(w,muons);
    else m_logger << INFO  << hist << " not found in mapCLhistMu" << LQLogger::endmsg;
  }
  else  m_logger << INFO  << "Type not set to muhist, is this a mistake?" << LQLogger::endmsg;

}


void AnalyzerCore::FillCLHist(histtype type, TString hist, vector<snu::KElectron> electrons,double w){

  if(type==elhist){
    map<TString, ElectronPlots*>::iterator elpit = mapCLhistEl.find(hist);
    if(elpit !=mapCLhistEl.end()) elpit->second->Fill(w,electrons);
    else m_logger << INFO  << hist << " not found in mapCLhistEl" <<LQLogger::endmsg;
  }
  else  m_logger << INFO  << "Type not set to elhist, is this a mistake?" << LQLogger::endmsg;
}

void AnalyzerCore::FillCLHist(histtype type, TString hist, vector<snu::KJet> jets, double w){
  
  if(type==jethist){
    map<TString, JetPlots*>::iterator jetpit = mapCLhistJet.find(hist);
    if(jetpit !=mapCLhistJet.end()) jetpit->second->Fill(w,jets);
    else m_logger << INFO  << hist << " not found in mapCLhistJet" <<LQLogger::endmsg;
  }
  else  m_logger << INFO  <<"Type not set to jethist, is this a mistake?" << LQLogger::endmsg;

}
void AnalyzerCore::FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,double w, int nbjet){
  if(type==hnpairmm){
    map<TString, HNpairPlotsMM*>::iterator HNpairmmit = mapCLhistHNpairMM.find(hist);
    if(HNpairmmit !=mapCLhistHNpairMM.end()) HNpairmmit->second->Fill(ev, muons, electrons, jets, w, nbjet);
    else {
      mapCLhistHNpairMM[hist] = new HNpairPlotsMM(hist);
      HNpairmmit = mapCLhistHNpairMM.find(hist);
      HNpairmmit->second->Fill(ev, muons, electrons, jets, w, nbjet);
    }
  }
}


void AnalyzerCore::FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,double w){
  
  vector<snu::KFatJet> fatjets;
  vector<snu::KJet> alljets;
  FillCLHist(type, hist, ev, muons,electrons ,jets, alljets,fatjets,w);
}

void AnalyzerCore::FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets,vector<snu::KJet> alljets,double w){

  vector<snu::KFatJet> fatjets;
  FillCLHist(type, hist, ev, muons,electrons ,jets, alljets, fatjets,w);
}


void AnalyzerCore::FillCLHist(histtype type, TString hist, snu::KEvent ev,vector<snu::KMuon> muons, vector<snu::KElectron> electrons, vector<snu::KJet> jets, vector<snu::KJet> alljets, vector<snu::KFatJet> fatjets,double w){

  if(type==trilephist){

    map<TString, TriLeptonPlots*>::iterator trilepit = mapCLhistTriLep.find(hist);
    if(trilepit !=mapCLhistTriLep.end()) trilepit->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistTriLep[hist] = new TriLeptonPlots(hist);
      trilepit = mapCLhistTriLep.find(hist);
      trilepit->second->Fill(ev, muons, electrons, jets,w);
    }
  }
  else if(type==hntrilephist){

    map<TString, HNTriLeptonPlots*>::iterator hntrilepit = mapCLhistHNTriLep.find(hist);
    if(hntrilepit !=mapCLhistHNTriLep.end()) hntrilepit->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistHNTriLep[hist] = new HNTriLeptonPlots(hist);
      hntrilepit = mapCLhistHNTriLep.find(hist);
      hntrilepit->second->Fill(ev, muons, electrons, jets,w);
    }
  }
 
  else if(type==sighist_ee){

    map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.find(hist);
    if(sigpit_ee !=mapCLhistSigEE.end()) sigpit_ee->second->Fill(ev, muons, electrons, jets,  alljets,fatjets,w);
    else {
      mapCLhistSigEE[hist] = new SignalPlotsEE(hist,2);
      sigpit_ee = mapCLhistSigEE.find(hist);
      sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    }
  }
  else if(type==sighist_eee){

    map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.find(hist);
    if(sigpit_ee !=mapCLhistSigEE.end()) sigpit_ee->second->Fill(ev, muons, electrons, jets,alljets, fatjets,w);
    else {
      mapCLhistSigEE[hist] = new SignalPlotsEE(hist, 3);
      sigpit_ee = mapCLhistSigEE.find(hist);
      sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets,fatjets,w);
    }
  }
  else if(type==sighist_eeee){

    map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.find(hist);
    if(sigpit_ee !=mapCLhistSigEE.end()) sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    else {
      mapCLhistSigEE[hist] = new SignalPlotsEE(hist,-1);
      sigpit_ee = mapCLhistSigEE.find(hist);
      sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    }
  }
  else if(type==sighist_e){

    map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.find(hist);
    if(sigpit_ee !=mapCLhistSigEE.end()) sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    else {
      mapCLhistSigEE[hist] = new SignalPlotsEE(hist,1);
      sigpit_ee = mapCLhistSigEE.find(hist);
      sigpit_ee->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    }
  }


  else if(type==sighist_m){

    map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.find(hist);
    if(sigpit_mm !=mapCLhistSigMM.end())  sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    else {
      mapCLhistSigMM[hist] = new SignalPlotsMM(hist,1);
      sigpit_mm = mapCLhistSigMM.find(hist);
      sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);

    }
  }
  else if(type==sighist_mm){

    map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.find(hist);
    if(sigpit_mm !=mapCLhistSigMM.end())  sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    else {
      mapCLhistSigMM[hist] = new SignalPlotsMM(hist,2);
      sigpit_mm = mapCLhistSigMM.find(hist);
      sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
      
    }
  }
  else if(type==sighist_mmm){

    map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.find(hist);
    if(sigpit_mm !=mapCLhistSigMM.end())  sigpit_mm->second->Fill(ev, muons, electrons, jets,alljets, fatjets,w);
    else {
      mapCLhistSigMM[hist] = new SignalPlotsMM(hist,3);
      sigpit_mm = mapCLhistSigMM.find(hist);
      sigpit_mm->second->Fill(ev, muons, electrons, jets,alljets, fatjets,w);
      
    }
  }
  else if(type==sighist_mmmm){

    map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.find(hist);
    if(sigpit_mm !=mapCLhistSigMM.end())  sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
    else {
      mapCLhistSigMM[hist] = new SignalPlotsMM(hist,-1);
      sigpit_mm = mapCLhistSigMM.find(hist);
      sigpit_mm->second->Fill(ev, muons, electrons, jets, alljets, fatjets,w);
      
    }
  }

  else if(type==sighist_em){

    map<TString, SignalPlotsEM*>::iterator sigpit_em = mapCLhistSigEM.find(hist);
    if(sigpit_em !=mapCLhistSigEM.end()) sigpit_em->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistSigEM[hist] = new SignalPlotsEM(hist);
      sigpit_em = mapCLhistSigEM.find(hist);
      sigpit_em->second->Fill(ev, muons, electrons, jets,w);
    }
  }
 else  m_logger << INFO  <<"Type not set to sighist, is this a mistake?" << LQLogger::endmsg;


}


void AnalyzerCore::WriteHistograms() throw (LQError){
  // This function is called after the cycle is ran. It wrues all histograms to the output file. This function is not used by user. But by the contrioller code.
  WriteHists();
  WriteCLHists();
  WriteNtp();
}

  
void AnalyzerCore::WriteCLHists(){

  for(map<TString, MuonPlots*>::iterator mupit = mapCLhistMu.begin(); mupit != mapCLhistMu.end(); mupit++){

    Dir = m_outputFile->mkdir(mupit->first);
    m_outputFile->cd( Dir->GetName() );
    mupit->second->Write();
    m_outputFile->cd();
  }

  for(map<TString, ElectronPlots*>::iterator elpit = mapCLhistEl.begin(); elpit != mapCLhistEl.end(); elpit++)\
    {

      Dir = m_outputFile->mkdir(elpit->first);
      m_outputFile->cd( Dir->GetName() );
      elpit->second->Write();
      m_outputFile->cd();
    }

  for(map<TString, JetPlots*>::iterator jetpit = mapCLhistJet.begin(); jetpit != mapCLhistJet.end(); jetpit++)\
    {
      
      Dir = m_outputFile->mkdir(jetpit->first);
      m_outputFile->cd( Dir->GetName() );
      jetpit->second->Write();
      m_outputFile->cd();
    }
  for(map<TString, SignalPlotsEE*>::iterator sigpit_ee = mapCLhistSigEE.begin(); sigpit_ee != mapCLhistSigEE.end(); sigpit_ee++){
    
    Dir = m_outputFile->mkdir(sigpit_ee->first);
    m_outputFile->cd( Dir->GetName() );
    sigpit_ee->second->Write();
    m_outputFile->cd();
  }
  for(map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.begin(); sigpit_mm != mapCLhistSigMM.end(); sigpit_mm++){

    Dir = m_outputFile->mkdir(sigpit_mm->first);
    m_outputFile->cd( Dir->GetName() );
    sigpit_mm->second->Write();
    m_outputFile->cd();
  }
  for(map<TString, SignalPlotsEM*>::iterator sigpit_em = mapCLhistSigEM.begin(); sigpit_em != mapCLhistSigEM.end(); sigpit_em++){

    Dir = m_outputFile->mkdir(sigpit_em->first);
    m_outputFile->cd( Dir->GetName() );
    sigpit_em->second->Write();
    m_outputFile->cd();
  }
  
  for(map<TString, TriLeptonPlots*>::iterator trilepit = mapCLhistTriLep.begin(); trilepit != mapCLhistTriLep.end(); trilepit++){

    Dir = m_outputFile->mkdir(trilepit->first);
    m_outputFile->cd( Dir->GetName() );
    trilepit->second->Write();
    m_outputFile->cd();
  }

  for(map<TString, HNpairPlotsMM*>::iterator HNpairmmit = mapCLhistHNpairMM.begin(); HNpairmmit != mapCLhistHNpairMM.end(); HNpairmmit++){
    
    Dir = m_outputFile->mkdir(HNpairmmit->first);
    m_outputFile->cd( Dir->GetName() );
    HNpairmmit->second->Write();
    m_outputFile->cd();
  }

  for(map<TString, HNTriLeptonPlots*>::iterator hntrilepit = mapCLhistHNTriLep.begin(); hntrilepit != mapCLhistHNTriLep.end(); hntrilepit++){

    //==== (jskim)I don't need director!
    //Dir = m_outputFile->mkdir(hntrilepit->first);
    //m_outputFile->cd( Dir->GetName() );
    hntrilepit->second->Write();
    m_outputFile->cd();
  }



  return;
}

void AnalyzerCore::WriteHists(){

  /// Open Output rootfile
  m_outputFile->cd();

  for(map<TString, TH1*>::iterator mapit = maphist.begin(); mapit != maphist.end(); mapit++){
    
    if(mapit->first.Contains("cutflow")){
      mapit->second->Write();
    }
    else if(mapit->first.Contains("Basic")){
      if(!m_outputFile->GetDirectory("Basic")){
	Dir = m_outputFile->mkdir("Basic");
	m_outputFile->cd(Dir->GetName());
        mapit -> second -> Write();
        m_outputFile -> cd();
      }
      else{
        m_outputFile->cd("Basic");
        mapit -> second -> Write();
        m_outputFile -> cd();
      }
    }
    else{
//      TDirectory *dir = m_outputFile->GetDirectory("Hists");
//   
//      if (dir) {
//	m_outputFile->cd("Hists");
//	mapit->second->Write();
//	m_outputFile->cd();
//      }
//      else{
//	Dir = m_outputFile->mkdir("Hists");
//	m_outputFile->cd( Dir->GetName() );
        mapit->second->Write();
//	m_outputFile->cd();
//      }
    }
  }
  for(map<TString, TH2*>::iterator mapit = maphist2D.begin(); mapit != maphist2D.end(); mapit++){
    
    TDirectory *dir = m_outputFile->GetDirectory("Hists2D");

    if (dir) {
      m_outputFile->cd("Hists2D");
      mapit->second->Write();
      m_outputFile->cd();
    }
    else{
      Dir = m_outputFile->mkdir("Hists2D");
      m_outputFile->cd( Dir->GetName() );
      mapit->second->Write();
      m_outputFile->cd();
    }
  }
  
  for(map<TString, TH3*>::iterator mapit = maphist3D.begin(); mapit != maphist3D.end(); mapit++){
    TDirectory *dir = m_outputFile->GetDirectory("Hists3D");

    if (dir) {
      m_outputFile->cd("Hists3D");
      mapit->second->Write();
      m_outputFile->cd();
    }
    else{
      Dir = m_outputFile->mkdir("Hists3D");
      m_outputFile->cd( Dir->GetName() );
      mapit->second->Write();
      m_outputFile->cd();
    }

  }

  //==== HN Gen Matching
  if((TString(getenv("USER")) == "jskim" || TString(getenv("USER")) =="shjeon")){
    m_HNgenmatch->WriteHNGenHists();
  }
  return;
}

TH1* AnalyzerCore::GetHist(TString hname){

  TH1* h = NULL;
  std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit->second;
  else m_logger << DEBUG  << hname << " was not found in map" << LQLogger::endmsg;

  return h;
}



TH2* AnalyzerCore::GetHist2D(TString hname){

  TH2* h = NULL;
  std::map<TString, TH2*>::iterator mapit = maphist2D.find(hname);
  if(mapit != maphist2D.end()) return mapit->second;
  else m_logger << DEBUG  << hname << " was not found in map" << LQLogger::endmsg;

  return h;
}


TH3* AnalyzerCore::GetHist3D(TString hname){

  TH3* h = NULL;
  std::map<TString, TH3*>::iterator mapit = maphist3D.find(hname);
  if(mapit != maphist3D.end()) return mapit->second;
  else m_logger << DEBUG  << hname << " was not found in map" << LQLogger::endmsg;

  return h;
}




void AnalyzerCore::FillCutFlow(TString cut, float weight){

  bool cut_exists=false;
  m_logger << DEBUG  << "FillCutFlow" <<  LQLogger::endmsg;
  for(unsigned int i=0; i < cutflow_list.size();  i++){
    if (cut == cutflow_list.at(i)) cut_exists=true;
  }

  if (!cut_exists){
    m_logger << DEBUG  << "FillCutFlow: " << cut <<  LQLogger::endmsg;

    n_cutflowcuts=n_cutflowcuts+1;
    if(GetHist("cutflow")) {
      m_logger << DEBUG  << "FillCutFlow: mod " << cut <<  LQLogger::endmsg;

      map<TString, TH1*>::iterator it = maphist.find("cutflow");
      vector<float> counters;
      for(int j = 0; j  < n_cutflowcuts ; j++){
	counters.push_back(it->second->GetBinContent(1+j));
      }
      delete it->second;
      
      AnalyzerCore::MakeHistograms("cutflow", n_cutflowcuts,0.,float(n_cutflowcuts));
      for(unsigned int i=0;i < cutflow_list.size();  i++){
	GetHist("cutflow")->GetXaxis()->SetBinLabel(i+1,cutflow_list.at(i));
	GetHist("cutflow")->Fill(cutflow_list.at(i), weight);
      }
      GetHist("cutflow")->GetXaxis()->SetBinLabel(n_cutflowcuts,cut);
      cutflow_list.push_back(cut);
      
    }
    else{
      m_logger << DEBUG  << "FillCutFlow:Fill " << cut <<  LQLogger::endmsg;
      AnalyzerCore::MakeHistograms("cutflow", n_cutflowcuts,0.,float(n_cutflowcuts));
      GetHist("cutflow")->GetXaxis()->SetBinLabel(1,cut);
      GetHist("cutflow")->Fill(cut,weight);
      cutflow_list.push_back(cut);
    }
  }
  else {
    if(GetHist("cutflow")) {
      GetHist("cutflow")->Fill(cut,weight);
    }
  }
}

////###############################################################################################
/// @@@@@@@@@@@@@@@@@@@@@@@@@ JETPILEUP   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



////############################################################################################### 
/// @@@@@@@@@@@@@@@@@@@@@@@@@ CUTS   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

bool AnalyzerCore::PassMETFilter(){

  bool pass (true);
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#MiniAOD_8011_ICHEP_dataset

  if (!eventbase->GetEvent().PassTightHalo2016Filter()) {
    pass = false;    m_logger << DEBUG << "Event Fails PassTightHalo2016Filter " << LQLogger::endmsg;
  }
  if (!eventbase->GetEvent().PassHBHENoiseFilter()) {
    pass = false;
    m_logger << DEBUG << "Event Fails PassHBHENoiseFilter " << LQLogger::endmsg;
  }
  if (!eventbase->GetEvent().PassHBHENoiseIsoFilter()) {
    pass = false;
    m_logger << DEBUG << "Event Fails PassHBHENoiseIsoFilter " << LQLogger::endmsg;
  }
  if(!eventbase->GetEvent().PassEcalDeadCellTriggerPrimitiveFilter()) {
    pass = false;
    m_logger << DEBUG << "Event Fails PassEcalDeadCellTriggerPrimitiveFilter" << LQLogger::endmsg;
  }
  
  if(!eventbase->GetEvent().PassBadChargedCandidateFilter()) {
    pass = false;
    m_logger << DEBUG << "Event Fails PassBadChargedCandidateFilterr" << LQLogger::endmsg;
  }

  if(!eventbase->GetEvent().PassBadPFMuonFilter()) {
    pass = false;
    m_logger << DEBUG << "Event Fails  PassBadPFMuonFilter" << LQLogger::endmsg;
  }

  if (isData){
    if(!eventbase->GetEvent().PassBadEESupercrystalFilter()) {
      pass = false;
      m_logger << DEBUG << "Event FailsPassBadEESupercrystalFilter" << LQLogger::endmsg;
    }
  }
  

  return pass;
}



bool AnalyzerCore::Zcandidate(std::vector<snu::KMuon> muons, float interval, bool require_os){
  
  if(muons.size()!=2) return false;
  if(require_os&&SameCharge(muons)) return false;
  
  KParticle Z = muons.at(0) + muons.at(1);
  if(fabs(Z.M() - 90.) <  interval) return true;
  else return false;
  
}
  
bool AnalyzerCore::SameCharge(std::vector<snu::KMuon> muons){
  
  if(muons.size()!=2) return false;
  if(muons.at(0).Charge() == muons.at(1).Charge()) return true;
  return false;
}


bool AnalyzerCore::Zcandidate(std::vector<snu::KElectron> electrons, float interval, bool require_os){

  if(electrons.size()!=2) return false;
  if(require_os&&SameCharge(electrons)) return false;

  KParticle Z = electrons.at(0) + electrons.at(1);
  if(fabs(Z.M() - 90.) <  interval) return true;
  else return false;

}

bool AnalyzerCore::SameCharge(std::vector<snu::KElectron> electrons, bool runningcf){

  if(electrons.size()!=2) return false;


  if(!runningcf){
    if(electrons.at(0).Charge() == electrons.at(1).Charge()) return true;
  }
  else     if(electrons.at(0).Charge() != electrons.at(1).Charge()) return true;

  return false;
}


bool AnalyzerCore::OppositeCharge(std::vector<snu::KElectron> electrons, std::vector<snu::KMuon> muons){
  
  if(electrons.size() != 1) return false;
  if(muons.size() != 1) return false;
  
  if(electrons[0].Charge() == muons[0].Charge()) return false;
  return true;
}

bool AnalyzerCore::OppositeCharge(std::vector<snu::KElectron> electrons, bool runningcf){

  if(electrons.size()!=2) return false;

  if(!runningcf){
    if(electrons.at(0).Charge() != electrons.at(1).Charge()) return true;
  }
  else     if(electrons.at(0).Charge() == electrons.at(1).Charge()) return true;

  return false;
}

std::vector<snu::KElectron> AnalyzerCore::ShiftElectronEnergy(std::vector<snu::KElectron> beforeshift, TString el_ID, bool applyshift){

  if(el_ID != "ELECTRON_HN_TIGHTv4") return beforeshift;
  if(!applyshift) return beforeshift;

  std::vector<snu::KElectron> aftershift;
  double shiftrate = -999.;
  if(beforeshift.size() == 1) shiftrate = (1-0.024);
  if(beforeshift.size() == 2) shiftrate = (1-0.015);
  if(beforeshift.size() > 2) shiftrate = (-999.);
   

   for(unsigned int i=0; i < beforeshift.size(); i++){
     beforeshift.at(i).SetPtEtaPhiM(beforeshift.at(i).Pt()*shiftrate, beforeshift.at(i).Eta(), beforeshift.at(i).Phi(), 0.511e-3);
     aftershift.push_back(beforeshift.at(i));
   }
   return aftershift;
 }

float AnalyzerCore::GetCFweight(int syst, std::vector<snu::KElectron> electrons, bool apply_sf, TString el_ID){

  if(el_ID != "ELECTRON_HN_TIGHTv4") return 0.;
  if(electrons.size() > 2) return 0.;

  std::vector<snu::KElectron> lep;
  for(int i=0; i<electrons.size(); i++){
    lep.push_back(electrons.at(i));
  }

  if(lep.size()==2){
    if(lep.at(0).Charge() == lep.at(1).Charge()) return 0.;
  }

  std::vector<double> CFrate, CFweight, sf;
  for(int i=0; i<lep.size(); i++){
    CFrate.push_back(GetCFRates(lep.at(i).Pt(), lep.at(i).SCEta(), el_ID));
    CFweight.push_back( (CFrate.at(i)/(1-CFrate.at(i))) );
  }

  for(int i=0; i<lep.size(); i++){
    if(apply_sf){
      if(fabs(lep.at(i).SCEta()) < 1.4442){
        sf.push_back(0.693589 + (syst*0.693589*0.11));
      }
      else{
        sf.push_back(0.684761 + (syst*0.684761*0.08));
      }
    }
    else{
      sf.push_back( 1 );
    }
  }

  double cfweight = 0.;
  for(int i=0; i<lep.size(); i++){
    cfweight += (sf.at(i)) * (CFweight.at(i));
  }
  return cfweight;
}


float AnalyzerCore::GetCFRates(double el_pt, double el_eta, TString el_ID){
  if(el_ID != "ELECTRON_HN_TIGHTv4") return 0.;

  el_eta = fabs(el_eta);
  if(el_eta > 1.4442 && el_eta < 1.556) return 0.;

  double invPt = 1./el_pt;
  double a = 999., b= 999.;
  if(el_eta < 0.9){
    if(invPt< 0.023){a=(-0.001381); b=(4.334e-05);}
    else{a=(0.001010); b=(-1.146e-05);}
  }
  else if(el_eta < 1.4442){
    if(invPt< 0.015){a=(-0.04296); b=(0.0008670);}
    else if(invPt< 0.023){a=(-0.01529); b=(0.0004522);}
    else{a=(-0.001546); b=(0.0001272);}
  }
  else{
    if(invPt< 0.012){a=(-0.4238); b=(0.006366);}
    else if(invPt< 0.020){a=(-0.1040); b=(0.002550);}
    else{a=(-0.01603); b=(0.0007672);}
  }

  double rate = (a)*invPt + (b);
  if(rate < 0) rate = 0.;
  return rate;

}


int AnalyzerCore::NBJet(std::vector<snu::KJet> jets,  KJet::Tagger tag, KJet::WORKING_POINT wp, int period){

  int nbjet=0;

  if(period == 0) {
    Message("period not set in AnalyzerCore::NBJet. Will assign mcperiod for you but this may not give correct behaviour", WARNING);
    period=GetPeriod();
  }

  TString btag_key_lf("") , btag_key_hf("");
  TString wp_string="";
  if(wp == snu::KJet::Loose)wp_string = "Loose";
  if(wp == snu::KJet::Medium)wp_string = "Medium";
  if(wp == snu::KJet::Tight)wp_string = "Tight";

  TString tag_string="";

  if(tag== snu::KJet::CSVv2) tag_string ="CSVv2Moriond17_2017_1_26";
  
  if(tag== snu::KJet::cMVAv2) tag_string ="cMVAv2Moriond17_2017_1_26";
   
  btag_key_lf = tag_string+"_"+wp_string+"_lf";
  btag_key_hf = tag_string+"_"+wp_string+"_hf";
  std::map<TString,BTagSFUtil*>::iterator it_lf = MapBTagSF.find(btag_key_lf);
  std::map<TString,BTagSFUtil*>::iterator it_hf = MapBTagSF.find(btag_key_hf);
  
  if(it_lf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}
  if(it_hf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}

  /// systematics allowed are +-1 and +-3 for HN analysis 
  if ( tag == snu::KJet::JETPROB) return -999;
  for(unsigned int ij=0; ij <jets.size(); ij++){

    if(IsBTagged(jets.at(ij), tag, wp)) nbjet++;
    continue;

    bool isBtag=false;
    if (isData) {

      if (it_lf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  -999999, jets.at(ij).Pt(), jets.at(ij).Eta(),period))
	isBtag=true;
    }
    else if (jets.at(ij).HadronFlavour() > 1){
      if (it_hf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  jets.at(ij).HadronFlavour(),jets.at(ij).Pt(), jets.at(ij).Eta(),period))
        isBtag=true;
    }
    else{
      if (it_lf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  jets.at(ij).HadronFlavour(),jets.at(ij).Pt(), jets.at(ij).Eta(),period))
	isBtag=true;
    }
    
    if(isBtag )nbjet++;
  }
  return nbjet;
}


bool AnalyzerCore::IsBTagged(snu::KJet jet,  KJet::Tagger tag, KJet::WORKING_POINT wp, int mcperiod, int syst){

  if(mcperiod == 0) {
    Message("mcperiod not set in AnalyzerCore::IsBTagged. Will assign mcperiod for you but this may not give correct behaviour", WARNING);      
    mcperiod=GetPeriod();
  }

  TString btag_key_lf("") , btag_key_hf("");
  TString wp_string="";
  if(wp == snu::KJet::Loose)wp_string = "Loose";
  if(wp == snu::KJet::Medium)wp_string = "Medium";
  if(wp == snu::KJet::Tight)wp_string = "Tight";

  TString tag_string="";
  if(tag== snu::KJet::CSVv2) tag_string ="CSVv2Moriond17_2017_1_26";
  if(tag== snu::KJet::cMVAv2) tag_string ="cMVAv2Moriond17_2017_1_26";



  /// Data applied no correction. So only mcperiod is set

  btag_key_lf = tag_string+"_"+wp_string+"_lf";
  btag_key_hf = tag_string+"_"+wp_string+"_hf";

  if(syst==0){

  }
  //==== Heavy (Eff) Up
  else if(syst==1){
    btag_key_hf += "_systup";
  }
  //==== Heavy (Eff) Down
  else if(syst==-1){
    btag_key_hf += "_systdown";
  }
  //==== Light (Miss) Up
  else if(syst==3){
    btag_key_lf += "_systup";
  }
  //==== Light (Miss) Down
  else if(syst==-3){
    btag_key_lf += "_systdown";
  }
  else{
    // wrong syst?
  }

  std::map<TString,BTagSFUtil*>::iterator it_lf = MapBTagSF.find(btag_key_lf);
  std::map<TString,BTagSFUtil*>::iterator it_hf = MapBTagSF.find(btag_key_hf);

  if(it_lf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}
  if(it_hf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}

  /// systematics allowed are +-1 and +-3 for HN analysis
  if ( tag == snu::KJet::JETPROB) return -999;
  
  bool isBtag=false;
  if (isData) {
    
    if (it_lf->second->IsTagged(jet.BJetTaggerValue(tag),  -999999, jet.Pt(), jet.Eta(), mcperiod))
      isBtag=true;
  }
    else if (jet.HadronFlavour() > 1){
      if (it_hf->second->IsTagged(jet.BJetTaggerValue(tag),  jet.HadronFlavour(),jet.Pt(), jet.Eta(),mcperiod))
        isBtag=true;
    }
    else{
      if (it_lf->second->IsTagged(jet.BJetTaggerValue(tag),  jet.HadronFlavour(),jet.Pt(), jet.Eta(),mcperiod))
        isBtag=true;
    }
  
  return isBtag;
}


float AnalyzerCore::BTagScaleFactor_1a(std::vector<snu::KJet> jetColl, KJet::Tagger tag, KJet::WORKING_POINT wp, int mcperiod, TString Option){

  //BTag SF from 1a method.
  //This is coded for H+->WA analysis. I'm fine with anybody else using this function, but be aware that HN analyses decided to use 2a method.
  //And I currently have no plan to use multiple WP. so I just coded to work only for single WP regime.

  if(isData) return 1.;

  if(mcperiod == 0) {
    Message("FYI : mcperiod not set in AnalyzerCore::BTagScaleFactor_1a: meaning auto-set", DEBUG);
    mcperiod=GetPeriod();
  }

  TString btag_key_lf(""), btag_key_hf("");
  TString wp_string="";
  if(wp == snu::KJet::Loose)  wp_string = "Loose";
  if(wp == snu::KJet::Medium) wp_string = "Medium";
  if(wp == snu::KJet::Tight)  wp_string = "Tight";

  TString tag_string="";
  if(tag== snu::KJet::CSVv2)  tag_string ="CSVv2Moriond17_2017_1_26";
  if(tag== snu::KJet::cMVAv2) tag_string ="cMVAv2Moriond17_2017_1_26";

  TString Str_SystDir_L="", Str_SystDir_BC="";
  if(Option.Contains("Syst")){
    if(Option.Contains("Up")){
      if     (Option.Contains("LTag"))  Str_SystDir_L ="_systup";
      else if(Option.Contains("BCTag")) Str_SystDir_BC="_systup";
    }
    else if(Option.Contains("Down")){
      if     (Option.Contains("LTag"))  Str_SystDir_L ="_systdown";
      else if(Option.Contains("BCTag")) Str_SystDir_BC="_systdown";
    }
  }

  btag_key_lf = tag_string+"_"+wp_string+"_lf"+Str_SystDir_L;
  btag_key_hf = tag_string+"_"+wp_string+"_hf"+Str_SystDir_BC;
  std::map<TString,BTagSFUtil*>::iterator it_lf = MapBTagSF.find(btag_key_lf);
  std::map<TString,BTagSFUtil*>::iterator it_hf = MapBTagSF.find(btag_key_hf);

  if(it_lf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}
  if(it_hf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}

  if ( tag == snu::KJet::JETPROB) return -999;

  float BTagSF=1.;
  for(unsigned int i=0; i<jetColl.size(); i++){
    if(jetColl.at(i).IsBTagged(tag, wp)){
      if(jetColl.at(i).HadronFlavour()==5 || jetColl.at(i).HadronFlavour()==4){
        BTagSF *= it_hf->second->GetJetSF(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta(),mcperiod);
      }
      else{
        BTagSF *= it_lf->second->GetJetSF(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta(),mcperiod);
      }
    }
    else{
      float SFj=1., Effj=1.;
      if(jetColl.at(i).HadronFlavour()==5 || jetColl.at(i).HadronFlavour()==4){
        SFj  = it_hf->second->GetJetSF(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta(),mcperiod);
        Effj = it_hf->second->JetTagEfficiency(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta());
      }
      else{
        SFj  = it_lf->second->GetJetSF(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta(),mcperiod);
        Effj = it_lf->second->JetTagEfficiency(jetColl.at(i).HadronFlavour(), jetColl.at(i).Pt(), jetColl.at(i).Eta());
      }

      if( (1.-Effj)==0. ) return 0.;
      BTagSF *= (1.-SFj*Effj)/(1.-Effj);
    }
  }

  return BTagSF;

}


double AnalyzerCore::MuonDYMassCorrection(std::vector<snu::KMuon> mu, double w){
  
  if(mu.size()< 2) return 0.;
  snu::KParticle Z = mu.at(0) + mu.at(1);
  
  double factor (1.);
  if(Z.M() > 90.){
    factor = 8.37401e-01 + 1.61277e-03*Z.M();
  }
  return w*factor;
}


bool AnalyzerCore::IsTight(snu::KMuon muon){
  /// ADD TIGHT BaseSelection::MUON REQUIREMENT
  float reliso= muon.RelIso04();

  if(( reliso >= 0.1)) return false;
  if(fabs(muon.dXY()) >= 0.05) return false; 
  return true;
}

bool AnalyzerCore::PassID(snu::KMuon mu, TString muid){


  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsMuon.find(muid);
  std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfMuon.find(muid);
  if(it== selectionIDMapsMuon.end()){
    cerr << "Muon ID ["+muid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {

    return eventbase->GetMuonSel()->PassUserID(muid,mu, it->second,fit->second);
  }
  return true;
}


bool AnalyzerCore::PassID(snu::KElectron el, TString elid){
  

  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsElectron.find(elid);
  std::map<TString, vector<pair<TString,float> > >::iterator fit = selectionIDMapfElectron.find(elid);
  if(it== selectionIDMapsElectron.end()){
    cerr << "Electron ID ["+elid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {

    return eventbase->GetElectronSel()->PassUserID(elid,el, it->second, fit->second);
  }
  return true;
}

bool AnalyzerCore::IsCF(snu::KElectron el){
  vector<snu::KTruth> truth =  eventbase->GetTruth();
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    if(fabs(eventbase->GetTruth().at(ig).PdgId()) == 11){
      if(fabs(eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()) == 23 ||
	 fabs(eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()) == 24){
	if(eventbase->GetTruth().at(ig).PdgId() * el.Charge() > 0 ) return true;
	else return false;
      }
    }
  }
  return false;
}

vector<snu::KElectron> AnalyzerCore::GetTruePrompt(vector<snu::KElectron> electrons, bool keep_chargeflip, bool keepfake){
  
  if(electrons.size() == 0)
    return electrons;
  

  vector<snu::KElectron> prompt_electrons;
  for(unsigned int i = 0; i < electrons.size(); i++){

    if(!k_isdata){
      if(k_running_taudecays){
	
	if(electrons.at(i).MCFromTau())  prompt_electrons.push_back(electrons.at(i));
      }
      else{
	
	bool ismatched = TruthMatched(electrons.at(i),  keep_chargeflip);                                                                                                                            
	if(keepfake&&keep_chargeflip) prompt_electrons.push_back(electrons.at(i));
	else if(keep_chargeflip&& ismatched) prompt_electrons.push_back(electrons.at(i));
	else if(keepfake&&! MCIsCF(electrons.at(i))) prompt_electrons.push_back(electrons.at(i)); 
	else if(ismatched && !MCIsCF(electrons.at(i))) prompt_electrons.push_back(electrons.at(i));
      }
    }// Data
    else prompt_electrons.push_back(electrons.at(i));
  }/// loop

  return prompt_electrons;


}
bool  AnalyzerCore::MCIsCF(snu::KElectron el){
  
  if(el.GetType() == 4) return true;
  if(el.GetType() == 5) return true;
  if(el.GetType() == 6) return true;
  if(el.GetType() == 13) return true;
  if(el.GetType() == 19) return true;
  if(el.GetType() == 20) return true;
  if(el.GetType() == 21) return true;
  
  return false;
}
vector<snu::KMuon> AnalyzerCore::GetTruePrompt(vector<snu::KMuon> muons, bool keepfake){
  if(muons.size()==0)return muons;

  vector<snu::KMuon> prompt_muons;

  for(unsigned int i = 0; i < muons.size(); i++){
    if(!k_isdata){
      
      if(keepfake) prompt_muons.push_back(muons.at(i));
      else if(TruthMatched(muons.at(i))) prompt_muons.push_back(muons.at(i));
    }
    // Data
    else prompt_muons.push_back(muons.at(i));
  }/// loop
  return prompt_muons;
  
}


void AnalyzerCore::CorrectMuonMomentum(vector<snu::KMuon>& k_muons){
  
  mcdata_correction->CorrectMuonMomentum(k_muons,eventbase->GetTruth());
  Message("END CorrectMuonMomentum",DEBUG);    
  return;
    
}

void AnalyzerCore::SetCorrectedMomentum(vector<snu::KMuon>& k_muons, vector<snu::KTruth> truth){

  for(std::vector<snu::KMuon>::iterator it = k_muons.begin(); it != k_muons.end(); it++){

    if(it->RochPt() < 0.){
      if(it->IsPF() && (it->IsGlobal()==1 || it->IsTracker() == 1)&& it->Pt() > 5. && fabs(it->Eta()) < 2.5){
	it->SetRochPt(mcdata_correction->GetCorrectedMuonMomentum(*it, truth));
      }
      else it->SetRochPt(it->Pt());
    }
  }

}


void AnalyzerCore::SetCorrectedMomentum(vector<snu::KMuon>& k_muons){
  
  for(std::vector<snu::KMuon>::iterator it = k_muons.begin(); it != k_muons.end(); it++){
    //if(k_classname=="SKTreeMaker" && (it->RochPt() >0.)) exit(EXIT_FAILURE);
    if(k_classname!="SKTreeMaker" && k_classname.Contains("SKTreeMaker")){
      
    if(it->RochPt() < 0.) {
	cerr << "Roch Pt wrongly set in dilep skim" << endl;
	exit(EXIT_FAILURE);
      }
    }
    
    if(it->RochPt() < 0.){
      /// If not loose muon then it can crash (this is safe unless using FLATCAT)
      if(it->IsPF() && (it->IsGlobal()==1 || it->IsTracker() == 1)&& it->Pt() > 5. && fabs(it->Eta()) < 2.5){
	it->SetRochPt(mcdata_correction->GetCorrectedMuonMomentum(*it, eventbase->GetTruth()));
      }
      else it->SetRochPt(it->Pt());
    }    
  }
  
}

void AnalyzerCore::MakeNtp(TString hname, TString myvar){

  mapntp[hname] =  new TNtupleD(hname.Data(),hname.Data(),myvar.Data());
}


void AnalyzerCore::FillNtp(TString hname, Double_t myinput[]){

  if (GetNtp(hname)) GetNtp(hname)->Fill(myinput);
  else m_logger << INFO << hname << " was NOT found. Check you ntp. " << LQLogger::endmsg;

  return;
}

// //
void AnalyzerCore::WriteNtp(){

  /// Open Output rootfile
  m_outputFile->cd();

  for(map<TString, TNtupleD*>::iterator mapit = mapntp.begin(); mapit != mapntp.end(); mapit++){
    mapit->second->Write();
  }

  return;
}
// //



// //
TNtupleD* AnalyzerCore::GetNtp(TString hname){

  TNtupleD* n = NULL;
  std::map<TString, TNtupleD*>::iterator mapit = mapntp.find(hname);
  if (mapit != mapntp.end()) return mapit->second;
  else m_logger << INFO << hname << " was not found in map" << LQLogger::endmsg;

  return n;
}

std::vector<snu::KMuon> AnalyzerCore::sort_muons_ptorder(std::vector<snu::KMuon> muons){

  std::vector<snu::KMuon> outmuon;
  while(outmuon.size() != muons.size()){
    double this_maxpt = 0.;
    int index(0);
    for(unsigned int i=0; i<muons.size(); i++){
      bool isthisused = std::find( outmuon.begin(), outmuon.end(), muons.at(i) ) != outmuon.end();
      if(isthisused) continue;
      if( muons.at(i).Pt() > this_maxpt ){
        index = i;
        this_maxpt = muons.at(i).Pt();
      }
    }
    outmuon.push_back( muons.at(index) );
  }
  return outmuon;
 

}

//------------------------------------------//
//         H+->WA Shared Tools              //
//------------------------------------------//

//--Gen-Matching Tools----------------------------------------------------------------------//
int AnalyzerCore::GenMatchedIdx(snu::KElectron El, std::vector<snu::KTruth>& truthColl){
  //Find Matched Index within dR01; if ambiguous closest dR one chosen(Resolution way better than dPtRel)
  //Seed from RecoLepton

  int MatchedIdx=-1;
  float dR=999., dRmax=0.1;
  
  for(int i=2; i<truthColl.size(); i++){
    if(truthColl.at(i).IndexMother()<0 )  continue;
    if(truthColl.at(i).GenStatus()!=1)    continue;
    if(fabs(truthColl.at(i).PdgId())!=11) continue;
    if(truthColl.at(i).DeltaR(El)>dRmax)  continue;

    if(truthColl.at(i).DeltaR(El)<dR){ dR=truthColl.at(i).DeltaR(El); MatchedIdx=i; }
  }

  return MatchedIdx;
}


int AnalyzerCore::GenMatchedIdx(snu::KMuon Mu, std::vector<snu::KTruth>& truthColl){
  //Find Matched Index within dR01; if ambiguous closest dR one chosen(Resolution way better than dPtRel)
  //Seed from RecoLepton

  int MatchedIdx=-1;
  float dR=999., dRmax=0.1;
  
  for(int i=2; i<truthColl.size(); i++){
    if(truthColl.at(i).IndexMother()<0 )  continue;
    if(truthColl.at(i).GenStatus()!=1)    continue;
    if(fabs(truthColl.at(i).PdgId())!=13) continue;
    if(truthColl.at(i).DeltaR(Mu)>dRmax)  continue;

    if(truthColl.at(i).DeltaR(Mu)<dR){ dR=truthColl.at(i).DeltaR(Mu); MatchedIdx=i; }
  }

  return MatchedIdx;
}


int AnalyzerCore::GetNearPhotonIdx(snu::KElectron Ele, std::vector<snu::KTruth>& TruthColl, TString Option){//1)

  int NearPhotonIdx=-1;
  bool OnlyHardPhoton=false;
    if(Option=="Hard") OnlyHardPhoton=true;
  float PTthreshold=10.;
  float dRmax=0.2;//2)
  float dRmin=999.;
  for(int i=2; i<TruthColl.size(); i++){
    if( TruthColl.at(i).IndexMother()<0   ) continue;
    if( !(TruthColl.at(i).PdgId()==22 && (TruthColl.at(i).GenStatus()==1 || TruthColl.at(i).GenStatus()==23)) ) continue;
    if( TruthColl.at(i).Pt()<PTthreshold  ) continue;
    if( !(Ele.Pt()/TruthColl.at(i).Pt()>0.8 && Ele.Pt()/TruthColl.at(i).Pt()<1.2) ) continue;//3)
    if( Ele.DeltaR(TruthColl.at(i))>dRmax ) continue;

    if( TruthColl.at(i).GenStatus()==23 && !IsFinalPhotonSt23(TruthColl) ) continue;//4)
    if( Ele.DeltaR(TruthColl.at(i))<dRmin ){ dRmin=Ele.DeltaR(TruthColl.at(i)); NearPhotonIdx=i; }
  }


  return NearPhotonIdx;
//footnote
//1) When checked with ZG sample with CBPOGT ElePt>25, Hardscattered photons, External conversion is only meaningful for electrons.
//   Electron External Conversion~3000 Muon external conversion. This is in agreement with theoretical calculation that xsec~M^{-2} in asymmetric limit.
//   ref. Arxiv:1110.1368v1
//2) When checked with ZG sample with CBPOGT ElePt>25GeV, HardScatter Photons, Matching Eff~99.7% even when tested with dR01 cone.
//   But just used conservative cone size.
//3) When checked with ZG sample with CBPOGT ElePt>25GeV, HardScatter Photons, Matching Eff~99.5%
//   External Conversion object's momentum is roughly symmetric with mother photon's momentum(similar fraction lower/upper than PT(G))
//4) In some cases hard scattered photon(GenSt23) has no daughter like final state particles.
}


int AnalyzerCore::GetNearPhotonIdx(snu::KMuon Mu, std::vector<snu::KTruth>& TruthColl, TString Option){//1)

  int NearPhotonIdx=-1;
  bool OnlyHardPhoton=false;
    if(Option=="Hard") OnlyHardPhoton=true;
  float PTthreshold=10.;
  float dRmax=0.2;//2)
  float dRmin=999.;
  for(int i=2; i<TruthColl.size(); i++){
    if( TruthColl.at(i).IndexMother()<0   ) continue;
    if( !(TruthColl.at(i).PdgId()==22 && (TruthColl.at(i).GenStatus()==1 || TruthColl.at(i).GenStatus()==23)) ) continue;
    if( TruthColl.at(i).Pt()<PTthreshold  ) continue;
    if( !(Mu.Pt()/TruthColl.at(i).Pt()>0.8 && Mu.Pt()/TruthColl.at(i).Pt()<1.2) ) continue;//3)
    if( Mu.DeltaR(TruthColl.at(i))>dRmax ) continue;

    if( TruthColl.at(i).GenStatus()==23 && !IsFinalPhotonSt23(TruthColl) ) continue;//4)
    if( Mu.DeltaR(TruthColl.at(i))<dRmin ){ dRmin=Mu.DeltaR(TruthColl.at(i)); NearPhotonIdx=i; }
  }

  return NearPhotonIdx;
//footnote
//1) When checked with ZG sample with CBPOGT ElePt>25, Hardscattered photons, External conversion is only meaningful for electrons.
//   Electron External Conversion~3000 Muon external conversion. This is in agreement with theoretical calculation that xsec~M^{-2} in asymmetric limit.
//   ref. Arxiv:1110.1368v1
//2) When checked with ZG sample with CBPOGT ElePt>25GeV, HardScatter Photons, Matching Eff~99.7% even when tested with dR01 cone.
//   But just used conservative cone size.
//3) When checked with ZG sample with CBPOGT ElePt>25GeV, HardScatter Photons, Matching Eff~99.5%
//   External Conversion object's momentum is roughly symmetric with mother photon's momentum(similar fraction lower/upper than PT(G))
//4) In some cases hard scattered photon(GenSt23) has no daughter like final state particles.
}


int AnalyzerCore::FirstNonSelfMotherIdx(int TruthIdx, std::vector<snu::KTruth>& TruthColl){

  if(TruthIdx<2) return -1;

  int pid=TruthColl.at(TruthIdx).PdgId(), midx=TruthIdx;
  while(TruthColl.at(midx).PdgId()==pid){
    midx=TruthColl.at(midx).IndexMother();  
    if(midx<0) break;
  }

  return midx;
}


int AnalyzerCore::LastSelfMotherIdx(int TruthIdx,std::vector<snu::KTruth>& TruthColl){

  if(TruthIdx<2) return TruthIdx;

  int pid=TruthColl.at(TruthIdx).PdgId(), midx=TruthIdx, currentidx=TruthIdx;
  while(TruthColl.at(midx).PdgId()==pid){
    currentidx=midx;
    midx=TruthColl.at(midx).IndexMother();  
    if(midx<0) break;
  }

  return currentidx;
}


bool AnalyzerCore::HasHadronicAncestor(int TruthIdx, std::vector<snu::KTruth>& TruthColl){
  //Returns true  if 1)has hadron mother, 2)has quark mother(!top) 3)Incident protons
  //        false if 1)is hardscattered truth, 2)EW/H/BSM/t daughter, 3)not above, 4)invalid input(e.g. unmatched case)
  
  if(TruthIdx<0) return false;
  if(TruthIdx<2) return true;

  bool HasPartonHadronAncestor=false;
  int  midx=TruthIdx, fmid=fabs(TruthColl.at(midx).PdgId()), MSt_orig=-1;
  int  St_orig=TruthColl.at(LastSelfMotherIdx(TruthIdx, TruthColl)).GenStatus();
  if( St_orig>20 && St_orig<30) return false;

  while( midx>=2 ){
    midx=FirstNonSelfMotherIdx(midx,TruthColl);
    MSt_orig=TruthColl.at(LastSelfMotherIdx(midx,TruthColl)).GenStatus();
    fmid=fabs(TruthColl.at(midx).PdgId());
    if(  fmid==23 || fmid==24 || fmid==25 || fmid==6 || fmid==36 || fmid==32 ){ HasPartonHadronAncestor=false; break; }
    if( (fmid==11 || fmid==13 || fmid==15 || fmid==22) && (MSt_orig>20 && MSt_orig<30)){ HasPartonHadronAncestor=false; break; }
    if( fmid>50 ) { HasPartonHadronAncestor=true; break; }
    if( (fmid>=1 && fmid<=5) || fmid==21 ){ HasPartonHadronAncestor=true; break; }
  }

  return HasPartonHadronAncestor;
}


bool AnalyzerCore::IsFinalPhotonSt23(std::vector<snu::KTruth> TruthColl){
//In Some XG proc events, it seems there is status 23 photon, yet no status 1 photon and no other genparticle is daughter of this photon.
//This is to check whether this is the case for the event.
//And this is designed only for 1 hard photon case as W+G or Z+G or TT+G

  bool IsFinalGammaStatus23 = false;
  bool HasStatus23Photon    = false;
  for(int i=2; i<TruthColl.size(); i++){
    int fpid  = fabs(TruthColl.at(i).PdgId());
    int GenSt = TruthColl.at(i).GenStatus();
    int MPID_direct= TruthColl.at(TruthColl.at(i).IndexMother()).PdgId();
    if( !((fpid!=22 && MPID_direct==22) || (fpid==22 && (GenSt==23||GenSt==1))) ) continue;

    int LastSelfIdx  = LastSelfMotherIdx(i,TruthColl);
    int LastSelfSt   = TruthColl.at(LastSelfIdx).GenStatus();
    int MotherIdx    = FirstNonSelfMotherIdx(i,TruthColl);
    int LastSelfMIdx=-1, MStatus_orig=-1;
    if(MotherIdx!=-1){
      LastSelfMIdx = LastSelfMotherIdx(MotherIdx,TruthColl);
      MStatus_orig = TruthColl.at(LastSelfMIdx).GenStatus();
    }

    if(fpid==22){
      if(GenSt==23) {HasStatus23Photon=true; IsFinalGammaStatus23=true;}
      else if(GenSt==1 && LastSelfSt==23) {IsFinalGammaStatus23=false; break;}//a)
    }
    else if( MPID_direct==22 && MStatus_orig==23 ){ IsFinalGammaStatus23=false; break;}//b)
  }

  if(!HasStatus23Photon) return false;
  
  return IsFinalGammaStatus23;

//**footnotes
//a) The status 1 photon is end of the history of status 23 photon.
//b) Some particle is daughter of status 23 photon.
}


int AnalyzerCore::GetLeptonType(int TruthIdx, std::vector<snu::KTruth>& TruthColl, TString Option){
//Type : 1:EWPrompt  /  2:Signal Daughter /  3:EWtau daughter / 4:Internal Conversion daughter from t/EWV/EWlep(Implicit,Explicit) / 5:Internal Conversion daughter from HardScatterPhoton
//      -1:Unmatched & not EW Conversion candidate / -2:Hadron daughter / -3:Daughter of tau from hadron or parton / -4:Internal conversion daughter(implicit,explicit) having hadronic origin / -5:External conversion candidate(Hard scattered photon) / -6:External conversion from t/EWV/EWlep
//      (-4:Daughter of Non-hard scattered photon & has parton or hadron ancestor OR implicit Conv from quark)
//       0:Error / >0: Non-fake: Non-hadronic origin / <0 : Fakes: Hadronic origin or external conversion


  //Only consider Status 1 lepton
  if(TruthIdx<2) return 0;
  if(TruthColl.at(TruthIdx).GenStatus()!=1) return 0;
  if( !(fabs(TruthColl.at(TruthIdx).PdgId())==11 || fabs(TruthColl.at(TruthIdx).PdgId())==13) ) return 0;

  int LeptonType=0;
  int LastSelfIdx     = LastSelfMotherIdx(TruthIdx,TruthColl);
  int MotherIdx       = FirstNonSelfMotherIdx(TruthIdx,TruthColl);
  int LastSelfMIdx    = LastSelfMotherIdx(MotherIdx,TruthColl);
  int GrMotherIdx     = FirstNonSelfMotherIdx(MotherIdx,TruthColl);
  int LastSelfGrMIdx  = LastSelfMotherIdx(GrMotherIdx,TruthColl);

  int MPID=0, GrMPID=0;
  int Status_orig=0, MStatus_orig=0, MStatus_last=0, GrMStatus_orig=0, GrMStatus_last=0;
  bool HadronicOrigin = false;
    if(    TruthIdx!=-1   ){ Status_orig    = TruthColl.at(LastSelfIdx).GenStatus();
                             HadronicOrigin = HasHadronicAncestor(TruthIdx, TruthColl);
                           }                           
    if(   MotherIdx!=-1   ){ MPID         = TruthColl.at(MotherIdx).PdgId();
                             MStatus_orig = TruthColl.at(LastSelfMIdx).GenStatus();
                             MStatus_last = TruthColl.at(MotherIdx).GenStatus();
                           }
    if(  GrMotherIdx!=-1  ){ GrMPID         = TruthColl.at(GrMotherIdx).PdgId();
                             GrMStatus_orig = TruthColl.at(LastSelfGrMIdx).GenStatus();
                             GrMStatus_last = TruthColl.at(GrMotherIdx).GenStatus();
                           }
 
  if     ( TruthIdx==-1 )                                       LeptonType= 0;
  else if( fabs(MPID)==23 || fabs(MPID)==24 || fabs(MPID)==25 ) LeptonType= 1;
  else if( fabs(MPID)==36 || fabs(MPID)==32 )                   LeptonType= 2;
  else if( Status_orig>20 && Status_orig<30 )                   LeptonType= 1;//1)
  else if( fabs(MPID)>50 )                                      LeptonType=-2;
  else if( fabs(MPID)==15 && MStatus_last==2 ){
           if     ( fabs(GrMPID)==23 || fabs(GrMPID)==24 || fabs(GrMPID)==25 ) LeptonType= 3;
           else if( MStatus_orig>20  && MStatus_orig<30  )                     LeptonType= 3;//1)
           else if( HadronicOrigin )                                           LeptonType=-3;//2-a)
           else if( fabs(GrMPID)==22  && GrMStatus_orig>20 && GrMStatus_orig<30 )                     LeptonType= 5;//2-b)
           else if( fabs(GrMPID)==22 )                                                                LeptonType= 4;//2-c)
           else if( (fabs(GrMPID)==11 || fabs(GrMPID)==13 || fabs(GrMPID)==15) && GrMStatus_last!=2 ) LeptonType= 4;//2-d)
           else                                                                                       LeptonType= 0;
         }
  else if( fabs(MPID)==22 ){
           if( MStatus_orig>20 && MStatus_orig<30 )                            LeptonType= 5;//3-a)
           else if( HadronicOrigin )                                           LeptonType=-4;//3-b)
           else if( fabs(GrMPID)==24 || fabs(GrMPID)==23 || fabs(GrMPID)==6  ) LeptonType= 4;//3-c)
           else if( fabs(GrMPID)==11 || fabs(GrMPID)==13 || fabs(GrMPID)==15 ) LeptonType= 4;//3-d)
           else                                                                LeptonType= 0;
         }
  else if( (fabs(MPID)==11 || fabs(MPID)==13 || fabs(MPID)==15) && MStatus_last!=2 && !HadronicOrigin ) LeptonType= 4;//4-a)
  else if( ((fabs(MPID)>=1 && fabs(MPID)<=5) || fabs(MPID)==21) && MStatus_last!=2 )                    LeptonType=-4;//4-b)
  else if( fabs(MPID)==6 ) LeptonType=4;//4-c)
  else LeptonType=0;


  return LeptonType;

//**footnote
//These are based on observation in DY,ZG,TT sample(DY,ZG:amcnlo+pythia, TT:powheg+pythia) for other PS generator, convention may differ.
//1) In amcnlo generator, output of ME level generation does not have specific guage field mother. e.g. u u~ > l+ l- -> fabs(MID)=1
//   This perhaps due to multiple field can interplay in production, and apparently it is not possible to distinguish them in any logic.
//   e.g. think about previous example. you cannot say whether this is from gamma or Z or H...
//   But in PS procedure, corrections on ME proc is done sometimes. In that case it seems mother is set Z for OS ll prod. W for lnu prod.
//   e.g. If pythia applies ISR process on input u u~, than it should affect momentum of all the consequent processes, or in case of lnu, W radiating gamma can be added.
//   You may think lnu case is obvious, but it may not like in case of pp > lllnu(You never know which one is from W>lnu and Z>ll)
//2-a) e.g. a)Had > ta+X, ta>l+2nu b) q>ta+X in jet fragmentation (ta is not hardscattered, since it is already considered prev. step)
//2-b) e.g. gamma>ta(+)+ta(-)+X, ta>lnu (St=2)
//2-c) e.g. " " " " " " " " " " " " " " " " " ", but soft gamma case. this is not observed in test sample but put here just in case.
//          (Non hadronic origin since such case already counted before, gamma should be from non-hadronic source)
//2-d) e.g. l>tata..+l.. , ta>l+2nu (Implicit tau conv. from non-hadronic lepton and decay) In implicit conv. GenStatus!=2
//3-a) e.g. hard gamma>ll
//3-b) e.g. a) Had>gamma+X, gamma>ll+X (in PS+Had stage intermediate process is omitted you see just Had>Nphoton+Mhadrons+..)
//          b) q>gamma+q, gamms>ll+X in jet fragmentation or radiations of tops.
//          c) gluon>Ngamma+Mhadrons in jet fragmentation (Actually observed in samples)
//3-c) e.g. W+>W+ gamma, or t>t+gamma, gamma>ll+X, not yet observed in test sample but possible (upto radiation is observed so far)
//3-d) e.g. ta>ta+gamma, gamma>ll+X, tau not from hadron(e.g. pp>tata)
//4-a) e.g. EW lep l, l>lll... just implicit conversion. 
//4-b) e.g. q or g> Nlepton +MHadrons... in parton shower history
//4-c) e.g. t>t+ll.. implicit conversion
}


int AnalyzerCore::GetLeptonType(snu::KElectron El, std::vector<snu::KTruth>& TruthColl, TString Option){
//Type : 1:EWPrompt  /  2:Signal Daughter /  3:EWtau daughter / 4:Internal Conversion daughter from t/EWV/EWlep(Implicit,Explicit) / 5:Internal Conversion daughter from HardScatterPhoton
//      -1:Unmatched & not EW Conversion candidate / -2:Hadron daughter / -3:Daughter of tau from hadron or parton / -4:Internal conversion daughter(implicit,explicit) having hadronic origin / -5:External conversion candidate(Hard scattered photon) / -6:External conversion from t/EWV/EWlep
//      (-4:Daughter of Non-hard scattered photon & has parton or hadron ancestor OR implicit Conv from quark)
//       0:Error / >0: Non-fake: Non-hadronic origin / <0 : Fakes: Hadronic origin or external conversion

  int LeptonType=0;
  int MatchedTruthIdx = GenMatchedIdx(El,TruthColl);
  int LastSelfIdx     = LastSelfMotherIdx(MatchedTruthIdx,TruthColl);
  int MotherIdx       = FirstNonSelfMotherIdx(MatchedTruthIdx,TruthColl);
  int LastSelfMIdx    = LastSelfMotherIdx(MotherIdx,TruthColl);
  int GrMotherIdx     = FirstNonSelfMotherIdx(MotherIdx,TruthColl);
  int LastSelfGrMIdx  = LastSelfMotherIdx(GrMotherIdx,TruthColl);

  int NearPhotonType=0, NearPhotonIdx=-1;
  int MPID=0, GrMPID=0;
  int Status_orig=0, MStatus_orig=0, MStatus_last=0, GrMStatus_orig=0, GrMStatus_last=0;
  bool HadronicOrigin = false;
    if(MatchedTruthIdx!=-1){ Status_orig    = TruthColl.at(LastSelfIdx).GenStatus();
                             HadronicOrigin = HasHadronicAncestor(MatchedTruthIdx, TruthColl);
                           }
    else                   { NearPhotonIdx  = GetNearPhotonIdx(El, TruthColl);
                             NearPhotonType = GetPhotonType(NearPhotonIdx, TruthColl);
                           }
    if(   MotherIdx!=-1   ){ MPID         = TruthColl.at(MotherIdx).PdgId();
                             MStatus_orig = TruthColl.at(LastSelfMIdx).GenStatus();
                             MStatus_last = TruthColl.at(MotherIdx).GenStatus();
                           }
    if(  GrMotherIdx!=-1  ){ GrMPID         = TruthColl.at(GrMotherIdx).PdgId();
                             GrMStatus_orig = TruthColl.at(LastSelfGrMIdx).GenStatus();
                             GrMStatus_last = TruthColl.at(GrMotherIdx).GenStatus();
                           }

  if     ( fabs(MPID)==23 || fabs(MPID)==24 || fabs(MPID)==25 ) LeptonType= 1;
  else if( fabs(MPID)==36 || fabs(MPID)==32 )                   LeptonType= 2;
  else if( Status_orig>20 && Status_orig<30 )                   LeptonType= 1;//1)
  else if( fabs(MPID)>50 )                                      LeptonType=-2;
  else if( fabs(MPID)==15 && MStatus_last==2 ){
           if     ( fabs(GrMPID)==23 || fabs(GrMPID)==24 || fabs(GrMPID)==25 ) LeptonType= 3;
           else if( MStatus_orig>20  && MStatus_orig<30  )                     LeptonType= 3;//1)
           else if( HadronicOrigin )                                           LeptonType=-3;//2-a)
           else if( fabs(GrMPID)==22  && GrMStatus_orig>20 && GrMStatus_orig<30 )                     LeptonType= 5;//2-b)
           else if( fabs(GrMPID)==22 )                                                                LeptonType= 4;//2-c)
           else if( (fabs(GrMPID)==11 || fabs(GrMPID)==13 || fabs(GrMPID)==15) && GrMStatus_last!=2 ) LeptonType= 4;//2-d)
           else                                                                                       LeptonType= 0;
         }
  else if( fabs(MPID)==22 ){
           if( MStatus_orig>20 && MStatus_orig<30 )                            LeptonType= 5;//3-a)
           else if( HadronicOrigin )                                           LeptonType=-4;//3-b)
           else if( fabs(GrMPID)==24 || fabs(GrMPID)==23 || fabs(GrMPID)==6  ) LeptonType= 4;//3-c)
           else if( fabs(GrMPID)==11 || fabs(GrMPID)==13 || fabs(GrMPID)==15 ) LeptonType= 4;//3-d)
           else                                                                LeptonType= 0;
         }
  else if( MatchedTruthIdx==-1 ){
           if     ( NearPhotonType<=0 ) LeptonType=-1;
           else if( NearPhotonType==1 ) LeptonType=-5;
           else if( NearPhotonType==2 ) LeptonType=-6;
         }
  else if( (fabs(MPID)==11 || fabs(MPID)==13 || fabs(MPID)==15) && MStatus_last!=2 && !HadronicOrigin ) LeptonType= 4;//4-a)
  else if( ((fabs(MPID)>=1 && fabs(MPID)<=5) || fabs(MPID)==21) && MStatus_last!=2 )                    LeptonType=-4;//4-b)
  else if( fabs(MPID)==6 ) LeptonType=4;//4-c)
  else LeptonType=0;

 
  return LeptonType;

//**footnote
//These are based on observation in DY,ZG,TT sample(DY,ZG:amcnlo+pythia, TT:powheg+pythia) for other PS generator, convention may differ.
//1) In amcnlo generator, output of ME level generation does not have specific guage field mother. e.g. u u~ > l+ l- -> fabs(MID)=1
//   This perhaps due to multiple field can interplay in production, and apparently it is not possible to distinguish them in any logic.
//   e.g. think about previous example. you cannot say whether this is from gamma or Z or H...
//   But in PS procedure, corrections on ME proc is done sometimes. In that case it seems mother is set Z for OS ll prod. W for lnu prod.
//   e.g. If pythia applies ISR process on input u u~, than it should affect momentum of all the consequent processes, or in case of lnu, W radiating gamma can be added.
//   You may think lnu case is obvious, but it may not like in case of pp > lllnu(You never know which one is from W>lnu and Z>ll)
//2-a) e.g. a)Had > ta+X, ta>l+2nu b) q>ta+X in jet fragmentation (ta is not hardscattered, since it is already considered prev. step)
//2-b) e.g. gamma>ta(+)+ta(-)+X, ta>lnu (St=2)
//2-c) e.g. " " " " " " " " " " " " " " " " " ", but soft gamma case. this is not observed in test sample but put here just in case.
//          (Non hadronic origin since such case already counted before, gamma should be from non-hadronic source)
//2-d) e.g. l>tata..+l.. , ta>l+2nu (Implicit tau conv. from non-hadronic lepton and decay) In implicit conv. GenStatus!=2
//3-a) e.g. hard gamma>ll
//3-b) e.g. a) Had>gamma+X, gamma>ll+X (in PS+Had stage intermediate process is omitted you see just Had>Nphoton+Mhadrons+..)
//          b) q>gamma+q, gamms>ll+X in jet fragmentation or radiations of tops.
//          c) gluon>Ngamma+Mhadrons in jet fragmentation (Actually observed in samples)
//3-c) e.g. W+>W+ gamma, or t>t+gamma, gamma>ll+X, not yet observed in test sample but possible (upto radiation is observed so far)
//3-d) e.g. ta>ta+gamma, gamma>ll+X, tau not from hadron(e.g. pp>tata)
//4-a) e.g. EW lep l, l>lll... just implicit conversion. 
//4-b) e.g. q or g> Nlepton +MHadrons... in parton shower history
//4-c) e.g. t>t+ll.. implicit conversion
}


int AnalyzerCore::GetLeptonType(snu::KMuon Mu, std::vector<snu::KTruth>& TruthColl, TString Option){
//Type : 1:EWPrompt  /  2:Signal Daughter /  3:EWtau daughter / 4:Internal Conversion daughter from t/EWV/EWlep(Implicit,Explicit) / 5:Internal Conversion daughter from HardScatterPhoton
//      -1:Unmatched & not EW Conversion candidate / -2:Hadron daughter / -3:Daughter of tau from hadron or parton / -4:Internal conversion daughter(implicit,explicit) having hadronic origin / -5:External conversion candidate(Hard scattered photon) / -6:External conversion from t/EWV/EWlep
//      (-4:Daughter of Non-hard scattered photon & has parton or hadron ancestor OR implicit Conv from quark)
//       0:Error / >0: Non-fake: Non-hadronic origin / <0 : Fakes: Hadronic origin or external conversion

  int LeptonType=0;
  int MatchedTruthIdx = GenMatchedIdx(Mu,TruthColl);
  int LastSelfIdx     = LastSelfMotherIdx(MatchedTruthIdx,TruthColl);
  int MotherIdx       = FirstNonSelfMotherIdx(MatchedTruthIdx,TruthColl);
  int LastSelfMIdx    = LastSelfMotherIdx(MotherIdx,TruthColl);
  int GrMotherIdx     = FirstNonSelfMotherIdx(MotherIdx,TruthColl);
  int LastSelfGrMIdx  = LastSelfMotherIdx(GrMotherIdx,TruthColl);

  int NearPhotonType=0, NearPhotonIdx=-1;
  int MPID=0, GrMPID=0;
  int Status_orig=0, MStatus_orig=0, MStatus_last=0, GrMStatus_orig=0, GrMStatus_last=0;
  bool HadronicOrigin = false;
    if(MatchedTruthIdx!=-1){ Status_orig    = TruthColl.at(LastSelfIdx).GenStatus();
                             HadronicOrigin = HasHadronicAncestor(MatchedTruthIdx, TruthColl);
                           }
    else                   { NearPhotonIdx  = GetNearPhotonIdx(Mu, TruthColl);
                             NearPhotonType = GetPhotonType(NearPhotonIdx, TruthColl);
                           }
    if(   MotherIdx!=-1   ){ MPID         = TruthColl.at(MotherIdx).PdgId();
                             MStatus_orig = TruthColl.at(LastSelfMIdx).GenStatus();
                             MStatus_last = TruthColl.at(MotherIdx).GenStatus();
                           }
    if(  GrMotherIdx!=-1  ){ GrMPID         = TruthColl.at(GrMotherIdx).PdgId();
                             GrMStatus_orig = TruthColl.at(LastSelfGrMIdx).GenStatus();
                             GrMStatus_last = TruthColl.at(GrMotherIdx).GenStatus();
                           }

  if     ( fabs(MPID)==23 || fabs(MPID)==24 || fabs(MPID)==25 ) LeptonType= 1;
  else if( fabs(MPID)==36 || fabs(MPID)==32 )                   LeptonType= 2;
  else if( Status_orig>20 && Status_orig<30 )                   LeptonType= 1;//1)
  else if( fabs(MPID)>50 )                                      LeptonType=-2;
  else if( fabs(MPID)==15 && MStatus_last==2 ){
           if     ( fabs(GrMPID)==23 || fabs(GrMPID)==24 || fabs(GrMPID)==25 ) LeptonType= 3;
           else if( MStatus_orig>20  && MStatus_orig<30  )                     LeptonType= 3;//1)
           else if( HadronicOrigin )                                           LeptonType=-3;//2-a)
           else if( fabs(GrMPID)==22  && GrMStatus_orig>20 && GrMStatus_orig<30 )                     LeptonType= 5;//2-b)
           else if( fabs(GrMPID)==22 )                                                                LeptonType= 4;//2-c)
           else if( (fabs(GrMPID)==11 || fabs(GrMPID)==13 || fabs(GrMPID)==15) && GrMStatus_last!=2 ) LeptonType= 4;//2-d)
           else                                                                                       LeptonType= 0;
         }
  else if( fabs(MPID)==22 ){
           if( MStatus_orig>20 && MStatus_orig<30 )                            LeptonType= 5;//3-a)
           else if( HadronicOrigin )                                           LeptonType=-4;//3-b)
           else if( fabs(GrMPID)==24 || fabs(GrMPID)==23 || fabs(GrMPID)==6  ) LeptonType= 4;//3-c)
           else if( fabs(GrMPID)==11 || fabs(GrMPID)==13 || fabs(GrMPID)==15 ) LeptonType= 4;//3-d)
           else                                                                LeptonType= 0;
         }
  else if( MatchedTruthIdx==-1 ){
           if     ( NearPhotonType<=0 ) LeptonType=-1;
           else if( NearPhotonType==1 ) LeptonType=-5;
           else if( NearPhotonType==2 ) LeptonType=-6;
         }
  else if( (fabs(MPID)==11 || fabs(MPID)==13 || fabs(MPID)==15) && MStatus_last!=2 && !HadronicOrigin ) LeptonType= 4;//4-a)
  else if( ((fabs(MPID)>=1 && fabs(MPID)<=5) || fabs(MPID)==21) && MStatus_last!=2 )                    LeptonType=-4;//4-b)
  else if( fabs(MPID)==6 ) LeptonType=4;//4-c)
  else LeptonType=0;

 
  return LeptonType;

//**footnote
//These are based on observation in DY,ZG,TT sample(DY,ZG:amcnlo+pythia, TT:powheg+pythia) for other PS generator, convention may differ.
//1) In amcnlo generator, output of ME level generation does not have specific guage field mother. e.g. u u~ > l+ l- -> fabs(MID)=1
//   This perhaps due to multiple field can interplay in production, and apparently it is not possible to distinguish them in any logic.
//   e.g. think about previous example. you cannot say whether this is from gamma or Z or H...
//   But in PS procedure, corrections on ME proc is done sometimes. In that case it seems mother is set Z for OS ll prod. W for lnu prod.
//   e.g. If pythia applies ISR process on input u u~, than it should affect momentum of all the consequent processes, or in case of lnu, W radiating gamma can be added.
//   You may think lnu case is obvious, but it may not like in case of pp > lllnu(You never know which one is from W>lnu and Z>ll)
//2-a) e.g. a)Had > ta+X, ta>l+2nu b) q>ta+X in jet fragmentation (ta is not hardscattered, since it is already considered prev. step)
//2-b) e.g. gamma>ta(+)+ta(-)+X, ta>lnu (St=2)
//2-c) e.g. " " " " " " " " " " " " " " " " " ", but soft gamma case. this is not observed in test sample but put here just in case.
//          (Non hadronic origin since such case already counted before, gamma should be from non-hadronic source)
//2-d) e.g. l>tata..+l.. , ta>l+2nu (Implicit tau conv. from non-hadronic lepton and decay) In implicit conv. GenStatus!=2
//3-a) e.g. hard gamma>ll
//3-b) e.g. a) Had>gamma+X, gamma>ll+X (in PS+Had stage intermediate process is omitted you see just Had>Nphoton+Mhadrons+..)
//          b) q>gamma+q, gamms>ll+X in jet fragmentation or radiations of tops.
//          c) gluon>Ngamma+Mhadrons in jet fragmentation (Actually observed in samples)
//3-c) e.g. W+>W+ gamma, or t>t+gamma, gamma>ll+X, not yet observed in test sample but possible (upto radiation is observed so far)
//3-d) e.g. ta>ta+gamma, gamma>ll+X, tau not from hadron(e.g. pp>tata)
//4-a) e.g. EW lep l, l>lll... just implicit conversion. 
//4-b) e.g. q or g> Nlepton +MHadrons... in parton shower history
//4-c) e.g. t>t+ll.. implicit conversion
}


int AnalyzerCore::GetPhotonType(int PhotonIdx, std::vector<snu::KTruth> TruthColl){
//Type : 
// 0: Invalid input or Error or HardScatter is input when hardscatter is not final state
// 1: HardScatter / 2: Else prompt daughter(l,V,t)
//-1: Reserved for unmatched(Not used now) / -2: Hadronic origin

  if( PhotonIdx<2 ) return 0;
  if( !(TruthColl.at(PhotonIdx).PdgId()==22 && (TruthColl.at(PhotonIdx).GenStatus()==1 || TruthColl.at(PhotonIdx).GenStatus()==23)) ) return 0;

  if(TruthColl.at(PhotonIdx).GenStatus()==23){
    if(IsFinalPhotonSt23(TruthColl)) return 1;
    else                             return 0;
  }//From this pt, only St1 Photon is treated.

  int PhotonType=0;
  int LastSelfIdx    = LastSelfMotherIdx(PhotonIdx,TruthColl);
  int MotherIdx      = FirstNonSelfMotherIdx(PhotonIdx,TruthColl);
  int fMPID=0, Status_orig=0;
  bool HadronicOrigin = false;
    if( PhotonIdx!=-1 ){ Status_orig    = TruthColl.at(LastSelfIdx).GenStatus();
                         HadronicOrigin = HasHadronicAncestor(PhotonIdx, TruthColl);
                       }                           
    if( MotherIdx!=-1 ){ fMPID          = fabs(TruthColl.at(MotherIdx).PdgId()); }


  if     (  Status_orig>20 && Status_orig<30   ) PhotonType= 1;//1)
  else if(       fMPID==23 || fMPID==25        ) PhotonType= 1;//2)
  else if( fMPID==24 || fMPID==6  || fMPID==37 ) PhotonType= 2;//3)
  else if(           HadronicOrigin            ) PhotonType=-2;//4)
  else if( fMPID==11 || fMPID==13 || fMPID==15 ) PhotonType= 2;//5)
  else                                           PhotonType= 0;
  
  return PhotonType;
//**footnote
//1) In case of hard scattered photon, they may have history; GenSt=23>...>GenSt1, And depending on generator, their mother can be explicitly
//  written in history as Z>GG St2 but sometimes their field mother is not designated for avoiding confusion of gauge symmetry.
//  e.g. qq>llG instead of qq>Z>llG
//  To cover all the case, first thing to check is original state of photon is hard scattered or not
//2) Sometimes, if there is no correction on gamma is applied on PS step, photon's final state is 1 before any history.
//   e.g. G;St=1, Mother=Z ; cannot find hard scatter Z / in such case only possibility is to check mother.
//   But in some case, it is not obvious. because in PS step, charged ptls can radiate photons including bosons.
//   So it is artificial to distinguish hard scattered photon and photon from PS step. And fraction of this case is not negilgible; very frequently observed.
//3) top and charged bosons radiate photons, and some case the photon is very energetic.
//4) This category does not include tops. Photons from hadrons and quarks. But predominantly, in most of the cases they are daughter of pi0.
//   But rarely other mesons as eta, B, or even some quarks can also radiate energetic photons.
//5) Photons radiated from lepton FSR, but sometimes they radiate quite energetic photons.
}

//------------------------------------------------------------------------------------------//


//Jihwan Bhyun Modification//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////


std::vector<snu::KJet> AnalyzerCore::SelBJets(std::vector<snu::KJet>& jetColl,TString level){
//Gets Jetcollection as argument and selects Btagged jet in the collection and returns the b tagged jet collection
//If argument jetcoll is in Pt order, then by the algorithm returned bjetcoll will also be in Pt order
   std::vector<snu::KJet> bjetColl;
   for(unsigned int i=0; i<jetColl.size(); i++){
     if(level.Contains("Loose")){ if(jetColl.at(i).BJetTaggerValue(KJet::CSVv2) > 0.5426) bjetColl.push_back(jetColl.at(i));}
     if(level.Contains("Medium")){ if(jetColl.at(i).BJetTaggerValue(KJet::CSVv2) > 0.8484) bjetColl.push_back(jetColl.at(i));}
     if(level.Contains("Tight")){ if(jetColl.at(i).BJetTaggerValue(KJet::CSVv2) > 0.9535) bjetColl.push_back(jetColl.at(i));}
   }
   return bjetColl;
}


std::vector<snu::KJet> AnalyzerCore::SelLightJets(std::vector<snu::KJet>& jetColl, TString level){

    std::vector<snu::KJet> ljetColl;
    for(vector<KJet>::iterator itj=jetColl.begin(); itj!=jetColl.end(); ++itj){
      if(level.Contains("Loose")){ if(itj->BJetTaggerValue(KJet::CSVv2)<0.5426) ljetColl.push_back((*itj));}
      if(level.Contains("Medium")){ if(itj->BJetTaggerValue(KJet::CSVv2)<0.8484) ljetColl.push_back((*itj));}
      if(level.Contains("Tight")){ if(itj->BJetTaggerValue(KJet::CSVv2)<0.9535) ljetColl.push_back((*itj));}
    }
    return ljetColl;
}

std::vector<int> AnalyzerCore::GetSFBJetIdx(std::vector<snu::KJet>& jetColl,TString level){
   std::vector<int> bIdxColl;
   for(unsigned int i=0; i<jetColl.size(); i++){
     if(level.Contains("Medium")){ if(IsBTagged(jetColl.at(i),snu::KJet::CSVv2,snu::KJet::Medium)){ bIdxColl.push_back(i); }}
   }
   for(unsigned int i=0; i<jetColl.size(); i++){
     if(level.Contains("Tight")){ if(IsBTagged(jetColl.at(i),snu::KJet::CSVv2,snu::KJet::Tight)){ bIdxColl.push_back(i); }}
   }

   return bIdxColl;
}


std::vector<int> AnalyzerCore::GetSFLJetIdx(std::vector<snu::KJet>& jetColl, std::vector<int>& bIdxColl,TString level){
   std::vector<int> ljIdxColl;
   if(jetColl.size()<bIdxColl.size()){ cout<<"[Error] Input Wrong to GetSFLJetIdx"<<endl; return ljIdxColl;}
   for(unsigned int i=0; i<jetColl.size(); i++){
     bool trig=true;
     for(int j=0; j<bIdxColl.size(); j++){ if(bIdxColl.at(j)==i) trig=false; }
     if(trig) ljIdxColl.push_back(i);
   }
   return ljIdxColl;
}

std::vector<snu::KElectron> AnalyzerCore::SelEndcapElectrons(std::vector<snu::KElectron>& electronColl){
//19. Jul. 2016, Used for data validation
//selects endcap electrons of the electron collection
   std::vector<snu::KElectron> eECColl;
   for(unsigned int i=0; i<electronColl.size(); i++){
     if(fabs(electronColl.at(i).Eta())>1.479) eECColl.push_back(electronColl.at(i));
   }
   return eECColl;
}


int AnalyzerCore::SumCharge(std::vector<snu::KMuon>& MuonColl){
    int Q=0;
    for(unsigned int i=0; i<MuonColl.size(); i++){Q+=MuonColl.at(i).Charge();}
    return Q;
}


int AnalyzerCore::TriMuChargeIndex(std::vector<snu::KMuon>& MuonColl, TString charge){
    //First Choose 2SS, 1OS muons(++- or --+) SS means 2 of them having same sign, OS means 1 of them having different sign from others.
    // charge="OS" will return the index of muon that having different charge from other 2,
    // charge="SS1" will return the index of first muon that having same sign
    // charge="SS2" will return the index of second muon that having same sign
    int totQ=0, n=0, ns=0;
    for(unsigned int i=0; i<MuonColl.size(); i++){totQ+=MuonColl.at(i).Charge();}

    if(charge.Contains("OS")){
      if(totQ==1) {
         for(unsigned int j=0; j<MuonColl.size(); j++) {if(MuonColl.at(j).Charge()==-1) n=j;}
      }
      if(totQ==-1){
         for(unsigned int j=0; j<MuonColl.size(); j++) {if(MuonColl.at(j).Charge()==1) n=j;}
      }
    }
    if(charge.Contains("SS1")){
      if(totQ==1){
         for(unsigned int j=0; j<MuonColl.size(); j++){
             if(MuonColl.at(j).Charge()==1) {ns+=1; if(ns==1) n=j;}
         }
      }
      if(totQ==-1){
         for(unsigned int j=0; j<MuonColl.size(); j++){
             if(MuonColl.at(j).Charge()==-1) {ns+=1; if(ns==1) n=j;}
         }
      }
    }
    if(charge.Contains("SS2")){
       if(totQ==1){
          for(unsigned int j=0; j<MuonColl.size(); j++){
             if(MuonColl.at(j).Charge()==1) {ns+=1; if(ns==2) n=j;}
          }
       }
       if(totQ==-1){
          for(unsigned int j=0; j<MuonColl.size(); j++){
             if(MuonColl.at(j).Charge()==-1) {ns+=1; if(ns==2) n=j;}
          }
       }
    }
    return n;
}

double AnalyzerCore::GetvPz(snu::KParticle v, snu::KElectron e, int pm){
     double RecoPz=0;
     double X=80.4*80.4/2+e.Px()*v.Px()+e.Py()*v.Py();
     double D=e.E()*e.E()*(X*X-v.Pt()*v.Pt()*e.Pt()*e.Pt());
     if(pm==1) {if(D>=0) {double RecoPzv1 = (e.Pz()*X+sqrt(D))/(e.Pt()*e.Pt()); RecoPz=RecoPzv1;}
                if(D<0)  {double RecoPzv1 = (e.Pz()*X)/(e.Pt()*e.Pt()); RecoPz=RecoPzv1;}
               }
     if(pm==2) {if(D>=0) {double RecoPzv2 = (e.Pz()*X-sqrt(D))/(e.Pt()*e.Pt()); RecoPz=RecoPzv2;}
                if(D<0)  {double RecoPzv2 = (e.Pz()*X)/(e.Pt()*e.Pt()); RecoPz=RecoPzv2;}
               }
     return RecoPz;
}

double AnalyzerCore::GetvPz(snu::KParticle v, snu::KMuon mu, int pm){
     double RecoPz=0;
     double X=80.4*80.4/2+mu.Px()*v.Px()+mu.Py()*v.Py();
     double D=mu.E()*mu.E()*(X*X-v.Pt()*v.Pt()*mu.Pt()*mu.Pt());
     if(pm==1) {if(D>=0) {double RecoPzv1 = (mu.Pz()*X+sqrt(D))/(mu.Pt()*mu.Pt()); RecoPz=RecoPzv1;}
                if(D<0)  {double RecoPzv1 = (mu.Pz()*X)/(mu.Pt()*mu.Pt()); RecoPz=RecoPzv1;}
               }
     if(pm==2) {if(D>=0) {double RecoPzv2 = (mu.Pz()*X-sqrt(D))/(mu.Pt()*mu.Pt()); RecoPz=RecoPzv2;}
                if(D<0)  {double RecoPzv2 = (mu.Pz()*X)/(mu.Pt()*mu.Pt()); RecoPz=RecoPzv2;}
               }
     return RecoPz;
}

double AnalyzerCore::GetAngle(snu::KElectron e, snu::KMuon mu){
    TVector3 v1,v2;
    v1.SetPtEtaPhi(e.Pt(), e.Eta(), e.Phi());  v2.SetPtEtaPhi(mu.Pt(), mu.Eta(), mu.Phi());
    double angle=v1.Angle(v2);
    return angle;
}

double AnalyzerCore::GetAngle(snu::KMuon mu1, snu::KMuon mu2){
    TVector3 v1,v2;
    v1.SetPtEtaPhi(mu1.Pt(), mu1.Eta(), mu1.Phi());  v2.SetPtEtaPhi(mu2.Pt(), mu2.Eta(), mu2.Phi());
    double angle=v1.Angle(v2);
    return angle;
}

double AnalyzerCore::GetAngle(snu::KElectron e, snu::KJet j){
    TVector3 v1,v2;
    v1.SetPtEtaPhi(e.Pt(), e.Eta(), e.Phi());  v2.SetPtEtaPhi(j.Pt(), j.Eta(), j.Phi());
    double angle=v1.Angle(v2);
    return angle;
}

double AnalyzerCore::GetAngle(snu::KMuon mu, snu::KJet j){
    TVector3 v1,v2;
    v1.SetPtEtaPhi(mu.Pt(), mu.Eta(), mu.Phi());  v2.SetPtEtaPhi(j.Pt(), j.Eta(), j.Phi());
    double angle=v1.Angle(v2);
    return angle;
}

double AnalyzerCore::GetAngle(snu::KJet j1, snu::KJet j2){
    TVector3 v1,v2;
    v1.SetPtEtaPhi(j1.Pt(), j1.Eta(), j1.Phi());  v2.SetPtEtaPhi(j2.Pt(), j2.Eta(), j2.Phi());
    double angle=v1.Angle(v2);
    return angle;
}


double AnalyzerCore::dPtRel(snu::KTruth T, snu::KElectron e){
    double X=fabs(e.Pt()-T.Pt())/T.Pt();
    return X;
}

double AnalyzerCore::dPtRel(snu::KTruth T, snu::KMuon mu){
    double X=fabs(mu.Pt()-T.Pt())/T.Pt();
    return X;
}

double AnalyzerCore::dPtRel(snu::KTruth T, snu::KJet j){
    double X=fabs(j.Pt()-T.Pt())/T.Pt();
    return X;
}

/*
double AnalyzerCore::RelPt(snu::KTruth T, snu::KJet j){
  if(j.Pt()==0) return -1;
  double
  
}*/


int AnalyzerCore::GenMatchedIdx(snu::KTruth T, std::vector<snu::KMuon>& MuonColl){
    //Usage: Returns index of the muon in Muoncoll that matches Truth object T
    //Algo: Make a collecton passing minimum cut : dRmax=0.05, dPtRelmax=0.2 
    //**Here I consider dR more imp. => dPtRel is relatively looser.
    //if none, return index -1, if 1, return the index
    //if ambiguous, select one with minimum dRRel relative to dRmaxcut + dPtRel relative to dPtRelmax
    std::vector<int> muCandIdxColl;
    int muFinalCandIdx;
    double tmp, score, dRoCut, dPtReloCut; //oCut:over Cut

    for( int i=0; i<MuonColl.size(); i++){
      if(T.DeltaR(MuonColl.at(i))<0.1) muCandIdxColl.push_back(i);
      //if((T.DeltaR(MuonColl.at(i))<0.1)&&(dPtRel(T,MuonColl.at(i))<0.2)) muCandIdxColl.push_back(i);
    }
    if(muCandIdxColl.size()==0) muFinalCandIdx=-1;
    else if(muCandIdxColl.size()==1) muFinalCandIdx=muCandIdxColl.at(0);
    else if(muCandIdxColl.size()>1){
       for( int i=0; i<muCandIdxColl.size(); i++){
          dRoCut=T.DeltaR(MuonColl.at(muCandIdxColl.at(i)))/0.05;
          dPtReloCut=0.;
          //dPtReloCut=dPtRel(T,MuonColl.at(muCandIdxColl.at(i)))/0.2;
          tmp=dRoCut+dPtReloCut;
          if(i==0) {muFinalCandIdx=muCandIdxColl.at(i); score=dRoCut+dPtReloCut;}
          else if(tmp<score) {muFinalCandIdx=muCandIdxColl.at(i); score=tmp;}
       }
    }

    return muFinalCandIdx;
}

int AnalyzerCore::GenMatchedIdx(snu::KTruth T, std::vector<snu::KElectron>& ElectronColl){
    //Usage: Returns index of the electtron in Electroncoll that matches Truth object T
    //Algo: Make a collecton passing minimum cut : dRmax=0.05, dPtRelmax=0.2 
    //**Here I consider dR more imp. => dPtRel is relatively looser.
    //if none, return index -1, if 1, return the index
    //if ambiguous, select one with minimum dRRel relative to dRmaxcut + dPtRel relative to dPtRelmax
    std::vector<int> eCandIdxColl;
    int eFinalCandIdx;
    double tmp, score, dRoCut, dPtReloCut; //oCut:over Cut

    for( int i=0; i<ElectronColl.size(); i++){
      if((T.DeltaR(ElectronColl.at(i))<0.1)&&(dPtRel(T,ElectronColl.at(i))<0.2)) eCandIdxColl.push_back(i);
    }
    if(eCandIdxColl.size()==0) eFinalCandIdx=-1;
    else if(eCandIdxColl.size()==1) eFinalCandIdx=eCandIdxColl.at(0);
    else if(eCandIdxColl.size()>1){
       for( int i=0; i<eCandIdxColl.size(); i++){
          dRoCut=T.DeltaR(ElectronColl.at(eCandIdxColl.at(i)))/0.05;
          dPtReloCut=dPtRel(T,ElectronColl.at(eCandIdxColl.at(i)))/0.2;
          tmp=dRoCut+dPtReloCut;
          if(i==0) {eFinalCandIdx=eCandIdxColl.at(i); score=dRoCut+dPtReloCut;}
          else if(tmp<score) {eFinalCandIdx=eCandIdxColl.at(i); score=tmp;}
       }
    }

    return eFinalCandIdx;
}

int AnalyzerCore::GenMatchedIdx(snu::KTruth T, std::vector<snu::KJet>& JetColl){
    //Usage: Returns index of the jet in Jetcoll that matches Truth object T
    //Algo: Make a collecton passing minimum cut : dRmax=0.5, dPtRelmax=0.2 
    //**Here I consider dR more imp. => dPtRel is relatively looser.
    //if none, return index -1, if 1, return the index
    //if ambiguous, select one with minimum dRRel relative to dRmaxcut + dPtRel relative to dPtRelmax
    std::vector<int> bjCandIdxColl;
    std::vector<int> ljCandIdxColl;
    int jFinalCandIdx;
    double tmp, score, dRoCut, dPtReloCut; //oCut:over Cut

    for( int i=0; i<JetColl.size(); i++){
      if((T.DeltaR(JetColl.at(i))<0.5)&&(dPtRel(T,JetColl.at(i))<0.7)){
         if(JetColl.at(i).IsBTagged(KJet::CSVv2,KJet::Medium)) bjCandIdxColl.push_back(i);
         else ljCandIdxColl.push_back(i);
      }
    }

    //B Matcher Algo
    if(fabs(T.PdgId())==5){
      if((bjCandIdxColl.size()==0)&&(ljCandIdxColl.size()==0)) jFinalCandIdx=-1;
      else if((bjCandIdxColl.size()==0)&&(ljCandIdxColl.size()==1)) jFinalCandIdx=ljCandIdxColl.at(0);
      else if((bjCandIdxColl.size()==0)){// 0 bjCandColl , >1 ljCandColl
         for( int i=0; i<ljCandIdxColl.size(); i++){
            dRoCut=T.DeltaR(JetColl.at(ljCandIdxColl.at(i)))/0.5;
            dPtReloCut=dPtRel(T,JetColl.at(ljCandIdxColl.at(i)))/0.7;
            tmp=dRoCut+dPtReloCut;
            if(i==0) {jFinalCandIdx=ljCandIdxColl.at(i); score=dRoCut+dPtReloCut;}
            else if(tmp<score) {jFinalCandIdx=ljCandIdxColl.at(i); score=tmp;}
         }
      }
      else if(bjCandIdxColl.size()==1) jFinalCandIdx=bjCandIdxColl.at(0);
      else {//>1 bjCand, & any # ljCand
         for( int i=0; i<bjCandIdxColl.size(); i++){
            dRoCut=T.DeltaR(JetColl.at(bjCandIdxColl.at(i)))/0.5;
            dPtReloCut=dPtRel(T,JetColl.at(bjCandIdxColl.at(i)))/0.7;
            tmp=dRoCut+dPtReloCut;
            if(i==0) {jFinalCandIdx=bjCandIdxColl.at(i); score=dRoCut+dPtReloCut;}
            else if(tmp<score) {jFinalCandIdx=bjCandIdxColl.at(i); score=tmp;}
         }
      }
   }//B matcher Algo
   //light jet matcher
   else{
     if(ljCandIdxColl.size()==0) jFinalCandIdx=-1;
     else if(ljCandIdxColl.size()==1) jFinalCandIdx=ljCandIdxColl.at(0);
     else{
         for( int i=0; i<ljCandIdxColl.size(); i++){
            dRoCut=T.DeltaR(JetColl.at(ljCandIdxColl.at(i)))/0.5;
            dPtReloCut=dPtRel(T,JetColl.at(ljCandIdxColl.at(i)))/0.7;
            tmp=dRoCut+dPtReloCut;
            if(i==0) {jFinalCandIdx=ljCandIdxColl.at(i); score=dRoCut+dPtReloCut;}
            else if(tmp<score) {jFinalCandIdx=ljCandIdxColl.at(i); score=tmp;}
         }
     }
   }//light jet matcher

    return jFinalCandIdx;
}


bool AnalyzerCore::GenDecayInfo(std::vector<snu::KTruth>& TruthColl, TString Option){
//Option: Decay3lv  ; return true if hc decays 3l+v(1e2mu or 3mu or 1ta2mu)
//        DecayPr3lv; return true if hc decays to 1e2mu or 3mu
//        Decayta3lv: return true if hc decays to 1ta2mu

  bool DecayPr3lv=false, Decayta3lv=false, Decay3lv=false, Decay2l2j=false, Decaytlv=false, Decaytatlv=false;
  bool TrigPr3lv = Option.Contains("DecayPr3lv"), Trigta3lv = Option.Contains("Decayta3lv");
  bool Trig3lv   = Option.Contains("Decay3lv"),   Trig2l2j  = Option.Contains("Decay2l2j"), Trig4l = Option.Contains("Decay4l");
  bool Trigtatlv = Option.Contains("Decaytalv");
  bool Trigger   = false;

  for(int i=2; i<(int) TruthColl.size(); i++){
    if(TruthColl.at(i).IndexMother()==-1) continue;
      int fpid=abs(TruthColl.at(i).PdgId());
      int GenSt=TruthColl.at(i).GenStatus();
    if( !(GenSt==1 || GenSt==2 || (GenSt>20 && GenSt<30)) ) continue;
    if( !(fpid==11 || fpid==13 || fpid==15 || fpid<6) ) continue;

    int MotherIdx       = FirstNonSelfMotherIdx(i,TruthColl);
    int GrMotherIdx     = FirstNonSelfMotherIdx(MotherIdx,TruthColl);
    int MPID=0, GrMPID=0;
    if( MotherIdx!=-1     ) MPID     = TruthColl.at(MotherIdx).PdgId();
    if( GrMotherIdx!=-1   ) GrMPID   = TruthColl.at(GrMotherIdx).PdgId();

    if( fpid<6 && abs(MPID)==24  && abs(GrMPID)==37                  )  Decay2l2j=true;
    if( (fpid==11 || fpid==13 )   && abs(MPID)==24 && abs(GrMPID)==37 ){ Decay3lv =true; DecayPr3lv=true; }
    if( fpid==15 && abs(MPID)==24 && abs(GrMPID)==37                  ){ Decay3lv =true; Decayta3lv=true; }
    if( (fpid==11 || fpid==13 || fpid==15) && abs(MPID)==24 && abs(GrMPID)==6 ) Decaytlv=true;
    if( fpid==15 && abs(MPID)==24 && abs(GrMPID)==6 ) Decaytatlv=true;
  }

  if     (Trig3lv  ) Trigger = Decay3lv;
  else if(TrigPr3lv) Trigger = DecayPr3lv;
  else if(Trigta3lv) Trigger = Decayta3lv;
  else if(Trig2l2j ) Trigger = Decay2l2j;
  else if(Trig4l   ) Trigger = Decay3lv and Decaytlv;
  else if(Trigtatlv) Trigger = Decaytatlv;

  return Trigger;
}


int AnalyzerCore::GetGenMatchedSigIndex(std::vector<snu::KTruth>& truthColl, std::vector<snu::KMuon>& muonColl, std::vector<snu::KElectron>& electronColl, std::vector<snu::KJet>& jetColl, TString PtlName, float weight){
// Possible PtlName list : "mum_A", "mup_A", "e_W" , "mu_W", "j1_W", "j2_W", "b_t" , "bx_tx" 

//////////////////////////////////////////////////////////////////////////////////////////////
////Gen Matching//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

//Truth Objects/////////
   int count=0;
   snu::KTruth b_t, bx_tx, mum_A, mup_A; // Doesn't dep. on decaymode
   snu::KTruth e_W, ve_W, mu_W, vmu_W, j1_W, j2_W;//dep. on decaymode
//   snu::KTruth tmpj1_W_tx, tmpj2_W_tx;
   bool Decay3lv=false, emumu_gen=false, trimu_gen=false;
   bool emumu=false, trimu=false;
   int mum_Ai, mup_Ai, e_Wi, mu_Wi, j1_Wi, j2_Wi, b_ti, bx_txi;

   if((muonColl.size()==2)&&(electronColl.size()==1)) emumu=true;
   else if((muonColl.size()==3)&&(electronColl.size()==0)) trimu=true;


   for(unsigned int i=0; i<truthColl.size(); i++){
     if(truthColl.at(i).IndexMother()==-1) continue;
     int pid=truthColl.at(i).PdgId();
     int MotherIdx=truthColl.at(i).IndexMother();
     int Motherpid=truthColl.at(MotherIdx).PdgId();
//     if(truthColl.at(MotherIdx).IndexMother()==-1) continue;
//     int GrmaIdx=truthColl.at(MotherIdx).IndexMother();
//     int Grmapid=truthColl.at(GrmaIdx).PdgId();
//     Actually ^ is not needed since we only generated H+ not H- > W+ from H+ && W- from tx
     //H+ side(t side); t> b H+ > b (A w+)
     if(Motherpid==6){
       if(pid==5) b_t=truthColl.at(i);
     }
     else if(Motherpid==36){
       if(pid==13) mum_A=truthColl.at(i);
       else if(pid==-13) mup_A=truthColl.at(i);
     }
     else if(Motherpid==24){
       if((fabs(pid)>10)&&(fabs(pid)<20)){
         Decay3lv=true;
         if(pid==-11) {e_W=truthColl.at(i); emumu_gen=true;}
         else if(pid==-13) {mu_W=truthColl.at(i); trimu_gen=true;}
         else if(pid==12) ve_W=truthColl.at(i);
         else if(pid==14) vmu_W=truthColl.at(i);
       }
       else if(fabs(pid)<10){
         Decay3lv=false;
         if(count==0) j1_W=truthColl.at(i);
         else j2_W=truthColl.at(i);
         count++;
       }
     }

     //tx side ; tx>bx w-
     else if(Motherpid==-6){
         if(pid==-5) bx_tx=truthColl.at(i);
     }
     else if(Motherpid==-24){
       if((fabs(pid)>10)&&(fabs(pid)<20)){
         if(pid==11) {e_W=truthColl.at(i); emumu_gen=true;}
         else if(pid==13) {mu_W=truthColl.at(i); trimu_gen=true;}
         else if(pid==-12) ve_W=truthColl.at(i);
         else if(pid==-14) vmu_W=truthColl.at(i);
       }
       else if(fabs(pid)<10){
         if(count==0) j1_W=truthColl.at(i);
         else j2_W=truthColl.at(i);
         count++;
       }
     }

   }


//Gen-Reco Matching/////////
   //if((emumu)&&(emumu_gen)){
   if(emumu){
      if(Decay3lv){
        //top-H+ side
        mum_Ai=GenMatchedIdx(mum_A, muonColl);  mup_Ai=GenMatchedIdx(mup_A, muonColl);
        if(emumu_gen){e_Wi=GenMatchedIdx(e_W, electronColl);} else e_Wi=-1;
        b_ti=GenMatchedIdx(b_t, jetColl);

        //tx side
        bx_txi=GenMatchedIdx(bx_tx, jetColl);
        j1_Wi=GenMatchedIdx(j1_W, jetColl);  j2_Wi=GenMatchedIdx(j2_W, jetColl);
     }
     else{
        //top-H+ side
        mum_Ai=GenMatchedIdx(mum_A, muonColl); mup_Ai=GenMatchedIdx(mup_A, muonColl);
        b_ti=GenMatchedIdx(b_t, jetColl);
        j1_Wi=GenMatchedIdx(j1_W, jetColl); j2_Wi=GenMatchedIdx(j2_W, jetColl);

        //tx side
        bx_txi=GenMatchedIdx(bx_tx, jetColl);
        if(emumu_gen){ e_Wi=GenMatchedIdx(e_W, electronColl);} else e_Wi=-1;
     }

     //Filtering and cutflows
     FillHist("Basic_GenMatchFlow_emu", 0, weight, 0., 6., 6);//Nocut
     if((mum_Ai!=-1)&&(mup_Ai!=-1)){ FillHist("Basic_GenMatchFlow_emu", 1, weight, 0., 6., 6);//step1
        if(e_Wi!=-1){ FillHist("Basic_GenMatchFlow_emu", 2, weight, 0., 6., 6);//step2
           if(mum_Ai!=mup_Ai){ FillHist("Basic_GenMatchFlow_emu", 3, weight, 0., 6., 6);
             if((j1_Wi!=-1)&&(j2_Wi!=-1)){ FillHist("Basic_GenMatchFlow_emu", 4, weight, 0., 6., 6);//step3
                if((bx_txi!=-1)&&(b_ti!=-1)) FillHist("Basic_GenMatchFlow_emu", 5, weight, 0., 6., 6);//step4
             }
           }
        }
     }

     if     (PtlName=="mum_A") return mum_Ai;
     else if(PtlName=="mup_A") return mup_Ai;
     else if(PtlName=="e_W")   return e_Wi;
     else if(PtlName=="mu_W")  return -1;
     else if(PtlName=="j1_W")  return j1_Wi;
     else if(PtlName=="j2_W")  return j2_Wi;
     else if(PtlName=="b_t")   return b_ti;
     else if(PtlName=="bx_tx") return bx_txi;
   }
   else if(trimu){
      if(Decay3lv){
        //top-H+ side
        mum_Ai=GenMatchedIdx(mum_A, muonColl); mup_Ai=GenMatchedIdx(mup_A, muonColl);
        if(trimu_gen){ mu_Wi=GenMatchedIdx(mu_W, muonColl);} else mu_Wi=-1;
        b_ti=GenMatchedIdx(b_t, jetColl);

        //tx side
        bx_txi=GenMatchedIdx(bx_tx, jetColl);
        j1_Wi=GenMatchedIdx(j1_W, jetColl); j2_Wi=GenMatchedIdx(j2_W, jetColl);
     }
     else{
        //top-H+ side
        mum_Ai=GenMatchedIdx(mum_A, muonColl); mup_Ai=GenMatchedIdx(mup_A, muonColl);
        b_ti=GenMatchedIdx(b_t, jetColl);
        j1_Wi=GenMatchedIdx(j1_W, jetColl); j2_Wi=GenMatchedIdx(j2_W, jetColl);

        //tx side
        bx_txi=GenMatchedIdx(bx_tx, jetColl);
        if(trimu_gen){ mu_Wi=GenMatchedIdx(mu_W, muonColl);} else mu_Wi=-1;
     }

     //Filtering and cutflows
     FillHist("Basic_GenMatchFlow_3mu", 0, weight, 0., 7., 7);//Nocut
     if((mum_Ai!=-1)&&(mup_Ai!=-1)){ FillHist("Basic_GenMatchFlow_3mu", 1, weight, 0., 7., 7);//step1
        if(mu_Wi!=-1){ FillHist("Basic_GenMatchFlow_3mu", 2, weight, 0., 7., 7);//step2
          if(mum_Ai!=mup_Ai){ FillHist("Basic_GenMatchFlow_3mu", 3, weight, 0., 7., 7);
            if((mu_Wi!=mum_Ai)&&(mu_Wi!=mup_Ai)){ FillHist("Basic_GenMatchFlow_3mu", 4, weight, 0., 7., 7);
               if((j1_Wi!=-1)&&(j2_Wi!=-1)){ FillHist("Basic_GenMatchFlow_3mu", 5, weight, 0., 7., 7);//step3
                  if((bx_txi!=-1)&&(b_ti!=-1)) FillHist("Basic_GenMatchFlow_3mu", 6, weight, 0., 7., 7);//step4
               }
            }
         }
       }
     }
     if     (PtlName=="mum_A") return mum_Ai;
     else if(PtlName=="mup_A") return mup_Ai;
     else if(PtlName=="e_W")   return -1;
     else if(PtlName=="mu_W")  return mu_Wi;
     else if(PtlName=="j1_W")  return j1_Wi;
     else if(PtlName=="j2_W")  return j2_Wi;
     else if(PtlName=="b_t")   return b_ti;
     else if(PtlName=="bx_tx") return bx_txi;
   }//

//////////////////////////////////////////////////////////////////////////////////////
return -1;
//End of function
}

int AnalyzerCore::NPromptLeptons(std::vector<snu::KTruth>& truthColl, TString Option){
  
  int NPromptLepton_Tot=0, NPromptLepton_EW=0, NPromptLepton_BSM=0;
  int NTauHard=0;
  bool InAcceptance=(Option.Contains("InAcceptance") || Option.Contains("_A"));
  bool InclTauLep  =(Option.Contains("TauLep")       || Option.Contains("_TL"));
  bool InclTau     =(Option.Contains("InclTau")      || Option.Contains("_T"));
  bool OnlyConv    =(Option.Contains("OnlyConv")     || Option.Contains("_C"));
    if(InclTau && InclTauLep) InclTauLep=false;//To avoid double counting;
  std::vector<int> PromptIdxColl;

  //Prompt Count
  for(unsigned int i=2; i<truthColl.size(); i++){
    if( truthColl.at(i).IndexMother()<0 )  continue;
     int fpid=fabs(truthColl.at(i).PdgId());
    if( !(fpid==11 || fpid==13 || fpid==15) ) continue;
    if( fpid==15 ){
       int mpid=fabs(truthColl.at(truthColl.at(i).IndexMother()).PdgId());
       if((mpid==23 || mpid==24) && truthColl.at(i).GenStatus()==2)              NTauHard++;
       else if(truthColl.at(i).GenStatus()>20 && truthColl.at(i).GenStatus()<30) NTauHard++;
    }
    if( !((fpid==11 || fpid==13) && truthColl.at(i).GenStatus()==1 ) )  continue;

    int LepType=GetLeptonType(i, truthColl);

    if(!InAcceptance){
      if(LepType==1) {NPromptLepton_EW++; PromptIdxColl.push_back(i);}
      if(LepType==2) {NPromptLepton_BSM++; PromptIdxColl.push_back(i);}
      if(InclTauLep && LepType==3) NPromptLepton_EW++;
    }
    else{
      if( (fpid==11 && fabs(truthColl.at(i).Eta())<2.5) 
         || (fpid==13 && fabs(truthColl.at(i).Eta())<2.4) ){
        if(LepType==1) NPromptLepton_EW++;
        if(LepType==2) NPromptLepton_BSM++;
        if(InclTauLep && LepType==3) NPromptLepton_EW++;
      }
    }
  }

  //Conversion Count
  std::vector<int> ConvBranchIdxColl;
  for(int i=0; i<PromptIdxColl.size(); i++){
    for(int j=i+1; j<PromptIdxColl.size(); j++){
      int pid1 =truthColl.at(PromptIdxColl.at(i)).PdgId();
      int midx1=truthColl.at(PromptIdxColl.at(i)).IndexMother(), mid1=truthColl.at(midx1).PdgId();
      while(mid1==pid1){
        int pid2 =truthColl.at(PromptIdxColl.at(j)).PdgId();
        int midx2=truthColl.at(PromptIdxColl.at(j)).IndexMother(), mid2=truthColl.at(midx2).PdgId();
        while(mid2==pid2){
          if(midx1==midx2) break;
          else {midx2=truthColl.at(midx2).IndexMother(); mid2=truthColl.at(midx2).PdgId();}
        }
        if(midx1==midx2){ ConvBranchIdxColl.push_back(midx1); break; }
        else {midx1=truthColl.at(midx1).IndexMother(); mid1=truthColl.at(midx1).PdgId();}
      }
    }
  }
  int NConv=ConvBranchIdxColl.size(), NMultiConv=0, SameBranchCount=0;
  for(int i=0; i<ConvBranchIdxColl.size(); i++){
    for(int j=i+1; j<ConvBranchIdxColl.size(); j++){
      if(ConvBranchIdxColl.at(i)==ConvBranchIdxColl.at(j)) SameBranchCount++;
    }
  }
  if(SameBranchCount==3) NMultiConv++;//Basically we should also count all the nC2 (n=3,4...) but higher order terms need not be considered till now yet.
  NConv-=NMultiConv;

  NPromptLepton_Tot=NPromptLepton_EW+NPromptLepton_BSM-NConv;
  if(InclTau) NPromptLepton_Tot+=NTauHard;

  if     (Option.Contains("EW"))  return (NPromptLepton_Tot-NPromptLepton_BSM);
  else if(Option.Contains("BSM")) return NPromptLepton_BSM;
  else if(OnlyConv)               return NConv;

  return NPromptLepton_Tot;
}

int AnalyzerCore::NLeptonicBosonDecay(std::vector<snu::KTruth>& TruthColl){

   int Counter=0;
   for(int i=2; i<TruthColl.size(); i++){
      int fpid=fabs(TruthColl.at(i).PdgId());     
     if( !(fpid==11 || fpid==13 || fpid==15) ) continue;
      int mfpid=fabs(TruthColl.at(TruthColl.at(i).IndexMother()).PdgId());
      int GenSt=TruthColl.at(i).GenStatus();
     if( GenSt>20 && GenSt<30 ) Counter++;
     else if( mfpid==23 || mfpid==24 || mfpid==25 ) Counter++;
   }

   return Counter;
}


int AnalyzerCore::GenMatchedIdx(snu::KJet Jet, std::vector<snu::KTruth>& TruthColl){
  //Seeding from jet, find parton matched to jet within dRmax;
  //Does not return matched index if not matched or matched multiple partons.(-1:Unmatched, -2:Multiple)

  int MatchedIdx=-1;
  float dR=999., dRmax=0.5;
  int CountMatched=0;
  for(int i=2; i<TruthColl.size(); i++){
    if( TruthColl.at(i).IndexMother()<0) continue;
    if( TruthColl.at(i).GenStatus()!=23 ) continue;
      int fpid=fabs(TruthColl.at(i).PdgId());
    if( !((fpid>0 && fpid<=5) || fpid==21) ) continue;
    if( TruthColl.at(i).DeltaR(Jet)>dRmax )  continue;
    if( CountMatched>1 ) return -2;

    if(TruthColl.at(i).DeltaR(Jet)<dR){ dR=TruthColl.at(i).DeltaR(Jet); MatchedIdx=i; CountMatched++;}
  }

  return MatchedIdx;
}


bool AnalyzerCore::IsHardPhotonConverted(std::vector<snu::KTruth> TruthColl){

  bool IsConverted=false;
  for(int i=2; i<TruthColl.size(); i++){
    int fpid  = fabs(TruthColl.at(i).PdgId());
    int MPID_direct= TruthColl.at(TruthColl.at(i).IndexMother()).PdgId();
    if( !(fpid!=22 && MPID_direct==22) ) continue;

    int MotherIdx  = FirstNonSelfMotherIdx(i,TruthColl);
    int LastSelfMIdx=-1, MStatus_orig=-1;
    if(MotherIdx!=-1){
      LastSelfMIdx = LastSelfMotherIdx(MotherIdx,TruthColl);
      MStatus_orig = TruthColl.at(LastSelfMIdx).GenStatus();
    }

    if( MPID_direct==22 && MStatus_orig==23 ){ IsConverted=true; break; }
  }

  return IsConverted;
}


bool AnalyzerCore::IsJetConsistentPartonHadronMatch(snu::KJet Jet, std::vector<snu::KTruth>& TruthColl, TString Option){

  bool IsConsistent=false;
  bool AllFlav=false, BFlav=false, HeavyFlav=false;
    if     (Option.Contains("AllFlav")) AllFlav  =true;
    else if(Option.Contains("BFlav"))   BFlav    =true;
    else if(Option.Contains("Heavy"))   HeavyFlav=true;
  int MatchedIdx=GenMatchedIdx(Jet,TruthColl);
  int MatchedPartonPID=0;
  int JetHadronFlav=Jet.HadronFlavour();
    if(MatchedIdx>=0) MatchedPartonPID=TruthColl.at(MatchedIdx).PdgId();

  if(AllFlav){
    if     ( JetHadronFlav==5 && fabs(MatchedPartonPID)==5 ) IsConsistent=true;
    else if( JetHadronFlav==4 && fabs(MatchedPartonPID)==4 ) IsConsistent=true;
    else if( JetHadronFlav==0 && ((fabs(MatchedPartonPID)<4 && MatchedPartonPID!=0) || MatchedPartonPID==21) ) IsConsistent=true;
  }
  else if(BFlav){
    IsConsistent=true;
    if     ( JetHadronFlav==5 && fabs(MatchedPartonPID)!=5 ) IsConsistent=false;
    else if( JetHadronFlav!=5 && fabs(MatchedPartonPID)==5 ) IsConsistent=false;
  }
  else if(HeavyFlav){
    IsConsistent=true;
    if     ( JetHadronFlav==5 &&   fabs(MatchedPartonPID)!=5 ) IsConsistent=false;
    else if( JetHadronFlav==4 && !(fabs(MatchedPartonPID)==5 || fabs(MatchedPartonPID)==4) ) IsConsistent=false;
    else if( JetHadronFlav==0 &&  (fabs(MatchedPartonPID)==5 || fabs(MatchedPartonPID)==4) ) IsConsistent=false;
  }

  return IsConsistent;
}

bool AnalyzerCore::HasEWLepInJet(snu::KJet Jet, std::vector<snu::KTruth>& TruthColl, TString Option){

  int HasEWLep=false;
  bool InclEWTau=Option.Contains("InclTau");
  for(int i=2; i<TruthColl.size(); i++){
    if(TruthColl.at(i).IndexMother()<0 )  continue;
    bool IsCand=false;
    int abspid=fabs(TruthColl.at(i).PdgId());
    int absmpid=fabs(TruthColl.at(TruthColl.at(i).IndexMother()).PdgId());
    if( TruthColl.at(i).GenStatus()==1 ){
      if( abspid==11 || abspid==13 ){
        int LepType=GetLeptonType(i, TruthColl);
        if( LepType==1 || LepType==2 || LepType==3 ) IsCand=true; 
      }
    }
    if( InclEWTau && abspid==15 && TruthColl.at(i).GenStatus()>20 && TruthColl.at(i).GenStatus()<30 ) IsCand=true;
    else if( InclEWTau && abspid==15 && (absmpid==23 || absmpid==24) && TruthColl.at(i).GenStatus()==2 ) IsCand=true;

    if( IsCand && TruthColl.at(i).DeltaR(Jet)<0.4) HasEWLep=true;
    if( HasEWLep ) break;
  }

  return HasEWLep;
}

bool AnalyzerCore::NearEWLep(snu::KElectron Ele, std::vector<snu::KTruth>& TruthColl, TString Option){

  bool NearEWLep=false;
  bool InclEWTau=Option.Contains("InclTau");
  float PTthreshold=10.; //Intended for suppressing conversion fake;
  for(int i=2; i<TruthColl.size(); i++){
    if(TruthColl.at(i).IndexMother()<0 ) continue;
    //if(TruthColl.at(i).Pt()<PTthreshold) continue;
    bool IsCand=false;
    int abspid=fabs(TruthColl.at(i).PdgId());
    int absmpid=fabs(TruthColl.at(TruthColl.at(i).IndexMother()).PdgId());
    if( TruthColl.at(i).GenStatus()==1 ){
      if( abspid==11 || abspid==13 ){
        int LepType=GetLeptonType(i, TruthColl);
        if( LepType==1 || LepType==2 || LepType==3 ) IsCand=true; 
      }
    }
    if( InclEWTau && abspid==15 && TruthColl.at(i).GenStatus()>20 && TruthColl.at(i).GenStatus()<30 ) IsCand=true;
    else if( InclEWTau && abspid==15 && (absmpid==23 || absmpid==24) && TruthColl.at(i).GenStatus()==2 ) IsCand=true;

    if( IsCand && TruthColl.at(i).DeltaR(Ele)<0.4) NearEWLep=true;
    if( NearEWLep ) break;
  }

  return NearEWLep;
}


bool AnalyzerCore::NearPhoton(snu::KElectron Ele, std::vector<snu::KTruth>& TruthColl, TString Option){//1)

  bool NearPhoton=false;
  bool OnlyHardPhoton=false;
    if(Option=="Hard") OnlyHardPhoton=true;
  float PTthreshold=10.;
  float dRmax=0.2;//2)

  for(int i=2; i<TruthColl.size(); i++){
    if( TruthColl.at(i).IndexMother()<0   ) continue;
    if( TruthColl.at(i).PdgId()!=22       ) continue;
    if( !(TruthColl.at(i).GenStatus()==1 || TruthColl.at(i).GenStatus()==23) ) continue;
    if( TruthColl.at(i).Pt()<PTthreshold  ) continue;
    if( !(Ele.Pt()/TruthColl.at(i).Pt()>0.8 && Ele.Pt()/TruthColl.at(i).Pt()<1.2) ) continue;//3)
    if( Ele.DeltaR(TruthColl.at(i))>dRmax ) continue;

    if( GetPhotonType(i,TruthColl)>0      ){ NearPhoton=true; break; }
  }

  return NearPhoton;
//footnote
//1) When checked with ZG sample with CBPOGT ElePt>25, Hardscattered photons, External conversion is only meaningful for electrons.
//   Electron External Conversion~3000 Muon external conversion. This is in agreement with theoretical calculation that xsec~M^{-2} in asymmetric limit.
//   ref. Arxiv:1110.1368v1
//2) When checked with ZG sample with CBPOGT ElePt>25GeV, HardScatter Photons, Matching Eff~99.7% even when tested with dR01 cone.
//   But just used conservative cone size.
//3) When checked with ZG sample with CBPOGT ElePt>25GeV, HardScatter Photons, Matching Eff~99.5%
//   External Conversion object's momentum is roughly symmetric with mother photon's momentum(similar fraction lower/upper than PT(G))
}


std::vector<snu::KElectron> AnalyzerCore::SkimLepColl(std::vector<snu::KElectron>& EleColl, std::vector<snu::KTruth>& TruthColl, TString Option){

  bool GetPrompt=false, GetHadFake=false, GetEWtau=false, GetNHIntConv=false, GetNHExtConv=false;
  if(Option.Contains("Prompt"))          GetPrompt    =true;
  if(Option.Contains("HFake"))           GetHadFake   =true;
  if(Option.Contains("EWtau"))           GetEWtau     =true;
  if(Option.Contains("NHConv"))         {GetNHIntConv =true; GetNHExtConv=true;}
  else{ if(Option.Contains("NHIntConv")) GetNHIntConv =true;
        if(Option.Contains("NHExtConv")) GetNHExtConv =true; }
  if(     Option=="Fake"     )          {GetHadFake   =true; GetNHExtConv=true;}


  std::vector<snu::KElectron> ReturnVec;
  for(int i=0; i<EleColl.size(); i++){
    //if(GetHadFake){
      //if(NearEWLep(EleColl.at(i), TruthColl, "InclTau")) continue;
      //if(NearPhoton(EleColl.at(i), TruthColl))           continue;
    //}

    int LepType=GetLeptonType(EleColl.at(i), TruthColl);
    if( GetPrompt    && (LepType==1 || LepType==2) ) ReturnVec.push_back(EleColl.at(i));
    if( GetHadFake   && (LepType<0 && LepType>=-4) ) ReturnVec.push_back(EleColl.at(i));
    if( GetEWtau     &&         LepType==3         ) ReturnVec.push_back(EleColl.at(i));
    if( GetNHIntConv &&         LepType>=4         ) ReturnVec.push_back(EleColl.at(i));
    if( GetNHExtConv &&         LepType<-4         ) ReturnVec.push_back(EleColl.at(i));
  }

  return ReturnVec;
}


std::vector<snu::KElectron> AnalyzerCore::SkimLepColl(std::vector<snu::KElectron>& EleColl, TString Option, float PTmin){
  
  std::vector<snu::KElectron> ReturnColl;
  bool Barrel1=false, Barrel2=false, Endcap=false, PtCut=false;
  if(Option.Contains("B1")) Barrel1=true;
  if(Option.Contains("B2")) Barrel2=true;
  if(Option.Contains("E"))  Endcap =true;
  if(Option.Contains("Pt")) PtCut  =true;

  for(int i=0; i<EleColl.size(); i++){
    bool PassSel=false;
    if( PtCut   && EleColl.at(i).Pt()<PTmin ) continue;
    if( Barrel1 && fabs(EleColl.at(i).Eta())<0.8 ) ReturnColl.push_back(EleColl.at(i));
    if( Barrel2 && fabs(EleColl.at(i).Eta())>=0.8 && fabs(EleColl.at(i).Eta())<1.479 ) ReturnColl.push_back(EleColl.at(i));
    if( Endcap  && fabs(EleColl.at(i).Eta())>=1.479 && fabs(EleColl.at(i).Eta())<2.5 ) ReturnColl.push_back(EleColl.at(i));
  }

  return ReturnColl;
}


std::vector<snu::KMuon> AnalyzerCore::SkimLepColl(std::vector<snu::KMuon>& MuColl, std::vector<snu::KTruth>& TruthColl, TString Option){

  bool GetPrompt=false, GetHadFake=false, GetEWtau=false, GetNHIntConv=false, GetNHExtConv=false;
  if(Option.Contains("Prompt"))          GetPrompt    =true;
  if(Option.Contains("HFake"))           GetHadFake   =true;
  if(Option.Contains("EWtau"))           GetEWtau     =true;
  if(Option.Contains("NHConv"))         {GetNHIntConv =true; GetNHExtConv=true;}
  else{ if(Option.Contains("NHIntConv")) GetNHIntConv =true;
        if(Option.Contains("NHExtConv")) GetNHExtConv =true; }
  if(     Option=="Fake"     )          {GetHadFake   =true; GetNHExtConv=true;}


  std::vector<snu::KMuon> ReturnVec;
  for(int i=0; i<MuColl.size(); i++){
    int LepType=GetLeptonType(MuColl.at(i), TruthColl);
    if( GetPrompt    && (LepType==1 || LepType==2) ) ReturnVec.push_back(MuColl.at(i));
    if( GetHadFake   && (LepType<0 && LepType>=-4) ) ReturnVec.push_back(MuColl.at(i));
    if( GetEWtau     &&         LepType==3         ) ReturnVec.push_back(MuColl.at(i));
    if( GetNHIntConv &&         LepType>=4         ) ReturnVec.push_back(MuColl.at(i));
    if( GetNHExtConv &&         LepType<-4         ) ReturnVec.push_back(MuColl.at(i));
  }

  return ReturnVec;
}

std::vector<snu::KJet> AnalyzerCore::SkimJetColl(std::vector<snu::KJet>& JetColl, std::vector<snu::KTruth>& TruthColl, TString Option){

  bool GetTrueJet=false, GetFakeJet=false, GetPrLepCleanJet=false, ExcTau=false;
  TString Criteria="";
  if(Option.Contains("True"))   GetTrueJet       =true;
  if(Option.Contains("Fake"))   GetFakeJet       =true;
  if(Option.Contains("NoPr"))   GetPrLepCleanJet =true;
  if(Option.Contains("NoTau")) {ExcTau           =true; Criteria="InclTau";}

  std::vector<snu::KJet> ReturnVec;
  for(int i=0; i<JetColl.size(); i++){
    bool HasEWLep=HasEWLepInJet(JetColl.at(i), TruthColl, Criteria);
    if( GetPrLepCleanJet && (!HasEWLep) ) ReturnVec.push_back(JetColl.at(i));
  }

  return ReturnVec;
}

std::vector<snu::KJet> AnalyzerCore::SkimJetColl(std::vector<snu::KJet>& JetColl, std::vector<snu::KElectron>& EleColl, std::vector<snu::KMuon>& MuColl, TString Option){

  std::vector<snu::KJet> ReturnVec;
  for(int i=0; i<JetColl.size(); i++){

    bool NearLep=false;
    for(int j=0; j<EleColl.size(); j++){
      if(JetColl.at(i).DeltaR(EleColl.at(j))<0.4) NearLep=true;
    }

    for(int j=0; j<MuColl.size(); j++){
      if(JetColl.at(i).DeltaR(MuColl.at(j))<0.4) NearLep=true;
    }

    if(!NearLep) ReturnVec.push_back(JetColl.at(i));
  }

  return ReturnVec;
}



int AnalyzerCore::GetNearJetIdx(snu::KElectron Ele, std::vector<snu::KJet> JetColl){

   int NearJetIdx=-1;
   float dRmin=999., dRmax=0.4;
   for(int i=0; i<JetColl.size(); i++){
     float dR=Ele.DeltaR(JetColl.at(i));
     if(dR>dRmax) continue;
     if(dR<dRmin){ dRmin=dR; NearJetIdx=i;}
   }
   
   return NearJetIdx;
}



//template <class Ptl> int AnalyzerCore::GetDaughterCandIdx(std::vector<Ptl> PtlColl, TString MotherPtl, float WindowWidth, TString Option){
int AnalyzerCore::GetDaughterCandIdx(std::vector<snu::KMuon> PtlColl, TString MotherPtl, float WindowWidth, TString Option){

  if(MotherPtl=="Z"){
    float mindM=9999.; int IdxLead=-1, IdxSubl=-1;
    for(unsigned int i=0; i<PtlColl.size(); i++){
      for(unsigned int j=i+1; j<PtlColl.size(); j++){
        float Mass=(PtlColl.at(i)+PtlColl.at(j)).M();
        if(PtlColl.at(i).Charge()==PtlColl.at(j).Charge()) continue;
        if(fabs(Mass-91.2)<mindM) {mindM=fabs(Mass-91.2); IdxLead=i; IdxSubl=j;}
      }
    }
    if(mindM<WindowWidth){
      if     (Option.Contains("Lead")) return IdxLead;
      else if(Option.Contains("Subl")) return IdxSubl;
    }
    else if(WindowWidth<0 || Option.Contains("NoLimit")){
      if     (Option.Contains("Lead")) return IdxLead;
      else if(Option.Contains("Subl")) return IdxSubl;
    }
  }

  return -1;
};

int AnalyzerCore::GetDaughterCandIdx(std::vector<snu::KJet> PtlColl, TString MotherPtl, float WindowWidth, TString Option){

  float MotherMass=-1., mindM=9999.; int IdxLead=-1, IdxSubl=-1;

  if     (MotherPtl=="Z") MotherMass=91.2; 
  else if(MotherPtl=="W") MotherMass=80.4;
  else return -1;

  for(unsigned int i=0; i<PtlColl.size(); i++){
    for(unsigned int j=i+1; j<PtlColl.size(); j++){
      float Mass=(PtlColl.at(i)+PtlColl.at(j)).M();
      if(fabs(Mass-MotherMass)<mindM) {mindM=fabs(Mass-MotherMass); IdxLead=i; IdxSubl=j;}
    }
  }
  if(mindM<WindowWidth){
    if     (Option.Contains("Lead")) return IdxLead;
    else if(Option.Contains("Subl")) return IdxSubl;
  }
  else if(WindowWidth<0 || Option.Contains("NoLimit")){
    if     (Option.Contains("Lead")) return IdxLead;
    else if(Option.Contains("Subl")) return IdxSubl;
  }

  return -1;
};



double AnalyzerCore::TopPTReweight(std::vector<snu::KTruth> TruthColl){
  //Reference: https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting
  //This is coded for H+>WA CR, but this can be generally used for any kinds of decay of SM top pair events.
  //Caution : This can be used for SM ttbar samples only. Though top pt disagreement is obseved in several generators,
  //          but the SF are derived from Powheg+Pythia sample(Run2). So it is safe to use this only for Powheg ttbar.
  //          And this MUST NOT be applied to single top sample and tops produced from BSM mechanism
  //Current normalisation factor version : v8-0-6

  double weight=1.;

  for(std::vector<snu::KTruth>::iterator it_truth= TruthColl.begin(); it_truth!=TruthColl.end(); it_truth++){
    if(fabs(it_truth->PdgId())==6 && fabs(it_truth->GenStatus())<30 && fabs(it_truth->GenStatus())>20){
      weight*=exp(0.0615-0.0005*it_truth->Pt());
    }
  }
  return sqrt(weight)*0.999777;

}

float AnalyzerCore::GenFilterEfficiency(TString SampleName){

  if( isData ) return 1.;
  else if( !(SampleName.Contains("3mu") || SampleName.Contains("1e2mu")) ) return 1.;

  if     (SampleName.Contains("MHc90_MZp2" )) return 0.309;
  else if(SampleName.Contains("MHc100_MZp2")) return 0.451;
  else if(SampleName.Contains("MHc110_MZp2")) return 0.557;
  else if(SampleName.Contains("MHc120_MZp2")) return 0.614;
  else if(SampleName.Contains("MHc130_MZp2")) return 0.642;
  else if(SampleName.Contains("MHc140_MZp2")) return 0.656;
  else if(SampleName.Contains("MHc150_MZp2")) return 0.628;
  else if(SampleName.Contains("MHc160_MZp2")) return 0.510;

  else if(SampleName.Contains("MHc90_MZp5" )) return 0.275;
  else if(SampleName.Contains("MHc100_MZp5")) return 0.445;
  else if(SampleName.Contains("MHc110_MZp5")) return 0.546;
  else if(SampleName.Contains("MHc120_MZp5")) return 0.603;
  else if(SampleName.Contains("MHc130_MZp5")) return 0.637;
  else if(SampleName.Contains("MHc140_MZp5")) return 0.650;
  else if(SampleName.Contains("MHc150_MZp5")) return 0.631;
  else if(SampleName.Contains("MHc160_MZp5")) return 0.504;

  else if(SampleName.Contains("MHc90_MZp8" )) return 0.268;
  else if(SampleName.Contains("MHc100_MZp8")) return 0.423;
  else if(SampleName.Contains("MHc110_MZp8")) return 0.533;
  else if(SampleName.Contains("MHc120_MZp8")) return 0.596;
  else if(SampleName.Contains("MHc130_MZp8")) return 0.629;
  else if(SampleName.Contains("MHc140_MZp8")) return 0.642;
  else if(SampleName.Contains("MHc150_MZp8")) return 0.616;
  else if(SampleName.Contains("MHc160_MZp8")) return 0.503;

  else if(SampleName.Contains("MHc90_MA10" )) return 0.373;
  else if(SampleName.Contains("MHc100_MA10")) return 0.380;
  else if(SampleName.Contains("MHc110_MA10")) return 0.457;
  else if(SampleName.Contains("MHc120_MA10")) return 0.501;
  else if(SampleName.Contains("MHc130_MA10")) return 0.531;
  else if(SampleName.Contains("MHc140_MA10")) return 0.550;
  else if(SampleName.Contains("MHc150_MA10")) return 0.535;
  else if(SampleName.Contains("MHc160_MA10")) return 0.436;

  else if(SampleName.Contains("MHc90_MA10" )) return 0.373;
  else if(SampleName.Contains("MHc100_MA10")) return 0.380;
  else if(SampleName.Contains("MHc110_MA10")) return 0.457;
  else if(SampleName.Contains("MHc120_MA10")) return 0.501;
  else if(SampleName.Contains("MHc130_MA10")) return 0.531;
  else if(SampleName.Contains("MHc140_MA10")) return 0.550;
  else if(SampleName.Contains("MHc150_MA10")) return 0.535;
  else if(SampleName.Contains("MHc160_MA10")) return 0.436;

  else if(SampleName.Contains("MHc100_MA15")) return 0.464;
  else if(SampleName.Contains("MHc110_MA15")) return 0.496;
  else if(SampleName.Contains("MHc120_MA15")) return 0.528;
  else if(SampleName.Contains("MHc130_MA15")) return 0.547;
  else if(SampleName.Contains("MHc140_MA15")) return 0.558;
  else if(SampleName.Contains("MHc150_MA15")) return 0.540;
  else if(SampleName.Contains("MHc160_MA15")) return 0.440;

  else if(SampleName.Contains("MHc100_MA20")) return 0.556;
  else if(SampleName.Contains("MHc110_MA20")) return 0.564;
  else if(SampleName.Contains("MHc120_MA20")) return 0.572;
  else if(SampleName.Contains("MHc130_MA20")) return 0.576;
  else if(SampleName.Contains("MHc140_MA20")) return 0.576;
  else if(SampleName.Contains("MHc150_MA20")) return 0.554;
  else if(SampleName.Contains("MHc160_MA20")) return 0.450;

  else if(SampleName.Contains("MHc110_MA25")) return 0.620;
  else if(SampleName.Contains("MHc120_MA25")) return 0.617;
  else if(SampleName.Contains("MHc130_MA25")) return 0.613;
  else if(SampleName.Contains("MHc140_MA25")) return 0.609;
  else if(SampleName.Contains("MHc150_MA25")) return 0.577;
  else if(SampleName.Contains("MHc160_MA25")) return 0.459;

  else if(SampleName.Contains("MHc110_MA30")) return 0.648;
  else if(SampleName.Contains("MHc120_MA30")) return 0.652;
  else if(SampleName.Contains("MHc130_MA30")) return 0.644;
  else if(SampleName.Contains("MHc140_MA30")) return 0.632;
  else if(SampleName.Contains("MHc150_MA30")) return 0.594;
  else if(SampleName.Contains("MHc160_MA30")) return 0.471;

  else if(SampleName.Contains("MHc120_MA35")) return 0.672;
  else if(SampleName.Contains("MHc130_MA35")) return 0.665;
  else if(SampleName.Contains("MHc140_MA35")) return 0.649;
  else if(SampleName.Contains("MHc150_MA35")) return 0.612;
  else if(SampleName.Contains("MHc160_MA35")) return 0.485;


  return 1.;

}


float AnalyzerCore::GetHiggsMass(TString SampleName, TString Option){

  bool AMass=Option.Contains("A");
  bool HcMass=Option.Contains("H+");
  if(AMass && HcMass) return -1.;
  if(isData) return -1.;

  float HiggsMass=0.;
  if(AMass){
    if     (SampleName.Contains("MA15")) HiggsMass= 15.;
    else if(SampleName.Contains("MA20")) HiggsMass= 20.;
    else if(SampleName.Contains("MA25")) HiggsMass= 25.;
    else if(SampleName.Contains("MA30")) HiggsMass= 30.;
    else if(SampleName.Contains("MA35")) HiggsMass= 35.;
  }
  else if(HcMass){
    if     (SampleName.Contains("MHc100")) HiggsMass= 100.;
    else if(SampleName.Contains("MHc110")) HiggsMass= 110.;
    else if(SampleName.Contains("MHc120")) HiggsMass= 120.;
    else if(SampleName.Contains("MHc130")) HiggsMass= 130.;
    else if(SampleName.Contains("MHc140")) HiggsMass= 140.;
    else if(SampleName.Contains("MHc150")) HiggsMass= 150.;
    else if(SampleName.Contains("MHc160")) HiggsMass= 160.;
  }

  return HiggsMass;
}


float AnalyzerCore::SignalNorm(TString SampleName, float Xsec){

  //Normalise xsec(tt)*[2*Br(t>bH+)*Br(H+>AW)*Br(A>mumu)*Br(t->bW)] to required value.

  if( isData ) return 1.;
  if( !(SampleName.Contains("TTToHcToWA") || SampleName.Contains("TTToHcToWZp")) ) return 1.;

  float weight=Xsec/20.;
  //20 is current default xsec value for all signal samples
  //At Xsec, we have N*2*[Xsec/(20*2/14%)] =N*Xsec/20*14%

  //Branching fraction weights to each channels
  if     (SampleName.Contains("1e2mu"))  weight*=0.1464;
  else if(SampleName.Contains("3mu"))    weight*=0.1464;
  else if(SampleName.Contains("1ta2mu")) weight*=0.1464;
  else if(SampleName.Contains("2l2mu"))  weight*=0.1061;

  return weight;

}



bool AnalyzerCore::PassIDCriteria(snu::KElectron Ele, TString ID, TString Option){

  bool PassID=true;
  if(ID=="POGCBTIP"){
    if(Ele.SNUID()<1000) PassID=false;
    if( !((fabs(Ele.Eta())<1.447 && Ele.PFRelIso(0.3)< 0.0588) || (fabs(Ele.Eta())>=1.447 && Ele.PFRelIso(0.3)<0.0571)) ) PassID=false;
    if( !(fabs(Ele.dxy())<0.05)  ) PassID=false;
    if( !(fabs(Ele.dz())<0.1)    ) PassID=false;
    if( !(fabs(Ele.dxySig())<3.) ) PassID=false;
  }
  else if(ID=="POGMVAMIP"){
    if( !(Ele.PassNotrigMVAMedium()) ) PassID=false;
    if( !(Ele.IsTrigMVAValid())      ) PassID=false;
    if( !(Ele.PassesConvVeto())      ) PassID=false;
    if( !(Ele.PFRelIso(0.3)<0.1)     ) PassID=false;
    if( !(fabs(Ele.dxy())<0.025)     ) PassID=false;
    if( !(fabs(Ele.dz())<0.05)       ) PassID=false;
    if( !(fabs(Ele.dxySig())<4.)     ) PassID=false;
  }
  else if(ID=="POGMVATIP"){
    if( !(Ele.PassNotrigMVATight())  ) PassID=false;
    if( !(Ele.IsTrigMVAValid())      ) PassID=false;
    if( !(Ele.PassesConvVeto())      ) PassID=false;
    if( !(Ele.PFRelIso(0.3)<0.1)     ) PassID=false;
    if( !(fabs(Ele.dxy())<0.025)     ) PassID=false;
    if( !(fabs(Ele.dz())<0.05)       ) PassID=false;
    if( !(fabs(Ele.dxySig())<4.)     ) PassID=false;
  }
  else if(ID=="POGWP90Isop1IPp025p05sig4"){
    if( !(Ele.PassNotrigMVAMedium()) ) PassID=false;
    if( !(Ele.IsTrigMVAValid())      ) PassID=false;
    if( !(Ele.PassesConvVeto())      ) PassID=false;
    if( !(Ele.PFRelIso(0.3)<0.1)     ) PassID=false;
    if( !(fabs(Ele.dxy())<0.025)     ) PassID=false;
    if( !(fabs(Ele.dz())<0.05)       ) PassID=false;
    if( !(fabs(Ele.dxySig())<4.)     ) PassID=false;
  }
  else if(ID=="POGWP90Isop08IPp025p05sig4"){
    if( !(Ele.PassNotrigMVAMedium()) ) PassID=false;
    if( !(Ele.IsTrigMVAValid())      ) PassID=false;
    if( !(Ele.PassesConvVeto())      ) PassID=false;
    if( !(Ele.PFRelIso(0.3)<0.08)    ) PassID=false;
    if( !(fabs(Ele.dxy())<0.025)     ) PassID=false;
    if( !(fabs(Ele.dz())<0.05)       ) PassID=false;
    if( !(fabs(Ele.dxySig())<4.)     ) PassID=false;
  }
  else if(ID=="POGWP90Isop06IPp025p05sig4"){
    if( !(Ele.PassNotrigMVAMedium()) ) PassID=false;
    if( !(Ele.IsTrigMVAValid())      ) PassID=false;
    if( !(Ele.PassesConvVeto())      ) PassID=false;
    if( !(Ele.PFRelIso(0.3)<0.06)    ) PassID=false;
    if( !(fabs(Ele.dxy())<0.025)     ) PassID=false;
    if( !(fabs(Ele.dz())<0.05)       ) PassID=false;
    if( !(fabs(Ele.dxySig())<4.)     ) PassID=false;
  }
  else if(ID=="LMVA06Isop4IPp025p05sig4"){
    if( !(Ele.IsTrigMVAValid())      ) PassID=false;
    if( !(Ele.PassesConvVeto())      ) PassID=false;
    if( !(Ele.PFRelIso(0.3)<0.4)     ) PassID=false;
    if( !(fabs(Ele.dxy())<0.025)     ) PassID=false;
    if( !(fabs(Ele.dz())<0.05)       ) PassID=false;
    if( !(fabs(Ele.dxySig())<4.)     ) PassID=false;

    if     ( fabs(Ele.Eta())<0.8   ){ if(Ele.MVA()<-0.92) PassID=false; }
    else if( fabs(Ele.Eta())<1.479 ){ if(Ele.MVA()<-0.85) PassID=false; }
    else                            { if(Ele.MVA()<-0.76) PassID=false; }
  }
  else if(ID=="LNoMVANoIsoIPp05p1sig4"){
    if( !(Ele.IsTrigMVAValid())      ) PassID=false;
    if( !(Ele.PassesConvVeto())      ) PassID=false;
    if( !(fabs(Ele.dxy())<0.05)     ) PassID=false;
    if( !(fabs(Ele.dz())<0.1)       ) PassID=false;
    if( !(fabs(Ele.dxySig())<4.)     ) PassID=false;
  }
  else if(ID=="POGMVAMIso"){
    if( !(Ele.PassNotrigMVAMedium()) ) PassID=false;
    if( !(Ele.IsTrigMVAValid())      ) PassID=false;
    if( !(Ele.PassesConvVeto())      ) PassID=false;
    if( !(Ele.PFRelIso(0.3)<0.1)     ) PassID=false;
  }
  else if(ID=="POGMVAM"){
    if( !(Ele.PassNotrigMVAMedium()) ) PassID=false;
    if( !(Ele.IsTrigMVAValid())      ) PassID=false;
    if( !(Ele.PassesConvVeto())      ) PassID=false;
  }
  else if(ID=="POGWP90Isop06IPp025p1sig4"){
    if( !(Ele.PassNotrigMVAMedium()) ) PassID=false;
    if( !(Ele.IsTrigMVAValid())      ) PassID=false;
    if( !(Ele.PassesConvVeto())      ) PassID=false;
    if( !(Ele.PFRelIso(0.3)<0.06)    ) PassID=false;
    if( !(fabs(Ele.dxy())<0.025)     ) PassID=false;
    if( !(fabs(Ele.dz())<0.1)        ) PassID=false;
    if( !(fabs(Ele.dxySig())<4.)     ) PassID=false;
  }
  else if(ID=="HctoWAFakeLoose"){
    if( !(Ele.IsTrigMVAValid())      ) PassID=false;
    if( !(Ele.PassesConvVeto())      ) PassID=false;
    if( !(Ele.PFRelIso(0.3)<0.4)     ) PassID=false;
    if( !(fabs(Ele.dxy())<0.025)     ) PassID=false;
    if( !(fabs(Ele.dz())<0.1   )     ) PassID=false;
    if( !(fabs(Ele.dxySig())<4.)     ) PassID=false;

    if     ( fabs(Ele.Eta())<0.8   ){ if(Ele.MVA()<-0.92) PassID=false; }
    else if( fabs(Ele.Eta())<1.479 ){ if(Ele.MVA()<-0.88) PassID=false; }
    else                            { if(Ele.MVA()<-0.78) PassID=false; }
  }
  else{ cout<<"Error: No Matched ID!"<<endl; return false;}

  return PassID;

}

bool AnalyzerCore::PassIDCriteria(snu::KMuon Mu, TString ID, TString Option){

  bool PassID=true;
  bool Debug=Option.Contains("Debug"); Debug=true; int IdxCase=0;
  bool RochCorr=Option.Contains("Roch");
  float fEta=fabs(Mu.Eta());
  if(ID=="HNTrilepFakeL2"){
    if( !(Mu.IsTight())      ) PassID=false;
    if     ( !RochCorr && !(Mu.RelIso04()<0.4) ) PassID=false;
    else if(  RochCorr && Mu.RochPt()>0 && !(Mu.RelIso04()*Mu.MiniAODPt()/Mu.RochPt()<0.4) ) PassID=false;
    if( !(fabs(Mu.dXY())<0.01) ) PassID=false;
    if( !(fabs(Mu.dZ()) <0.1 ) ) PassID=false;
    if( !(fabs(Mu.dXYSig())<3.)) PassID=false;
  }
  else if(ID=="POGLNoIso"){
    if( !(Mu.IsLoose()) ) PassID=false;
  }
  else if(ID=="POGTIsop4IPp2p1"){
    if     ( !RochCorr && !(Mu.RelIso04()<0.4)     ) PassID=false;
    else if(  RochCorr && !(RochIso(Mu,"0.4")<0.4) ) PassID=false;
    if     ( !(Mu.IsTight())                       ) PassID=false;
    if     ( !(fabs(Mu.dXY())<0.2)      ) PassID=false;
    if     ( !(fabs(Mu.dZ())<0.1)       ) PassID=false;
  }
  else if(ID=="POGTIsop6IPp2p1NoChi"){
    if     ( !RochCorr && !(Mu.RelIso04()<0.6)     ) PassID=false;
    else if(  RochCorr && !(RochIso(Mu,"0.4")<0.6) ) PassID=false;
    if     ( !(Mu.IsPF() && Mu.IsGlobal())         ) PassID=false;
    if     ( !( Mu.validHits()    >0
             && Mu.validPixHits() >0
             && Mu.validStations()>1
             && Mu.ActiveLayer()  >5 )  ) PassID=false;
    if     ( !(fabs(Mu.dXY())<0.2)      ) PassID=false;
    if     ( !(fabs(Mu.dZ())<0.1)       ) PassID=false;
  }
  else if(ID=="POGTIsop4IPp2p1sig4NoChi"){
    if     ( !RochCorr && !(Mu.RelIso04()<0.4)     ) PassID=false;
    else if(  RochCorr && !(RochIso(Mu,"0.4")<0.4) ) PassID=false;
    if     ( !(Mu.IsPF() && Mu.IsGlobal())         ) PassID=false;
    if     ( !( Mu.validHits()    >0
             && Mu.validPixHits() >0
             && Mu.validStations()>1
             && Mu.ActiveLayer()  >5 )  ) PassID=false;
    if     ( !(fabs(Mu.dXY())<0.2)      ) PassID=false;
    if     ( !(fabs(Mu.dZ())<0.1)       ) PassID=false;
    if     ( !(fabs(Mu.dXYSig())<4.)    ) PassID=false;
  }
  else if(ID=="POGTIsop6IPp2p1sig4"){
    if     ( !RochCorr && !(Mu.RelIso04()<0.6)     ) PassID=false;
    else if(  RochCorr && !(RochIso(Mu,"0.4")<0.6) ) PassID=false;
    if     ( !(Mu.IsTight())                       ) PassID=false;
    if     ( !(fabs(Mu.dXY())<0.2)      ) PassID=false;
    if     ( !(fabs(Mu.dZ())<0.1)       ) PassID=false;
    if     ( !(fabs(Mu.dXYSig())<4.)    ) PassID=false;
  }
  else if(ID=="POGTIsop6IPp2p1"){
    if     ( !RochCorr && !(Mu.RelIso04()<0.6)     ) PassID=false;
    else if(  RochCorr && !(RochIso(Mu,"0.4")<0.6) ) PassID=false;
    if     ( !(Mu.IsTight())                       ) PassID=false;
    if     ( !(fabs(Mu.dXY())<0.2)      ) PassID=false;
    if     ( !(fabs(Mu.dZ())<0.1)       ) PassID=false;
  }
  else if(ID=="POGTIsop4IPp01p05sig4Chi4"){
    if     ( !RochCorr && !(Mu.RelIso04()<0.4)     ) PassID=false;
    else if(  RochCorr && !(RochIso(Mu,"0.4")<0.4) ) PassID=false;
    if     ( !(Mu.IsTight())                       ) PassID=false;
    if     ( !(fabs(Mu.dXY())<0.01)     ) PassID=false;
    if     ( !(fabs(Mu.dZ() )<0.05)     ) PassID=false;
    if     ( !(fabs(Mu.dXYSig())<4.)    ) PassID=false;
    if     ( !(fabs(Mu.GlobalChi2())<4.)) PassID=false;
  }
  else if(ID=="POGTIsop6IPp01p05sig4Chi4"){
    if     ( !RochCorr && !(Mu.RelIso04()<0.6)     ) PassID=false;
    else if(  RochCorr && !(RochIso(Mu,"0.4")<0.6) ) PassID=false;
    if     ( !(Mu.IsTight())                       ) PassID=false;
    if     ( !(fabs(Mu.dXY())<0.01)     ) PassID=false;
    if     ( !(fabs(Mu.dZ() )<0.05)     ) PassID=false;
    if     ( !(fabs(Mu.dXYSig())<4.)    ) PassID=false;
    if     ( !(fabs(Mu.GlobalChi2())<4.)) PassID=false;
  }
  else if(ID=="POGLIsop4IPp5p1Chi100"){
    if     ( !RochCorr && !(Mu.RelIso04()<0.4)     ) PassID=false;
    else if(  RochCorr && !(RochIso(Mu,"0.4")<0.4) ) PassID=false;
    if     ( !(Mu.IsLoose())            ) PassID=false;
    if     ( !(fabs(Mu.dXY())<0.5)      ) PassID=false;
    if     ( !(fabs(Mu.dZ())<0.1)       ) PassID=false;
    if     ( !(fabs(Mu.GlobalChi2())<100.)) PassID=false;
  }
  else if(ID=="POGTIsop20IPp01p05sig4Chi4"){
    if     ( !RochCorr && !(Mu.RelIso04()    <0.20) ) PassID=false;
    else if(  RochCorr && !(RochIso(Mu,"0.4")<0.20) ) PassID=false;
    if     ( !(Mu.IsTight())            ) PassID=false;
    if     ( !(fabs(Mu.dXY())<0.01)     ) PassID=false;
    if     ( !(fabs(Mu.dZ() )<0.05)    ) PassID=false;
    if     ( !(fabs(Mu.dXYSig())<4.)    ) PassID=false;
    if     ( !(fabs(Mu.GlobalChi2())<4.)) PassID=false;
  }
  else{ cout<<"Error: No Matched ID!"<<endl; return false;}

  return PassID;
}



bool AnalyzerCore::HasStrangeEleCand(std::vector<snu::KElectron> EleColl){
//Recommended to check whether weird electrons affect each of analysis or not.
//Ref : https://hypernews.cern.ch/HyperNews/CMS/get/physics-announcements/4803.html

  bool HasIt=false;
  for(int i=0; i<EleColl.size(); i++){
    if(CouldBeStrangeEleCand(EleColl.at(i))){ HasIt=true; break; }
  }

  return HasIt;
}

bool AnalyzerCore::CouldBeStrangeEleCand(snu::KElectron Ele){
//Recommended to check whether weird electrons affect each of analysis or not.
//Ref : https://hypernews.cern.ch/HyperNews/CMS/get/physics-announcements/4803.html

  return fabs(Ele.Eta()-Ele.SCEta())>0.2;

}

float AnalyzerCore::RochIso(snu::KMuon Mu, TString ConeSize){

  if     ( !Mu.IsRochesterCorrected() ) return 0.;
  else if(  Mu.RochPt()==0.           ) return 0.;

  float PFRelIsodbeta = ConeSize=="0.4"? Mu.RelIso04():Mu.RelIso03();
  return PFRelIsodbeta*Mu.MiniAODPt()/Mu.RochPt();

}


int AnalyzerCore::GetSigGenPtlIdx(vector<snu::KTruth>& TruthColl, TString PtlName){

  int IdxPtl=-1;

  if(PtlName=="mup_A"){
    for(int i=2; i<TruthColl.size(); i++){
      if(TruthColl.at(i).GenStatus()!=1) continue;
      if(TruthColl.at(i).PdgId()!=-13) continue;
      if(TruthColl.at(FirstNonSelfMotherIdx(i, TruthColl)).PdgId()!=36) continue;
      IdxPtl=i; break;
    }
  }
  else if(PtlName=="mum_A"){
    for(int i=2; i<TruthColl.size(); i++){
      if(TruthColl.at(i).GenStatus()!=1) continue;
      if(TruthColl.at(i).PdgId()!=13) continue;
      if(TruthColl.at(FirstNonSelfMotherIdx(i, TruthColl)).PdgId()!=36) continue;
      IdxPtl=i; break;
    }
  }
  else if(PtlName=="mu_W_hc" || PtlName=="l_W_hc"){
    for(int i=2; i<TruthColl.size(); i++){
      if(TruthColl.at(i).GenStatus()!=1) continue;
      if(TruthColl.at(i).PdgId()!=-13) continue;
      if(TruthColl.at(FirstNonSelfMotherIdx(i, TruthColl)).PdgId()!=24) continue;
      IdxPtl=i; break;
    }
  }
  else if(PtlName=="mu_W_tx" || PtlName=="l_W_tx"){
    for(int i=2; i<TruthColl.size(); i++){
      if(TruthColl.at(i).GenStatus()!=1) continue;
      if(TruthColl.at(i).PdgId()!=13) continue;
      if(TruthColl.at(FirstNonSelfMotherIdx(i, TruthColl)).PdgId()!=-24) continue;
      IdxPtl=i; break;
    }
  }
  else if(PtlName=="e_W_hc" || PtlName=="l_W_hc"){
    for(int i=2; i<TruthColl.size(); i++){
      if(TruthColl.at(i).GenStatus()!=1) continue;
      if(TruthColl.at(i).PdgId()!=-11) continue;
      if(TruthColl.at(FirstNonSelfMotherIdx(i, TruthColl)).PdgId()!=24) continue;
      IdxPtl=i; break;
    }
  }
  else if(PtlName=="e_W_tx" || PtlName=="l_W_tx"){
    for(int i=2; i<TruthColl.size(); i++){
      if(TruthColl.at(i).GenStatus()!=1) continue;
      if(TruthColl.at(i).PdgId()!=11) continue;
      if(TruthColl.at(FirstNonSelfMotherIdx(i, TruthColl)).PdgId()!=-24) continue;
      IdxPtl=i; break;
    }
  }
  else if(PtlName=="mu_ta_W_hc" || PtlName=="l_ta_W_hc"){
    for(int i=2; i<TruthColl.size(); i++){
      if(TruthColl.at(i).GenStatus()!=1) continue;
      if(TruthColl.at(i).PdgId()!=-13) continue;
      int MotherIdx=FirstNonSelfMotherIdx(i, TruthColl);
      if(TruthColl.at(MotherIdx).PdgId()!=-15) continue;
      if(TruthColl.at(FirstNonSelfMotherIdx(MotherIdx, TruthColl)).PdgId()!=24) continue;
      IdxPtl=i; break;
    }
  }
  else if(PtlName=="mu_ta_W_tx" || PtlName=="l_ta_W_tx"){
    for(int i=2; i<TruthColl.size(); i++){
      if(TruthColl.at(i).GenStatus()!=1) continue;
      if(TruthColl.at(i).PdgId()!=13) continue;
      int MotherIdx=FirstNonSelfMotherIdx(i, TruthColl);
      if(TruthColl.at(MotherIdx).PdgId()!=15) continue;
      if(TruthColl.at(FirstNonSelfMotherIdx(MotherIdx, TruthColl)).PdgId()!=-24) continue;
      IdxPtl=i; break;
    }
  }
  else if(PtlName=="e_ta_W_hc" || PtlName=="l_ta_W_hc"){
    for(int i=2; i<TruthColl.size(); i++){
      if(TruthColl.at(i).GenStatus()!=1) continue;
      if(TruthColl.at(i).PdgId()!=-11) continue;
      int MotherIdx=FirstNonSelfMotherIdx(i, TruthColl);
      if(TruthColl.at(MotherIdx).PdgId()!=-15) continue;
      if(TruthColl.at(FirstNonSelfMotherIdx(MotherIdx, TruthColl)).PdgId()!=24) continue;
      IdxPtl=i; break;
    }
  }
  else if(PtlName=="e_ta_W_tx" || PtlName=="l_ta_W_tx"){
    for(int i=2; i<TruthColl.size(); i++){
      if(TruthColl.at(i).GenStatus()!=1) continue;
      if(TruthColl.at(i).PdgId()!=11) continue;
      int MotherIdx=FirstNonSelfMotherIdx(i, TruthColl);
      if(TruthColl.at(MotherIdx).PdgId()!=15) continue;
      if(TruthColl.at(FirstNonSelfMotherIdx(MotherIdx, TruthColl)).PdgId()!=-24) continue;
      IdxPtl=i; break;
    }
  }

  return IdxPtl;

}


void AnalyzerCore::SetMuonResCorrection(vector<snu::KMuon>& MuColl, TString Option){

  bool SystVar=false;
  float SystDir=0.;
  if(Option.Contains("Syst")){
    if     (Option.Contains("Up"))  {SystVar=true; SystDir= 1.;}
    else if(Option.Contains("Down")){SystVar=true; SystDir=-1.;}
  }

  for(int i=0; i<(int) MuColl.size(); i++){
    double dRelPt = mcdata_correction->GetRochesterMomentumWidth(MuColl.at(i));
    MuColl.at(i) *= (1.+SystDir*dRelPt);
  }

  return;
}


float AnalyzerCore::GetXsecUncertainty(TString SampleName, TString Option){

  int ReturnType=1; //1:Total 2:Q2 3:PDF 4:as
  if     (Option.Contains("Q2") ) ReturnType=2;
  else if(Option.Contains("PDF")) ReturnType=3;
  else if(Option.Contains("as") ) ReturnType=4;

  float Uncertainty=0.;
  if     (SampleName.Contains("WZTo3LNu_powheg")){
    if(ReturnType==1) Uncertainty=0.041;//NNLO
  }
  else if(SampleName.Contains("ZZTo4L_powheg")){
    if(ReturnType==1) Uncertainty=0.08;//NNLO
  }
  else if(SampleName.Contains("ttW")){
    if(ReturnType==1) Uncertainty=0.133;
  }
  else if(SampleName.Contains("ttZ")){
    if(ReturnType==1) Uncertainty=0.12;
  }
  else if(SampleName.Contains("ttH_nonbb")){
    if(ReturnType==1) Uncertainty=0.128;
  }
  else if(SampleName.Contains("tZq")){
    if(ReturnType==1) Uncertainty=0.5;
  }

  return Uncertainty;
}
