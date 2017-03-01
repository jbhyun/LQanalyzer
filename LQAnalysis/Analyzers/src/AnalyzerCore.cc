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


AnalyzerCore::AnalyzerCore() : LQCycleBase(), n_cutflowcuts(0), MCweight(-999.),reset_lumi_mask(false),changed_target_lumi(false), k_reset_period(false), a_mcperiod(-1) {

  k_debugmode=false;
  
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

  mcdata_correction = new MCDataCorrections();
  
  string lqdir = getenv("LQANALYZER_DIR");

  /// DataDrivenBackgrounds class is used to get data driven fake + charge flip backgrounds
  m_datadriven_bkg = new DataDrivenBackgrounds();

  /// Currently only have csvv2 or cMVAv2 btaggers: In HN we use csvv2 
  /// List of taggers
  std::vector<TString> vtaggers;
  vtaggers.push_back("CSVv2Moriond17_2017_1_26_BtoF");
  vtaggers.push_back("CSVv2Moriond17_2017_1_26_GtoH");
  vtaggers.push_back("cMVAv2Moriond17_2017_1_26_BtoF");
  vtaggers.push_back("cMVAv2Moriond17_2017_1_26_GtoH");
  /// Will add DeepCSV in 805

  if( getenv("CATDEBUG") == "True") k_debugmode=true;
  cout << "Setting up 2016 selection " << endl;
  SetupSelectionMuon(lqdir + "/CATConfig/SelectionConfig/muons.sel");
  SetupSelectionMuon(lqdir + "/CATConfig/SelectionConfig/user_muons.sel");
  
  SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/electrons.sel");
  SetupSelectionElectron(lqdir + "/CATConfig/SelectionConfig/user_electrons.sel");
  
  SetupSelectionJet(lqdir + "/CATConfig/SelectionConfig/jets.sel");
  SetupSelectionJet(lqdir + "/CATConfig/SelectionConfig/user_jets.sel");
  
  if(k_debugmode){
    for( map<TString,vector<pair<TString,float> > >::iterator it=  selectionIDMapfMuon.begin() ; it !=  selectionIDMapfMuon.end(); it++){
      cout << it->first << endl;
      for (unsigned int i=0 ; i < it->second.size(); i++){
	cout << it->second.at(i).first << " " << it->second.at(i).second << endl;
      }
    }
  }
  Message("SetupSelection DONE", DEBUG); 

  
  // List of working points
  std::vector<TString> v_wps;
  v_wps.push_back("Loose");
  v_wps.push_back("Medium");
  v_wps.push_back("Tight");
  MapBTagSF = SetupBTagger(vtaggers,v_wps);
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
  

  //==== HN Gen Matching Class
  m_HNgenmatch = new HNGenMatching();
  
  cout <<  "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;

}

void AnalyzerCore::SetupLuminosityMap(bool initialsetup, TString forceperiod){

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
    else if(singleperiod=="D")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_277772_278808.txt";
    else if(singleperiod=="E")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_276831_277420.txt";
    else if(singleperiod=="F")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_276315_276811.txt";
    else if(singleperiod=="G")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_280919_284044.txt";
    else if(singleperiod=="H")lumitriggerpath=lqdir + "/data/Luminosity/"+getenv("yeartag")+"/triggers_catversion_" + getenv("CATVERSION")+"_278820_280385.txt";
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



float AnalyzerCore::CorrectedMETRochester(TString muid_formet, bool update_met){

  /// function returns corrected met + can be used to set event met to corrected met

  float met_x =eventbase->GetEvent().PFMETx(); 
  float met_y =eventbase->GetEvent().PFMETy();
  std::vector<snu::KMuon> muall = GetMuons(muid_formet);
  
  float px_orig(0.), py_orig(0.),px_corrected(0.), py_corrected(0.);
  for(unsigned int im=0; im < muall.size() ; im++){
      
      px_orig+= muall.at(im).MiniAODPt()*TMath::Cos(muall.at(im).Phi());
      py_orig+= muall.at(im).MiniAODPt()*TMath::Sin(muall.at(im).Phi());
      px_corrected += muall.at(im).Px();
      py_corrected += muall.at(im).Py();
      
  }
  met_x = met_x + px_orig - px_corrected;	
  met_y = met_y + py_orig - py_corrected;	
  
  if(update_met){
    if(!eventbase->GetEvent().PropagatedRochesterToMET()){
      snu::KEvent tempev = eventbase->GetEvent();
      tempev.SetMET(snu::KEvent::pfmet,  sqrt(met_x*met_x + met_y*met_y), eventbase->GetEvent().METPhi(), eventbase->GetEvent().SumET());
      tempev.SetPFMETx(met_x);
      tempev.SetPFMETy(met_y);
      tempev.SetPropagatedRochesterToMET(true);
      eventbase->SetEventBase(tempev);
    }
  }
  return sqrt(met_x*met_x + met_y*met_y);
}   





float AnalyzerCore::CorrectedMETElectron(TString elid_formet, int sys){

  float met_x =eventbase->GetEvent().PFMETx();
  float met_y =eventbase->GetEvent().PFMETy();
  std::vector<snu::KElectron> elall = GetElectrons(elid_formet);

  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int iel=0; iel < elall.size() ; iel++){


    px_orig+= elall.at(iel).Px();
    py_orig+= elall.at(iel).Py();
    if(sys==1){
      px_shifted += elall.at(iel).Px()*elall.at(iel).PtShiftedUp();
    }
    if(sys==-1){
      px_shifted += elall.at(iel).Px()*elall.at(iel).PtShiftedDown();
    }


  }
  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;


  return sqrt(met_x*met_x + met_y*met_y);

}

float AnalyzerCore::CorrectedMETMuon(TString muid_formet, int sys){
  
  float met_x =eventbase->GetEvent().PFMETx();
  float met_y =eventbase->GetEvent().PFMETy();
  std::vector<snu::KMuon> muall = GetMuons(muid_formet);
  
  float px_orig(0.), py_orig(0.),px_shifted(0.), py_shifted(0.);
  for(unsigned int imu=0; imu < muall.size() ; imu++){
    
    px_orig+= muall.at(imu).Px();
    py_orig+= muall.at(imu).Py();
    if(sys==1){
      px_shifted += muall.at(imu).Px()*muall.at(imu).PtShiftedUp();
    }
    if(sys==-1){
      px_shifted += muall.at(imu).Px()*muall.at(imu).PtShiftedDown();
    }  
  }
  met_x = met_x + px_orig - px_shifted;
  met_y = met_y + py_orig - py_shifted;
  
  
  return sqrt(met_x*met_x + met_y*met_y);
  
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

int AnalyzerCore::GetMCPeriod(){
  /// This function returns a period B-H for MC events. 
  /// It uses a random number and retrunds a period based on the luminosity of each period
  /// It assumes the trigger used is unprescaled
  if(isData) return -1;
  if(!k_reset_period) return a_mcperiod;
  k_reset_period=false;
  

  if(k_mcperiod > 0) {
    a_mcperiod = k_mcperiod;
    return a_mcperiod;
  }
    

  //  double r = ((double) rand() / (RAND_MAX));
  double r =gRandom->Rndm(); /// random number between 0 and 1
  
  
  /// values obtained from cattuple googledoc
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

        }
        else if ( x ==1) {is >> idlabel;}
        else {
          is >> tmpf;
          float_jetsel.push_back(make_pair(cutnames.at(x),tmpf));
        }
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

        if (x > 6 && x < 18){
          is >> tmp;
          string_jetsel.push_back(make_pair(cutnames.at(x),tmp) );

        }
        else if ( x ==1) {is >> idlabel;}
        else {
          is >> tmpf;
          float_jetsel.push_back(make_pair(cutnames.at(x),tmpf));
        }
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

	if (x > 10 && x < 16){
	  is >> tmp;
	  string_muonsel.push_back(make_pair(cutnames.at(x),tmp) );

	}
	else if ( x ==1) {is >> idlabel;}
	else if (x > 26 && x < 28){
          is >> tmp;
          string_muonsel.push_back(make_pair(cutnames.at(x),tmp) );
        }
	else {
	  is >> tmpf;
	  float_muonsel.push_back(make_pair(cutnames.at(x),tmpf));
	}
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
          string_elsel.push_back(make_pair(cutnames.at(x),tmp) );
        }
	else  if (x > 26 && x < 30){
          is >> tmp;
          string_elsel.push_back(make_pair(cutnames.at(x),tmp) );
        }
        else if ( x ==1) {is >> idlabel;}
        else {
          is >> tmpf;
          float_elsel.push_back(make_pair(cutnames.at(x),tmpf));
        }
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
	tmpmap[*it + "_" + *it2 + "_lf"]= new BTagSFUtil("incl", it->Data(), it2->Data());
	tmpmap[*it +  "_" + *it2 + "_hf"]= new BTagSFUtil("mujets", it->Data(), it2->Data());
	// tmpmap[*it +  "_" + *it2 + "_hfcomb"]= new BTagSFUtil("comb", it->Data(), it2->Data());                /// SWITCH ON IF USER NEEDS THIS METHOD  
	// tmpmap[*it +  "_" + *it2 + "iterativefit"]= new BTagSFUtil("iterativefit", it->Data(), it2->Data());   /// SWITCH ON IF USER NEEDS THIS METHOD
      }
      if (it->Contains("cMVA")){
	tmpmap[*it + "_" + *it2 + "_lf"]= new BTagSFUtil("incl", it->Data(), it2->Data());
	tmpmap[*it +  "_" + *it2 + "_hf"]= new BTagSFUtil("ttbar", it->Data(), it2->Data());
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
  return GetElectrons(true, true,GetStringID(electronid), ptcut, etacut);
}

std::vector<snu::KElectron> AnalyzerCore::GetElectrons(bool keepcf, bool keepfake, BaseSelection::ID electronid, float ptcut, float etacut){
  return GetElectrons(keepcf, keepfake,GetStringID(electronid), ptcut, etacut);
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
    
    if (muontag.Contains("NONE") && eltag.Contains("NONE"))  eventbase->GetJetSel()->SelectJets(jetColl,jetid, ptcut,etacut);
    else if (muontag.Contains("NONE") || eltag.Contains("NONE")) {    cerr << "cannot choose to remove jets near only one lepton" << endl; exit(EXIT_FAILURE);}
    else eventbase->GetJetSel()->SelectJets(jetColl, GetMuons(muontag), GetElectrons(eltag) ,jetid, ptcut,etacut);
  }

  return jetColl;
  
}


std::vector<snu::KFatJet> AnalyzerCore::GetFatJets(TString fatjetid, float ptcut, float etacut){

  std::vector<snu::KFatJet> fatjetColl;

  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsFatJet.find(fatjetid);
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
    if (muontag.Contains("NONE") && eltag.Contains("NONE"))  eventbase->GetFatJetSel()->SelectFatJets(fatjetColl,fatjetid, ptcut,etacut);
    else if (muontag.Contains("NONE") || eltag.Contains("NONE")) {    cerr << "cannot choose to remove jets near only one lepton" << endl; exit(EXIT_FAILURE);}
    else eventbase->GetFatJetSel()->SelectFatJets(fatjetColl, GetMuons(muontag), GetElectrons(eltag) ,fatjetid, ptcut,etacut);
    
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
  if(it== selectionIDMapsMuon.end()){
    cerr << "Muon ID ["+muid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {
    if (ptcut == -999.)  eventbase->GetMuonSel()->SelectMuons(muonColl,muid);
    else eventbase->GetMuonSel()->SelectMuons(muonColl,muid, ptcut, etacut);
  }
  //if(k_running_nonprompt) eventbase->GetMuonSel()->SelectMuons(muonColl, BaseSelection::MUON_HN_FAKELOOSE, 15., 2.4);
  //else eventbase->GetMuonSel()->SelectMuons(muonColl, BaseSelection::MUON_HN_TIGHT, 15., 2.4);
  
  return  GetTruePrompt(muonColl, keepfakes);
  
}


std::vector<snu::KElectron> AnalyzerCore::GetElectrons(TString elid,float ptcut, float etacut){
  return GetElectrons( true,  true, elid, ptcut, etacut);
}

std::vector<snu::KElectron> AnalyzerCore::GetElectrons(bool keepcf, bool keepfake, TString elid,float ptcut, float etacut){
  
  std::vector<snu::KElectron> electronColl;
  
  if(elid.Contains("NONE")) return electronColl;


  std::map<TString, vector<pair<TString,TString> > >::iterator it = selectionIDMapsElectron.find(elid);
  if(it== selectionIDMapsElectron.end()){
    cerr << "Electron ID ["+elid+"] not found" << endl; exit(EXIT_FAILURE);
  }
  else {
    if (ptcut == -999.)  eventbase->GetElectronSel()->SelectElectrons(electronColl,elid);
    else eventbase->GetElectronSel()->SelectElectrons(electronColl,elid, ptcut, etacut);
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

  if(period < 0) {
    Message("period not set in AnalyzerCore::HasCloseBJet. Will assign mcperiod for you but this may not give correct behaviour", WARNING);
    period=GetMCPeriod();
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
    if(triggername.Contains(mit->first)) return (mit->second / tlumi);
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

  for(map<TString, TH1*>::iterator it = maphist.begin(); it!= maphist.end(); it++){
    delete it->second;
  }
  maphist.clear();

  for(map<TString, TH2*>::iterator it = maphist2D.begin(); it!= maphist2D.end(); it++){
    delete it->second;
  }
  maphist2D.clear();



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


  for(std::map<TString,BTagSFUtil*>::iterator it = MapBTagSF.begin(); it!= MapBTagSF.end(); it++){
    delete it->second;
  }
  MapBTagSF.clear();

  for(map<TString, HNTriLeptonPlots*>::iterator it = mapCLhistHNTriLep.begin(); it != mapCLhistHNTriLep.end(); it++){
    delete it->second;
  }
  mapCLhistHNTriLep.clear();

  //// New class functions for databkg+corrections
  delete mcdata_correction;
  delete m_datadriven_bkg;

  //==== HN Gen Matching
  delete m_HNgenmatch;

}

//###
//###   IMPORTANT BASE FUNCTION: SETS UP EVENT FOR ALL CYCLES
//###

void AnalyzerCore::SetUpEvent(Long64_t entry, float ev_weight) throw( LQError ) {
  
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
  std::vector<snu::KFatJet> skfatjets= GetAllFatJets();
  std::vector<snu::KGenJet> skgenjets=GetAllGenJets();
  
   
  /// np == numberof particles you want to store at truth info. 30 is default unless running nocut sktree OR signal
  int np =  AssignnNumberOfTruth();
  
  LQEvent lqevent(GetAllMuons(), GetAllElectrons(), GetAllPhotons(), skjets,skfatjets, skgenjets,GetTruthParticles(np), triggerinfo,eventinfo);
  
  //  eventbase is master class to use in analysis 
  //
  
  eventbase = new EventBase(lqevent);

  eventbase->GetElectronSel()->SetIDSMap(selectionIDMapsElectron);
  eventbase->GetElectronSel()->SetIDFMap(selectionIDMapfElectron);
  eventbase->GetMuonSel()->SetIDSMap(selectionIDMapsMuon);
  eventbase->GetMuonSel()->SetIDFMap(selectionIDMapfMuon);
  eventbase->GetJetSel()->SetIDSMap(selectionIDMapsJet);
  eventbase->GetJetSel()->SetIDFMap(selectionIDMapfJet);
  eventbase->GetFatJetSel()->SetIDSMap(selectionIDMapsFatJet);
  eventbase->GetFatJetSel()->SetIDFMap(selectionIDMapfFatJet);

  m_datadriven_bkg->SetEventBase(eventbase);
  
  if(!k_isdata){
    if(!changed_target_lumi){
      changed_target_lumi=true;
    }
  }
  

  /// Setup correction class
  k_reset_period=true;
  
  mcdata_correction->SetMCPeriod(GetMCPeriod());
  mcdata_correction->SetIsData(isData);

  
  
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

int AnalyzerCore::AssignnNumberOfTruth(){
  int np = 1000;
  if(k_classname.Contains("SKTreeMaker")) np = 1000;
  if(k_classname.Contains("SKTreeMakerDiLep")) np = 0;
  if(k_classname.Contains("SKTreeMakerTriLep")) np = 0;

  if(k_classname.Contains("SKTreeMaker")){
    if(k_sample_name.Contains("QCD") && !k_sample_name.Contains("mad")) np = 0;
  }

  /// List of signal samples
  /// G.Yu needs to add signal here
  
  if(IsSignal()) np = 1000;
  
  return np;
}



bool AnalyzerCore::IsSignal(){
  
  if(k_sample_name.Contains("Majornana")) return true;
  if(k_sample_name.Contains("HN")) return true;
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
void AnalyzerCore::TruthPrintOut(){
  if(isData) return;
  m_logger << INFO<< "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  cout << "Particle Index |  PdgId  | GenStatus   | Mother PdgId |  Part_Eta | Part_Pt | Part_Phi | Mother Index |   " << endl;
  for(unsigned int ig=0; ig < eventbase->GetTruth().size(); ig++){
    
    if(eventbase->GetTruth().at(ig).IndexMother() <= 0)continue;
    if(eventbase->GetTruth().at(ig).IndexMother() >= int(eventbase->GetTruth().size()))continue;
    if (eventbase->GetTruth().at(ig).PdgId() == 2212)  cout << ig << " | " << eventbase->GetTruth().at(ig).PdgId() << "  |               |         |        |         |        |         |" << endl;

    cout << ig << " |  " <<  eventbase->GetTruth().at(ig).PdgId() << " |  " << eventbase->GetTruth().at(ig).GenStatus() << " |  " << eventbase->GetTruth().at(eventbase->GetTruth().at(ig).IndexMother()).PdgId()<< " |   " << eventbase->GetTruth().at(ig).Eta() << " | " << eventbase->GetTruth().at(ig).Pt() << " | " << eventbase->GetTruth().at(ig).Phi() << " |   " << eventbase->GetTruth().at(ig).IndexMother()  << endl;
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
  
  delete eventbase;                                                                                                            

}
  

void AnalyzerCore::ListTriggersAvailable(){
  cout << "Set of triggers you can use are: " << endl;
  for(unsigned int i=0; i < eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames().size(); i++){
    cout << eventbase->GetTrigger().GetHLTInsideDatasetTriggerNames().at(i)<< " has prescale " << eventbase->GetTrigger().GetHLTInsideDatasetTriggerPrescales().at(i)<< endl;
  }
  return;
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
  cout << "index" << '\t' << "pdgid" << '\t' << "mother" << '\t' << "mother pid" << endl;
  for(int i=2; i<truthColl.size(); i++){
    cout << i << '\t' << truthColl.at(i).PdgId() << '\t' << truthColl.at(i).IndexMother() << '\t' << truthColl.at( truthColl.at(i).IndexMother() ).PdgId() << endl;
  }

}


void AnalyzerCore::Message(TString message, LQMsgType type){
  m_logger <<  type << message << LQLogger::endmsg;
}



////###############################################################################################                                                                           
 /// @@@@@@@@@@@@@@@@@@@@@@@@@ HIST   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                                                                             

void AnalyzerCore::MakeCleverHistograms(histtype type, TString clhistname ){
  
  //// ELECTRON PLOTs                                                                                          
  if(type==elhist) mapCLhistEl[clhistname] = new ElectronPlots(clhistname);
  //// BaseSelection::MUON PLOTs                                                                                              
  if(type==muhist) mapCLhistMu[clhistname] = new MuonPlots(clhistname);
  /// JET PLOTs                                                                                                
  if(type==jethist) mapCLhistJet[clhistname] = new JetPlots(clhistname);
  /// Signal plots                                                                                             
  if(type==sighist_ee)  mapCLhistSigEE[clhistname] = new SignalPlotsEE(clhistname);
  if(type==sighist_mm)  mapCLhistSigMM[clhistname] = new SignalPlotsMM(clhistname);
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

void AnalyzerCore::MakeHistograms2D(TString hname, int nbinsx, float xmin, float xmax, int nbinsy, float ymin, float ymax, TString label) {

  maphist2D[hname] =  new TH2D(hname.Data(),hname.Data(),nbinsx,xmin,xmax, nbinsy,ymin,ymax);
  maphist2D[hname]->GetXaxis()->SetTitle(label);
}

void AnalyzerCore::MakeHistograms2D(TString hname, int nbinsx,  float xbins[], int nbinsy,  float ybins[], TString label) {

  maphist2D[hname] =  new TH2D(hname.Data(),hname.Data(),nbinsx , xbins, nbinsy,ybins);
  maphist2D[hname]->GetXaxis()->SetTitle(label);
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

void AnalyzerCore::FillHist(TString histname, float value1, float value2, float w, float xmin, float xmax, int nbinsx, float ymin, float ymax, int nbinsy , TString label){

  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist2D(histname)) GetHist2D(histname)->Fill(value1,value2, w);
  else{
    if (nbinsx < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms2D(histname, nbinsx, xmin, xmax,nbinsy, ymin, ymax , label);
    if(GetHist2D(histname)) GetHist2D(histname)->GetXaxis()->SetTitle(label);
    if(GetHist2D(histname)) GetHist2D(histname)->Fill(value1,value2, w);
  }

}

void AnalyzerCore::FillHist(TString histname, float valuex, float valuey, float w, float xbins[], int nxbins, float ybins[], int nybins , TString label){
  m_logger << DEBUG << "FillHist : " << histname << LQLogger::endmsg;
  if(GetHist2D(histname)) GetHist2D(histname)->Fill(valuex,valuey, w);

  else{
    if (nxbins < 0) {
      m_logger << ERROR << histname << " was NOT found. Nbins was not set also... please configure histogram maker correctly" << LQLogger::endmsg;
      exit(0);
    }
    m_logger << DEBUG << "Making the histogram" << LQLogger::endmsg;
    MakeHistograms2D(histname, nxbins, xbins, nybins, ybins , label);
    if(GetHist2D(histname)) GetHist2D(histname)->GetXaxis()->SetTitle(label);
    
    if(GetHist2D(histname)) GetHist2D(histname)->Fill(valuex, valuey, w);
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
    if(sigpit_ee !=mapCLhistSigEE.end()) sigpit_ee->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistSigEE[hist] = new SignalPlotsEE(hist);
      sigpit_ee = mapCLhistSigEE.find(hist);
      sigpit_ee->second->Fill(ev, muons, electrons, jets,w);
    }
  }
  else if(type==sighist_mm){

    map<TString, SignalPlotsMM*>::iterator sigpit_mm = mapCLhistSigMM.find(hist);
    if(sigpit_mm !=mapCLhistSigMM.end()) sigpit_mm->second->Fill(ev, muons, electrons, jets,w);
    else {
      mapCLhistSigMM[hist] = new SignalPlotsMM(hist);
      sigpit_mm = mapCLhistSigMM.find(hist);
      sigpit_mm->second->Fill(ev, muons, electrons, jets,w);
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
    
    
    
    if(mapit->first.Contains("Basic")){
      if(!m_outputFile->GetDirectory( "Basic" )){
	Dir = m_outputFile->mkdir("Basic");
	m_outputFile->cd( Dir->GetName() );
      }
      else  m_outputFile->cd("Basic");
      mapit->second->Write();
      m_outputFile->cd();
    }

    
    
    else {
      mapit->second->Write();
    }
  }
  
  for(map<TString, TH2*>::iterator mapit = maphist2D.begin(); mapit != maphist2D.end(); mapit++){
    mapit->second->Write();
  }

  //==== HN Gen Matching
  m_HNgenmatch->WriteHNGenHists();

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
  
  if(electrons.size() > 2){
    int p_charge=0;
    int m_charge=0;
    for(unsigned int iel = 0 ; iel < electrons.size() ; iel++){
      if(electrons.at(iel).Charge() < 0 ) m_charge++;
      if(electrons.at(iel).Charge() > 0 ) p_charge++;
    }
    if(p_charge > 1) return true;
    if(m_charge > 1) return true;
  }
  if(electrons.size()!=2) return false;


  if(!runningcf){
    if(electrons.at(0).Charge() == electrons.at(1).Charge()) return true;
  }
  else     if(electrons.at(0).Charge() != electrons.at(1).Charge()) return true;

  return false;
}


int AnalyzerCore::NBJet(std::vector<snu::KJet> jets,  KJet::Tagger tag, KJet::WORKING_POINT wp, int period){

  int nbjet=0;

  if(period < 0) {
    Message("period not set in AnalyzerCore::NBJet. Will assign mcperiod for you but this may not give correct behaviour", WARNING);
    period=GetMCPeriod();
  }

  TString btag_key_lf("") , btag_key_hf("");
  TString wp_string="";
  if(wp == snu::KJet::Loose)wp_string = "Loose";
  if(wp == snu::KJet::Medium)wp_string = "Medium";
  if(wp == snu::KJet::Tight)wp_string = "Tight";

  TString tag_string="";

  if(period < 6){
    if(tag== snu::KJet::CSVv2) tag_string ="CSVv2Moriond17_2017_1_26_BtoF";
  }
  else
    if(tag== snu::KJet::CSVv2) tag_string ="CSVv2Moriond17_2017_1_26_GtoH";

  if(period < 6){
    if(tag== snu::KJet::cMVAv2) tag_string ="cMVAv2Moriond17_2017_1_26_BtoF";
  }
  else
    if(tag== snu::KJet::cMVAv2) tag_string ="cMVAv2Moriond17_2017_1_26_GtoH";

   
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

      if (it_lf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  -999999, jets.at(ij).Pt(), jets.at(ij).Eta()))
	isBtag=true;
    }
    else if (jets.at(ij).HadronFlavour() > 1){
      if (it_hf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  jets.at(ij).HadronFlavour(),jets.at(ij).Pt(), jets.at(ij).Eta()))
        isBtag=true;
    }
    else{
      if (it_lf->second->IsTagged(jets.at(ij).BJetTaggerValue(tag),  jets.at(ij).HadronFlavour(),jets.at(ij).Pt(), jets.at(ij).Eta()))
	isBtag=true;
    }
    
    if(isBtag )nbjet++;
  }
  return nbjet;
}


bool AnalyzerCore::IsBTagged(snu::KJet jet,  KJet::Tagger tag, KJet::WORKING_POINT wp, int mcperiod){

  if(mcperiod < 0) {
    Message("mcperiod not set in AnalyzerCore::IsBTagged. Will assign mcperiod for you but this may not give correct behaviour", WARNING);      
    mcperiod=GetMCPeriod();
  }

  TString btag_key_lf("") , btag_key_hf("");
  TString wp_string="";
  if(wp == snu::KJet::Loose)wp_string = "Loose";
  if(wp == snu::KJet::Medium)wp_string = "Medium";
  if(wp == snu::KJet::Tight)wp_string = "Tight";

  TString tag_string="";
  if(mcperiod < 6){
    if(tag== snu::KJet::CSVv2) tag_string ="CSVv2Moriond17_2017_1_26_BtoF";
  }
  else
    if(tag== snu::KJet::CSVv2) tag_string ="CSVv2Moriond17_2017_1_26_GtoH";

  if(mcperiod < 6){
    if(tag== snu::KJet::cMVAv2) tag_string ="cMVAv2Moriond17_2017_1_26_BtoF";
  }
  else
    if(tag== snu::KJet::cMVAv2) tag_string ="cMVAv2Moriond17_2017_1_26_GtoH";




  btag_key_lf = tag_string+"_"+wp_string+"_lf";
  btag_key_hf = tag_string+"_"+wp_string+"_hf";
  std::map<TString,BTagSFUtil*>::iterator it_lf = MapBTagSF.find(btag_key_lf);
  std::map<TString,BTagSFUtil*>::iterator it_hf = MapBTagSF.find(btag_key_hf);

  if(it_lf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}
  if(it_hf == MapBTagSF.end()){   Message("Combination of btagger and working point is not allowed. Check configation of MapBTagSF", ERROR);  exit(EXIT_FAILURE);}

  /// systematics allowed are +-1 and +-3 for HN analysis
  if ( tag == snu::KJet::JETPROB) return -999;
  
  bool isBtag=false;
  if (isData) {
    
    if (it_lf->second->IsTagged(jet.BJetTaggerValue(tag),  -999999, jet.Pt(), jet.Eta()))
      isBtag=true;
  }
    else if (jet.HadronFlavour() > 1){
      if (it_hf->second->IsTagged(jet.BJetTaggerValue(tag),  jet.HadronFlavour(),jet.Pt(), jet.Eta()))
        isBtag=true;
    }
    else{
      if (it_lf->second->IsTagged(jet.BJetTaggerValue(tag),  jet.HadronFlavour(),jet.Pt(), jet.Eta()))
        isBtag=true;
    }
  
  return isBtag;
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


bool AnalyzerCore::IsTight(snu::KElectron el){
  
  return eventbase->GetElectronSel()->PassUserID("ELECTRON_HN_TIGHT",el);

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
      if(keepfake&&keep_chargeflip) prompt_electrons.push_back(electrons.at(i));
      else if(keep_chargeflip&&electrons.at(i).MCMatched()) prompt_electrons.push_back(electrons.at(i));
      else if(keepfake&&! electrons.at(i).MCIsCF()) prompt_electrons.push_back(electrons.at(i)); 
      else if(electrons.at(i).MCMatched() && !electrons.at(i).MCIsCF()) prompt_electrons.push_back(electrons.at(i));
    }// Data
    else prompt_electrons.push_back(electrons.at(i));
  }/// loop

  return prompt_electrons;


}

vector<snu::KMuon> AnalyzerCore::GetTruePrompt(vector<snu::KMuon> muons, bool keepfake){
  if(muons.size()==0)return muons;

  vector<snu::KMuon> prompt_muons;

  for(unsigned int i = 0; i < muons.size(); i++){
    if(!k_isdata){

      if(keepfake) prompt_muons.push_back(muons.at(i));
      else if(muons.at(i).MCMatched()) prompt_muons.push_back(muons.at(i));
    }// Data
    else prompt_muons.push_back(muons.at(i));
  }/// loop
  return prompt_muons;

}


void AnalyzerCore::CorrectMuonMomentum(vector<snu::KMuon>& k_muons){
  
  mcdata_correction->CorrectMuonMomentum(k_muons,eventbase->GetTruth());
  Message("END CorrectMuonMomentum",DEBUG);    
  return;
    
}


void AnalyzerCore::SetCorrectedMomentum(vector<snu::KMuon>& k_muons){
  
  for(std::vector<snu::KMuon>::iterator it = k_muons.begin(); it != k_muons.end(); it++){
    it->SetRochPt(mcdata_correction->GetCorrectedMuonMomentum(*it, eventbase->GetTruth()));
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


//std::vector<snu::KJet> AnalyzerCore::GetTopBasicJets(/*std::vector<snu::KMuon>& muonColl, std::vector<snu::KElectron>& electronColl*/){
//    //// Jet Selection for Topfitter Unclustered energy estimation
//    //// Jets with minimum cuts are selected. PT should be>10GeV : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools
//    //// Lepton Veto part in Dr.Yu's function is deleted. since we are estimating unclustered energy. We don't need to abandon jets because of leptons
//    //// Unclustered energy should be Detector level definition not User analysis level definition; i.e. unclE=-sumjet-sumlep-met but the jets, leps are not user defines jet, leps but detector level(meaning datasetlevel) reconstructed jets, leps(wonder I need to include tau) 
//    //// 
//
//    std::vector<snu::KJet> jetColl;
//
//    eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
//    eventbase->GetJetSel()->SetPt(10.);
//    eventbase->GetJetSel()->SetEta(3.5);//Basically we need till eta=5 but this is the limit of SKtree
//    eventbase->GetJetSel()->Selection(jetColl);
////    eventbase->GetJetSel()->JetSelectionLeptonVeto(jetColl, muonColl, electronColl);
//
//    return jetColl;
//}


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
      if((T.DeltaR(MuonColl.at(i))<0.1)&&(dPtRel(T,MuonColl.at(i))<0.2)) muCandIdxColl.push_back(i);
    }
    if(muCandIdxColl.size()==0) muFinalCandIdx=-1;
    else if(muCandIdxColl.size()==1) muFinalCandIdx=muCandIdxColl.at(0);
    else if(muCandIdxColl.size()>1){
       for( int i=0; i<muCandIdxColl.size(); i++){
          dRoCut=T.DeltaR(MuonColl.at(muCandIdxColl.at(i)))/0.05;
          dPtReloCut=dPtRel(T,MuonColl.at(muCandIdxColl.at(i)))/0.2;
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


bool AnalyzerCore::GenDecayInfo(std::vector<snu::KTruth>& truthColl, TString Decaymode){
//Usage: TString: 3lv ; return true if hc decays 3l+v, false if hc decays 2mu+2j (in genlevel)
//                trimu: return true if final state is 3mu+MET+4j (in genlevel)
//                emumu: return true if final state is e+2mu+MET+4j (in genlevel)
//Algorithm: Based on information that we set t>hc b, tx>b W- => w+ is always from hc, so if we check whether w+ decay leptonically or hadronically then we know whether hc decay leptonically or hadronically
//          +I checked that for this sample there is no additional W which is not originated from signal process but NLO effect. all the W+ are from hc. Checked that with code for all events.

   bool trigger=false;

   if(Decaymode=="3lv"){
      for(int i=truthColl.size()-1; i>=0; i--){
         if(truthColl.at(i).IndexMother()==-1) continue;
           int pid=truthColl.at(i).PdgId();
           int MotherIdx=truthColl.at(i).IndexMother();
           int Motherpid=truthColl.at(MotherIdx).PdgId();
//         cout<<i<<" PID: "<<truthColl.at(i).PdgId()<<" GenStatus: "<<truthColl.at(i).GenStatus()<<" MotherIdx: "<<MotherIdx<<" MotherPID: "<<Motherpid<<" Pt: "<<truthColl.at(i).Pt()<<endl;

           //For algorithm, read description on headline.
           if(Motherpid==24){
              if(fabs(pid)<10) trigger=false;
              else if(fabs(pid)<20) trigger=true;//else if->fabs(pid)>10 is auto-included
           }
      }
   }

   if(Decaymode=="trimu"){
      for(int i=truthColl.size()-1; i>=0; i--){
         if(truthColl.at(i).IndexMother()==-1) continue;
           int pid=truthColl.at(i).PdgId();
           int MotherIdx=truthColl.at(i).IndexMother();
           int Motherpid=truthColl.at(MotherIdx).PdgId();

           if((fabs(Motherpid)==24)&&(fabs(pid)==13)) trigger=true;
      }
   }

   if(Decaymode=="emumu"){
      for(int i=truthColl.size()-1; i>=0; i--){
         if(truthColl.at(i).IndexMother()==-1) continue;
           int pid=truthColl.at(i).PdgId();
           int MotherIdx=truthColl.at(i).IndexMother();
           int Motherpid=truthColl.at(MotherIdx).PdgId();

           if((fabs(Motherpid)==24)&&(fabs(pid)==11)) trigger=true;
      }
   }

  return trigger;
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
  bool InAcceptance=Option.Contains("InAcceptance");
  for(unsigned int i=0; i<truthColl.size(); i++){
    if(truthColl.at(i).IndexMother() < 0 ) continue;
    int pid=truthColl.at(i).PdgId(), mpid=truthColl.at(truthColl.at(i).IndexMother()).PdgId();
   
    if(fabs(pid)==11 || fabs(pid)==13){
      if(!InAcceptance){
        if     (fabs(mpid)==23 || fabs(mpid)==24) NPromptLepton_EW++;
        else if(fabs(mpid)==32 || fabs(mpid)==36) NPromptLepton_BSM++;//For H+>AW analysis
      }
      else{
        if(truthColl.at(i).Pt()<10) continue;
        if((fabs(truthColl.at(i).Eta())>2.4 && fabs(pid)==13) || (fabs(truthColl.at(i).Eta())>2.5 && fabs(pid)==11)){
          continue;
        }
        if     (fabs(mpid)==23 || fabs(mpid)==24) NPromptLepton_EW++;
        else if(fabs(mpid)==32 || fabs(mpid)==36) NPromptLepton_BSM++;//For H+>AW analysis
      }
    }
  }
  NPromptLepton_Tot=NPromptLepton_EW+NPromptLepton_BSM;

  if     (Option.Contains("EW"))  return NPromptLepton_EW;
  else if(Option.Contains("BSM")) return NPromptLepton_BSM;

  return NPromptLepton_Tot;
}
