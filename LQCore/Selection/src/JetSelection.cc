#include "JetSelection.h"
#include <iostream>


using namespace snu;

JetSelection::JetSelection(LQEvent ev) :BaseSelection() {
  k_lqevent = ev;
}

JetSelection::~JetSelection() {}

//// This code is used to make selection cuts to vectors of KJets


void JetSelection::BasicSelection(std::vector<KJet>& jetColl) {
  
  //// This is a basic set of cuts on jets

  std::vector<KJet> alljets = k_lqevent.GetJets();

  for (std::vector<KJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){
  
    if ( (jit->Pt() >= pt_cut_min) &&  (fabs(jit->Eta()) < eta_cut)){
      if ( PassUserID(PFJET_LOOSE, *jit) &&    (jit->Pt() >= pt_cut_min) &&  (fabs(jit->Eta()) < eta_cut))  jetColl.push_back(*jit);
    }
       
  }
}

  

void JetSelection::Selection(std::vector<KJet>& jetColl) {
  
  //// This is a basic set of cuts on jets
 
  std::vector<KJet> alljets = k_lqevent.GetJets();
  
  for (std::vector<KJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){
    
    bool pileupjet=false;
    if(applypileuptool) pileupjet =  ( !jit->PileupJetIDLoose());  ///---> CHECK THIS
    pileupjet=false;


    if(apply_ID) {
      if ( jit->Pt() >= pt_cut_min && jit->Pt() < pt_cut_max &&
	   fabs(jit->Eta()) < eta_cut
	   &&PassUserID(k_id, *jit) && !pileupjet)  jetColl.push_back(*jit);
    }
    else{
      if ( jit->Pt() >= pt_cut_min && jit->Pt() < pt_cut_max && 
	   fabs(jit->Eta()) < eta_cut
	   && PassUserID(PFJET_LOOSE, *jit)&& !pileupjet)  jetColl.push_back(*jit);
    }
  }
  BaseSelection::reset();
  return;
  
}
void JetSelection::JetHNSelection(std::vector<KJet>& jetColl, std::vector<KMuon> muonColl, std::vector<KElectron> electronColl){
  return JetHNSelection (jetColl, muonColl, electronColl, 20., 2.5, true, "Loose");
} 

void JetSelection::JetHNSelection(std::vector<KJet>& jetColl, std::vector<KMuon> muonColl, std::vector<KElectron> electronColl, float ptcut, float etacut, bool pileupID, TString ID ) {
  
  
  std::vector<KJet> pre_jetColl; 
  std::vector<KJet> alljets = k_lqevent.GetJets();
  
  for (std::vector<KJet>::iterator jit = alljets.begin(); jit!=alljets.end(); jit++){
    
    bool pass_pileupID = jit->PileupJetIDLoose();
    if(!pileupID) pass_pileupID = true; 
    
    if(ID.Contains("Loose")){
      if ( (jit->Pt() >= ptcut) && fabs(jit->Eta()) < etacut   && jit->PassLooseID()&&pass_pileupID)  pre_jetColl.push_back(*jit);
    }
    else if(ID.Contains("TightLepVeto")){
      if ( (jit->Pt() >= ptcut) && fabs(jit->Eta()) < etacut   && jit->PassTightLepVetoID()&&pass_pileupID)  pre_jetColl.push_back(*jit);
    } 
    else if(ID.Contains("Tight")){
      if ( (jit->Pt() >= ptcut) && fabs(jit->Eta()) < etacut   && jit->PassTightID()&&pass_pileupID)  pre_jetColl.push_back(*jit);
    } 
    else {cout << "Jet ID " << ID << " not found" << endl; exit(EXIT_FAILURE);}
  }
  

  //cout << "Number of loose jets = " << pre_jetColl.size() << endl;
  //cout << "Number of electrons = " << electronColl.size() << endl;
  //cout << "Number of muons = " << muonColl.size() << endl;
  
  for (UInt_t ijet = 0; ijet < pre_jetColl.size(); ijet++) {
    jetIsOK = true;
    for (UInt_t ilep = 0; ilep < muonColl.size(); ilep++) {
      if (muonColl[ilep].DeltaR( pre_jetColl[ijet] ) < 0.4) {
        jetIsOK = false;
	//cout << "Muon eta/phi = " << muonColl[ilep].Eta() << " " << muonColl[ilep].Phi() << endl;
	//cout << "Jet eta/phi = " <<  pre_jetColl[ijet].Eta() << " " <<  pre_jetColl[ijet].Phi() << endl;
	ilep = muonColl.size();
      }
    }/// End of muon loop
    for (UInt_t ilep = 0; ilep < electronColl.size(); ilep++) {
      if (electronColl[ilep].DeltaR( pre_jetColl[ijet] ) < 0.4 ) {
        jetIsOK = false;
	//cout <<"electron eta/phi =" << electronColl[ilep].Eta() << " " << electronColl[ilep].Phi() << endl;
	//cout <<"Jet eta/phi = " <<  pre_jetColl[ijet].Eta() <<" " <<pre_jetColl[ijet].Phi() << endl;
        ilep = electronColl.size();
      }
    }/// End of electron loop
    
    if (jetIsOK) jetColl.push_back( pre_jetColl[ijet] );
  }/// End of Jet loop
  
}

void JetSelection::JetSelectionLeptonVeto(std::vector<KJet>& jetColl, std::vector<KMuon> muonColl, std::vector<KElectron> electronColl) {
  
  //// This is a basic set of cuts on jets
  ///  + the jets are removed that are close to leptons
  
  std::vector<KJet> pre_jetColl;
  Selection(pre_jetColl);
  
  for (UInt_t ijet = 0; ijet < pre_jetColl.size(); ijet++) {
    jetIsOK = true;
    for (UInt_t ilep = 0; ilep < muonColl.size(); ilep++) {
      if (muonColl[ilep].DeltaR( pre_jetColl[ijet] ) < 0.4) {
	jetIsOK = false;
	ilep = muonColl.size();
      }
    }/// End of muon loop
    
    
    for (UInt_t ilep = 0; ilep < electronColl.size(); ilep++) {
      if (electronColl[ilep].DeltaR( pre_jetColl[ijet] ) < 0.4 ) {
	jetIsOK = false;
	ilep = electronColl.size();
      }
    }/// End of electron loop
    
    if (jetIsOK) jetColl.push_back( pre_jetColl[ijet] );
  } /// End of Jet loop
  
}



bool JetSelection::PassUserID (ID id, snu::KJet jet){ 
  if      ( id == PFJET_LOOSE  ) return PassUserID_PFJetLoose  (jet);
  else if ( id == PFJET_TIGHT  ) return PassUserID_PFJetTight  (jet);
  else return false;
}


bool JetSelection::PassUserID_PFJetLoose ( snu::KJet jet){
  
  return jet.PassLooseID();
}


bool JetSelection::PassUserID_PFJetTight ( snu::KJet jet)
{
  return jet.PassTightID();
}


JetSelection& JetSelection::operator= (const JetSelection& ms) {
  if(this != &ms){    
    BaseSelection::operator = (ms);
    k_lqevent = ms.k_lqevent;  
  }
  return *this;
};

JetSelection::JetSelection(const JetSelection& ms):
  BaseSelection(ms)
{
  k_lqevent = ms.k_lqevent;  
};



