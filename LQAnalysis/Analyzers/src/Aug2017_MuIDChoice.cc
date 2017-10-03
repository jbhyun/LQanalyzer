/***************************************************************************
 * @Project: Aug2017_MuIDChoice 
 * @Package: LQanalyzer, ROOT-based analysis framework for Korea SNU / Main package author: John Almond
 *
 * @Code Author: Jihwan Bhyun 
 *
 ***************************************************************************/

/// Local includes
 #include "Aug2017_MuIDChoice.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Aug2017_MuIDChoice);

 Aug2017_MuIDChoice::Aug2017_MuIDChoice() : AnalyzerCore(), out_muons(0) {

   SetLogName("Aug2017_MuIDChoice");
   Message("In Aug2017_MuIDChoice constructor", INFO);
   InitialiseAnalysis();
 }


 void Aug2017_MuIDChoice::InitialiseAnalysis() throw( LQError ) {
   
   /// Initialise histograms
   MakeHistograms();  
   // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
   // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
   Message("Feb2016, HwA analysis", INFO);
   return;
 }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Loop///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Aug2017_MuIDChoice::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   if(!isData) weight*=MCweight;
   FillHist("GenWeight", MCweight, 1, -2., 2., 4);

   /// Acts on data to remove bad reconstructed event 
   //if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;

   m_logger<<DEBUG<<"RunNumber/Event Number = "<<eventbase->GetEvent().RunNumber()<<" : "<<eventbase->GetEvent().EventNumber()<<LQLogger::endmsg;
   m_logger<<DEBUG<<"isData = "<<isData<<LQLogger::endmsg;


   //Numbet of Vertex NoPUreweight
   FillHist("Nvtx_nocut_NoRW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);

   //Pileup Reweight
   float pileup_reweight=1.;
   if (!k_isdata) { pileup_reweight=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central); }

   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight*pileup_reweight, 0., 50., 50);


   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


  
   bool PunziComp=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("PunziComp"))      PunziComp      = true;
   }

    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //ListTriggersAvailable();

   //Normalisation Lumi Setting
   bool Pass_Trigger=true;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   //if     ( PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")   ) Pass_Trigger=true;
   //else if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;
   if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", TargetLumi);
   FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 1., 100);
   weight*=trigger_ps_weight*trigger_period_weight;


   //Trigger Cut
   if(!Pass_Trigger) return;
   FillCutFlow("TriggerCut", weight*pileup_reweight);
   /**********************************************************************************/

   //METFilterCut : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
   if(!PassMETFilter()) return;  FillCutFlow("EventCut", weight*pileup_reweight);

   //Vertex Cut
   //(vtx.ndof>4&&maxAbsZ<=0)||abs(vtx.z)<= 24)&&((maxd0 <=0)||abs(vtx.position.rho)<=2)&&!(vtx.isFake))
   if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight*pileup_reweight);
   


///////Objects in Analysis/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  
   /**PreSelCut***********************************************************************************************/
   //Intended for Code speed boosting up.
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
   //std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, false);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
   if(PunziComp){ if(!(muonPreColl.size()>=3)) return; }
   /**********************************************************************************************************/

   //For Fake Study
   //Muon ID's to Test
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
   std::vector<snu::KMuon> muonHN2FakeLColl; eventbase->GetMuonSel()->Selection(muonHN2FakeLColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
   std::vector<snu::KMuon> muonHN2FakeTColl; eventbase->GetMuonSel()->Selection(muonHN2FakeTColl, true);

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);
   std::vector<snu::KMuon> muonPOGTIsop15Coll; eventbase->GetMuonSel()->Selection(muonPOGTIsop15Coll, true);
   std::vector<snu::KMuon> muonPOGLIsop6IPp5p1Chi30Coll;
   std::vector<snu::KMuon> muonPOGTIsop15IPp01p1Chi3Coll, muonPOGTIsop20IPp01p1Chi3Coll, muonPOGTIsop25IPp01p1Chi3Coll;
     for(int i=0; i<muonPreColl.size(); i++){
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGTIsop15IPp01p1Chi3","Roch")) muonPOGTIsop15IPp01p1Chi3Coll.push_back(muonPreColl.at(i));
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGTIsop20IPp01p1Chi3","Roch")) muonPOGTIsop20IPp01p1Chi3Coll.push_back(muonPreColl.at(i));
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGTIsop25IPp01p1Chi3","Roch")) muonPOGTIsop25IPp01p1Chi3Coll.push_back(muonPreColl.at(i));
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop6IPp5p1Chi30","Roch"))  muonPOGLIsop6IPp5p1Chi30Coll.push_back(muonPreColl.at(i));
     }


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
     bool LeptonVeto=false;
   std::vector<snu::KJet> jetNoVetoColl; eventbase->GetJetSel()->Selection(jetNoVetoColl, LeptonVeto, muonHN2FakeLColl, electronLooseColl);
   std::vector<snu::KJet> bjetNoVetoColl = SelBJets(jetNoVetoColl, "Medium");

   std::vector<snu::KJet> jetColl  = SkimJetColl(jetNoVetoColl,  electronLooseColl, muonPOGLIsop6IPp5p1Chi30Coll, "EleMuVeto");
   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");


   float  met    = eventbase->GetEvent().MET(snu::KEvent::pfmet);
   float  metphi = eventbase->GetEvent().METPhi();
   double met_x  = eventbase->GetEvent().PFMETx();
   double met_y  = eventbase->GetEvent().PFMETy();
   double Pzv,Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));

   int Nvtx=eventbase->GetEvent().nVertices();


   /*****************************************************
   **Scale Factors
   *****************************************************/
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
   float trigger_sf=1.;
   float fake_weight=1.; bool EventCand=false;
   float geneff_weight=1., gennorm_weight=1.;

   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
//   if     (PunziComp) { if(muonHN2FakeLColl.size()==3)     EventCand=true; }
//   else if(SigRegFake && EMuMu) { if(muonHN2FakeLColl.size()==2 && electronLooseColl.size()==1) EventCand=true; }
//   else if(DiMuon_analysis){ if(muonLooseColl.size()==2)     EventCand=true; }
//   else if(DiEle_analysis) { if(electronLooseColl.size()==2) EventCand=true; }

   if(true){
   //if(SigRegFake){
     if(!isData){
      //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
  
      //id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
      //reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronFakeLColl);
  
      //id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
      //iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
      //trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
  
      //btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

       geneff_weight   = GenFilterEfficiency(k_sample_name);
       gennorm_weight  = SignalNorm(k_sample_name, 200.);
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
       //if(k_running_nonprompt){
         //fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", electronLooseColl.size());
       //}
     }
   }

   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf*geneff_weight*gennorm_weight;
   /***************************************************************************************************/



   //=============================================================================//
   // Main Analysis Code     -----------------------------------------------------//
   //=============================================================================//


   if(PunziComp){

     ComparePunziDisc(muonHN2FakeTColl, muonPOGLIsop6IPp5p1Chi30Coll,
       electronTightColl, electronLooseColl, jetColl, bjetColl, weight, "_HNTriT", "");
     ComparePunziDisc(muonPOGTIsop15Coll, muonPOGLIsop6IPp5p1Chi30Coll,
       electronTightColl, electronLooseColl, jetColl, bjetColl, weight, "_StdPOGT", "");
     ComparePunziDisc(muonPOGTIsop15IPp01p1Chi3Coll, muonPOGLIsop6IPp5p1Chi30Coll,
       electronTightColl, electronLooseColl, jetColl, bjetColl, weight, "_TIsop15", "");
     ComparePunziDisc(muonPOGTIsop20IPp01p1Chi3Coll, muonPOGLIsop6IPp5p1Chi30Coll,
       electronTightColl, electronLooseColl, jetColl, bjetColl, weight, "_TIsop20", "");
     ComparePunziDisc(muonPOGTIsop25IPp01p1Chi3Coll, muonPOGLIsop6IPp5p1Chi30Coll,
       electronTightColl, electronLooseColl, jetColl, bjetColl, weight, "_TIsop25", "");

   }


return;
}// End of execute event loop
  

void Aug2017_MuIDChoice::ComparePunziDisc(std::vector<snu::KMuon> muonColl, std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KElectron> electronColl, std::vector<snu::KElectron> electronLooseColl, std::vector<snu::KJet> jetColl, std::vector<snu::KJet> bjetColl, float weight, TString Label, TString Option){
  //Input : jetVetoColl, bjetVetoColl
  //Purpose : Compare punzi discriminants for each of mu ID, and choose best ID for the analysis

  if( !(muonLooseColl.size()==3 && electronLooseColl.size()==0) ) return;
  if( !(muonColl.size()==3 && electronColl.size()==0) ) return;
  if( fabs(SumCharge(muonColl))!=1 ) return;

  int IdxOS  = TriMuChargeIndex(muonColl, "OS");
  int IdxSS1 = TriMuChargeIndex(muonColl, "SS1");
  int IdxSS2 = TriMuChargeIndex(muonColl, "SS2");

  float MmumuOSSS1=(muonColl.at(IdxOS)+muonColl.at(IdxSS1)).M();
  float MmumuOSSS2=(muonColl.at(IdxOS)+muonColl.at(IdxSS2)).M();

  if( !(MmumuOSSS1>12. && MmumuOSSS2>12.) ) return;
  
  bool DiMuTrSel=false, TriMuTrSel1=false, TriMuTrSel2=false, DiMuTrig=false, TriMuTrig=false;
  if( muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>10 ) DiMuTrSel=true;
  if( muonColl.at(0).Pt()>15 && muonColl.at(1).Pt()>12 && muonColl.at(2).Pt()>7  ) TriMuTrSel1=true;
  if( muonColl.at(0).Pt()>20 && muonColl.at(1).Pt()>12 && muonColl.at(2).Pt()>10 ) TriMuTrSel2=true;
  if( PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")
     || PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) DiMuTrig=true;
  if( PassTrigger("HLT_TripleMu_12_10_5_v") ) TriMuTrig=true;
 
  bool JetSel=false, BJetSel=false, MmuSR=false;
  if(  jetColl.size()>=3 ) JetSel  = true;
  if( bjetColl.size()>=1 ) BJetSel = true;
  if(    MmumuOSSS2<40   ) MmuSR   = true;

  if(DiMuTrig && DiMuTrSel){
    bool PassSel=true;
    //Step1: 3lOS above trig cut & passtrig
    FillHist("CutFlow_DiMuTr"+Label, 0., weight, 0., 10., 10);
    //Step2: 3lOS+geq1b
    if(PassSel && BJetSel){
      FillHist("CutFlow_DiMuTr"+Label, 1., weight, 0., 10., 10);} else PassSel=false;
    //Step3: 3lOS+geq3j+geq1b
    if(PassSel && BJetSel && JetSel){
      FillHist("CutFlow_DiMuTr"+Label, 2., weight, 0., 10., 10);} else PassSel=false;
    //Step4: Off Z
    if(PassSel && BJetSel && JetSel && MmuSR ){
      FillHist("CutFlow_DiMuTr"+Label, 3., weight, 0., 10., 10);} else PassSel=false;
  }
  if(TriMuTrig && TriMuTrSel1){
    bool PassSel=true;
    //Step1: 3lOS above trig cut & passtrig
    FillHist("CutFlow_TriMuTr1"+Label, 0., weight, 0., 10., 10);
    //Step2: 3lOS+geq1b
    if(PassSel && BJetSel){
      FillHist("CutFlow_TriMuTr1"+Label, 1., weight, 0., 10., 10);} else PassSel=false;
    //Step3: 3lOS+geq3j+geq1b
    if(PassSel && BJetSel && JetSel){
      FillHist("CutFlow_TriMuTr1"+Label, 2., weight, 0., 10., 10);} else PassSel=false;
    //Step4: Off Z
    if(PassSel && BJetSel && JetSel && MmuSR ){
      FillHist("CutFlow_TriMuTr1"+Label, 3., weight, 0., 10., 10);} else PassSel=false;
  }
  if(TriMuTrig && TriMuTrSel2){
    bool PassSel=true;
    //Step1: 3lOS above trig cut & passtrig
    FillHist("CutFlow_TriMuTr2"+Label, 0., weight, 0., 10., 10);
    //Step2: 3lOS+geq1b
    if(PassSel && BJetSel){
      FillHist("CutFlow_TriMuTr2"+Label, 1., weight, 0., 10., 10);} else PassSel=false;
    //Step3: 3lOS+geq3j+geq1b
    if(PassSel && BJetSel && JetSel){
      FillHist("CutFlow_TriMuTr2"+Label, 2., weight, 0., 10., 10);} else PassSel=false;
    //Step4: Off Z
    if(PassSel && BJetSel && JetSel && MmuSR ){
      FillHist("CutFlow_TriMuTr2"+Label, 3., weight, 0., 10., 10);} else PassSel=false;
  }


}



void Aug2017_MuIDChoice::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Aug2017_MuIDChoice::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
//  string analysisdir = getenv("FILEDIR");  
//  if(!k_isdata) reweightPU = new Reweight((analysisdir + "SNUCAT_Pileup.root").c_str());

  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
//  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
//  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

Aug2017_MuIDChoice::~Aug2017_MuIDChoice() {
  
  Message("In Aug2017_MuIDChoice Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}




float Aug2017_MuIDChoice::FakeRateMC(snu::KElectron Ele, TString Option){

  float FR=0.;

  float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  float fEta=fabs(Ele.Eta());
  if(Option=="TT_powheg"){
    if(fEta<0.8){
      if(PTCorr<35) FR=0.1459;
      else if(PTCorr<50) FR=0.1138;
      else if(PTCorr<70) FR=0.09918;
      else if(PTCorr<100) FR=0.0831;
      else FR=0.09295;
    }
    else if(fEta<1.479){
      if(PTCorr<35) FR=0.1676;
      else if(PTCorr<50) FR=0.1313;
      else if(PTCorr<70) FR=0.1235;
      else if(PTCorr<100) FR=0.1129;
      else FR=0.1289;
    }
    else{
      if(PTCorr<35) FR=0.2567;
      else if(PTCorr<50) FR=0.2371;
      else if(PTCorr<70) FR=0.2257;
      else if(PTCorr<100) FR=0.2275;
      else FR=0.2407;
    }
 }

  return FR;
}


int Aug2017_MuIDChoice::GetFakeLepSrcIdx(snu::KMuon Mu, std::vector<snu::KTruth> TruthColl){

  int MatchedIdx=-1;
  int LepType=GetLeptonType(Mu,TruthColl);
  if(LepType==-1){
    float mindR=999., maxdR=0.4; int IdxMindR=-1;
    for(int j=2; j<TruthColl.size(); j++){
      if( TruthColl.at(j).GenStatus()!=1 ) continue;
      if( fabs(TruthColl.at(j).PdgId())<50 ) continue;
      if( !(TruthColl.at(j).Pt()>5 && fabs(TruthColl.at(j).Eta())<2.5) ) continue;
      if( Mu.DeltaR(TruthColl.at(j))<mindR ){
        mindR=Mu.DeltaR(TruthColl.at(j));
        IdxMindR=j;
      }
    }
    if(mindR<maxdR){ MatchedIdx=IdxMindR; }
  }
  if(LepType==-2 || LepType==-3){
    int MatchedTruthIdx = GenMatchedIdx(Mu,TruthColl);
    int MotherIdx       = FirstNonSelfMotherIdx(MatchedTruthIdx,TruthColl);
    int GrMotherIdx     = FirstNonSelfMotherIdx(MotherIdx,TruthColl);
    
    if(LepType==-2) MatchedIdx=MotherIdx;
    if(LepType==-3) MatchedIdx=GrMotherIdx;    
  }

  return MatchedIdx;
}



int Aug2017_MuIDChoice::GetFakeLepJetSrcType(snu::KMuon Mu, std::vector<snu::KJet> JetColl){
  //Type: -1: Unmatched, 1:L, 2:C, 3:B
  int SrcType=-1;
  bool NearB=false, NearC=false, NearL=false;
  for(int i=0; i<JetColl.size(); i++){
    if(Mu.DeltaR(JetColl.at(i))<0.4){
      if(JetColl.at(i).HadronFlavour()==5){ NearB=true; break; }//1)
      else if(JetColl.at(i).HadronFlavour()==4){ NearC=true; }
      else if(JetColl.at(i).HadronFlavour()==0){ NearL=true; }
    }
  }

  if     (NearB) SrcType=3;
  else if(NearC) SrcType=2;
  else if(NearL) SrcType=1;

  return SrcType;
//1) Higher Priority to B. if there's multiple near jets, then b-jet has higher priority
}

float Aug2017_MuIDChoice::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  //float PTCorr=Ele.Pt()*(1+min(0,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}


int Aug2017_MuIDChoice::StepPassed(std::vector<snu::KMuon> MuColl, std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, TString Option){

   int PassedSteps=0;
   bool EMuMu=false, TriMu=false;
   if     (Option.Contains("EMuMu")) EMuMu=true;
   else if(Option.Contains("TriMu")) TriMu=true;

   if(EMuMu){
     //Step1
     if(EleColl.size()==1 && MuColl.size()==2){
        if(EleColl.at(0).Pt()>25 && MuColl.at(0).Pt()>15 && MuColl.at(1).Pt()>10){PassedSteps++;} else return PassedSteps;
     }else return PassedSteps;
     

     //Step2
     if(MuColl.at(0).Charge()!=MuColl.at(1).Charge()) {PassedSteps++;} else return PassedSteps;
     //Step3
     if((MuColl.at(0)+MuColl.at(1)).M()>12)           {PassedSteps++;} else return PassedSteps;
     //Step4
     if(BJetColl.size()>=1)                           {PassedSteps++;} else return PassedSteps;
     //Step5
     if(JetColl.size()>=3)                            {PassedSteps++;} else return PassedSteps;
     //Step6
     if((MuColl.at(0)+MuColl.at(1)).M()>12 && (MuColl.at(0)+MuColl.at(1)).M()<40) {PassedSteps++;} else return PassedSteps;
   }
   else if(TriMu){
     //Step1
     if(EleColl.size()==0 && MuColl.size()==3){
       if(MuColl.at(0).Pt()>20 && MuColl.at(1).Pt()>10 && MuColl.at(2).Pt()>10 ){PassedSteps++;} else return PassedSteps;
     }else return PassedSteps;

     //Step2
     if(fabs(SumCharge(MuColl))==1){ PassedSteps++;} else return PassedSteps;

     int Idx_OS=TriMuChargeIndex(MuColl,"OS");
     int Idx_SS1=TriMuChargeIndex(MuColl,"SS1");
     int Idx_SS2=TriMuChargeIndex(MuColl,"SS2");
     float MOSSS1=(MuColl.at(Idx_OS)+MuColl.at(Idx_SS1)).M();
     float MOSSS2=(MuColl.at(Idx_OS)+MuColl.at(Idx_SS2)).M();
     //Step3
     if(MOSSS1>12 && MOSSS2>12){PassedSteps++;} else return PassedSteps;
     //Step4
     if(BJetColl.size()>=1){PassedSteps++;} else return PassedSteps;
     //Step5
     if(JetColl.size()>=3){PassedSteps++;} else return PassedSteps;
    
   }

   return PassedSteps;
}


int Aug2017_MuIDChoice::NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   if(Option.Contains("EWPrompt")) ReturnPrompt=true;
   if(Option.Contains("HFake"))    ReturnFake=true;

   for(int i=0; i<MuColl.size(); i++){
     int MuType=GetLeptonType(MuColl.at(i), truthColl);
     //if     (MuType>0 ) Nprompt++;

     if     (MuType>0 && MuType<4) Nprompt++;
     else if(MuType<0 && MuType>-5 ) Nfake++;
//     cout<<i<<" MuType "<<MuType<<endl;
   }

//   cout<<"Tot "<<Nprompt+Nfake<<" NPr "<<Nprompt<<" Nfake "<<Nfake<<endl;
   if(ReturnPrompt && ReturnFake) return Nprompt+Nfake;
   else if(ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}

int Aug2017_MuIDChoice::NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false, ReturnConv=false;
   if(Option.Contains("EWPrompt")) ReturnPrompt=true;
   if(Option.Contains("HFake"))    ReturnFake=true;
   if(Option.Contains("Conv"))     ReturnConv=true;

   for(int i=0; i<EleColl.size(); i++){
     int EleType=GetLeptonType(EleColl.at(i), truthColl);
     //if     (EleType>0 )  Nprompt++;

     if     (EleType>0 && EleType<4)  Nprompt++;
     else if(EleType<0 && EleType>-5) Nfake++;
     //else if(EleType<-4)              Nfake++;

//     cout<<i<<" ElType "<<EleType<<endl;
   }


//   cout<<"Tot "<<Nprompt+Nfake<<" NPr "<<Nprompt<<" Nfake "<<Nfake<<endl;
   if(ReturnPrompt && ReturnFake) return Nprompt+Nfake;
   else if(ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}



bool Aug2017_MuIDChoice::IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option){

  bool IsConversionCandidate=false;
  int EleType=GetLeptonType(Ele, TruthColl);
  int HardPhotonIdx=-1;
  if(EleType==1 || fabs(EleType)==2 || fabs(EleType)==3 ) return false; 

  for(int i=2; i<TruthColl.size(); i++){

    int pid=TruthColl.at(i).PdgId();
    int midx=TruthColl.at(i).IndexMother();
    int mpid=TruthColl.at(midx).PdgId();
    int GenSt=TruthColl.at(i).GenStatus();

    if( fabs(pid)==22 && ((GenSt>20 && GenSt<30) || GenSt==1 ) ) HardPhotonIdx=i;
  }
  if(HardPhotonIdx<2) return false;
  
  if(TruthColl.at(HardPhotonIdx).DeltaR(Ele)<0.4) IsConversionCandidate=true;

  return IsConversionCandidate;

}



void Aug2017_MuIDChoice::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow_W") && GetHist("cutflow_N")){
    GetHist("cutflow_W")->Fill(cut,weight);
    GetHist("cutflow_N")->Fill(cut,1);
  }
  else{
    if(!GetHist("cutflow_W")){
      AnalyzerCore::MakeHistograms("cutflow_W", 11, 0., 11.);
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(5,"3lCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(6,"M(#mu#mu)>4");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(7,"#geq1b");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(8,"2j");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(9,"3j");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
    if(!GetHist("cutflow_N")){
      AnalyzerCore::MakeHistograms("cutflow_N", 11, 0., 11.);
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(5,"3lCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(6,"M(#mu#mu)>4");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(7,"#geq1b");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(8,"2j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(9,"3j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
  }
}



void Aug2017_MuIDChoice::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Aug2017_MuIDChoice::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);



  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Aug2017_MuIDChoiceCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Aug2017_MuIDChoice::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
