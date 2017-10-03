// $Id: Jun2017_MCFakeStudy.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQJun2017_MCFakeStudy Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Jun2017_MCFakeStudy.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Jun2017_MCFakeStudy);

 Jun2017_MCFakeStudy::Jun2017_MCFakeStudy() : AnalyzerCore(), out_muons(0) {

   SetLogName("Jun2017_MCFakeStudy");
   Message("In Jun2017_MCFakeStudy constructor", INFO);
   InitialiseAnalysis();
 }


 void Jun2017_MCFakeStudy::InitialiseAnalysis() throw( LQError ) {
   
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

void Jun2017_MCFakeStudy::ExecuteEvents()throw( LQError ){

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


   bool EleFakeIDOpt=false, EleFakeParam=false, ClosureTest=false, FineTune=false, EleOnlyClosure=false, CompositionCheck=false;
   bool SelBiasTest=false, TrigBiasTest=false;
   bool NewClosureEMuMu=false;
  
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EleFakeIDOpt"))     EleFakeIDOpt=true;
     else if(k_flags.at(i).Contains("EleFakeParam"))     EleFakeParam=true;
     else if(k_flags.at(i).Contains("ClosureTest"))      ClosureTest=true;
     else if(k_flags.at(i).Contains("FineTune"))         FineTune=true;
     else if(k_flags.at(i).Contains("EleOnlyClosure"))   EleOnlyClosure=true;
     else if(k_flags.at(i).Contains("CompositionCheck")) CompositionCheck=true;
     else if(k_flags.at(i).Contains("TrigBiasTest"))     TrigBiasTest=true;
     else if(k_flags.at(i).Contains("SelBiasTest"))      SelBiasTest=true;
     else if(k_flags.at(i).Contains("NewClosureEMuMu"))  NewClosureEMuMu=true;
   }

    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //ListTriggersAvailable();

   //Normalisation Lumi Setting
   float trigger_ps_weight=1.;
   if(!isData) trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
   weight*=trigger_ps_weight;
   FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 1., 100);


   //Trigger Cut
   //if(!Pass_Trigger) return;
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
   if(ClosureTest || NewClosureEMuMu){ if( !((muonPreColl.size()>=2 && electronPreColl.size()>=1) || muonPreColl.size()>=3) ) return; }
   /**********************************************************************************************************/

   //For Fake Study
   //Muon ID's to Test
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
   std::vector<snu::KMuon> muonVetoColl; eventbase->GetMuonSel()->Selection(muonVetoColl, true);
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


   //Electron ID's to Test
   //
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
   std::vector<snu::KElectron> electronVetoColl; eventbase->GetElectronSel()->Selection(electronVetoColl);

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_TIGHT);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
   std::vector<snu::KElectron> electronCBFakeLColl; eventbase->GetElectronSel()->Selection(electronCBFakeLColl);


     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronFakeLPreOptColl; eventbase->GetElectronSel()->Selection(electronFakeLPreOptColl);

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronFakeLColl; eventbase->GetElectronSel()->Selection(electronFakeLColl);

     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronMVAMColl; eventbase->GetElectronSel()->Selection(electronMVAMColl);

   std::vector<snu::KElectron> electronNull;

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
   //  //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
   //  bool LeptonVeto=true;
     bool LeptonVeto=false;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonHN2FakeLColl, electronFakeLColl);


   //Method to apply 2a method SF
   //std::vector<int> bIdxColl  = GetSFBJetIdx(jetColl,"Medium");
   //std::vector<int> ljIdxColl = GetSFLJetIdx(jetColl, bIdxColl, "Medium");
   //std::vector<snu::KJet> bjetColl2a; for(int i=0; i<bIdxColl.size(); i++) {bjetColl2a.push_back(jetColl.at(bIdxColl.at(i)));}
   //std::vector<snu::KJet> ljetColl2a; for(int i=0; i<ljIdxColl.size(); i++){ljetColl2a.push_back(jetColl.at(ljIdxColl.at(i)));}

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   //std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

   double met = eventbase->GetEvent().PFMETType1();
   double met_x = eventbase->GetEvent().PFMETType1x();
   double met_y = eventbase->GetEvent().PFMETType1y();
   //double met = eventbase->GetEvent().MET();
   double Pzv,Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
   //snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
   //                     v[1].SetPx(met_x); v[1].SetPy(met_y);

   //int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();//njets-nbjets;//number of light jets
   int Nvtx=eventbase->GetEvent().nVertices();


   /*****************************************************
   **Scale Factors
   *****************************************************/
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
   float trigger_sf=1.;
   float fake_weight=1.; bool EventCand=false;
   float geneff_weight=1., gennorm_weight=1.;

   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
//   if     (SigRegFake && TriMu) { if(muonHN2FakeLColl.size()==3)     EventCand=true; }
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
       //  fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", electronLooseColl.size());
       //}
     }
   }

   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf*geneff_weight*gennorm_weight;
   /***************************************************************************************************/



/************************************************************************************/
//////////////////////////////////////////////////////////////////////////////////////
/////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************/


   if(EleFakeIDOpt){
     std::vector<snu::KElectron> FakeColl     = SkimLepColl(electronFakeLPreOptColl, truthColl, "HFake");
//     std::vector<snu::KElectron> PromptColl   = SkimLepColl(electronFakeLPreOptColl, truthColl, "Prompt");
     std::vector<snu::KElectron> FakeLColl     = SkimLepColl(electronFakeLColl, truthColl, "HFake");
     std::vector<snu::KJet>      JetCleanColl = SkimJetColl(jetColl, truthColl, "NoPrNoTau");
     //int NPromptGenLepAll = NPromptLeptons(truthColl,"InclTau");
     //int NConversion      = NPromptLeptons(truthColl,"OnlyConv");
     float HT=0; for(int i=0; i<JetCleanColl.size(); i++){ HT+=JetCleanColl.at(i).Pt(); }

     
     const int NIsoCuts=11, NMVACuts=21;
     float IsoCuts[NIsoCuts]={0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
     float MVACuts[NMVACuts]={-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};

     const int NPtEdges=7;
     float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 200.};


     //2D Scanning MVA-Iso
     for(int i=0; i<FakeColl.size(); i++){
       FillHist("TrigCutCheck_Iso_SumW", FakeColl.at(i).PFRelIso(0.3), weight, 0., 1., 100);
       FillHist("TrigCutCheck_MVA_SumW", FakeColl.at(i).MVA(), weight, 0., 1., 100);
       if(FakeColl.at(i).TriggerMatched("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v")){
         FillHist("TrigCutCheck_Iso_TrigSumW", FakeColl.at(i).PFRelIso(0.3), weight, 0., 1., 100);
         FillHist("TrigCutCheck_MVA_TrigSumW", FakeColl.at(i).MVA(), weight, 0., 1., 100);
       }

       FillHist("NFakeHasNearJet", 0., weight, 0., 5., 5);
       for(int j=0; j<JetCleanColl.size(); j++){
         if(FakeColl.at(i).DeltaR(JetCleanColl.at(j))<0.4){
           FillHist("NFakeHasNearJet", 1., weight, 0., 5., 5);

           for(int it_iso=1; it_iso<NIsoCuts; it_iso++){
             for(int it_mva=0; it_mva<NMVACuts-1; it_mva++){
               if(FakeColl.at(i).PFRelIso(0.3)<IsoCuts[it_iso] && FakeColl.at(i).MVA()>MVACuts[it_mva]){
                 FillHist("EleSumW_jMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                 if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleIDSumW_jMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);

                 //Eta Inclusive
                 if(JetCleanColl.at(j).HadronFlavour()==5){
                   FillHist("NFakeHasNearBJet", 2., weight, 0., 5., 5);
                   FillHist("EleSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                   if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleIDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                 }//End of B-nearjet loop
                 if(JetCleanColl.at(j).HadronFlavour()==4){
                   FillHist("NFakeHasNearCJet", 3., weight, 0., 5., 5);
                   FillHist("EleSumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                   if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleIDSumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                 }//End of C-near jet loop
                 if(JetCleanColl.at(j).HadronFlavour()==0){
                   FillHist("NFakeHasNearLJet", 4., weight, 0., 5., 5);
                   FillHist("EleSumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                   if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleIDSumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                 }//End of L-near jet loop
    
                 //Barrel1
                 if(fabs(FakeColl.at(i).Eta())<0.8){
                   if(JetCleanColl.at(j).HadronFlavour()==5){
                     FillHist("EleB1SumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                     if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleB1IDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                   }//End of B-nearjet loop
                   if(JetCleanColl.at(j).HadronFlavour()==4){
                     FillHist("EleB1SumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                     if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleB1IDSumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                   }//End of C-near jet loop
                   if(JetCleanColl.at(j).HadronFlavour()==0){
                     FillHist("EleB1SumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                     if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleB1IDSumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                   }//End of L-near jet loop
                 }
                 else if(fabs(FakeColl.at(i).Eta())<1.479){//Barrel2
                   if(JetCleanColl.at(j).HadronFlavour()==5){
                     FillHist("EleB2SumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                     if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleB2IDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                   }//End of B-nearjet loop
                   if(JetCleanColl.at(j).HadronFlavour()==4){
                     FillHist("EleB2SumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                     if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleB2IDSumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                   }//End of C-near jet loop
                   if(JetCleanColl.at(j).HadronFlavour()==0){
                     FillHist("EleB2SumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                     if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleB2IDSumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                   }//End of L-near jet loop
                 }
                 else if(fabs(FakeColl.at(i).Eta())<2.5){//EndCap
                   if(JetCleanColl.at(j).HadronFlavour()==5){
                     FillHist("EleESumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                     if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleEIDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                   }//End of B-nearjet loop
                   if(JetCleanColl.at(j).HadronFlavour()==4){
                     FillHist("EleESumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                     if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleEIDSumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                   }//End of C-near jet loop
                   if(JetCleanColl.at(j).HadronFlavour()==0){
                     FillHist("EleESumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
                     if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("EleEIDSumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
                   }//End of L-near jet loop

                 }
               }//End of if Iso<Cut && MVA<Cut
             }//End of MVACut loop
           }//End of IsoCut Loop

         }//End of Fake Matched Near Jet
       }//End of Jet Loop
     }//End of Fake Lepton Loop 

     //So Cut values : (Though Need Fine Tune)
     //(0.4  FR) - B1 : -0.2, 0.15 / B2 : -0.1  0.14 / E  : -0.3  0.25


     for(int i=0; i<FakeLColl.size(); i++){
       float PTCorr=ConeCorrectedPT(FakeLColl.at(i), 0.1);
       int FakeSrcType=GetFakeLepSrcType(FakeLColl.at(i), JetCleanColl);

       if(PTCorr<25) continue;

       FillHist("FakeSrcType", FakeSrcType, weight, -1., 4., 5);
       if( fabs(FakeLColl.at(i).Eta())<0.8 ){
         //AllFakes
         FillHist("EleB1SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
         if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleB1IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);

         //PerSources(Only for matched ones)
         if     (FakeSrcType==3){
           FillHist("EleB1SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
           if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleB1IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else if(FakeSrcType==2){
           FillHist("EleB1SumW_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
           if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleB1IDSumW_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else if(FakeSrcType==1){
           FillHist("EleB1SumW_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
           if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleB1IDSumW_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
       }
       else if( fabs(FakeLColl.at(i).Eta())<1.479 ){
         //AllFakes
         FillHist("EleB2SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
         if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleB2IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);

         //PerSources(Only for matched ones)
         if     (FakeSrcType==3){
           FillHist("EleB2SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
           if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleB2IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else if(FakeSrcType==2){
           FillHist("EleB2SumW_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
           if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleB2IDSumW_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else if(FakeSrcType==1){
           FillHist("EleB2SumW_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
           if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleB2IDSumW_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
       }
       else if( fabs(FakeLColl.at(i).Eta())<2.5 ){
         //AllFakes
         FillHist("EleESumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
         if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleEIDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);

         //PerSources(Only for matched ones)
         if     (FakeSrcType==3){
           FillHist("EleESumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
           if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleEIDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else if(FakeSrcType==2){
           FillHist("EleESumW_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
           if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleEIDSumW_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else if(FakeSrcType==1){
           FillHist("EleESumW_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
           if(PassIDCriteria(FakeLColl.at(i), "POGMVAMIP")) FillHist("EleEIDSumW_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }

       }

     }//End of Fake Lepton Loop 

     
   }//End of Ele Fake WP Scan
   if(FineTune){
     std::vector<snu::KElectron> FakeColl     = SkimLepColl(electronFakeLPreOptColl, truthColl, "HFake");
     std::vector<snu::KElectron> FakeB1Coll   = SkimLepColl(FakeColl, "B1");
     std::vector<snu::KElectron> FakeB2Coll   = SkimLepColl(FakeColl, "B2");
     std::vector<snu::KElectron> FakeEColl    = SkimLepColl(FakeColl, "E");
     //std::vector<snu::KElectron> FakeLColl    = SkimLepColl(electronFakeLColl, truthColl, "HFake");
     std::vector<snu::KJet>      JetCleanColl = SkimJetColl(jetColl, truthColl, "NoPrNoTau");
     float HT=0; for(int i=0; i<JetCleanColl.size(); i++){ HT+=JetCleanColl.at(i).Pt(); }

     //const int NIsoCuts=3, NMVACuts=3;
     //float IsoCuts[NIsoCuts]={0.3, 0.4, 0.5};
     //float MVACuts[NMVACuts]={-0.4, -0.2, 0.};


     const int NPtEdges=7;
     float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 200.};

     bool FullScanning=true, ShortScanning=false;
     if(FullScanning){
       //Full Range Scanning
       const int NIsoCuts=11, NMVACuts=21;
       float IsoCuts[NIsoCuts]={0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
       float MVACuts[NMVACuts]={-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
  
  
       for(int it_iso=0; it_iso<NIsoCuts; it_iso++){
         for(int it_mva=0; it_mva<NMVACuts; it_mva++){
           Draw1DFakePlot(FakeColl, JetCleanColl, MVACuts[it_mva], IsoCuts[it_iso], NPtEdges, PtEdges, "All");
           Draw1DFakePlot(FakeB1Coll, JetCleanColl, MVACuts[it_mva], IsoCuts[it_iso], NPtEdges, PtEdges, "B1");
           Draw1DFakePlot(FakeB2Coll, JetCleanColl, MVACuts[it_mva], IsoCuts[it_iso], NPtEdges, PtEdges, "B2");
           Draw1DFakePlot(FakeEColl, JetCleanColl, MVACuts[it_mva], IsoCuts[it_iso], NPtEdges, PtEdges, "E");
         }
       }
     }
     if(ShortScanning){
       Draw1DFakePlot(FakeColl, JetCleanColl, -0.5, -0.4, -0.5, 0.3, NPtEdges, PtEdges, "All");
       Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.5, -0.4, -0.5, 0.3, NPtEdges, PtEdges, "B1");
       Draw1DFakePlot(FakeB2Coll, JetCleanColl, -0.5, -0.4, -0.5, 0.3, NPtEdges, PtEdges, "B2");
       Draw1DFakePlot(FakeEColl, JetCleanColl, -0.5, -0.4, -0.5, 0.3, NPtEdges, PtEdges, "E");

       Draw1DFakePlot(FakeColl, JetCleanColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "All");
       Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "B1");
       Draw1DFakePlot(FakeB2Coll, JetCleanColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "B2");
       Draw1DFakePlot(FakeEColl, JetCleanColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "E");

       Draw1DFakePlot(FakeColl, JetCleanColl, -0.7, -0.6, -0.6, 0.5, NPtEdges, PtEdges, "All");
       Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.7, -0.6, -0.6, 0.5, NPtEdges, PtEdges, "B1");
       Draw1DFakePlot(FakeB2Coll, JetCleanColl, -0.7, -0.6, -0.6, 0.5, NPtEdges, PtEdges, "B2");
       Draw1DFakePlot(FakeEColl, JetCleanColl, -0.7, -0.6, -0.6, 0.5, NPtEdges, PtEdges, "E");
     }


   }
   if(CompositionCheck){
     
     std::vector<snu::KElectron> electronFakeL1Coll;
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMFakeLIso04Opt1")) electronFakeL1Coll.push_back(electronFakeLPreOptColl.at(i));
       }
     std::vector<snu::KElectron> electronFakeL2Coll;
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMFakeLIso04Opt2")) electronFakeL2Coll.push_back(electronFakeLPreOptColl.at(i));
       }


     for(int i=0; i<electronFakeL1Coll.size(); i++){
       if(k_sample_name.Contains("DY") || k_sample_name.Contains("TT")){
         if(PassIDCriteria(electronFakeL1Coll.at(i), "POGMVAMIP")) continue;
         int LepType=GetLeptonType(electronFakeL1Coll.at(i), truthColl);
         if(LepType>0 || LepType<-4) continue;
       }
       int FakeSrcType=GetFakeLepSrcType(electronFakeL1Coll.at(i), jetColl);
       float PTCorr=ConeCorrectedPT(electronFakeL1Coll.at(i), 0.1);
       FillHist("FakeSrcType_IDOpt1", FakeSrcType, weight, -5., 5., 10);
       FillHist("JetPtSrc_IDOpt1", PTCorr, weight, 0., 500., 25); 
     }
     for(int i=0; i<electronFakeL2Coll.size(); i++){
       if(k_sample_name.Contains("DY") || k_sample_name.Contains("TT")){
         if(PassIDCriteria(electronFakeL2Coll.at(i), "POGMVAMIP")) continue;
         int LepType=GetLeptonType(electronFakeL2Coll.at(i), truthColl);
         if(LepType>0 || LepType<-4) continue;
       }
       int FakeSrcType=GetFakeLepSrcType(electronFakeL2Coll.at(i), jetColl);
       float PTCorr=ConeCorrectedPT(electronFakeL2Coll.at(i), 0.1);
       FillHist("FakeSrcType_IDOpt2", FakeSrcType, weight, -5., 5., 10);
       FillHist("JetPtSrc_IDOpt2", PTCorr, weight, 0., 500., 25); 
     }
     for(int i=0; i<electronCBFakeLColl.size(); i++){
       if(k_sample_name.Contains("DY") || k_sample_name.Contains("TT")){
         if(PassIDCriteria(electronCBFakeLColl.at(i), "POGCBTIP")) continue;
         int LepType=GetLeptonType(electronCBFakeLColl.at(i), truthColl);
         if(LepType>0 || LepType<-4) continue;
       }
       int FakeSrcType=GetFakeLepSrcType(electronCBFakeLColl.at(i), jetColl);
       float PTCorr=ConeCorrectedPT(electronCBFakeLColl.at(i), 0.1);
       FillHist("FakeSrcType_IDCBL", FakeSrcType, weight, -5., 5., 10);
       FillHist("JetPtSrc_IDCBL", PTCorr, weight, 0., 500., 25); 
     }

   }
   if(EleOnlyClosure){

     std::vector<snu::KJet> JetCleanColl  = SkimJetColl(jetColl, truthColl, "NoPrNoTau");
     std::vector<snu::KJet> BJetCleanColl = SelBJets(JetCleanColl, "Medium");

     Draw1DClosurePlot(muonHN2FakeTColl, electronFakeLPreOptColl, JetCleanColl, BJetCleanColl, truthColl, met, -0.5, -0.4, -0.5, 0.3, "");
     Draw1DClosurePlot(muonHN2FakeTColl, electronFakeLPreOptColl, JetCleanColl, BJetCleanColl, truthColl, met, -0.7, -0.5, -0.6, 0.4, "");
     Draw1DClosurePlot(muonHN2FakeTColl, electronFakeLPreOptColl, JetCleanColl, BJetCleanColl, truthColl, met, -0.7, -0.6, -0.6, 0.5, "");

   }
   if(TrigBiasTest){
     std::vector<snu::KElectron> FakeColl     = SkimLepColl(electronFakeLPreOptColl, truthColl, "HFake");
     std::vector<snu::KElectron> FakeB1Coll   = SkimLepColl(FakeColl, "B1");
     std::vector<snu::KElectron> FakeB2Coll   = SkimLepColl(FakeColl, "B2");
     std::vector<snu::KElectron> FakeEColl    = SkimLepColl(FakeColl, "E");
     //std::vector<snu::KElectron> FakeLColl    = SkimLepColl(electronFakeLColl, truthColl, "HFake");
     std::vector<snu::KJet>      JetCleanColl = SkimJetColl(jetColl, truthColl, "NoPrNoTau");
     float HT=0; for(int i=0; i<JetCleanColl.size(); i++){ HT+=JetCleanColl.at(i).Pt(); }

     //const int NIsoCuts=3, NMVACuts=3;
     //float IsoCuts[NIsoCuts]={0.3, 0.4, 0.5};
     //float MVACuts[NMVACuts]={-0.4, -0.2, 0.};


     const int NPtEdges=8;
     float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 150., 200.};


     Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.2, -0.1, -0.3, 0.4, NPtEdges, PtEdges, "B1Opt1");
     Draw1DFakePlot(FakeB2Coll, JetCleanColl, -0.2, -0.1, -0.3, 0.4, NPtEdges, PtEdges, "B2Opt1");
     Draw1DFakePlot(FakeEColl, JetCleanColl, -0.2, -0.1, -0.3, 0.4, NPtEdges, PtEdges, "EOpt1");

     Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "B1Opt2");
     Draw1DFakePlot(FakeB2Coll, JetCleanColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "B2Opt2");
     Draw1DFakePlot(FakeEColl, JetCleanColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "EOpt2");

     if(PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")){
       Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.2, -0.1, -0.3, 0.4, NPtEdges, PtEdges, "B1Opt1Trig");
       Draw1DFakePlot(FakeB2Coll, JetCleanColl, -0.2, -0.1, -0.3, 0.4, NPtEdges, PtEdges, "B2Opt1Trig");
       Draw1DFakePlot(FakeEColl, JetCleanColl, -0.2, -0.1, -0.3, 0.4, NPtEdges, PtEdges, "EOpt1Trig");
  
       Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "B1Opt2Trig");
       Draw1DFakePlot(FakeB2Coll, JetCleanColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "B2Opt2Trig");
       Draw1DFakePlot(FakeEColl, JetCleanColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "EOpt2Trig");
     }

   }
   if(SelBiasTest){
     std::vector<snu::KElectron> electronFakeL1Coll;
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMFakeLIso04Opt1")) electronFakeL1Coll.push_back(electronFakeLPreOptColl.at(i));
       }
     std::vector<snu::KElectron> electronFakeL2Coll;
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMFakeLIso04Opt2")) electronFakeL2Coll.push_back(electronFakeLPreOptColl.at(i));
       }
     std::vector<snu::KJet> jetVeto1Coll= SkimJetColl(jetColl, electronFakeL1Coll, muonHN2FakeLColl, "EleMuVeto");
     std::vector<snu::KJet> jetVeto2Coll= SkimJetColl(jetColl, electronFakeL2Coll, muonHN2FakeLColl, "EleMuVeto");

     const int NPtEdges=8;
     float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 150., 200.};

     //ID1
     bool PassSel=true;
     if(electronFakeL1Coll.size()!=1          ) PassSel=false;//To remove DY Cont. && Ambiguity in meas.
     if(PassSel && jetVeto1Coll.size()==0     ) PassSel=false;//Should fire trigger
     if(PassSel && jetVeto1Coll.at(0).Pt()<40 ) PassSel=false;//Should fire trigger
     if(PassSel && jetVeto1Coll.at(0).DeltaR(electronFakeL1Coll.at(0))<1.0) PassSel=false;//the lepton is independent of this jet
     if(PassSel && met>30 ) PassSel=false;
     if(PassSel){
       float MTW=sqrt(2)*sqrt(met*electronFakeL1Coll.at(0).Pt()-met_x*electronFakeL1Coll.at(0).Px()-met_y*electronFakeL1Coll.at(0).Py());
       if(MTW>30) PassSel=false;
     } 
     if(PassSel){
       if(PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")){
         if(fabs(electronFakeL1Coll.at(0).Eta())<0.8){
           Draw1DFakePlot(electronFakeL1Coll, jetColl, -0.2, -0.1, -0.3, 0.4, NPtEdges, PtEdges, "B1Opt1Trig");
         }
         else if(fabs(electronFakeL1Coll.at(0).Eta())<1.479){
           Draw1DFakePlot(electronFakeL1Coll, jetColl, -0.2, -0.1, -0.3, 0.4, NPtEdges, PtEdges, "B2Opt1Trig");
         }
         else{
           Draw1DFakePlot(electronFakeL1Coll, jetColl, -0.2, -0.1, -0.3, 0.4, NPtEdges, PtEdges, "EOpt1Trig");
         }
       }
     }
     //ID2
     PassSel=true;
     if(electronFakeL2Coll.size()!=1          ) PassSel=false;//To remove DY Cont. && Ambiguity in meas.
     if(PassSel && jetVeto2Coll.size()==0     ) PassSel=false;//Should fire trigger
     if(PassSel && jetVeto2Coll.at(0).Pt()<40 ) PassSel=false;//Should fire trigger
     if(PassSel && jetVeto2Coll.at(0).DeltaR(electronFakeL2Coll.at(0))<1.0) PassSel=false;//the lepton is independent of this jet
     if(PassSel && met>30 ) PassSel=false;
     if(PassSel){
       float MTW=sqrt(2)*sqrt(met*electronFakeL2Coll.at(0).Pt()-met_x*electronFakeL2Coll.at(0).Px()-met_y*electronFakeL2Coll.at(0).Py());
       if(MTW>30) PassSel=false;
     } 
     if(PassSel){
       if(PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")){
         if(fabs(electronFakeL2Coll.at(0).Eta())<0.8){
           Draw1DFakePlot(electronFakeL2Coll, jetColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "B1Opt2Trig");
         }
         else if(fabs(electronFakeL2Coll.at(0).Eta())<1.479){
           Draw1DFakePlot(electronFakeL2Coll, jetColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "B2Opt2Trig");
         }
         else{
           Draw1DFakePlot(electronFakeL2Coll, jetColl, -0.7, -0.5, -0.6, 0.4, NPtEdges, PtEdges, "EOpt2Trig");
         }
       }
     }


   }
   if(ClosureTest){


     bool EMuMuClosure=true, TriMuClosure=true;
     bool IsCand=false;
     //m_datadriven_bkg->GetFakeObj()->SetUseQCDFake(true);
     float fakeweight=-1.;

     int NValidEleL = NPromptFake_Ele(electronFakeLColl, truthColl, "EWPromptHFake");
     int NValidMuL  = NPromptFake_Mu(muonHN2FakeLColl, truthColl, "EWPromptHFake");

     //Fake Composition Test
     for(int i=0; i<electronFakeLColl.size(); i++){
       int EleType=GetLeptonType(electronFakeLColl.at(i),truthColl);
       if(EleType>0 && EleType<4) continue;
       FillHist("EleFakeType", EleType, weight, -10., 10., 20);
       if(!PassIDCriteria(electronFakeLColl.at(i), "POGMVAMIP")){
         FillHist("EleFakeType_TfLp", EleType, weight, -10., 10., 20);
       }
     }
     for(int i=0; i<muonHN2FakeLColl.size(); i++){
       int MuType=GetLeptonType(muonHN2FakeLColl.at(i), truthColl);
       if(MuType>0 && MuType<4) continue;
       FillHist("MuFakeType", MuType, weight, -10., 10., 20);
       if(muonHN2FakeLColl.at(i).RelIso04()>0.1){
         FillHist("MuFakeType_TfLp", MuType, weight, -10., 10., 20);
       }
     }


     //Use only 3 valid leptons 
     if(EMuMuClosure){ if( NValidMuL==2 && NValidEleL==1 ) IsCand=true; }
     if(TriMuClosure){ if( NValidMuL==3 ) IsCand=true; }
       
     if(!IsCand) return;
     

     for(int i=0; i<muonHN2FakeLColl.size(); i++){
       if(muonHN2FakeLColl.at(i).RelIso04()>0.1){
 //       float FR=m_datadriven_bkg->GetFakeObj()->getTrilepFakeRate_muon(false, muonHN2FakeLColl.at(i).Pt(), muonHN2FakeLColl.at(i).Eta());
        float FR=0.;
        fakeweight*=-FR/(1-FR);
       }
     }
     for(int i=0; i<electronFakeLColl.size(); i++){
       if(!PassIDCriteria(electronFakeLColl.at(i), "POGMVAMIP")){
         //float FR=FakeRateMC(electronFakeLColl.at(i), "QCD_EGM");
         //float FR=FakeRateMC(electronFakeLColl.at(i), "DY");
         float FR=FakeRateMC(electronFakeLColl.at(i), "TT_powheg");
         fakeweight*=-FR/(1-FR);
       }
     }
     if(EMuMuClosure && electronFakeLColl.size()==1 && electronMVAMColl.size()==1 
                     && muonHN2FakeLColl.size()==2 && muonHN2FakeTColl.size()==2 ) fakeweight=0;
     if(TriMuClosure && muonHN2FakeLColl.size()==3 && muonHN2FakeTColl.size()==3 ) fakeweight=0;


     if(TriMuClosure){
       if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") || PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") ){
         if(muonHN2FakeLColl.size()==3 && electronFakeLColl.size()==0){
           int NStepPassed_exp=StepPassed(muonHN2FakeLColl, electronFakeLColl, jetColl, bjetColl, met, "TriMu");
           if(NStepPassed_exp>=5) FillHist("NStepPassed_3mu_exp", 5., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=4) FillHist("NStepPassed_3mu_exp", 4., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=3) FillHist("NStepPassed_3mu_exp", 3., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=2) FillHist("NStepPassed_3mu_exp", 2., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=1) FillHist("NStepPassed_3mu_exp", 1., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=0) FillHist("NStepPassed_3mu_exp", 0., weight*fakeweight, 0., 10., 10);
         }
         if(muonHN2FakeTColl.size()==3 && electronFakeLColl.size()==0 && NValidMuL==3 ){
           int NStepPassed_obs=StepPassed(muonHN2FakeTColl, electronMVAMColl, jetColl, bjetColl, met, "TriMu");
           if(NStepPassed_obs>=5) FillHist("NStepPassed_3mu_obs", 5., weight, 0., 10., 10);
           if(NStepPassed_obs>=4) FillHist("NStepPassed_3mu_obs", 4., weight, 0., 10., 10);
           if(NStepPassed_obs>=3) FillHist("NStepPassed_3mu_obs", 3., weight, 0., 10., 10);
           if(NStepPassed_obs>=2) FillHist("NStepPassed_3mu_obs", 2., weight, 0., 10., 10);
           if(NStepPassed_obs>=1) FillHist("NStepPassed_3mu_obs", 1., weight, 0., 10., 10);
           if(NStepPassed_obs>=0) FillHist("NStepPassed_3mu_obs", 0., weight, 0., 10., 10);
         }
       }
     }
     if(EMuMuClosure){
       int Pass_Trigger1=0, Pass_Trigger2=0;
       if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
       if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;
  
       float trigger_period_weight=1.;
       bool Pass_Trigger=false;
       if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_Trigger=true;
       trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;


       if(Pass_Trigger){
         weight *= trigger_period_weight;
         if(muonHN2FakeLColl.size()==2 && electronFakeLColl.size()==1 ){
           int NStepPassed_exp=StepPassed(muonHN2FakeLColl, electronFakeLColl, jetColl, bjetColl, met, "EMuMu");
           if(NStepPassed_exp>=6) FillHist("NStepPassed_1e2mu_exp", 6., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=5) FillHist("NStepPassed_1e2mu_exp", 5., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=4) FillHist("NStepPassed_1e2mu_exp", 4., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=3) FillHist("NStepPassed_1e2mu_exp", 3., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=2) FillHist("NStepPassed_1e2mu_exp", 2., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=1) FillHist("NStepPassed_1e2mu_exp", 1., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=0) FillHist("NStepPassed_1e2mu_exp", 0., weight*fakeweight, 0., 10., 10);
         }
         if(muonHN2FakeTColl.size()==2 && electronMVAMColl.size()==1 ){
           int NStepPassed_obs=StepPassed(muonHN2FakeTColl, electronMVAMColl, jetColl, bjetColl, met, "EMuMu");
           if(NStepPassed_obs>=6) FillHist("NStepPassed_1e2mu_obs", 6., weight, 0., 10., 10);
           if(NStepPassed_obs>=5) FillHist("NStepPassed_1e2mu_obs", 5., weight, 0., 10., 10);
           if(NStepPassed_obs>=4) FillHist("NStepPassed_1e2mu_obs", 4., weight, 0., 10., 10);
           if(NStepPassed_obs>=3) FillHist("NStepPassed_1e2mu_obs", 3., weight, 0., 10., 10);
           if(NStepPassed_obs>=2) FillHist("NStepPassed_1e2mu_obs", 2., weight, 0., 10., 10);
           if(NStepPassed_obs>=1) FillHist("NStepPassed_1e2mu_obs", 1., weight, 0., 10., 10);
           if(NStepPassed_obs>=0) FillHist("NStepPassed_1e2mu_obs", 0., weight, 0., 10., 10);
         }
       }
     }

   }
   if(NewClosureEMuMu){
     bool EMuMuClosure=true, TriMuClosure=true;
     bool IsCand=false;
     //m_datadriven_bkg->GetFakeObj()->SetUseQCDFake(true);
     float fakeweight=-1.;

     int NValidEleL = NPromptFake_Ele(electronFakeLColl, truthColl, "EWPromptHFake");
     int NValidMuL  = NPromptFake_Mu (muonHN2FakeLColl , truthColl, "EWPromptHFake");
     int NValidEleT = NPromptFake_Ele(electronMVAMColl , truthColl, "EWPromptHFake");
     int NValidMuT  = NPromptFake_Mu (muonHN2FakeTColl , truthColl, "EWPromptHFake");
     int NHFakeEleL = NPromptFake_Ele(electronFakeLColl, truthColl, "HFake");
     int NHFakeEleT = NPromptFake_Ele(electronMVAMColl , truthColl, "HFake");
     int NHFakeMuL  = NPromptFake_Mu (muonHN2FakeLColl , truthColl, "HFake");
     int NHFakeMuT  = NPromptFake_Mu (muonHN2FakeTColl , truthColl, "HFake");
     int NHFakeL = NHFakeEleL+NHFakeMuL;
     int NHFakeT = NHFakeEleT+NHFakeMuT;

     //Fake Composition Test
     for(int i=0; i<electronFakeLColl.size(); i++){
       int EleType=GetLeptonType(electronFakeLColl.at(i),truthColl);
       if(EleType>0 && EleType<4) continue;
       FillHist("EleFakeType", EleType, weight, -10., 10., 20);
       if(!PassIDCriteria(electronFakeLColl.at(i), "POGMVAMIP")){
         FillHist("EleFakeType_TfLp", EleType, weight, -10., 10., 20);
       }
     }
     for(int i=0; i<muonHN2FakeLColl.size(); i++){
       int MuType=GetLeptonType(muonHN2FakeLColl.at(i), truthColl);
       if(MuType>0 && MuType<4) continue;
       FillHist("MuFakeType", MuType, weight, -10., 10., 20);
       if(muonHN2FakeLColl.at(i).RelIso04()>0.1){
         FillHist("MuFakeType_TfLp", MuType, weight, -10., 10., 20);
       }
     }



     for(int i=0; i<muonHN2FakeLColl.size(); i++){
       if(muonHN2FakeLColl.at(i).RelIso04()>0.1){
        //float FR=m_datadriven_bkg->GetFakeObj()->getTrilepFakeRate_muon(false, muonHN2FakeLColl.at(i).Pt(), muonHN2FakeLColl.at(i).Eta());
        float FR=0.;
        fakeweight*=-FR/(1-FR);
       }
     }
     for(int i=0; i<electronFakeLColl.size(); i++){
       if(!PassIDCriteria(electronFakeLColl.at(i), "POGMVAMIP")){
         float FR=FakeRateMC(electronFakeLColl.at(i), "QCD_EGM");
         //float FR=FakeRateMC(electronFakeLColl.at(i), "DY");
         //float FR=FakeRateMC(electronFakeLColl.at(i), "TT_powheg");
         fakeweight*=-FR/(1-FR);
       }
     }
     if(EMuMuClosure && electronFakeLColl.size()==1 && electronMVAMColl.size()==1 
                     && muonHN2FakeLColl.size()==2 && muonHN2FakeTColl.size()==2 ) fakeweight=0;
     if(TriMuClosure && muonHN2FakeLColl.size()==3 && muonHN2FakeTColl.size()==3 ) fakeweight=0;


     if(TriMuClosure){
       if( PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") || PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") ){
         if(muonHN2FakeLColl.size()==3 && electronFakeLColl.size()==0 ){
         //if(muonHN2FakeLColl.size()==3 && electronFakeLColl.size()==0 && NValidMuL==3 ){
           int NStepPassed_exp=StepPassed(muonHN2FakeLColl, electronFakeLColl, jetColl, bjetColl, met, "TriMu");
           if(NStepPassed_exp>=5) FillHist("NStepPassed_3mu_exp", 5., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=4) FillHist("NStepPassed_3mu_exp", 4., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=3) FillHist("NStepPassed_3mu_exp", 3., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=2) FillHist("NStepPassed_3mu_exp", 2., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=1) FillHist("NStepPassed_3mu_exp", 1., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=0) FillHist("NStepPassed_3mu_exp", 0., weight*fakeweight, 0., 10., 10);
         }
         if(muonHN2FakeTColl.size()==3 && electronFakeLColl.size()==0 && NHFakeT>0){
         //if(muonHN2FakeTColl.size()==3 && electronFakeLColl.size()==0 && NValidMuT==3 ){
           int NStepPassed_obs=StepPassed(muonHN2FakeTColl, electronMVAMColl, jetColl, bjetColl, met, "TriMu");
           if(NStepPassed_obs>=5) FillHist("NStepPassed_3mu_obs", 5., weight, 0., 10., 10);
           if(NStepPassed_obs>=4) FillHist("NStepPassed_3mu_obs", 4., weight, 0., 10., 10);
           if(NStepPassed_obs>=3) FillHist("NStepPassed_3mu_obs", 3., weight, 0., 10., 10);
           if(NStepPassed_obs>=2) FillHist("NStepPassed_3mu_obs", 2., weight, 0., 10., 10);
           if(NStepPassed_obs>=1) FillHist("NStepPassed_3mu_obs", 1., weight, 0., 10., 10);
           if(NStepPassed_obs>=0) FillHist("NStepPassed_3mu_obs", 0., weight, 0., 10., 10);
         }
       }
     }
     if(EMuMuClosure){
       int Pass_Trigger1=0, Pass_Trigger2=0;
       if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
       if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;
  
       float trigger_period_weight=1.;
       bool Pass_Trigger=false;
       if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_Trigger=true;
       trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;


       if(Pass_Trigger){
         weight *= trigger_period_weight;
         //if(muonHN2FakeLColl.size()==2 && electronFakeLColl.size()==1 && NValidEleL==1 && NValidMuL==2){
         if(muonHN2FakeLColl.size()==2 && electronFakeLColl.size()==1 ){
           int NStepPassed_exp=StepPassed(muonHN2FakeLColl, electronFakeLColl, jetColl, bjetColl, met, "EMuMu");
           if(NStepPassed_exp>=6) FillHist("NStepPassed_1e2mu_exp", 6., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=5) FillHist("NStepPassed_1e2mu_exp", 5., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=4) FillHist("NStepPassed_1e2mu_exp", 4., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=3) FillHist("NStepPassed_1e2mu_exp", 3., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=2) FillHist("NStepPassed_1e2mu_exp", 2., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=1) FillHist("NStepPassed_1e2mu_exp", 1., weight*fakeweight, 0., 10., 10);
           if(NStepPassed_exp>=0) FillHist("NStepPassed_1e2mu_exp", 0., weight*fakeweight, 0., 10., 10);
         }
         //if(muonHN2FakeTColl.size()==2 && electronMVAMColl.size()==1 && NValidEleT==1 && NValidMuT==2 ){
         if(muonHN2FakeTColl.size()==2 && electronMVAMColl.size()==1 && NHFakeT>0){
           int NStepPassed_obs=StepPassed(muonHN2FakeTColl, electronMVAMColl, jetColl, bjetColl, met, "EMuMu");
           if(NStepPassed_obs>=6) FillHist("NStepPassed_1e2mu_obs", 6., weight, 0., 10., 10);
           if(NStepPassed_obs>=5) FillHist("NStepPassed_1e2mu_obs", 5., weight, 0., 10., 10);
           if(NStepPassed_obs>=4) FillHist("NStepPassed_1e2mu_obs", 4., weight, 0., 10., 10);
           if(NStepPassed_obs>=3) FillHist("NStepPassed_1e2mu_obs", 3., weight, 0., 10., 10);
           if(NStepPassed_obs>=2) FillHist("NStepPassed_1e2mu_obs", 2., weight, 0., 10., 10);
           if(NStepPassed_obs>=1) FillHist("NStepPassed_1e2mu_obs", 1., weight, 0., 10., 10);
           if(NStepPassed_obs>=0) FillHist("NStepPassed_1e2mu_obs", 0., weight, 0., 10., 10);
         }
       }
     }

   }
 
/////////////////////////////////////////////////////////////////////////////////// 

return;
}// End of execute event loop
  


void Jun2017_MCFakeStudy::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Jun2017_MCFakeStudy::BeginCycle() throw( LQError ){
  
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

Jun2017_MCFakeStudy::~Jun2017_MCFakeStudy() {
  
  Message("In Jun2017_MCFakeStudy Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


void Jun2017_MCFakeStudy::Draw1DClosurePlot(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> ElePreColl, std::vector<snu::KJet> jetColl, std::vector<snu::KJet> bjetColl, std::vector<snu::KTruth> truthColl, float met, float MVAB1Cut, float MVAB2Cut, float MVAECut, float IsoCut, TString Option){

  std::vector<snu::KElectron> EleLColl;
  std::vector<snu::KElectron> EleTColl;
    for(int i=0; i<ElePreColl.size(); i++){
      if( ElePreColl.at(i).PFRelIso(0.3)>IsoCut ) continue;

      if     ( fabs(ElePreColl.at(i).Eta())<0.8   ){ if( ElePreColl.at(i).MVA()<MVAB1Cut ) continue; }
      else if( fabs(ElePreColl.at(i).Eta())<1.479 ){ if( ElePreColl.at(i).MVA()<MVAB2Cut ) continue; }
      else                                         { if( ElePreColl.at(i).MVA()<MVAECut  ) continue; }
     
      EleLColl.push_back(ElePreColl.at(i));
      if(PassIDCriteria(ElePreColl.at(i), "POGMVAMIP")) EleTColl.push_back(ElePreColl.at(i)); 
    }

 
  //Checking 2PromptMu+1Fk Ele
  float fakeweight = -1.;
  int   NValidEleL = NPromptFake_Ele(EleLColl, truthColl, "EWPromptHFake");
  int   NValidMuL  = NPromptFake_Mu(MuPreColl, truthColl, "EWPrompt");
 
  //Use only 3 valid leptons 
  bool IsCand = false;
  if( NValidMuL==2 && NValidEleL==1 ) IsCand=true;
  if(!IsCand) return;

  std::ostringstream s1; s1<<IsoCut;
  TString Str_IsoCut=s1.str(); Str_IsoCut.ReplaceAll(".","p");


  for(int i=0; i<EleLColl.size(); i++){
    if(!PassIDCriteria(EleLColl.at(i), "POGMVAMIP")){
      //float FR=1.;
      //if(Str_IsoCut=="0p3") FR=FakeRateMC(EleLColl.at(i), "QCD_Iso03");
      //else if(Str_IsoCut=="0p4") FR=FakeRateMC(EleLColl.at(i), "QCD_Iso04");
      //else if(Str_IsoCut=="0p5") FR=FakeRateMC(EleLColl.at(i), "QCD_Iso05");
      float FR=FakeRateMC(EleLColl.at(i), "QCD_EGM");
      //float FR=FakeRateMC(EleLColl.at(i), "DY_Iso04");
      //float FR=FakeRateMC(EleLColl.at(i), "TT_Iso04");
      //float FR=FakeRateMC(EleLColl.at(i), "TT_powheg");
      //cout<<"FR "<<FR<<endl;
      fakeweight*=-FR/(1-FR);
    }
  }

  if( EleLColl.size()==1 && EleTColl.size()==1 ) fakeweight=0;

  int Pass_Trigger1=0, Pass_Trigger2=0;
  if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
  if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;
  
  float trigger_period_weight=1.;
  bool Pass_Trigger=false;
  if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_Trigger=true;
  trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;

  int EleType=GetLeptonType(EleLColl.at(0),truthColl);


  if(Pass_Trigger){
    weight *= trigger_period_weight;
    if(NValidMuL==2 && EleLColl.size()==1 ){
      int NStepPassed_exp=StepPassed(MuPreColl, EleLColl, jetColl, bjetColl, met, "EMuMu");
      if(NStepPassed_exp>=6) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_exp", 6., weight*fakeweight, 0., 10., 10);
      if(NStepPassed_exp>=5) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_exp", 5., weight*fakeweight, 0., 10., 10);
      if(NStepPassed_exp>=4) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_exp", 4., weight*fakeweight, 0., 10., 10);
      if(NStepPassed_exp>=3) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_exp", 3., weight*fakeweight, 0., 10., 10);
      if(NStepPassed_exp>=2) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_exp", 2., weight*fakeweight, 0., 10., 10);
      if(NStepPassed_exp>=1) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_exp", 1., weight*fakeweight, 0., 10., 10);
      if(NStepPassed_exp>=0) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_exp", 0., weight*fakeweight, 0., 10., 10);
    }
    if(NValidMuL==2 && EleTColl.size()==1 && EleType>=-4 && EleType<=-1 ){
      int NStepPassed_obs=StepPassed(MuPreColl, EleTColl, jetColl, bjetColl, met, "EMuMu");
      if(NStepPassed_obs>=6) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_obs", 6., weight, 0., 10., 10);
      if(NStepPassed_obs>=5) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_obs", 5., weight, 0., 10., 10);
      if(NStepPassed_obs>=4) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_obs", 4., weight, 0., 10., 10);
      if(NStepPassed_obs>=3) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_obs", 3., weight, 0., 10., 10);
      if(NStepPassed_obs>=2) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_obs", 2., weight, 0., 10., 10);
      if(NStepPassed_obs>=1) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_obs", 1., weight, 0., 10., 10);
      if(NStepPassed_obs>=0) FillHist("NStepPassed_1e2mu_iso"+Str_IsoCut+"_obs", 0., weight, 0., 10., 10);
    }
  }
  

}


void Jun2017_MCFakeStudy::Draw1DFakePlot(std::vector<snu::KElectron> FakeColl, std::vector<snu::KJet> JetColl, float MVACut, float IsoCut, int NPtEdges, float PtEdges[], TString Option){

  std::ostringstream s1,s2;
  s1<<MVACut; s2<<IsoCut;
  TString Str_MVACut=s1.str(),    Str_IsoCut=s2.str();
  Str_MVACut.ReplaceAll(".","p"); Str_IsoCut.ReplaceAll(".","p");
  Str_MVACut.ReplaceAll("-","m");
  
  TString EtaReg="";
  if(Option.Contains("B1")) EtaReg="B1";
  else if(Option.Contains("B2")) EtaReg="B2";
  else if(Option.Contains("E")) EtaReg="E";
  else if(Option.Contains("All")) EtaReg="";

  for(int i=0; i<FakeColl.size(); i++){
    if(FakeColl.at(i).MVA()<MVACut) continue;
    if(FakeColl.at(i).PFRelIso(0.3)>IsoCut) continue;

    float PTCorr=ConeCorrectedPT(FakeColl.at(i), 0.1);
    int FakeSrcType=GetFakeLepSrcType(FakeColl.at(i), JetColl);

    if(PTCorr<25) continue;
    //AllFakes
    FillHist("Ele"+EtaReg+"SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
    if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("Ele"+EtaReg+"IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);

    //PerSources(Only for matched ones)
    if     (FakeSrcType==3){
      FillHist("Ele"+EtaReg+"SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("Ele"+EtaReg+"IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else if(FakeSrcType==2){
      FillHist("Ele"+EtaReg+"SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("Ele"+EtaReg+"IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else if(FakeSrcType==1){
      FillHist("Ele"+EtaReg+"SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("Ele"+EtaReg+"IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }
  }//End of Fake Lepton Loop 


}


void Jun2017_MCFakeStudy::Draw1DFakePlot(std::vector<snu::KElectron> FakeColl, std::vector<snu::KJet> JetColl, float MVAB1Cut, float MVAB2Cut, float MVAECut, float IsoCut, int NPtEdges, float PtEdges[], TString Option){

  std::ostringstream s2; s2<<IsoCut;
  TString Str_IsoCut=s2.str();
  Str_IsoCut.ReplaceAll(".","p");
  
  TString EtaReg="";
  if(Option.Contains("B1")) EtaReg="B1";
  else if(Option.Contains("B2")) EtaReg="B2";
  else if(Option.Contains("E")) EtaReg="E";
  else if(Option.Contains("All")) EtaReg="";


  for(int i=0; i<FakeColl.size(); i++){
    if( EtaReg=="B1" && fabs(FakeColl.at(i).Eta())>=0.8 ) continue;
    if( EtaReg=="B2" && (fabs(FakeColl.at(i).Eta())<0.8 || fabs(FakeColl.at(i).Eta())>=1.479 ) ) continue;
    if( EtaReg=="E"  && fabs(FakeColl.at(i).Eta())<1.479 ) continue;

    if     ( fabs(FakeColl.at(i).Eta())<0.8   ){ if(FakeColl.at(i).MVA()<MVAB1Cut) continue; }
    else if( fabs(FakeColl.at(i).Eta())<1.479 ){ if(FakeColl.at(i).MVA()<MVAB2Cut) continue; }
    else                                       { if(FakeColl.at(i).MVA()<MVAECut)  continue; }

    if(FakeColl.at(i).PFRelIso(0.3)>IsoCut) continue;

    float PTCorr=ConeCorrectedPT(FakeColl.at(i), 0.1);
    int FakeSrcType=GetFakeLepSrcType(FakeColl.at(i), JetColl);

    if(PTCorr<25) continue;

    //AllFakes
    FillHist("Ele"+Option+"SumW_iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
    if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("Ele"+Option+"IDSumW_iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);

    //PerSources(Only for matched ones)
    if     (FakeSrcType==3){
      FillHist("Ele"+Option+"SumW_iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("Ele"+Option+"IDSumW_iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else if(FakeSrcType==2){
      FillHist("Ele"+Option+"SumW_iso"+Str_IsoCut+"_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("Ele"+Option+"IDSumW_iso"+Str_IsoCut+"_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else if(FakeSrcType==1){
      FillHist("Ele"+Option+"SumW_iso"+Str_IsoCut+"_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), "POGMVAMIP")) FillHist("Ele"+Option+"IDSumW_iso"+Str_IsoCut+"_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }

  }//End of Fake Lepton Loop 


}


float Jun2017_MCFakeStudy::FakeRateMC(snu::KElectron Ele, TString Option){

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
  else if(Option=="TT_Iso04"){
    if(fEta<0.8){
      if(PTCorr<35) FR=0.217223;
      else if(PTCorr<50) FR=0.097823;
      else if(PTCorr<70) FR=0.0765813;
      else if(PTCorr<100) FR=0.0651159;
      else FR=0.06603;
    }
    else if(fEta<1.479){
      if(PTCorr<35) FR=0.243832;
      else if(PTCorr<50) FR=0.115953;
      else if(PTCorr<70) FR=0.106106;
      else if(PTCorr<100) FR=0.0941181;
      else FR=0.0965011;
    }
    else{
      if(PTCorr<35) FR=0.314817;
      else if(PTCorr<50) FR=0.193874;
      else if(PTCorr<70) FR=0.177807;
      else if(PTCorr<100) FR=0.172306;
      else FR=0.171792;
    }
  }
  else if(Option=="QCD_EGM"){
    if(fEta<0.8){
      //if(PTCorr<35) FR=0.09618;//@PT10
      if(PTCorr<35) FR=0.1636;//@PT25
      else if(PTCorr<50) FR=0.1187;
      else if(PTCorr<70) FR=0.1699;
      else if(PTCorr<100) FR=0.1384;
      else FR=0.1430;
    }
    else if(fEta<1.479){
      //if(PTCorr<35) FR=0.1297;//@PT10
      if(PTCorr<35) FR=0.2244;//@PT25
      else if(PTCorr<50) FR=0.1582;
      else if(PTCorr<70) FR=0.158;
      else if(PTCorr<100) FR=0.1148;
      else FR=0.2209;
    }
    else{
      //if(PTCorr<35) FR=0.2479;//@PT10
      if(PTCorr<35) FR=0.3798;//@PT25
      else if(PTCorr<50) FR=0.2129;
      else if(PTCorr<70) FR=0.2157;
      else if(PTCorr<100) FR=0.2727;
      else FR=0.2477;
    }
  }
  else if(Option=="QCD_Iso03"){
    if(fEta<0.8){
      if(PTCorr<35) FR=0.1389;
      else if(PTCorr<50) FR=0.1136;
      else if(PTCorr<70) FR=0.1549;
      else if(PTCorr<100) FR=0.1158;
      else FR=0.1187;
    }
    else if(fEta<1.479){
      if(PTCorr<35) FR=0.1911;
      else if(PTCorr<50) FR=0.1478;
      else if(PTCorr<70) FR=0.1435;
      else if(PTCorr<100) FR=0.09653;
      else FR=0.1797;
    }
    else{
      if(PTCorr<35) FR=0.3239;//@PT25
      else if(PTCorr<50) FR=0.1998;
      else if(PTCorr<70) FR=0.1908;
      else if(PTCorr<100) FR=0.2346;
      else FR=0.2070;
    }
  }
  else if(Option=="QCD_Iso04"){
    if(fEta<0.8){
      if(PTCorr<35) FR=0.1185;//@PT25
      else if(PTCorr<50) FR=0.08259;
      else if(PTCorr<70) FR=0.1231;
      else if(PTCorr<100) FR=0.08276;
      else FR=0.08602;
    }
    else if(fEta<1.479){
      if(PTCorr<35) FR=0.1737;
      else if(PTCorr<50) FR=0.1199;
      else if(PTCorr<70) FR=0.1153;
      else if(PTCorr<100) FR=0.07766;
      else FR=0.1483;
    }
    else{
      if(PTCorr<35) FR=0.2912;//@PT25
      else if(PTCorr<50) FR=0.1568;
      else if(PTCorr<70) FR=0.15865;
      else if(PTCorr<100) FR=0.193951;
      else FR=0.166254;
    }
  }
  else if(Option=="QCD_Iso05"){
    if(fEta<0.8){
      if(PTCorr<35) FR=0.1185;//@PT25
      else if(PTCorr<50) FR=0.07862;
      else if(PTCorr<70) FR=0.1182;
      else if(PTCorr<100) FR=0.0796;
      else FR=0.08313;
    }
    else if(fEta<1.479){
      if(PTCorr<35) FR=0.1593;//@PT25
      else if(PTCorr<50) FR=0.1043;
      else if(PTCorr<70) FR=0.09615;
      else if(PTCorr<100) FR=0.07045;
      else FR=0.1297;
    }
    else{
      if(PTCorr<35) FR=0.291249;//@PT25
      else if(PTCorr<50) FR=0.149709;
      else if(PTCorr<70) FR=0.149482;
      else if(PTCorr<100) FR=0.185653;
      else FR=0.16068;
    }
  }
  else if(Option=="DY"){
    if(fEta<0.8){
      if(PTCorr<35) FR=0.222;//@PT25
      else if(PTCorr<50) FR=0.215;
      else if(PTCorr<70) FR=0.183;
      else if(PTCorr<100) FR=0.219;
      else FR=0.2143;
    }
    else if(fEta<1.479){
      if(PTCorr<35) FR=0.2266;//@PT25
      else if(PTCorr<50) FR=0.2278;
      else if(PTCorr<70) FR=0.2705;
      else if(PTCorr<100) FR=0.2280;
      else FR=0.4105;
    }
    else{
      if(PTCorr<35) FR=0.3865;//@PT25
      else if(PTCorr<50) FR=0.3448;
      else if(PTCorr<70) FR=0.3467;
      else if(PTCorr<100) FR=0.3912;
      else FR=0.3640;
    }
  }
  else if(Option=="DY_Iso04"){
    if(fEta<0.8){
      if(PTCorr<35) FR=0.148342;//@PT25
      else if(PTCorr<50) FR=0.144186;
      else if(PTCorr<70) FR=0.129305;
      else if(PTCorr<100) FR=0.138778;
      else FR=0.126477;
    }
    else if(fEta<1.479){
      if(PTCorr<35) FR=0.153721;//@PT25
      else if(PTCorr<50) FR=0.17511;
      else if(PTCorr<70) FR=0.180541;
      else if(PTCorr<100) FR=0.171512;
      else FR=0.294216;
    }
    else{
      if(PTCorr<35) FR=0.285454;//@PT25
      else if(PTCorr<50) FR=0.256927;
      else if(PTCorr<70) FR=0.262726;
      else if(PTCorr<100) FR=0.295326;
      else FR=0.212645;
    }
  }

  return FR;
}

int Jun2017_MCFakeStudy::GetFakeLepSrcType(snu::KElectron Ele, std::vector<snu::KJet> JetColl){
  //Type: -1: Unmatched, 1:L, 2:C, 3:B
  int SrcType=-1;
  bool NearB=false, NearC=false, NearL=false;
  for(int i=0; i<JetColl.size(); i++){
    if(Ele.DeltaR(JetColl.at(i))<0.4){
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

float Jun2017_MCFakeStudy::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  //float PTCorr=Ele.Pt()*(1+min(0,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}


int Jun2017_MCFakeStudy::StepPassed(std::vector<snu::KMuon> MuColl, std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, TString Option){

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


int Jun2017_MCFakeStudy::NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   if(Option.Contains("EWPrompt")) ReturnPrompt=true;
   if(Option.Contains("HFake"))    ReturnFake=true;

   for(int i=0; i<MuColl.size(); i++){
     int MuType=GetLeptonType(MuColl.at(i), truthColl);
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

int Jun2017_MCFakeStudy::NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   if(Option.Contains("EWPrompt")) ReturnPrompt=true;
   if(Option.Contains("HFake"))    ReturnFake=true;

   for(int i=0; i<EleColl.size(); i++){
     int EleType=GetLeptonType(EleColl.at(i), truthColl);
     if     (EleType>0 && EleType<4)  Nprompt++;
     else if(EleType<0 && EleType>-5) Nfake++;

//     cout<<i<<" ElType "<<EleType<<endl;
   }


//   cout<<"Tot "<<Nprompt+Nfake<<" NPr "<<Nprompt<<" Nfake "<<Nfake<<endl;
   if(ReturnPrompt && ReturnFake) return Nprompt+Nfake;
   else if(ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}



bool Jun2017_MCFakeStudy::IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option){

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



void Jun2017_MCFakeStudy::FillCutFlow(TString cut, float weight){
  
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



void Jun2017_MCFakeStudy::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Jun2017_MCFakeStudy::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);



  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Jun2017_MCFakeStudyCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Jun2017_MCFakeStudy::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
