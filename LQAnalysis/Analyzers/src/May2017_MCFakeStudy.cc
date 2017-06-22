// $Id: May2017_MCFakeStudy.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQMay2017_MCFakeStudy Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "May2017_MCFakeStudy.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (May2017_MCFakeStudy);

 May2017_MCFakeStudy::May2017_MCFakeStudy() : AnalyzerCore(), out_muons(0) {

   SetLogName("May2017_MCFakeStudy");
   Message("In May2017_MCFakeStudy constructor", INFO);
   InitialiseAnalysis();
 }


 void May2017_MCFakeStudy::InitialiseAnalysis() throw( LQError ) {
   
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

void May2017_MCFakeStudy::ExecuteEvents()throw( LQError ){

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


   bool EleFake=false, BasicCompositionStudy=false, SigRegFake=false, EMuMu=false, TriMu=false, MatchingTest=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EleFake"))    EleFake=true;
     else if(k_flags.at(i).Contains("BasicCompositionStudy")) BasicCompositionStudy=true;
     else if(k_flags.at(i).Contains("SigRegFake")) SigRegFake=true;
     else if(k_flags.at(i).Contains("EMuMu"))      EMuMu=true;
     else if(k_flags.at(i).Contains("TriMu"))      TriMu=true;
     else if(k_flags.at(i).Contains("MatchingTest")) MatchingTest=true;
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
   //  if     (MuEff ) { if( !(muonPreColl.size()>=1)     ) return; }
   //  else if(EleEff) { if( !(electronPreColl.size()>=1) ) return; }
     if     (SigRegFake && EMuMu){ if( !(muonPreColl.size()>=2 && electronPreColl.size()>=1) ) return; }
     else if(SigRegFake && TriMu){ if( !(muonPreColl.size()>=3) ) return; }
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
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(3.);
   std::vector<snu::KElectron> electronFakeLColl; eventbase->GetElectronSel()->Selection(electronFakeLColl);//For POG Type Study

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
     //bool LeptonVeto=true;
     bool LeptonVeto=false;
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonHN2FakeLColl, electronVetoColl);


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

   //////Basic Objects Check//////////////////////////////////////////////////////////
/*   FillHist("Basic_Nj_orig", njets, weight, 0., 10., 10);
   FillHist("Basic_Nb_orig", nbjets, weight, 0., 10., 10);//Bjet def CVSincV2>0.89;pfKJet::CSVv2IVF_discriminator 
   FillHist("Basic_NmuT_orig", muonColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NmuL_orig", muonLooseColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeT_orig", electronColl.size(), weight, 0., 10., 10);
   FillHist("Basic_NeL_orig", electronLooseColl.size(), weight, 0., 10., 10);
   if(electronColl.size()==1) FillHist("Basic_NmuT_weT", muonColl.size(), weight, 0., 10., 10);
   if(TriMu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_3mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_3mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>2) FillHist("Basic_Ptmu3_3mu_orig", muonColl.at(2).Pt(), weight, 0., 200., 200);
   }
   if(EMuMu_analysis){
     if(muonColl.size()>0) FillHist("Basic_Ptmu1_1e2mu_orig", muonColl.at(0).Pt(), weight, 0., 200., 200);
     if(muonColl.size()>1) FillHist("Basic_Ptmu2_1e2mu_orig", muonColl.at(1).Pt(), weight, 0., 200., 200);
     if(electronColl.size()>0) FillHist("Basic_Pte_1e2mu_orig", electronColl.at(0).Pt(), weight, 0., 200., 200);
   }
   if(njets>0)  FillHist("Basic_j1_Et_orig", jetColl.at(0).Et(), weight, 0., 200., 200);
   if(njets>1)  FillHist("Basic_j2_Et_orig", jetColl.at(1).Et(), weight, 0., 200., 200);
   if(njets>2)  FillHist("Basic_j3_Et_orig", jetColl.at(2).Et(), weight, 0., 200., 200);
   if(njets>3)  FillHist("Basic_j4_Et_orig", jetColl.at(3).Et(), weight, 0., 200., 200);
   if(nbjets>0) FillHist("Basic_b1_Et_orig", bjetColl.at(0).Et(), weight, 0, 200., 200);
   if(nbjets>1) FillHist("Basic_b2_Et_orig", bjetColl.at(1).Et(), weight, 0., 200., 200);
   FillHist("Basic_METdist_orig", met, weight, 0., 200., 100);*/
   ///////////////////////////////////////////////////////////////////////////////////


/************************************************************************************/
//////////////////////////////////////////////////////////////////////////////////////
/////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************/

   //ToCheck
   //1. Iso Dist of fake - DY vs. TT vs. QCD
   //  --Incl & Eta & HT
   //  --Also investigate source
   //  --near jet bjet PT,eta
   //  --RelIso, AbsIso (Is it possible to get only PV track CH component?)
   //  --In Case of tops check depending on decay channel(0,1,2lep); More the leps, less the HT
   //
   //2. HT per Nfake 
   //3. PT per source PT
   //4. dR from near jet - kind of jet?

   if(EleFake){
     std::vector<snu::KElectron> FakeColl     = SkimLepColl(electronFakeLColl, truthColl, "Fake");
     std::vector<snu::KElectron> PromptColl   = SkimLepColl(electronFakeLColl, truthColl, "Prompt");
     std::vector<snu::KJet>      JetCleanColl = SkimJetColl(jetColl, truthColl, "NoPrNoTau");
     int NPromptGenLepAll = NPromptLeptons(truthColl,"InclTau");
     int NConversion      = NPromptLeptons(truthColl,"OnlyConv");
     FillHist("NConv", NConversion, 1., -10., 10., 20);
     FillHist("NprL", NPromptGenLepAll, 1., -10., 10., 20);

     //PrintTruth(); 

     //Not limited to Acc & Tau daughter because hadronic activity hugely depends on the number of W decaying hadronically
     float HT=0; for(int i=0; i<JetCleanColl.size(); i++){ HT+=JetCleanColl.at(i).Pt(); }
     
     //Nnear Jets
     std::vector<int> NJetNearFakeLep;
     for(int i=0; i<FakeColl.size(); i++){
       int JetCount=0;
       for(int j=0; j<JetCleanColl.size(); j++){
         if(FakeColl.at(i).DeltaR(JetCleanColl.at(j))<0.4) JetCount++;
       }
       NJetNearFakeLep.push_back(JetCount);
       FillHist("NJetNearFakeLep", JetCount, weight, 0., 10., 10);
     }

     if(NPromptGenLepAll==0){
       FillHist("HT_0l", HT, weight, 0., 1000., 50);
       for(int i=0; i<PromptColl.size(); i++){
         FillHist("PromptElRelIso_0l", PromptColl.at(i).PFRelIso(0.3), weight, 0., 1., 20);

         FillHist("PromptElAvgRelIso_HT_1D_0l", HT, weight*PromptColl.at(i).PFRelIso(0.3), 0., 1000., 10);
         FillHist("PromptElSumW_HT_1D_0l", HT, weight, 0., 1000., 10);
         FillHist("PromptElAvgRelIso_fEta_1D_0l", fabs(PromptColl.at(i).Eta()), weight*PromptColl.at(i).PFRelIso(0.3), 0., 2.5, 5);
         FillHist("PromptElSumW_fEta_1D_0l", fabs(PromptColl.at(i).Eta()), weight, 0., 2.5, 5);

         FillHist("PromptPT_PromptfEta_1D_0l", fabs(PromptColl.at(i).Eta()), PromptColl.at(i).Pt()*weight, 0., 2.5, 5);
         FillHist("PromptPT_PromptfEta_2D_0l", PromptColl.at(i).Pt(), fabs(PromptColl.at(i).Eta()), weight, 0., 300., 60, 0., 2.5, 5);
       }
       for(int i=0; i<FakeColl.size(); i++){
         FillHist("NCount_HasCloseJet_0l", 0., weight, 0., 3., 3);
         if(NearEWLep(FakeColl.at(i), truthColl, "InclTau")) FillHist("NCount_HasCloseJet_0l", 1., weight, 0., 3., 3);

         FillHist("FakeElSumW_HT_1D_0l", HT, weight, 0., 1000., 10);
         FillHist("FakeElSumW_fEta_1D_0l", fabs(FakeColl.at(i).Eta()), weight, 0., 2.5, 5);
         FillHist("FakeElSumW_Pt_1D_0l", FakeColl.at(i).Pt(), weight, 0., 200., 40);
         FillHist("FakeElSumW_HT_fEta_2D_0l", HT, fabs(FakeColl.at(i).Eta()), weight, 0., 1000., 10, 0., 2.5, 5);

         FillHist("FakeElRelIso_0l", FakeColl.at(i).PFRelIso(0.3), weight, 0., 1., 20);
         FillHist("FakeElAvgRelIso_HT_1D_0l", HT, weight*FakeColl.at(i).PFRelIso(0.3), 0., 1000., 10);
         FillHist("FakeElAvgRelIso_fEta_1D_0l", fabs(FakeColl.at(i).Eta()), weight*FakeColl.at(i).PFRelIso(0.3), 0., 2.5, 5);

         FillHist("FakePT_FakefEta_1D_0l", fabs(FakeColl.at(i).Eta()), FakeColl.at(i).Pt()*weight, 0., 2.5, 5);
         FillHist("FakePT_FakefEta_2D_0l", FakeColl.at(i).Pt(), fabs(FakeColl.at(i).Eta()), weight, 0., 200, 40, 0., 2.5, 5);

         if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
           FillHist("FakeIDElSumW_HT_1D_0l", HT, weight, 0., 1000., 10);
           FillHist("FakeIDElSumW_fEta_1D_0l", fabs(FakeColl.at(i).Eta()), weight, 0., 2.5, 5);
           FillHist("FakeIDElSumW_Pt_1D_0l", FakeColl.at(i).Pt(), weight, 0., 200., 40);
           FillHist("FakeIDElSumW_HT_fEta_2D_0l", HT, fabs(FakeColl.at(i).Eta()), weight, 0., 1000., 10, 0., 2.5, 5);
         }


         for(int j=0; j<JetCleanColl.size(); j++){
           if(FakeColl.at(i).DeltaR(JetCleanColl.at(j))<0.4){
             FillHist("NCount_HasCloseJet_0l", 2., weight, 0., 3., 3);
             FillHist("FakeSumW_JetPt_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
             FillHist("FakeSumW_JetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);

             FillHist("dRej_fj_0l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
             FillHist("FakeRelPT_0l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
             FillHist("FakeRelPT_JetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
             FillHist("FakeRelPT_JetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
             FillHist("FakePT_JetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
             FillHist("FakePT_JetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);

             FillHist("FakeElAvgRelIso_JetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
             FillHist("FakeElAvgRelIso_JetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);

             if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
               FillHist("dRIDej_fj_0l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeIDSumW_JetPt_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeIDSumW_JetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
               FillHist("FakeIDPT_JetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
             }
             
             //Bjet
             if(JetCleanColl.at(j).HadronFlavour()==5){
             //if(JetCleanColl.at(j).BJetTaggerValue(snu::KJet::CSVv2)>0.8484){
               FillHist("FakeSumW_BJetPt_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeSumW_BJetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
  
               FillHist("dReBj_fBj_0l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeRelPT_relB_0l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
               FillHist("FakeRelPT_BJetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
               FillHist("FakeRelPT_BJetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
               FillHist("FakePT_BJetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               FillHist("FakePT_BJetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);
  
               FillHist("FakeElAvgRelIso_BJetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
               FillHist("FakeElAvgRelIso_BJetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);
  
               if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
                 FillHist("dRIDeBj_fBj_0l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
                 FillHist("FakeIDSumW_BJetPt_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("FakeIDSumW_BJetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
                 FillHist("FakeIDPT_BJetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               }
             }
             else if(JetCleanColl.at(j).HadronFlavour()==4){
               FillHist("FakeSumW_CJetPt_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeSumW_CJetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
  
               FillHist("dReCj_fCj_0l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeRelPT_relC_0l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
               FillHist("FakeRelPT_CJetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
               FillHist("FakeRelPT_CJetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
               FillHist("FakePT_CJetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               FillHist("FakePT_CJetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);
  
               FillHist("FakeElAvgRelIso_CJetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
               FillHist("FakeElAvgRelIso_CJetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);
  
               if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
                 FillHist("dRIDeCj_fCj_0l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
                 FillHist("FakeIDSumW_CJetPt_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("FakeIDSumW_CJetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
                 FillHist("FakeIDPT_CJetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               }
             }
             else if(JetCleanColl.at(j).HadronFlavour()==0){
               FillHist("FakeSumW_LJetPt_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeSumW_LJetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
  
               FillHist("dReLj_fLj_0l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeRelPT_relL_0l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
               FillHist("FakeRelPT_LJetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
               FillHist("FakeRelPT_LJetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
               FillHist("FakePT_LJetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               FillHist("FakePT_LJetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);
  
               FillHist("FakeElAvgRelIso_LJetPt_1D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
               FillHist("FakeElAvgRelIso_LJetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);
  

               //Parametrisation in HadFlav+PartonPt(+JetPt)
               int MatchedPartIdx=GenMatchedIdx(JetCleanColl.at(i),truthColl);
               bool PartHadConsist=IsJetConsistentPartonHadronMatch(JetCleanColl.at(i), truthColl, "Heavy");
               if(MatchedPartIdx>=0 && PartHadConsist){
                 FillHist("FakeSumW_LParPt_1D_0l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                 if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("FakeSumW_LJetPt_LParPt25_50_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("FakeSumW_LJetPt_LParPt50_100_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("FakeSumW_LJetPt_LParPt100_150_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("FakeSumW_LJetPt_LParPt150_200_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else                                           FillHist("FakeSumW_LJetPt_LParPt200_Inf_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               }

               //Parametrisation in HadFlav+NnearJet(+JetPt)
               FillHist("FakeSumW_LjMatch_NnearJet_1D_0l", NJetNearFakeLep.at(i), weight, 0., 10., 10);
               if(NJetNearFakeLep.at(i)==1) FillHist("FakeSumW_LJetPt_1nearj_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(NJetNearFakeLep.at(i)==2) FillHist("FakeSumW_LJetPt_2nearj_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(NJetNearFakeLep.at(i)==3) FillHist("FakeSumW_LJetPt_3nearj_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

               //Parametrisation in HadFlav+HT(+JetPt)
               FillHist("FakeSumW_LjMatch_HT_1D_0l", HT, weight, 0., 1000., 20);
               if(HT<100)      FillHist("FakeSumW_LJetPt_HT0_100_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(HT<200) FillHist("FakeSumW_LJetPt_HT100_200_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(HT<500) FillHist("FakeSumW_LJetPt_HT200_500_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else            FillHist("FakeSumW_LJetPt_HT500_Inf_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

                
               //Parametrisation in LeptonSource
               int LepType=GetLeptonType(FakeColl.at(i),truthColl);
               FillHist("FakeSumW_LjMatch_EleType_1D_0l", LepType, weight, -10., 10., 20);
               if(JetCleanColl.at(j).Pt()<50) FillHist("FakeSumW_LjMatch_EleType_LjPt25_50_1D_0l", LepType, weight, -10., 10., 20);
               if(LepType<0){
                 FillHist("HadFakeSumW_LJetPt_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("HadFakeSumW_LJetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);

                 if(MatchedPartIdx>=0 && PartHadConsist){
                   FillHist("HadFakeSumW_LParPt_1D_0l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                   if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("HadFakeSumW_LJetPt_LParPt25_50_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("HadFakeSumW_LJetPt_LParPt50_100_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("HadFakeSumW_LJetPt_LParPt100_150_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("HadFakeSumW_LJetPt_LParPt150_200_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else                                           FillHist("HadFakeSumW_LJetPt_LParPt200_Inf_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 } 

                 FillHist("HadFakeSumW_LjMatch_HT_1D_0l", HT, weight, 0., 1000., 20);
                 if(HT<100)      FillHist("HadFakeSumW_LJetPt_HT0_100_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<200) FillHist("HadFakeSumW_LJetPt_HT100_200_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<500) FillHist("HadFakeSumW_LJetPt_HT200_500_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else            FillHist("HadFakeSumW_LJetPt_HT500_Inf_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

               }


               if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
                 FillHist("dRIDeLj_fLj_0l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
                 FillHist("FakeIDSumW_LJetPt_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("FakeIDSumW_LJetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
                 FillHist("FakeIDPT_LJetPt_2D_0l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);


                 if(MatchedPartIdx>=0 && PartHadConsist){
                   FillHist("FakeIDSumW_LParPt_1D_0l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                   if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("FakeIDSumW_LJetPt_LParPt25_50_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("FakeIDSumW_LJetPt_LParPt50_100_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("FakeIDSumW_LJetPt_LParPt100_150_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("FakeIDSumW_LJetPt_LParPt150_200_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else                                           FillHist("FakeIDSumW_LJetPt_LParPt200_Inf_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 } 

                 //Parametrisation in HadFlav+NnearJet(+JetPt) -> Doesn't work!
                 FillHist("FakeIDSumW_LjMatch_NnearJet_1D_0l", NJetNearFakeLep.at(i), weight, 0., 10., 10);
                 if(NJetNearFakeLep.at(i)==1) FillHist("FakeIDSumW_LJetPt_1nearj_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(NJetNearFakeLep.at(i)==2) FillHist("FakeIDSumW_LJetPt_2nearj_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(NJetNearFakeLep.at(i)==3) FillHist("FakeIDSumW_LJetPt_3nearj_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
  
                 //Parametrisation in HadFlav+HT(+JetPt) ->Doesn't work! (0l,1l similar but cannot explain 2l!)
                 FillHist("FakeIDSumW_LjMatch_HT_1D_0l", HT, weight, 0., 1000., 20);
                 if(HT<100)      FillHist("FakeIDSumW_LJetPt_HT0_100_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<200) FillHist("FakeIDSumW_LJetPt_HT100_200_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<500) FillHist("FakeIDSumW_LJetPt_HT200_500_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else            FillHist("FakeIDSumW_LJetPt_HT500_Inf_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);


                 //Parametrisation in LeptonSource
                 FillHist("FakeIDSumW_LjMatch_EleType_1D_0l", LepType, weight, -10., 10., 20);
                 if(JetCleanColl.at(j).Pt()<50) FillHist("FakeIDSumW_LjMatch_EleType_LjPt25_50_1D_0l", LepType, weight, -10., 10., 20);

                 if(LepType<0){
                   FillHist("HadFakeIDSumW_LJetPt_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   FillHist("HadFakeIDSumW_LJetfEta_1D_0l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);

                   if(MatchedPartIdx>=0 && PartHadConsist){
                     FillHist("HadFakeIDSumW_LParPt_1D_0l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                     if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("HadFakeIDSumW_LJetPt_LParPt25_50_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("HadFakeIDSumW_LJetPt_LParPt50_100_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("HadFakeIDSumW_LJetPt_LParPt100_150_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("HadFakeIDSumW_LJetPt_LParPt150_200_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else                                           FillHist("HadFakeIDSumW_LJetPt_LParPt200_Inf_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   } 
  
                   FillHist("HadFakeIDSumW_LjMatch_HT_1D_0l", HT, weight, 0., 1000., 20);
                   if(HT<100)      FillHist("HadFakeIDSumW_LJetPt_HT0_100_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(HT<200) FillHist("HadFakeIDSumW_LJetPt_HT100_200_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(HT<500) FillHist("HadFakeIDSumW_LJetPt_HT200_500_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else            FillHist("HadFakeIDSumW_LJetPt_HT500_Inf_1D_0l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

                 }

               }//End of ID Pass
             }//End of Light HadronFlav Jet
           }//End of Jet-Fake Matched Case
         }//End of JetLoop
       }//End of FakeLoop
     }//End of Gen 0lep case
     else if(NPromptGenLepAll==1){
       FillHist("HT_1l", HT, weight, 0., 1000., 50);
       for(int i=0; i<PromptColl.size(); i++){
         FillHist("PromptElRelIso_1l", PromptColl.at(i).PFRelIso(0.3), weight, 0., 1., 20);

         FillHist("PromptElAvgRelIso_HT_1D_1l", HT, weight*PromptColl.at(i).PFRelIso(0.3), 0., 1000., 10);
         FillHist("PromptElSumW_HT_1D_1l", HT, weight, 0., 1000., 10);
         FillHist("PromptElAvgRelIso_fEta_1D_1l", fabs(PromptColl.at(i).Eta()), weight*PromptColl.at(i).PFRelIso(0.3), 0., 2.5, 5);
         FillHist("PromptElSumW_fEta_1D_1l", fabs(PromptColl.at(i).Eta()), weight, 0., 2.5, 5);

         FillHist("PromptPT_PromptfEta_1D_1l", fabs(PromptColl.at(i).Eta()), PromptColl.at(i).Pt()*weight, 0., 2.5, 5);
         FillHist("PromptPT_PromptfEta_2D_1l", PromptColl.at(i).Pt(), fabs(PromptColl.at(i).Eta()), weight, 0., 300., 60, 0., 2.5, 5);
       }
       for(int i=0; i<FakeColl.size(); i++){
         FillHist("NCount_HasCloseJet_1l", 0., weight, 0., 3., 3);
         if(NearEWLep(FakeColl.at(i), truthColl, "InclTau")) FillHist("NCount_HasCloseJet_1l", 1., weight, 0., 3., 3);

         FillHist("FakeElSumW_HT_1D_1l", HT, weight, 0., 1000., 10);
         FillHist("FakeElSumW_fEta_1D_1l", fabs(FakeColl.at(i).Eta()), weight, 0., 2.5, 5);
         FillHist("FakeElSumW_Pt_1D_1l", FakeColl.at(i).Pt(), weight, 0., 200., 40);
         FillHist("FakeElSumW_HT_fEta_2D_1l", HT, fabs(FakeColl.at(i).Eta()), weight, 0., 1000., 10, 0., 2.5, 5);

         FillHist("FakeElRelIso_1l", FakeColl.at(i).PFRelIso(0.3), weight, 0., 1., 20);
         FillHist("FakeElAvgRelIso_HT_1D_1l", HT, weight*FakeColl.at(i).PFRelIso(0.3), 0., 1000., 10);
         FillHist("FakeElAvgRelIso_fEta_1D_1l", fabs(FakeColl.at(i).Eta()), weight*FakeColl.at(i).PFRelIso(0.3), 0., 2.5, 5);

         FillHist("FakePT_FakefEta_1D_1l", fabs(FakeColl.at(i).Eta()), FakeColl.at(i).Pt()*weight, 0., 2.5, 5);
         FillHist("FakePT_FakefEta_2D_1l", FakeColl.at(i).Pt(), fabs(FakeColl.at(i).Eta()), weight, 0., 200, 40, 0., 2.5, 5);

         if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
           FillHist("FakeIDElSumW_HT_1D_1l", HT, weight, 0., 1000., 10);
           FillHist("FakeIDElSumW_fEta_1D_1l", fabs(FakeColl.at(i).Eta()), weight, 0., 2.5, 5);
           FillHist("FakeIDElSumW_Pt_1D_1l", FakeColl.at(i).Pt(), weight, 0., 200., 40);
           FillHist("FakeIDElSumW_HT_fEta_2D_1l", HT, fabs(FakeColl.at(i).Eta()), weight, 0., 1000., 10, 0., 2.5, 5);
         }


         for(int j=0; j<JetCleanColl.size(); j++){
           if(FakeColl.at(i).DeltaR(JetCleanColl.at(j))<0.4){
             FillHist("NCount_HasCloseJet_1l", 2., weight, 0., 3., 3);
             FillHist("FakeSumW_JetPt_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
             FillHist("FakeSumW_JetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);

             FillHist("dRej_fj_1l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
             FillHist("FakeRelPT_1l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
             FillHist("FakeRelPT_JetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
             FillHist("FakeRelPT_JetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
             FillHist("FakePT_JetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
             FillHist("FakePT_JetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);

             FillHist("FakeElAvgRelIso_JetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
             FillHist("FakeElAvgRelIso_JetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);

             if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
               FillHist("dRIDej_fj_1l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeIDSumW_JetPt_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeIDSumW_JetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
               FillHist("FakeIDPT_JetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
             }
             
             //Bjet
             if(JetCleanColl.at(j).HadronFlavour()==5){
             //if(JetCleanColl.at(j).BJetTaggerValue(snu::KJet::CSVv2)>0.8484){
               FillHist("FakeSumW_BJetPt_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeSumW_BJetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
  
               FillHist("dReBj_fBj_1l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeRelPT_relB_1l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
               FillHist("FakeRelPT_BJetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
               FillHist("FakeRelPT_BJetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
               FillHist("FakePT_BJetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               FillHist("FakePT_BJetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);
  
               FillHist("FakeElAvgRelIso_BJetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
               FillHist("FakeElAvgRelIso_BJetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);
  
               if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
                 FillHist("dRIDeBj_fBj_1l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
                 FillHist("FakeIDSumW_BJetPt_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("FakeIDSumW_BJetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
                 FillHist("FakeIDPT_BJetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               }
             }
             else if(JetCleanColl.at(j).HadronFlavour()==4){
               FillHist("FakeSumW_CJetPt_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeSumW_CJetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
  
               FillHist("dReCj_fCj_1l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeRelPT_relC_1l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
               FillHist("FakeRelPT_CJetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
               FillHist("FakeRelPT_CJetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
               FillHist("FakePT_CJetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               FillHist("FakePT_CJetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);
  
               FillHist("FakeElAvgRelIso_CJetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
               FillHist("FakeElAvgRelIso_CJetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);
  
               if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
                 FillHist("dRIDeCj_fCj_1l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
                 FillHist("FakeIDSumW_CJetPt_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("FakeIDSumW_CJetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
                 FillHist("FakeIDPT_CJetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               }
             }
             else if(JetCleanColl.at(j).HadronFlavour()==0){
               FillHist("FakeSumW_LJetPt_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeSumW_LJetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
  
               FillHist("dReLj_fLj_1l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeRelPT_relL_1l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
               FillHist("FakeRelPT_LJetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
               FillHist("FakeRelPT_LJetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
               FillHist("FakePT_LJetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               FillHist("FakePT_LJetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);
  
               FillHist("FakeElAvgRelIso_LJetPt_1D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
               FillHist("FakeElAvgRelIso_LJetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);
  

               //Parametrisation in HadFlav+PartonPt(+JetPt)
               int MatchedPartIdx=GenMatchedIdx(JetCleanColl.at(i),truthColl);
               bool PartHadConsist=IsJetConsistentPartonHadronMatch(JetCleanColl.at(i), truthColl, "Heavy");
               if(MatchedPartIdx>=0 && PartHadConsist){
                 FillHist("FakeSumW_LParPt_1D_1l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                 if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("FakeSumW_LJetPt_LParPt25_50_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("FakeSumW_LJetPt_LParPt50_100_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("FakeSumW_LJetPt_LParPt100_150_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("FakeSumW_LJetPt_LParPt150_200_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else                                           FillHist("FakeSumW_LJetPt_LParPt200_Inf_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               }

               //Parametrisation in HadFlav+NnearJet(+JetPt)
               FillHist("FakeSumW_LjMatch_NnearJet_1D_1l", NJetNearFakeLep.at(i), weight, 0., 10., 10);
               if(NJetNearFakeLep.at(i)==1) FillHist("FakeSumW_LJetPt_1nearj_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(NJetNearFakeLep.at(i)==2) FillHist("FakeSumW_LJetPt_2nearj_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(NJetNearFakeLep.at(i)==3) FillHist("FakeSumW_LJetPt_3nearj_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

               //Parametrisation in HadFlav+HT(+JetPt)
               FillHist("FakeSumW_LjMatch_HT_1D_1l", HT, weight, 0., 1000., 20);
               if(HT<100)      FillHist("FakeSumW_LJetPt_HT0_100_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(HT<200) FillHist("FakeSumW_LJetPt_HT100_200_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(HT<500) FillHist("FakeSumW_LJetPt_HT200_500_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else            FillHist("FakeSumW_LJetPt_HT500_Inf_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

                
               //Parametrisation in LeptonSource
               int LepType=GetLeptonType(FakeColl.at(i),truthColl);
               FillHist("FakeSumW_LjMatch_EleType_1D_1l", LepType, weight, -10., 10., 20);
               if(JetCleanColl.at(j).Pt()<50) FillHist("FakeSumW_LjMatch_EleType_LjPt25_50_1D_1l", LepType, weight, -10., 10., 20);
               if(LepType<0){
                 FillHist("HadFakeSumW_LJetPt_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("HadFakeSumW_LJetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);

                 if(MatchedPartIdx>=0 && PartHadConsist){
                   FillHist("HadFakeSumW_LParPt_1D_1l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                   if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("HadFakeSumW_LJetPt_LParPt25_50_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("HadFakeSumW_LJetPt_LParPt50_100_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("HadFakeSumW_LJetPt_LParPt100_150_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("HadFakeSumW_LJetPt_LParPt150_200_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else                                           FillHist("HadFakeSumW_LJetPt_LParPt200_Inf_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 } 

                 FillHist("HadFakeSumW_LjMatch_HT_1D_1l", HT, weight, 0., 1000., 20);
                 if(HT<100)      FillHist("HadFakeSumW_LJetPt_HT0_100_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<200) FillHist("HadFakeSumW_LJetPt_HT100_200_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<500) FillHist("HadFakeSumW_LJetPt_HT200_500_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else            FillHist("HadFakeSumW_LJetPt_HT500_Inf_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

               }


               if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
                 FillHist("dRIDeLj_fLj_1l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
                 FillHist("FakeIDSumW_LJetPt_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("FakeIDSumW_LJetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
                 FillHist("FakeIDPT_LJetPt_2D_1l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);


                 if(MatchedPartIdx>=0 && PartHadConsist){
                   FillHist("FakeIDSumW_LParPt_1D_1l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                   if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("FakeIDSumW_LJetPt_LParPt25_50_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("FakeIDSumW_LJetPt_LParPt50_100_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("FakeIDSumW_LJetPt_LParPt100_150_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("FakeIDSumW_LJetPt_LParPt150_200_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else                                           FillHist("FakeIDSumW_LJetPt_LParPt200_Inf_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 } 

                 //Parametrisation in HadFlav+NnearJet(+JetPt) -> Doesn't work!
                 FillHist("FakeIDSumW_LjMatch_NnearJet_1D_1l", NJetNearFakeLep.at(i), weight, 0., 10., 10);
                 if(NJetNearFakeLep.at(i)==1) FillHist("FakeIDSumW_LJetPt_1nearj_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(NJetNearFakeLep.at(i)==2) FillHist("FakeIDSumW_LJetPt_2nearj_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(NJetNearFakeLep.at(i)==3) FillHist("FakeIDSumW_LJetPt_3nearj_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
  
                 //Parametrisation in HadFlav+HT(+JetPt) ->Doesn't work! (1l,1l similar but cannot explain 2l!)
                 FillHist("FakeIDSumW_LjMatch_HT_1D_1l", HT, weight, 0., 1000., 20);
                 if(HT<100)      FillHist("FakeIDSumW_LJetPt_HT0_100_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<200) FillHist("FakeIDSumW_LJetPt_HT100_200_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<500) FillHist("FakeIDSumW_LJetPt_HT200_500_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else            FillHist("FakeIDSumW_LJetPt_HT500_Inf_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);


                 //Parametrisation in LeptonSource
                 FillHist("FakeIDSumW_LjMatch_EleType_1D_1l", LepType, weight, -10., 10., 20);
                 if(JetCleanColl.at(j).Pt()<50) FillHist("FakeIDSumW_LjMatch_EleType_LjPt25_50_1D_1l", LepType, weight, -10., 10., 20);

                 if(LepType<0){
                   FillHist("HadFakeIDSumW_LJetPt_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   FillHist("HadFakeIDSumW_LJetfEta_1D_1l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);

                   if(MatchedPartIdx>=0 && PartHadConsist){
                     FillHist("HadFakeIDSumW_LParPt_1D_1l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                     if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("HadFakeIDSumW_LJetPt_LParPt25_50_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("HadFakeIDSumW_LJetPt_LParPt50_100_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("HadFakeIDSumW_LJetPt_LParPt100_150_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("HadFakeIDSumW_LJetPt_LParPt150_200_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else                                           FillHist("HadFakeIDSumW_LJetPt_LParPt200_Inf_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   } 
  
                   FillHist("HadFakeIDSumW_LjMatch_HT_1D_1l", HT, weight, 0., 1000., 20);
                   if(HT<100)      FillHist("HadFakeIDSumW_LJetPt_HT0_100_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(HT<200) FillHist("HadFakeIDSumW_LJetPt_HT100_200_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(HT<500) FillHist("HadFakeIDSumW_LJetPt_HT200_500_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else            FillHist("HadFakeIDSumW_LJetPt_HT500_Inf_1D_1l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

                 }

               }//End of ID Pass
             }//End of Light HadronFlav Jet
           }//End of Jet-Fake Matched Case
         }//End of JetLoop
       }//End of FakeLoop

     }//End of Gen 1l case
     else if(NPromptGenLepAll==2){
       FillHist("HT_2l", HT, weight, 0., 1000., 50);
       for(int i=0; i<PromptColl.size(); i++){
         FillHist("PromptElRelIso_2l", PromptColl.at(i).PFRelIso(0.3), weight, 0., 1., 20);

         FillHist("PromptElAvgRelIso_HT_1D_2l", HT, weight*PromptColl.at(i).PFRelIso(0.3), 0., 1000., 10);
         FillHist("PromptElSumW_HT_1D_2l", HT, weight, 0., 1000., 10);
         FillHist("PromptElAvgRelIso_fEta_1D_2l", fabs(PromptColl.at(i).Eta()), weight*PromptColl.at(i).PFRelIso(0.3), 0., 2.5, 5);
         FillHist("PromptElSumW_fEta_1D_2l", fabs(PromptColl.at(i).Eta()), weight, 0., 2.5, 5);

         FillHist("PromptPT_PromptfEta_1D_2l", fabs(PromptColl.at(i).Eta()), PromptColl.at(i).Pt()*weight, 0., 2.5, 5);
         FillHist("PromptPT_PromptfEta_2D_2l", PromptColl.at(i).Pt(), fabs(PromptColl.at(i).Eta()), weight, 0., 300., 60, 0., 2.5, 5);
       }
       for(int i=0; i<FakeColl.size(); i++){
         FillHist("NCount_HasCloseJet_2l", 0., weight, 0., 3., 3);
         if(NearEWLep(FakeColl.at(i), truthColl, "InclTau")) FillHist("NCount_HasCloseJet_2l", 1., weight, 0., 3., 3);

         FillHist("FakeElSumW_HT_1D_2l", HT, weight, 0., 1000., 10);
         FillHist("FakeElSumW_fEta_1D_2l", fabs(FakeColl.at(i).Eta()), weight, 0., 2.5, 5);
         FillHist("FakeElSumW_Pt_1D_2l", FakeColl.at(i).Pt(), weight, 0., 200., 40);
         FillHist("FakeElSumW_HT_fEta_2D_2l", HT, fabs(FakeColl.at(i).Eta()), weight, 0., 1000., 10, 0., 2.5, 5);

         FillHist("FakeElRelIso_2l", FakeColl.at(i).PFRelIso(0.3), weight, 0., 1., 20);
         FillHist("FakeElAvgRelIso_HT_1D_2l", HT, weight*FakeColl.at(i).PFRelIso(0.3), 0., 1000., 10);
         FillHist("FakeElAvgRelIso_fEta_1D_2l", fabs(FakeColl.at(i).Eta()), weight*FakeColl.at(i).PFRelIso(0.3), 0., 2.5, 5);

         FillHist("FakePT_FakefEta_1D_2l", fabs(FakeColl.at(i).Eta()), FakeColl.at(i).Pt()*weight, 0., 2.5, 5);
         FillHist("FakePT_FakefEta_2D_2l", FakeColl.at(i).Pt(), fabs(FakeColl.at(i).Eta()), weight, 0., 200, 40, 0., 2.5, 5);

         if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
           FillHist("FakeIDElSumW_HT_1D_2l", HT, weight, 0., 1000., 10);
           FillHist("FakeIDElSumW_fEta_1D_2l", fabs(FakeColl.at(i).Eta()), weight, 0., 2.5, 5);
           FillHist("FakeIDElSumW_Pt_1D_2l", FakeColl.at(i).Pt(), weight, 0., 200., 40);
           FillHist("FakeIDElSumW_HT_fEta_2D_2l", HT, fabs(FakeColl.at(i).Eta()), weight, 0., 1000., 10, 0., 2.5, 5);
         }


         for(int j=0; j<JetCleanColl.size(); j++){
           if(FakeColl.at(i).DeltaR(JetCleanColl.at(j))<0.4){
             FillHist("NCount_HasCloseJet_2l", 2., weight, 0., 3., 3);
             FillHist("FakeSumW_JetPt_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
             FillHist("FakeSumW_JetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);

             FillHist("dRej_fj_2l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
             FillHist("FakeRelPT_2l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
             FillHist("FakeRelPT_JetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
             FillHist("FakeRelPT_JetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
             FillHist("FakePT_JetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
             FillHist("FakePT_JetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);

             FillHist("FakeElAvgRelIso_JetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
             FillHist("FakeElAvgRelIso_JetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);

             if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
               FillHist("dRIDej_fj_2l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeIDSumW_JetPt_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeIDSumW_JetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
               FillHist("FakeIDPT_JetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
             }
             
             //Bjet
             if(JetCleanColl.at(j).HadronFlavour()==5){
             //if(JetCleanColl.at(j).BJetTaggerValue(snu::KJet::CSVv2)>0.8484){
               FillHist("FakeSumW_BJetPt_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeSumW_BJetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
  
               FillHist("dReBj_fBj_2l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeRelPT_relB_2l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
               FillHist("FakeRelPT_BJetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
               FillHist("FakeRelPT_BJetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
               FillHist("FakePT_BJetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               FillHist("FakePT_BJetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);
  
               FillHist("FakeElAvgRelIso_BJetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
               FillHist("FakeElAvgRelIso_BJetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);
  
               if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
                 FillHist("dRIDeBj_fBj_2l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
                 FillHist("FakeIDSumW_BJetPt_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("FakeIDSumW_BJetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
                 FillHist("FakeIDPT_BJetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               }
             }
             else if(JetCleanColl.at(j).HadronFlavour()==4){
               FillHist("FakeSumW_CJetPt_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeSumW_CJetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
  
               FillHist("dReCj_fCj_2l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeRelPT_relC_2l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
               FillHist("FakeRelPT_CJetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
               FillHist("FakeRelPT_CJetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
               FillHist("FakePT_CJetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               FillHist("FakePT_CJetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);
  
               FillHist("FakeElAvgRelIso_CJetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
               FillHist("FakeElAvgRelIso_CJetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);
  
               if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
                 FillHist("dRIDeCj_fCj_2l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
                 FillHist("FakeIDSumW_CJetPt_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("FakeIDSumW_CJetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
                 FillHist("FakeIDPT_CJetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               }
             }
             else if(JetCleanColl.at(j).HadronFlavour()==0){
               FillHist("FakeSumW_LJetPt_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               FillHist("FakeSumW_LJetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
  
               FillHist("dReLj_fLj_2l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
               FillHist("FakeRelPT_relL_2l", FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 1., 20);
               FillHist("FakeRelPT_LJetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt()*weight, 0., 500., 20);
               FillHist("FakeRelPT_LJetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()/JetCleanColl.at(j).Pt(), weight, 0., 500., 20, 0., 1., 20);
               FillHist("FakePT_LJetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);
               FillHist("FakePT_LJetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt()*weight, 0., 500., 20);
  
               FillHist("FakeElAvgRelIso_LJetPt_1D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 500., 20);
               FillHist("FakeElAvgRelIso_LJetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), FakeColl.at(i).PFRelIso(0.3)*weight, 0., 2.5, 5);
  

               //Parametrisation in HadFlav+PartonPt(+JetPt)
               int MatchedPartIdx=GenMatchedIdx(JetCleanColl.at(i),truthColl);
               bool PartHadConsist=IsJetConsistentPartonHadronMatch(JetCleanColl.at(i), truthColl, "Heavy");
               if(MatchedPartIdx>=0 && PartHadConsist){
                 FillHist("FakeSumW_LParPt_1D_2l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                 if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("FakeSumW_LJetPt_LParPt25_50_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("FakeSumW_LJetPt_LParPt50_100_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("FakeSumW_LJetPt_LParPt100_150_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("FakeSumW_LJetPt_LParPt150_200_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else                                           FillHist("FakeSumW_LJetPt_LParPt200_Inf_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               }

               //Parametrisation in HadFlav+NnearJet(+JetPt)
               FillHist("FakeSumW_LjMatch_NnearJet_1D_2l", NJetNearFakeLep.at(i), weight, 0., 10., 10);
               if(NJetNearFakeLep.at(i)==1) FillHist("FakeSumW_LJetPt_1nearj_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(NJetNearFakeLep.at(i)==2) FillHist("FakeSumW_LJetPt_2nearj_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(NJetNearFakeLep.at(i)==3) FillHist("FakeSumW_LJetPt_3nearj_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

               //Parametrisation in HadFlav+HT(+JetPt)
               FillHist("FakeSumW_LjMatch_HT_1D_2l", HT, weight, 0., 1000., 20);
               if(HT<100)      FillHist("FakeSumW_LJetPt_HT0_100_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(HT<200) FillHist("FakeSumW_LJetPt_HT100_200_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else if(HT<500) FillHist("FakeSumW_LJetPt_HT200_500_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
               else            FillHist("FakeSumW_LJetPt_HT500_Inf_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

                
               //Parametrisation in LeptonSource
               int LepType=GetLeptonType(FakeColl.at(i),truthColl);
               FillHist("FakeSumW_LjMatch_EleType_1D_2l", LepType, weight, -10., 10., 20);
               if(JetCleanColl.at(j).Pt()<50) FillHist("FakeSumW_LjMatch_EleType_LjPt25_50_1D_2l", LepType, weight, -10., 10., 20);
               if(LepType<0){
                 FillHist("HadFakeSumW_LJetPt_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("HadFakeSumW_LJetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);

                 if(MatchedPartIdx>=0 && PartHadConsist){
                   FillHist("HadFakeSumW_LParPt_1D_2l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                   if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("HadFakeSumW_LJetPt_LParPt25_50_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("HadFakeSumW_LJetPt_LParPt50_100_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("HadFakeSumW_LJetPt_LParPt100_150_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("HadFakeSumW_LJetPt_LParPt150_200_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else                                           FillHist("HadFakeSumW_LJetPt_LParPt200_Inf_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 } 

                 FillHist("HadFakeSumW_LjMatch_HT_1D_2l", HT, weight, 0., 1000., 20);
                 if(HT<100)      FillHist("HadFakeSumW_LJetPt_HT0_100_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<200) FillHist("HadFakeSumW_LJetPt_HT100_200_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<500) FillHist("HadFakeSumW_LJetPt_HT200_500_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else            FillHist("HadFakeSumW_LJetPt_HT500_Inf_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

               }


               if(FakeColl.at(i).PFRelIso(0.3)<0.0588){
                 FillHist("dRIDeLj_fLj_2l", FakeColl.at(i).DeltaR(JetCleanColl.at(j)), weight, 0., 0.4, 40);
                 FillHist("FakeIDSumW_LJetPt_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 FillHist("FakeIDSumW_LJetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);
                 FillHist("FakeIDPT_LJetPt_2D_2l", JetCleanColl.at(j).Pt(), FakeColl.at(i).Pt(), weight, 0., 500., 20, 0., 200., 20);


                 if(MatchedPartIdx>=0 && PartHadConsist){
                   FillHist("FakeIDSumW_LParPt_1D_2l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                   if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("FakeIDSumW_LJetPt_LParPt25_50_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("FakeIDSumW_LJetPt_LParPt50_100_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("FakeIDSumW_LJetPt_LParPt100_150_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("FakeIDSumW_LJetPt_LParPt150_200_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else                                           FillHist("FakeIDSumW_LJetPt_LParPt200_Inf_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 } 

                 //Parametrisation in HadFlav+NnearJet(+JetPt) -> Doesn't work!
                 FillHist("FakeIDSumW_LjMatch_NnearJet_1D_2l", NJetNearFakeLep.at(i), weight, 0., 10., 10);
                 if(NJetNearFakeLep.at(i)==1) FillHist("FakeIDSumW_LJetPt_1nearj_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(NJetNearFakeLep.at(i)==2) FillHist("FakeIDSumW_LJetPt_2nearj_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(NJetNearFakeLep.at(i)==3) FillHist("FakeIDSumW_LJetPt_3nearj_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
  
                 //Parametrisation in HadFlav+HT(+JetPt) ->Doesn't work! (2l,1l similar but cannot explain 2l!)
                 FillHist("FakeIDSumW_LjMatch_HT_1D_2l", HT, weight, 0., 1000., 20);
                 if(HT<100)      FillHist("FakeIDSumW_LJetPt_HT0_100_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<200) FillHist("FakeIDSumW_LJetPt_HT100_200_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else if(HT<500) FillHist("FakeIDSumW_LJetPt_HT200_500_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                 else            FillHist("FakeIDSumW_LJetPt_HT500_Inf_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);


                 //Parametrisation in LeptonSource
                 FillHist("FakeIDSumW_LjMatch_EleType_1D_2l", LepType, weight, -10., 10., 20);
                 if(JetCleanColl.at(j).Pt()<50) FillHist("FakeIDSumW_LjMatch_EleType_LjPt25_50_1D_2l", LepType, weight, -10., 10., 20);

                 if(LepType<0){
                   FillHist("HadFakeIDSumW_LJetPt_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   FillHist("HadFakeIDSumW_LJetfEta_1D_2l", fabs(JetCleanColl.at(j).Eta()), weight, 0., 2.5, 5);

                   if(MatchedPartIdx>=0 && PartHadConsist){
                     FillHist("HadFakeIDSumW_LParPt_1D_2l", truthColl.at(MatchedPartIdx).Pt(), weight, 0., 500., 20);
                     if(truthColl.at(MatchedPartIdx).Pt()<50)       FillHist("HadFakeIDSumW_LJetPt_LParPt25_50_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else if(truthColl.at(MatchedPartIdx).Pt()<100) FillHist("HadFakeIDSumW_LJetPt_LParPt50_100_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else if(truthColl.at(MatchedPartIdx).Pt()<150) FillHist("HadFakeIDSumW_LJetPt_LParPt100_150_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else if(truthColl.at(MatchedPartIdx).Pt()<200) FillHist("HadFakeIDSumW_LJetPt_LParPt150_200_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                     else                                           FillHist("HadFakeIDSumW_LJetPt_LParPt200_Inf_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   } 
  
                   FillHist("HadFakeIDSumW_LjMatch_HT_1D_2l", HT, weight, 0., 1000., 20);
                   if(HT<100)      FillHist("HadFakeIDSumW_LJetPt_HT0_100_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(HT<200) FillHist("HadFakeIDSumW_LJetPt_HT100_200_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else if(HT<500) FillHist("HadFakeIDSumW_LJetPt_HT200_500_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);
                   else            FillHist("HadFakeIDSumW_LJetPt_HT500_Inf_1D_2l", JetCleanColl.at(j).Pt(), weight, 0., 500., 20);

                 }

               }//End of ID Pass
             }//End of Light HadronFlav Jet
           }//End of Jet-Fake Matched Case
         }//End of JetLoop
       }//End of FakeLoop

     }//End of Gen 2l case
   }//End of EleFakeStudy
   else if(BasicCompositionStudy){
     //std::vector<snu::KElectron> FakeColl     = SkimLepColl(electronFakeLColl, truthColl, "Fake");
     //std::vector<snu::KElectron> PromptColl   = SkimLepColl(electronFakeLColl, truthColl, "Prompt");
     std::vector<snu::KJet>      JetCleanColl = SkimJetColl(jetColl, truthColl, "NoPrNoTau");
     //int NPromptGenLepAll = NPromptLeptons(truthColl,"InclTau");
     //int NConversion      = NPromptLeptons(truthColl,"OnlyConv");

    //What To Check
    //1) How many jets have PFMu/GsfEle & POGLoose/POGLooseEle?
    //   - What happens if lep included jets are not corrected?
    //2) if jets are close to lepton, how many leptons do they have in average?(Simplest/Loose/Tight)
    //3) Does it have dependence on flavour? or source?(LepInclJetFr & AvgNLepInJet)
    //4) Lepton Type of Fake leptons in Jet?
    //   (Is it true that light jet fakes are midid and heavy jet fakes are secondary lepton?)
     for(int i=0; i<JetCleanColl.size(); i++){
       
       bool IsBJet=JetCleanColl.at(i).HadronFlavour()==5? true:false;
       bool IsCJet=JetCleanColl.at(i).HadronFlavour()==4? true:false;
       bool IsLJet=JetCleanColl.at(i).HadronFlavour()==0? true:false;

       bool EleInJet=false,  POGEleInJet=false,  MuInJet=false,  POGMuInJet=false;
       bool EleInBJet=false, POGEleInBJet=false, MuInBJet=false, POGMuInBJet=false;
       bool EleInCJet=false, POGEleInCJet=false, MuInCJet=false, POGMuInCJet=false;
       bool EleInLJet=false, POGEleInLJet=false, MuInLJet=false, POGMuInLJet=false;

       int NEleInThisJet=0,  NPOGEleInThisJet=0,  NMuInThisJet=0,  NPOGMuInThisJet=0;
       int NEleInThisBJet=0, NPOGEleInThisBJet=0, NMuInThisBJet=0, NPOGMuInThisBJet=0;
       int NEleInThisCJet=0, NPOGEleInThisCJet=0, NMuInThisCJet=0, NPOGMuInThisCJet=0;
       int NEleInThisLJet=0, NPOGEleInThisLJet=0, NMuInThisLJet=0, NPOGMuInThisLJet=0;

       int EleType=0, MuType=0;
       for(int j=0; j<electronPreColl.size(); j++){

         if(JetCleanColl.at(i).DeltaR(electronPreColl.at(j))<0.4){
           EleType=GetLeptonType(electronPreColl.at(j), truthColl);
           EleInJet=true;  NEleInThisJet++;
           FillHist("EleTypeInJet", EleType, weight, -5., 5., 10);
           if(electronPreColl.at(j).PassLoose() && electronPreColl.at(j).PFRelIso(0.3)<0.1){
             POGEleInJet=true; NPOGEleInThisJet++;
             FillHist("POGEleTypeInJet", EleType, weight, -5., 5., 10);
           }

           if(IsBJet){
             EleInBJet=true; NEleInThisBJet++;
             FillHist("EleTypeInBJet", EleType, weight, -5., 5., 10);
             if(electronPreColl.at(j).PassLoose() && electronPreColl.at(j).PFRelIso(0.3)<0.1){
               POGEleInBJet=true; NPOGEleInThisBJet++;
               FillHist("POGEleTypeInBJet", EleType, weight, -5., 5., 10);
             }
           }
           else if(IsCJet){
             EleInCJet=true; NEleInThisCJet++;
             FillHist("EleTypeInCJet", EleType, weight, -5., 5., 10);
             if(electronPreColl.at(j).PassLoose() && electronPreColl.at(j).PFRelIso(0.3)<0.1){
               POGEleInCJet=true; NPOGEleInThisCJet++;
               FillHist("POGEleTypeInCJet", EleType, weight, -5., 5., 10);
             }
           }
           else if(IsLJet){
             EleInLJet=true; NEleInThisLJet++;
             FillHist("EleTypeInLJet", EleType, weight, -5., 5., 10);
             if(electronPreColl.at(j).PassLoose() && electronPreColl.at(j).PFRelIso(0.3)<0.1){
               POGEleInLJet=true; NPOGEleInThisLJet++;
               FillHist("POGEleTypeInLJet", EleType, weight, -5., 5., 10);
             }
           }

         }//Ele-Jet Matched dR04
       }//End of Ele Loop


       for(int j=0; j<muonPreColl.size(); j++){

         if(JetCleanColl.at(i).DeltaR(muonPreColl.at(j))<0.4){
           MuType=GetLeptonType(muonPreColl.at(j), truthColl);
           MuInJet=true;  NMuInThisJet++;
           FillHist("MuTypeInJet", MuType, weight, -5., 5., 10);
           if(muonPreColl.at(j).IsLoose() && muonPreColl.at(j).RelIso04()<0.25){
             POGMuInJet=true; NPOGMuInThisJet++;
             FillHist("POGMuTypeInJet", MuType, weight, -5., 5., 10);
           }

           if(IsBJet){
             MuInBJet=true; NMuInThisBJet++;
             FillHist("MuTypeInBJet", MuType, weight, -5., 5., 10);
             if(muonPreColl.at(j).IsLoose() && muonPreColl.at(j).RelIso04()<0.25){
               POGMuInBJet=true; NPOGMuInThisBJet++;
               FillHist("POGMuTypeInBJet", MuType, weight, -5., 5., 10);
             }
           }
           else if(IsCJet){
             MuInCJet=true; NMuInThisCJet++;
             FillHist("MuTypeInCJet", MuType, weight, -5., 5., 10);
             if(muonPreColl.at(j).IsLoose() && muonPreColl.at(j).RelIso04()<0.25){
               POGMuInCJet=true; NPOGMuInThisCJet++;
               FillHist("POGMuTypeInCJet", MuType, weight, -5., 5., 10);
             }
           }
           else if(IsLJet){
             MuInLJet=true; NMuInThisLJet++;
             FillHist("MuTypeInLJet", MuType, weight, -5., 5., 10);
             if(muonPreColl.at(j).IsLoose() && muonPreColl.at(j).RelIso04()<0.25){
               POGMuInLJet=true; NPOGMuInThisLJet++;
               FillHist("POGMuTypeInLJet", MuType, weight, -5., 5., 10);
             }
           }

         }//Mu-Jet Matched dR04
       }//End of Mu Loop


       //Inclusive
       FillHist("LepInJetFr", 0., weight, 0., 7., 7);
       if(MuInJet || EleInJet)       FillHist("LepInJetFr", 1., weight, 0., 7., 7);
       if(EleInJet)                  FillHist("LepInJetFr", 2., weight, 0., 7., 7);
       if(MuInJet)                   FillHist("LepInJetFr", 3., weight, 0., 7., 7);
       if(POGEleInJet || POGMuInJet) FillHist("LepInJetFr", 4., weight, 0., 7., 7);
       if(POGEleInJet)               FillHist("LepInJetFr", 5., weight, 0., 7., 7);
       if(POGMuInJet)                FillHist("LepInJetFr", 6., weight, 0., 7., 7);

       FillHist("NEleInJet", NEleInThisJet, weight, 0., 10., 10);
       FillHist("NMuInJet", NMuInThisJet, weight, 0., 10., 10);
       FillHist("NPOGEleInJet", NPOGEleInThisJet, weight, 0., 10., 10);
       FillHist("NPOGMuInJet", NPOGMuInThisJet, weight, 0., 10., 10);

       //BJet
       if(IsBJet){
         FillHist("LepInBJetFr", 0., weight, 0., 7., 7);
         if(MuInBJet || EleInBJet)       FillHist("LepInBJetFr", 1., weight, 0., 7., 7);
         if(EleInBJet)                   FillHist("LepInBJetFr", 2., weight, 0., 7., 7);
         if(MuInBJet)                    FillHist("LepInBJetFr", 3., weight, 0., 7., 7);
         if(POGEleInBJet || POGMuInBJet) FillHist("LepInBJetFr", 4., weight, 0., 7., 7);
         if(POGEleInBJet)                FillHist("LepInBJetFr", 5., weight, 0., 7., 7);
         if(POGMuInBJet)                 FillHist("LepInBJetFr", 6., weight, 0., 7., 7);
  
         FillHist("NEleInBJet", NEleInThisBJet, weight, 0., 10., 10);
         FillHist("NMuInBJet", NMuInThisBJet, weight, 0., 10., 10);
         FillHist("NPOGEleInBJet", NPOGEleInThisBJet, weight, 0., 10., 10);
         FillHist("NPOGMuInBJet", NPOGMuInThisBJet, weight, 0., 10., 10);
       }
       //CJet
       else if(IsCJet){
         FillHist("LepInCJetFr", 0., weight, 0., 7., 7);
         if(MuInCJet || EleInCJet)       FillHist("LepInCJetFr", 1., weight, 0., 7., 7);
         if(EleInCJet)                   FillHist("LepInCJetFr", 2., weight, 0., 7., 7);
         if(MuInCJet)                    FillHist("LepInCJetFr", 3., weight, 0., 7., 7);
         if(POGEleInCJet || POGMuInCJet) FillHist("LepInCJetFr", 4., weight, 0., 7., 7);
         if(POGEleInCJet)                FillHist("LepInCJetFr", 5., weight, 0., 7., 7);
         if(POGMuInCJet)                 FillHist("LepInCJetFr", 6., weight, 0., 7., 7);
  
         FillHist("NEleInCJet", NEleInThisCJet, weight, 0., 10., 10);
         FillHist("NMuInCJet", NMuInThisCJet, weight, 0., 10., 10);
         FillHist("NPOGEleInCJet", NPOGEleInThisCJet, weight, 0., 10., 10);
         FillHist("NPOGMuInCJet", NPOGMuInThisCJet, weight, 0., 10., 10);
       }
       //LJet
       else if(IsLJet){
         FillHist("LepInLJetFr", 0., weight, 0., 7., 7);
         if(MuInLJet || EleInLJet)       FillHist("LepInLJetFr", 1., weight, 0., 7., 7);
         if(EleInLJet)                   FillHist("LepInLJetFr", 2., weight, 0., 7., 7);
         if(MuInLJet)                    FillHist("LepInLJetFr", 3., weight, 0., 7., 7);
         if(POGEleInLJet || POGMuInLJet) FillHist("LepInLJetFr", 4., weight, 0., 7., 7);
         if(POGEleInLJet)                FillHist("LepInLJetFr", 5., weight, 0., 7., 7);
         if(POGMuInLJet)                 FillHist("LepInLJetFr", 6., weight, 0., 7., 7);
  
         FillHist("NEleInLJet", NEleInThisLJet, weight, 0., 10., 10);
         FillHist("NMuInLJet", NMuInThisLJet, weight, 0., 10., 10);
         FillHist("NPOGEleInLJet", NPOGEleInThisLJet, weight, 0., 10., 10);
         FillHist("NPOGMuInLJet", NPOGMuInThisLJet, weight, 0., 10., 10);
       }


     }//End of Jet Loop
   
   }


 
/////////////////////////////////////////////////////////////////////////////////// 

return;
}// End of execute event loop
  


void May2017_MCFakeStudy::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void May2017_MCFakeStudy::BeginCycle() throw( LQError ){
  
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

May2017_MCFakeStudy::~May2017_MCFakeStudy() {
  
  Message("In May2017_MCFakeStudy Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


int May2017_MCFakeStudy::StepPassed(std::vector<snu::KMuon> MuColl, std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, TString Option){

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

   return PassedSteps;
}


int May2017_MCFakeStudy::NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   bool DrawHist=false;
   if     (Option.Contains("Prompt")) ReturnPrompt=true;
   else if(Option.Contains("Fake"))   ReturnFake=true;

   for(int i=0; i<MuColl.size(); i++){
     int MuType=GetLeptonType(MuColl.at(i), truthColl);

     if     (MuType==1) Nprompt++;
     else if(MuType==2) Nprompt++;
     else if(MuType==3) Nprompt++;
     else if(MuType<0 ) Nfake++;
     else if(MuType>=4) Nfake++;
   }

   if     (ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}

int May2017_MCFakeStudy::NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   bool DrawHist=false;
   if     (Option.Contains("Prompt")) ReturnPrompt=true;
   else if(Option.Contains("Fake"))   ReturnFake=true;

   for(int i=0; i<EleColl.size(); i++){
     int EleType=GetLeptonType(EleColl.at(i), truthColl);

     if     (EleType==1) Nprompt++;
     else if(EleType==2) Nprompt++;
     else if(EleType==3) Nprompt++;
     else if(EleType<0 ) Nfake++;
     else if(EleType>=4) Nfake++;
   }

   if     (ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}



bool May2017_MCFakeStudy::IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option){

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



void May2017_MCFakeStudy::FillCutFlow(TString cut, float weight){
  
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



void May2017_MCFakeStudy::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void May2017_MCFakeStudy::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);



  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this May2017_MCFakeStudyCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void May2017_MCFakeStudy::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
