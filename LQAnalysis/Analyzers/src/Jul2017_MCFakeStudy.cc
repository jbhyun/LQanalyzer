/***************************************************************************
 * @Project: Jul2017_MCFakeStudy 
 * @Package: LQanalyzer, ROOT-based analysis framework for Korea SNU / Main package author: John Almond
 *
 * @Code Author: Jihwan Bhyun 
 *
 ***************************************************************************/

/// Local includes
 #include "Jul2017_MCFakeStudy.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Jul2017_MCFakeStudy);

 Jul2017_MCFakeStudy::Jul2017_MCFakeStudy() : AnalyzerCore(), out_muons(0) {

   SetLogName("Jul2017_MCFakeStudy");
   Message("In Jul2017_MCFakeStudy constructor", INFO);
   InitialiseAnalysis();
 }


 void Jul2017_MCFakeStudy::InitialiseAnalysis() throw( LQError ) {
   
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

void Jul2017_MCFakeStudy::ExecuteEvents()throw( LQError ){

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


   bool Scan_AvgFR_MVAIso_2D=false, FineTune=false, EleOnlyClosure=false, CompositionCheck=false;
   bool SelBiasTest=false, TrigBiasTest=false, VetoImpactStudy=false, VJetReductionStudy=false, METMTWCutStudy=false;
   bool AppRegCompStudy=false, IPMissHitImpactStudy=false, HighdXYApproachStudy=false, IPOptStudy=false;
   bool FREmul=false, DetailClosure=false;
  
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("Scan_AvgFR_MVAIso_2D")) Scan_AvgFR_MVAIso_2D =true;
     else if(k_flags.at(i).Contains("FineTune"))             FineTune             =true;
     else if(k_flags.at(i).Contains("FREmul"))               FREmul               =true;
     else if(k_flags.at(i).Contains("DetailClosure"))        DetailClosure        =true;
     else if(k_flags.at(i).Contains("CompositionCheck"))     CompositionCheck     =true;
     else if(k_flags.at(i).Contains("EleOnlyClosure"))       EleOnlyClosure       =true;
     else if(k_flags.at(i).Contains("TrigBiasTest"))         TrigBiasTest         =true;
     else if(k_flags.at(i).Contains("SelBiasTest"))          SelBiasTest          =true;
     else if(k_flags.at(i).Contains("VetoImpactStudy"))      VetoImpactStudy      =true;
     else if(k_flags.at(i).Contains("METMTWCutStudy"))       METMTWCutStudy       =true;
     else if(k_flags.at(i).Contains("VJetReductionStudy"))   VJetReductionStudy   =true;
     else if(k_flags.at(i).Contains("AppRegCompStudy"))      AppRegCompStudy      =true;
     else if(k_flags.at(i).Contains("IPMissHitImpactStudy")) IPMissHitImpactStudy =true;
     else if(k_flags.at(i).Contains("HighdXYApproachStudy")) HighdXYApproachStudy =true;
     else if(k_flags.at(i).Contains("IPOptStudy"))           IPOptStudy           =true;
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
   if(AppRegCompStudy){ if( !((muonPreColl.size()>=2 && electronPreColl.size()>=1) || muonPreColl.size()>=3) ) return; }
   else if(DetailClosure || EleOnlyClosure){ if( !(muonPreColl.size()>=2 && electronPreColl.size()>=1) ) return; }
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

   std::vector<snu::KMuon> muonLooseColl;
     for(int i=0; i<muonPreColl.size(); i++){
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop4IPp5p1Chi100","Roch")) muonLooseColl.push_back(muonPreColl.at(i));
       //if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop6IPp5p1","Roch")) muonLooseColl.push_back(muonPreColl.at(i));
     }



   //Electron ID's to Test
   //
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_LOOSE);
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.5, 0.5);//2016 80X tuned WP
     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.1);    eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.1);//Not in ID, but additional safe WP
     eventbase->GetElectronSel()->SetdxySigMax(3.);
   std::vector<snu::KElectron> electronVetoColl; eventbase->GetElectronSel()->Selection(electronVetoColl);


     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronFakeLPreOptColl; eventbase->GetElectronSel()->Selection(electronFakeLPreOptColl);


//     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
//     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
//     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
//     eventbase->GetElectronSel()->SetBETrRegIncl(false);
//     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
//     eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
//     eventbase->GetElectronSel()->SetdxySigMax(3.);
//     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronFakeLColl; //eventbase->GetElectronSel()->Selection(electronFakeLColl);
     for(int i=0; i<electronFakeLPreOptColl.size(); i++){
       if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMTIsop06IPp025p05Sig4FakeLIsop4")) electronFakeLColl.push_back(electronFakeLPreOptColl.at(i));
     }



     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronMVAMColl; eventbase->GetElectronSel()->Selection(electronMVAMColl);

//     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP80);
//     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
//     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
//     eventbase->GetElectronSel()->SetBETrRegIncl(false);
//     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
//     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
//     eventbase->GetElectronSel()->SetdxySigMax(4.);
//     eventbase->GetElectronSel()->SetApplyConvVeto(true);
//   std::vector<snu::KElectron> electronMVATColl; eventbase->GetElectronSel()->Selection(electronMVATColl);
//   std::vector<snu::KElectron> electronMVAMColl=electronMVATColl;


   std::vector<snu::KElectron> electronNull;

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
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
         //fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", electronLooseColl.size());
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


   if(Scan_AvgFR_MVAIso_2D){
     std::vector<snu::KElectron> FakeColl     = SkimLepColl(electronFakeLPreOptColl, truthColl, "HFake");
     std::vector<snu::KJet>      JetCleanColl = SkimJetColl(jetColl, truthColl, "NoPrNoTau");
     
     const int NIsoCuts=11, NMVACuts=21;
     float IsoCuts[NIsoCuts]={0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};
     float MVACuts[NMVACuts]={-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};

     const int NPtEdges=8;
     float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 150., 200.};


     //2D Scanning MVA-Iso
     //1. Check ID var range unbiased from trigger requirement---------------------------------//
     for(int i=0; i<FakeColl.size(); i++){
       FillHist("TrigCutCheck_Iso_SumW", FakeColl.at(i).PFRelIso(0.3), weight, 0., 1., 100);
       FillHist("TrigCutCheck_MVA_SumW", FakeColl.at(i).MVA(), weight, -1., 1., 20);
       if(FakeColl.at(i).TriggerMatched("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v")){
         FillHist("TrigCutCheck_Iso_TrigSumW", FakeColl.at(i).PFRelIso(0.3), weight, 0., 1., 100);
         FillHist("TrigCutCheck_MVA_TrigSumW", FakeColl.at(i).MVA(), weight, -1., 1., 20);
       }
     }
     //----------------------------------------------------------------------------------------//


     //2. 2D Scan based on avg FR--------------------------------------------------------------//
     Draw2DFakePlot(FakeColl, JetCleanColl, NMVACuts, MVACuts, NIsoCuts, IsoCuts, NPtEdges, PtEdges, "ConePt");
     //----------------------------------------------------------------------------------------//
   
   }//End of Ele Fake WP Scan in MVA-Iso 2D Plane
   if(FineTune){
     std::vector<snu::KElectron> FakeColl     = SkimLepColl(electronFakeLPreOptColl, truthColl, "HFake");
     std::vector<snu::KElectron> FakeB1Coll   = SkimLepColl(FakeColl, "B1");
     std::vector<snu::KElectron> FakeB2Coll   = SkimLepColl(FakeColl, "B2");
     std::vector<snu::KElectron> FakeEColl    = SkimLepColl(FakeColl, "E");
     std::vector<snu::KJet>      JetCleanColl = SkimJetColl(jetColl, truthColl, "NoPrNoTau");


     const int NPtEdges=7;
     float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 200.};

     bool FullScanning=false, ScanSpecificWPs=true;
     if(FullScanning){
       //Full Range Scanning
       //const int NIsoCuts=6, NMVACuts=21;
       //float IsoCuts[NIsoCuts]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
       //float MVACuts[NMVACuts]={-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
  
       //For Scanning with Fine Steps
       const int NIsoCuts=1, NMVACuts=100;
       float IsoCuts[NIsoCuts]={0.4};
       float MVACuts[NMVACuts]={-0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93, -0.92, -0.91, -0.9,
                                -0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83, -0.82, -0.81, -0.8,
                                -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73, -0.72, -0.71, -0.7,
                                -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63, -0.62, -0.61, -0.6,
                                -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53, -0.52, -0.51, -0.5,
                                -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43, -0.42, -0.41, -0.4,
                                -0.39, -0.38, -0.37, -0.36, -0.35, -0.34, -0.33, -0.32, -0.31, -0.3,
                                -0.29, -0.28, -0.27, -0.26, -0.25, -0.24, -0.23, -0.22, -0.21, -0.2,
                                -0.19, -0.18, -0.17, -0.16, -0.15, -0.14, -0.13, -0.12, -0.11, -0.1,
                                -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.0};
                                

  
       for(int it_iso=0; it_iso<NIsoCuts; it_iso++){
         for(int it_mva=0; it_mva<NMVACuts; it_mva++){

           //TIsop06
           ScanFakeRate(FakeColl  , JetCleanColl, MVACuts[it_mva], IsoCuts[it_iso], NPtEdges, PtEdges,
                        "LNoMVANoIsoIPp025p1sig4", "POGWP90Isop06IPp025p1sig4", "_POGWP90Isop06IPp025p1sig4", "All");
           ScanFakeRate(FakeB1Coll, JetCleanColl, MVACuts[it_mva], IsoCuts[it_iso], NPtEdges, PtEdges,
                        "LNoMVANoIsoIPp025p1sig4", "POGWP90Isop06IPp025p1sig4", "_POGWP90Isop06IPp025p1sig4", "B1");
           ScanFakeRate(FakeB2Coll, JetCleanColl, MVACuts[it_mva], IsoCuts[it_iso], NPtEdges, PtEdges,
                        "LNoMVANoIsoIPp025p1sig4", "POGWP90Isop06IPp025p1sig4", "_POGWP90Isop06IPp025p1sig4", "B2");
           ScanFakeRate(FakeEColl , JetCleanColl, MVACuts[it_mva], IsoCuts[it_iso], NPtEdges, PtEdges,
                        "LNoMVANoIsoIPp025p1sig4", "POGWP90Isop06IPp025p1sig4", "_POGWP90Isop06IPp025p1sig4", "E");

         }
       }
     }
     if(ScanSpecificWPs){
       Draw1DFakePlot(FakeColl  , JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "All");
       Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "B1" );
       Draw1DFakePlot(FakeB2Coll, JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "B2" );
       Draw1DFakePlot(FakeEColl , JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "E"  );
     }
   }
   if(FREmul){
     const int NPtEdges=7; float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 200.};

//     EmulateFRMeasurement(electronPreColl, muonLooseColl, jetColl, bjetColl, met, met_x, met_y, weight, NPtEdges, PtEdges,
//                         "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "_WP90Isop06_L928576Isop4", "");


     EmulateFRMeasurement(electronPreColl, muonLooseColl, jetColl, bjetColl, met, met_x, met_y, weight, NPtEdges, PtEdges,
                         "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "_WP90Isop06IPp025p1sig4_L928878Isop4", "");


   }
   if(CompositionCheck){
     //Check Src Flav Composition in 1) Meas Reg(QCD) 2) Application Reg(DY,TT)(But except TTT evts)
     
     std::vector<snu::KElectron> electronFakeLColl;
       for(int i=0; i<(int) electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "HctoWAFakeLoose")) electronFakeLColl.push_back(electronFakeLPreOptColl.at(i));
       }


     for(int i=0; i<(int) electronFakeLColl.size(); i++){
       if(k_sample_name.Contains("DY") || k_sample_name.Contains("TT")){
         if(PassIDCriteria(electronFakeLColl.at(i), "POGWP90Isop06IPp025p1sig4")) continue;
         int LepType=GetLeptonType(electronFakeLColl.at(i), truthColl);
         if(LepType>0 || LepType<-4) continue;
         FillHist("FakeLepType", LepType, weight, -5., 5., 10);
       }
       int FakeSrcType=GetFakeLepSrcType(electronFakeLColl.at(i), jetColl);
       FillHist("FakeSrcType", FakeSrcType, weight, -5., 5., 10);
     }

   }
   if(EleOnlyClosure){

     std::vector<snu::KJet> JetCleanColl  = SkimJetColl(jetColl, truthColl, "NoPrNoTau");
     std::vector<snu::KJet> BJetCleanColl = SelBJets(JetCleanColl, "Medium");

     //Draw1DClosurePlot(muonHN2FakeTColl, electronFakeLPreOptColl, JetCleanColl, BJetCleanColl, truthColl, met, -0.7, -0.5, -0.6, 0.4, "");
     //Draw1DClosurePlot(muonHN2FakeTColl, electronFakeLPreOptColl, JetCleanColl, BJetCleanColl, truthColl, met, -0.9, -0.9, -0.8, 0.4, "");
     Draw1DClosurePlot(muonHN2FakeTColl, electronFakeLPreOptColl, JetCleanColl, BJetCleanColl, truthColl, met, 0., 0.1, -0.1, 0.4, "");

   }
   if(TrigBiasTest){
     std::vector<snu::KElectron> FakeColl     = SkimLepColl(electronFakeLPreOptColl, truthColl, "HFake");
     std::vector<snu::KElectron> FakeB1Coll   = SkimLepColl(FakeColl, "B1");
     std::vector<snu::KElectron> FakeB2Coll   = SkimLepColl(FakeColl, "B2");
     std::vector<snu::KElectron> FakeEColl    = SkimLepColl(FakeColl, "E");
     std::vector<snu::KJet>      JetCleanColl = SkimJetColl(jetColl, truthColl, "NoPrNoTau");

     const int NPtEdges=6;
     float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100.};

     Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "All");
     Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "B1");
     Draw1DFakePlot(FakeB2Coll, JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "B2");
     Draw1DFakePlot(FakeEColl,  JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "E");
     if(PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")){
       Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "AllTrig");
       Draw1DFakePlot(FakeB1Coll, JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "B1Trig");
       Draw1DFakePlot(FakeB2Coll, JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "B2Trig");
       Draw1DFakePlot(FakeEColl,  JetCleanColl, -0.92, -0.88, -0.78, 0.4, NPtEdges, PtEdges, "ETrig");
     }
   }
   if(SelBiasTest){
     std::vector<snu::KElectron> electronFakeL1Coll;
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMFakeLIso04Opt1")) electronFakeL1Coll.push_back(electronFakeLPreOptColl.at(i));
       }
     std::vector<snu::KJet> jetVeto1Coll= SkimJetColl(jetColl, electronFakeL1Coll, muonHN2FakeLColl, "EleMuVeto");

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

   }
   if(VetoImpactStudy){
     //To check whether we need additional veto point

     std::vector<snu::KElectron> electronFakeL1Coll;
       for(int i=0; i<electronPreColl.size(); i++){
         if(PassIDCriteria(electronPreColl.at(i), "POGMVAMFakeLIso04Opt1")) electronFakeL1Coll.push_back(electronPreColl.at(i));
       }
     std::vector<snu::KElectron> electron25FakeL1Coll = SkimLepColl(electronFakeL1Coll, "PtB1B2E", 25.);
     std::vector<snu::KJet> jetVeto1Coll= SkimJetColl(jetColl, electronFakeL1Coll, muonHN2FakeLColl, "EleMuVeto");


     std::vector<snu::KElectron> electronFakeL2Coll;
       for(int i=0; i<electronPreColl.size(); i++){
         if(PassIDCriteria(electronPreColl.at(i), "POGMVAMFakeLIso04Opt2")) electronFakeL2Coll.push_back(electronPreColl.at(i));
       }
     std::vector<snu::KElectron> electron25FakeL2Coll=SkimLepColl(electronFakeL2Coll, "PtB1B2E", 25.);
     std::vector<snu::KJet> jetVeto2Coll= SkimJetColl(jetColl, electronFakeL2Coll, muonHN2FakeLColl, "EleMuVeto");

     std::vector<snu::KElectron> electronFakeL2NoTrigColl;
       for(int i=0; i<electronPreColl.size(); i++){
         if(PassIDCriteria(electronPreColl.at(i), "POGMVAMFakeLIso04Opt2NoTrig")) electronFakeL2NoTrigColl.push_back(electronPreColl.at(i));
       }
     std::vector<snu::KElectron> electron25FakeL2NoTrigColl=SkimLepColl(electronFakeL2NoTrigColl, "PtB1B2E", 25.);


     //1. N(2l)/N(1l) in fake sample
     if(electronFakeL1Coll.size()>0 && jetVeto1Coll.size()>0){
       if(electronFakeL1Coll.at(0).Pt()>25 && jetVeto1Coll.at(0).Pt()>40){
         FillHist("NEleLOpt1_25_10", electronFakeL1Coll.size(), weight, 0., 10., 10);
         FillHist("NEleLOpt1_25_25", electron25FakeL1Coll.size(), weight, 0., 10., 10);
       }
     }

     if(electronFakeL2Coll.size()>0 && jetVeto2Coll.size()>0){
       if(electronFakeL2Coll.at(0).Pt()>25 && jetVeto2Coll.at(0).Pt()>40){
         FillHist("NEleLOpt2_25_10", electronFakeL2Coll.size(), weight, 0., 10., 10);
         FillHist("NEleLOpt2_25_25", electron25FakeL2Coll.size(), weight, 0., 10., 10);
       }
     }

     //2. Efficiency Check
     std::vector<snu::KElectron> electron25PreColl          = SkimLepColl(electronPreColl, "PtB1B2E", 25.);
     std::vector<snu::KElectron> electronPreFakeColl        = SkimLepColl(electronPreColl, truthColl, "HFake");
     std::vector<snu::KElectron> electronPrePromptColl      = SkimLepColl(electronPreColl, truthColl, "Prompt");
     std::vector<snu::KElectron> electron25PreFakeColl      = SkimLepColl(electron25PreColl, truthColl, "HFake");
     std::vector<snu::KElectron> electron25PrePromptColl    = SkimLepColl(electron25PreColl, truthColl, "Prompt");


     std::vector<snu::KElectron> electronFakeL1FakeColl     = SkimLepColl(electronFakeL1Coll, truthColl, "HFake");
     std::vector<snu::KElectron> electronFakeL1PromptColl   = SkimLepColl(electronFakeL1Coll, truthColl, "Prompt");
     std::vector<snu::KElectron> electron25FakeL1FakeColl   = SkimLepColl(electron25FakeL1Coll, truthColl, "HFake");
     std::vector<snu::KElectron> electron25FakeL1PromptColl = SkimLepColl(electron25FakeL1Coll, truthColl, "Prompt");

     std::vector<snu::KElectron> electronFakeL2FakeColl     = SkimLepColl(electronFakeL2Coll, truthColl, "HFake");
     std::vector<snu::KElectron> electronFakeL2PromptColl   = SkimLepColl(electronFakeL2Coll, truthColl, "Prompt");
     std::vector<snu::KElectron> electron25FakeL2FakeColl   = SkimLepColl(electron25FakeL2Coll, truthColl, "HFake");
     std::vector<snu::KElectron> electron25FakeL2PromptColl = SkimLepColl(electron25FakeL2Coll, truthColl, "Prompt");

     std::vector<snu::KElectron> electronFakeL2NoTrigFakeColl     = SkimLepColl(electronFakeL2NoTrigColl, truthColl, "HFake");
     std::vector<snu::KElectron> electronFakeL2NoTrigPromptColl   = SkimLepColl(electronFakeL2NoTrigColl, truthColl, "Prompt");
     std::vector<snu::KElectron> electron25FakeL2NoTrigFakeColl   = SkimLepColl(electron25FakeL2NoTrigColl, truthColl, "HFake");
     std::vector<snu::KElectron> electron25FakeL2NoTrigPromptColl = SkimLepColl(electron25FakeL2NoTrigColl, truthColl, "Prompt");


     FillHist("IDeff_NElepf_10", 0., electronPrePromptColl.size()*weight, 0., 2., 2);
     FillHist("IDeff_NElepf_10", 1., electronPreFakeColl.size()*weight, 0., 2., 2);

     FillHist("IDeff_NEleIDp_10", 0., electronFakeL1PromptColl.size()*weight, 0., 5., 5);
     FillHist("IDeff_NEleIDp_10", 1., electronFakeL2PromptColl.size()*weight, 0., 5., 5);
     FillHist("IDeff_NEleIDp_10", 2., electronFakeL2NoTrigPromptColl.size()*weight, 0., 5., 5);
     FillHist("IDeff_NEleIDf_10", 0., electronFakeL1FakeColl.size()*weight, 0., 5., 5);
     FillHist("IDeff_NEleIDf_10", 1., electronFakeL2FakeColl.size()*weight, 0., 5., 5);
     FillHist("IDeff_NEleIDf_10", 2., electronFakeL2NoTrigFakeColl.size()*weight, 0., 5., 5);

     FillHist("IDeff_NElepf_25", 0., electron25PrePromptColl.size()*weight, 0., 2., 2);
     FillHist("IDeff_NElepf_25", 1., electron25PreFakeColl.size()*weight, 0., 2., 2);

     FillHist("IDeff_NEleIDp_25", 0., electron25FakeL1PromptColl.size()*weight, 0., 5., 5);
     FillHist("IDeff_NEleIDp_25", 1., electron25FakeL2PromptColl.size()*weight, 0., 5., 5);
     FillHist("IDeff_NEleIDp_25", 2., electron25FakeL2NoTrigPromptColl.size()*weight, 0., 5., 5);
     FillHist("IDeff_NEleIDf_25", 0., electron25FakeL1FakeColl.size()*weight, 0., 5., 5);
     FillHist("IDeff_NEleIDf_25", 1., electron25FakeL2FakeColl.size()*weight, 0., 5., 5);
     FillHist("IDeff_NEleIDf_25", 2., electron25FakeL2NoTrigFakeColl.size()*weight, 0., 5., 5);

   }
   if(VJetReductionStudy){
     std::vector<snu::KElectron> electronFakeL2Coll;
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMFakeLIso04Opt2")) electronFakeL2Coll.push_back(electronFakeLPreOptColl.at(i));
       }
     std::vector<snu::KJet> jetVeto2Coll= SkimJetColl(jetColl, electronFakeL2Coll, muonHN2FakeLColl, "EleMuVeto");

     bool PassSel=true;
     if( !(electronMVAMColl.size()==1 && electronFakeL2Coll.size()==1) ) PassSel=false;//To remove DY Cont. && Ambiguity in meas.
     if(PassSel && jetVeto2Coll.size()==0     ) PassSel=false;//Should fire trigger
     if(PassSel && jetVeto2Coll.at(0).Pt()<40 ) PassSel=false;//Should fire trigger
     if(PassSel && jetVeto2Coll.at(0).DeltaR(electronFakeL2Coll.at(0))<1.0) PassSel=false;//the lepton is independent of this jet
     if(PassSel){

       float MTW=sqrt(2)*sqrt(met*electronFakeL2Coll.at(0).Pt()-met_x*electronFakeL2Coll.at(0).Px()-met_y*electronFakeL2Coll.at(0).Py());
       FillHist("NEvt_MET_1D", met, weight, 0., 200., 200);
       FillHist("NEvt_MTW_1D", MTW, weight, 0., 200., 200);
       FillHist("NEvt_METMTW_2D", met, MTW, weight, 0., 200., 40, 0., 200., 40); 
       for(int i=1; i<=40; i++){
         if(met<i*5.0){
           for(int j=1; j<=40; j++){
             if(MTW<j*5.0) FillHist("NEvt_ltMETltMTW_2D", (i-1)*5.0, (j-1)*5.0, weight, 0., 200., 40, 0., 200., 40); 
             else continue;
           }
         }
         else continue;
       }

       if( !(met<25 && MTW<35) ) PassSel=false;
     }
     if(PassSel){
       FillHist("RelPt_jol", jetVeto2Coll.at(0).Pt()/electronMVAMColl.at(0).Pt(), weight, 0., 10., 100);
       FillHist("CHEMFr_j", jetVeto2Coll.at(0).ChargedEMEnergyFraction(), weight, 0., 1., 100); 
     }

     
   }
   if(METMTWCutStudy){

     std::vector<snu::KElectron> electronFakeL2Coll;
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMFakeLIso04Opt2")) electronFakeL2Coll.push_back(electronFakeLPreOptColl.at(i));
       }
     std::vector<snu::KJet> jetVeto2Coll= SkimJetColl(jetColl, electronFakeL2Coll, muonHN2FakeLColl, "EleMuVeto");

     //const int NMETEdges=5, NMTWEdges=5;
     //float METEdges[NMETEdges]={0., 50., 100., 150., 200.};
     //float MTWEdges[NMTWEdges]={0., 50., 100., 150., 200.};

     bool PassSel=true;
     if(electronFakeL2Coll.size()!=1          ) PassSel=false;//To remove DY Cont. && Ambiguity in meas.
     if(PassSel && jetVeto2Coll.size()==0     ) PassSel=false;//Should fire trigger
     if(PassSel && jetVeto2Coll.at(0).Pt()<40 ) PassSel=false;//Should fire trigger
     if(PassSel && jetVeto2Coll.at(0).DeltaR(electronFakeL2Coll.at(0))<1.0) PassSel=false;//the lepton is independent of this jet
     if(PassSel){
       float MTW=sqrt(2)*sqrt(met*electronFakeL2Coll.at(0).Pt()-met_x*electronFakeL2Coll.at(0).Px()-met_y*electronFakeL2Coll.at(0).Py());
       FillHist("NEvt_MET_1D", met, weight, 0., 200., 200);
       FillHist("NEvt_MTW_1D", MTW, weight, 0., 200., 200);
       FillHist("NEvt_METMTW_2D", met, MTW, weight, 0., 200., 40, 0., 200., 40); 
       for(int i=1; i<=40; i++){
         if(met<i*5.0){
           for(int j=1; j<=40; j++){
             if(MTW<j*5.0) FillHist("NEvt_ltMETltMTW_2D", (i-1)*5.0, (j-1)*5.0, weight, 0., 200., 40, 0., 200., 40); 
             else continue;
           }
         }
         else continue;
       }
       //FillHist("NEvt_METMTW_2D", met, MTW, weight, METEdges, NMETEdges-1, MTWEdges, NMTWEdges-1); 
     }

   }
   if(AppRegCompStudy){

     //Trig Part---------------------------------------------------------------------------//
     int Pass_Trigger1=0, Pass_Trigger2=0;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
     if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;
  
     float trigger_period_weight=1.;
     bool Pass_Trigger=false;
     if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_Trigger=true;
     trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;
     if(!Pass_Trigger) return;
     weight*=trigger_period_weight;
     //-------------------------------------------------------------------------------------//


     if(k_running_nonprompt){} //m_datadriven_bkg->GetFakeObj()->SetUseQCDFake(true);

     std::vector<snu::KElectron> electronFakeL2Coll;
       for(int i=0; i<electronPreColl.size(); i++){
         if(PassIDCriteria(electronPreColl.at(i), "POGMVAMFakeLIso04Opt2")) electronFakeL2Coll.push_back(electronPreColl.at(i));
       }


     if( !(electronFakeL2Coll.size()==1      && muonHN2FakeLColl.size()==2) ) return;
     if( !(electronFakeL2Coll.at(0).Pt()>25. && muonHN2FakeLColl.at(0).Pt()>15. && muonHN2FakeLColl.at(1).Pt()>10.) ) return;
     if( !(muonHN2FakeLColl.at(0).Charge()!=muonHN2FakeLColl.at(1).Charge()) ) return;

     float Mmumu=(muonHN2FakeLColl.at(0)+muonHN2FakeLColl.at(1)).M();
     if( Mmumu<12. ) return;

     bool OnZ=false;
     if( fabs(Mmumu-91.2)<10. ) OnZ=true;


     int NEleT=0, NEleLnT=0, NMuT=0, NMuLnT=0;
     int EleType=0, Mu1Type=0, Mu2Type=0;
     int EleTypeNew=0, Mu1TypeNew=0, Mu2TypeNew=0;

     float fakeweight=-1.;
     if(PassIDCriteria(electronFakeL2Coll.at(0), "POGMVAMIP")) NEleT++;
     else {NEleLnT++; 
           float FR=FakeRateMC(electronFakeL2Coll.at(0), "QCD_EGM");
           fakeweight*=-FR/(1.-FR);}

     if(muonHN2FakeLColl.at(0).RelIso04()<0.1) NMuT++;
     else {NMuLnT++;
           float FR=0.;
           //float FR=m_datadriven_bkg->GetFakeObj()->getTrilepFakeRate_muon(false, muonHN2FakeLColl.at(0).Pt(), muonHN2FakeLColl.at(0).Eta());
           fakeweight*=-FR/(1.-FR);}


     if(muonHN2FakeLColl.at(1).RelIso04()<0.1) NMuT++;
     else{ NMuLnT++;
           float FR=0.;
           //float FR=m_datadriven_bkg->GetFakeObj()->getTrilepFakeRate_muon(false, muonHN2FakeLColl.at(1).Pt(), muonHN2FakeLColl.at(1).Eta());
           fakeweight*=-FR/(1.-FR);}

     if( (NEleT+NMuT)==3 ) fakeweight==0;


     EleType=GetLeptonType(electronFakeL2Coll.at(0), truthColl);
     Mu1Type=GetLeptonType(muonHN2FakeLColl.at(0)  , truthColl);
     Mu2Type=GetLeptonType(muonHN2FakeLColl.at(1)  , truthColl);

     
     if     (EleType==3)              EleTypeNew= 2;
     else if(EleType>3)               EleTypeNew= 3;
     else if(EleType<0 && EleType>-5) EleTypeNew=-1;
     else if(EleType<-4)              EleTypeNew=-2;

     if     (Mu1Type==3)              Mu1TypeNew= 2;   
     else if(Mu1Type>3)               Mu1TypeNew= 3;   
     else if(Mu1Type<0 && Mu1Type>-5) Mu1TypeNew=-1;   
     else if(Mu1Type<-4)              Mu1TypeNew=-2;   

     if     (Mu2Type==3)              Mu2TypeNew= 2;
     else if(Mu2Type>3)               Mu2TypeNew= 3;
     else if(Mu2Type<0 && Mu2Type>-5) Mu2TypeNew=-1;
     else if(Mu2Type<-4)              Mu2TypeNew=-2;     
     
     if( !(NEleT==1 && NMuT==2) ){

       float EleFR=FakeRateMC(electronFakeL2Coll.at(0), "QCD_EGM");

       FillHist("NTightLoose", NEleT+NMuT, weight, 0., 4., 4);
       FillHist("NTightLoose_Ele", NEleT, weight, 0., 2., 2);
       FillHist("NTightLoose_Mu", NMuT, weight, 0., 3., 3);

       if(NEleLnT==1){
         FillHist("EleLnTType", EleType, weight, -10., 10., 20);
         FillHist("EleLnTType_Summary", EleTypeNew, weight, -5., 5., 10);
         FillHist("Closure_3lOS_Ele_exp", 0., weight*EleFR, 0., 1., 1);
         if(EleType<0 && EleType>-5) FillHist("Closure_3lOS_HadEle_exp", 0., weight*EleFR, 0., 1., 1);
       } 
       else{
         FillHist("Mu1LnTType", Mu1Type, weight, -10., 10., 20);
         FillHist("Mu1LnTType_Summary", Mu1TypeNew, weight, -10., 10., 20);
         FillHist("Mu2LnTType", Mu2Type, weight, -10., 10., 20);
         FillHist("Mu2LnTType_Summary", Mu2TypeNew, weight, -10., 10., 20);
       }
       FillHist("Closure_3lOS_exp", 0., weight*fakeweight, 0., 1., 1);


       if(OnZ){
         FillHist("NTightLoose_OnZ", NEleT+NMuT, weight, 0., 4., 4);
         FillHist("NTightLoose_Ele_OnZ", NEleT, weight, 0., 2., 2);
         FillHist("NTightLoose_Mu_OnZ", NMuT, weight, 0., 3., 3);
  
         if(NEleLnT==1){
           FillHist("EleLnTType_OnZ", EleType, weight, -10., 10., 20);
           FillHist("EleLnTType_Summary_OnZ", EleTypeNew, weight, -5., 5., 10);
           FillHist("Closure_3lOSOnZ_Ele_exp", 0., weight*EleFR, 0., 1., 1);
         } 
         else{
           FillHist("Mu1LnTType_OnZ", Mu1Type, weight, -10., 10., 20);
           FillHist("Mu1LnTType_Summary_OnZ", Mu1TypeNew, weight, -10., 10., 20);
           FillHist("Mu2LnTType_OnZ", Mu2Type, weight, -10., 10., 20);
           FillHist("Mu2LnTType_Summary_OnZ", Mu2TypeNew, weight, -10., 10., 20);
         }
         FillHist("Closure_3lOSOnZ_exp", 0., weight*fakeweight, 0., 1., 1);

       }
     }//At least 1 loose
     else{

       FillHist("Closure_3lOS_Ele_obs", 0., weight, 0., 1., 1);
       FillHist("Closure_3lOS_obs", 0., weight, 0., 1., 1);
       if(OnZ){ FillHist("Closure_3lOSOnZ_Ele_obs", 0., weight, 0., 1., 1);
                FillHist("Closure_3lOSOnZ_obs", 0., weight, 0., 1., 1);}
     }//AllTight Ends

   }//End of AppRegCompStudy
   if(IPMissHitImpactStudy){
     std::vector<snu::KElectron> electronFakeL2Coll;
       for(int i=0; i<electronPreColl.size(); i++){
         if(PassIDCriteria(electronPreColl.at(i), "POGMVAMFakeLIso04Opt2")) electronFakeL2Coll.push_back(electronPreColl.at(i));
       }
     std::vector<snu::KElectron> electron25FakeL2Coll=SkimLepColl(electronFakeL2Coll, "PtB1B2E", 25.);
     std::vector<snu::KJet> jetVeto2Coll= SkimJetColl(jetColl, electronFakeL2Coll, muonHN2FakeLColl, "EleMuVeto");

     for(int i=0; i<electron25FakeL2Coll.size(); i++){
       int EleType=GetLeptonType(electron25FakeL2Coll.at(i), truthColl);
       FillHist("EleType", EleType, weight, -10., 10., 20);
       float d0=fabs(electron25FakeL2Coll.at(i).dxy()), dz=fabs(electron25FakeL2Coll.at(i).dz()), d0Sig=fabs(electron25FakeL2Coll.at(i).dxySig());
       int NMissHit=electron25FakeL2Coll.at(i).MissingHits();
       
       
       //d0
       if(EleType==1)                      FillHist("d0_Type1", d0, weight, 0., 0.2, 200); 
       else if(EleType==3)                 FillHist("d0_Type3", d0, weight, 0., 0.2, 200); 
       else if(EleType>3 )                 FillHist("d0_Type45_IntConv", d0, weight, 0., 0.2, 200); 
       else if(EleType==-1)                FillHist("d0_Typem1", d0, weight, 0., 0.2, 200);
       else if(EleType<=-2 && EleType>=-4) FillHist("d0_Typem234", d0, weight, 0., 0.2, 200);
       else if(EleType<-4)                 FillHist("d0_Typem56_ExtConv", d0, weight, 0., 0.2, 200);

       //dz
       if     (EleType==1)                 FillHist("dz_Type1", dz, weight, 0., 0.2, 200);
       else if(EleType==3)                 FillHist("dz_Type3", dz, weight, 0., 0.2, 200);
       else if(EleType>3)                  FillHist("dz_Type45_IntConv", dz, weight, 0., 0.2, 200);
       else if(EleType==-1)                FillHist("dz_Typem1", dz, weight, 0., 0.2, 200);
       else if(EleType<=-2 && EleType>=-4) FillHist("dz_Typem234", dz, weight, 0., 0.2, 200);
       else if(EleType<-4)                 FillHist("dz_Typem56_ExtConv", dz, weight, 0., 0.2, 200);

       //d0Sig
       if     (EleType==1)                 FillHist("d0Sig_Type1", d0Sig, weight, 0., 10, 100);
       else if(EleType==3)                 FillHist("d0Sig_Type3", d0Sig, weight, 0., 10, 100);
       else if(EleType>3)                  FillHist("d0Sig_Type45_IntConv", d0Sig, weight, 0., 10, 100);
       else if(EleType==-1)                FillHist("d0Sig_Typem1", d0Sig, weight, 0., 10, 100);
       else if(EleType<=-2 && EleType>=-4) FillHist("d0Sig_Typem234", d0Sig, weight, 0., 10, 100);
       else if(EleType<-4)                 FillHist("d0Sig_Typem56_ExtConv", d0Sig, weight, 0., 10, 100);

       //MissHits
       if     (EleType==1)                 FillHist("NMissHit_Type1", NMissHit, weight, 0., 10., 10);
       else if(EleType==3)                 FillHist("NMissHit_Type3", NMissHit, weight, 0., 10., 10);
       else if(EleType>3)                  FillHist("NMissHit_Type45_IntConv", NMissHit, weight, 0., 10., 10);
       else if(EleType==-1)                FillHist("NMissHit_Typem1", NMissHit, weight, 0., 10., 10);
       else if(EleType<=-2 && EleType>=-4) FillHist("NMissHit_Typem234", NMissHit, weight, 0., 10., 10);
       else if(EleType<-4)                 FillHist("NMissHit_Typem56_ExtConv", NMissHit, weight, 0., 10., 10);

       if(PassIDCriteria(electron25FakeL2Coll.at(i), "POGMVAMIP")){
         FillHist("EleType_ID", EleType, weight, -10., 10., 20);

         //d0
         if(EleType==1)                      FillHist("d0_Type1_ID", d0, weight, 0., 0.2, 200); 
         else if(EleType==3)                 FillHist("d0_Type3_ID", d0, weight, 0., 0.2, 200); 
         else if(EleType>3 )                 FillHist("d0_Type45_IntConv_ID", d0, weight, 0., 0.2, 200); 
         else if(EleType==-1)                FillHist("d0_Typem1_ID", d0, weight, 0., 0.2, 200);
         else if(EleType<=-2 && EleType>=-4) FillHist("d0_Typem234_ID", d0, weight, 0., 0.2, 200);
         else if(EleType<-4)                 FillHist("d0_Typem56_ExtConv_ID", d0, weight, 0., 0.2, 200);
  
         //dz
         if     (EleType==1)                 FillHist("dz_Type1_ID", dz, weight, 0., 0.2, 200);
         else if(EleType==3)                 FillHist("dz_Type3_ID", dz, weight, 0., 0.2, 200);
         else if(EleType>3)                  FillHist("dz_Type45_IntConv_ID", dz, weight, 0., 0.2, 200);
         else if(EleType==-1)                FillHist("dz_Typem1_ID", dz, weight, 0., 0.2, 200);
         else if(EleType<=-2 && EleType>=-4) FillHist("dz_Typem234_ID", dz, weight, 0., 0.2, 200);
         else if(EleType<-4)                 FillHist("dz_Typem56_ExtConv_ID", dz, weight, 0., 0.2, 200);
  
         //d0Sig
         if     (EleType==1)                 FillHist("d0Sig_Type1_ID", d0Sig, weight, 0., 10, 100);
         else if(EleType==3)                 FillHist("d0Sig_Type3_ID", d0Sig, weight, 0., 10, 100);
         else if(EleType>3)                  FillHist("d0Sig_Type45_IntConv_ID", d0Sig, weight, 0., 10, 100);
         else if(EleType==-1)                FillHist("d0Sig_Typem1_ID", d0Sig, weight, 0., 10, 100);
         else if(EleType<=-2 && EleType>=-4) FillHist("d0Sig_Typem234_ID", d0Sig, weight, 0., 10, 100);
         else if(EleType<-4)                 FillHist("d0Sig_Typem56_ExtConv_ID", d0Sig, weight, 0., 10, 100);
  
         //MissHits
         if     (EleType==1)                 FillHist("NMissHit_Type1_ID", NMissHit, weight, 0., 10., 10);
         else if(EleType==3)                 FillHist("NMissHit_Type3_ID", NMissHit, weight, 0., 10., 10);
         else if(EleType>3)                  FillHist("NMissHit_Type45_IntConv_ID", NMissHit, weight, 0., 10., 10);
         else if(EleType==-1)                FillHist("NMissHit_Typem1_ID", NMissHit, weight, 0., 10., 10);
         else if(EleType<=-2 && EleType>=-4) FillHist("NMissHit_Typem234_ID", NMissHit, weight, 0., 10., 10);
         else if(EleType<-4)                 FillHist("NMissHit_Typem56_ExtConv_ID", NMissHit, weight, 0., 10., 10);
       }
     }
   }
   if(HighdXYApproachStudy){
     //Questions
     //1. Does high d0 approach have source bias?
     // 1-1) What kinds of fakes have high dxy? And what cut of high dxy is proper for those fakes?
     // 1-2) Does high PT b-jet fakes have high IP?
     //2. Is this approach worth studying?
     // 2-1) In high PT l+jet in high dxy, what fraction is prompt?(how much prompt reduction can be expected from high dxy method?)
     //3. Does MC based SF have bias from wrong IP modelling?
     // 3-1) does near b fakes(mostly b,c had daughter) have similar shape as MC b-fakes?

       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> electronFakeLNod0PreOptColl; eventbase->GetElectronSel()->Selection(electronFakeLNod0PreOptColl);

     std::vector<snu::KElectron> electronFakeL2Coll;
       for(int i=0; i<electronFakeLNod0PreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLNod0PreOptColl.at(i), "POGMVAMFakeLIso04Opt2Nod0")) electronFakeL2Coll.push_back(electronFakeLNod0PreOptColl.at(i));
       }
     std::vector<snu::KJet>      JetCleanColl = SkimJetColl(jetColl, truthColl, "NoPrNoTau");


     for(int i=0; i<electronFakeL2Coll.size(); i++){
       int EleType=GetLeptonType(electronFakeL2Coll.at(i), truthColl);
       int FakeSrcType = GetFakeLepSrcType(electronFakeL2Coll.at(i), JetCleanColl);

       float d0=fabs(electronFakeL2Coll.at(i).dxy()), d0Sig=fabs(electronFakeL2Coll.at(i).dxySig());

       FillHist("EleSumW_All_EleType", EleType, weight, -10., 10., 20);
       FillHist("EleSumWd0_All_EleType", EleType, weight*d0, -10., 10., 20);
       FillHist("EleSumWd0Sig_All_EleType", EleType, weight*d0Sig, -10., 10., 20);

       if(EleType<0){
         FillHist("d0_All", d0, weight, 0., 1., 1000);
         FillHist("d0Sig_All", d0Sig, weight, 0., 50., 50);
       }
       if(EleType<0 && FakeSrcType==3){
         FillHist("EleSumW_BjMatch_EleType", EleType, weight, -10., 10., 20);
         FillHist("EleSumWd0_BjMatch_EleType", EleType, weight*d0, -10., 10., 20);
         FillHist("EleSumWd0Sig_BjMatch_EleType", EleType, weight*d0Sig, -10., 10., 20);

         FillHist("d0_BjMatch", d0, weight, 0., 1., 1000);
         FillHist("d0Sig_BjMatch", d0Sig, weight, 0., 50., 50);
       }
     }

   }
   if(IPOptStudy){

       eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
     std::vector<snu::KElectron> electronPre25Coll; eventbase->GetElectronSel()->Selection(electronPre25Coll);
     std::vector<snu::KElectron> electronPOGMVAMColl; 
       for(int i=0; i<electronPre25Coll.size(); i++){
         if(PassIDCriteria(electronPre25Coll.at(i), "POGMVAM")) electronPOGMVAMColl.push_back(electronPre25Coll.at(i));
       }
     std::vector<snu::KElectron> electronPOGMVAMPrColl=SkimLepColl(electronPOGMVAMColl, truthColl, "Prompt");
     std::vector<snu::KElectron> electronPOGMVAMFkColl=SkimLepColl(electronPOGMVAMColl, truthColl, "HFake");


     //IsoWP Study
     int NIsoWP=40;
     for(int i=0; i<electronPOGMVAMPrColl.size(); i++){
       FillHist("Pr_Iso", electronPOGMVAMPrColl.at(i).PFRelIso(0.3), weight, 0., 0.3, 300);
       for(int it_iso=1; it_iso<=NIsoWP; it_iso++){
         FillHist("NElePrSumW_Iso", (it_iso-1.)*0.00501, weight, 0., 0.2, 40);
         if(electronPOGMVAMPrColl.at(i).PFRelIso(0.3)<it_iso*0.005){
           FillHist("NElePrIDSumW_Iso", (it_iso-1.)*0.00501, weight, 0., 0.2, 40);
         }
       }
     }
     for(int i=0; i<electronPOGMVAMFkColl.size(); i++){
       FillHist("Fk_Iso", electronPOGMVAMFkColl.at(i).PFRelIso(0.3), weight, 0., 0.3, 300);
       for(int it_iso=1; it_iso<=NIsoWP; it_iso++){
         FillHist("NEleFkSumW_Iso", (it_iso-1.)*0.00501, weight, 0., 0.2, 40);
         if(electronPOGMVAMFkColl.at(i).PFRelIso(0.3)<it_iso*0.005){
           FillHist("NEleFkIDSumW_Iso", (it_iso-1.)*0.00501, weight, 0., 0.2, 40);
         }
       }
     }

     //IPStudy
     int Nd0WP=10, NdzWP=20;
     for(int i=0; i<electronPOGMVAMPrColl.size(); i++){
       float d0=fabs(electronPOGMVAMPrColl.at(i).dxy()), dz=fabs(electronPOGMVAMPrColl.at(i).dz());
       FillHist("Pr_d0", d0, weight, 0., 0.1, 20);
       FillHist("Pr_dz", dz, weight, 0., 0.2, 40);
       for(int it_d0=1; it_d0<=Nd0WP; it_d0++){
         for(int it_dz=1; it_dz<=NdzWP; it_dz++){
           if(electronPOGMVAMPrColl.at(i).PFRelIso(0.3)<0.06){
             FillHist("NElePrSumW_d0dz_iso0p06", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             if(d0<it_d0*0.005 && dz< it_dz*0.005){
               FillHist("NElePrIDSumW_d0dz_iso0p06", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             }
           }
           if(electronPOGMVAMPrColl.at(i).PFRelIso(0.3)<0.08){
             FillHist("NElePrSumW_d0dz_iso0p08", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             if(d0<it_d0*0.005 && dz< it_dz*0.005){
               FillHist("NElePrIDSumW_d0dz_iso0p08", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             }
           }
           if(electronPOGMVAMPrColl.at(i).PFRelIso(0.3)<0.1){
             FillHist("NElePrSumW_d0dz_iso0p1", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             if(d0<it_d0*0.005 && dz< it_dz*0.005){
               FillHist("NElePrIDSumW_d0dz_iso0p1", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             }
           }
         }
       }
     }
     for(int i=0; i<electronPOGMVAMFkColl.size(); i++){
       float d0=fabs(electronPOGMVAMFkColl.at(i).dxy()), dz=fabs(electronPOGMVAMFkColl.at(i).dz());
       FillHist("Fk_d0", d0, weight, 0., 0.1, 20);
       FillHist("Fk_dz", dz, weight, 0., 0.2, 40);
       for(int it_d0=1; it_d0<=Nd0WP; it_d0++){
         for(int it_dz=1; it_dz<=NdzWP; it_dz++){
           if(electronPOGMVAMFkColl.at(i).PFRelIso(0.3)<0.06){
             FillHist("NEleFkSumW_d0dz_iso0p06", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             if(d0<it_d0*0.005 && dz< it_dz*0.005){
               FillHist("NEleFkIDSumW_d0dz_iso0p06", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             }
           }
           if(electronPOGMVAMFkColl.at(i).PFRelIso(0.3)<0.08){
             FillHist("NEleFkSumW_d0dz_iso0p08", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             if(d0<it_d0*0.005 && dz< it_dz*0.005){
               FillHist("NEleFkIDSumW_d0dz_iso0p08", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             }
           }
           if(electronPOGMVAMFkColl.at(i).PFRelIso(0.3)<0.1){
             FillHist("NEleFkSumW_d0dz_iso0p1", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             if(d0<it_d0*0.005 && dz< it_dz*0.005){
               FillHist("NEleFkIDSumW_d0dz_iso0p1", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             }
           }
         }
       }
     }
   }//End of IP optimisation
   if(DetailClosure){

     CheckMCClosure(muonPreColl, electronPreColl, jetColl, bjetColl, truthColl,
       "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4",
       "_WP90Isop06LMVA928576Isop4", "TrkIsoVVLConeSUSY");

//     CheckMCClosure(muonPreColl, electronPreColl, jetColl, bjetColl, truthColl,
//            "LMVA06Isop4IPp025p05sig4", "POGWP90Isop06IPp025p05sig4", "POGLIsop4IPp5p1Chi100", "POGTIsop20IPp01p05sig4Chi4", 
//            "_WP90Isop06LMVA928576Isop4", "TrkIsoVVLConeSUSY");
     
     CheckMCClosure(muonPreColl, electronPreColl, jetColl, bjetColl, truthColl,
       "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "POGLIsop4IPp5p1Chi100", "POGTIsop20IPp01p05sig4Chi4", 
       "_WP90Isop06IPp025p1sig4_L928878Isop4", "TrkIsoVVLConeSUSY");

   }
 
/////////////////////////////////////////////////////////////////////////////////// 

return;
}// End of execute event loop
  


void Jul2017_MCFakeStudy::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Jul2017_MCFakeStudy::BeginCycle() throw( LQError ){
  
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

Jul2017_MCFakeStudy::~Jul2017_MCFakeStudy() {
  
  Message("In Jul2017_MCFakeStudy Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


void Jul2017_MCFakeStudy::CheckMCClosure(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> ElePreColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, std::vector<snu::KTruth> TruthColl, TString EleLID, TString EleTID, TString MuLID, TString MuTID, TString Label, TString Option){


  TString FilterInfo="", ConeMethod="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";

  std::vector<snu::KMuon> MuLColl, MuTColl;
    for(int i=0; i<MuPreColl.size(); i++){
      if(PassIDCriteria(MuPreColl.at(i), MuLID)) MuLColl.push_back(MuPreColl.at(i));
      if(PassIDCriteria(MuPreColl.at(i), MuTID)) MuTColl.push_back(MuPreColl.at(i));
    }
  std::vector<snu::KElectron> EleLColl, EleTColl;
    for(int i=0; i<ElePreColl.size(); i++){
      if(PassIDCriteria(ElePreColl.at(i), EleLID)) EleLColl.push_back(ElePreColl.at(i));
      if(PassIDCriteria(ElePreColl.at(i), EleTID)) EleTColl.push_back(ElePreColl.at(i));
    }

  std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");

 
  //Checking 2PromptMu+1Fk Ele
  float fakeweight = -1.;
  int   NHFakeMu   = NPromptFake_Mu(MuLColl, TruthColl, "HFake");
  int   NHFakeEle  = NPromptFake_Ele(EleLColl, TruthColl, "HFake");
  int   NHFakeLep  = NHFakeMu+NHFakeEle;
  int   NLooseLep  = MuLColl.size()+EleLColl.size();
  int   NTightLep  = MuTColl.size()+EleTColl.size();
  if(NLooseLep!=3) return;
  //cout<<"NHFakeMu "<<NHFakeMu<<" NHFakeEle "<<NHFakeEle<<endl;
  if(NHFakeLep==0) return;//Include at least 1 fake in the sample
  

  for(int i=0; i<MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID)){
      float FR=FakeRateMC(MuLColl.at(i), "QCD_"+MuTID.ReplaceAll("Test_","")+"_"+MuLID.ReplaceAll("Test_","")+"_"+FilterInfo+ConeMethod);
      fakeweight*=-FR/(1-FR);
  //    cout<<"MuFR"<<FR<<endl;
    }
  }
  for(int i=0; i<EleLColl.size(); i++){
    if(!PassIDCriteria(EleLColl.at(i), EleTID)){
      float FR=FakeRateMC(EleLColl.at(i), "QCD_"+EleTID.ReplaceAll("Test_","")+"_"+EleLID.ReplaceAll("Test_","")+"_"+ConeMethod);
      fakeweight*=-FR/(1-FR);
   //   cout<<"EleFR"<<FR<<endl;
    }
  }
  //if( NLooseLep==NTightLep ) fakeweight=0;
  if( NLooseLep==3 && NTightLep==3 ) fakeweight=0.;

  int Pass_Trigger1=0, Pass_Trigger2=0;
  if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
  if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;
  
  float trigger_period_weight=1.;
  bool Pass_Trigger=false;
  if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_Trigger=true;
  trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;
  //weight*=trigger_period_weight;


  if(Pass_Trigger){
    if(MuLColl.size()==2 && EleLColl.size()==1 && NTightLep!=3 ){
      bool PassSel=true;
      float Mmumu=(MuLColl.at(0)+MuLColl.at(1)).M();
      if( !(EleLColl.at(0).Pt()>25 && MuLColl.at(0).Pt()>10 && MuLColl.at(1).Pt()>10) ) PassSel=false;
      if( PassSel && fabs(SumCharge(MuLColl))!=0) PassSel=false;
      if( PassSel && !(Mmumu>12) ) PassSel=false;
      if(PassSel){
        bool JetSelPass=JetColl.size()>=2, BJetSelPass=BJetColl.size()>=1;
        FillHist("CutFlow_exp"+Label, 0., weight*fakeweight, 0., 10., 10);
        if(JetSelPass) FillHist("CutFlow_exp"+Label, 1., weight*fakeweight, 0., 10., 10);
        if(JetSelPass && BJetSelPass) FillHist("CutFlow_exp"+Label, 2., weight*fakeweight, 0., 10., 10);

        FillHist("PTmu1_exp"+Label, MuLColl.at(0).Pt(), weight*fakeweight, 0., 200., 40);
        FillHist("PTmu2_exp"+Label, MuLColl.at(1).Pt(), weight*fakeweight, 0., 200., 40);
        FillHist("PTe_exp"+Label, EleLColl.at(0).Pt(), weight*fakeweight, 0., 200., 40);
        FillHist("Mmumu_exp"+Label, (MuLColl.at(0)+MuLColl.at(1)).M(), weight*fakeweight, 0., 200., 40);
        FillHist("Etamu1_exp"+Label, MuLColl.at(0).Eta(), weight*fakeweight, -5., 5., 20);
        FillHist("Etamu2_exp"+Label, MuLColl.at(1).Eta(), weight*fakeweight, -5., 5., 20);
        FillHist("Etae_exp"+Label, EleLColl.at(0).Eta(), weight*fakeweight, -5., 5., 20);
        FillHist("Nj_exp"+Label, JetColl.size(), weight*fakeweight, 0., 10., 10);
        FillHist("Nb_exp"+Label, BJetColl.size(), weight*fakeweight, 0., 10., 10);
      }
    }
    if(MuTColl.size()==2 && EleTColl.size()==1){
      bool PassSel=true;
      float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
      if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) PassSel=false;
      if( PassSel && fabs(SumCharge(MuTColl))!=0) PassSel=false;
      if( PassSel && !(Mmumu>12) ) PassSel=false;
      if(PassSel){
        bool JetSelPass=JetColl.size()>=2, BJetSelPass=BJetColl.size()>=1;
        FillHist("CutFlow_obs"+Label, 0., weight, 0., 10., 10);
        if(JetSelPass) FillHist("CutFlow_obs"+Label, 1., weight, 0., 10., 10);
        if(JetSelPass && BJetSelPass) FillHist("CutFlow_obs"+Label, 2., weight, 0., 10., 10);

        FillHist("PTmu1_obs"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
        FillHist("PTmu2_obs"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
        FillHist("PTe_obs"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
        FillHist("Mmumu_obs"+Label, (MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 200., 40);
        FillHist("Etamu1_obs"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
        FillHist("Etamu2_obs"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
        FillHist("Etae_obs"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
        FillHist("Nj_obs"+Label, JetColl.size(), weight, 0., 10., 10);
        FillHist("Nb_obs"+Label, BJetColl.size(), weight, 0., 10., 10);
      }
    }
  } 

}


void Jul2017_MCFakeStudy::Draw1DClosurePlot(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> ElePreColl, std::vector<snu::KJet> jetColl, std::vector<snu::KJet> bjetColl, std::vector<snu::KTruth> truthColl, float met, float MVAB1Cut, float MVAB2Cut, float MVAECut, float IsoCut, TString Option){

  std::vector<snu::KElectron> EleLColl;
  std::vector<snu::KElectron> EleTColl;
    for(int i=0; i<ElePreColl.size(); i++){
      if( ElePreColl.at(i).PFRelIso(0.3)>IsoCut ) continue;

      if     ( fabs(ElePreColl.at(i).Eta())<0.8   ){ if( ElePreColl.at(i).MVA()<MVAB1Cut ) continue; }
      else if( fabs(ElePreColl.at(i).Eta())<1.479 ){ if( ElePreColl.at(i).MVA()<MVAB2Cut ) continue; }
      else                                         { if( ElePreColl.at(i).MVA()<MVAECut  ) continue; }
     
      EleLColl.push_back(ElePreColl.at(i));
      //if(PassIDCriteria(ElePreColl.at(i), "POGMVAMIP")) EleTColl.push_back(ElePreColl.at(i)); 
      if(PassIDCriteria(ElePreColl.at(i), "POGMVATIP")) EleTColl.push_back(ElePreColl.at(i)); 
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
    //if(!PassIDCriteria(EleLColl.at(i), "POGMVAMIP")){
    if(!PassIDCriteria(EleLColl.at(i), "POGMVATIP")){
      //float FR=1.;
      //if(Str_IsoCut=="0p3") FR=FakeRateMC(EleLColl.at(i), "QCD_Iso03");
      //else if(Str_IsoCut=="0p4") FR=FakeRateMC(EleLColl.at(i), "QCD_Iso04");
      //else if(Str_IsoCut=="0p5") FR=FakeRateMC(EleLColl.at(i), "QCD_Iso05");

      //float FR=FakeRateMC(EleLColl.at(i), "QCD_EGM");
      //float FR=FakeRateMC(EleLColl.at(i), "QCD_WP90TIsop06IPp025p05sig4_LIsop4");
      float FR=FakeRateMC(EleLColl.at(i), "QCD_WP80TIsop1IPp025p05sig4_LIsop4");
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



void Jul2017_MCFakeStudy::Draw2DFakePlot(std::vector<snu::KElectron> FakeColl, std::vector<snu::KJet> JetColl, int NMVACuts, float MVACuts[], int NIsoCuts, float IsoCuts[], int NPtEdges, float PtEdges[], TString Option){
//Input only fake leptons

  bool JetPtParam  = Option.Contains("JetPt");
  bool ConePtParam = Option.Contains("ConePt");
  if(JetPtParam && ConePtParam) return;
  else if( !JetPtParam && !ConePtParam ) return;

  for(int i=0; i<FakeColl.size(); i++){
       
    float MotherPt    = ConeCorrectedPT(FakeColl.at(i), 0.06);
    int   FakeSrcType = GetFakeLepSrcType(FakeColl.at(i), JetColl);
    if(JetPtParam){
      int NearJetIdx  = GetNearJetIdx(FakeColl.at(i), JetColl);
      if(NearJetIdx==-1) return;
      MotherPt=JetColl.at(NearJetIdx).Pt();
    }

    for(int it_iso=1; it_iso<NIsoCuts; it_iso++){
      for(int it_mva=0; it_mva<NMVACuts-1; it_mva++){
        //Check the only region that is not biased by LepPt Threshold.
        if(MotherPt<(25.+FakeColl.at(i).PFRelIso(0.3)*IsoCuts[it_iso])) continue;

        if(FakeColl.at(i).PFRelIso(0.3)<IsoCuts[it_iso] && FakeColl.at(i).MVA()>MVACuts[it_mva]){

          if(fabs(FakeColl.at(i).Eta())<0.8){//Barrel1
            //All fake regardless of near jets.
            FillHist("EleB1SumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
            if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleB1IDSumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);

            if(FakeSrcType==3){
              FillHist("EleB1SumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleB1IDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }//End of B-nearjet loop
            if(FakeSrcType==2){
              FillHist("EleB1SumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleB1IDSumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }//End of C-near jet loop
            if(FakeSrcType==1){
              FillHist("EleB1SumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleB1IDSumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }//End of L-near jet loop
          }
          else if(fabs(FakeColl.at(i).Eta())<1.479){//Barrel2
            //All fake regardless of near jets.
            FillHist("EleB2SumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
            if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleB2IDSumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);

            if(FakeSrcType==3){
              FillHist("EleB2SumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleB2IDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }//End of B-nearjet loop
            if(FakeSrcType==2){
              FillHist("EleB2SumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleB2IDSumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }//End of C-near jet loop
            if(FakeSrcType==1){
              FillHist("EleB2SumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleB2IDSumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }//End of L-near jet loop
          }
          else if(fabs(FakeColl.at(i).Eta())<2.5){//EndCap
            //All fake regardless of near jets.
            FillHist("EleESumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
            if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleEIDSumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);

            if(FakeSrcType==3){
              FillHist("EleESumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleEIDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }//End of B-nearjet loop
            if(FakeSrcType==2){
              FillHist("EleESumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleEIDSumW_CjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }//End of C-near jet loop
            if(FakeSrcType==1){
              FillHist("EleESumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("EleEIDSumW_LjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }//End of L-near jet loop
          }
        }//End of if Iso<Cut && MVA<Cut
      }//End of MVACut loop
    }//End of IsoCut Loop
  }//End of Fake Lepton Loop 

}


void Jul2017_MCFakeStudy::ScanFakeRate(std::vector<snu::KElectron> FakePreColl, std::vector<snu::KJet> JetColl, float MVACut, float IsoCut, int NPtEdges, float PtEdges[], TString PreID, TString TightID, TString Label, TString Option){

  std::vector<snu::KElectron> FakeColl;
    for(int i=0; i<FakePreColl.size(); i++){ if(PassIDCriteria(FakePreColl.at(i),PreID)) FakeColl.push_back(FakePreColl.at(i)); }

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

  float TightIso=0.;
  if     (TightID.Contains("Isop1"))  TightIso=0.1;
  else if(TightID.Contains("Isop08")) TightIso=0.08;
  else if(TightID.Contains("Isop06")) TightIso=0.06;


  for(int i=0; i<FakeColl.size(); i++){
    if(FakeColl.at(i).MVA()<MVACut) continue;
    if(FakeColl.at(i).PFRelIso(0.3)>IsoCut) continue;

    float PTCorr=ConeCorrectedPT(FakeColl.at(i), TightIso);
    int FakeSrcType=GetFakeLepSrcType(FakeColl.at(i), JetColl);

    if(PTCorr<25) continue;
    //AllFakes
    FillHist("Ele"+EtaReg+"SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);   
    if(PassIDCriteria(FakeColl.at(i), TightID)) FillHist("Ele"+EtaReg+"IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

    //PerSources(Only for matched ones)
    if     (FakeSrcType==3){
      FillHist("Ele"+EtaReg+"SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), TightID)) FillHist("Ele"+EtaReg+"IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else if(FakeSrcType==2){
      FillHist("Ele"+EtaReg+"SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_CjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), TightID)) FillHist("Ele"+EtaReg+"IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_CjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else if(FakeSrcType==1){
      FillHist("Ele"+EtaReg+"SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_LjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), TightID)) FillHist("Ele"+EtaReg+"IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_LjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    }
  }//End of Fake Lepton Loop 


}


void Jul2017_MCFakeStudy::Draw1DFakePlot(std::vector<snu::KElectron> FakeColl, std::vector<snu::KJet> JetColl, float MVAB1Cut, float MVAB2Cut, float MVAECut, float IsoCut, int NPtEdges, float PtEdges[], TString Option){

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

    float PTCorr=ConeCorrectedPT(FakeColl.at(i), 0.06);
    int FakeSrcType=GetFakeLepSrcType(FakeColl.at(i), JetColl);

    if(PTCorr<25) continue;

    //AllFakes
    FillHist("Ele"+Option+"SumW_iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
    if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("Ele"+Option+"IDSumW_iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);

    //PerSources(Only for matched ones)
    if     (FakeSrcType==3){
      FillHist("Ele"+Option+"SumW_iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("Ele"+Option+"IDSumW_iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else if(FakeSrcType==2){
      FillHist("Ele"+Option+"SumW_iso"+Str_IsoCut+"_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("Ele"+Option+"IDSumW_iso"+Str_IsoCut+"_CjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else if(FakeSrcType==1){
      FillHist("Ele"+Option+"SumW_iso"+Str_IsoCut+"_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);   
      if(PassIDCriteria(FakeColl.at(i), "POGWP90Isop06IPp025p1sig4")) FillHist("Ele"+Option+"IDSumW_iso"+Str_IsoCut+"_LjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
    }

  }//End of Fake Lepton Loop 


}


void Jul2017_MCFakeStudy::EmulateFRMeasurement(std::vector<snu::KElectron> ElePreColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight, int NPtEdges, float PtEdges[], TString LooseID, TString TightID, TString Label, TString Option){
  //BJet : NoVetoColl,  Jet: No matter whether lepton vetoed or not, it will veto in code any way.


  float MinElePt=25.;

  std::vector<snu::KElectron> EleLColl;
    for(int i=0; i<ElePreColl.size(); i++){ if(PassIDCriteria(ElePreColl.at(i),LooseID)) EleLColl.push_back(ElePreColl.at(i)); }
  std::vector<snu::KJet> JetVetoColl = SkimJetColl(JetColl, EleLColl, MuLColl, "EleMuVeto");


  if(EleLColl.size()==1 && MuLColl.size()==0 && JetVetoColl.size()>0){
    bool PassJetReq=false;

    float PTCorr=ConeCorrectedPT(EleLColl.at(0), 0.06), fEta=fabs(EleLColl.at(0).Eta());

    //Selection Requirement--------------------------------------------------------------//
    if(!PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")) return;
    if( EleLColl.at(0).Pt()<MinElePt ) return;
    if( PTCorr<25. ) return;
    for(int j=0; j<JetVetoColl.size(); j++){
      if(JetVetoColl.at(j).Pt()<40) continue;
      if(JetVetoColl.at(j).DeltaR(EleLColl.at(0))<1.0) continue;
      PassJetReq=true; break;
    }
    if(!PassJetReq) return;
    float MTW = sqrt(2)*sqrt(MET*EleLColl.at(0).Pt()-METx*EleLColl.at(0).Px()-METy*EleLColl.at(0).Py());

    //if( !(MET<25 && MTW<25) ) return;
    //-----------------------------------------------------------------------------------//

    bool IsNearB = IsNearBJet(EleLColl.at(0), BJetColl);
    if(fEta<0.8){
      FillHist("EleB1SumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
   
      if(IsNearB){ FillHist("EleB1SumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB1IDSumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
    }
    else if(fEta<1.479){
      FillHist("EleB2SumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
   
      if(IsNearB){ FillHist("EleB2SumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleB2IDSumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
    }
    else{
      FillHist("EleESumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_All_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
   
      if(IsNearB){ FillHist("EleESumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
                   if(PassIDCriteria(EleLColl.at(0), TightID)) FillHist("EleEIDSumW_BjMatch_PT_FR1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1); } 
    }
  }//End of 1e+0mu+geq1j
 
}


float Jul2017_MCFakeStudy::FakeRateMC(snu::KMuon Mu, TString Option){

  float FR=0., TightIsoCut=0.;
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("TIsop15")) TightIsoCut = 0.15;
    else if(Option.Contains("TIsop20")) TightIsoCut = 0.2;
    else if(Option.Contains("TIsop25")) TightIsoCut = 0.25;
  }
  float PTCorr=ConeCorrectedPT(Mu, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Mu.Pt();
  float fEta=fabs(Mu.Eta());

  if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLFOPt"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.282776;
      else if(PTCorr<20)  FR=0.254416;
      else if(PTCorr<25)  FR=0.228421;
      else if(PTCorr<35)  FR=0.212552;
      else                FR=0.197174;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.299608;
      else if(PTCorr<20)  FR=0.281139;
      else if(PTCorr<25)  FR=0.269316;
      else if(PTCorr<35)  FR=0.249154;
      else                FR=0.258796;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.341341;
      else if(PTCorr<20)  FR=0.328367;
      else if(PTCorr<25)  FR=0.321994;
      else if(PTCorr<35)  FR=0.318312;
      else                FR=0.292036;
    }
    else{
      if     (PTCorr<15)  FR=0.341761;
      else if(PTCorr<20)  FR=0.328543;
      else if(PTCorr<25)  FR=0.320256;
      else if(PTCorr<35)  FR=0.310155;
      else                FR=0.305548;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.31426 ;
      else if(PTCorr<20)  FR=0.222439;
      else if(PTCorr<25)  FR=0.201205;
      else if(PTCorr<35)  FR=0.174831;
      else                FR=0.144662;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.32857 ;
      else if(PTCorr<20)  FR=0.250185;
      else if(PTCorr<25)  FR=0.239277;
      else if(PTCorr<35)  FR=0.204876;
      else                FR=0.200625;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.369924;
      else if(PTCorr<20)  FR=0.298439;
      else if(PTCorr<25)  FR=0.282484;
      else if(PTCorr<35)  FR=0.27146 ;
      else                FR=0.233521;
    }
    else{
      if     (PTCorr<15)  FR=0.368448;
      else if(PTCorr<20)  FR=0.298501;
      else if(PTCorr<25)  FR=0.280768;
      else if(PTCorr<35)  FR=0.27085 ;
      else                FR=0.248833;
    }
  }


  return FR;
}


float Jul2017_MCFakeStudy::FakeRateMC(snu::KElectron Ele, TString Option){

  float FR=0., TightIsoCut=0.;
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("Isop06")) TightIsoCut = 0.06;
    else if(Option.Contains("Isop08")) TightIsoCut = 0.08;
    else if(Option.Contains("Isop1"))  TightIsoCut = 0.1;
  }
  float PTCorr=ConeCorrectedPT(Ele, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Ele.Pt();
  float fEta=fabs(Ele.Eta());

  if(Option=="QCD_POGWP90Isop06IPp025p1sig4_HctoWAFakeLoose_ConeSUSY"){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.0927136;
      else if(PTCorr<50)  FR=0.0720359;
      else if(PTCorr<70)  FR=0.0711901;
      else                FR=0.0624766;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.149674 ;
      else if(PTCorr<50)  FR=0.0613648;
      else if(PTCorr<70)  FR=0.0846761;
      else                FR=0.0534003;
    }
    else{
      if     (PTCorr<35)  FR=0.159323;
      else if(PTCorr<50)  FR=0.104854;
      else if(PTCorr<70)  FR=0.109523;
      else                FR=0.120316;
    }
  }
  else if(Option=="QCD_POGWP90Isop06IPp025p05sig4_LMVA06Isop4IPp025p05sig4_ConeSUSY"){
    if(fEta<0.8){
      if     (PTCorr<35)  FR=0.0933608;
      else if(PTCorr<50)  FR=0.0721879;
      else if(PTCorr<70)  FR=0.0714611;
      else                FR=0.0627619;
    }
    else if(fEta<1.479){
      if     (PTCorr<35)  FR=0.156215 ;
      else if(PTCorr<50)  FR=0.0664036;
      else if(PTCorr<70)  FR=0.0919452;
      else                FR=0.0562619;
    }
    else{
      if     (PTCorr<35)  FR=0.161256;
      else if(PTCorr<50)  FR=0.109329;
      else if(PTCorr<70)  FR=0.114383;
      else                FR=0.126938;
    }
  }

  return FR;
}

int Jul2017_MCFakeStudy::GetFakeLepSrcType(snu::KElectron Ele, std::vector<snu::KJet> JetColl){
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


float Jul2017_MCFakeStudy::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){


  //float PTCorr=Ele.Pt();
  //float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  float PTCorr=Ele.Pt()*(1.+max(0.,Ele.PFRelIso(0.3)-TightIsoCut));

  return PTCorr;
}


float Jul2017_MCFakeStudy::ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut){

  //float PTCorr=Mu.Pt();
  //float PTCorr=Mu.Pt()*(1.+RochIso(Mu,"0.4"));
  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}


int Jul2017_MCFakeStudy::StepPassed(std::vector<snu::KMuon> MuColl, std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, TString Option){

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


int Jul2017_MCFakeStudy::NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> truthColl, TString Option){

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

int Jul2017_MCFakeStudy::NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> truthColl, TString Option){

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



bool Jul2017_MCFakeStudy::IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option){

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


bool Jul2017_MCFakeStudy::IsNearBJet(snu::KElectron Ele, std::vector<snu::KJet> bjetNoVetoColl){

  bool IsNearB=false;
  for(int i=0; i<bjetNoVetoColl.size(); i++){
    if(Ele.DeltaR(bjetNoVetoColl.at(i))<0.4){
      IsNearB=true; break;
    }
  }

  return IsNearB;
}



void Jul2017_MCFakeStudy::FillCutFlow(TString cut, float weight){
  
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



void Jul2017_MCFakeStudy::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Jul2017_MCFakeStudy::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);



  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Jul2017_MCFakeStudyCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Jul2017_MCFakeStudy::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
