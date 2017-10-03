/***************************************************************************
 * @Project: Aug2017_MuFakeMCStudy 
 * @Package: LQanalyzer, ROOT-based analysis framework for Korea SNU / Main package author: John Almond
 *
 * @Code Author: Jihwan Bhyun 
 *
 ***************************************************************************/

/// Local includes
 #include "Aug2017_MuFakeMCStudy.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Aug2017_MuFakeMCStudy);

 Aug2017_MuFakeMCStudy::Aug2017_MuFakeMCStudy() : AnalyzerCore(), out_muons(0) {

   SetLogName("Aug2017_MuFakeMCStudy");
   Message("In Aug2017_MuFakeMCStudy constructor", INFO);
   InitialiseAnalysis();
 }


 void Aug2017_MuFakeMCStudy::InitialiseAnalysis() throw( LQError ) {
   
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

void Aug2017_MuFakeMCStudy::ExecuteEvents()throw( LQError ){

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


  
   bool FakeCompCheck=false, IDEffCheck=false, IDVarSensitivity=false, FRInspection=false, SiglWP=false, ScanFR=false;
   bool TrigBiasCheck=false, FRMeasEmul=false, Closure=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("FakeCompCheck"))    FakeCompCheck    = true;
     else if(k_flags.at(i).Contains("IDEffCheck"))       IDEffCheck       = true;
     else if(k_flags.at(i).Contains("IDVarSensitivity")) IDVarSensitivity = true;
     else if(k_flags.at(i).Contains("FRInspection"))     FRInspection     = true;
     else if(k_flags.at(i).Contains("SiglWP"))           SiglWP           = true;
     else if(k_flags.at(i).Contains("ScanFR"))           ScanFR           = true;
     else if(k_flags.at(i).Contains("TrigBiasCheck"))    TrigBiasCheck    = true;
     else if(k_flags.at(i).Contains("FRMeasEmul"))       FRMeasEmul       = true;
     else if(k_flags.at(i).Contains("Closure"))          Closure          = true;
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

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
   std::vector<snu::KMuon> muonPOGLIsoVLColl; eventbase->GetMuonSel()->Selection(muonPOGLIsoVLColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_MEDIUM);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPOGMColl; eventbase->GetMuonSel()->Selection(muonPOGMColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPOGTColl; eventbase->GetMuonSel()->Selection(muonPOGTColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.15);
   std::vector<snu::KMuon> muonPOGTIsoColl; eventbase->GetMuonSel()->Selection(muonPOGTIsoColl, true);


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
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

   std::vector<snu::KJet> jetColl  = SkimJetColl(jetNoVetoColl,  electronLooseColl, muonHN2FakeLColl, "EleMuVeto");
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



   //=============================================================================//
   // Main Analysis Code     -----------------------------------------------------//
   //=============================================================================//

   if(FakeCompCheck){
     //Purpose : Fake source and types of 1) signal region(3l+b+3j) & trilep(nob)
     //          2) TT and DY in general, 3) QCD measurement region
     //Check Types: 1) b/c/ljet, 2)LeptonType, 3) SourceHadron(Including unmatched)

     CheckFakeSources(muonPreColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_NoIDPFMu", "");
     CheckFakeSources(muonPOGLIsoVLColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_POGLIso", "");
     CheckFakeSources(muonPOGMColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_POGM", "");
     CheckFakeSources(muonPOGTColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_POGT", "");
     CheckFakeSources(muonPOGTIsoColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_POGTIso", "");

   }
   if(IDEffCheck){
     //Purpose : Check overall ID performance

     CheckIDEfficiency(muonPreColl, truthColl, weight, "", "");

   }
   if(IDVarSensitivity){
     //Purpose : Which ID variable is effective to which kinds of fakes.

     std::vector<snu::KMuon> muonPOGLIsop6IPpXp1;
     for(int i=0; i<muonPreColl.size(); i++){
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop6IPpXp1","Roch")) muonPOGLIsop6IPpXp1.push_back(muonPreColl.at(i));
     }

     CheckIDVarSensitivity(muonPOGLIsop6IPpXp1, jetNoVetoColl, truthColl, weight, "", "");

     std::vector<snu::KMuon> muon10to20Coll, muon20to50Coll, muon50to100Coll;
     for(int i=0; i<muonPOGLIsop6IPpXp1.size(); i++){
       if     (muonPOGLIsop6IPpXp1.at(i).Pt()<20)  muon10to20Coll.push_back(muonPOGLIsop6IPpXp1.at(i));
       else if(muonPOGLIsop6IPpXp1.at(i).Pt()<50)  muon20to50Coll.push_back(muonPOGLIsop6IPpXp1.at(i));
       else if(muonPOGLIsop6IPpXp1.at(i).Pt()<100) muon50to100Coll.push_back(muonPOGLIsop6IPpXp1.at(i));
     }

     CheckIDVarSensitivity(muon10to20Coll , jetNoVetoColl, truthColl, weight, "_10to20" , "");
     CheckIDVarSensitivity(muon20to50Coll , jetNoVetoColl, truthColl, weight, "_20to50" , "");
     CheckIDVarSensitivity(muon50to100Coll, jetNoVetoColl, truthColl, weight, "_50to100", "");
   }
   if(FRInspection){
     //Purpose : Check FR of various scenarios.

     if(SiglWP){
       InspectFakeRate(muonPreColl, jetNoVetoColl, truthColl,
         "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGLIsop4IPp5p1Chi100", "");
       InspectFakeRate(muonPreColl, jetNoVetoColl, truthColl,
         "Test_POGLIsop4IPp5p1", "Test_POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGLIsop4IPp5p1", "");
       InspectFakeRate(muonPreColl, jetNoVetoColl, truthColl,
         "Test_POGLIsop6IPp5p1", "Test_POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGLIsop6IPp5p1", "");
     }
     if(ScanFR){
       //const int Nd0Cuts  =12; float d0Cuts[Nd0Cuts]    ={0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.1, 0.2, 0.5, 1., 10., 100.};
       //const int Nd0Cuts   =1;  float d0Cuts[Nd0Cuts]    ={0.5};
       const int NChi2Cuts =8; float Chi2Cuts[NChi2Cuts]  ={4., 10., 20., 30., 50., 80., 100., 999.};
       const int Nd0SigCuts=11; float d0SigCuts[Nd0SigCuts]={4., 5., 7., 8., 9., 10., 20., 30., 50., 100., 999.};

       ScanFakeRate(muonPreColl, electronLooseColl, jetNoVetoColl, truthColl, "Test_POGLIsop6IPp5p1", "Test_POGTIsop20IPp01p05sig4Chi4", Nd0SigCuts, d0SigCuts, NChi2Cuts, Chi2Cuts, weight, "_TIsop20IPp01p05sig4Chi4LIsop6", "");
       ScanFakeRate(muonPreColl, electronLooseColl, jetNoVetoColl, truthColl, "Test_POGLIsop4IPp5p1", "Test_POGTIsop20IPp01p05sig4Chi4", Nd0SigCuts, d0SigCuts, NChi2Cuts, Chi2Cuts, weight, "_TIsop20IPp01p05sig4Chi4LIsop4", "");

     }

   }
   if(TrigBiasCheck){
     //Purpose : Check trigger bias. 1) Maximum unbiased extrapolation region : how much can I loosen ID variable?
     //                              2) Does trigger requirement change FR even in trigger safe range?

     //CheckTriggerBias(muonPreColl, jetNoVetoColl, truthColl, "Test_POGLIsop5IPp5p1Chi", "Test_POGTIsop20IPp01p1Chi3", weight, "", "");
     CheckTriggerBias(muonPreColl, jetNoVetoColl, truthColl, "POGLNoIso", "Test_POGTIsop20IPp01p1Chi3", weight, "", "");
   }
   if(FRMeasEmul){
     //Purpose : Fully emulate FR measurement for QCD. Result will be used for closure to DY, TT

     //std::vector<snu::KMuon> muonPOGLIsop5IPp5p1ChiColl, muonPOGLIsop6IPp5p1Chi30Coll;
     std::vector<snu::KMuon> muonPOGLIsop4IPp5p1Coll, muonPOGLIsop4IPp5p1Chi100Coll, muonPOGLIsop6IPp5p1Coll;
     for(int i=0; i<muonPreColl.size(); i++){
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop4IPp5p1","Roch")) muonPOGLIsop4IPp5p1Coll.push_back(muonPreColl.at(i));
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop4IPp5p1Chi100","Roch")) muonPOGLIsop4IPp5p1Chi100Coll.push_back(muonPreColl.at(i));
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop6IPp5p1","Roch")) muonPOGLIsop6IPp5p1Coll.push_back(muonPreColl.at(i));
       //if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop6IPp5p1Chi30","Roch")) muonPOGLIsop6IPp5p1Chi30Coll.push_back(muonPreColl.at(i));
       //if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop5IPp5p1Chi"  ,"Roch")) muonPOGLIsop5IPp5p1ChiColl.push_back(muonPreColl.at(i));
     }


     EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p05sig4Chi4", weight, "_TIsop20IPp01p05sig4Chi4LIsop4IPp5p1", "TrkIsoVVL");
     EmulateFRMeasurement(muonPOGLIsop4IPp5p1Chi100Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p05sig4Chi4", weight, "_TIsop20IPp01p05sig4Chi4LIsop4IPp5p1Chi100", "TrkIsoVVL");
     EmulateFRMeasurement(muonPOGLIsop6IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p05sig4Chi4", weight, "_TIsop20IPp01p05sig4Chi4Lp6", "TrkIsoVVL");

     //Round2
     //2-1
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p05sig4Chi3", weight, "_TIsop20IPp01p05sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop6IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p05sig4Chi3", weight, "_TIsop20IPp01p05sig4Lp6", "TrkIsoVVL");
     //2-2
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p025sig4Chi4", weight, "_TIsop20IPp01p025sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p05sig4Chi4" , weight, "_TIsop20IPp01p05sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p075sig4Chi4", weight, "_TIsop20IPp01p075sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p1sig4Chi4"  , weight, "_TIsop20IPp01p1sig4Lp4", "TrkIsoVVL");

     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp0125p04sig4Chi4", weight, "_TIsop20IPp0125p04sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp0125p05sig4Chi4", weight, "_TIsop20IPp0125p05sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp0125p06sig4Chi4", weight, "_TIsop20IPp0125p06sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp0125p07sig4Chi4", weight, "_TIsop20IPp0125p07sig4Lp4", "TrkIsoVVL");

     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp015p04sig4Chi4", weight, "_TIsop20IPp015p04sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp015p05sig4Chi4", weight, "_TIsop20IPp015p05sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp015p06sig4Chi4", weight, "_TIsop20IPp015p06sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp015p07sig4Chi4", weight, "_TIsop20IPp015p07sig4Lp4", "TrkIsoVVL");

     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp0175p04sig4Chi4", weight, "_TIsop20IPp0175p04sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp0175p05sig4Chi4", weight, "_TIsop20IPp0175p05sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp0175p06sig4Chi4", weight, "_TIsop20IPp0175p06sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp0175p07sig4Chi4", weight, "_TIsop20IPp0175p07sig4Lp4", "TrkIsoVVL");

     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp02p04sig4Chi4", weight, "_TIsop20IPp02p04sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp02p05sig4Chi4", weight, "_TIsop20IPp02p05sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp02p06sig4Chi4", weight, "_TIsop20IPp02p06sig4Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp02p07sig4Chi4", weight, "_TIsop20IPp02p07sig4Lp4", "TrkIsoVVL");


     //2-3; Chi impact on slope
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p05sig4" , weight, "_TIsop20IPp01p05sig4Chi10Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp0125p05sig4", weight, "_TIsop20IPp0125p05sig4Chi10Lp4", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop4IPp5p1Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp015p05sig4", weight, "_TIsop20IPp015p05sig4Chi10Lp4", "TrkIsoVVL");


     //Round 1
     //EmulateFRMeasurement(muonPOGLIsop6IPp5p1Chi30Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p1Chi3", weight, "_TIsop20Lp6", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop5IPp5p1ChiColl  , electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop20IPp01p1Chi3", weight, "_TIsop20Lp5", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop6IPp5p1Chi30Coll, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop15IPp01p1Chi3", weight, "_TIsop15Lp6", "TrkIsoVVL");
     //EmulateFRMeasurement(muonPOGLIsop5IPp5p1ChiColl  , electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "Test_POGTIsop15IPp01p1Chi3", weight, "_TIsop15Lp5", "TrkIsoVVL");
     //EmulateFRMeasurement(muonHN2FakeLColl            , electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, "HNTrilepTight2", weight, "_HNTrilep", "TrkIsoVVL");

   }
   if(Closure){
     //Purpose : Description of yields in different sample expected from QCD measured FR

     CheckMCClosure(muonPreColl, electronLooseColl, jetNoVetoColl, truthColl, "Test_POGLIsop4IPp5p1Chi100", "Test_POGTIsop20IPp01p05sig4Chi4", weight, "_TIsop20IPp01p05sig4Chi4LIsop4IPp5p1Chi100", "TrkIsoVVLConeSUSY");   
     

     //#"TrkIsoVVLFOPt" #"TrkIsoVVLConeSUSY" #"TrkIsoVVLConeE" #"NoFilterFOPt" #"NoFilterConeSUSY" #"NoFilterConeE"

   }

return;
}// End of execute event loop
  

void Aug2017_MuFakeMCStudy::CheckFakeSources(std::vector<snu::KMuon> muonColl, std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KElectron> electronLooseColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, std::vector<snu::KTruth> truthColl, float weight, TString Label, TString Option){
  //Input JetColl, BJetColl as novetocoll since I may use nearby jet.
  //Purpose : Fake source and types of 1) signal region(3l+b+3j) & trilep(nob)
  //          2) TT and DY in general, 3) QCD measurement region
  //Check Types: 1) b/c/ljet, 2)LeptonType, 3) SourceHadron(Including unmatched)

  //std::vector<snu::KElectron> FakeColl     = SkimLepColl(muonLooseColl, truthColl, "HFake");
  //std::vector<snu::KJet>      JetCleanColl = SkimJetColl(jetColl, truthColl, "NoPrNoTau");

  std::vector<snu::KJet> jetVetoColl  = SkimJetColl(JetColl,  electronLooseColl, muonLooseColl, "EleMuVeto");
  std::vector<snu::KJet> bjetVetoColl = SkimJetColl(BJetColl, electronLooseColl, muonLooseColl, "EleMuVeto");
  bool TrilepNob=false, TrilepHasB=false, SRSel=false, FRMeasReg=false;
  bool TrilepNob_L=false, TrilepHasB_L=false, SRSel_L=false, FRMeasReg_L=false;
  bool PassMeasJetReq=false;

  for(int i=0; i<muonColl.size(); i++){
    int LepType    = GetLeptonType(muonColl.at(i),truthColl);
    int SrcJetType = GetFakeLepJetSrcType(muonColl.at(i), JetColl);
    int SrcIdx     = GetFakeLepSrcIdx(muonColl.at(i), truthColl);
    int SrcfPID    = SrcIdx==-1 ? 0 : fabs(truthColl.at(SrcIdx).PdgId());

    if(LepType>0 && LepType<4) continue;

    //Overall LepTypes
    FillHist("LepType"+Label, LepType, weight, -10., 10., 20);

    //SrcJetType Composition and its dependence on flavour //e)
    FillHist("SrcJetType"+Label, SrcJetType, weight, -1., 4., 5);
    if(SrcJetType==1) FillHist("LepType_Lj"+Label, LepType, weight, -10., 10., 20);   
    if(SrcJetType==2) FillHist("LepType_Cj"+Label, LepType, weight, -10., 10., 20);   
    if(SrcJetType==3) FillHist("LepType_Bj"+Label, LepType, weight, -10., 10., 20);   

    //Source hadrons for each fake type of leptons
    if(LepType==-1){
      FillHist("LepTypem1_MatchFr"+Label, 0., weight, 0., 2., 2);
      if(SrcIdx!=-1){
        FillHist("LepTypem1_MatchFr"+Label, 1., weight, 0., 2., 2);//a)

        //b)
        float dRmufo=muonColl.at(i).DeltaR(truthColl.at(SrcIdx));
        float dPtRel=fabs(muonColl.at(i).Pt()-truthColl.at(SrcIdx).Pt())/muonColl.at(i).Pt(); 
        FillHist("dRmufo_Typem1"+Label, dRmufo, weight, 0., 0.5, 100);
        FillHist("SrcPID_Typem1"+Label, SrcfPID, weight, 0., 6000., 6000);
        if(SrcfPID==211){ FillHist("dRmufo_Typem1_pid211"+Label, dRmufo, weight, 0., 0.5, 100);
                         FillHist("dPtRelmufo_Typem1_pid211"+Label, dPtRel, weight, 0., 10., 100); }
        if(SrcfPID==310){ FillHist("dRmufo_Typem1_pid310"+Label, dRmufo, weight, 0., 0.5, 100);
                         FillHist("dPtRelmufo_Typem1_pid310"+Label, dPtRel, weight, 0., 10., 100); }
        if(SrcfPID==2212){ FillHist("dRmufo_Typem1_pid2212"+Label, dRmufo, weight, 0., 0.5, 100);
                          FillHist("dPtRelmufo_Typem1_pid2212"+Label, dPtRel, weight, 0., 10., 100); }
      }
    }
    else if(LepType==-2 || LepType==-3){//c)//d)
      if(SrcIdx!=-1){
        FillHist("SrcPID_Typem2"+Label, SrcfPID, weight, 0., 6000., 6000);
      }
    }
//a)How many unmatched fakes have nearby hadrons? ->Only 30%.. Can be a little biased study for TT
//                                                  but TightIso DY case is nearly ~60% so it can be meaningful
//b)It is observed that for non-ID PFMuon, Neutral Kaon(70%), Proton(20%), Charged Pion(10%) are most dominant fake sources
//  for unmatched cases. But isolated fakes are mostly pions regardless of ID-tightness.
//  And it is also observed that pion momentum is nearly the same with muon in dR, Pt
//  But kaon and proton case, momentum is hugely different suggests muon reco due to misreconstruction
//c)Mother of Type-2,-3 are observed to be similar, so merged
//d)B Mesons are dominant for TT(~60%), but still ~30% comes from charm meson as D0, D+..
//e)hadron decays are main source for b,c jet, but still unmatched is the main source for light jets
//  and in DY light jet is the biggest source of fakes just as electron. especially pions in light jet.

  } 

}



void Aug2017_MuFakeMCStudy::CheckIDEfficiency(std::vector<snu::KMuon> muonColl, std::vector<snu::KTruth> truthColl, float weight, TString Label, TString Option){

  std::vector<snu::KMuon> PromptColl   = SkimLepColl(muonColl, truthColl, "Prompt");
  std::vector<snu::KMuon> FakeColl     = SkimLepColl(muonColl, truthColl, "HFake");


  for(int i=0; i<PromptColl.size(); i++){
    //Inclusive
    FillHist("NMuSumW_Prompt"+Label, 0., weight, 0., 10., 10);
    FillHist("NMuSumW_Prompt"+Label, 1., weight, 0., 10., 10);
    FillHist("NMuSumW_Prompt"+Label, 2., weight, 0., 10., 10);
    FillHist("NMuSumW_Prompt"+Label, 3., weight, 0., 10., 10);
    FillHist("NMuSumW_Prompt"+Label, 4., weight, 0., 10., 10);
    FillHist("NMuSumW_Prompt"+Label, 5., weight, 0., 10., 10);
    FillHist("NMuSumW_Prompt"+Label, 6., weight, 0., 10., 10);
    FillHist("NMuSumW_Prompt"+Label, 7., weight, 0., 10., 10);
    FillHist("NMuSumW_Prompt"+Label, 8., weight, 0., 10., 10);
    FillHist("NMuSumW_Prompt"+Label, 9., weight, 0., 10., 10);
    if(PromptColl.at(i).IsMedium()) FillHist("NMuIDSumW_Prompt"+Label, 0., weight, 0., 10., 10);
    if(PromptColl.at(i).IsTight())  FillHist("NMuIDSumW_Prompt"+Label, 1., weight, 0., 10., 10);
    if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.6 ) FillHist("NMuIDSumW_Prompt"+Label, 2., weight, 0., 10., 10);
    if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.15) FillHist("NMuIDSumW_Prompt"+Label, 3., weight, 0., 10., 10);
    if(PassIDCriteria(PromptColl.at(i), "HNTrilepFakeL2"))         FillHist("NMuIDSumW_Prompt"+Label, 4., weight, 0., 10., 10);
    if(PassIDCriteria(PromptColl.at(i), "HNTrilepTight2"))         FillHist("NMuIDSumW_Prompt"+Label, 5., weight, 0., 10., 10);
    if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop6IPpXp1"))   FillHist("NMuIDSumW_Prompt"+Label, 6., weight, 0., 10., 10);
    if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop4IPpXp1"))   FillHist("NMuIDSumW_Prompt"+Label, 7., weight, 0., 10., 10);
    if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop15IPp01p1")) FillHist("NMuIDSumW_Prompt"+Label, 8., weight, 0., 10., 10);
    if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop25IPp01p1")) FillHist("NMuIDSumW_Prompt"+Label, 9., weight, 0., 10., 10);



    //LowPt Efficiency
    if(PromptColl.at(i).Pt()<15){
      FillHist("NMuSumW_PromptPt10to15"+Label, 0., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt10to15"+Label, 1., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt10to15"+Label, 2., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt10to15"+Label, 3., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt10to15"+Label, 4., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt10to15"+Label, 5., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt10to15"+Label, 6., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt10to15"+Label, 7., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt10to15"+Label, 8., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt10to15"+Label, 9., weight, 0., 10., 10);
      if(PromptColl.at(i).IsMedium()) FillHist("NMuIDSumW_PromptPt10to15"+Label, 0., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight())  FillHist("NMuIDSumW_PromptPt10to15"+Label, 1., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.6 ) FillHist("NMuIDSumW_PromptPt10to15"+Label, 2., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.15) FillHist("NMuIDSumW_PromptPt10to15"+Label, 3., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "HNTrilepFakeL2"))         FillHist("NMuIDSumW_PromptPt10to15"+Label, 4., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "HNTrilepTight2"))         FillHist("NMuIDSumW_PromptPt10to15"+Label, 5., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop6IPpXp1"))   FillHist("NMuIDSumW_PromptPt10to15"+Label, 6., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop4IPpXp1"))   FillHist("NMuIDSumW_PromptPt10to15"+Label, 7., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop15IPp01p1")) FillHist("NMuIDSumW_PromptPt10to15"+Label, 8., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop25IPp01p1")) FillHist("NMuIDSumW_PromptPt10to15"+Label, 9., weight, 0., 10., 10);
    }
    else if(PromptColl.at(i).Pt()<20){
      FillHist("NMuSumW_PromptPt15to20"+Label, 0., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt15to20"+Label, 1., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt15to20"+Label, 2., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt15to20"+Label, 3., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt15to20"+Label, 4., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt15to20"+Label, 5., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt15to20"+Label, 6., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt15to20"+Label, 7., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt15to20"+Label, 8., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt15to20"+Label, 9., weight, 0., 10., 10);
      if(PromptColl.at(i).IsMedium()) FillHist("NMuIDSumW_PromptPt15to20"+Label, 0., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight())  FillHist("NMuIDSumW_PromptPt15to20"+Label, 1., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.6 ) FillHist("NMuIDSumW_PromptPt15to20"+Label, 2., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.15) FillHist("NMuIDSumW_PromptPt15to20"+Label, 3., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "HNTrilepFakeL2"))         FillHist("NMuIDSumW_PromptPt15to20"+Label, 4., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "HNTrilepTight2"))         FillHist("NMuIDSumW_PromptPt15to20"+Label, 5., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop6IPpXp1"))   FillHist("NMuIDSumW_PromptPt15to20"+Label, 6., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop4IPpXp1"))   FillHist("NMuIDSumW_PromptPt15to20"+Label, 7., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop15IPp01p1")) FillHist("NMuIDSumW_PromptPt15to20"+Label, 8., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop25IPp01p1")) FillHist("NMuIDSumW_PromptPt15to20"+Label, 9., weight, 0., 10., 10);
    }
    else if(PromptColl.at(i).Pt()<30){
      FillHist("NMuSumW_PromptPt20to30"+Label, 0., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt20to30"+Label, 1., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt20to30"+Label, 2., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt20to30"+Label, 3., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt20to30"+Label, 4., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt20to30"+Label, 5., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt20to30"+Label, 6., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt20to30"+Label, 7., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt20to30"+Label, 8., weight, 0., 10., 10);
      FillHist("NMuSumW_PromptPt20to30"+Label, 9., weight, 0., 10., 10);
      if(PromptColl.at(i).IsMedium()) FillHist("NMuIDSumW_PromptPt20to30"+Label, 0., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight())  FillHist("NMuIDSumW_PromptPt20to30"+Label, 1., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.6 ) FillHist("NMuIDSumW_PromptPt20to30"+Label, 2., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.15) FillHist("NMuIDSumW_PromptPt20to30"+Label, 3., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "HNTrilepFakeL2"))         FillHist("NMuIDSumW_PromptPt20to30"+Label, 4., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "HNTrilepTight2"))         FillHist("NMuIDSumW_PromptPt20to30"+Label, 5., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop6IPpXp1"))   FillHist("NMuIDSumW_PromptPt20to30"+Label, 6., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop4IPpXp1"))   FillHist("NMuIDSumW_PromptPt20to30"+Label, 7., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop15IPp01p1")) FillHist("NMuIDSumW_PromptPt20to30"+Label, 8., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop25IPp01p1")) FillHist("NMuIDSumW_PromptPt20to30"+Label, 9., weight, 0., 10., 10);
    }

    if(GetLeptonType(PromptColl.at(i),truthColl)==2){
      //Inclusive
      FillHist("NMuSumW_BSM"+Label, 0., weight, 0., 10., 10);
      FillHist("NMuSumW_BSM"+Label, 1., weight, 0., 10., 10);
      FillHist("NMuSumW_BSM"+Label, 2., weight, 0., 10., 10);
      FillHist("NMuSumW_BSM"+Label, 3., weight, 0., 10., 10);
      FillHist("NMuSumW_BSM"+Label, 4., weight, 0., 10., 10);
      FillHist("NMuSumW_BSM"+Label, 5., weight, 0., 10., 10);
      FillHist("NMuSumW_BSM"+Label, 6., weight, 0., 10., 10);
      FillHist("NMuSumW_BSM"+Label, 7., weight, 0., 10., 10);
      FillHist("NMuSumW_BSM"+Label, 8., weight, 0., 10., 10);
      FillHist("NMuSumW_BSM"+Label, 9., weight, 0., 10., 10);
      if(PromptColl.at(i).IsMedium()) FillHist("NMuIDSumW_BSM"+Label, 0., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight())  FillHist("NMuIDSumW_BSM"+Label, 1., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.6 ) FillHist("NMuIDSumW_BSM"+Label, 2., weight, 0., 10., 10);
      if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.15) FillHist("NMuIDSumW_BSM"+Label, 3., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "HNTrilepFakeL2"))         FillHist("NMuIDSumW_BSM"+Label, 4., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "HNTrilepTight2"))         FillHist("NMuIDSumW_BSM"+Label, 5., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop6IPpXp1"))   FillHist("NMuIDSumW_BSM"+Label, 6., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop4IPpXp1"))   FillHist("NMuIDSumW_BSM"+Label, 7., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop15IPp01p1")) FillHist("NMuIDSumW_BSM"+Label, 8., weight, 0., 10., 10);
      if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop25IPp01p1")) FillHist("NMuIDSumW_BSM"+Label, 9., weight, 0., 10., 10);
  
      //LowPt Efficiency
      if(PromptColl.at(i).Pt()<15){
        FillHist("NMuSumW_BSM10to15"+Label, 0., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM10to15"+Label, 1., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM10to15"+Label, 2., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM10to15"+Label, 3., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM10to15"+Label, 4., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM10to15"+Label, 5., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM10to15"+Label, 6., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM10to15"+Label, 7., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM10to15"+Label, 8., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM10to15"+Label, 9., weight, 0., 10., 10);
        if(PromptColl.at(i).IsMedium()) FillHist("NMuIDSumW_BSM10to15"+Label, 0., weight, 0., 10., 10);
        if(PromptColl.at(i).IsTight())  FillHist("NMuIDSumW_BSM10to15"+Label, 1., weight, 0., 10., 10);
        if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.6 ) FillHist("NMuIDSumW_BSM10to15"+Label, 2., weight, 0., 10., 10);
        if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.15) FillHist("NMuIDSumW_BSM10to15"+Label, 3., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "HNTrilepFakeL2"))         FillHist("NMuIDSumW_BSM10to15"+Label, 4., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "HNTrilepTight2"))         FillHist("NMuIDSumW_BSM10to15"+Label, 5., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop6IPpXp1"))   FillHist("NMuIDSumW_BSM10to15"+Label, 6., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop4IPpXp1"))   FillHist("NMuIDSumW_BSM10to15"+Label, 7., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop15IPp01p1")) FillHist("NMuIDSumW_BSM10to15"+Label, 8., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop25IPp01p1")) FillHist("NMuIDSumW_BSM10to15"+Label, 9., weight, 0., 10., 10);
      }
      else if(PromptColl.at(i).Pt()<20){
        FillHist("NMuSumW_BSM15to20"+Label, 0., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM15to20"+Label, 1., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM15to20"+Label, 2., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM15to20"+Label, 3., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM15to20"+Label, 4., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM15to20"+Label, 5., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM15to20"+Label, 6., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM15to20"+Label, 7., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM15to20"+Label, 8., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM15to20"+Label, 9., weight, 0., 10., 10);
        if(PromptColl.at(i).IsMedium()) FillHist("NMuIDSumW_BSM15to20"+Label, 0., weight, 0., 10., 10);
        if(PromptColl.at(i).IsTight())  FillHist("NMuIDSumW_BSM15to20"+Label, 1., weight, 0., 10., 10);
        if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.6 ) FillHist("NMuIDSumW_BSM15to20"+Label, 2., weight, 0., 10., 10);
        if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.15) FillHist("NMuIDSumW_BSM15to20"+Label, 3., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "HNTrilepFakeL2"))         FillHist("NMuIDSumW_BSM15to20"+Label, 4., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "HNTrilepTight2"))         FillHist("NMuIDSumW_BSM15to20"+Label, 5., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop6IPpXp1"))   FillHist("NMuIDSumW_BSM15to20"+Label, 6., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop4IPpXp1"))   FillHist("NMuIDSumW_BSM15to20"+Label, 7., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop15IPp01p1")) FillHist("NMuIDSumW_BSM15to20"+Label, 8., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop25IPp01p1")) FillHist("NMuIDSumW_BSM15to20"+Label, 9., weight, 0., 10., 10);
      }
      else if(PromptColl.at(i).Pt()<30){
        FillHist("NMuSumW_BSM20to30"+Label, 0., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM20to30"+Label, 1., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM20to30"+Label, 2., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM20to30"+Label, 3., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM20to30"+Label, 4., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM20to30"+Label, 5., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM20to30"+Label, 6., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM20to30"+Label, 7., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM20to30"+Label, 8., weight, 0., 10., 10);
        FillHist("NMuSumW_BSM20to30"+Label, 9., weight, 0., 10., 10);
        if(PromptColl.at(i).IsMedium()) FillHist("NMuIDSumW_BSM20to30"+Label, 0., weight, 0., 10., 10);
        if(PromptColl.at(i).IsTight())  FillHist("NMuIDSumW_BSM20to30"+Label, 1., weight, 0., 10., 10);
        if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.6 ) FillHist("NMuIDSumW_BSM20to30"+Label, 2., weight, 0., 10., 10);
        if(PromptColl.at(i).IsTight() && PromptColl.at(i).RelIso04()<0.15) FillHist("NMuIDSumW_BSM20to30"+Label, 3., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "HNTrilepFakeL2"))         FillHist("NMuIDSumW_BSM20to30"+Label, 4., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "HNTrilepTight2"))         FillHist("NMuIDSumW_BSM20to30"+Label, 5., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop6IPpXp1"))   FillHist("NMuIDSumW_BSM20to30"+Label, 6., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGLIsop4IPpXp1"))   FillHist("NMuIDSumW_BSM20to30"+Label, 7., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop15IPp01p1")) FillHist("NMuIDSumW_BSM20to30"+Label, 8., weight, 0., 10., 10);
        if(PassIDCriteria(PromptColl.at(i), "Test_POGTIsop25IPp01p1")) FillHist("NMuIDSumW_BSM20to30"+Label, 9., weight, 0., 10., 10);
      }
    }//End of BSM Prompt
  }//End of PromptColl
  for(int i=0; i<FakeColl.size(); i++){
    FillHist("NMuSumW_Fake"+Label, 0., weight, 0., 10., 10);
    FillHist("NMuSumW_Fake"+Label, 1., weight, 0., 10., 10);
    FillHist("NMuSumW_Fake"+Label, 2., weight, 0., 10., 10);
    FillHist("NMuSumW_Fake"+Label, 3., weight, 0., 10., 10);
    FillHist("NMuSumW_Fake"+Label, 4., weight, 0., 10., 10);
    FillHist("NMuSumW_Fake"+Label, 5., weight, 0., 10., 10);
    FillHist("NMuSumW_Fake"+Label, 6., weight, 0., 10., 10);
    FillHist("NMuSumW_Fake"+Label, 7., weight, 0., 10., 10);
    FillHist("NMuSumW_Fake"+Label, 8., weight, 0., 10., 10);
    FillHist("NMuSumW_Fake"+Label, 9., weight, 0., 10., 10);
    if(FakeColl.at(i).IsMedium()) FillHist("NMuIDSumW_Fake"+Label, 0., weight, 0., 10., 10);
    if(FakeColl.at(i).IsTight())  FillHist("NMuIDSumW_Fake"+Label, 1., weight, 0., 10., 10);
    if(FakeColl.at(i).IsTight() && FakeColl.at(i).RelIso04()<0.6 ) FillHist("NMuIDSumW_Fake"+Label, 2., weight, 0., 10., 10);
    if(FakeColl.at(i).IsTight() && FakeColl.at(i).RelIso04()<0.15) FillHist("NMuIDSumW_Fake"+Label, 3., weight, 0., 10., 10);
    if(PassIDCriteria(FakeColl.at(i), "HNTrilepFakeL2"))         FillHist("NMuIDSumW_Fake"+Label, 4., weight, 0., 10., 10);
    if(PassIDCriteria(FakeColl.at(i), "HNTrilepTight2"))         FillHist("NMuIDSumW_Fake"+Label, 5., weight, 0., 10., 10);
    if(PassIDCriteria(FakeColl.at(i), "Test_POGLIsop6IPpXp1"))   FillHist("NMuIDSumW_Fake"+Label, 6., weight, 0., 10., 10);
    if(PassIDCriteria(FakeColl.at(i), "Test_POGLIsop4IPpXp1"))   FillHist("NMuIDSumW_Fake"+Label, 7., weight, 0., 10., 10);
    if(PassIDCriteria(FakeColl.at(i), "Test_POGTIsop15IPp01p1")) FillHist("NMuIDSumW_Fake"+Label, 8., weight, 0., 10., 10);
    if(PassIDCriteria(FakeColl.at(i), "Test_POGTIsop25IPp01p1")) FillHist("NMuIDSumW_Fake"+Label, 9., weight, 0., 10., 10);
  }

//Conclusion
//1. Low PT ID efficiency of signal muon is same as DY. (55~60%: PT10-15, ~70%: PT15-20, ~80%: PT20-30 @POGT+iso04)
//   It seems it's just my signal muon efficiency is low because pt's low. Tight Iso is critical for both DY, my signal
//   (But my signal eff is lower by ~10% of DY(i.e. 5~8% depending on PT) for the same cond. So evt signature matters to some extent)
//   It's best if we can use lower iso. (but it can be fatal due to increase in fake contribution and reduction of extrapolation space for fake est.)
//2. Tight ID is better than Medium ID.
//   POGT prompt eff only smaller by ~1% than medium, but fake rate(DY; usually light fake) is smaller by 20%
//   Cuts as Nstation are actually quite much important for fakes as pions
//3. ID efficiency is ~100% for secondary fakes just as prompts. It can only be controlled with Iso, IP
//4. Iso Eff. works both for light fakes(unmatched) and secondary leptons well. 
//   But iso fake rate is ~2 times higher in l-fake(~10% vs. 20%)
//   Instead I observe ID efficiency is ~50% for light fakes. So if IP can be well tuned, fakes can be controlled with them.
}


void Aug2017_MuFakeMCStudy::CheckIDVarSensitivity(std::vector<snu::KMuon> muonColl, std::vector<snu::KJet> JetColl, std::vector<snu::KTruth> truthColl, float weight, TString Label, TString Option){

  for(int i=0; i<muonColl.size(); i++){
    int LepType    = GetLeptonType(muonColl.at(i),truthColl);
    int SrcJetType = GetFakeLepJetSrcType(muonColl.at(i), JetColl);
    int SrcIdx     = GetFakeLepSrcIdx(muonColl.at(i), truthColl);
    int SrcfPID    = SrcIdx==-1 ? 0 : fabs(truthColl.at(SrcIdx).PdgId());
    bool IsHFake=false, IsPrompt=false;

    if( LepType<0 && LepType>-5 ) IsHFake=true;
    else if(LepType>0 && LepType<4) IsPrompt=true;
    else continue;
    if( !muonColl.at(i).IsLoose() || muonColl.at(i).RelIso04()>0.6) continue;
   
    if(IsHFake){
      FillHist("IsGlobal_HFk"+Label, 0., weight, 0., 2., 2);
      if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_HFk"+Label, 1., weight, 0., 2., 2); 
      FillHist("NValidHits_HFk"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
      FillHist("NValidPixHits_HFk"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
      FillHist("NValidStations_HFk"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
      FillHist("NActiveLayer_HFk"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
      FillHist("GlobalChi2_HFk"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 50., 50);
      FillHist("D0_HFk"+Label, muonColl.at(i).dXY(), weight, 0., 1., 1000);
      FillHist("DZ_HFk"+Label, muonColl.at(i).dZ(), weight, 0., 1., 1000); 
      FillHist("RelIso04_HFk"+Label, muonColl.at(i).RelIso04(), weight, 0., 0.5, 50); 
    }
    if(IsPrompt){
      FillHist("IsGlobal_Pr"+Label, 0., weight, 0., 2., 2);
      if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_Pr"+Label, 1., weight, 0., 2., 2); 
      FillHist("NValidHits_Pr"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
      FillHist("NValidPixHits_Pr"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
      FillHist("NValidStations_Pr"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
      FillHist("NActiveLayer_Pr"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
      FillHist("GlobalChi2_Pr"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 50., 50);
      FillHist("D0_Pr"+Label, muonColl.at(i).dXY(), weight, 0., 1., 1000);
      FillHist("DZ_Pr"+Label, muonColl.at(i).dZ(), weight, 0., 1., 1000); 
      FillHist("RelIso04_Pr"+Label, muonColl.at(i).RelIso04(), weight, 0., 0.5, 50); 
    }
    if(LepType<-1){
      FillHist("IsGlobal_Typem234"+Label, 0., weight, 0., 2., 2);
      if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_Typem234"+Label, 1., weight, 0., 2., 2); 
      FillHist("NValidHits_Typem234"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
      FillHist("NValidPixHits_Typem234"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
      FillHist("NValidStations_Typem234"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
      FillHist("NActiveLayer_Typem234"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
      FillHist("GlobalChi2_Typem234"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 50., 50);
      FillHist("D0_Typem234"+Label, muonColl.at(i).dXY(), weight, 0., 1., 1000);
      FillHist("DZ_Typem234"+Label, muonColl.at(i).dZ(), weight, 0., 1., 1000); 
      FillHist("RelIso04_Typem234"+Label, muonColl.at(i).RelIso04(), weight, 0., 0.5, 50); 
    }
    if(LepType==2){
      FillHist("IsGlobal_Type2"+Label, 0., weight, 0., 2., 2);
      if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_Type2"+Label, 1., weight, 0., 2., 2); 
      FillHist("NValidHits_Type2"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
      FillHist("NValidPixHits_Type2"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
      FillHist("NValidStations_Type2"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
      FillHist("NActiveLayer_Type2"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
      FillHist("GlobalChi2_Type2"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 50., 50);
      FillHist("D0_Type2"+Label, muonColl.at(i).dXY(), weight, 0., 1., 1000);
      FillHist("DZ_Type2"+Label, muonColl.at(i).dZ(), weight, 0., 1., 1000); 
      FillHist("RelIso04_Type2"+Label, muonColl.at(i).RelIso04(), weight, 0., 0.5, 50); 
    }
    if(LepType==-1){
      FillHist("IsGlobal_Typem1"+Label, 0., weight, 0., 2., 2);
      if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_Typem1"+Label, 1., weight, 0., 2., 2); 
      FillHist("NValidHits_Typem1"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
      FillHist("NValidPixHits_Typem1"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
      FillHist("NValidStations_Typem1"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
      FillHist("NActiveLayer_Typem1"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
      FillHist("GlobalChi2_Typem1"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 50., 50);
      FillHist("D0_Typem1"+Label, muonColl.at(i).dXY(), weight, 0., 1., 1000);
      FillHist("DZ_Typem1"+Label, muonColl.at(i).dZ(), weight, 0., 1., 1000); 
      FillHist("RelIso04_Typem1"+Label, muonColl.at(i).RelIso04(), weight, 0., 0.5, 50); 

      if(SrcfPID==211){
        FillHist("IsGlobal_Typem1PID211"+Label, 0., weight, 0., 2., 2);
        if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_Typem1PID211"+Label, 1., weight, 0., 2., 2); 
        FillHist("NValidHits_Typem1PID211"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
        FillHist("NValidPixHits_Typem1PID211"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
        FillHist("NValidStations_Typem1PID211"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
        FillHist("NActiveLayer_Typem1PID211"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
        FillHist("GlobalChi2_Typem1PID211"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 50., 50);
        FillHist("D0_Typem1PID211"+Label, muonColl.at(i).dXY(), weight, 0., 1., 1000);
        FillHist("DZ_Typem1PID211"+Label, muonColl.at(i).dZ(), weight, 0., 1., 1000); 
        FillHist("RelIso04_Typem1PID211"+Label, muonColl.at(i).RelIso04(), weight, 0., 0.5, 50); 
      }
      if(SrcfPID==310){
        FillHist("IsGlobal_Typem1PID310"+Label, 0., weight, 0., 2., 2);
        if(muonColl.at(i).IsGlobal()) FillHist("IsGlobal_Typem1PID310"+Label, 1., weight, 0., 2., 2); 
        FillHist("NValidHits_Typem1PID310"+Label, muonColl.at(i).validHits(), weight, 0., 50., 50);
        FillHist("NValidPixHits_Typem1PID310"+Label, muonColl.at(i).validPixHits(), weight, 0., 20., 20);
        FillHist("NValidStations_Typem1PID310"+Label, muonColl.at(i).validStations(), weight, 0., 20., 20);
        FillHist("NActiveLayer_Typem1PID310"+Label, muonColl.at(i).ActiveLayer(), weight, 0., 20., 20);
        FillHist("GlobalChi2_Typem1PID310"+Label, muonColl.at(i).GlobalChi2(), weight, 0., 50., 50);
        FillHist("D0_Typem1PID310"+Label, muonColl.at(i).dXY(), weight, 0., 1., 1000);
        FillHist("DZ_Typem1PID310"+Label, muonColl.at(i).dZ(), weight, 0., 1., 1000); 
        FillHist("RelIso04_Typem1PID310"+Label, muonColl.at(i).RelIso04(), weight, 0., 0.5, 50); 
      }
    }//End of LepType -1
  }//End of Muon Loop
//Summary
//1)Typem1 fakes have some parts with low pixel hits, low tracker hits, low stations, low valid hits, high chi2.
//  IP is generally large, but around 50-70% smaller than Typem2 fakes; less sensitive to IP, but still responds to IP.
//  1-a) Neutral Kaons have some portions with low pixel hit, low active layer=>decayed in the tracker
//                                             low stations and high chi2(distorted at yoke and hcal)
//  1-b) Pions have prompt like tracker hits(active layer), pixel hits =>kept its momentum in tracker
//                              but low stations, high chi2 =>track distorted in hcal or chamber
//2)At first glance, active global hits and chi2 seems to be possible to retuned before applying all cuts,
//  but after applying all cuts except IP, chi2, global Nhits, global Nhit pattern becomes similar between other fakes(Typem2) and prompts.
//  So global Nhit is not worth putting effort to reoptimise.
//  But chi2 can work. it can be reduced to 3 if one wants to reduce light fakes.
}


void Aug2017_MuFakeMCStudy::InspectFakeRate(std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KJet> JetColl, std::vector<snu::KTruth> truthColl, TString LooseID, TString TightID, float weight, TString Label, TString Option){
  //Input : Very Loose MuColl; Looser than Loose ID, jetNoVetoColl.

  const int NPtEdges=9;
  float PtEdges[NPtEdges] = {0., 10., 15., 20., 25., 35., 50., 70., 100.};
  //float EtaEdges[NEtaEdges]={0., 0.9, 1.6, 2.1, 2.4};
  //float EtaEdges[NEtaEdges]={0., 1.2, 2.1, 2.4};

  

  for(int i=0; i<muonLooseColl.size(); i++){
    if(!PassIDCriteria(muonLooseColl.at(i), LooseID)) continue;
    float PTCorr     = muonLooseColl.at(i).Pt()*(1+RochIso(muonLooseColl.at(i), "0.4"));
    float fEta       = fabs(muonLooseColl.at(i).Eta());
    int   LepType    = GetLeptonType(muonLooseColl.at(i),truthColl);
    int   SrcJetType = GetFakeLepJetSrcType(muonLooseColl.at(i), JetColl);

    if(LepType<-4 || LepType>0) continue;


    FillHist("NMuSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
    if     (LepType==-1)    FillHist("NMuSumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    else if(LepType>-4 )    FillHist("NMuSumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    if     (SrcJetType==3)  FillHist("NMuSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    else if(SrcJetType==2)  FillHist("NMuSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    else if(SrcJetType==1)  FillHist("NMuSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

    if(PassIDCriteria(muonLooseColl.at(i), TightID)){
      FillHist("NMuIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (LepType==-1)    FillHist("NMuIDSumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(LepType>-4 )    FillHist("NMuIDSumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("NMuIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("NMuIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("NMuIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    }

    if(fEta<0.9){
      FillHist("NMuBSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (LepType==-1)    FillHist("NMuBSumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(LepType>-4 )    FillHist("NMuBSumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("NMuBSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("NMuBSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("NMuBSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(PassIDCriteria(muonLooseColl.at(i), TightID)){
        FillHist("NMuBIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (LepType==-1)    FillHist("NMuBIDSumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(LepType>-4 )    FillHist("NMuBIDSumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("NMuBIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuBIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuBIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }
    else if(fEta<1.6){
      FillHist("NMuBESumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (LepType==-1)    FillHist("NMuBESumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(LepType>-4 )    FillHist("NMuBESumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("NMuBESumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("NMuBESumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("NMuBESumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(PassIDCriteria(muonLooseColl.at(i), TightID)){
        FillHist("NMuBEIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (LepType==-1)    FillHist("NMuBEIDSumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(LepType>-4 )    FillHist("NMuBEIDSumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("NMuBEIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuBEIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuBEIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }
    else if(fEta<2.1){
      FillHist("NMuE1SumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (LepType==-1)    FillHist("NMuE1SumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(LepType>-4 )    FillHist("NMuE1SumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("NMuE1SumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("NMuE1SumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("NMuE1SumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(PassIDCriteria(muonLooseColl.at(i), TightID)){
        FillHist("NMuE1IDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (LepType==-1)    FillHist("NMuE1IDSumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(LepType>-4 )    FillHist("NMuE1IDSumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("NMuE1IDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuE1IDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuE1IDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }
    else{
      FillHist("NMuE2SumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (LepType==-1)    FillHist("NMuE2SumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(LepType>-4 )    FillHist("NMuE2SumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("NMuE2SumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("NMuE2SumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("NMuE2SumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(PassIDCriteria(muonLooseColl.at(i), TightID)){
        FillHist("NMuE2IDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (LepType==-1)    FillHist("NMuE2IDSumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(LepType>-4 )    FillHist("NMuE2IDSumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("NMuE2IDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuE2IDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuE2IDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }

  }//End of MuLoop
  
}



void Aug2017_MuFakeMCStudy::ScanFakeRate(std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KTruth> truthColl, TString LooseID, TString TightID, int Nd0SigCuts, float d0SigCuts[], int NChi2Cuts, float Chi2Cuts[], float weight, TString Label, TString Option){
  //Input : Loose ID: ID criteria for loose except d0Sig, chi2 cut.

  const int NPtEdges=9;
  float PtEdges[NPtEdges] = {0., 10., 15., 20., 25., 35., 50., 70., 100.};
  
  float TightIsoCut=0.;
  if     (TightID.Contains("TIsop15")) TightIsoCut = 0.15;
  else if(TightID.Contains("TIsop20")) TightIsoCut = 0.2;
  else if(TightID.Contains("TIsop25")) TightIsoCut = 0.25;

  bool LowHT=false, HighHT=false;
  float HT=0.;
//  std::vector<snu::KJet> JetVetoColl=SkimJetColl(JetColl, EleLColl, muonLooseColl, "EleMuVeto");
//  if(JetVetoColl.size()==0) return;
//  if(JetVetoColl.at(0).Pt()<40) return;
  //for(int i=0; i<JetVetoColl.size(); i++){ HT+=JetVetoColl.at(i).Pt(); }
  //if(HT<50){ LowHT=true;} else{ HighHT=true; }
 

  for(int it_d0Sig=0; it_d0Sig<Nd0SigCuts; it_d0Sig++){
    for(int it_chi2=0; it_chi2<NChi2Cuts; it_chi2++){
      //float d0SigCut=d0Sig_init+d0SigStep*it_d0Sig, Chi2Cut=Chi2_init+Chi2Step*it_chi2;
      float d0SigCut=d0SigCuts[it_d0Sig], Chi2Cut=Chi2Cuts[it_chi2];
      std::ostringstream s1; s1<<d0SigCut;
      TString Str_d0SigCut=s1.str();
      Str_d0SigCut.ReplaceAll("0.","p");   Str_d0SigCut.ReplaceAll(".","p");

      std::ostringstream s2; s2<<Chi2Cut;
      TString Str_Chi2Cut=s2.str();
      Str_Chi2Cut.ReplaceAll("0.","p"); Str_Chi2Cut.ReplaceAll(".","p");


      for(int i=0; i<muonLooseColl.size(); i++){
        if(!PassIDCriteria(muonLooseColl.at(i), LooseID) ) continue;
        if(fabs(muonLooseColl.at(i).dXYSig())>d0SigCut   ) continue;
        if(fabs(muonLooseColl.at(i).GlobalChi2())>Chi2Cut) continue;

        float PTCorr     = ConeCorrectedPT(muonLooseColl.at(i), TightIsoCut);
        float fEta       = fabs(muonLooseColl.at(i).Eta());
        int   LepType    = GetLeptonType(muonLooseColl.at(i),truthColl);
        int   SrcJetType = GetFakeLepJetSrcType(muonLooseColl.at(i), JetColl);
    
        if( LepType<-4 || LepType>0 ) continue;

        //if( !(muonLooseColl.at(i).TriggerMatched("HLT_Mu8_TrkIsoVVL_v")||muonLooseColl.at(i).TriggerMatched("HLT_Mu17_TrkIsoVVL_v")) ) continue;

        if(fEta<0.9){
          FillHist("NMuBSumW_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuBSumW_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuBSumW_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuBSumW_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
          if(PassIDCriteria(muonLooseColl.at(i), TightID)){
            FillHist("NMuBIDSumW_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuBIDSumW_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuBIDSumW_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuBIDSumW_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          }
        }
        else if(fEta<1.6){
          FillHist("NMuBESumW_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuBESumW_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuBESumW_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuBESumW_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
          if(PassIDCriteria(muonLooseColl.at(i), TightID)){
            FillHist("NMuBEIDSumW_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuBEIDSumW_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuBEIDSumW_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuBEIDSumW_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          }
        }
        else if(fEta<2.1){
          FillHist("NMuE1SumW_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuE1SumW_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuE1SumW_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuE1SumW_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
          if(PassIDCriteria(muonLooseColl.at(i), TightID)){
            FillHist("NMuE1IDSumW_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuE1IDSumW_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuE1IDSumW_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuE1IDSumW_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          }
        }
        else{
          FillHist("NMuE2SumW_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuE2SumW_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuE2SumW_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuE2SumW_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
          if(PassIDCriteria(muonLooseColl.at(i), TightID)){
            FillHist("NMuE2IDSumW_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuE2IDSumW_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuE2IDSumW_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuE2IDSumW_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          }
        }


        if(LowHT){
    
          if(fEta<0.9){
            FillHist("NMuBSumW_LowHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuBSumW_LowHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuBSumW_LowHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuBSumW_LowHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
            if(PassIDCriteria(muonLooseColl.at(i), TightID)){
              FillHist("NMuBIDSumW_LowHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
              if     (SrcJetType==3)  FillHist("NMuBIDSumW_LowHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==2)  FillHist("NMuBIDSumW_LowHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==1)  FillHist("NMuBIDSumW_LowHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
          else if(fEta<1.6){
            FillHist("NMuBESumW_LowHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuBESumW_LowHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuBESumW_LowHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuBESumW_LowHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
            if(PassIDCriteria(muonLooseColl.at(i), TightID)){
              FillHist("NMuBEIDSumW_LowHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
              if     (SrcJetType==3)  FillHist("NMuBEIDSumW_LowHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==2)  FillHist("NMuBEIDSumW_LowHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==1)  FillHist("NMuBEIDSumW_LowHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
          else if(fEta<2.1){
            FillHist("NMuE1SumW_LowHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuE1SumW_LowHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuE1SumW_LowHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuE1SumW_LowHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
            if(PassIDCriteria(muonLooseColl.at(i), TightID)){
              FillHist("NMuE1IDSumW_LowHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
              if     (SrcJetType==3)  FillHist("NMuE1IDSumW_LowHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==2)  FillHist("NMuE1IDSumW_LowHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==1)  FillHist("NMuE1IDSumW_LowHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
          else{
            FillHist("NMuE2SumW_LowHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuE2SumW_LowHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuE2SumW_LowHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuE2SumW_LowHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
            if(PassIDCriteria(muonLooseColl.at(i), TightID)){
              FillHist("NMuE2IDSumW_LowHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
              if     (SrcJetType==3)  FillHist("NMuE2IDSumW_LowHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==2)  FillHist("NMuE2IDSumW_LowHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==1)  FillHist("NMuE2IDSumW_LowHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
        }//End of Low HT

        if(HighHT){
    
          if(fEta<0.9){
            FillHist("NMuBSumW_HighHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuBSumW_HighHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuBSumW_HighHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuBSumW_HighHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
            if(PassIDCriteria(muonLooseColl.at(i), TightID)){
              FillHist("NMuBIDSumW_HighHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
              if     (SrcJetType==3)  FillHist("NMuBIDSumW_HighHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==2)  FillHist("NMuBIDSumW_HighHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==1)  FillHist("NMuBIDSumW_HighHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
          else if(fEta<1.6){
            FillHist("NMuBESumW_HighHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuBESumW_HighHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuBESumW_HighHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuBESumW_HighHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
            if(PassIDCriteria(muonLooseColl.at(i), TightID)){
              FillHist("NMuBEIDSumW_HighHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
              if     (SrcJetType==3)  FillHist("NMuBEIDSumW_HighHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==2)  FillHist("NMuBEIDSumW_HighHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==1)  FillHist("NMuBEIDSumW_HighHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
          else if(fEta<2.1){
            FillHist("NMuE1SumW_HighHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuE1SumW_HighHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuE1SumW_HighHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuE1SumW_HighHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
            if(PassIDCriteria(muonLooseColl.at(i), TightID)){
              FillHist("NMuE1IDSumW_HighHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
              if     (SrcJetType==3)  FillHist("NMuE1IDSumW_HighHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==2)  FillHist("NMuE1IDSumW_HighHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==1)  FillHist("NMuE1IDSumW_HighHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
          else{
            FillHist("NMuE2SumW_HighHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
            if     (SrcJetType==3)  FillHist("NMuE2SumW_HighHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==2)  FillHist("NMuE2SumW_HighHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            else if(SrcJetType==1)  FillHist("NMuE2SumW_HighHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      
            if(PassIDCriteria(muonLooseColl.at(i), TightID)){
              FillHist("NMuE2IDSumW_HighHT_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
              if     (SrcJetType==3)  FillHist("NMuE2IDSumW_HighHT_BjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==2)  FillHist("NMuE2IDSumW_HighHT_CjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
              else if(SrcJetType==1)  FillHist("NMuE2IDSumW_HighHT_LjMatch_d0Sig"+Str_d0SigCut+"chi"+Str_Chi2Cut+"_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
        }//End of High HT

      }//End of MuLoop
    }//End of chi2loop
  }//End of d0Sig loop


}


void Aug2017_MuFakeMCStudy::CheckTriggerBias(std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KJet> JetColl, std::vector<snu::KTruth> truthColl, TString LooseID, TString TightID, float weight, TString Label, TString Option){
//Purpose : Check trigger bias. 1) Maximum unbiased extrapolation region : how much can I loosen ID variable?
//                              2) Does trigger requirement change FR even in trigger safe range?

  const int NPtEdges = 9;
  float PtEdges[NPtEdges] = {0., 10., 15., 20., 25., 35., 50., 70., 100.};

  for(int i=0; i<muonLooseColl.size(); i++){
    if(!PassIDCriteria(muonLooseColl.at(i), LooseID)) continue;
    float PTCorr     = muonLooseColl.at(i).Pt()*(1+RochIso(muonLooseColl.at(i), "0.4"));
    float PT         = muonLooseColl.at(i).Pt();
    float fEta       = fabs(muonLooseColl.at(i).Eta());
    int   LepType    = GetLeptonType(muonLooseColl.at(i),truthColl);
    int   SrcJetType = GetFakeLepJetSrcType(muonLooseColl.at(i), JetColl);

    if(LepType<-4 || LepType>0) continue;

    //Trigger Efficiency Cut-Off Check
    FillHist("MuSumW_RelIso04", RochIso(muonLooseColl.at(i), "0.4"), weight, 0., 1., 100);
    if(muonLooseColl.at(i).TriggerMatched("HLT_Mu8_TrkIsoVVL_v")) FillHist("MuTrigSumW_RelIso04", RochIso(muonLooseColl.at(i), "0.4"), weight, 0., 1., 100);


    //Check the FR after passing trigger.
    if(muonLooseColl.at(i).TriggerMatched("HLT_Mu8_TrkIsoVVL_v")){
    
      if(fEta<0.9){
        FillHist("NMuBTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if    (SrcJetType==3 )  FillHist("NMuBTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuBTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuBTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuBIDTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuBIDTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuBIDTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuBIDTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else if(fEta<1.6){
        FillHist("NMuBETrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("NMuBETrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuBETrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuBETrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuBEIDTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuBEIDTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuBEIDTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuBEIDTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else if(fEta<2.1){
        FillHist("NMuE1TrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("NMuE1TrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuE1TrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuE1TrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuE1IDTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuE1IDTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuE1IDTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuE1IDTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else{
        FillHist("NMuE2TrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("NMuE2TrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuE2TrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuE2TrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuE2IDTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuE2IDTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuE2IDTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuE2IDTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
  
    }


    //Selection Division for trigger choices are not giving bias
    if(PTCorr>35 && PT>20 && muonLooseColl.at(i).TriggerMatched("HLT_Mu17_TrkIsoVVL_v")){
      float trigger_ps_reweight=WeightByTrigger("HLT_Mu17_TrkIsoVVL_v", TargetLumi)/WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
      weight*=trigger_ps_reweight;

      if(fEta<0.9){
        FillHist("NMuBMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if    (SrcJetType==3 )  FillHist("NMuBMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuBMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuBMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuBIDMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuBIDMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuBIDMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuBIDMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else if(fEta<1.6){
        FillHist("NMuBEMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("NMuBEMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuBEMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuBEMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuBEIDMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuBEIDMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuBEIDMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuBEIDMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else if(fEta<2.1){
        FillHist("NMuE1MultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("NMuE1MultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuE1MultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuE1MultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuE1IDMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuE1IDMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuE1IDMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuE1IDMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else{
        FillHist("NMuE2MultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("NMuE2MultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuE2MultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuE2MultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuE2IDMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuE2IDMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuE2IDMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuE2IDMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }

    }//PT>20 && PTCorr35-
    else if(PTCorr<35 && PT>10 && muonLooseColl.at(i).TriggerMatched("HLT_Mu8_TrkIsoVVL_v")){
      float trigger_ps_reweight=WeightByTrigger("HLT_Mu8_TrkIsoVVL_v", TargetLumi)/WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
      weight*=trigger_ps_reweight;

      if(fEta<0.9){
        FillHist("NMuBMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if    (SrcJetType==3 )  FillHist("NMuBMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuBMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuBMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuBIDMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuBIDMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuBIDMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuBIDMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else if(fEta<1.6){
        FillHist("NMuBEMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("NMuBEMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuBEMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuBEMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuBEIDMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuBEIDMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuBEIDMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuBEIDMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else if(fEta<2.1){
        FillHist("NMuE1MultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("NMuE1MultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuE1MultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuE1MultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuE1IDMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuE1IDMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuE1IDMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuE1IDMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else{
        FillHist("NMuE2MultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("NMuE2MultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuE2MultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuE2MultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
        if(PassIDCriteria(muonLooseColl.at(i), TightID)){
          FillHist("NMuE2IDMultTrigSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("NMuE2IDMultTrigSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("NMuE2IDMultTrigSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("NMuE2IDMultTrigSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }

    }//End of PT>10 & PTCorr10-35

  }//End of MuLoop

}

void Aug2017_MuFakeMCStudy::EmulateFRMeasurement(std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, std::vector<snu::KTruth> truthColl, TString TightID, float weight, TString Label, TString Option){
  //Input : Very Loose MuColl; Looser than Loose ID, jetNoVetoColl.


  TString AddFilterLabel="";
  if(Option.Contains("TrkIsoVVL")) AddFilterLabel="_TrkIsoVVL";

  const int NPtEdges=9;
  float PtEdges[NPtEdges] = {0., 10., 15., 20., 25., 35., 50., 70., 100.};

  bool PassJetReq=false, PassMETMTWReq=false;
  bool TrigSel=false;
  if( !(MuLColl.size()==1 && EleLColl.size()==0) ) return;

  for(int i=0; i<JetColl.size(); i++){
    if(JetColl.at(i).Pt()<40) continue;
    if(MuLColl.at(0).DeltaR(JetColl.at(i))>0.4) PassJetReq=true;
  }
  if(!PassJetReq) return;

  float MTW = sqrt(2)*sqrt(MET*MuLColl.at(0).Pt()-METx*MuLColl.at(0).Px()-METy*MuLColl.at(0).Py());
  if( (MET<25 && MTW<35) ) PassMETMTWReq=true;
  
  float PT         = MuLColl.at(0).Pt();
  float TightIsoCut=0.;
  if     (TightID.Contains("TIsop15")) TightIsoCut = 0.15;
  else if(TightID.Contains("TIsop20")) TightIsoCut = 0.2;
  else if(TightID.Contains("TIsop25")) TightIsoCut = 0.25;
  //float PTCorr     = MuLColl.at(0).Pt()*(1+RochIso(MuLColl.at(0), "0.4"));
  float PTCorr     = ConeCorrectedPT(MuLColl.at(0), TightIsoCut);
  float fEta       = fabs(MuLColl.at(0).Eta());
  int   LepType    = GetLeptonType(MuLColl.at(0),truthColl);
  int   SrcJetType = GetFakeLepJetSrcType(MuLColl.at(0), JetColl);

  if(LepType<-4 || LepType>0) return;
 

  if     (PassTrigger("HLT_Mu17"+AddFilterLabel+"_v") && PT>20 && PTCorr>35){
    TrigSel=true; weight*=WeightByTrigger("HLT_Mu17"+AddFilterLabel+"_v", TargetLumi)/WeightByTrigger("HLT_IsoMu24_v", TargetLumi);}
  else if(PassTrigger("HLT_Mu8"+AddFilterLabel+"_v")  && PT>10 && PTCorr<35){
    TrigSel=true; weight*=WeightByTrigger("HLT_Mu8"+AddFilterLabel+"_v" , TargetLumi)/WeightByTrigger("HLT_IsoMu24_v", TargetLumi)*1.33;}


  if(TrigSel){
    if(fEta<0.9){
      FillHist("MuBSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (SrcJetType==3)  FillHist("MuBSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("MuBSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("MuBSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
      if(PassIDCriteria(MuLColl.at(0), TightID)){
        FillHist("MuBIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("MuBIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuBIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuBIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }
    else if(fEta<1.6){
      FillHist("MuBESumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (SrcJetType==3)  FillHist("MuBESumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("MuBESumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("MuBESumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
      if(PassIDCriteria(MuLColl.at(0), TightID)){
        FillHist("MuBEIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("MuBEIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuBEIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuBEIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }
    else if(fEta<2.1){
      FillHist("MuE1SumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (SrcJetType==3)  FillHist("MuE1SumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("MuE1SumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("MuE1SumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
      if(PassIDCriteria(MuLColl.at(0), TightID)){
        FillHist("MuE1IDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("MuE1IDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuE1IDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuE1IDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }
    else{
      FillHist("MuE2SumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (SrcJetType==3)  FillHist("MuE2SumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("MuE2SumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("MuE2SumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
      if(PassIDCriteria(MuLColl.at(0), TightID)){
        FillHist("MuE2IDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("MuE2IDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuE2IDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuE2IDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }

    if(PassMETMTWReq){
      if(fEta<0.9){
        FillHist("MuBMETMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("MuBMETMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuBMETMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuBMETMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    
        if(PassIDCriteria(MuLColl.at(0), TightID)){
          FillHist("MuBMETMTWIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("MuBMETMTWIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("MuBMETMTWIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("MuBMETMTWIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else if(fEta<1.6){
        FillHist("MuBEMETMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("MuBEMETMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuBEMETMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuBEMETMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    
        if(PassIDCriteria(MuLColl.at(0), TightID)){
          FillHist("MuBEMETMTWIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("MuBEMETMTWIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("MuBEMETMTWIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("MuBEMETMTWIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else if(fEta<2.1){
        FillHist("MuE1METMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("MuE1METMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuE1METMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuE1METMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    
        if(PassIDCriteria(MuLColl.at(0), TightID)){
          FillHist("MuE1METMTWIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("MuE1METMTWIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("MuE1METMTWIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("MuE1METMTWIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else{
        FillHist("MuE2METMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (SrcJetType==3)  FillHist("MuE2METMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuE2METMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuE2METMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    
        if(PassIDCriteria(MuLColl.at(0), TightID)){
          FillHist("MuE2METMTWIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if     (SrcJetType==3)  FillHist("MuE2METMTWIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("MuE2METMTWIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("MuE2METMTWIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
    }
  }


}




void Aug2017_MuFakeMCStudy::CheckMCClosure(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KTruth> TruthColl, TString LooseID, TString TightID, float weight, TString Label, TString Option){

  TString FilterInfo="", ConeMethod="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";

  std::vector<snu::KMuon> MuLColl, MuTColl;
    for(int i=0; i<MuPreColl.size(); i++){
      if(PassIDCriteria(MuPreColl.at(i), LooseID)) MuLColl.push_back(MuPreColl.at(i));
      if(PassIDCriteria(MuPreColl.at(i), TightID)) MuTColl.push_back(MuPreColl.at(i));
    }
  std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");

 
  //Checking 2PromptMu+1Fk Ele
  float fakeweight = -1.;
  int   NHFakeMu   = NPromptFake_Mu(MuLColl, TruthColl, "HFake");
  int   NHFakeEle  = NPromptFake_Ele(EleLColl, TruthColl, "HFake");
  int   NHFakeLep  = NHFakeMu+NHFakeEle;
  int   NLooseLep  = MuLColl.size()+EleLColl.size();

  //Use only 3 valid leptons 
  bool IsCand = false;
  if( NLooseLep!=3 ) return;
  if( NHFakeLep==0 ) return;//Include at least 1 fake in the sample


  for(int i=0; i<MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), TightID)){
      float FR=FakeRateMC(MuLColl.at(i), "QCD_"+TightID.ReplaceAll("Test_","")+"_"+LooseID.ReplaceAll("Test_","")+"_"+FilterInfo+ConeMethod);
      fakeweight*=-FR/(1-FR);
      //cout<<FR<<endl;
    }
  }
  if( MuLColl.size()==3 && MuTColl.size()==3 ) fakeweight=0;

  bool Pass_Trigger=false;
  if( PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")
     || PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_Trigger=true;
  

  if(Pass_Trigger){
    if(MuLColl.size()==3 && EleLColl.size()==0){
      bool PassSel=true;
      int IdxOS=-1, IdxSS1=-1, IdxSS2=-1; float MOSSS1=-1, MOSSS2=-1;
      if( !(MuLColl.at(0).Pt()>20 && MuLColl.at(1).Pt()>10 && MuLColl.at(2).Pt()>10) ) PassSel=false;
      if( PassSel && fabs(SumCharge(MuLColl))!=1) PassSel=false;
      if( PassSel ){ IdxOS=TriMuChargeIndex(MuLColl,"OS"); IdxSS1=TriMuChargeIndex(MuLColl,"SS1"); IdxSS2=TriMuChargeIndex(MuLColl,"SS2"); 
                     MOSSS1=(MuLColl.at(IdxOS)+MuLColl.at(IdxSS1)).M(), MOSSS2=(MuLColl.at(IdxOS)+MuLColl.at(IdxSS2)).M();}
      if( PassSel && !(MOSSS1>12 && MOSSS2>12) ) PassSel=false;
      if(PassSel){
        bool JetSelPass=JetColl.size()>=3, BJetSelPass=BJetColl.size()>=1;
        FillHist("CutFlow_exp"+Label, 0., weight*fakeweight, 0., 10., 10);
        if(JetSelPass) FillHist("CutFlow_exp"+Label, 1., weight*fakeweight, 0., 10., 10);
        if(JetSelPass && BJetSelPass) FillHist("CutFlow_exp"+Label, 2., weight*fakeweight, 0., 10., 10);

        FillHist("PTmu1_exp"+Label, MuLColl.at(0).Pt(), weight*fakeweight, 0., 200., 40);
        FillHist("PTmu2_exp"+Label, MuLColl.at(1).Pt(), weight*fakeweight, 0., 200., 40);
        FillHist("PTmu3_exp"+Label, MuLColl.at(2).Pt(), weight*fakeweight, 0., 200., 40);
        FillHist("MOSSS1_exp"+Label, MOSSS1, weight*fakeweight, 0., 300., 60);
        FillHist("MOSSS2_exp"+Label, MOSSS2, weight*fakeweight, 0., 200., 40);
        FillHist("Etamu1_exp"+Label, MuLColl.at(0).Eta(), weight*fakeweight, -5., 5., 20);
        FillHist("Etamu2_exp"+Label, MuLColl.at(1).Eta(), weight*fakeweight, -5., 5., 20);
        FillHist("Etamu3_exp"+Label, MuLColl.at(2).Eta(), weight*fakeweight, -5., 5., 20);
        FillHist("Nj_exp"+Label, JetColl.size(), weight*fakeweight, 0., 10., 10);
        FillHist("Nb_exp"+Label, BJetColl.size(), weight*fakeweight, 0., 10., 10);
      }
    }
    if(MuLColl.size()==3 && MuTColl.size()==3 && EleLColl.size()==0){
      bool PassSel=true;
      int IdxOS=-1, IdxSS1=-1, IdxSS2=-1; float MOSSS1=-1, MOSSS2=-1;
      if( !(MuLColl.at(0).Pt()>20 && MuLColl.at(1).Pt()>10 && MuLColl.at(2).Pt()>10) ) PassSel=false;
      if( PassSel && fabs(SumCharge(MuLColl))!=1 ) PassSel=false;
      if( PassSel ){ IdxOS =TriMuChargeIndex(MuLColl,"OS"); IdxSS1=TriMuChargeIndex(MuLColl,"SS1"); IdxSS2=TriMuChargeIndex(MuLColl,"SS2"); 
                     MOSSS1=(MuLColl.at(IdxOS)+MuLColl.at(IdxSS1)).M(), MOSSS2=(MuLColl.at(IdxOS)+MuLColl.at(IdxSS2)).M();}
      if( PassSel && !(MOSSS1>12 && MOSSS2>12) ) PassSel=false;
      if( PassSel ){
        bool JetSelPass=JetColl.size()>=3, BJetSelPass=BJetColl.size()>=1;
        FillHist("CutFlow_obs"+Label, 0., weight, 0., 10., 10);
        if(JetSelPass) FillHist("CutFlow_obs"+Label, 1., weight, 0., 10., 10);
        if(JetSelPass && BJetSelPass) FillHist("CutFlow_obs"+Label, 2., weight, 0., 10., 10);

        FillHist("PTmu1_obs"+Label, MuLColl.at(0).Pt(), weight, 0., 200., 40);
        FillHist("PTmu2_obs"+Label, MuLColl.at(1).Pt(), weight, 0., 200., 40);
        FillHist("PTmu3_obs"+Label, MuLColl.at(2).Pt(), weight, 0., 200., 40);
        FillHist("MOSSS1_obs"+Label, MOSSS1, weight, 0., 300., 60);
        FillHist("MOSSS2_obs"+Label, MOSSS2, weight, 0., 200., 40);
        FillHist("Etamu1_obs"+Label, MuLColl.at(0).Eta(), weight, -5., 5., 20);
        FillHist("Etamu2_obs"+Label, MuLColl.at(1).Eta(), weight, -5., 5., 20);
        FillHist("Etamu3_obs"+Label, MuLColl.at(2).Eta(), weight, -5., 5., 20);
        FillHist("Nj_obs"+Label, JetColl.size(), weight, 0., 10., 10);
        FillHist("Nb_obs"+Label, BJetColl.size(), weight, 0., 10., 10);
      }
    }
  } 

}


void Aug2017_MuFakeMCStudy::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Aug2017_MuFakeMCStudy::BeginCycle() throw( LQError ){

  Message("In begin Cycle", INFO);

  return;

}

Aug2017_MuFakeMCStudy::~Aug2017_MuFakeMCStudy() {
  
  Message("In Aug2017_MuFakeMCStudy Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}




float Aug2017_MuFakeMCStudy::FakeRateMC(snu::KElectron Ele, TString Option){

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


float Aug2017_MuFakeMCStudy::FakeRateMC(snu::KMuon Mu, TString Option){

  float FR=0., TightIsoCut=0.;
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("TIsop15")) TightIsoCut = 0.15;
    else if(Option.Contains("TIsop20")) TightIsoCut = 0.2;
    else if(Option.Contains("TIsop25")) TightIsoCut = 0.25;
  }
  float PTCorr=ConeCorrectedPT(Mu, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Mu.Pt();
  float fEta=fabs(Mu.Eta());

  if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.314222;
      else if(PTCorr<20)  FR=0.222494;
      else if(PTCorr<25)  FR=0.201252;
      else if(PTCorr<35)  FR=0.17477 ;
      else if(PTCorr<50)  FR=0.144774;
      else                FR=0.137049;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.328538;
      else if(PTCorr<20)  FR=0.250277;
      else if(PTCorr<25)  FR=0.239337;
      else if(PTCorr<35)  FR=0.204935;
      else if(PTCorr<50)  FR=0.200698;
      else                FR=0.18832;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.36993 ;
      else if(PTCorr<20)  FR=0.298505;
      else if(PTCorr<25)  FR=0.282521;
      else if(PTCorr<35)  FR=0.271262;
      else if(PTCorr<50)  FR=0.233547;
      else                FR=0.228264;
    }
    else{
      if     (PTCorr<15)  FR=0.368468;
      else if(PTCorr<20)  FR=0.298542;
      else if(PTCorr<25)  FR=0.28087 ;
      else if(PTCorr<35)  FR=0.270855;
      else if(PTCorr<50)  FR=0.248846;
      else                FR=0.24872;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGLIsop6IPp5p1sig4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.370499;
      else if(PTCorr<20)  FR=0.192715;
      else if(PTCorr<25)  FR=0.162086;
      else if(PTCorr<35)  FR=0.138789;
      else                FR=0.109161;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.396116;
      else if(PTCorr<20)  FR=0.228319;
      else if(PTCorr<25)  FR=0.210723;
      else if(PTCorr<35)  FR=0.173594;
      else                FR=0.16023 ;
    }
    else if(fEta<2.1){
      if     (PTCorr<15)  FR=0.437458;
      else if(PTCorr<20)  FR=0.279085;
      else if(PTCorr<25)  FR=0.260136;
      else if(PTCorr<35)  FR=0.249217;
      else                FR=0.203574;
    }
    else{
      if     (PTCorr<15)  FR=0.410967;
      else if(PTCorr<20)  FR=0.268164;
      else if(PTCorr<25)  FR=0.25555 ;
      else if(PTCorr<35)  FR=0.238569;
      else                FR=0.209493;
    }
  }


  return FR;
}



int Aug2017_MuFakeMCStudy::GetFakeLepSrcIdx(snu::KMuon Mu, std::vector<snu::KTruth> TruthColl){

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



int Aug2017_MuFakeMCStudy::GetFakeLepJetSrcType(snu::KMuon Mu, std::vector<snu::KJet> JetColl){
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


float Aug2017_MuFakeMCStudy::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  //float PTCorr=Ele.Pt()*(1+min(0,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}


float Aug2017_MuFakeMCStudy::ConeCorrectedPT(snu::KMuon Mu, double TightIsoCut){

  //float PTCorr=Mu.Pt();
  //float PTCorr=Mu.Pt()*(1+RochIso(Mu,"0.4"));
  float PTCorr=Mu.Pt()*(1.+max(0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}


int Aug2017_MuFakeMCStudy::StepPassed(std::vector<snu::KMuon> MuColl, std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, TString Option){

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


int Aug2017_MuFakeMCStudy::NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false;
   if(Option.Contains("EWPrompt")) ReturnPrompt=true;
   if(Option.Contains("HFake"))    ReturnFake=true;

   for(int i=0; i<MuColl.size(); i++){
     int MuType=GetLeptonType(MuColl.at(i), truthColl);

     if     (MuType>0 && MuType<4) Nprompt++;
     else if(MuType<0 && MuType>-5 ) Nfake++;
   }

   if(ReturnPrompt && ReturnFake) return Nprompt+Nfake;
   else if(ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;
  
}

int Aug2017_MuFakeMCStudy::NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> truthColl, TString Option){

   int Nprompt=0, Nfake=0;

   bool ReturnPrompt=false, ReturnFake=false, ReturnConv=false;
   if(Option.Contains("EWPrompt")) ReturnPrompt=true;
   if(Option.Contains("HFake"))    ReturnFake=true;
   if(Option.Contains("Conv"))     ReturnConv=true;

   for(int i=0; i<EleColl.size(); i++){
     int EleType=GetLeptonType(EleColl.at(i), truthColl);

     if     (EleType>0 && EleType<4)  Nprompt++;
     else if(EleType<0 && EleType>-5) Nfake++;
   }

   if(ReturnPrompt && ReturnFake) return Nprompt+Nfake;
   else if(ReturnPrompt) return Nprompt;
   else if(ReturnFake)   return Nfake;

   return 0.;

}



bool Aug2017_MuFakeMCStudy::IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option){

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



void Aug2017_MuFakeMCStudy::FillCutFlow(TString cut, float weight){
  
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



void Aug2017_MuFakeMCStudy::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Aug2017_MuFakeMCStudy::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);



  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Aug2017_MuFakeMCStudyCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Aug2017_MuFakeMCStudy::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
