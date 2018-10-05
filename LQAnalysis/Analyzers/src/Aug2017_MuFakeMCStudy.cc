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
   bool AltCRTest=false;
   bool ConvenerAsk  =false;
   for(int i=0; i<(int) k_flags.size(); i++){
     if     (k_flags.at(i).Contains("FakeCompCheck"))    FakeCompCheck    = true;
     else if(k_flags.at(i).Contains("IDEffCheck"))       IDEffCheck       = true;
     else if(k_flags.at(i).Contains("IDVarSensitivity")) IDVarSensitivity = true;
     else if(k_flags.at(i).Contains("FRInspection"))     FRInspection     = true;
     else if(k_flags.at(i).Contains("SiglWP"))           SiglWP           = true;
     else if(k_flags.at(i).Contains("ScanFR"))           ScanFR           = true;
     else if(k_flags.at(i).Contains("TrigBiasCheck"))    TrigBiasCheck    = true;
     else if(k_flags.at(i).Contains("FRMeasEmul"))       FRMeasEmul       = true;
     else if(k_flags.at(i).Contains("Closure"))          Closure          = true;
     else if(k_flags.at(i).Contains("ConvenerAsk"))      ConvenerAsk      = true;
     else if(k_flags.at(i).Contains("AltCRTest"))        AltCRTest        = true;
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
   if(Closure){ if( !(muonPreColl.size()>=3) ) return; }
   /**********************************************************************************************************/

   //For Fake Study
   //Muon ID's to Test
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
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronLooseColl; eventbase->GetElectronSel()->Selection(electronLooseColl);
     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     //eventbase->GetJetSel()->SetPileUpJetID(true,"Loose");
     bool LeptonVeto=false;
   std::vector<snu::KJet> jetNoVetoColl; eventbase->GetJetSel()->Selection(jetNoVetoColl, LeptonVeto, muonPOGLIsoVLColl, electronLooseColl);
   std::vector<snu::KJet> bjetNoVetoColl = SelBJets(jetNoVetoColl, "Medium");

   std::vector<snu::KJet> jetColl  = SkimJetColl(jetNoVetoColl,  electronLooseColl, muonPOGLIsoVLColl, "EleMuVeto");
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

   if(false){
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

     //CheckFakeSources(muonPreColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_NoIDPFMu", "");
     //CheckFakeSources(muonPOGLIsoVLColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_POGLIso", "");
     //CheckFakeSources(muonPOGMColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_POGM", "");
     //CheckFakeSources(muonPOGTColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_POGT", "");
     //CheckFakeSources(muonPOGTIsoColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_POGTIso", "");
     std::vector<snu::KMuon> muonHctoWATColl;
     for(int i=0; i<(int) muonPreColl.size(); i++){
       if(PassIDCriteria(muonPreColl.at(i),"POGTIsop20IPp01p05sig4Chi4","Roch")) muonHctoWATColl.push_back(muonPreColl.at(i));
     }
     bool Mgt12=true;
     for(int i=0; i<(int) muonHctoWATColl.size(); i++){
       for(int j=i+1; j<(int) muonHctoWATColl.size(); j++){
         if( (muonHctoWATColl.at(i)+muonHctoWATColl.at(j)).M()<12 ){ Mgt12=false; break; }
       }
     }
     if(Mgt12 && muonHctoWATColl.size()>2 && muonHctoWATColl.at(0).Pt()>20){
       CheckFakeSources(muonHctoWATColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_HctoWAPreSel", "");
       if(bjetNoVetoColl.size()>0){
         CheckFakeSources(muonHctoWATColl, muonPreColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl, weight, "_HctoWAPreSelHasB", "");
       }
     }
   }
   if(IDEffCheck){
     //Purpose : Check overall ID performance

     CheckIDEfficiency(muonPreColl, truthColl, weight, "", "");

   }
   if(IDVarSensitivity){
     //Purpose : Which ID variable is effective to which kinds of fakes.

     std::vector<snu::KMuon> muonPOGLIsop6IPpXp1;
     for(int i=0; i<(int) muonPreColl.size(); i++){
       if(PassIDCriteria(muonPreColl.at(i),"Test_POGLIsop6IPpXp1","Roch")) muonPOGLIsop6IPpXp1.push_back(muonPreColl.at(i));
     }

     CheckIDVarSensitivity(muonPOGLIsop6IPpXp1, jetNoVetoColl, truthColl, weight, "", "");

     std::vector<snu::KMuon> muon10to20Coll, muon20to50Coll, muon50to100Coll;
     for(int i=0; i<(int) muonPOGLIsop6IPpXp1.size(); i++){
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

     std::vector<snu::KMuon> muonTrigPreColl;
     for(int i=0; i<(int) muonPreColl.size(); i++){
       if( muonPreColl.at(i).TriggerMatched("HLT_Mu8_TrkIsoVVL_v") ) muonTrigPreColl.push_back(muonPreColl.at(i));
     }

     if(SiglWP){
       InspectFakeRate(muonPreColl, jetNoVetoColl, truthColl,
         "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1sig4", "");

       InspectFakeRate(muonTrigPreColl, jetNoVetoColl, truthColl,
         "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1sig4Trig", "");
     }
     if(ScanFR){
       //const int Nd0Cuts  =12; float d0Cuts[Nd0Cuts]    ={0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.1, 0.2, 0.5, 1., 10., 100.};
       //const int Nd0Cuts   =1;  float d0Cuts[Nd0Cuts]    ={0.5};
       const int NChi2Cuts =8;  float Chi2Cuts[NChi2Cuts]  ={4., 10., 20., 30., 50., 80., 100., 999.};
       const int Nd0SigCuts=11; float d0SigCuts[Nd0SigCuts]={4., 5., 7., 8., 9., 10., 20., 30., 50., 100., 999.};

       ScanFakeRate(muonPreColl, electronLooseColl, jetNoVetoColl, truthColl, "POGTIsop4IPp2p1NoChi", "POGTIsop20IPp01p05sig4Chi4", Nd0SigCuts, d0SigCuts, NChi2Cuts, Chi2Cuts, weight, "_TIsop20IPp01p05sig4Chi4LIsop4", "");
       ScanFakeRate(muonPreColl, electronLooseColl, jetNoVetoColl, truthColl, "POGTIsop6IPp2p1NoChi", "POGTIsop20IPp01p05sig4Chi4", Nd0SigCuts, d0SigCuts, NChi2Cuts, Chi2Cuts, weight, "_TIsop20IPp01p05sig4Chi4LIsop6", "");
       ScanFakeRate(muonTrigPreColl, electronLooseColl, jetNoVetoColl, truthColl, "POGTIsop4IPp2p1NoChi", "POGTIsop20IPp01p05sig4Chi4", Nd0SigCuts, d0SigCuts, NChi2Cuts, Chi2Cuts, weight, "_TIsop20IPp01p05sig4Chi4LIsop4Trig", "");
       ScanFakeRate(muonTrigPreColl, electronLooseColl, jetNoVetoColl, truthColl, "POGTIsop6IPp2p1NoChi", "POGTIsop20IPp01p05sig4Chi4", Nd0SigCuts, d0SigCuts, NChi2Cuts, Chi2Cuts, weight, "_TIsop20IPp01p05sig4Chi4LIsop6Trig", "");
     }

   }
   if(TrigBiasCheck){
     //Purpose : Check trigger bias. 1) Maximum unbiased extrapolation region : how much can I loosen ID variable?
     //                              2) Does trigger requirement change FR even in trigger safe range?

     std::vector<snu::KMuon> muonTrigDiPreColl, muonTrigSiglPreColl;
     for(int i=0; i<(int) muonPreColl.size(); i++){
       if( muonPreColl.at(i).TriggerMatched("HLT_Mu8_TrkIsoVVL_v") ) muonTrigDiPreColl.push_back(muonPreColl.at(i));
       if(    muonPreColl.at(i).TriggerMatched("HLT_IsoMu24_v")
           or muonPreColl.at(i).TriggerMatched("HLT_IsoTkMu24_v")  ) muonTrigSiglPreColl.push_back(muonPreColl.at(i));
     }

     InspectFakeRate(muonPreColl, jetNoVetoColl, truthColl,
       "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1sig4", "");

     InspectFakeRate(muonTrigDiPreColl, jetNoVetoColl, truthColl,
       "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1sig4TrigDi", "");
     InspectFakeRate(muonTrigSiglPreColl, jetNoVetoColl, truthColl,
       "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1sig4TrigSigl", "");


     //CheckTriggerBias(muonPreColl, jetNoVetoColl, truthColl, "Test_POGLIsop5IPp5p1Chi", "Test_POGTIsop20IPp01p1Chi3", weight, "", "");
     //CheckTriggerBias(muonPreColl, jetNoVetoColl, truthColl, "POGLIsop4IPp5p1Chi100", "POGTIsop20IPp01p05sig4Chi4", weight, "", "");
   }
   if(FRMeasEmul){
     //Purpose : Fully emulate FR measurement for QCD. Result will be used for closure to DY, TT

//     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
//       "POGLIsop4IPp5p1Chi100", "POGTIsop20IPp01p05sig4Chi4", weight, "_TIsop20IPp01p05sig4Chi4LIsop4IPp5p1Chi100", "TrkIsoVVL");
//     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
//       "POGTIsop4IPp2p1NoChi", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop4IPp2p1NoChi", "TrkIsoVVL");
//     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
//       "POGTIsop6IPp2p1NoChi", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1NoChi", "TrkIsoVVL");
//     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
//       "POGTIsop4IPp01p05sig4Chi4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4LIsop4", "TrkIsoVVL");
//     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
//       "POGTIsop6IPp01p05sig4Chi4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4LIsop6", "TrkIsoVVL");
//     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
//       "POGTIsop4IPp2p1sig4NoChi", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop4IPp2p1sig4NoChi", "TrkIsoVVL");
//     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
//       "POGTIsop6IPp2p1sig4NoChi", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1sig4NoChi", "TrkIsoVVL");
//     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
//       "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1sig4", "TrkIsoVVL");

     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
       "POGLIsop4IPp5p1Chi100", "POGTIsop20IPp01p05sig4Chi4", weight, "_TIsop20IPp01p05sig4Chi4LIsop4IPp5p1Chi100", "");
     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
       "POGTIsop4IPp2p1", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop4IPp2p1", "");
     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
       "POGTIsop6IPp2p1", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1", "");
     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
       "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4POGTIsop6IPp2p1sig4", "");
     EmulateFRMeasurement(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
       "POGTIsop6IPp01p05sig4Chi4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4LIsop6", "");

   }
   if(Closure){
     //Purpose : Description of yields in different sample expected from QCD measured FR
     CheckMCClosure(muonPreColl, electronLooseColl, jetNoVetoColl, truthColl,
       "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", weight, "_TIsop20IPp01p05sig4Chi4LIsop6IPp2p1sig4", "TrkIsoVVLConeSUSY"); 
     CheckMCClosure(muonPreColl, electronLooseColl, jetNoVetoColl, truthColl,
       "POGTIsop6IPp01p05sig4Chi4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4LIsop6", "TrkIsoVVLConeSUSY");   

     CheckMCClosure(muonPreColl, electronLooseColl, jetNoVetoColl, truthColl,
       "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", weight, "_TIsop20IPp01p05sig4Chi4LIsop6IPp2p1sig4_NoFilter", "NoFilterConeSUSY"); 
     CheckMCClosure(muonPreColl, electronLooseColl, jetNoVetoColl, truthColl,
       "POGTIsop6IPp01p05sig4Chi4", "POGTIsop20IPp01p05sig4Chi4", weight, "_POGTIsop20IPp01p05sig4Chi4LIsop6_NoFilter", "NoFilterConeSUSY");   

     //#"TrkIsoVVLFOPt" #"TrkIsoVVLConeSUSY" #"TrkIsoVVLConeE" #"NoFilterFOPt" #"NoFilterConeSUSY" #"NoFilterConeE"
   }
   if(ConvenerAsk){
       eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
       eventbase->GetJetSel()->SetPt(10.);                     eventbase->GetJetSel()->SetEta(2.4);
       bool LeptonVeto=false;
     std::vector<snu::KJet> jetNoVeto10Coll; eventbase->GetJetSel()->Selection(jetNoVeto10Coll, LeptonVeto, muonPOGLIsoVLColl, electronLooseColl);
     std::vector<snu::KJet> bjetNoVeto10Coll = SelBJets(jetNoVeto10Coll, "Medium");

     std::vector<snu::KGenJet> genjetColl;  eventbase->GetGenJetSel()->Selection(genjetColl);//PT>8 at ntuple level skim

     std::vector<snu::KMuon> muonTightNoIPIsoColl;
     for(int i=0; i<muonPreColl.size(); i++){ if(muonPreColl.at(i).IsTight() && muonPreColl.at(i).GlobalChi2()<4.) muonTightNoIPIsoColl.push_back(muonPreColl.at(i)); }

     CheckMotherDaugherRelationship(muonTightNoIPIsoColl, muonPreColl, electronLooseColl, jetNoVeto10Coll, bjetNoVeto10Coll, met, met*cos(metphi), met*sin(metphi), truthColl, genjetColl,
       weight, "", "");
   } 
   if(AltCRTest){

     CheckAltCRAvailability(muonPreColl, electronLooseColl, jetNoVetoColl, met, met*cos(metphi), met*sin(metphi), truthColl,
       "POGTIsop6IPp2p1sig4", "POGTIsop20IPp01p05sig4Chi4", weight, "", ""); 
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

  for(int i=0; i<(int) muonColl.size(); i++){
    int LepType    = GetLeptonType(muonColl.at(i),truthColl);
    int SrcJetType = GetFakeLepJetSrcType(muonColl.at(i), JetColl);
    int SrcIdx     = GetFakeLepSrcIdx(muonColl.at(i), truthColl);
    int SrcfPID    = SrcIdx==-1 ? 0 : fabs(truthColl.at(SrcIdx).PdgId());

    if( LepType<-4 || LepType>-1 ) continue;

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


  for(int i=0; i<(int) PromptColl.size(); i++){
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
  for(int i=0; i<(int) FakeColl.size(); i++){
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

  for(int i=0; i<(int) muonColl.size(); i++){
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

  float TightIsoCut=0.;
  if     (TightID.Contains("TIsop15")) TightIsoCut = 0.15;
  else if(TightID.Contains("TIsop20")) TightIsoCut = 0.2;
  else if(TightID.Contains("TIsop25")) TightIsoCut = 0.25;

  

  for(int i=0; i<(int) muonLooseColl.size(); i++){
    if(!PassIDCriteria(muonLooseColl.at(i), LooseID)) continue;
    float PTCorr     = ConeCorrectedPT(muonLooseColl.at(i), TightIsoCut);
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
    else{
      FillHist("NMuESumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (LepType==-1)    FillHist("NMuESumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(LepType>-4 )    FillHist("NMuESumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("NMuESumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("NMuESumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("NMuESumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);

      if(PassIDCriteria(muonLooseColl.at(i), TightID)){
        FillHist("NMuEIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (LepType==-1)    FillHist("NMuEIDSumW_Typem1_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(LepType>-4 )    FillHist("NMuEIDSumW_Typem23_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("NMuEIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("NMuEIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("NMuEIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
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


      for(int i=0; i<(int) muonLooseColl.size(); i++){
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

  for(int i=0; i<(int) muonLooseColl.size(); i++){
    float PTCorr     = muonLooseColl.at(i).Pt()*(1+RochIso(muonLooseColl.at(i), "0.4"));
    float PT         = muonLooseColl.at(i).Pt();
    float fEta       = fabs(muonLooseColl.at(i).Eta());
    int   LepType    = GetLeptonType(muonLooseColl.at(i),truthColl);
    int   SrcJetType = GetFakeLepJetSrcType(muonLooseColl.at(i), JetColl);

    if(LepType<-4 || LepType>0) continue;

    //Trigger Efficiency Cut-Off Check
    int IdxID=0.;
     if(muonLooseColl.at(i).IsLoose()) IdxID++;
     if(muonLooseColl.at(i).IsTight()) IdxID++;
    FillHist("MuSumW_RelIso04"      +Label, RochIso(muonLooseColl.at(i), "0.4")   , weight, 0., 1.  , 100);
    FillHist("MuSumW_d0"            +Label, fabs(muonLooseColl.at(i).dXY())       , weight, 0., 1.  , 1000);
    FillHist("MuSumW_dz"            +Label, fabs(muonLooseColl.at(i).dZ())        , weight, 0., 1.  , 1000);
    FillHist("MuSumW_Chi2"          +Label, fabs(muonLooseColl.at(i).GlobalChi2()), weight, 0., 100., 1000);
    FillHist("MuSumW_NValidHits"    +Label, muonLooseColl.at(i).validHits()       , weight, 0., 50. , 50);
    FillHist("MuSumW_NValidPixHits" +Label, muonLooseColl.at(i).validPixHits()    , weight, 0., 20. , 20);
    FillHist("MuSumW_NValidStations"+Label, muonLooseColl.at(i).validStations()   , weight, 0., 20. , 20);
    FillHist("MuSumW_NActiveLayer"  +Label, muonLooseColl.at(i).ActiveLayer()     , weight, 0., 20. , 20);
    FillHist("MuSumW_ID"            +Label, IdxID                                 , weight, 0., 20. , 20);
    if(muonLooseColl.at(i).TriggerMatched("HLT_Mu8_TrkIsoVVL_v")){
      FillHist("MuTrigSumW_RelIso04"      +Label, RochIso(muonLooseColl.at(i), "0.4")   , weight, 0., 1.  , 100);
      FillHist("MuTrigSumW_d0"            +Label, fabs(muonLooseColl.at(i).dXY())       , weight, 0., 1.  , 1000);
      FillHist("MuTrigSumW_dz"            +Label, fabs(muonLooseColl.at(i).dZ())        , weight, 0., 1.  , 1000);
      FillHist("MuTrigSumW_Chi2"          +Label, fabs(muonLooseColl.at(i).GlobalChi2()), weight, 0., 100., 1000);
      FillHist("MuTrigSumW_NValidHits"    +Label, muonLooseColl.at(i).validHits()       , weight, 0., 50. , 50);
      FillHist("MuTrigSumW_NValidPixHits" +Label, muonLooseColl.at(i).validPixHits()    , weight, 0., 20. , 20);
      FillHist("MuTrigSumW_NValidStations"+Label, muonLooseColl.at(i).validStations()   , weight, 0., 20. , 20);
      FillHist("MuTrigSumW_NActiveLayer"  +Label, muonLooseColl.at(i).ActiveLayer()     , weight, 0., 20. , 20);
      FillHist("MuTrigSumW_ID"            +Label, IdxID                                 , weight, 0., 20. , 20);
    }
    if(LepType==-1){
      FillHist("MuTypem1SumW_RelIso04"      +Label, RochIso(muonLooseColl.at(i), "0.4")   , weight, 0., 1.  , 100);
      FillHist("MuTypem1SumW_d0"            +Label, fabs(muonLooseColl.at(i).dXY())       , weight, 0., 1.  , 1000);
      FillHist("MuTypem1SumW_dz"            +Label, fabs(muonLooseColl.at(i).dZ())        , weight, 0., 1.  , 1000);
      FillHist("MuTypem1SumW_Chi2"          +Label, fabs(muonLooseColl.at(i).GlobalChi2()), weight, 0., 100., 1000);
      FillHist("MuTypem1SumW_NValidHits"    +Label, muonLooseColl.at(i).validHits()       , weight, 0., 50. , 50);
      FillHist("MuTypem1SumW_NValidPixHits" +Label, muonLooseColl.at(i).validPixHits()    , weight, 0., 20. , 20);
      FillHist("MuTypem1SumW_NValidStations"+Label, muonLooseColl.at(i).validStations()   , weight, 0., 20. , 20);
      FillHist("MuTypem1SumW_NActiveLayer"  +Label, muonLooseColl.at(i).ActiveLayer()     , weight, 0., 20. , 20);
      FillHist("MuTypem1SumW_ID"            +Label, IdxID                                 , weight, 0., 20. , 20);
      if(muonLooseColl.at(i).TriggerMatched("HLT_Mu8_TrkIsoVVL_v")){
        FillHist("MuTypem1TrigSumW_RelIso04"      +Label, RochIso(muonLooseColl.at(i), "0.4")   , weight, 0., 1.  , 100);
        FillHist("MuTypem1TrigSumW_d0"            +Label, fabs(muonLooseColl.at(i).dXY())       , weight, 0., 1.  , 1000);
        FillHist("MuTypem1TrigSumW_dz"            +Label, fabs(muonLooseColl.at(i).dZ())        , weight, 0., 1.  , 1000);
        FillHist("MuTypem1TrigSumW_Chi2"          +Label, fabs(muonLooseColl.at(i).GlobalChi2()), weight, 0., 100., 1000);
        FillHist("MuTypem1TrigSumW_NValidHits"    +Label, muonLooseColl.at(i).validHits()       , weight, 0., 50. , 50);
        FillHist("MuTypem1TrigSumW_NValidPixHits" +Label, muonLooseColl.at(i).validPixHits()    , weight, 0., 20. , 20);
        FillHist("MuTypem1TrigSumW_NValidStations"+Label, muonLooseColl.at(i).validStations()   , weight, 0., 20. , 20);
        FillHist("MuTypem1TrigSumW_NActiveLayer"  +Label, muonLooseColl.at(i).ActiveLayer()     , weight, 0., 20. , 20);
        FillHist("MuTypem1TrigSumW_ID"            +Label, IdxID                                 , weight, 0., 20. , 20);
      }
    }
    else if(LepType>=-4 && LepType<-1){
      FillHist("MuTypem234SumW_RelIso04"      +Label, RochIso(muonLooseColl.at(i), "0.4")   , weight, 0., 1.  , 100);
      FillHist("MuTypem234SumW_d0"            +Label, fabs(muonLooseColl.at(i).dXY())       , weight, 0., 1.  , 1000);
      FillHist("MuTypem234SumW_dz"            +Label, fabs(muonLooseColl.at(i).dZ())        , weight, 0., 1.  , 1000);
      FillHist("MuTypem234SumW_Chi2"          +Label, fabs(muonLooseColl.at(i).GlobalChi2()), weight, 0., 100., 1000);
      FillHist("MuTypem234SumW_NValidHits"    +Label, muonLooseColl.at(i).validHits()       , weight, 0., 50. , 50);
      FillHist("MuTypem234SumW_NValidPixHits" +Label, muonLooseColl.at(i).validPixHits()    , weight, 0., 20. , 20);
      FillHist("MuTypem234SumW_NValidStations"+Label, muonLooseColl.at(i).validStations()   , weight, 0., 20. , 20);
      FillHist("MuTypem234SumW_NActiveLayer"  +Label, muonLooseColl.at(i).ActiveLayer()     , weight, 0., 20. , 20);
      FillHist("MuTypem234SumW_ID"            +Label, IdxID                                 , weight, 0., 20. , 20);
      if(muonLooseColl.at(i).TriggerMatched("HLT_Mu8_TrkIsoVVL_v")){
        FillHist("MuTypem234TrigSumW_RelIso04"      +Label, RochIso(muonLooseColl.at(i), "0.4")   , weight, 0., 1.  , 100);
        FillHist("MuTypem234TrigSumW_d0"            +Label, fabs(muonLooseColl.at(i).dXY())       , weight, 0., 1.  , 1000);
        FillHist("MuTypem234TrigSumW_dz"            +Label, fabs(muonLooseColl.at(i).dZ())        , weight, 0., 1.  , 1000);
        FillHist("MuTypem234TrigSumW_Chi2"          +Label, fabs(muonLooseColl.at(i).GlobalChi2()), weight, 0., 100., 1000);
        FillHist("MuTypem234TrigSumW_NValidHits"    +Label, muonLooseColl.at(i).validHits()       , weight, 0., 50. , 50);
        FillHist("MuTypem234TrigSumW_NValidPixHits" +Label, muonLooseColl.at(i).validPixHits()    , weight, 0., 20. , 20);
        FillHist("MuTypem234TrigSumW_NValidStations"+Label, muonLooseColl.at(i).validStations()   , weight, 0., 20. , 20);
        FillHist("MuTypem234TrigSumW_NActiveLayer"  +Label, muonLooseColl.at(i).ActiveLayer()     , weight, 0., 20. , 20);
        FillHist("MuTypem234TrigSumW_ID"            +Label, IdxID                                 , weight, 0., 20. , 20);
      }
    }


    //For Here only loose ID muon used.
    if(!PassIDCriteria(muonLooseColl.at(i), LooseID)) continue;

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

void Aug2017_MuFakeMCStudy::EmulateFRMeasurement(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, float MET, float METx, float METy, std::vector<snu::KTruth> truthColl, TString LooseID, TString TightID, float weight, TString Label, TString Option){
  //Input : Very Loose MuColl; Looser than Loose ID, jetNoVetoColl.


  TString AddFilterLabel="";
  if(Option.Contains("TrkIsoVVL")) AddFilterLabel="_TrkIsoVVL";

  std::vector<snu::KMuon> MuLColl;
   for(int i=0; i<(int) MuPreColl.size(); i++){
     if(PassIDCriteria(MuPreColl.at(i),LooseID,"Roch")) MuLColl.push_back(MuPreColl.at(i));
   }


  std::vector<snu::KJet> JetVetoColl    = SkimJetColl(JetNoVetoColl, EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetVetoColl   = SelBJets(JetVetoColl, "Medium");
  std::vector<snu::KJet> BJetNoVetoColl = SelBJets(JetNoVetoColl, "Medium");
  //std::vector<snu::KJet> BJetVetoMColl = SelBJets(JetVetoColl, "Medium");
  //std::vector<snu::KJet> BJetMColl     = SelBJets(JetNoVetoColl, "Medium");
  //std::vector<snu::KJet> BJetVetoTColl = SelBJets(JetVetoColl, "Tight");
  //std::vector<snu::KJet> BJetTColl     = SelBJets(JetNoVetoColl, "Tight");
  //std::vector<snu::KJet> BJetVetoLColl = SelBJets(JetVetoColl, "Loose");
  //std::vector<snu::KJet> BJetLColl     = SelBJets(JetNoVetoColl, "Loose");

//  int NBJets=100;
//  if(Option.Contains("Veto")){
//    if(Option.Contains("Loose"))  NBJets=BJetVetoLColl.size();
//    if(Option.Contains("Medium")) NBJets=BJetVetoMColl.size();
//    if(Option.Contains("Tight"))  NBJets=BJetVetoTColl.size();
//  }
//  else{
//    if(Option.Contains("Loose"))  NBJets=BJetLColl.size();
//    if(Option.Contains("Medium")) NBJets=BJetMColl.size();
//    if(Option.Contains("Tight"))  NBJets=BJetTColl.size();
//  }
  int NBJets = BJetVetoColl.size();
  bool IsNearB=false;

  const int NPtEdges=9;
  float PtEdges[NPtEdges] = {0., 10., 15., 20., 25., 35., 50., 70., 100.};

  bool PassJetReq=false, PassMETMTWReq=false;
  bool TrigSel=false;
  if( !(MuLColl.size()==1 && EleLColl.size()==0) ) return;

  for(int i=0; i<(int) JetNoVetoColl.size(); i++){
    if(JetNoVetoColl.at(i).Pt()<40) continue;
    if(MuLColl.at(0).DeltaR(JetNoVetoColl.at(i))>0.4) PassJetReq=true;
  }
  if(!PassJetReq) return;

   for(int i=0; i<(int) BJetNoVetoColl.size(); i++){
     if(BJetNoVetoColl.at(i).DeltaR(MuLColl.at(0))<0.4) IsNearB=true;
   }


  float MTW = sqrt(2)*sqrt(MET*MuLColl.at(0).Pt()-METx*MuLColl.at(0).Px()-METy*MuLColl.at(0).Py());
  if( (MET<25 && MTW<25) ) PassMETMTWReq=true;
  
  float PT         = MuLColl.at(0).Pt();
  float TightIsoCut=0.;
  if     (TightID.Contains("TIsop15")) TightIsoCut = 0.15;
  else if(TightID.Contains("TIsop20")) TightIsoCut = 0.2;
  else if(TightID.Contains("TIsop25")) TightIsoCut = 0.25;

  float PTCorr     = ConeCorrectedPT(MuLColl.at(0), TightIsoCut);
  float fEta       = fabs(MuLColl.at(0).Eta());
  int   LepType    = GetLeptonType(MuLColl.at(0),truthColl);
  int   SrcJetType = GetFakeLepJetSrcType(MuLColl.at(0), JetNoVetoColl);

  if(LepType<-4 || LepType>0) return;
 

  if     (PassTrigger("HLT_Mu17"+AddFilterLabel+"_v") && PT>20 && PTCorr>35){
    TrigSel=true; weight*=WeightByTrigger("HLT_Mu17"+AddFilterLabel+"_v", TargetLumi)/WeightByTrigger("HLT_IsoMu24_v", TargetLumi);}
  else if(PassTrigger("HLT_Mu8"+AddFilterLabel+"_v")  && PT>10 && PTCorr<35){
    TrigSel=true; weight*=WeightByTrigger("HLT_Mu8"+AddFilterLabel+"_v" , TargetLumi)/WeightByTrigger("HLT_IsoMu24_v", TargetLumi)*1.33;}


  if(TrigSel){
    //For Composition Check==============================================//
    FillHist("MuSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
    if     (SrcJetType==3)  FillHist("MuSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    else if(SrcJetType==2)  FillHist("MuSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    else if(SrcJetType==1)  FillHist("MuSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    if(NBJets>0){
      FillHist("MuSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (IsNearB      )  FillHist("MuSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else                    FillHist("MuSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("MuSumW_BjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("MuSumW_CjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("MuSumW_LjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    }
    else{
      FillHist("MuSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("MuSumW_BjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("MuSumW_CjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("MuSumW_LjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    }
    if(PassIDCriteria(MuLColl.at(0), TightID)){
      FillHist("MuIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (IsNearB      )  FillHist("MuIDSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else                    FillHist("MuIDSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("MuIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("MuIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("MuIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(NBJets>0){
        FillHist("MuIDSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("MuIDSumW_BjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuIDSumW_CjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuIDSumW_LjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("MuIDSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("MuIDSumW_BjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuIDSumW_CjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuIDSumW_LjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }
    //======================================================================//

    if(fEta<0.9){
      FillHist("MuBSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if(NBJets>0) FillHist("MuBSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      else         FillHist("MuBSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (IsNearB      )  FillHist("MuBSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else                    FillHist("MuBSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("MuBSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("MuBSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("MuBSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
      if(PassIDCriteria(MuLColl.at(0), TightID)){
        FillHist("MuBIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if(NBJets>0) FillHist("MuBIDSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        else         FillHist("MuBIDSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (IsNearB      )  FillHist("MuBIDSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else                    FillHist("MuBIDSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("MuBIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuBIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuBIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }
    else if(fEta<1.6){
      FillHist("MuBESumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if(NBJets>0) FillHist("MuBESumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      else         FillHist("MuBESumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (IsNearB      )  FillHist("MuBESumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else                    FillHist("MuBESumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("MuBESumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("MuBESumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("MuBESumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
      if(PassIDCriteria(MuLColl.at(0), TightID)){
        FillHist("MuBEIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if(NBJets>0) FillHist("MuBEIDSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        else         FillHist("MuBEIDSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (IsNearB      )  FillHist("MuBEIDSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else                    FillHist("MuBEIDSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("MuBEIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuBEIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuBEIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }
    else{
      FillHist("MuESumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if(NBJets>0) FillHist("MuESumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      else         FillHist("MuESumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (IsNearB      )  FillHist("MuESumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else                    FillHist("MuESumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if     (SrcJetType==3)  FillHist("MuESumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("MuESumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("MuESumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
  
      if(PassIDCriteria(MuLColl.at(0), TightID)){
        FillHist("MuEIDSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if(NBJets>0) FillHist("MuEIDSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        else         FillHist("MuEIDSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (IsNearB      )  FillHist("MuEIDSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else                    FillHist("MuEIDSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("MuEIDSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuEIDSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuEIDSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
    }

    if(PassMETMTWReq){
      //For Composition Check==============================================//
      FillHist("MuMETMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
      if     (SrcJetType==3)  FillHist("MuMETMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==2)  FillHist("MuMETMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      else if(SrcJetType==1)  FillHist("MuMETMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      if(NBJets>0){
        FillHist("MuMETMTWSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (IsNearB      )  FillHist("MuMETMTWSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else                    FillHist("MuMETMTWSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("MuMETMTWSumW_BjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuMETMTWSumW_CjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuMETMTWSumW_LjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      else{
        FillHist("MuMETMTWSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("MuMETMTWSumW_BjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuMETMTWSumW_CjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuMETMTWSumW_LjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
      }
      if(PassIDCriteria(MuLColl.at(0), TightID)){
        FillHist("MuIDMETMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if     (IsNearB      )  FillHist("MuIDMETMTWSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else                    FillHist("MuIDMETMTWSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("MuIDMETMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuIDMETMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuIDMETMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if(NBJets>0){
          FillHist("MuIDMETMTWSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          if     (SrcJetType==3)  FillHist("MuIDMETMTWSumW_BjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("MuIDMETMTWSumW_CjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("MuIDMETMTWSumW_LjMatchHasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
        else{
          FillHist("MuIDMETMTWSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          if     (SrcJetType==3)  FillHist("MuIDMETMTWSumW_BjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("MuIDMETMTWSumW_CjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("MuIDMETMTWSumW_LjMatchNoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      //======================================================================//
  
      if(fEta<0.9){
        FillHist("MuBMETMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if(NBJets>0) FillHist("MuBMETMTWSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        else         FillHist("MuBMETMTWSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (IsNearB      )  FillHist("MuBMETMTWSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else                    FillHist("MuBMETMTWSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("MuBMETMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuBMETMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuBMETMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    
        if(PassIDCriteria(MuLColl.at(0), TightID)){
          FillHist("MuBIDMETMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if(NBJets>0) FillHist("MuBIDMETMTWSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          else         FillHist("MuBIDMETMTWSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          if     (IsNearB      )  FillHist("MuBIDMETMTWSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else                    FillHist("MuBIDMETMTWSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          if     (SrcJetType==3)  FillHist("MuBIDMETMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("MuBIDMETMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("MuBIDMETMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else if(fEta<1.6){
        FillHist("MuBEMETMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if(NBJets>0) FillHist("MuBEMETMTWSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        else         FillHist("MuBEMETMTWSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (IsNearB      )  FillHist("MuBEMETMTWSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else                    FillHist("MuBEMETMTWSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("MuBEMETMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuBEMETMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuBEMETMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    
        if(PassIDCriteria(MuLColl.at(0), TightID)){
          FillHist("MuBEIDMETMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if(NBJets>0) FillHist("MuBEIDMETMTWSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          else         FillHist("MuBEIDMETMTWSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          if     (IsNearB      )  FillHist("MuBEIDMETMTWSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else                    FillHist("MuBEIDMETMTWSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          if     (SrcJetType==3)  FillHist("MuBEIDMETMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("MuBEIDMETMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("MuBEIDMETMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
      else{
        FillHist("MuEMETMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        if(NBJets>0) FillHist("MuEMETMTWSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
        else         FillHist("MuEMETMTWSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (IsNearB      )  FillHist("MuEMETMTWSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else                    FillHist("MuEMETMTWSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        if     (SrcJetType==3)  FillHist("MuEMETMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==2)  FillHist("MuEMETMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        else if(SrcJetType==1)  FillHist("MuEMETMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
    
        if(PassIDCriteria(MuLColl.at(0), TightID)){
          FillHist("MuEIDMETMTWSumW_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          if(NBJets>0) FillHist("MuEIDMETMTWSumW_HasB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);  
          else         FillHist("MuEIDMETMTWSumW_NoB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          if     (IsNearB      )  FillHist("MuEIDMETMTWSumW_NearB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else                    FillHist("MuEIDMETMTWSumW_AwayB_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          if     (SrcJetType==3)  FillHist("MuEIDMETMTWSumW_BjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==2)  FillHist("MuEIDMETMTWSumW_CjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
          else if(SrcJetType==1)  FillHist("MuEIDMETMTWSumW_LjMatch_PT1D"+Label, PTCorr, weight, PtEdges, NPtEdges-1);
        }
      }
    }//End of METMTW
  }//End of TrigSel


}




void Aug2017_MuFakeMCStudy::CheckMCClosure(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KTruth> TruthColl, TString LooseID, TString TightID, float weight, TString Label, TString Option){

  TString FilterInfo="", ConeMethod="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";

  std::vector<snu::KMuon> MuLColl, MuTColl;
    for(int i=0; i<(int) MuPreColl.size(); i++){
      if(PassIDCriteria(MuPreColl.at(i), LooseID, "Roch")) MuLColl.push_back(MuPreColl.at(i));
      if(PassIDCriteria(MuPreColl.at(i), TightID, "Roch")) MuTColl.push_back(MuPreColl.at(i));
    }
  std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");
  std::vector<snu::KJet> BJetNoVetoColl = SelBJets(JetNoVetoColl, "Medium");
  float MET    = eventbase->GetEvent().MET(snu::KEvent::pfmet);
  float metphi = eventbase->GetEvent().METPhi();
  float METx=MET*cos(metphi), METy=MET*sin(metphi);
  TString FlavOpt="";
   
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


  for(int i=0; i<(int) MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), TightID)){
      if(Option.Contains("ParB")){ 
        FlavOpt="_AwayB";
        for(int j=0; j<(int) BJetNoVetoColl.size(); j++){
          if(BJetNoVetoColl.at(j).DeltaR(MuLColl.at(i))<0.4){ FlavOpt="_HasB"; break; }
        }
      }
      if(Option.Contains("GenB")){
        bool BjMatch=false, CjMatch=false; FlavOpt="";
        for(int j=0; j<(int) JetNoVetoColl.size(); j++){ 
          if(JetNoVetoColl.at(j).DeltaR(MuLColl.at(i))<0.4){
            int FlavIdx=JetNoVetoColl.at(j).HadronFlavour();
            if(FlavIdx==5                        ){ FlavOpt="_BjMatch"; break; }
            if(!BjMatch && FlavIdx==4            ){ FlavOpt="_CjMatch"; }
            if(!BjMatch && !CjMatch && FlavIdx==0){ FlavOpt="_LjMatch"; }
          }
        }
        if(k_sample_name.Contains("TT") && FlavOpt=="") FlavOpt="_BjMatch";
      }
      float FR=FakeRateMC(MuLColl.at(i), "QCD"+FlavOpt+"_"+TightID+"_"+LooseID+"_"+FilterInfo+ConeMethod);
      fakeweight*=-FR/(1-FR);
      //cout<<"FR: "<<FR<<" Weight: "<<fakeweight<<endl;
    }
  }
  if( MuLColl.size()==3 && MuTColl.size()==3 ) fakeweight=0;

  bool Pass_Trigger=false, Pass_TriggerBG=false, Pass_TriggerH=false;
  float trigger_ps_weight=1., trigger_period_weight=1.;
  float LumiBG=27.257618, LumiH=8.605696, LumiBH=35.863314;
  if(  PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")
     ||PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") )    Pass_TriggerBG=true;
  if(  PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")
     ||PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_TriggerH =true;

  if( Pass_TriggerBG || Pass_TriggerH ) Pass_Trigger=true;
  trigger_period_weight=( (Pass_TriggerBG? 27.257618:0.)+(Pass_TriggerH? 8.605696:0.) )/35.863314;
  weight*=trigger_period_weight;

  

  if(Pass_Trigger){
    if(MuLColl.size()==3 && EleLColl.size()==0){
      bool PassSel=true;
      int IdxOS=-1, IdxSS1=-1, IdxSS2=-1, IdxSSW=-1, IdxSSA=-1;
      float MOSSS1=-1, MOSSS2=-1, Mmumu=-1, MTW=-1.;
      if( !(MuLColl.at(0).Pt()>20 && MuLColl.at(1).Pt()>10 && MuLColl.at(2).Pt()>10) ) PassSel=false;
      if( PassSel && fabs(SumCharge(MuLColl))!=1) PassSel=false;
      if( PassSel ){
        IdxOS  = TriMuChargeIndex(MuLColl,"OS");
        IdxSS1 = TriMuChargeIndex(MuLColl,"SS1");
        IdxSS2 = TriMuChargeIndex(MuLColl,"SS2"); 
        IdxSSW = TriMuChargeIndex(MuLColl, MET, METx, METy, "SSW");
        IdxSSA = TriMuChargeIndex(MuLColl, MET, METx, METy, "SSA");
        MOSSS1 = (MuLColl.at(IdxOS)+MuLColl.at(IdxSS1)).M();
        MOSSS2 = (MuLColl.at(IdxOS)+MuLColl.at(IdxSS2)).M();
        Mmumu  = (MuLColl.at(IdxOS)+MuLColl.at(IdxSSA)).M();
        MTW    = sqrt(2)*sqrt(MET*MuLColl.at(IdxSSW).Pt()-METx*MuLColl.at(IdxSSW).Px()-METy*MuLColl.at(IdxSSW).Py());
      }
      if( PassSel && !(MOSSS1>12 && MOSSS2>12) ) PassSel=false;
      if(PassSel){
        bool JetSelPass=JetColl.size()>=2, BJetSelPass=BJetColl.size()>=1;
        bool OffZ=fabs(MOSSS1-91.2)>10. && fabs(MOSSS2-91.2)>10.;
        FillHist("CutFlow_exp"+Label, 0., weight*fakeweight, 0., 10., 10);
        if(OffZ && JetSelPass) FillHist("CutFlow_exp"+Label, 1., weight*fakeweight, 0., 10., 10);
        if(OffZ && JetSelPass && BJetSelPass) FillHist("CutFlow_exp"+Label, 2., weight*fakeweight, 0., 10., 10);

        FillHist("PTmu1_exp"+Label, MuLColl.at(0).Pt(), weight*fakeweight, 0., 200., 40);
        FillHist("PTmu2_exp"+Label, MuLColl.at(1).Pt(), weight*fakeweight, 0., 200., 40);
        FillHist("PTmu3_exp"+Label, MuLColl.at(2).Pt(), weight*fakeweight, 0., 200., 40);
        FillHist("MOSSS1_exp"+Label, MOSSS1, weight*fakeweight, 0., 300., 60);
        FillHist("MOSSS2_exp"+Label, MOSSS2, weight*fakeweight, 0., 200., 40);
        FillHist("Mmumu_exp" +Label, Mmumu , weight*fakeweight, 0., 200., 40);
        FillHist("Etamu1_exp"+Label, MuLColl.at(0).Eta(), weight*fakeweight, -5., 5., 20);
        FillHist("Etamu2_exp"+Label, MuLColl.at(1).Eta(), weight*fakeweight, -5., 5., 20);
        FillHist("Etamu3_exp"+Label, MuLColl.at(2).Eta(), weight*fakeweight, -5., 5., 20);
        FillHist("Nj_exp"+Label, JetColl.size(), weight*fakeweight, 0., 10., 10);
        FillHist("Nb_exp"+Label, BJetColl.size(), weight*fakeweight, 0., 10., 10);
        FillHist("MET_exp"+Label, MET, weight, 0., 300., 30);
        FillHist("MTW_exp"+Label, MTW, weight, 0., 200., 20);
        if(JetSelPass && BJetSelPass && OffZ){
          FillHist("PTmu1_1b2j_exp"+Label, MuLColl.at(0).Pt(), weight*fakeweight, 0., 200., 40);
          FillHist("PTmu2_1b2j_exp"+Label, MuLColl.at(1).Pt(), weight*fakeweight, 0., 200., 40);
          FillHist("PTmu3_1b2j_exp"+Label, MuLColl.at(2).Pt(), weight*fakeweight, 0., 200., 40);
          FillHist("MOSSS1_1b2j_exp"+Label, MOSSS1, weight*fakeweight, 0., 300., 60);
          FillHist("MOSSS2_1b2j_exp"+Label, MOSSS2, weight*fakeweight, 0., 200., 40);
          FillHist("Mmumu_1b2j_exp" +Label, Mmumu , weight*fakeweight, 0., 200., 40);
          FillHist("Etamu1_1b2j_exp"+Label, MuLColl.at(0).Eta(), weight*fakeweight, -5., 5., 20);
          FillHist("Etamu2_1b2j_exp"+Label, MuLColl.at(1).Eta(), weight*fakeweight, -5., 5., 20);
          FillHist("Etamu3_1b2j_exp"+Label, MuLColl.at(2).Eta(), weight*fakeweight, -5., 5., 20);
          FillHist("Nj_1b2j_exp"+Label, JetColl.size(), weight*fakeweight, 0., 10., 10);
          FillHist("Nb_1b2j_exp"+Label, BJetColl.size(), weight*fakeweight, 0., 10., 10);
          FillHist("MET_1b2j_exp"+Label, MET, weight, 0., 300., 30);
          FillHist("MTW_1b2j_exp"+Label, MTW, weight, 0., 200., 20);
        }
      }
    }
    if(MuLColl.size()==3 && MuTColl.size()==3 && EleLColl.size()==0){
      bool PassSel=true;
      int IdxOS=-1, IdxSS1=-1, IdxSS2=-1, IdxSSA=-1, IdxSSW=-1; float MOSSS1=-1, MOSSS2=-1, Mmumu=-1, MTW=-1;
      if( !(MuLColl.at(0).Pt()>20 && MuLColl.at(1).Pt()>10 && MuLColl.at(2).Pt()>10) ) PassSel=false;
      if( PassSel && fabs(SumCharge(MuLColl))!=1 ) PassSel=false;
      if( PassSel ){
        IdxOS  = TriMuChargeIndex(MuLColl,"OS");
        IdxSS1 = TriMuChargeIndex(MuLColl,"SS1");
        IdxSS2 = TriMuChargeIndex(MuLColl,"SS2"); 
        IdxSSW = TriMuChargeIndex(MuLColl, MET, METx, METy, "SSW");
        IdxSSA = TriMuChargeIndex(MuLColl, MET, METx, METy, "SSA");
        MOSSS1 = (MuLColl.at(IdxOS)+MuLColl.at(IdxSS1)).M();
        MOSSS2 = (MuLColl.at(IdxOS)+MuLColl.at(IdxSS2)).M();
        Mmumu  = (MuLColl.at(IdxOS)+MuLColl.at(IdxSSA)).M();
        MTW    = sqrt(2)*sqrt(MET*MuLColl.at(IdxSSW).Pt()-METx*MuLColl.at(IdxSSW).Px()-METy*MuLColl.at(IdxSSW).Py());
      }
      if( PassSel && !(MOSSS1>12 && MOSSS2>12) ) PassSel=false;
      if( PassSel ){
        bool JetSelPass=JetColl.size()>=2, BJetSelPass=BJetColl.size()>=1;
        bool OffZ=fabs(MOSSS1-91.2)>10. && fabs(MOSSS2-91.2)>10.;
        FillHist("CutFlow_obs"+Label, 0., weight, 0., 10., 10);
        if(OffZ && JetSelPass) FillHist("CutFlow_obs"+Label, 1., weight, 0., 10., 10);
        if(OffZ && JetSelPass && BJetSelPass) FillHist("CutFlow_obs"+Label, 2., weight, 0., 10., 10);

        int LepTypeOS  = GetLeptonType(MuLColl.at(IdxOS),TruthColl);
        int LepTypeSSW = GetLeptonType(MuLColl.at(IdxSSW),TruthColl);
        int LepTypeSSA = GetLeptonType(MuLColl.at(IdxSSA),TruthColl);

        if(OffZ){
          if     ( LepTypeOS>0 && LepTypeSSA>0 ) FillHist("MuMuType"+Label, 0., weight, 0., 10., 10);
          else if( LepTypeOS>0 && LepTypeSSA<0 ) FillHist("MuMuType"+Label, 1., weight, 0., 10., 10);
          else                                   FillHist("MuMuType"+Label, 2., weight, 0., 10., 10);
        }

        FillHist("PTmu1_obs"+Label, MuLColl.at(0).Pt(), weight, 0., 200., 40);
        FillHist("PTmu2_obs"+Label, MuLColl.at(1).Pt(), weight, 0., 200., 40);
        FillHist("PTmu3_obs"+Label, MuLColl.at(2).Pt(), weight, 0., 200., 40);
        FillHist("MOSSS1_obs"+Label, MOSSS1, weight, 0., 300., 60);
        FillHist("MOSSS2_obs"+Label, MOSSS2, weight, 0., 200., 40);
        FillHist("Mmumu_obs" +Label, Mmumu , weight, 0., 200., 40);
        FillHist("Etamu1_obs"+Label, MuLColl.at(0).Eta(), weight, -5., 5., 20);
        FillHist("Etamu2_obs"+Label, MuLColl.at(1).Eta(), weight, -5., 5., 20);
        FillHist("Etamu3_obs"+Label, MuLColl.at(2).Eta(), weight, -5., 5., 20);
        FillHist("Nj_obs"+Label, JetColl.size(), weight, 0., 10., 10);
        FillHist("Nb_obs"+Label, BJetColl.size(), weight, 0., 10., 10);
        FillHist("MET_obs"+Label, MET, weight, 0., 300., 30);
        FillHist("MTW_obs"+Label, MTW, weight, 0., 200., 20);
        if(JetSelPass && BJetSelPass && OffZ){
          FillHist("PTmu1_1b2j_obs"+Label, MuLColl.at(0).Pt(), weight, 0., 200., 40);
          FillHist("PTmu2_1b2j_obs"+Label, MuLColl.at(1).Pt(), weight, 0., 200., 40);
          FillHist("PTmu3_1b2j_obs"+Label, MuLColl.at(2).Pt(), weight, 0., 200., 40);
          FillHist("MOSSS1_1b2j_obs"+Label, MOSSS1, weight, 0., 300., 60);
          FillHist("MOSSS2_1b2j_obs"+Label, MOSSS2, weight, 0., 200., 40);
          FillHist("Mmumu_1b2j_obs" +Label, MOSSS2, weight, 0., 200., 40);
          FillHist("Etamu1_1b2j_obs"+Label, MuLColl.at(0).Eta(), weight, -5., 5., 20);
          FillHist("Etamu2_1b2j_obs"+Label, MuLColl.at(1).Eta(), weight, -5., 5., 20);
          FillHist("Etamu3_1b2j_obs"+Label, MuLColl.at(2).Eta(), weight, -5., 5., 20);
          FillHist("Nj_1b2j_obs"+Label, JetColl.size(), weight, 0., 10., 10);
          FillHist("Nb_1b2j_obs"+Label, BJetColl.size(), weight, 0., 10., 10);
          FillHist("MET_1b2j_obs"+Label, MET, weight, 0., 300., 30);
          FillHist("MTW_1b2j_obs"+Label, MTW, weight, 0., 200., 20);
        }
      }
    }
  } 
}



void Aug2017_MuFakeMCStudy::CheckMotherDaugherRelationship(std::vector<snu::KMuon> muonColl, std::vector<snu::KMuon> muonLooseColl, std::vector<snu::KElectron> electronLooseColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, std::vector<snu::KTruth> truthColl, std::vector<snu::KGenJet> GenJetColl, float weight, TString Label, TString Option){

  for(int i=0; i<(int) muonLooseColl.size(); i++){
    float PTCorr     = muonLooseColl.at(i).Pt()*(1+RochIso(muonLooseColl.at(i), "0.4"));
    float PT         = muonLooseColl.at(i).Pt();
    float fEta       = fabs(muonLooseColl.at(i).Eta());
    int   LepType    = GetLeptonType(muonLooseColl.at(i),truthColl);
    int   SrcJetType = GetFakeLepJetSrcType(muonLooseColl.at(i), JetColl);
    int   SrcJetIdx04   = GetFakeLepJetSrcIdx(muonLooseColl.at(i), JetColl, "DR04");
    int   SrcJetIdxNear = GetFakeLepJetSrcIdx(muonLooseColl.at(i), JetColl, "Near");
    int   SrcGenJetIdx04   = GetFakeLepGenJetSrcIdx(muonLooseColl.at(i), GenJetColl, "DR04");
    int   SrcGenJetIdxNear = GetFakeLepGenJetSrcIdx(muonLooseColl.at(i), GenJetColl, "Near");


    if(LepType<-4 || LepType>0) continue;
    
    FillHist("NMatchedJet", -1., weight, -1., 2., 3);
    FillHist("Fk_RelIso", muonLooseColl.at(i).RelIso04(), weight, 0., 10., 500);
    if(SrcJetIdxNear!=-1){
      FillHist("dRlj_MatchNear"       , muonLooseColl.at(i).DeltaR(JetColl.at(SrcJetIdxNear)), weight, 0., 5., 100);
    }
    if(SrcJetIdx04!=-1){
      FillHist("NMatchedJet", 1., weight, -1., 2., 3);
      FillHist("Fk_RelIso_Match04"  , muonLooseColl.at(i).RelIso04(), weight, 0., 10., 500);
      FillHist("Fk_RelPtJet_Match04", muonLooseColl.at(i).Pt()/JetColl.at(SrcJetIdx04).Pt(), weight, 0., 2., 40);
      FillHist("dRlj_Match04"       , muonLooseColl.at(i).DeltaR(JetColl.at(SrcJetIdx04)), weight, 0., 5., 100);
      FillHist("PTMu_Match04"       , muonLooseColl.at(i).Pt(), weight, 0., 200., 40);
      FillHist("PTJet_Match04"      , JetColl.at(SrcJetIdx04).Pt(), weight, 0., 500., 100);
      FillHist("EtaMu_Match04"      , muonLooseColl.at(i).Eta()    , weight, -5., 5., 20);
      FillHist("EtaJet_Match04"     , JetColl.at(SrcJetIdx04).Eta(), weight, -5., 5., 20);
    }
    else{
      FillHist("NMatchedJet", 0., weight, -1., 2., 3);
    }

    if(SrcGenJetIdxNear!=-1){
      FillHist("dRlj_GenMatchNear"       , muonLooseColl.at(i).DeltaR(GenJetColl.at(SrcGenJetIdxNear)), weight, 0., 5., 100);
    }
    if(SrcGenJetIdx04!=-1){
      FillHist("NMatchedGenJet", 1., weight, -1., 2., 3);
      FillHist("Fk_RelIso_GenMatch04"  , muonLooseColl.at(i).RelIso04(), weight, 0., 10., 500);
      FillHist("Fk_RelPtGenJet_GenMatch04", muonLooseColl.at(i).Pt()/GenJetColl.at(SrcGenJetIdx04).Pt(), weight, 0., 2., 40);
      FillHist("dRlj_GenMatch04"       , muonLooseColl.at(i).DeltaR(GenJetColl.at(SrcGenJetIdx04)), weight, 0., 5., 100);
    }
    else{
      FillHist("NMatchedGenJet", 0., weight, -1., 2., 3);
    }


    if(PassIDCriteria(muonLooseColl.at(i), "POGTIsop20IPp01p05sig4Chi4")){
      FillHist("NMatchedJet_ID", -1., weight, -1., 2., 3);
      FillHist("Fk_RelIso_ID", muonLooseColl.at(i).RelIso04(), weight, 0., 10., 500);
      if(SrcJetIdxNear!=-1){
        FillHist("dRlj_MatchNear_ID"       , muonLooseColl.at(i).DeltaR(JetColl.at(SrcJetIdxNear)), weight, 0., 5., 100);
      }
      if(SrcJetIdx04!=-1){
        FillHist("NMatchedJet_ID", 1., weight, -1., 2., 3);
        FillHist("Fk_RelIso_Match04_ID"  , muonLooseColl.at(i).RelIso04(), weight, 0., 10., 500);
        FillHist("Fk_RelPtJet_Match04_ID", muonLooseColl.at(i).Pt()/JetColl.at(SrcJetIdx04).Pt(), weight, 0., 2., 40);
        FillHist("dRlj_Match04_ID"       , muonLooseColl.at(i).DeltaR(JetColl.at(SrcJetIdx04)), weight, 0., 5., 100);
        FillHist("PTMu_Match04_ID"       , muonLooseColl.at(i).Pt(), weight, 0., 200., 40);
        FillHist("PTJet_Match04_ID"      , JetColl.at(SrcJetIdx04).Pt(), weight, 0., 500., 100);
        FillHist("EtaMu_Match04_ID"      , muonLooseColl.at(i).Eta()    , weight, -5., 5., 20);
        FillHist("EtaJet_Match04_ID"     , JetColl.at(SrcJetIdx04).Eta(), weight, -5., 5., 20);
      }
      else{
        FillHist("NMatchedJet_ID", 0., weight, -1., 2., 3);
      }


      if(SrcGenJetIdxNear!=-1){
        FillHist("dRlj_GenMatchNear_ID"       , muonLooseColl.at(i).DeltaR(GenJetColl.at(SrcGenJetIdxNear)), weight, 0., 5., 100);
      }
      if(SrcGenJetIdx04!=-1){
        FillHist("NMatchedGenJet_ID", 1., weight, -1., 2., 3);
        FillHist("Fk_RelIso_GenMatch04_ID"  , muonLooseColl.at(i).RelIso04(), weight, 0., 10., 500);
        FillHist("Fk_RelPtGenJet_GenMatch04_ID", muonLooseColl.at(i).Pt()/GenJetColl.at(SrcGenJetIdx04).Pt(), weight, 0., 2., 40);
        FillHist("dRlj_GenMatch04_ID"       , muonLooseColl.at(i).DeltaR(GenJetColl.at(SrcGenJetIdx04)), weight, 0., 5., 100);
      }
      else{
        FillHist("NMatchedGenJet_ID", 0., weight, -1., 2., 3);
      }

    }
  }//End of ID Variable

}

void Aug2017_MuFakeMCStudy::CheckAltCRAvailability(std::vector<snu::KMuon> MuPreColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, float MET, float METx, float METy, std::vector<snu::KTruth> TruthColl, TString LooseID, TString TightID, float weight, TString Label, TString Option){


  std::vector<snu::KMuon> MuLColl, MuTColl;
    for(int i=0; i<(int) MuPreColl.size(); i++){
      if(PassIDCriteria(MuPreColl.at(i), LooseID, "Roch")) MuLColl.push_back(MuPreColl.at(i));
      if(PassIDCriteria(MuPreColl.at(i), TightID, "Roch")) MuTColl.push_back(MuPreColl.at(i));
    }
  std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");
//  std::vector<snu::KJet> BJetNoVetoColl = SelBJets(JetNoVetoColl, "Medium");
   
  //Checking 2PromptMu+1Fk Ele
  float fakeweight = -1.;
  int   NHFakeMu   = NPromptFake_Mu(MuLColl, TruthColl, "HFake");
  int   NHFakeEle  = NPromptFake_Ele(EleLColl, TruthColl, "HFake");
  int   NHFakeLep  = NHFakeMu+NHFakeEle;
  int   NLooseLep  = MuLColl.size()+EleLColl.size();

  //Use only 3 valid leptons 
  bool DiLCand = NLooseLep==2;
  bool TriLCand = NLooseLep==3;
  if( !(DiLCand || TriLCand) ) return;//Include at least 1 fake in the sample

  if(TriLCand && MuTColl.size()==3){
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    int   IdxSSW = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSW");
    int   IdxSSA = TriMuChargeIndex(MuTColl, MET, METx, METy, "SSA");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    float Mmumu  = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSA)).M();

    if( !(MOSSS1>12 && MOSSS2>12) ) return;
  
      FillHist("Mmumu_3lOS" +Label, Mmumu , weight, 0., 200., 200);
      FillHist("MOSSS1_3lOS"+Label, MOSSS1, weight, 0., 200., 200);
      FillHist("MOSSS2_3lOS"+Label, MOSSS2, weight, 0., 200., 200);
      if( fabs(MOSSS1-91.2)>10. ) FillHist("MOSSS2_OSSS1OffZ_3lOS"+Label, MOSSS2, weight, 0., 200., 200);
    if( fabs(MOSSS1-91.2)<10 || fabs(MOSSS2-91.2)<10 ) return;
      FillHist("Mmumu_3lOSOffZ" +Label, Mmumu , weight, 0., 200., 200);
      FillHist("MOSSS1_3lOSOffZ"+Label, MOSSS1, weight, 0., 200., 200);
      FillHist("MOSSS2_3lOSOffZ"+Label, MOSSS2, weight, 0., 200., 200);
  }
  if(DiLCand && MuTColl.size()==2){
    float Mmumu = (MuTColl.at(0)+MuTColl.at(1)).M();
    float PTMu1 = MuTColl.at(0).Pt(), PTMu2 = MuTColl.at(1).Pt();
    int njets= JetColl.size(), nbjets= BJetColl.size();
    if(Mmumu<12.) return;
    if(NHFakeLep==0){
      if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PP2010_Incl"+Label, Mmumu, weight, 0., 200., 200);
      if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PP1010_Incl"+Label, Mmumu, weight, 0., 200., 200);
    }
    if(NHFakeLep==1){
      if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PF2010_Incl"+Label, Mmumu, weight, 0., 200., 200);
      if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PF1010_Incl"+Label, Mmumu, weight, 0., 200., 200);
    }
    if(njets==0){
      if(NHFakeLep==0){
        if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PP2010_0j"+Label, Mmumu, weight, 0., 200., 200);
        if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PP1010_0j"+Label, Mmumu, weight, 0., 200., 200);
      }
      if(NHFakeLep==1){
        if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PF2010_0j"+Label, Mmumu, weight, 0., 200., 200);
        if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PF1010_0j"+Label, Mmumu, weight, 0., 200., 200);
      }
    }
    if(njets==1){
      if(NHFakeLep==0){
        if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PP2010_1j"+Label, Mmumu, weight, 0., 200., 200);
        if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PP1010_1j"+Label, Mmumu, weight, 0., 200., 200);
      }
      if(NHFakeLep==1){
        if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PF2010_1j"+Label, Mmumu, weight, 0., 200., 200);
        if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PF1010_1j"+Label, Mmumu, weight, 0., 200., 200);
      }
    }
    if(njets==2){
      if(NHFakeLep==0){
        if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PP2010_2j"+Label, Mmumu, weight, 0., 200., 200);
        if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PP1010_2j"+Label, Mmumu, weight, 0., 200., 200);
      }
      if(NHFakeLep==1){
        if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PF2010_2j"+Label, Mmumu, weight, 0., 200., 200);
        if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PF1010_2j"+Label, Mmumu, weight, 0., 200., 200);
      }
    }
    if(njets==3){
      if(NHFakeLep==0){
        if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PP2010_3j"+Label, Mmumu, weight, 0., 200., 200);
        if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PP1010_3j"+Label, Mmumu, weight, 0., 200., 200);
      }
      if(NHFakeLep==1){
        if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PF2010_3j"+Label, Mmumu, weight, 0., 200., 200);
        if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PF1010_3j"+Label, Mmumu, weight, 0., 200., 200);
      }
    }
    if(nbjets==1){
      if(NHFakeLep==0){
        if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PP2010_1b"+Label, Mmumu, weight, 0., 200., 200);
        if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PP1010_1b"+Label, Mmumu, weight, 0., 200., 200);
      }
      if(NHFakeLep==1){
        if(PTMu1>20. && PTMu2>10.) FillHist("Mmumu_DiL_PF2010_1b"+Label, Mmumu, weight, 0., 200., 200);
        if(PTMu1>10. && PTMu2>10.) FillHist("Mmumu_DiL_PF1010_1b"+Label, Mmumu, weight, 0., 200., 200);
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
      if     (PTCorr<15)  FR=0.31422 ;
      else if(PTCorr<20)  FR=0.2225  ;
      else if(PTCorr<25)  FR=0.201252;
      else if(PTCorr<35)  FR=0.174787;
      else if(PTCorr<50)  FR=0.144774;
      else                FR=0.137056;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.328501;
      else if(PTCorr<20)  FR=0.250262;
      else if(PTCorr<25)  FR=0.239327;
      else if(PTCorr<35)  FR=0.204919;
      else if(PTCorr<50)  FR=0.200699;
      else                FR=0.188323;
    }
    else{
      if     (PTCorr<15)  FR=0.369477;
      else if(PTCorr<20)  FR=0.29852 ;
      else if(PTCorr<25)  FR=0.28201 ;
      else if(PTCorr<35)  FR=0.271146;
      else if(PTCorr<50)  FR=0.238203;
      else                FR=0.234892;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.319891;
      else if(PTCorr<20)  FR=0.246482;
      else if(PTCorr<25)  FR=0.249877;
      else if(PTCorr<35)  FR=0.223699;
      else if(PTCorr<50)  FR=0.195308;
      else                FR=0.18359 ;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.338585;
      else if(PTCorr<20)  FR=0.286107;
      else if(PTCorr<25)  FR=0.297202;
      else if(PTCorr<35)  FR=0.262352;
      else if(PTCorr<50)  FR=0.257737;
      else                FR=0.252453;
    }
    else{
      if     (PTCorr<15)  FR=0.385327;
      else if(PTCorr<20)  FR=0.345415;
      else if(PTCorr<25)  FR=0.336095;
      else if(PTCorr<35)  FR=0.333724;
      else if(PTCorr<50)  FR=0.298049;
      else                FR=0.288971;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
//    if(fEta<0.9){//HasB
//      if     (PTCorr<15)  FR=0.300815;
//      else if(PTCorr<20)  FR=0.224423;
//      else if(PTCorr<25)  FR=0.195447;
//      else if(PTCorr<35)  FR=0.170102;
//      else if(PTCorr<50)  FR=0.137973;
//      else                FR=0.121946;
//    }
//    else if(fEta<1.6){
//      if     (PTCorr<15)  FR=0.339909;
//      else if(PTCorr<20)  FR=0.247718;
//      else if(PTCorr<25)  FR=0.238802;
//      else if(PTCorr<35)  FR=0.203254;
//      else if(PTCorr<50)  FR=0.209858;
//      else                FR=0.190979;
//    }
//    else{
//      if     (PTCorr<15)  FR=0.363888;
//      else if(PTCorr<20)  FR=0.293596;
//      else if(PTCorr<25)  FR=0.291579;
//      else if(PTCorr<35)  FR=0.285735;
//      else if(PTCorr<50)  FR=0.252648;
//      else                FR=0.241945;
//    }
    if(fEta<0.9){//NearB
      if     (PTCorr<15)  FR=0.082816 ;
      else if(PTCorr<20)  FR=0.0835979;
      else if(PTCorr<25)  FR=0.096017 ;
      else if(PTCorr<35)  FR=0.0864538;
      else if(PTCorr<50)  FR=0.0610422;
      else                FR=0.0662293;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.0954808;
      else if(PTCorr<20)  FR=0.100042 ;
      else if(PTCorr<25)  FR=0.108263 ;
      else if(PTCorr<35)  FR=0.0848722;
      else if(PTCorr<50)  FR=0.120748 ;
      else                FR=0.0643254;
    }
    else{
      if     (PTCorr<15)  FR=0.0803636;
      else if(PTCorr<20)  FR=0.0874016;
      else if(PTCorr<25)  FR=0.113497 ;
      else if(PTCorr<35)  FR=0.114078 ;
      else if(PTCorr<50)  FR=0.103653 ;
      else                FR=0.167258 ;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.329892;
      else if(PTCorr<20)  FR=0.231877;
      else if(PTCorr<25)  FR=0.210858;
      else if(PTCorr<35)  FR=0.18123 ;
      else if(PTCorr<50)  FR=0.150159;
      else                FR=0.142378;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.343916;
      else if(PTCorr<20)  FR=0.260858;
      else if(PTCorr<25)  FR=0.246079;
      else if(PTCorr<35)  FR=0.210801;
      else if(PTCorr<50)  FR=0.20615 ;
      else                FR=0.193187;
    }
    else{
      if     (PTCorr<15)  FR=0.372854;
      else if(PTCorr<20)  FR=0.301522;
      else if(PTCorr<25)  FR=0.287546;
      else if(PTCorr<35)  FR=0.274437;
      else if(PTCorr<50)  FR=0.239801;
      else                FR=0.237417;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.336205;
      else if(PTCorr<20)  FR=0.258356;
      else if(PTCorr<25)  FR=0.266948;
      else if(PTCorr<35)  FR=0.236031;
      else if(PTCorr<50)  FR=0.206611;
      else                FR=0.195049;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.354954;
      else if(PTCorr<20)  FR=0.300419;
      else if(PTCorr<25)  FR=0.308447;
      else if(PTCorr<35)  FR=0.273022;
      else if(PTCorr<50)  FR=0.268583;
      else                FR=0.261656;
    }
    else{
      if     (PTCorr<15)  FR=0.388919;
      else if(PTCorr<20)  FR=0.349613;
      else if(PTCorr<25)  FR=0.344924;
      else if(PTCorr<35)  FR=0.33913 ;
      else if(PTCorr<50)  FR=0.300837;
      else                FR=0.293124;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1NoChi_TrkIsoVVLConeSUSY"){
//    if(fEta<0.9){//HasB
//      if     (PTCorr<15)  FR=0.305642;
//      else if(PTCorr<20)  FR=0.227389;
//      else if(PTCorr<25)  FR=0.198893;
//      else if(PTCorr<35)  FR=0.172341;
//      else if(PTCorr<50)  FR=0.14026 ;
//      else                FR=0.123841;
//    }
//    else if(fEta<1.6){
//      if     (PTCorr<15)  FR=0.346488;
//      else if(PTCorr<20)  FR=0.252103;
//      else if(PTCorr<25)  FR=0.241488;
//      else if(PTCorr<35)  FR=0.20549 ;
//      else if(PTCorr<50)  FR=0.211249;
//      else                FR=0.192561;
//    }
//    else{
//      if     (PTCorr<15)  FR=0.366058;
//      else if(PTCorr<20)  FR=0.294558;
//      else if(PTCorr<25)  FR=0.292205;
//      else if(PTCorr<35)  FR=0.286885;
//      else if(PTCorr<50)  FR=0.253524;
//      else                FR=0.242219;
//    }

    if(fEta<0.9){//NearB
      if     (PTCorr<15)  FR=0.0836733;
      else if(PTCorr<20)  FR=0.0839263;
      else if(PTCorr<25)  FR=0.0964101;
      else if(PTCorr<35)  FR=0.087111 ;
      else if(PTCorr<50)  FR=0.0614359;
      else                FR=0.0665797;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.098713 ;
      else if(PTCorr<20)  FR=0.10115  ;
      else if(PTCorr<25)  FR=0.108608 ;
      else if(PTCorr<35)  FR=0.0854134;
      else if(PTCorr<50)  FR=0.121014 ;
      else                FR=0.0648354;
    }
    else{
      if     (PTCorr<15)  FR=0.08079  ;
      else if(PTCorr<20)  FR=0.0875968;
      else if(PTCorr<25)  FR=0.11368  ;
      else if(PTCorr<35)  FR=0.114311 ;
      else if(PTCorr<50)  FR=0.103623 ;
      else                FR=0.167314 ;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.271221 ;
      else if(PTCorr<20)  FR=0.134619 ;
      else if(PTCorr<25)  FR=0.110795 ;
      else if(PTCorr<35)  FR=0.0920508;
      else if(PTCorr<50)  FR=0.0693127;
      else                FR=0.0571533;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.291793 ;
      else if(PTCorr<20)  FR=0.161946 ;
      else if(PTCorr<25)  FR=0.145827 ;
      else if(PTCorr<35)  FR=0.117371 ;
      else if(PTCorr<50)  FR=0.103942 ;
      else                FR=0.0856054;
    }
    else{
      if     (PTCorr<15)  FR=0.330455;
      else if(PTCorr<20)  FR=0.204726;
      else if(PTCorr<25)  FR=0.18639 ;
      else if(PTCorr<35)  FR=0.172681;
      else if(PTCorr<50)  FR=0.141964;
      else                FR=0.126728;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.277055 ;
      else if(PTCorr<20)  FR=0.149739 ;
      else if(PTCorr<25)  FR=0.145405 ;
      else if(PTCorr<35)  FR=0.124283 ;
      else if(PTCorr<50)  FR=0.0978479;
      else                FR=0.0833598;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.301852;
      else if(PTCorr<20)  FR=0.185629;
      else if(PTCorr<25)  FR=0.185122;
      else if(PTCorr<35)  FR=0.156936;
      else if(PTCorr<50)  FR=0.137329;
      else                FR=0.119177;
    }
    else{
      if     (PTCorr<15)  FR=0.345734;
      else if(PTCorr<20)  FR=0.23771 ;
      else if(PTCorr<25)  FR=0.227932;
      else if(PTCorr<35)  FR=0.216493;
      else if(PTCorr<50)  FR=0.182448;
      else                FR=0.162061;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1NoChi_TrkIsoVVLConeSUSY"){
//    if(fEta<0.9){//HasB
//      if     (PTCorr<15)  FR=0.253064 ;
//      else if(PTCorr<20)  FR=0.132948 ;
//      else if(PTCorr<25)  FR=0.106011 ;
//      else if(PTCorr<35)  FR=0.0867408;
//      else if(PTCorr<50)  FR=0.0583407;
//      else                FR=0.0515706;
//    }
//    else if(fEta<1.6){
//      if     (PTCorr<15)  FR=0.299641 ;
//      else if(PTCorr<20)  FR=0.164955 ;
//      else if(PTCorr<25)  FR=0.147266 ;
//      else if(PTCorr<35)  FR=0.116065 ;
//      else if(PTCorr<50)  FR=0.112543 ;
//      else                FR=0.0944257;
//    }
//    else{
//      if     (PTCorr<15)  FR=0.325495;
//      else if(PTCorr<20)  FR=0.207752;
//      else if(PTCorr<25)  FR=0.201676;
//      else if(PTCorr<35)  FR=0.185382;
//      else if(PTCorr<50)  FR=0.155966;
//      else                FR=0.134836;
//    }
    if(fEta<0.9){//NearB
      if     (PTCorr<15)  FR=0.0623887;
      else if(PTCorr<20)  FR=0.0482613;
      else if(PTCorr<25)  FR=0.0448971;
      else if(PTCorr<35)  FR=0.0395998;
      else if(PTCorr<50)  FR=0.0245074;
      else                FR=0.0247607;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.0779777;
      else if(PTCorr<20)  FR=0.0610952;
      else if(PTCorr<25)  FR=0.0610987;
      else if(PTCorr<35)  FR=0.0423056;
      else if(PTCorr<50)  FR=0.0533732;
      else                FR=0.024341 ;
    }
    else{
      if     (PTCorr<15)  FR=0.0651505;
      else if(PTCorr<20)  FR=0.0577962;
      else if(PTCorr<25)  FR=0.0640869;
      else if(PTCorr<35)  FR=0.0668677;
      else if(PTCorr<50)  FR=0.0568092;
      else                FR=0.0828605;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp01p05sig4Chi4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.509318;
      else if(PTCorr<20)  FR=0.359136;
      else if(PTCorr<25)  FR=0.327135;
      else if(PTCorr<35)  FR=0.282876;
      else if(PTCorr<50)  FR=0.245291;
      else                FR=0.230259;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.54671 ;
      else if(PTCorr<20)  FR=0.413858;
      else if(PTCorr<25)  FR=0.380258;
      else if(PTCorr<35)  FR=0.32436 ;
      else if(PTCorr<50)  FR=0.328264;
      else                FR=0.313919;
    }
    else{
      if     (PTCorr<15)  FR=0.613294;
      else if(PTCorr<20)  FR=0.489137;
      else if(PTCorr<25)  FR=0.456957;
      else if(PTCorr<35)  FR=0.445719;
      else if(PTCorr<50)  FR=0.387913;
      else                FR=0.394796;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp01p05sig4Chi4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.509747;
      else if(PTCorr<20)  FR=0.36514 ;
      else if(PTCorr<25)  FR=0.343122;
      else if(PTCorr<35)  FR=0.303122;
      else if(PTCorr<50)  FR=0.272945;
      else                FR=0.255955;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.547486;
      else if(PTCorr<20)  FR=0.420919;
      else if(PTCorr<25)  FR=0.392663;
      else if(PTCorr<35)  FR=0.344241;
      else if(PTCorr<50)  FR=0.34337 ;
      else                FR=0.3426  ;
    }
    else{
      if     (PTCorr<15)  FR=0.614017;
      else if(PTCorr<20)  FR=0.494134;
      else if(PTCorr<25)  FR=0.460323;
      else if(PTCorr<35)  FR=0.457226;
      else if(PTCorr<50)  FR=0.407284;
      else                FR=0.414266;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp01p05sig4Chi4_TrkIsoVVLConeSUSY"){
//    if(fEta<0.9){//HasB
//      if     (PTCorr<15)  FR=0.498126;
//      else if(PTCorr<20)  FR=0.358516;
//      else if(PTCorr<25)  FR=0.316968;
//      else if(PTCorr<35)  FR=0.277669;
//      else if(PTCorr<50)  FR=0.229729;
//      else                FR=0.20088 ;
//    }
//    else if(fEta<1.6){
//      if     (PTCorr<15)  FR=0.566062;
//      else if(PTCorr<20)  FR=0.418663;
//      else if(PTCorr<25)  FR=0.379795;
//      else if(PTCorr<35)  FR=0.334079;
//      else if(PTCorr<50)  FR=0.336888;
//      else                FR=0.304086;
//    }
//    else{
//      if     (PTCorr<15)  FR=0.637943;
//      else if(PTCorr<20)  FR=0.49632 ;
//      else if(PTCorr<25)  FR=0.492386;
//      else if(PTCorr<35)  FR=0.465106;
//      else if(PTCorr<50)  FR=0.412026;
//      else                FR=0.438496;
//    }

    //NearB
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.383756;
      else if(PTCorr<20)  FR=0.279563;
      else if(PTCorr<25)  FR=0.256573;
      else if(PTCorr<35)  FR=0.210008;
      else if(PTCorr<50)  FR=0.149417;
      else                FR=0.169035;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.420788;
      else if(PTCorr<20)  FR=0.349621;
      else if(PTCorr<25)  FR=0.309854;
      else if(PTCorr<35)  FR=0.23608 ;
      else if(PTCorr<50)  FR=0.32618 ;
      else                FR=0.181229;
    }
    else{
      if     (PTCorr<15)  FR=0.499587;
      else if(PTCorr<20)  FR=0.378005;
      else if(PTCorr<25)  FR=0.406647;
      else if(PTCorr<35)  FR=0.345562;
      else if(PTCorr<50)  FR=0.284941;
      else                FR=0.406716;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.420458;
      else if(PTCorr<20)  FR=0.210949;
      else if(PTCorr<25)  FR=0.174416;
      else if(PTCorr<35)  FR=0.146076;
      else if(PTCorr<50)  FR=0.113188;
      else                FR=0.097384;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.462609;
      else if(PTCorr<20)  FR=0.256771;
      else if(PTCorr<25)  FR=0.228671;
      else if(PTCorr<35)  FR=0.185761;
      else if(PTCorr<50)  FR=0.167516;
      else                FR=0.142433;
    }
    else{
      if     (PTCorr<15)  FR=0.542713;
      else if(PTCorr<20)  FR=0.330838;
      else if(PTCorr<25)  FR=0.300315;
      else if(PTCorr<35)  FR=0.282436;
      else if(PTCorr<50)  FR=0.234402;
      else                FR=0.210434;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.421107;
      else if(PTCorr<20)  FR=0.214535;
      else if(PTCorr<25)  FR=0.186623;
      else if(PTCorr<35)  FR=0.157792;
      else if(PTCorr<50)  FR=0.124924;
      else                FR=0.110091;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.463734;
      else if(PTCorr<20)  FR=0.260781;
      else if(PTCorr<25)  FR=0.236924;
      else if(PTCorr<35)  FR=0.198074;
      else if(PTCorr<50)  FR=0.173928;
      else                FR=0.15409 ;
    }
    else{
      if     (PTCorr<15)  FR=0.543834;
      else if(PTCorr<20)  FR=0.335507;
      else if(PTCorr<25)  FR=0.303724;
      else if(PTCorr<35)  FR=0.288788;
      else if(PTCorr<50)  FR=0.245314;
      else                FR=0.22206 ;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_TrkIsoVVLConeSUSY"){
//    if(fEta<0.9){//HasB
//      if     (PTCorr<15)  FR=0.412987;
//      else if(PTCorr<20)  FR=0.214109;
//      else if(PTCorr<25)  FR=0.173707;
//      else if(PTCorr<35)  FR=0.142558;
//      else if(PTCorr<50)  FR=0.092292;
//      else                FR=0.087386;
//    }
//    else if(fEta<1.6){
//      if     (PTCorr<15)  FR=0.489415;
//      else if(PTCorr<20)  FR=0.276985;
//      else if(PTCorr<25)  FR=0.236499;
//      else if(PTCorr<35)  FR=0.198413;
//      else if(PTCorr<50)  FR=0.19158 ;
//      else                FR=0.159617;
//    }
//    else{
//      if     (PTCorr<15)  FR=0.568274;
//      else if(PTCorr<20)  FR=0.353558;
//      else if(PTCorr<25)  FR=0.345517;
//      else if(PTCorr<35)  FR=0.303058;
//      else if(PTCorr<50)  FR=0.258099;
//      else                FR=0.231955;
//    }

    if(fEta<0.9){//NearB
      if     (PTCorr<15)  FR=0.297745 ;
      else if(PTCorr<20)  FR=0.15785  ;
      else if(PTCorr<25)  FR=0.123827 ;
      else if(PTCorr<35)  FR=0.100511 ;
      else if(PTCorr<50)  FR=0.06543  ;
      else                FR=0.0693239;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.324575 ;
      else if(PTCorr<20)  FR=0.207702 ;
      else if(PTCorr<25)  FR=0.176821 ;
      else if(PTCorr<35)  FR=0.119094 ;
      else if(PTCorr<50)  FR=0.160515 ;
      else                FR=0.0776659;
    }
    else{
      if     (PTCorr<15)  FR=0.397288;
      else if(PTCorr<20)  FR=0.221381;
      else if(PTCorr<25)  FR=0.232421;
      else if(PTCorr<35)  FR=0.208876;
      else if(PTCorr<50)  FR=0.167898;
      else                FR=0.214372;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1sig4NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.477702;
      else if(PTCorr<20)  FR=0.352577;
      else if(PTCorr<25)  FR=0.336932;
      else if(PTCorr<35)  FR=0.298713;
      else if(PTCorr<50)  FR=0.269651;
      else                FR=0.250289;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.494674;
      else if(PTCorr<20)  FR=0.395798;
      else if(PTCorr<25)  FR=0.380622;
      else if(PTCorr<35)  FR=0.334755;
      else if(PTCorr<50)  FR=0.336507;
      else                FR=0.332443;
    }
    else{
      if     (PTCorr<15)  FR=0.49395 ;
      else if(PTCorr<20)  FR=0.427594;
      else if(PTCorr<25)  FR=0.413441;
      else if(PTCorr<35)  FR=0.411078;
      else if(PTCorr<50)  FR=0.370971;
      else                FR=0.375592;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1sig4NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){//NearB
      if     (PTCorr<15)  FR=0.388324;
      else if(PTCorr<20)  FR=0.245076;
      else if(PTCorr<25)  FR=0.23741 ;
      else if(PTCorr<35)  FR=0.206138;
      else if(PTCorr<50)  FR=0.161407;
      else                FR=0.156359;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.372285;
      else if(PTCorr<20)  FR=0.279095;
      else if(PTCorr<25)  FR=0.279538;
      else if(PTCorr<35)  FR=0.221924;
      else if(PTCorr<50)  FR=0.252549;
      else                FR=0.199088;
    }
    else{
      if     (PTCorr<15)  FR=0.280805;
      else if(PTCorr<20)  FR=0.242246;
      else if(PTCorr<25)  FR=0.306331;
      else if(PTCorr<35)  FR=0.28124 ;
      else if(PTCorr<50)  FR=0.233069;
      else                FR=0.251672;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1sig4NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.47658 ;
      else if(PTCorr<20)  FR=0.344532;
      else if(PTCorr<25)  FR=0.317605;
      else if(PTCorr<35)  FR=0.276683;
      else if(PTCorr<50)  FR=0.240431;
      else                FR=0.225033;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.491739;
      else if(PTCorr<20)  FR=0.383535;
      else if(PTCorr<25)  FR=0.361495;
      else if(PTCorr<35)  FR=0.309984;
      else if(PTCorr<50)  FR=0.316713;
      else                FR=0.300761;
    }
    else{
      if     (PTCorr<15)  FR=0.488566;
      else if(PTCorr<20)  FR=0.408808;
      else if(PTCorr<25)  FR=0.396667;
      else if(PTCorr<35)  FR=0.387424;
      else if(PTCorr<50)  FR=0.343767;
      else                FR=0.348801;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.395574;
      else if(PTCorr<20)  FR=0.207515;
      else if(PTCorr<25)  FR=0.183348;
      else if(PTCorr<35)  FR=0.155668;
      else if(PTCorr<50)  FR=0.123627;
      else                FR=0.107321;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.419747;
      else if(PTCorr<20)  FR=0.244993;
      else if(PTCorr<25)  FR=0.229699;
      else if(PTCorr<35)  FR=0.192354;
      else if(PTCorr<50)  FR=0.170429;
      else                FR=0.14975 ;
    }
    else{
      if     (PTCorr<15)  FR=0.439288;
      else if(PTCorr<20)  FR=0.289097;
      else if(PTCorr<25)  FR=0.271926;
      else if(PTCorr<35)  FR=0.26034 ;
      else if(PTCorr<50)  FR=0.223463;
      else                FR=0.202654;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){//NearB
      if     (PTCorr<15)  FR=0.299092 ;
      else if(PTCorr<20)  FR=0.143535 ;
      else if(PTCorr<25)  FR=0.117148 ;
      else if(PTCorr<35)  FR=0.103306 ;
      else if(PTCorr<50)  FR=0.0765811;
      else                FR=0.0630313;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.304042;
      else if(PTCorr<20)  FR=0.176138;
      else if(PTCorr<25)  FR=0.164007;
      else if(PTCorr<35)  FR=0.123931;
      else if(PTCorr<50)  FR=0.131938;
      else                FR=0.093126;
    }
    else{
      if     (PTCorr<15)  FR=0.23905 ;
      else if(PTCorr<20)  FR=0.16028 ;
      else if(PTCorr<25)  FR=0.19258 ;
      else if(PTCorr<35)  FR=0.176143;
      else if(PTCorr<50)  FR=0.137479;
      else                FR=0.130829;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.396487;
      else if(PTCorr<20)  FR=0.207776;
      else if(PTCorr<25)  FR=0.183522;
      else if(PTCorr<35)  FR=0.155841;
      else if(PTCorr<50)  FR=0.123879;
      else                FR=0.107774;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.421205;
      else if(PTCorr<20)  FR=0.245923;
      else if(PTCorr<25)  FR=0.230051;
      else if(PTCorr<35)  FR=0.192691;
      else if(PTCorr<50)  FR=0.170713;
      else                FR=0.150102;
    }
    else{
      if     (PTCorr<15)  FR=0.440847;
      else if(PTCorr<20)  FR=0.289948;
      else if(PTCorr<25)  FR=0.272952;
      else if(PTCorr<35)  FR=0.261306;
      else if(PTCorr<50)  FR=0.224457;
      else                FR=0.20334 ;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.300225 ;
      else if(PTCorr<20)  FR=0.143649 ;
      else if(PTCorr<25)  FR=0.117231 ;
      else if(PTCorr<35)  FR=0.103437 ;
      else if(PTCorr<50)  FR=0.0767057;
      else                FR=0.0631149;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.304083 ;
      else if(PTCorr<20)  FR=0.176359 ;
      else if(PTCorr<25)  FR=0.164183 ;
      else if(PTCorr<35)  FR=0.12428  ;
      else if(PTCorr<50)  FR=0.13231  ;
      else                FR=0.0932271;
    }
    else{
      if     (PTCorr<15)  FR=0.239239;
      else if(PTCorr<20)  FR=0.160474;
      else if(PTCorr<25)  FR=0.192764;
      else if(PTCorr<35)  FR=0.176604;
      else if(PTCorr<50)  FR=0.137952;
      else                FR=0.131124;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.395195 ;
      else if(PTCorr<20)  FR=0.202967 ;
      else if(PTCorr<25)  FR=0.169647 ;
      else if(PTCorr<35)  FR=0.143029 ;
      else if(PTCorr<50)  FR=0.111431 ;
      else                FR=0.0951533;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.418279;
      else if(PTCorr<20)  FR=0.238722;
      else if(PTCorr<25)  FR=0.217297;
      else if(PTCorr<35)  FR=0.177343;
      else if(PTCorr<50)  FR=0.161958;
      else                FR=0.136925;
    }
    else{
      if     (PTCorr<15)  FR=0.435522;
      else if(PTCorr<20)  FR=0.276537;
      else if(PTCorr<25)  FR=0.259684;
      else if(PTCorr<35)  FR=0.245707;
      else if(PTCorr<50)  FR=0.207092;
      else                FR=0.187308;
    }
  }
  else if(Option=="QCD_BjMatch_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.212887;
      else if(PTCorr<20)  FR=0.170592;
      else if(PTCorr<25)  FR=0.16553 ;
      else if(PTCorr<35)  FR=0.153397;
      else if(PTCorr<50)  FR=0.123731;
      else                FR=0.104981;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.239631;
      else if(PTCorr<20)  FR=0.206327;
      else if(PTCorr<25)  FR=0.216909;
      else if(PTCorr<35)  FR=0.182379;
      else if(PTCorr<50)  FR=0.170383;
      else                FR=0.156924;
    }
    else{
      if     (PTCorr<15)  FR=0.22021 ;
      else if(PTCorr<20)  FR=0.228938;
      else if(PTCorr<25)  FR=0.24682 ;
      else if(PTCorr<35)  FR=0.241839;
      else if(PTCorr<50)  FR=0.204816;
      else                FR=0.19224 ;
    }
  }
  else if(Option=="QCD_CjMatch_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.280125;
      else if(PTCorr<20)  FR=0.206693;
      else if(PTCorr<25)  FR=0.237403;
      else if(PTCorr<35)  FR=0.232245;
      else if(PTCorr<50)  FR=0.222692;
      else                FR=0.233209;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.285224;
      else if(PTCorr<20)  FR=0.263604;
      else if(PTCorr<25)  FR=0.263886;
      else if(PTCorr<35)  FR=0.261374;
      else if(PTCorr<50)  FR=0.276096;
      else                FR=0.250051;
    }
    else{
      if     (PTCorr<15)  FR=0.316214;
      else if(PTCorr<20)  FR=0.275994;
      else if(PTCorr<25)  FR=0.315642;
      else if(PTCorr<35)  FR=0.333266;
      else if(PTCorr<50)  FR=0.308791;
      else                FR=0.331183;
    }
  }
  else if(Option=="QCD_LjMatch_POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.235959;
      else if(PTCorr<20)  FR=0.163328;
      else if(PTCorr<25)  FR=0.255066;
      else if(PTCorr<35)  FR=0.220489;
      else if(PTCorr<50)  FR=0.150643;
      else                FR=0.217929;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.290873;
      else if(PTCorr<20)  FR=0.200495;
      else if(PTCorr<25)  FR=0.300621;
      else if(PTCorr<35)  FR=0.262121;
      else if(PTCorr<50)  FR=0.316673;
      else                FR=0.371495;
    }
    else{
      if     (PTCorr<15)  FR=0.389803;
      else if(PTCorr<20)  FR=0.364795;
      else if(PTCorr<25)  FR=0.398313;
      else if(PTCorr<35)  FR=0.382796;
      else if(PTCorr<50)  FR=0.383164;
      else                FR=0.421589;
    }
  }

  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.378176 ;
      else if(PTCorr<20)  FR=0.184572 ;
      else if(PTCorr<25)  FR=0.160799 ;
      else if(PTCorr<35)  FR=0.135949 ;
      else if(PTCorr<50)  FR=0.107135 ;
      else                FR=0.0933122;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.398919;
      else if(PTCorr<20)  FR=0.213567;
      else if(PTCorr<25)  FR=0.196568;
      else if(PTCorr<35)  FR=0.164957;
      else if(PTCorr<50)  FR=0.145089;
      else                FR=0.130676;
    }
    else{
      if     (PTCorr<15)  FR=0.415817;
      else if(PTCorr<20)  FR=0.250131;
      else if(PTCorr<25)  FR=0.231834;
      else if(PTCorr<35)  FR=0.221139;
      else if(PTCorr<50)  FR=0.187845;
      else                FR=0.173545;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.28293  ;
      else if(PTCorr<20)  FR=0.13529  ;
      else if(PTCorr<25)  FR=0.106368 ;
      else if(PTCorr<35)  FR=0.0874014;
      else if(PTCorr<50)  FR=0.0660323;
      else                FR=0.0541946;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.276023 ;
      else if(PTCorr<20)  FR=0.154409 ;
      else if(PTCorr<25)  FR=0.140863 ;
      else if(PTCorr<35)  FR=0.110194 ;
      else if(PTCorr<50)  FR=0.105941 ;
      else                FR=0.0796622;
    }
    else{
      if     (PTCorr<15)  FR=0.255897;
      else if(PTCorr<20)  FR=0.14699 ;
      else if(PTCorr<25)  FR=0.163857;
      else if(PTCorr<35)  FR=0.149301;
      else if(PTCorr<50)  FR=0.117876;
      else                FR=0.111363;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.376251 ;
      else if(PTCorr<20)  FR=0.180005 ;
      else if(PTCorr<25)  FR=0.148464 ;
      else if(PTCorr<35)  FR=0.122854 ;
      else if(PTCorr<50)  FR=0.0952236;
      else                FR=0.081313 ;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.394537;
      else if(PTCorr<20)  FR=0.206124;
      else if(PTCorr<25)  FR=0.18436 ;
      else if(PTCorr<35)  FR=0.151379;
      else if(PTCorr<50)  FR=0.135037;
      else                FR=0.117474;
    }
    else{
      if     (PTCorr<15)  FR=0.409673;
      else if(PTCorr<20)  FR=0.23768 ;
      else if(PTCorr<25)  FR=0.21949 ;
      else if(PTCorr<35)  FR=0.206409;
      else if(PTCorr<50)  FR=0.172383;
      else                FR=0.158458;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.264291;
      else if(PTCorr<20)  FR=0.133615;
      else if(PTCorr<25)  FR=0.128545;
      else if(PTCorr<35)  FR=0.10907 ;
      else if(PTCorr<50)  FR=0.085681;
      else                FR=0.073323;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.287358;
      else if(PTCorr<20)  FR=0.163346;
      else if(PTCorr<25)  FR=0.160179;
      else if(PTCorr<35)  FR=0.135796;
      else if(PTCorr<50)  FR=0.117967;
      else                FR=0.104968;
    }
    else{
      if     (PTCorr<15)  FR=0.328588;
      else if(PTCorr<20)  FR=0.206919;
      else if(PTCorr<25)  FR=0.195996;
      else if(PTCorr<35)  FR=0.185807;
      else if(PTCorr<50)  FR=0.155233;
      else                FR=0.140276;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.0846029;
      else if(PTCorr<20)  FR=0.0465758;
      else if(PTCorr<25)  FR=0.041396 ;
      else if(PTCorr<35)  FR=0.0355063;
      else if(PTCorr<50)  FR=0.0262166;
      else                FR=0.0210795;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.0993301;
      else if(PTCorr<20)  FR=0.058441 ;
      else if(PTCorr<25)  FR=0.0579149;
      else if(PTCorr<35)  FR=0.0449768;
      else if(PTCorr<50)  FR=0.041935 ;
      else                FR=0.0301986;
    }
    else{
      if     (PTCorr<15)  FR=0.0943277;
      else if(PTCorr<20)  FR=0.0611011;
      else if(PTCorr<25)  FR=0.0708035;
      else if(PTCorr<35)  FR=0.0657832;
      else if(PTCorr<50)  FR=0.0512415;
      else                FR=0.0505063;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.25635  ;
      else if(PTCorr<20)  FR=0.118343 ;
      else if(PTCorr<25)  FR=0.0960862;
      else if(PTCorr<35)  FR=0.0781842;
      else if(PTCorr<50)  FR=0.0592189;
      else                FR=0.0490953;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.274516 ;
      else if(PTCorr<20)  FR=0.139824 ;
      else if(PTCorr<25)  FR=0.12374  ;
      else if(PTCorr<35)  FR=0.0997158;
      else if(PTCorr<50)  FR=0.0866521;
      else                FR=0.0734869;
    }
    else{
      if     (PTCorr<15)  FR=0.31102 ;
      else if(PTCorr<20)  FR=0.175795;
      else if(PTCorr<25)  FR=0.158154;
      else if(PTCorr<35)  FR=0.146317;
      else if(PTCorr<50)  FR=0.119155;
      else                FR=0.107794;
    }
  }
  else if(Option=="QCD_AwayB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.402142 ;
      else if(PTCorr<20)  FR=0.190637 ;
      else if(PTCorr<25)  FR=0.163378 ;
      else if(PTCorr<35)  FR=0.137635 ;
      else if(PTCorr<50)  FR=0.108139 ;
      else                FR=0.0951788;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.439664;
      else if(PTCorr<20)  FR=0.225892;
      else if(PTCorr<25)  FR=0.202145;
      else if(PTCorr<35)  FR=0.169413;
      else if(PTCorr<50)  FR=0.147694;
      else                FR=0.134003;
    }
    else{
      if     (PTCorr<15)  FR=0.51159 ;
      else if(PTCorr<20)  FR=0.289179;
      else if(PTCorr<25)  FR=0.256916;
      else if(PTCorr<35)  FR=0.243065;
      else if(PTCorr<50)  FR=0.204101;
      else                FR=0.188659;
    }
  }
  else if(Option=="QCD_HasB_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.351364 ;
      else if(PTCorr<20)  FR=0.153304 ;
      else if(PTCorr<25)  FR=0.114185 ;
      else if(PTCorr<35)  FR=0.0917333;
      else if(PTCorr<50)  FR=0.068622 ;
      else                FR=0.0556129;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.379511;
      else if(PTCorr<20)  FR=0.192976;
      else if(PTCorr<25)  FR=0.163234;
      else if(PTCorr<35)  FR=0.122893;
      else if(PTCorr<50)  FR=0.115132;
      else                FR=0.086039;
    }
    else{
      if     (PTCorr<15)  FR=0.463989;
      else if(PTCorr<20)  FR=0.238416;
      else if(PTCorr<25)  FR=0.232298;
      else if(PTCorr<35)  FR=0.203586;
      else if(PTCorr<50)  FR=0.1504  ;
      else                FR=0.136749;
    }
  }
  else if(Option=="QCD_POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.401259 ;
      else if(PTCorr<20)  FR=0.187452 ;
      else if(PTCorr<25)  FR=0.152679 ;
      else if(PTCorr<35)  FR=0.125571 ;
      else if(PTCorr<50)  FR=0.0968874;
      else                FR=0.0831219;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.437923;
      else if(PTCorr<20)  FR=0.222317;
      else if(PTCorr<25)  FR=0.194375;
      else if(PTCorr<35)  FR=0.158588;
      else if(PTCorr<50)  FR=0.139727;
      else                FR=0.122123;
    }
    else{
      if     (PTCorr<15)  FR=0.510331;
      else if(PTCorr<20)  FR=0.284652;
      else if(PTCorr<25)  FR=0.253272;
      else if(PTCorr<35)  FR=0.236261;
      else if(PTCorr<50)  FR=0.193669;
      else                FR=0.177171;
    }
  }
  else { cout<<"RateMC Error"<<endl; }

  return FR;
}



int Aug2017_MuFakeMCStudy::GetFakeLepSrcIdx(snu::KMuon Mu, std::vector<snu::KTruth> TruthColl){

  int MatchedIdx=-1;
  int LepType=GetLeptonType(Mu,TruthColl);
  if(LepType==-1){
    float mindR=999., maxdR=0.4; int IdxMindR=-1;
    for(int j=2; j<(int) TruthColl.size(); j++){
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

int Aug2017_MuFakeMCStudy::GetFakeLepGenJetSrcIdx(snu::KMuon Mu, std::vector<snu::KGenJet> GenJetColl, TString Option){

  int SrcIdx=-1;
  float dRmin=999.;
  bool DR04=Option.Contains("DR04"), Nearest=Option.Contains("Near");
  for(int i=0; i<(int) GenJetColl.size(); i++){
    if     ( DR04  ){ if(Mu.DeltaR(GenJetColl.at(i))<0.4){ SrcIdx=i; break; } }
    else if(Nearest){ if(Mu.DeltaR(GenJetColl.at(i))<dRmin){ dRmin=Mu.DeltaR(GenJetColl.at(i)); SrcIdx=i; } }
  }

  return SrcIdx;
}


int Aug2017_MuFakeMCStudy::GetFakeLepJetSrcIdx(snu::KMuon Mu, std::vector<snu::KJet> JetColl, TString Option){

  int SrcIdx=-1;
  float dRmin=999.;
  bool DR04=Option.Contains("DR04"), Nearest=Option.Contains("Near");
  for(int i=0; i<(int) JetColl.size(); i++){
    if     ( DR04  ){ if(Mu.DeltaR(JetColl.at(i))<0.4){ SrcIdx=i; break; } }
    else if(Nearest){ if(Mu.DeltaR(JetColl.at(i))<dRmin){ dRmin=Mu.DeltaR(JetColl.at(i)); SrcIdx=i; } }
  }

  return SrcIdx;
}


int Aug2017_MuFakeMCStudy::GetFakeLepJetSrcType(snu::KMuon Mu, std::vector<snu::KJet> JetColl){
  //Type: -1: Unmatched, 1:L, 2:C, 3:B
  int SrcType=-1;
  bool NearB=false, NearC=false, NearL=false;
  for(int i=0; i<(int) JetColl.size(); i++){
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

   for(int i=0; i<(int) MuColl.size(); i++){
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

   for(int i=0; i<(int) EleColl.size(); i++){
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

  for(int i=2; i<(int) TruthColl.size(); i++){

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
