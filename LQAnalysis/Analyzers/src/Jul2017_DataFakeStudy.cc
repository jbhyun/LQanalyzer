// $Id: Jul2017_DataFakeStudy.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQJul2017_DataFakeStudy Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
 #include "Jul2017_DataFakeStudy.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Jul2017_DataFakeStudy);

 Jul2017_DataFakeStudy::Jul2017_DataFakeStudy() : AnalyzerCore(), out_muons(0) {

   SetLogName("Jul2017_DataFakeStudy");
   Message("In Jul2017_DataFakeStudy constructor", INFO);
   InitialiseAnalysis();
 }


 void Jul2017_DataFakeStudy::InitialiseAnalysis() throw( LQError ) {
   
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

void Jul2017_DataFakeStudy::ExecuteEvents()throw( LQError ){

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
   float pileup_reweight=1., pileup_reweight_systdown=1., pileup_reweight_systup=1.;
   if(!k_isdata){ pileup_reweight          = eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central); 
                  pileup_reweight_systup   = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::up)  /pileup_reweight :0. ;
                  pileup_reweight_systdown = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down)/pileup_reweight :0. ;
   }

   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight*pileup_reweight, 0., 50., 50);


   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);


   bool EleFR=false, TrigSel=false, NormCheck=false, UnPreTrig=false, PreTrig=false, SiglPreTrig=false, MultPreTrig=false, Closure=false;
   bool FRScan=false, IDValidation=false;
   bool HighdXYCheck=false, IsoIPOpt=false;
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("EleFR"))        EleFR        = true;
     else if(k_flags.at(i).Contains("TrigSel"))      TrigSel      = true;
     else if(k_flags.at(i).Contains("NormCheck"))    NormCheck    = true;
     else if(k_flags.at(i).Contains("UnPreTrig"))    UnPreTrig    = true;
     else if(k_flags.at(i).Contains("SiglPreTrig"))  SiglPreTrig  = true;
     else if(k_flags.at(i).Contains("MultPreTrig"))  MultPreTrig  = true;
     else if(k_flags.at(i).Contains("Closure"))      Closure      = true;
     else if(k_flags.at(i).Contains("FRScan"))       FRScan       = true;
     else if(k_flags.at(i).Contains("HighdXYCheck")) HighdXYCheck = true;
     else if(k_flags.at(i).Contains("IsoIPOpt"))     IsoIPOpt     = true;
     else if(k_flags.at(i).Contains("IDValidation")) IDValidation = true;
   }

   if(SiglPreTrig || MultPreTrig) PreTrig=true;
    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //ListTriggersAvailable();

   //Trigger Path of Analysis
   bool Pass_Trigger=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   if(EleFR || NormCheck || FRScan || HighdXYCheck){
     if(UnPreTrig){
       if( PassTrigger("HLT_Ele27_WPTight_Gsf_v") ) Pass_Trigger=true;
       if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Ele27_WPTight_Gsf_v", TargetLumi);
     }
     if(SiglPreTrig){
       if( PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) Pass_Trigger=true;
       if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
     }
     if(MultPreTrig){
//         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
//         eventbase->GetElectronSel()->SetBETrRegIncl(false);
//       std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
//       std::vector<snu::KElectron> electronFakeLColl;
//         for(int i=0; i<electronPreColl.size(); i++){
//           if(PassIDCriteria(electronPreColl.at(i), "POGMVAMFakeLIso04Opt2")) electronFakeLColl.push_back(electronPreColl.at(i));
//         }

         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
         eventbase->GetElectronSel()->SetdxySigMax(3.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> electronFakeLColl; eventbase->GetElectronSel()->Selection(electronFakeLColl);


       int Case=0;
       if( electronFakeLColl.size()==1 && electronFakeLColl.at(0).Pt()>25
           && PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) {Pass_Trigger=true; Case=1;}
       if( electronFakeLColl.size()==1 && electronFakeLColl.at(0).Pt()<25 && electronFakeLColl.at(0).Pt()>20 
           && PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) {Pass_Trigger=true; Case=2;}
       if( electronFakeLColl.size()==1 && electronFakeLColl.at(0).Pt()<20 && electronFakeLColl.at(0).Pt()>15 
           && PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") ) {Pass_Trigger=true; Case=3;}


       if     (!isData && Case==1) trigger_ps_weight=WeightByTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
       else if(!isData && Case==2) trigger_ps_weight=WeightByTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
       else if(!isData && Case==3) trigger_ps_weight=WeightByTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v", TargetLumi);
     }

   }
   if(Closure){
     int Pass_Trigger1=0, Pass_Trigger2=0;
     if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") ) Pass_Trigger1++;
     if( PassTrigger("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger2++;

     if(isData){
       int DataPeriod=GetDataPeriod();
       if( DataPeriod>0 && DataPeriod<7 && Pass_Trigger1==1 ) Pass_Trigger=true;
       else if( DataPeriod==7 && Pass_Trigger2==1 ) Pass_Trigger=true;
     }
     else{
       if( Pass_Trigger1>0 || Pass_Trigger2>0 ) Pass_Trigger=true;
       trigger_period_weight=(Pass_Trigger1*27.257618+Pass_Trigger2*8.605696)/35.863314;
       trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
     }
   }
   if(IsoIPOpt || IDValidation){
     if( PassTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_Trigger=true;
     if(!isData) trigger_ps_weight=WeightByTrigger("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", TargetLumi);
   }
   FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 1., 100);
   weight*=trigger_ps_weight*trigger_period_weight;


     //HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v   6.992 0.1M
     //HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v          6.162 0.2M
     //HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v 14.888 0.25M
     //HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v         30.397 1M
     //HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v 58.896 0.9M
     //HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v         16.43
     //HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v 63.046

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
   //std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  
   /**PreSelCut***********************************************************************************************/
   //Intended for Code speed boosting up.
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_LOOSE);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
     if(EleFR || TrigSel || NormCheck || FRScan || HighdXYCheck ){ if( !(electronPreColl.size()>=1) ) return; }
     else if(Closure)                   { if( !(electronPreColl.size()>=1 && muonPreColl.size()>=2) ) return; }
     else if(IsoIPOpt)                  { if( !(electronPreColl.size()>=2) ) return; }
     FillCutFlow("PreSel", weight);
   /**********************************************************************************************************/

     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.4);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");  eventbase->GetMuonSel()->SetRelIso(0.1);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);
   std::vector<snu::KMuon> muonColl;      if(k_running_nonprompt){ muonColl=muonLooseColl;} else{ muonColl=muonTightColl;}


     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronFakeLPreOptColl; eventbase->GetElectronSel()->Selection(electronFakeLPreOptColl);


     std::vector<snu::KElectron> electronFakeL2Coll;
     //WP90
//       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
//         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMTIsop1IPp025p05Sig4FakeLIsop4")) electronFakeL2Coll.push_back(electronFakeLPreOptColl.at(i));
//       }
//       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
//         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVAMTIsop06IPp025p05Sig4FakeLIsop4")) electronFakeL2Coll.push_back(electronFakeLPreOptColl.at(i));
//       }
//       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
//         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVATTIsop1IPp025p05Sig4FakeLIsop4")) electronFakeL2Coll.push_back(electronFakeLPreOptColl.at(i));
//       }
//       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
//         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "POGMVATTIsop06IPp025p05Sig4FakeLIsop4")) electronFakeL2Coll.push_back(electronFakeLPreOptColl.at(i));
//       }
       for(int i=0; i<electronFakeLPreOptColl.size(); i++){
         if(PassIDCriteria(electronFakeLPreOptColl.at(i), "HctoWAFakeLoose")) electronFakeL2Coll.push_back(electronFakeLPreOptColl.at(i));
       }
     


     eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
     eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.05, 0.05);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
   std::vector<snu::KElectron> electronColl;  if(!k_running_nonprompt) electronColl=electronTightColl;

     bool LeptonVeto=false;
     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
   std::vector<snu::KJet> jetColl; eventbase->GetJetSel()->Selection(jetColl, LeptonVeto, muonLooseColl, electronFakeLPreOptColl);



   //Method to apply 2a method SF
   //std::vector<int> bIdxColl  = GetSFBJetIdx(jetColl,"Medium");
   //std::vector<int> ljIdxColl = GetSFLJetIdx(jetColl, bIdxColl, "Medium");
   //std::vector<snu::KJet> bjetColl2a; for(int i=0; i<bIdxColl.size(); i++) {bjetColl2a.push_back(jetColl.at(bIdxColl.at(i)));}
   //std::vector<snu::KJet> ljetColl2a; for(int i=0; i<ljIdxColl.size(); i++){ljetColl2a.push_back(jetColl.at(ljIdxColl.at(i)));}

   std::vector<snu::KJet> bjetColl = SelBJets(jetColl, "Medium");
   std::vector<snu::KJet> ljetColl = SelLightJets(jetColl, "Medium");

   float PFMETType1    = eventbase->GetEvent().PFMETType1();
   float met           = eventbase->GetEvent().MET(snu::KEvent::pfmet);
   float met_JESup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::up);
   float met_JERup     = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::up);
   float met_Unclup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::up);
   float met_ElEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::up);
   float met_MuEnup    = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::up);
   float met_JESdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetEn       , snu::KEvent::down);
   float met_JERdown   = eventbase->GetEvent().PFMETShifted(snu::KEvent::JetRes      , snu::KEvent::down);
   float met_Uncldown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::Unclustered , snu::KEvent::down);
   float met_ElEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::ElectronEn  , snu::KEvent::down);
   float met_MuEndown  = eventbase->GetEvent().PFMETShifted(snu::KEvent::MuonEn      , snu::KEvent::down);

   float metphi = eventbase->GetEvent().METPhi();
   double met_x = eventbase->GetEvent().PFMETx();
   double met_y = eventbase->GetEvent().PFMETy();
   double Pzv,Pzv1, Pzv2;
   snu::KParticle v; v.SetPxPyPzE(met_x, met_y, 0, sqrt(met_x*met_x+met_y*met_y));
   //snu::KParticle v[4]; v[0].SetPx(met_x); v[0].SetPy(met_y);
   //                     v[1].SetPx(met_x); v[1].SetPy(met_y);

   int nbjets=bjetColl.size(); int njets=jetColl.size(); const int nljets=ljetColl.size();//njets-nbjets;//number of light jets
   int Nvtx=eventbase->GetEvent().nVertices();


   /*****************************************************
   **Scale Factors
   *****************************************************/
   float id_weight_ele=1., reco_weight_ele=1., trk_weight_mu=1., id_weight_mu=1., iso_weight_mu=1., btag_sf=1.;
   float trigger_sf=1.;
   float fake_weight=1.; bool EventCand=false;
   float geneff_weight=1., gennorm_weight=1.;
//
   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/

   if(true){
     if(!isData){
       if(NormCheck || Closure || IDValidation){
       //trigger_sf      = mcdata_correction->TriggerScaleFactor( electronColl, muonColl, "HLT_IsoMu24_v" );
  
       //id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_80", electronColl);
       id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_MVA_90", electronColl);
       //id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_POG_TIGHT", electronColl);
       reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
  
       id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HN_TRI_TIGHT", muonColl);
       //iso_weight_mu   = mcdata_correction->MuonISOScaleFactor("MUON_POG_TIGHT", muonColl);
       trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
       //btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

       //geneff_weight   = GenFilterEfficiency(k_sample_name);
       //gennorm_weight  = SignalNorm(k_sample_name, 200.);
       }
       if(FRScan){
         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
       }
       if(HighdXYCheck && (k_sample_name.Contains("qcd") ||k_sample_name.Contains("QCD"))) pileup_reweight=1.;
     }
     else{
       //Perfect Prompt Ratio Approximation applied.(p=1), Applicable to generic number, combination of leptons under premise of high prompt rate.
   //    if(k_running_nonprompt){
   //      fake_weight = m_datadriven_bkg->Get_DataDrivenWeight(false, muonLooseColl, "MUON_HN_TRI_TIGHT", muonLooseColl.size(), electronLooseColl, "ELECTRON_HN_LOWDXY_TIGHT", electronLooseColl.size());
//       }
     }
   }

   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*pileup_reweight*fake_weight*btag_sf*geneff_weight*gennorm_weight;
   //-----------------------------------------------------------------------------------------//


//----------------------------------------------------------------------------------//
//////////////////////////////////////////////////////////////////////////////////////
/////// MAIN ANALYSIS CODE ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//----------------------------------------------------------------------------------//

   if(TrigSel){
     if(electronColl.size()!=1) return;
     if(electronColl.at(0).Pt()<25) return;
     float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());

     if(PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v")){
       FillHist("Count_e17", 0., weight, 0., 10., 10);//1: Ele25

       if(jetColl.size()==1){
         if( jetColl.at(0).Pt()>30 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17", 1., weight, 0., 10., 10);//2: 1j away 30
           if( met<30 && MTW<30 ) FillHist("Count_e17", 2., weight, 0., 10., 10);//FakeEnrich 1j30
         }
         if( jetColl.at(0).Pt()>40 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17", 3., weight, 0., 10., 10);//2: 1j away 40
           if( met<30 && MTW<30 ) FillHist("Count_e17", 4., weight, 0., 10., 10);//FakeEnrich 1j40
         }
         if( jetColl.at(0).Pt()>60 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17", 5., weight, 0., 10., 10);//2: 1j away 60
           if( met<30 && MTW<30 ) FillHist("Count_e17", 6., weight, 0., 10., 10);//FakeEnrich 1j60
         }
       }//1j 

     }//End of Ele17 
     if(PassTrigger("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v")){
       FillHist("Count_e17j30", 0., weight, 0., 10., 10);//1: Ele25

       if(jetColl.size()==1){
         if( jetColl.at(0).Pt()>30 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17j30", 1., weight, 0., 10., 10);//2: 1j away 30
           if( met<30 && MTW<30 ) FillHist("Count_e17j30", 2., weight, 0., 10., 10);//FakeEnrich 1j30
         }
         if( jetColl.at(0).Pt()>40 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17j30", 3., weight, 0., 10., 10);//2: 1j away 40
           if( met<30 && MTW<30 ) FillHist("Count_e17j30", 4., weight, 0., 10., 10);//FakeEnrich 1j40
         }
         if( jetColl.at(0).Pt()>60 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e17j30", 5., weight, 0., 10., 10);//2: 1j away 60
           if( met<30 && MTW<30 ) FillHist("Count_e17j30", 6., weight, 0., 10., 10);//FakeEnrich 1j60
         }
       }//1j 

     }//End of Ele17Jet30 
     if(PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v")){
       FillHist("Count_e23", 0., weight, 0., 10., 10);//1: Ele25

       if(jetColl.size()==1){
         if( jetColl.at(0).Pt()>30 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23", 1., weight, 0., 10., 10);//2: 1j away 30
           if( met<30 && MTW<30 ) FillHist("Count_e23", 2., weight, 0., 10., 10);//FakeEnrich 1j30
         }
         if( jetColl.at(0).Pt()>40 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23", 3., weight, 0., 10., 10);//2: 1j away 40
           if( met<30 && MTW<30 ) FillHist("Count_e23", 4., weight, 0., 10., 10);//FakeEnrich 1j40
         }
         if( jetColl.at(0).Pt()>60 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23", 5., weight, 0., 10., 10);//2: 1j away 60
           if( met<30 && MTW<30 ) FillHist("Count_e23", 6., weight, 0., 10., 10);//FakeEnrich 1j60
         }
       }//1j 

     }//End of Ele23 
     if(PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v")){
       FillHist("Count_e23j30", 0., weight, 0., 10., 10);//1: Ele25

       if(jetColl.size()==1){
         if( jetColl.at(0).Pt()>30 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23j30", 1., weight, 0., 10., 10);//2: 1j away 30
           if( met<30 && MTW<30 ) FillHist("Count_e23j30", 2., weight, 0., 10., 10);//FakeEnrich 1j30
         }
         if( jetColl.at(0).Pt()>40 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23j30", 3., weight, 0., 10., 10);//2: 1j away 40
           if( met<30 && MTW<30 ) FillHist("Count_e23j30", 4., weight, 0., 10., 10);//FakeEnrich 1j40
         }
         if( jetColl.at(0).Pt()>60 && jetColl.at(0).DeltaR(electronColl.at(0))>1.0 ){
           FillHist("Count_e23j30", 5., weight, 0., 10., 10);//2: 1j away 60
           if( met<30 && MTW<30 ) FillHist("Count_e23j30", 6., weight, 0., 10., 10);//FakeEnrich 1j60
         }
       }//1j 

     }//End of Ele23Jet30

   //170602 - Conclusion -Ele23Jet30 is the best.
   //though for single lepton events other triggers as Ele17 may give better events, but for dijet topology lepton&jet events
   //Ele23 Jet30 gives the most events. So I chose this trigger.
    
   }//End of TriggerChoice

   if(NormCheck){
     //1) High MTW region(MTW>70, MET>50) mostly includes W related events, no fakes
     //   -a) For Norm check with 2)
     //   -b) For NLO xsec check.
     //2) DY+Jets events ; very clear for prescale check
     
     bool SystRun=true;
     if(SystRun){
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESUpColl; eventbase->GetJetSel()->Selection(jetJESUpColl, LeptonVeto, muonLooseColl, electronFakeLPreOptColl, "SystUpJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERUpColl; eventbase->GetJetSel()->Selection(jetJERUpColl, LeptonVeto, muonLooseColl, electronFakeLPreOptColl, "SystUpJER");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESDownColl; eventbase->GetJetSel()->Selection(jetJESDownColl, LeptonVeto, muonLooseColl, electronFakeLPreOptColl, "SystDownJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERDownColl; eventbase->GetJetSel()->Selection(jetJERDownColl, LeptonVeto, muonLooseColl, electronFakeLPreOptColl, "SystDownJER");



       DoSystRun("NormCheck", "PreTrig", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight);
       if(!isData){
         DoSystRun("NormCheck", "PreTrigSystUpPU", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight*pileup_reweight_systup);
         DoSystRun("NormCheck", "PreTrigSystUpUncl", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi), weight);
         DoSystRun("NormCheck", "PreTrigSystUpJES", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetJESUpColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi), weight);
         DoSystRun("NormCheck", "PreTrigSystUpJER", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetJERUpColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi), weight);
         DoSystRun("NormCheck", "PreTrigSystDownPU", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight*pileup_reweight_systdown);
         DoSystRun("NormCheck", "PreTrigSystDownUncl", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi), weight);
         DoSystRun("NormCheck", "PreTrigSystDownJES", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetJESDownColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi), weight);
         DoSystRun("NormCheck", "PreTrigSystDownJER", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetJERDownColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi), weight);
       }

     }
     if(!SystRun){
     
     float nvtx_reweight=1., purw_down=1.;      
     std::vector<snu::KJet> jetVeto2Coll = SkimJetColl(jetColl, electronFakeL2Coll, muonLooseColl, "EleMuVeto");

     if(!isData){
        nvtx_reweight=mcdata_correction->UserPileupWeight(eventbase->GetEvent(), jetVeto2Coll.size()); //Measured in DY UnPreTrig
        //nvtx_reweight=GetPreTrigPURW(Nvtx);//Measured in PreTrig
        if(pileup_reweight!=0) purw_down=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down)/pileup_reweight;
     }


     //Cut1| 1: 1l(Pass Trigger turn-on) / 2: 2l(Pass Trigger turn-on) / 0 : not 1, -1 case (0l or more than 2l)
     //   2| 0: 0j / 1: 1j / 2: geq2j (>0 means geq1j)
     //   3| 1: PT(j1)>40 / 0: fail cut (requires Cut2>0)
     //   4| 1: MET>30 / 0: fail cut
     //   5| 1: MTW>50 / 0: fail cut (implicitly requires Cut1==1)


     if(electronColl.size()==1 && electronFakeL2Coll.size()==1 && muonLooseColl.size()==0){
       if( UnPreTrig && electronColl.at(0).Pt()<30) return;
       if( SiglPreTrig   && electronColl.at(0).Pt()<25) return;
       float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());

       //Nvtx
       if(met>30){
         FillHist("Nvtx_met30", Nvtx, weight, 0., 50., 50);
         FillHist("Nvtx_met30_PUdown", Nvtx, weight*purw_down, 0., 50., 50);
         FillHist("Nvtx_met30_NvtxRW", Nvtx, weight*nvtx_reweight, 0., 50., 50);
         if(jetVeto2Coll.size()==0) {
           FillHist("Nvtx_0jmet30", Nvtx, weight, 0., 50., 50);
           FillHist("Nvtx_0jmet30_NvtxRW", Nvtx, weight*nvtx_reweight, 0., 50., 50);
           FillHist("NvtxRW_0jmet30", Nvtx, nvtx_reweight, 0., 50., 50);
         }
         if(jetVeto2Coll.size()>0 && jetVeto2Coll.at(0).Pt()>40){
           if(jetVeto2Coll.size()==1) {
             FillHist("Nvtx_1jmet30", Nvtx, weight, 0., 50., 50);
             FillHist("Nvtx_1jmet30_NvtxRW", Nvtx, weight*nvtx_reweight, 0., 50., 50);
             FillHist("NvtxRW_1jmet30", Nvtx, nvtx_reweight, 0., 50., 50);
           }
           FillHist("Nvtx_gt1jmet30", Nvtx, weight, 0., 50., 50);
           FillHist("Nvtx_gt1jmet30_PUdown", Nvtx, weight*purw_down, 0., 50., 50);
           FillHist("Nvtx_gt1jmet30_NvtxRW", Nvtx, weight*nvtx_reweight, 0., 50., 50);
         }         
       }

       //For Validation in unprescaled path
       FillHist("PTe", electronColl.at(0).Pt(), weight, 0., 200., 200);
       FillHist("Etae", electronColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Nj", jetVeto2Coll.size(), weight, 0., 10., 10);
       if(jetVeto2Coll.size()>0){
         FillHist("PTj1", jetVeto2Coll.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etaj1", jetVeto2Coll.at(0).Eta(), weight, -5., 5., 100);
       }
       FillHist("MET", met, weight, 0., 200., 200);
       FillHist("MET_noxy", PFMETType1, weight, 0., 200., 200);
       FillHist("MTW", MTW, weight, 0., 200., 200);
       FillHist("MTW_noxy", MTW, weight, 0., 200., 200);
       FillHist("MET_NvtxRW", met, weight*nvtx_reweight, 0., 200., 200);
       FillHist("MTW_NvtxRW", MTW, weight*nvtx_reweight, 0., 200., 200);
       if(met>30){
         FillHist("PTe_met30", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_met30", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Nj_met30", jetVeto2Coll.size(), weight, 0., 10., 10);
         if(jetVeto2Coll.size()>0){
           FillHist("PTj1_met30", jetVeto2Coll.at(0).Pt(), weight, 0., 200., 200);
           FillHist("Etaj1_met30", jetVeto2Coll.at(0).Eta(), weight, -5., 5., 100);
         }
         FillHist("MET_met30", met, weight, 0., 200., 200);
         if(muonLooseColl.size()==0) FillHist("MET_met30MuLVeto", met, weight, 0., 200., 200);
         FillHist("MTW_met30", MTW, weight, 0., 200., 200);



         //NvtxRWeighted
         FillHist("PTe_met30_NvtxRW", electronColl.at(0).Pt(), weight*nvtx_reweight, 0., 200., 200);
         FillHist("Etae_met30_NvtxRW", electronColl.at(0).Eta(), weight*nvtx_reweight, -5., 5., 100);
         FillHist("Nj_met30_NvtxRW", jetVeto2Coll.size(), weight*nvtx_reweight, 0., 10., 10);
         if(jetVeto2Coll.size()>0){
           FillHist("PTj1_met30_NvtxRW", jetVeto2Coll.at(0).Pt(), weight*nvtx_reweight, 0., 200., 200);
           FillHist("Etaj1_met30_NvtxRW", jetVeto2Coll.at(0).Eta(), weight*nvtx_reweight, -5., 5., 100);
         }
         FillHist("MET_met30_NvtxRW", met, weight*nvtx_reweight, 0., 200., 200);
         FillHist("MTW_met30_NvtxRW", MTW, weight*nvtx_reweight, 0., 200., 200);
       }
       //For both prescaled & unprescaled
       if(jetVeto2Coll.size()>0 && jetVeto2Coll.at(0).Pt()>40){
         FillHist("PTe_gt1j40", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_gt1j40", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Nj_gt1j40", jetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("PTj1_gt1j40", jetVeto2Coll.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etaj1_gt1j40", jetVeto2Coll.at(0).Eta(), weight, -5., 5., 100);
         FillHist("MET_gt1j40", met, weight, 0., 200., 200);
         FillHist("MTW_gt1j40", MTW, weight, 0., 200., 200);
         if(met>30){
           FillHist("PTe_gt1j40met30", electronColl.at(0).Pt(), weight, 0., 200., 200);
           FillHist("Etae_gt1j40met30", electronColl.at(0).Eta(), weight, -5., 5., 100);
           FillHist("Nj_gt1j40met30", jetVeto2Coll.size(), weight, 0., 10., 10);
           FillHist("PTj1_gt1j40met30", jetVeto2Coll.at(0).Pt(), weight, 0., 200., 200);
           FillHist("Etaj1_gt1j40met30", jetVeto2Coll.at(0).Eta(), weight, -5., 5., 100);

           FillHist("MET_gt1j40met30", met, weight, 0., 200., 200);
           FillHist("MTW_gt1j40met30", MTW, weight, 0., 200., 200);
         }

         //NvtxRW
         FillHist("PTe_gt1j40_NvtxRW", electronColl.at(0).Pt(), weight*nvtx_reweight, 0., 200., 200);
         FillHist("Etae_gt1j40_NvtxRW", electronColl.at(0).Eta(), weight*nvtx_reweight, -5., 5., 100);
         FillHist("Nj_gt1j40_NvtxRW", jetVeto2Coll.size(), weight*nvtx_reweight, 0., 10., 10);
         FillHist("PTj1_gt1j40_NvtxRW", jetVeto2Coll.at(0).Pt(), weight*nvtx_reweight, 0., 200., 200);
         FillHist("Etaj1_gt1j40_NvtxRW", jetVeto2Coll.at(0).Eta(), weight*nvtx_reweight, -5., 5., 100);
         FillHist("MET_gt1j40_NvtxRW", met, weight*nvtx_reweight, 0., 200., 200);
         FillHist("MTW_gt1j40_NvtxRW", MTW, weight*nvtx_reweight, 0., 200., 200);
         if(met>30){
           FillHist("MTW_gt1j40met30_NvtxRW", MTW, weight*nvtx_reweight, 0., 200., 200);
         }
       }

       //NormCheck
       if(met>30 && MTW>70){
         FillHist("Count_NormCR", 0., weight, 0., 5., 5);
         if(jetVeto2Coll.size()>0){
            if(jetVeto2Coll.at(0).Pt()>40){
              FillHist("Count_NormCR", 1., weight, 0., 5., 5);
              FillHist("Count_NormCR", 2., weight*nvtx_reweight, 0., 5., 5);
            }
         }
       }
     }//End of 1 electron 

     if(electronColl.size()==2 && electronFakeL2Coll.size()==2 && muonLooseColl.size()==0){
       bool Execute=false;
       if( SiglPreTrig && electronColl.at(0).Charge()!=electronColl.at(1).Charge()
                       && electronColl.at(0).Pt()>25. && electronColl.at(1).Pt()>10. 
                       && fabs((electronColl.at(0)+electronColl.at(1)).M()-91.2)<15.   ){ Execute=true; }
       if( UnPreTrig && electronColl.at(0).Charge()!=electronColl.at(1).Charge()
                     && electronColl.at(0).Pt()>30. && electronColl.at(1).Pt()>10. 
                     && fabs((electronColl.at(0)+electronColl.at(1)).M()-91.2)<15. ){ Execute=true;}

       if(Execute){
         FillHist("Mee_e25e10", (electronColl.at(0)+electronColl.at(1)).M(), weight, 60., 120., 60);
         FillHist("Count_NormCR", 3., weight, 0., 5., 5);
         if(jetVeto2Coll.size()>0 && jetVeto2Coll.at(0).Pt()>40){
           FillHist("Mee_e25e10gt1j", (electronColl.at(0)+electronColl.at(1)).M(), weight, 60., 120., 60);
           FillHist("Count_NormCR", 4., weight, 0., 5., 5);
           FillHist("MET_Zwin_e25e10gt1j", met, weight, 0., 200., 200);
           FillHist("MET_Zwin_e25e10gt1j_NvtxRW", met, weight*nvtx_reweight, 0., 200., 200);
           
           FillHist("Nvtx_2OSEle_PURW", Nvtx, weight, 0., 50., 50);
           if(isData){
             FillHist("Nvtx_nocut_data", Nvtx, weight, 0., 50., 50);
           }
           else{
             FillHist("Nvtx_nocut_mc", Nvtx, weight, 0., 50., 50);
           }

         }
       }
     }//End of DY Count 

     }//End of not SystRun
     
   }
   if(EleFR){

     std::vector<snu::KJet> jetVeto2Coll = SkimJetColl(jetColl, electronFakeL2Coll, muonLooseColl, "EleMuVeto");
     //float nvtx_reweight=mcdata_correction->UserPileupWeight(eventbase->GetEvent(), jetVeto2Coll.size());
     //float purw_down=0.; if(pileup_reweight!=0) purw_down=eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down)/pileup_reweight; 
     //pileup_reweight electronFakeVetoColl


     if(electronFakeL2Coll.size()==1 && jetVeto2Coll.size()>0){
       bool PassSel=true; float MTW=0.;
       if( MultPreTrig && electronFakeL2Coll.at(0).Pt()<15 ) PassSel=false;
       if( SiglPreTrig && electronFakeL2Coll.at(0).Pt()<25 ) PassSel=false;
       if( PassSel && !(jetVeto2Coll.at(0).Pt()>40) ) PassSel=false;
       if( PassSel && !(electronFakeL2Coll.at(0).DeltaR(jetVeto2Coll.at(0))>1.0) ) PassSel=false;
       if( PassSel ){
         MTW = sqrt(2)*sqrt(met*electronFakeL2Coll.at(0).Pt()-met_x*electronFakeL2Coll.at(0).Px()-met_y*electronFakeL2Coll.at(0).Py());
         
         for(int i=1; i<=20; i++){
           if(met<i*5.){
             for(int j=1; j<=20; j++){
               if(MTW<j*5.){
                 FillHist("NEvt_ltMETltMTW_2D", (i-1)*5., (j-1)*5., weight, 0., 100., 20, 0., 100., 20); 
                 if(PassIDCriteria(electronFakeL2Coll.at(0),"POGMVAMIP")){
                   FillHist("NIDEvt_ltMETltMTW_2D", (i-1)*5., (j-1)*5., weight, 0., 100., 20, 0., 100., 20); 
                 }
               }
             }
           }
         }

         if( !(met<25 && MTW<35) ) PassSel=false;
       }
       if( PassSel ){
         float PTCorr=electronFakeL2Coll.at(0).Pt()*(1.+electronFakeL2Coll.at(0).PFRelIso(0.3));
         float fEta=fabs(electronFakeL2Coll.at(0).Eta());
         const int NPtEdges=7, NEtaEdges=4;
         float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 200.};
         float EtaEdges[NEtaEdges]={0., 0.8, 1.479, 2.5};
  
       //  if(PTCorr<25.) return;

         if( fEta<1.479 ) FillHist("EleEBOpt2SumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
         else             FillHist("EleEEOpt2SumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
         FillHist("EleOpt2SumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
         FillHist("EleOpt2SumW_PTEta_2D", PTCorr, fEta, weight, PtEdges, NPtEdges-1, EtaEdges, NEtaEdges-1);
         if(PassIDCriteria(electronFakeL2Coll.at(0),"POGMVAMIP")){
           if( fEta<1.479 ) FillHist("EleEBOpt2IDSumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
           else             FillHist("EleEEOpt2IDSumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
           FillHist("EleOpt2IDSumW_PT_1D", PTCorr, weight, PtEdges, NPtEdges-1);
           FillHist("EleOpt2IDSumW_PTEta_2D", PTCorr, fEta, weight, PtEdges, NPtEdges-1, EtaEdges, NEtaEdges-1);
         }

         //Validation
         FillHist("dRej1", electronFakeL2Coll.at(0).DeltaR(jetVeto2Coll.at(0)),weight, 0., 5., 50);
         FillHist("PTe1", electronFakeL2Coll.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae1", electronFakeL2Coll.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTj1", jetVeto2Coll.at(0).Pt(), weight, 0., 500., 500);
         FillHist("MET", met, weight, 0., 200., 200);
         FillHist("MTW", MTW, weight, 0., 200., 200);
         FillHist("Nj", jetVeto2Coll.size(), weight, 0., 10., 10);
         
       }
     }//Endof Opt2

   }
   if(FRScan){
     
     //Full Scan
//     const int NIsoCuts=7, NMVACuts=21;
//     float IsoCuts[NIsoCuts]={0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
//     float MVACuts[NMVACuts]={-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};

     //Short Scan
     //const int NIsoCuts=2, NMVACuts=4;
     const int NIsoCuts=2, NMVACuts=4;
     float IsoCuts[NIsoCuts]={0., 0.4};
     //float MVACuts[NMVACuts]={-0.1, 0., 0.1, 1.};
     float MVACuts[NMVACuts]={-0.92, -0.85, -0.78, 1.};

     const int NPtEdges=7;
     float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 200.};

     //ScanFakeRate(electronFakeLPreOptColl, muonLooseColl, jetColl, bjetColl, met, met_x, met_y, NMVACuts, MVACuts, NIsoCuts, IsoCuts, NPtEdges, PtEdges, "SiglPreTrig2D1D");
     ScanFakeRate(electronFakeLPreOptColl, muonLooseColl, jetColl, bjetColl, met, met_x, met_y, NMVACuts, MVACuts, NIsoCuts, IsoCuts, NPtEdges, PtEdges, "SiglPreTrig1D");
   }
   if(Closure){
      

     bool SystRun=false;
     if(!SystRun){
       if( k_running_nonprompt )  electronColl=electronFakeL2Coll;
       if( !(electronFakeL2Coll.size()==1 && muonLooseColl.size()==2) ) return;
       if( !(electronColl.size()==1 && muonColl.size()==2) ) return;
       if( !(muonColl.at(0).Charge()!=muonColl.at(1).Charge()) ) return;
       if( !(electronColl.at(0).Pt()>25 && muonColl.at(0).Pt()>15 && muonColl.at(1).Pt()>10) ) return;
       if( fabs(electronColl.at(0).Eta())>2.5 ) return;
       if( k_running_nonprompt ){
         float fakeweight=-1.; int NLooseNotTight=0;
         for(int i=0; i<muonLooseColl.size(); i++){
           if(muonLooseColl.at(i).RelIso04()>0.1){
           //if(muonLooseColl.at(i).RelIso04()*muonLooseColl.at(i).Pt()/muonLooseColl.at(i).RochPt()>0.1){
             float FR=m_datadriven_bkg->GetFakeObj()->getTrilepFakeRate_muon(false, muonLooseColl.at(i).Pt(), muonLooseColl.at(i).Eta());
             fakeweight*=-FR/(1-FR);
             NLooseNotTight++;
           }
         }
         for(int i=0; i<electronFakeL2Coll.size(); i++){
           if(!PassIDCriteria(electronFakeL2Coll.at(i), "POGMVAMIP")){
             //float FR=FakeRateData(electronFakeL2Coll.at(i), "POGMVAMIPOpt2Iso04");
             //float FR=FakeRateData(electronFakeL2Coll.at(i), "POGMVAMTIsop1IPp025p05Sig4FakeLIsop4");
             float FR=FakeRateData(electronFakeL2Coll.at(i), "POGMVAMTIsop06IPp025p05Sig4FakeLIsop4");
             //float FR=FakeRateData(electronFakeL2Coll.at(i), "POGMVATTIsop1IPp025p05Sig4FakeLIsop4");
             //float FR=FakeRateData(electronFakeL2Coll.at(i), "POGMVATTIsop06IPp025p05Sig4FakeLIsop4");
            
             fakeweight*=-FR/(1-FR);
             NLooseNotTight++;
           }
         }
         if(NLooseNotTight==0) return;
         weight*=fakeweight;
       }
  
       float Mmumu=(muonColl.at(0)+muonColl.at(1)).M();
       float MTW = sqrt(2)*sqrt(met*electronColl.at(0).Pt()-met_x*electronColl.at(0).Px()-met_y*electronColl.at(0).Py());
       float M3l=(electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M();
       if(Mmumu<12) return;
  
  
       std::vector<snu::KJet> jetVeto2Coll  = SkimJetColl(jetColl, electronFakeL2Coll, muonLooseColl, "EleMuVeto");
       std::vector<snu::KJet> bjetVeto2Coll = SelBJets(jetVeto2Coll, "Medium");
       if(!isData)  weight *= BTagScaleFactor_1a(jetVeto2Coll, snu::KJet::CSVv2, snu::KJet::Medium);
  
       bool HasBJet=false, OnZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false, ZbbSel=false, AN_Sideband=false;
       if( bjetVeto2Coll.size()!=0 )                      HasBJet = true;
       if( fabs(Mmumu-91.2)<10     )                      OnZ     =true;
       if( fabs(M3l-91.2)<10       )                      OnZG    =true;
       if( OnZ && M3l>101.2 && met>50 )                   WZSel   =true;
       if( OnZG && Mmumu<81.2 && met<50 )                 ZGSel   =true;
       if( OnZ && HasBJet && jetVeto2Coll.size()>2)       ttZSel  =true;
       if( OnZ && HasBJet && jetVeto2Coll.size()<3)       ZbbSel  =true;
       if( Mmumu>40 && HasBJet && jetVeto2Coll.size()>2 ) AN_Sideband=true;
        
       if(isData){
         bool HasWeirdEle=HasStrangeEleCand(electronColl);
         if(HasWeirdEle){
           if(!k_running_nonprompt) FillHist("NEvtHasWeirdEleT", 1., 1, 0., 2., 2);
           else                     FillHist("NEvtHasWeirdEleL", 1., 1, 0., 2., 2);
         }
         else{
           if(!k_running_nonprompt) FillHist("NEvtHasWeirdEleT", 0., 1, 0., 2., 2);
           else                     FillHist("NEvtHasWeirdEleL", 0., 1, 0., 2., 2);
         }

         if(HasBJet && jetVeto2Coll.size()>2){
           if(HasWeirdEle){
             if(!k_running_nonprompt) FillHist("NEvtHasWeirdEleT_SigReg", 1., 1, 0., 2., 2);
             else                     FillHist("NEvtHasWeirdEleL_SigReg", 1., 1, 0., 2., 2);
           }
           else{
             if(!k_running_nonprompt) FillHist("NEvtHasWeirdEleT_SigReg", 0., 1, 0., 2., 2);
             else                     FillHist("NEvtHasWeirdEleL_SigReg", 0., 1, 0., 2., 2);
           }
         }
       }  

  
       //General 3l Selection
       if(!HasBJet){
         FillHist("PTe_3lOS", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_3lOS", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_3lOS", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOS", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOS", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOS", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_3lOS", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_3lOS", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_3lOS", jetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("Nb_3lOS", bjetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("MET_3lOS", met, weight, 0., 200., 200);
         FillHist("MTW_3lOS", MTW, weight, 0., 200., 200);
       }
       //3l OnZ ; Fake CR
       if(!HasBJet && OnZ){
         FillHist("PTe_3lOSOnZ", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_3lOSOnZ", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_3lOSOnZ", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOSOnZ", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOSOnZ", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOSOnZ", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_3lOSOnZ", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_3lOSOnZ", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_3lOSOnZ", jetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("Nb_3lOSOnZ", bjetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("MET_3lOSOnZ", met, weight, 0., 200., 200);
         FillHist("MTW_3lOSOnZ", MTW, weight, 0., 200., 200);
         if(!OnZG){
           FillHist("PTe_3lOSOnZOffZG", electronColl.at(0).Pt(), weight, 0., 200., 200);
           FillHist("Etae_3lOSOnZOffZG", electronColl.at(0).Eta(), weight, -5., 5., 100);
           FillHist("PTmu1_3lOSOnZOffZG", muonColl.at(0).Pt(), weight, 0., 200., 200);
           FillHist("PTmu2_3lOSOnZOffZG", muonColl.at(1).Pt(), weight, 0., 200., 200);
           FillHist("Etamu1_3lOSOnZOffZG", muonColl.at(0).Eta(), weight, -5., 5., 100);
           FillHist("Etamu2_3lOSOnZOffZG", muonColl.at(1).Eta(), weight, -5., 5., 100);
    
           FillHist("Mmumu_3lOSOnZOffZG", Mmumu, weight, 0., 200., 200);
           FillHist("M3l_3lOSOnZOffZG", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
           FillHist("Nj_3lOSOnZOffZG", jetVeto2Coll.size(), weight, 0., 10., 10);
           FillHist("Nb_3lOSOnZOffZG", bjetVeto2Coll.size(), weight, 0., 10., 10);
           FillHist("MET_3lOSOnZOffZG", met, weight, 0., 200., 200);
           FillHist("MTW_3lOSOnZOffZG", MTW, weight, 0., 200., 200);
         }
       }
       //3l OnZ & HasBJet
       if(HasBJet && OnZ){
         FillHist("PTe_3lOSOnZHasB", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_3lOSOnZHasB", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_3lOSOnZHasB", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_3lOSOnZHasB", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_3lOSOnZHasB", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_3lOSOnZHasB", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_3lOSOnZHasB", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_3lOSOnZHasB", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_3lOSOnZHasB", jetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("Nb_3lOSOnZHasB", bjetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("MET_3lOSOnZHasB", met, weight, 0., 200., 200);
         FillHist("MTW_3lOSOnZHasB", MTW, weight, 0., 200., 200);
       }
       if(WZSel){
         FillHist("PTe_WZSel", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_WZSel", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_WZSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_WZSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_WZSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_WZSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_WZSel", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_WZSel", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_WZSel", jetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("Nb_WZSel", bjetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("MET_WZSel", met, weight, 0., 200., 200);
         FillHist("MTW_WZSel", MTW, weight, 0., 200., 200);
       }
       if(ZGSel){
         FillHist("PTe_ZGSel", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_ZGSel", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_ZGSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_ZGSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_ZGSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_ZGSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_ZGSel", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_ZGSel", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_ZGSel", jetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("Nb_ZGSel", bjetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("MET_ZGSel", met, weight, 0., 200., 200);
         FillHist("MTW_ZGSel", MTW, weight, 0., 200., 200);
       }
       if(ttZSel){
         FillHist("PTe_ttZSel", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_ttZSel", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_ttZSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_ttZSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_ttZSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_ttZSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_ttZSel", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_ttZSel", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_ttZSel", jetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("Nb_ttZSel", bjetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("MET_ttZSel", met, weight, 0., 200., 200);
         FillHist("MTW_ttZSel", MTW, weight, 0., 200., 200);
       }
       if(ZbbSel){
         FillHist("PTe_ZbbSel", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_ZbbSel", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_ZbbSel", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_ZbbSel", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_ZbbSel", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_ZbbSel", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_ZbbSel", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_ZbbSel", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_ZbbSel", jetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("Nb_ZbbSel", bjetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("MET_ZbbSel", met, weight, 0., 200., 200);
         FillHist("MTW_ZbbSel", MTW, weight, 0., 200., 200);
       }
       if(AN_Sideband){
         FillHist("PTe_ANSideband", electronColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("Etae_ANSideband", electronColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("PTmu1_ANSideband", muonColl.at(0).Pt(), weight, 0., 200., 200);
         FillHist("PTmu2_ANSideband", muonColl.at(1).Pt(), weight, 0., 200., 200);
         FillHist("Etamu1_ANSideband", muonColl.at(0).Eta(), weight, -5., 5., 100);
         FillHist("Etamu2_ANSideband", muonColl.at(1).Eta(), weight, -5., 5., 100);
  
         FillHist("Mmumu_ANSideband", Mmumu, weight, 0., 200., 200);
         FillHist("M3l_ANSideband", (electronColl.at(0)+muonColl.at(0)+muonColl.at(1)).M(), weight, 0., 500., 500);
         FillHist("Nj_ANSideband", jetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("Nb_ANSideband", bjetVeto2Coll.size(), weight, 0., 10., 10);
         FillHist("MET_ANSideband", met, weight, 0., 200., 200);
         FillHist("MTW_ANSideband", MTW, weight, 0., 200., 200);
       }
     }//End of Not SystRun
     if(SystRun){


         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESUpColl; eventbase->GetJetSel()->Selection(jetJESUpColl, LeptonVeto, muonLooseColl, electronFakeLPreOptColl, "SystUpJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERUpColl; eventbase->GetJetSel()->Selection(jetJERUpColl, LeptonVeto, muonLooseColl, electronFakeLPreOptColl, "SystUpJER");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESDownColl; eventbase->GetJetSel()->Selection(jetJESDownColl, LeptonVeto, muonLooseColl, electronFakeLPreOptColl, "SystDownJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERDownColl; eventbase->GetJetSel()->Selection(jetJERDownColl, LeptonVeto, muonLooseColl, electronFakeLPreOptColl, "SystDownJER");

         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
         eventbase->GetMuonSel()->SetBSdz(0.1);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
       std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,true, "SystUpMuEn");
         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(3.);
         eventbase->GetMuonSel()->SetBSdz(0.1);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.1);
       std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,true, "SystDownMuEn");

         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
         eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
         eventbase->GetElectronSel()->SetdxySigMax(3.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");

         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
         eventbase->GetElectronSel()->SetdxyBEMax(0.05, 0.05);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
         eventbase->GetElectronSel()->SetdxySigMax(3.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnDownColl; eventbase->GetElectronSel()->Selection(EleTElEnDownColl, "SystDownElEn");



       DoSystRun("Closure", "", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight);
       if(!isData){
         DoSystRun("Closure", "SystUpPU", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight*pileup_reweight_systup);
         DoSystRun("Closure", "SystUpUncl", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi), weight);
         DoSystRun("Closure", "SystUpJES", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetJESUpColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi), weight);
         DoSystRun("Closure", "SystUpJER", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetJERUpColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi), weight);
         DoSystRun("Closure", "SystUpBTag_L", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight);
         DoSystRun("Closure", "SystUpBTag_BC", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight);
         DoSystRun("Closure", "SystUpElEn", EleTElEnUpColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi), weight);
         DoSystRun("Closure", "SystUpMuEn", electronColl, electronFakeL2Coll, MuTMuEnUpColl, muonLooseColl, jetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi), weight);


         DoSystRun("Closure", "SystDownPU", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight*pileup_reweight_systdown);
         DoSystRun("Closure", "SystDownUncl", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi), weight);
         DoSystRun("Closure", "SystDownJES", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetJESDownColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi), weight);
         DoSystRun("Closure", "SystDownJER", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetJERDownColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi), weight);
         DoSystRun("Closure", "SystDownBTag_L", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight);
         DoSystRun("Closure", "SystDownBTag_BC", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight);
         DoSystRun("Closure", "SystDownElEn", EleTElEnDownColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi), weight);
         DoSystRun("Closure", "SystDownMuEn", electronColl, electronFakeL2Coll, MuTMuEnDownColl, muonLooseColl, jetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi), weight);
       }
       else if( isData && k_running_nonprompt ){

         DoSystRun("Closure", "", electronFakeL2Coll, electronFakeL2Coll, muonLooseColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight);
//
//         DoSystRun("Closure", "SystUpEleFR", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight);
//         DoSystRun("Closure", "SystDownEleFR", electronColl, electronFakeL2Coll, muonColl, muonLooseColl, jetColl, met, met*cos(metphi), met*sin(metphi), weight);
       }

     }//End of SystRun
   }//End of Closure
   if(HighdXYCheck){
     std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);

       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> EleHighdxyTColl; eventbase->GetElectronSel()->Selection(EleHighdxyTColl);

       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_HctoWA_FAKELOOSE);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.4, 0.4);
       eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> EleHighdxyLColl; eventbase->GetElectronSel()->Selection(EleHighdxyLColl);

     std::vector<snu::KJet> jetVetoColl    = SkimJetColl(jetColl, EleHighdxyLColl   , muonLooseColl, "EleMuVeto");
     std::vector<snu::KJet> jetLowVetoColl = SkimJetColl(jetColl, electronFakeL2Coll, muonLooseColl, "EleMuVeto");


     //Series of Cuts for FR measurement--------------------------------------------------------------------------------//
     bool PassSelHighd0 = true, PassSelLowd0=true; float MTW=0.;
     if( EleHighdxyLColl.size()!=1 ) PassSelHighd0=false;
     if( PassSelHighd0 && !(EleHighdxyLColl.at(0).Pt()>25. && fabs(EleHighdxyLColl.at(0).dxySig())>4.) ) PassSelHighd0=false;
     if( PassSelHighd0 && !(jetVetoColl.size()>0 && jetVetoColl.at(0).Pt()>40) ) PassSelHighd0=false;
     if( PassSelHighd0 ){
       bool HasAway40Jet=false;
       for( int i=0; i<jetVetoColl.size(); i++){
         if(jetVetoColl.at(i).Pt()<40) continue;
         if(jetVetoColl.at(i).DeltaR(EleHighdxyLColl.at(0))>1.0) {HasAway40Jet=true; break;}
       }
       if(!HasAway40Jet) PassSelHighd0=false;
     }
     if( !isData && PassSelHighd0 && !(k_sample_name.Contains("qcd") || k_sample_name.Contains("QCD")) ){
       int EleType=isData? 0 : GetLeptonType(EleHighdxyLColl.at(0),truthColl);
       if( EleType<0 && EleType>-5 ) PassSelHighd0=false;// for MC, we only subtract that originate from prompts. Fakes are estimated in data.
     }


     if( electronFakeL2Coll.size()!=1 ) PassSelLowd0=false;
     if( PassSelLowd0 && !(electronFakeL2Coll.at(0).Pt()>25.) ) PassSelLowd0=false;
     if( PassSelLowd0 && !(jetLowVetoColl.size()>0 && jetLowVetoColl.at(0).Pt()>40) ) PassSelLowd0=false;
     if( PassSelLowd0 ){
       bool HasAway40Jet=false;
       for( int i=0; i<jetLowVetoColl.size(); i++){
         if(jetLowVetoColl.at(i).Pt()<40) continue;
         if(jetLowVetoColl.at(i).DeltaR(electronFakeL2Coll.at(0))>1.0) {HasAway40Jet=true; break;}
       }
       if(!HasAway40Jet) PassSelLowd0=false;
     }
     if( !isData && PassSelLowd0 && !(k_sample_name.Contains("qcd") || k_sample_name.Contains("QCD")) ){
       int EleType=isData? 0 : GetLeptonType(electronFakeL2Coll.at(0),truthColl);
       if( EleType<0 && EleType>-5 ) PassSelLowd0=false;// for MC, we only subtract that originate from prompts. Fakes are estimated in data.
     }


     //-----------------------------------------------------------------------------------------------------------------//
     
     if( PassSelHighd0 ){
       
       float PTCorr=EleHighdxyLColl.at(0).Pt()*(1.+EleHighdxyLColl.at(0).PFRelIso(0.3));
       float fEta=fabs(EleHighdxyLColl.at(0).Eta());
       const int NPtEdges=8; float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 150., 200.};
       bool  IsNearB = IsNearBJet(EleHighdxyLColl.at(0), bjetColl);
       float MTW = sqrt(2)*sqrt(met*EleHighdxyLColl.at(0).Pt()-met_x*EleHighdxyLColl.at(0).Px()-met_y*EleHighdxyLColl.at(0).Py());
       bool PassMETMTW=(met<25 && MTW<35);


       //Pass Loose
       if(fEta<0.8){
         FillHist("EleHighd0B1SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         if(IsNearB) FillHist("EleHighd0B1SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
       }
       else if(fEta<1.479){
         FillHist("EleHighd0B2SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         if(IsNearB) FillHist("EleHighd0B2SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
       }
       else{
         FillHist("EleHighd0ESumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         if(IsNearB) FillHist("EleHighd0ESumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
       }

       if(PassMETMTW){
         if(fEta<0.8){
           FillHist("EleHighd0METMTWB1SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleHighd0METMTWB1SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else if(fEta<1.479){
           FillHist("EleHighd0METMTWB2SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleHighd0METMTWB2SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else{
           FillHist("EleHighd0METMTWESumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleHighd0METMTWESumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
       }

       //Pass IDCriteria
       if(EleHighdxyTColl.size()==1){

         if(fEta<0.8){
           FillHist("EleHighd0B1IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleHighd0B1IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else if(fEta<1.479){
           FillHist("EleHighd0B2IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleHighd0B2IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else{
           FillHist("EleHighd0EIDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleHighd0EIDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
  
         if(PassMETMTW){
           if(fEta<0.8){
             FillHist("EleHighd0METMTWB1IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
             if(IsNearB) FillHist("EleHighd0METMTWB1IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           }
           else if(fEta<1.479){
             FillHist("EleHighd0METMTWB2IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
             if(IsNearB) FillHist("EleHighd0METMTWB2IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           }
           else{
             FillHist("EleHighd0METMTWEIDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
             if(IsNearB) FillHist("EleHighd0METMTWEIDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           }
         }
       }//End of highd0 Tight
     }//End of High d0 MeasReg

     //QCDMC Highd0 to Lowd0 SF Measurement-----------------------------------------------------------------------------//
     if( !isData && (k_sample_name.Contains("qcd") || k_sample_name.Contains("QCD") || k_sample_name.Contains("TT_powheg")) && PassSelLowd0 ){

       float PTCorr=electronFakeL2Coll.at(0).Pt()*(1.+electronFakeL2Coll.at(0).PFRelIso(0.3));
       float fEta=fabs(electronFakeL2Coll.at(0).Eta());
       const int NPtEdges=8; float PtEdges[NPtEdges]={0., 25., 35., 50., 70., 100., 150., 200.};
       bool  IsNearB = IsNearBJet(electronFakeL2Coll.at(0), bjetColl);
       float MTW = sqrt(2)*sqrt(met*electronFakeL2Coll.at(0).Pt()-met_x*electronFakeL2Coll.at(0).Px()-met_y*electronFakeL2Coll.at(0).Py());
       bool PassMETMTW=(met<25 && MTW<35);


       //Pass Loose
       if(fEta<0.8){
         FillHist("EleLowd0B1SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         if(IsNearB) FillHist("EleLowd0B1SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
       }
       else if(fEta<1.479){
         FillHist("EleLowd0B2SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         if(IsNearB) FillHist("EleLowd0B2SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
       }
       else{
         FillHist("EleLowd0ESumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         if(IsNearB) FillHist("EleLowd0ESumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
       }

       if(PassMETMTW){
         if(fEta<0.8){
           FillHist("EleLowd0METMTWB1SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleLowd0METMTWB1SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else if(fEta<1.479){
           FillHist("EleLowd0METMTWB2SumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleLowd0METMTWB2SumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else{
           FillHist("EleLowd0METMTWESumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleLowd0METMTWESumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
       }

       //Pass IDCriteria
       if(electronColl.size()==1){

         if(fEta<0.8){
           FillHist("EleLowd0B1IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleLowd0B1IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else if(fEta<1.479){
           FillHist("EleLowd0B2IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleLowd0B2IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
         else{
           FillHist("EleLowd0EIDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           if(IsNearB) FillHist("EleLowd0EIDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
         }
  
         if(PassMETMTW){
           if(fEta<0.8){
             FillHist("EleLowd0METMTWB1IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
             if(IsNearB) FillHist("EleLowd0METMTWB1IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           }
           else if(fEta<1.479){
             FillHist("EleLowd0METMTWB2IDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
             if(IsNearB) FillHist("EleLowd0METMTWB2IDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           }
           else{
             FillHist("EleLowd0METMTWEIDSumW_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
             if(IsNearB) FillHist("EleLowd0METMTWEIDSumW_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
           }
         }
       }//End of highd0 Tight
     }//End of QCDMC Lowd0 MeasReg
   }//End of HighdXYCheck
   if(IsoIPOpt){
     std::vector<snu::KElectron> electronPOGMVAMColl;
       for(int i=0; i<electronPreColl.size(); i++){
         if(electronPreColl.at(i).Pt()<25) continue;
         if(PassIDCriteria(electronPreColl.at(i), "POGMVAM")) electronPOGMVAMColl.push_back(electronPreColl.at(i));
       }
     if(electronPOGMVAMColl.size()!=2) return;
     if(fabs((electronPOGMVAMColl.at(0)+electronPOGMVAMColl.at(1)).M()-91.2)>15) return;

     for(int i=0; i<electronPOGMVAMColl.size(); i++){
       float Iso=electronPOGMVAMColl.at(i).PFRelIso(0.3);
       float d0=fabs(electronPOGMVAMColl.at(i).dxy()), dz=fabs(electronPOGMVAMColl.at(i).dz());

       FillHist("Data_Iso", electronPOGMVAMColl.at(i).PFRelIso(0.3), weight, 0., 0.3, 300);
       FillHist("Data_d0", d0, weight, 0., 0.1, 20);
       FillHist("Data_dz", dz, weight, 0., 0.2, 40);

       //IsoWP Study
       int NIsoWP=40;
       for(int it_iso=1; it_iso<=NIsoWP; it_iso++){
         FillHist("NEleSumW_Iso", (it_iso-1.)*0.00501, weight, 0., 0.2, 40);
         if(electronPOGMVAMColl.at(i).PFRelIso(0.3)<it_iso*0.005){
           FillHist("NEleIDSumW_Iso", (it_iso-1.)*0.00501, weight, 0., 0.2, 40);
         }
       }


       int Nd0WP=10, NdzWP=20;
       for(int it_d0=1; it_d0<=Nd0WP; it_d0++){
         for(int it_dz=1; it_dz<=NdzWP; it_dz++){
           if(electronPOGMVAMColl.at(i).PFRelIso(0.3)<0.06){
             FillHist("NEleSumW_d0dz_iso0p06", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             if(d0<it_d0*0.005 && dz< it_dz*0.005){
               FillHist("NEleIDSumW_d0dz_iso0p06", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             }
           }
           if(electronPOGMVAMColl.at(i).PFRelIso(0.3)<0.08){
             FillHist("NEleSumW_d0dz_iso0p08", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             if(d0<it_d0*0.005 && dz< it_dz*0.005){
               FillHist("NEleIDSumW_d0dz_iso0p08", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             }
           }
           if(electronPOGMVAMColl.at(i).PFRelIso(0.3)<0.1){
             FillHist("NEleSumW_d0dz_iso0p1", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             if(d0<it_d0*0.005 && dz< it_dz*0.005){
               FillHist("NEleIDSumW_d0dz_iso0p1", (it_d0-1)*0.00501, (it_dz-1)*0.00501, weight, 0., 0.05, Nd0WP, 0., 0.1, NdzWP);
             }
           }
         }
       }
     }//End of EleLoop
   }//End of IsoIPOpt
   if(IDValidation){
     
     bool PassSel=true; float Mee=0.;
     if( electronColl.size()!=2 ) PassSel=false;
     if( PassSel && !(electronColl.at(0).Pt()>25 && electronColl.at(1).Pt()>15) ) PassSel=false;
     if( PassSel &&   electronColl.at(0).Charge()== electronColl.at(1).Charge() ) PassSel=false;
     if( PassSel ) Mee=(electronColl.at(0)+electronColl.at(1)).M();
     if( PassSel && fabs(Mee-91.2)>15 ) PassSel=false;
     
     if(PassSel){
       FillHist("PTe1" , electronColl.at(0).Pt(),  weight, 0., 200., 200);
       FillHist("PTe2" , electronColl.at(1).Pt(),  weight, 0., 200., 200);
       FillHist("Etae1", electronColl.at(0).Eta(), weight, -5., 5., 100);
       FillHist("Etae2", electronColl.at(1).Eta(), weight, -5., 5., 100);
       FillHist("Mee"  , Mee, weight, 60., 120., 60);
     }
     
     //Check Efficiency
     if( true ){
     //if( !isData && (k_sample_name.Contains("DY") || k_sample_name.Contains("TT")) ){
     
       bool PassLoose=true;
       if( electronPreColl.size()!=2 ) PassLoose=false;
       if( PassLoose && electronPreColl.at(0).Charge()==electronPreColl.at(1).Charge() ) PassLoose=false;
       if( PassLoose && !(electronPreColl.at(0).Pt()>25 && electronPreColl.at(1).Pt()>15) ) PassLoose=false;
  
       if(PassLoose){
         for(int i=0; i<electronPreColl.size(); i++){
           FillHist("NEleSumW_PT" , electronPreColl.at(i).Pt(),  weight, 0., 200., 20);
           FillHist("NEleSumW_Eta", electronPreColl.at(i).Eta(), weight, -5., 5., 20);
           FillHist("NEleSumW_PTEta", electronPreColl.at(i).Pt(), fabs(electronPreColl.at(i).Eta()), weight, 0., 200., 20, 0, 2.5, 5);
           if(PassIDCriteria(electronPreColl.at(i), "POGMVATIP")){
             FillHist("NEleIDSumW_PT" , electronPreColl.at(i).Pt(),  weight, 0., 200., 20);
             FillHist("NEleIDSumW_Eta", electronPreColl.at(i).Eta(), weight, -5., 5., 20);
             FillHist("NEleIDSumW_PTEta", electronPreColl.at(i).Pt(), fabs(electronPreColl.at(i).Eta()), weight, 0., 200., 20, 0., 2.5, 5);
           } 
         }
       }
     }
     
   }//End of validation

/////////////////////////////////////////////////////////////////////////////////// 

return;

}// End of execute event loop
  


void Jul2017_DataFakeStudy::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Jul2017_DataFakeStudy::BeginCycle() throw( LQError ){
  
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

Jul2017_DataFakeStudy::~Jul2017_DataFakeStudy() {
  
  Message("In Jul2017_DataFakeStudy Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}


void Jul2017_DataFakeStudy::ScanFakeRate(std::vector<snu::KElectron> EleColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, int NMVACuts, float MVACuts[], int NIsoCuts, float IsoCuts[], int NPtEdges, float PtEdges[], TString Option){
//Input only fake leptons

  float MinElePt=0.;
  bool Scan1D=false, Scan2D=false;
  if(Option.Contains("SiglPreTrig")) MinElePt=25.;
  if(Option.Contains("MultPreTrig")) MinElePt=15.;
  if(Option.Contains("2D"))          Scan2D  =true;
  if(Option.Contains("1D"))          Scan1D  =true;

      
  for(int it_iso=1; it_iso<NIsoCuts; it_iso++){
    for(int it_mva=0; it_mva<NMVACuts-1; it_mva++){

      //Selection Requirement--------------------------------------------------------------//
      int NLooseEle=0;
      for(int i_ele=0; i_ele<EleColl.size(); i_ele++){
        if(EleColl.at(i_ele).PFRelIso(0.3)<IsoCuts[it_iso] && EleColl.at(i_ele).MVA()>MVACuts[it_mva]) NLooseEle++;
      }
      if(NLooseEle!=1) continue;
      if(MuLColl.size()!=0) continue;
      if(EleColl.at(0).Pt()<MinElePt) continue;
      bool PassJetReq=false; 
      for(int j=0; j<JetColl.size(); j++){
        if(JetColl.at(j).Pt()<40) continue;
        if(JetColl.at(j).DeltaR(EleColl.at(0))<1.0) continue;
        PassJetReq=true;
      }
      if(!PassJetReq) continue;
      float MTW = sqrt(2)*sqrt(MET*EleColl.at(0).Pt()-METx*EleColl.at(0).Px()-METy*EleColl.at(0).Py());
      if( !(MET<25 && MTW<35) ) continue;
      //-----------------------------------------------------------------------------------//

      float PTCorr = ConeCorrectedPT(EleColl.at(0), 0.1);
      bool  IsNearB  = IsNearBJet(EleColl.at(0), BJetColl);

      std::ostringstream s2; s2<<IsoCuts[it_iso];
      TString Str_IsoCut=s2.str();
      Str_IsoCut.ReplaceAll(".","p");

      std::ostringstream s1; s1<<MVACuts[it_mva];
      TString Str_MVACut=s1.str();
      Str_MVACut.ReplaceAll(".","p");  Str_MVACut.ReplaceAll("-","m");

      //Check the only region that is not biased by LepPt Threshold.
      if(PTCorr<(MinElePt+EleColl.at(0).PFRelIso(0.3)*IsoCuts[it_iso])) continue;

      if(EleColl.at(0).PFRelIso(0.3)<IsoCuts[it_iso] && EleColl.at(0).MVA()>MVACuts[it_mva]){

        if(fabs(EleColl.at(0).Eta())<0.8){//Barrel1
          if(Scan2D){
            FillHist("EleB1SumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
            if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleB1IDSumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
  
            if(IsNearB){
              FillHist("EleB1SumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleB1IDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }
          }
          if(Scan1D){
            FillHist("EleB1SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleB1IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
  
            if(IsNearB){
              FillHist("EleB1SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
              if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleB1IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
        }
        else if(fabs(EleColl.at(0).Eta())<1.479){//Barrel2
          if(Scan2D){
            FillHist("EleB2SumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
            if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleB2IDSumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
  
            if(IsNearB){
              FillHist("EleB2SumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleB2IDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }
          }
          if(Scan1D){
            FillHist("EleB2SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleB2IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
  
            if(IsNearB){
              FillHist("EleB2SumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
              if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleB2IDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
        }
        else if(fabs(EleColl.at(0).Eta())<2.5){//EndCap
          if(Scan2D){
            FillHist("EleESumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
            if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleEIDSumW_All_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
  
            if(IsNearB){
              FillHist("EleESumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);   
              if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleEIDSumW_BjMatch_ltIsoltMVA_FR2D", IsoCuts[it_iso-1], MVACuts[it_mva], weight, IsoCuts, NIsoCuts-1, MVACuts, NMVACuts-1);
            }
          }
          if(Scan1D){
            FillHist("EleESumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleEIDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_All_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
  
            if(IsNearB){
              FillHist("EleESumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
              if(PassIDCriteria(EleColl.at(0), "POGMVAMIP")) FillHist("EleEIDSumW_mva"+Str_MVACut+"iso"+Str_IsoCut+"_BjMatch_PT_FR1D", PTCorr, weight, PtEdges, NPtEdges-1);
            }
          }
        }
      }//End of if Iso<Cut && MVA<Cut
    }//End of MVACut loop
  }//End of IsoCut Loop


}

bool Jul2017_DataFakeStudy::IsNearBJet(snu::KElectron Ele, std::vector<snu::KJet> bjetNoVetoColl){

  bool IsNearB=false;
  for(int i=0; i<bjetNoVetoColl.size(); i++){
    if(Ele.DeltaR(bjetNoVetoColl.at(i))<0.4){
      IsNearB=true; break;
    }
  }

  return IsNearB;
}

float Jul2017_DataFakeStudy::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  //float PTCorr=Ele.Pt()*(1+min(0,Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}


/*
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
*/


void Jul2017_DataFakeStudy::DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, float MET, float METx, float METy, float weight){

  int SystDir=0; TString SystKindLabel="", SystDirLabel="";
  if(Mode.Contains("Syst")){
    if     (isData)                   return;
    if     (Mode.Contains("Up"))      {SystDir= 1; SystDirLabel="_systup";}
    else if(Mode.Contains("Down"))    {SystDir=-1; SystDirLabel="_systdown";}
    if     (Mode.Contains("PU"))      SystKindLabel="_PU";
    else if(Mode.Contains("JES"))     SystKindLabel="_JES";
    else if(Mode.Contains("JER"))     SystKindLabel="_JER";
    else if(Mode.Contains("Uncl"))    SystKindLabel="_Uncl";
    else if(Mode.Contains("ElEn"))    SystKindLabel="_ElEn";
    else if(Mode.Contains("MuEn"))    SystKindLabel="_MuEn";
    else if(Mode.Contains("BTag_L"))  SystKindLabel="_BTag_L";
    else if(Mode.Contains("BTag_BC")) SystKindLabel="_BTag_BC";
  }

  if(Cycle=="NormCheck"){

    bool UnPreTrig=false, PreTrig=false;
    if     (Mode.Contains("UnPreTrig")) UnPreTrig=true;
    else if(Mode.Contains("PreTrig"))   PreTrig  =true;

    std::vector<snu::KJet> JetVetoColl = SkimJetColl(JetColl, EleLColl, MuLColl, "EleMuVeto");

    if(EleColl.size()==1 && EleLColl.size()==1 && MuLColl.size()==0){
      if( UnPreTrig && EleColl.at(0).Pt()<30.) return;
      if( PreTrig   && EleColl.at(0).Pt()<25.) return;

      float MTW = sqrt(2)*sqrt(MET*EleColl.at(0).Pt()-METx*EleColl.at(0).Px()-METy*EleColl.at(0).Py());
      if(MET>30.){
        FillHist("MET_met30"+SystDirLabel+SystKindLabel, MET, weight, 0., 200., 200);
        FillHist("MTW_met30"+SystDirLabel+SystKindLabel, MTW, weight, 0., 200., 200);

        if(MTW>70.) FillHist("Count_NormCR"+SystDirLabel+SystKindLabel, 0., weight, 0., 10., 10);
        if(JetVetoColl.size()>0 && JetVetoColl.at(0).Pt()>40.){
          FillHist("MET_gt1j40met30"+SystDirLabel+SystKindLabel, MET, weight, 0., 200., 200);
          FillHist("MTW_gt1j40met30"+SystDirLabel+SystKindLabel, MTW, weight, 0., 200., 200);
          if(MTW>70.) FillHist("Count_NormCR"+SystDirLabel+SystKindLabel, 1., weight, 0., 10., 10); 
        }

      }

      if(JetVetoColl.size()>0 && JetVetoColl.at(0).Pt()>40.){
        FillHist("MET_gt1j40"+SystDirLabel+SystKindLabel, MET, weight, 0., 200., 200);
        FillHist("MTW_gt1j40"+SystDirLabel+SystKindLabel, MTW, weight, 0., 200., 200);
      }

    }//End of 1l

    if(EleColl.size()==2 && EleLColl.size()==2 && MuLColl.size()==0){
      if( UnPreTrig && !(EleColl.at(0).Pt()>30. && EleColl.at(1).Pt()>10.) ) return;
      if( PreTrig   && !(EleColl.at(0).Pt()>25. && EleColl.at(1).Pt()>10.) ) return;
      if( EleColl.at(0).Charge()==EleColl.at(1).Charge()  ) return;
      if( fabs((EleColl.at(0)+EleColl.at(1)).M()-91.2)>15 ) return;

      FillHist("Mee_e25e10"+SystDirLabel+SystKindLabel, (EleColl.at(0)+EleColl.at(1)).M(), weight, 60., 120., 60);
      FillHist("Count_NormCR"+SystDirLabel+SystKindLabel, 2., weight, 0., 10., 10);
      if(JetVetoColl.size()>0 && JetVetoColl.at(0).Pt()>40.){
        FillHist("Mee_e25e10gt1j"+SystDirLabel+SystKindLabel, (EleColl.at(0)+EleColl.at(1)).M(), weight, 60., 120., 60);
        FillHist("Count_NormCR"+SystDirLabel+SystKindLabel, 3., weight, 0., 10., 10);
      }
    }//End of 2l

  }//End of NormCheck
  else if(Cycle=="Closure"){

    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleColl.size() ==1 && MuColl.size() ==2) ) return;
    if( !(MuColl.at(0).Charge()!=MuColl.at(1).Charge()) ) return;
    if( !(EleColl.at(0).Pt()>25 && MuColl.at(0).Pt()>15 && MuColl.at(1).Pt()>10) ) return;
    if( fabs(EleColl.at(0).Eta())>2.5 ) return;
    if( k_running_nonprompt ){
      float fakeweight=-1.; int NLooseNotTight=0;
      for(int i=0; i<MuLColl.size(); i++){
        if(MuLColl.at(i).RelIso04()>0.1){
          float FR=m_datadriven_bkg->GetFakeObj()->getTrilepFakeRate_muon(false, MuLColl.at(i).Pt(), MuLColl.at(i).Eta());
          fakeweight*=-FR/(1-FR);
          NLooseNotTight++;
        }
      }
      for(int i=0; i<EleLColl.size(); i++){
        if(!PassIDCriteria(EleLColl.at(i), "POGMVAMIP")){
          //float FR=FakeRateData(EleLColl.at(i), "POGMVAMIPOpt2Iso04");
          float FR=FakeRateData(EleLColl.at(i), "POGMVAMTIsop1IPp025p05Sig4FakeLIsop4");
          fakeweight*=-FR/(1-FR);
          NLooseNotTight++;
        }
      }
      if(NLooseNotTight==0) return;
      weight*=fakeweight;
    }

    float Mmumu=(MuColl.at(0)+MuColl.at(1)).M();
    float MTW = sqrt(2)*sqrt(MET*EleColl.at(0).Pt()-METx*EleColl.at(0).Px()-METy*EleColl.at(0).Py());
    float M3l=(EleColl.at(0)+MuColl.at(0)+MuColl.at(1)).M();
    if(Mmumu<12) return;


    std::vector<snu::KJet> JetVetoColl  = SkimJetColl(JetColl, EleLColl, MuLColl, "EleMuVeto");
    std::vector<snu::KJet> BJetVetoColl = SelBJets(JetVetoColl, "Medium");
    if(!isData){
      if     (SystDir==0)                             weight *= BTagScaleFactor_1a(JetVetoColl, snu::KJet::CSVv2, snu::KJet::Medium);
      else if(SystKindLabel=="_BTag_L"  && SystDir>0) weight *= BTagScaleFactor_1a(JetVetoColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpLTag");
      else if(SystKindLabel=="_BTag_BC" && SystDir>0) weight *= BTagScaleFactor_1a(JetVetoColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpBCTag");
      else if(SystKindLabel=="_BTag_L"  && SystDir<0) weight *= BTagScaleFactor_1a(JetVetoColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownLTag");
      else if(SystKindLabel=="_BTag_BC" && SystDir<0) weight *= BTagScaleFactor_1a(JetVetoColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownBCTag");
    }

    bool HasBJet=false, OnZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false, ZbbSel=false;
    if( BJetVetoColl.size()!=0       )          HasBJet =true;
    if( fabs(Mmumu-91.2)<10          )          OnZ     =true;
    if( fabs(M3l-91.2)<10            )          OnZG    =true;
    if( OnZ && M3l>101.2   && MET>50 )          WZSel   =true;
    if( OnZG && Mmumu<81.2 && MET<50 )          ZGSel   =true;
    if( OnZ && HasBJet && JetVetoColl.size()>2) ttZSel  =true;
    if( OnZ && HasBJet && JetVetoColl.size()<3) ZbbSel  =true;
     


    //General 3l Selection
    if(!HasBJet){
      FillHist("PTe_3lOS"+SystDirLabel+SystKindLabel, EleColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_3lOS"+SystDirLabel+SystKindLabel, EleColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_3lOS"+SystDirLabel+SystKindLabel, MuColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_3lOS"+SystDirLabel+SystKindLabel, MuColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_3lOS"+SystDirLabel+SystKindLabel, MuColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_3lOS"+SystDirLabel+SystKindLabel, MuColl.at(1).Eta(), weight, -5., 5., 100);

      FillHist("Mmumu_3lOS"+SystDirLabel+SystKindLabel, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_3lOS"+SystDirLabel+SystKindLabel, (EleColl.at(0)+MuColl.at(0)+MuColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_3lOS"+SystDirLabel+SystKindLabel, JetVetoColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOS"+SystDirLabel+SystKindLabel, BJetVetoColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOS"+SystDirLabel+SystKindLabel, MET, weight, 0., 200., 200);
      FillHist("MTW_3lOS"+SystDirLabel+SystKindLabel, MTW, weight, 0., 200., 200);
    }
    //3l OnZ ; Fake CR
    if(!HasBJet && OnZ){
      FillHist("PTe_3lOSOnZ"+SystDirLabel+SystKindLabel, EleColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_3lOSOnZ"+SystDirLabel+SystKindLabel, EleColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_3lOSOnZ"+SystDirLabel+SystKindLabel, MuColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_3lOSOnZ"+SystDirLabel+SystKindLabel, MuColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_3lOSOnZ"+SystDirLabel+SystKindLabel, MuColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_3lOSOnZ"+SystDirLabel+SystKindLabel, MuColl.at(1).Eta(), weight, -5., 5., 100);

      FillHist("Mmumu_3lOSOnZ"+SystDirLabel+SystKindLabel, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_3lOSOnZ"+SystDirLabel+SystKindLabel, (EleColl.at(0)+MuColl.at(0)+MuColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_3lOSOnZ"+SystDirLabel+SystKindLabel, JetVetoColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOSOnZ"+SystDirLabel+SystKindLabel, BJetVetoColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOSOnZ"+SystDirLabel+SystKindLabel, MET, weight, 0., 200., 200);
      FillHist("MTW_3lOSOnZ"+SystDirLabel+SystKindLabel, MTW, weight, 0., 200., 200);
      if(!OnZG){
        FillHist("PTe_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, EleColl.at(0).Pt(), weight, 0., 200., 200);
        FillHist("Etae_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, EleColl.at(0).Eta(), weight, -5., 5., 100);
        FillHist("PTmu1_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, MuColl.at(0).Pt(), weight, 0., 200., 200);
        FillHist("PTmu2_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, MuColl.at(1).Pt(), weight, 0., 200., 200);
        FillHist("Etamu1_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, MuColl.at(0).Eta(), weight, -5., 5., 100);
        FillHist("Etamu2_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, MuColl.at(1).Eta(), weight, -5., 5., 100);
  
        FillHist("Mmumu_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, Mmumu, weight, 0., 200., 200);
        FillHist("M3l_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, (EleColl.at(0)+MuColl.at(0)+MuColl.at(1)).M(), weight, 0., 500., 500);
        FillHist("Nj_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, JetVetoColl.size(), weight, 0., 10., 10);
        FillHist("Nb_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, BJetVetoColl.size(), weight, 0., 10., 10);
        FillHist("MET_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, MET, weight, 0., 200., 200);
        FillHist("MTW_3lOSOnZOffZG"+SystDirLabel+SystKindLabel, MTW, weight, 0., 200., 200);
      }
    }
    //3l OnZ & HasBJet
    if(HasBJet && OnZ){
      FillHist("PTe_3lOSOnZHasB"+SystDirLabel+SystKindLabel, EleColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_3lOSOnZHasB"+SystDirLabel+SystKindLabel, EleColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_3lOSOnZHasB"+SystDirLabel+SystKindLabel, MuColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_3lOSOnZHasB"+SystDirLabel+SystKindLabel, MuColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_3lOSOnZHasB"+SystDirLabel+SystKindLabel, MuColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_3lOSOnZHasB"+SystDirLabel+SystKindLabel, MuColl.at(1).Eta(), weight, -5., 5., 100);

      FillHist("Mmumu_3lOSOnZHasB"+SystDirLabel+SystKindLabel, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_3lOSOnZHasB"+SystDirLabel+SystKindLabel, (EleColl.at(0)+MuColl.at(0)+MuColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_3lOSOnZHasB"+SystDirLabel+SystKindLabel, JetVetoColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOSOnZHasB"+SystDirLabel+SystKindLabel, BJetVetoColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOSOnZHasB"+SystDirLabel+SystKindLabel, MET, weight, 0., 200., 200);
      FillHist("MTW_3lOSOnZHasB"+SystDirLabel+SystKindLabel, MTW, weight, 0., 200., 200);
    }
    if(WZSel){
      FillHist("PTe_WZSel"+SystDirLabel+SystKindLabel, EleColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_WZSel"+SystDirLabel+SystKindLabel, EleColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_WZSel"+SystDirLabel+SystKindLabel, MuColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_WZSel"+SystDirLabel+SystKindLabel, MuColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_WZSel"+SystDirLabel+SystKindLabel, MuColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_WZSel"+SystDirLabel+SystKindLabel, MuColl.at(1).Eta(), weight, -5., 5., 100);

      FillHist("Mmumu_WZSel"+SystDirLabel+SystKindLabel, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_WZSel"+SystDirLabel+SystKindLabel, (EleColl.at(0)+MuColl.at(0)+MuColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_WZSel"+SystDirLabel+SystKindLabel, JetVetoColl.size(), weight, 0., 10., 10);
      FillHist("Nb_WZSel"+SystDirLabel+SystKindLabel, BJetVetoColl.size(), weight, 0., 10., 10);
      FillHist("MET_WZSel"+SystDirLabel+SystKindLabel, MET, weight, 0., 200., 200);
      FillHist("MTW_WZSel"+SystDirLabel+SystKindLabel, MTW, weight, 0., 200., 200);
    }
    if(ZGSel){
      FillHist("PTe_ZGSel"+SystDirLabel+SystKindLabel, EleColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_ZGSel"+SystDirLabel+SystKindLabel, EleColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_ZGSel"+SystDirLabel+SystKindLabel, MuColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_ZGSel"+SystDirLabel+SystKindLabel, MuColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_ZGSel"+SystDirLabel+SystKindLabel, MuColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_ZGSel"+SystDirLabel+SystKindLabel, MuColl.at(1).Eta(), weight, -5., 5., 100);

      FillHist("Mmumu_ZGSel"+SystDirLabel+SystKindLabel, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_ZGSel"+SystDirLabel+SystKindLabel, (EleColl.at(0)+MuColl.at(0)+MuColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_ZGSel"+SystDirLabel+SystKindLabel, JetVetoColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ZGSel"+SystDirLabel+SystKindLabel, BJetVetoColl.size(), weight, 0., 10., 10);
      FillHist("MET_ZGSel"+SystDirLabel+SystKindLabel, MET, weight, 0., 200., 200);
      FillHist("MTW_ZGSel"+SystDirLabel+SystKindLabel, MTW, weight, 0., 200., 200);
    }
    if(ttZSel){
      FillHist("PTe_ttZSel"+SystDirLabel+SystKindLabel, EleColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_ttZSel"+SystDirLabel+SystKindLabel, EleColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_ttZSel"+SystDirLabel+SystKindLabel, MuColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_ttZSel"+SystDirLabel+SystKindLabel, MuColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_ttZSel"+SystDirLabel+SystKindLabel, MuColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_ttZSel"+SystDirLabel+SystKindLabel, MuColl.at(1).Eta(), weight, -5., 5., 100);

      FillHist("Mmumu_ttZSel"+SystDirLabel+SystKindLabel, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_ttZSel"+SystDirLabel+SystKindLabel, (EleColl.at(0)+MuColl.at(0)+MuColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_ttZSel"+SystDirLabel+SystKindLabel, JetVetoColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ttZSel"+SystDirLabel+SystKindLabel, BJetVetoColl.size(), weight, 0., 10., 10);
      FillHist("MET_ttZSel"+SystDirLabel+SystKindLabel, MET, weight, 0., 200., 200);
      FillHist("MTW_ttZSel"+SystDirLabel+SystKindLabel, MTW, weight, 0., 200., 200);
    }
    if(ZbbSel){
      FillHist("PTe_ZbbSel"+SystDirLabel+SystKindLabel, EleColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("Etae_ZbbSel"+SystDirLabel+SystKindLabel, EleColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("PTmu1_ZbbSel"+SystDirLabel+SystKindLabel, MuColl.at(0).Pt(), weight, 0., 200., 200);
      FillHist("PTmu2_ZbbSel"+SystDirLabel+SystKindLabel, MuColl.at(1).Pt(), weight, 0., 200., 200);
      FillHist("Etamu1_ZbbSel"+SystDirLabel+SystKindLabel, MuColl.at(0).Eta(), weight, -5., 5., 100);
      FillHist("Etamu2_ZbbSel"+SystDirLabel+SystKindLabel, MuColl.at(1).Eta(), weight, -5., 5., 100);

      FillHist("Mmumu_ZbbSel"+SystDirLabel+SystKindLabel, Mmumu, weight, 0., 200., 200);
      FillHist("M3l_ZbbSel"+SystDirLabel+SystKindLabel, (EleColl.at(0)+MuColl.at(0)+MuColl.at(1)).M(), weight, 0., 500., 500);
      FillHist("Nj_ZbbSel"+SystDirLabel+SystKindLabel, JetVetoColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ZbbSel"+SystDirLabel+SystKindLabel, BJetVetoColl.size(), weight, 0., 10., 10);
      FillHist("MET_ZbbSel"+SystDirLabel+SystKindLabel, MET, weight, 0., 200., 200);
      FillHist("MTW_ZbbSel"+SystDirLabel+SystKindLabel, MTW, weight, 0., 200., 200);
    }

  }//End of Closure Cycle

  return;
}

float Jul2017_DataFakeStudy::FakeRateData(snu::KElectron Ele, TString ID){

  float FR=0., PTCorr=Ele.Pt()*(1.+Ele.PFRelIso(0.3)), fEta=fabs(Ele.Eta());
  if(ID=="POGMVAMIPOpt2Iso04"){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.193686 ;
      else if(PTCorr<50 ) FR= 0.115564 ;
      else                FR= 0.11773  ;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.233713;
      else if(PTCorr<50 ) FR= 0.127104;
      else                FR= 0.139068;
    }
    else{
      if     (PTCorr<35 ) FR= 0.261925;
      else if(PTCorr<50 ) FR= 0.168822;
      else                FR= 0.193664;
    }
  }
  else if(ID=="POGMVAMTIsop06IPp025p05Sig4FakeLIsop4"){
    //OLD(MVA:(-0.9,-0.9,-0.8) ; Coarse Optimisation
//    if(fEta<0.8){
//      if     (PTCorr<35 ) FR= 0.122623 ;
//      else if(PTCorr<50 ) FR= 0.0679368;
//      else                FR= 0.0707609;
//    }
//    else if(fEta<1.479){
//      if     (PTCorr<35 ) FR= 0.12563  ;
//      else if(PTCorr<50 ) FR= 0.0646067;
//      else                FR= 0.0715022;
//    }
//    else{
//      if     (PTCorr<35 ) FR= 0.156385 ;
//      else if(PTCorr<50 ) FR= 0.0942703;
//      else                FR= 0.104913 ;
//    }
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.119907 ;
      else if(PTCorr<50 ) FR= 0.065943 ;
      else                FR= 0.0682414;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.132106 ;
      else if(PTCorr<50 ) FR= 0.0682482;
      else                FR= 0.077325 ;
    }
    else{
      if     (PTCorr<35 ) FR= 0.160416 ;
      else if(PTCorr<50 ) FR= 0.0976898;
      else                FR= 0.109458 ;
    }
  }
  else if(ID=="POGMVAMTIsop1IPp025p05Sig4FakeLIsop4"){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.204863;
      else if(PTCorr<50 ) FR= 0.120178;
      else                FR= 0.123658;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.236692;
      else if(PTCorr<50 ) FR= 0.128829;
      else                FR= 0.141241;
    }
    else{
      if     (PTCorr<35 ) FR= 0.294857;
      else if(PTCorr<50 ) FR= 0.193047;
      else                FR= 0.219773;
    }
  }
  else if(ID=="POGMVATTIsop06IPp025p05Sig4FakeLIsop4"){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.0813384;
      else if(PTCorr<50 ) FR= 0.0473129;
      else                FR= 0.0490471;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.0949602;
      else if(PTCorr<50 ) FR= 0.0517796;
      else                FR= 0.0630693;
    }
    else{
      if     (PTCorr<35 ) FR= 0.109435 ;
      else if(PTCorr<50 ) FR= 0.0670226;
      else                FR= 0.079485 ;
    }
  }
  else if(ID=="POGMVATTIsop1IPp025p05Sig4FakeLIsop4"){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.146427 ;
      else if(PTCorr<50 ) FR= 0.0908473;
      else                FR= 0.0960944;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.18927 ;
      else if(PTCorr<50 ) FR= 0.108929;
      else                FR= 0.126564;
    }
    else{
      if     (PTCorr<35 ) FR= 0.217664;
      else if(PTCorr<50 ) FR= 0.141156;
      else                FR= 0.173934;
    }
  }
  else if(ID=="POGMVAMIPOpt2Iso06"){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.193733;
      else if(PTCorr<50 ) FR= 0.103705;
      else                FR= 0.104228;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.233731;
      else if(PTCorr<50 ) FR= 0.115073;
      else                FR= 0.12269 ;
    }
    else{
      if     (PTCorr<35 ) FR= 0.261928;
      else if(PTCorr<50 ) FR= 0.15727 ;
      else                FR= 0.179709;
    }
  }
  else if(ID=="POGMVAMFakeLIso04NoMVA"){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.134663 ;
      else if(PTCorr<50 ) FR= 0.0726663;
      else                FR= 0.0687263;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.139805 ;
      else if(PTCorr<50 ) FR= 0.0696455;
      else                FR= 0.0688009;
    }
    else{
      if     (PTCorr<35 ) FR= 0.0977045;
      else if(PTCorr<50 ) FR= 0.0474763;
      else                FR= 0.0413024;
    }
  }
  else if(ID=="POGMVAMFakeLIso04Opt2Data"){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.169891 ;
      else if(PTCorr<50 ) FR= 0.0980297;
      else                FR= 0.0972565;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.196989;
      else if(PTCorr<50 ) FR= 0.103387;
      else                FR= 0.109215;
    }
    else{
      if     (PTCorr<35 ) FR= 0.1734  ;
      else if(PTCorr<50 ) FR= 0.102403;
      else                FR= 0.108606;
    }
  }
  else if(ID=="POGMVAMFakeLIso04Opt1"){
    if(fEta<0.8){
      if     (PTCorr<35 ) FR= 0.232943;
      else if(PTCorr<50 ) FR= 0.143792;
      else                FR= 0.154709;
    }
    else if(fEta<1.479){
      if     (PTCorr<35 ) FR= 0.277706;
      else if(PTCorr<50 ) FR= 0.156332;
      else                FR= 0.177498;
    }
    else{
      if     (PTCorr<35 ) FR= 0.336682;
      else if(PTCorr<50 ) FR= 0.226226;
      else                FR= 0.266821;
    }
  }
  return FR;
} 


float Jul2017_DataFakeStudy::GetPreTrigPURW(int Nvtx){

   float weight=1.;
   if(Nvtx<2) weight=1.22087;
   else if(Nvtx<3) weight=4.79564;
   else if(Nvtx<4) weight=3.84689;
   else if(Nvtx<5) weight=2.47179;
   else if(Nvtx<6) weight=3.74121;
   else if(Nvtx<7) weight=3.16705;
   else if(Nvtx<8) weight=3.04533;
   else if(Nvtx<9) weight=2.85291;
   else if(Nvtx<10) weight=2.50907;
   else if(Nvtx<11) weight=2.29848;
   else if(Nvtx<12) weight=2.25329;
   else if(Nvtx<13) weight=1.9864;
   else if(Nvtx<14) weight=1.73597;
   else if(Nvtx<15) weight=1.57862;
   else if(Nvtx<16) weight=1.38849;
   else if(Nvtx<17) weight=1.21789;
   else if(Nvtx<18) weight=1.05361;
   else if(Nvtx<19) weight=0.909483;
   else if(Nvtx<20) weight=0.787018;
   else if(Nvtx<21) weight=0.699577;
   else if(Nvtx<22) weight=0.606739;
   else if(Nvtx<23) weight=0.533779;
   else if(Nvtx<24) weight=0.47826;
   else if(Nvtx<25) weight=0.430963;
   else if(Nvtx<26) weight=0.38291;
   else if(Nvtx<27) weight=0.368938;
   else if(Nvtx<28) weight=0.306386;
   else if(Nvtx<29) weight=0.321641;
   else if(Nvtx<30) weight=0.267767;
   else if(Nvtx<31) weight=0.303042;
   else if(Nvtx<32) weight=0.220255;
   else if(Nvtx<33) weight=0.237769;
   else if(Nvtx<34) weight=0.214341;
   else if(Nvtx<35) weight=0.254868;
   else if(Nvtx<36) weight=0.318919;
   else if(Nvtx<37) weight=0.290442;
   else if(Nvtx<38) weight=0.204838;
   else if(Nvtx<39) weight=0.389061;
   else if(Nvtx<40) weight=0.309147;
   else if(Nvtx<41) weight=0.286387;
   else if(Nvtx<42) weight=0.266532;
   else if(Nvtx<43) weight=0.327612;
   else if(Nvtx<44) weight=0.575947;
   else if(Nvtx<45) weight=0.23432;
   else if(Nvtx<46) weight=0.572198;
   else if(Nvtx<47) weight=0.541966;
   else if(Nvtx<48) weight=0.384433;
   else if(Nvtx<49) weight=1.57278;
   else weight=0.389116;

   return weight;
}



void Jul2017_DataFakeStudy::FillCutFlow(TString cut, float weight){
  
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
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(5,"PreSel");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(6,"1e");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(7,"#geq1j");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(8,"1j");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(9,"met20");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
    if(!GetHist("cutflow_N")){
      AnalyzerCore::MakeHistograms("cutflow_N", 11, 0., 11.);
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(2,"TriggerCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(3,"EventCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(5,"PreSel");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(6,"1e");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(7,"#geq1j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(8,"1j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(9,"met20");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
  }
}



void Jul2017_DataFakeStudy::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Jul2017_DataFakeStudy::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  //Basic Hists/////////////////////////////////////////////////
  //Data properties before selection  
  AnalyzerCore::MakeHistograms("Basic_Nj_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_Nb_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NmuL_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NeL_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NmuT_orig",10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_NeT_orig", 10, 0., 10.);
  AnalyzerCore::MakeHistograms("Basic_METdist_orig", 100, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_NmuT_weT", 10, 0., 10.);

  AnalyzerCore::MakeHistograms("Basic_Ptmu1_3mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu2_3mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Ptmu1_1e2mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_Pte_1e2mu_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j1_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j2_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j3_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_j4_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_b1_Et_orig", 200, 0., 200.);
  AnalyzerCore::MakeHistograms("Basic_b2_Et_orig", 200, 0., 200.);

  //After Nlcut
  AnalyzerCore::MakeHistograms("Basic_Nvtx_NoRW_wNlOScut", 50, 0., 50.);
  AnalyzerCore::MakeHistograms("Basic_Nvtx_PURW_wNlOScut", 50, 0., 50.);
  AnalyzerCore::MakeHistograms("Basic_MET_wNlOScut", 200, 0., 200.);
 


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Jul2017_DataFakeStudyCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Jul2017_DataFakeStudy::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
