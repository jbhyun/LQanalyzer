/***************************************************************************
 * @Project: Mar2018_TrilepCR 
 * @Package: LQanalyzer, ROOT-based analysis framework for Korea SNU / Main package author: John Almond
 *
 * @Code Author: Jihwan Bhyun 
 *
 ***************************************************************************/

/// Local includes
 #include "Mar2018_TrilepCR.h"

 //Core includes
 #include "Reweight.h"
 #include "EventBase.h"                                                                                                                           
 #include "BaseSelection.h"
 #include "AnalyzerCore.h"
 //// Needed to allow inheritance for use in LQCore/core classes
 ClassImp (Mar2018_TrilepCR);

 Mar2018_TrilepCR::Mar2018_TrilepCR() : AnalyzerCore(), out_muons(0) {

   SetLogName("Mar2018_TrilepCR");
   Message("In Mar2018_TrilepCR constructor", INFO);
   InitialiseAnalysis();
 }


 void Mar2018_TrilepCR::InitialiseAnalysis() throw( LQError ) {
   
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

void Mar2018_TrilepCR::ExecuteEvents()throw( LQError ){

////Basic Infos///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Apply the gen weight 
   float geneff_weight=1., gennorm_weight=1., k_factor_weight=1.;
   float pileup_reweight=1., pileup_reweight_systdown=1., pileup_reweight_systup=1.;
   if(!isData){
     weight*=MCweight;
     geneff_weight            = GenFilterEfficiency(k_sample_name);
     gennorm_weight           = SignalNorm(k_sample_name, 200.);
     k_factor_weight          = GetKFactor();
     pileup_reweight          = eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::central);                                                    
     pileup_reweight_systup   = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::up)  /pileup_reweight :0. ;
     pileup_reweight_systdown = pileup_reweight>0.? eventbase->GetEvent().PileUpWeight_Gold(snu::KEvent::down)/pileup_reweight :0. ;
    
     weight *= geneff_weight*gennorm_weight*k_factor_weight*pileup_reweight;
   }

   FillHist("GenWeight", MCweight, 1, -2., 2., 4);

   /// Acts on data to remove bad reconstructed event 
   //if(isData&& (! eventbase->GetEvent().LumiMask(lumimask))) return;

   m_logger<<DEBUG<<"RunNumber/Event Number = "<<eventbase->GetEvent().RunNumber()<<" : "<<eventbase->GetEvent().EventNumber()<<LQLogger::endmsg;
   m_logger<<DEBUG<<"isData = "<<isData<<LQLogger::endmsg;


   //Numbet of Vertex NoPUreweight
   FillHist("Nvtx_nocut_NoRW", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);

   //Numbet of Vertex PUreweight
   FillHist("Nvtx_nocut_PURW", eventbase->GetEvent().nVertices(), weight*pileup_reweight, 0., 50., 50);

   //Total Event(MCweight+PUrw)////////////////
   FillCutFlow("NoCut", weight*pileup_reweight);

   bool TrilepCR=false, SSDilepCR=false, TTZAnomaly=false, SystRun=false;
   bool DoubleMuon=false, ElectronMuon=false, TriMu=false, EMuMu=false;
   TString Cycle="";
   for(int i=0; i<k_flags.size(); i++){
     if     (k_flags.at(i).Contains("TrilepCR"))      {TrilepCR       = true; Cycle="TrilepCR";}
     else if(k_flags.at(i).Contains("SSDilepCR"))     {SSDilepCR      = true; Cycle="SSDilepCR";}
     else if(k_flags.at(i).Contains("TTZAnomaly"))    {TTZAnomaly     = true; Cycle="TTZAnomaly";}
     else if(k_flags.at(i).Contains("DoubleMuon"))     DoubleMuon     = true;
     else if(k_flags.at(i).Contains("ElectronMuon"))   ElectronMuon   = true;
     else if(k_flags.at(i).Contains("TriMu"))          TriMu          = true;
     else if(k_flags.at(i).Contains("EMuMu"))          EMuMu          = true;
     else if(k_flags.at(i).Contains("SystRun"))        SystRun        = true;
   }

    
   /***************************************************************************************
   **Trigger Treatment
   ***************************************************************************************/
   //ListTriggersAvailable();

   //Normalisation Lumi Setting
   bool Pass_Trigger=false, Pass_TriggerBG=false, Pass_TriggerH=false;
   float trigger_ps_weight=1., trigger_period_weight=1.;
   float LumiBG=27.257618, LumiH=8.605696, LumiBH=35.863314;

   if( TrilepCR || SSDilepCR || TTZAnomaly ){
     if(DoubleMuon || TriMu){
       if(  PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v")
          ||PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") )    Pass_TriggerBG=true;
       if(  PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")
          ||PassTrigger("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") ) Pass_TriggerH =true;
  
       if(isData){
         int DataPeriod=GetDataPeriod();
         if( DataPeriod>0 && DataPeriod<7 && Pass_TriggerBG ) Pass_Trigger=true;
         else if( DataPeriod==7 && Pass_TriggerH )            Pass_Trigger=true;
       }
       else{
         if( Pass_TriggerBG || Pass_TriggerH ) Pass_Trigger=true;
         trigger_period_weight=( (Pass_TriggerBG? 27.257618:0.)+(Pass_TriggerH? 8.605696:0.) )/35.863314;
         trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
       }
     }
     if(ElectronMuon || EMuMu){
       if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") )    Pass_TriggerBG=true;
       if( PassTrigger("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") ) Pass_TriggerH =true;
  
       if(isData){
         int DataPeriod=GetDataPeriod();
         if( DataPeriod>0 && DataPeriod<7 && Pass_TriggerBG ) Pass_Trigger=true;
         else if( DataPeriod==7 && Pass_TriggerH )            Pass_Trigger=true;
       }
       else{
         if( Pass_TriggerBG || Pass_TriggerH ) Pass_Trigger=true;
         trigger_period_weight=( (Pass_TriggerBG? 27.257618:0.)+(Pass_TriggerH? 8.605696:0.) )/35.863314;
         trigger_ps_weight=WeightByTrigger("HLT_IsoMu24_v", TargetLumi);
       }
     }
   }
   FillHist("TriggerPSWeight", trigger_ps_weight, 1., 0., 1., 100);
   weight*=trigger_ps_weight*trigger_period_weight;

   //Trigger Cut
   if(!Pass_Trigger) return;
   FillCutFlow("Trigger", weight*pileup_reweight);

   /**********************************************************************************/

   //METFilterCut : https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFilters
   if(!PassMETFilter()) return;  FillCutFlow("METFilter", weight*pileup_reweight);

   //Vertex Cut
   //(vtx.ndof>4&&maxAbsZ<=0)||abs(vtx.z)<= 24)&&((maxd0 <=0)||abs(vtx.position.rho)<=2)&&!(vtx.isFake))
   if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; FillCutFlow("VertexCut", weight*pileup_reweight);
   


///////Objects in Analysis/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   //Primary Object Collection
   std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  
   /**PreSelCut***********************************************************************************************/
   //Intended for Code speed boosting up.
     eventbase->GetMuonSel()->SetPt(5.);                     eventbase->GetMuonSel()->SetEta(2.4);
   std::vector<snu::KMuon> muonPreColl; eventbase->GetMuonSel()->Selection(muonPreColl, true);
     eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
   std::vector<snu::KElectron> electronPreColl; eventbase->GetElectronSel()->Selection(electronPreColl);
   if     (TrilepCR ){ if( !(muonPreColl.size()>=3 || (muonPreColl.size()>=2 && electronPreColl.size()>=1)) ) return; }
   else if(SSDilepCR){ bool Pass=false;
     int Nlep=muonPreColl.size()+electronPreColl.size();
     if(Nlep>=3) Pass=true;
     else if(Nlep==2){ int SumCharge=0; for(int i=0; i<(int) muonPreColl.size(); i++){ SumCharge+=muonPreColl.at(i).Charge();}
                                        for(int i=0; i<(int) electronPreColl.size(); i++){ SumCharge+=electronPreColl.at(i).Charge(); }  
                       if(SumCharge!=0) Pass=true; }
     if(!Pass) return;
   }
   /**********************************************************************************************************/


   FillCutFlow("PreSel", weight*pileup_reweight);
   //For Fake Study
   //Muon ID's to Test
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.2);                 eventbase->GetMuonSel()->SetdxySigMax(4.);
     eventbase->GetMuonSel()->SetBSdz(0.1);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.6);
   std::vector<snu::KMuon> muonLooseColl; eventbase->GetMuonSel()->Selection(muonLooseColl, true);
     eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
     eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
     eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetdxySigMax(4.);
     eventbase->GetMuonSel()->SetBSdz(0.05);
     eventbase->GetMuonSel()->SetChiNdof(4.);
     eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.2);
   std::vector<snu::KMuon> muonTightColl; eventbase->GetMuonSel()->Selection(muonTightColl,true);

   std::vector<snu::KMuon> muonColl;
     if(!k_running_nonprompt){ muonColl=muonTightColl; }else{ muonColl=muonLooseColl; }


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
     eventbase->GetElectronSel()->SetPt(25.);                eventbase->GetElectronSel()->SetEta(2.5);
     eventbase->GetElectronSel()->SetBETrRegIncl(false);
     eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
     eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
     eventbase->GetElectronSel()->SetdxySigMax(4.);
     eventbase->GetElectronSel()->SetApplyConvVeto(true);
   std::vector<snu::KElectron> electronTightColl; eventbase->GetElectronSel()->Selection(electronTightColl);
   std::vector<snu::KElectron> electronColl;
     if(!k_running_nonprompt){ electronColl=electronTightColl; }else{ electronColl=electronLooseColl; }

     eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
     eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     bool LeptonVeto=false;
   std::vector<snu::KJet> jetNoVetoColl; eventbase->GetJetSel()->Selection(jetNoVetoColl, LeptonVeto, muonLooseColl, electronLooseColl);
   std::vector<snu::KJet> bjetNoVetoColl = SelBJets(jetNoVetoColl, "Medium");

   std::vector<snu::KJet> jetColl  = SkimJetColl(jetNoVetoColl,  electronLooseColl, muonLooseColl, "EleMuVeto");
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
   float nvtx_reweight=1.; 

   /*This part is for boosting up speed.. SF part takes rather longer time than expected*/
   if(!SystRun){
     if     (TrilepCR    ){ if(muonLooseColl.size()>=3 or (muonLooseColl.size()>=2 && electronLooseColl.size()>=1)) EventCand=true; }
     else if(SSDilepCR   ){ if(muonLooseColl.size()>=2 or (muonLooseColl.size()>=1 && electronLooseColl.size()>=1)) EventCand=true; }

     if(!isData){
       if(EventCand && ( TrilepCR || SSDilepCR )){
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
      
         //trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl); - No longer recommendation
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);

         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);
         if(DoubleMuon || TriMu){
           float trigger_sf1 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v");
           float trigger_sf2 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v");
           trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         }
         if(ElectronMuon || EMuMu){
           float trigger_sf1 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
           float trigger_sf2 = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
           trigger_sf    = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         }
       }
     }
     else if(k_running_nonprompt && EventCand ){
       fake_weight = GetFakeWeight_Data(muonLooseColl, electronLooseColl,
                     "POGTIsop6IPp2p1sig4" , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "TrkIsoVVLConeSUSY");
     }
   }
   weight *= id_weight_ele*reco_weight_ele*id_weight_mu*iso_weight_mu*trk_weight_mu*trigger_sf*fake_weight*btag_sf;
   /***************************************************************************************************/



   //=============================================================================//
   // Main Analysis Code     -----------------------------------------------------//
   //=============================================================================//


   if(SystRun){
     // Syst Sel. and Syst Corr. -------------------------------------------------------------------//
       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
       eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
       eventbase->GetMuonSel()->SetdxySigMax(4.);
       eventbase->GetMuonSel()->SetChiNdof(4.);
     std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,true);
       SetMuonResCorrection(MuTMuEnUpColl, "SystUp");

       eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
       eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
       eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
       eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
       eventbase->GetMuonSel()->SetdxySigMax(4.);
       eventbase->GetMuonSel()->SetChiNdof(4.);
     std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,true);
       SetMuonResCorrection(MuTMuEnDownColl, "SystDown");


       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
       eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdxySigMax(4.);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");

       eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
       eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
       eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
       eventbase->GetElectronSel()->SetBETrRegIncl(false);
       eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
       eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
       eventbase->GetElectronSel()->SetdxySigMax(4.);
       eventbase->GetElectronSel()->SetApplyConvVeto(true);
     std::vector<snu::KElectron> EleTElEnDownColl; eventbase->GetElectronSel()->Selection(EleTElEnDownColl, "SystDownElEn");

     LeptonVeto=true;
       eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
       eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     std::vector<snu::KJet> jetJESUpColl; eventbase->GetJetSel()->Selection(jetJESUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJES");
       eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
       eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     std::vector<snu::KJet> jetJERUpColl; eventbase->GetJetSel()->Selection(jetJERUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJER");
       eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
       eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     std::vector<snu::KJet> jetJESDownColl; eventbase->GetJetSel()->Selection(jetJESDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJES");
       eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
       eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
     std::vector<snu::KJet> jetJERDownColl; eventbase->GetJetSel()->Selection(jetJERDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJER");

     std::vector<snu::KJet> bjetJESUpColl   = SelBJets(jetJESUpColl  , "Medium");
     std::vector<snu::KJet> bjetJERUpColl   = SelBJets(jetJERUpColl  , "Medium");
     std::vector<snu::KJet> bjetJESDownColl = SelBJets(jetJESDownColl, "Medium");
     std::vector<snu::KJet> bjetJERDownColl = SelBJets(jetJERDownColl, "Medium");

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

     //Scale Factors--------------------------------------------------------------------------------//      
     float id_weight_ele_sfup=1.    , id_weight_ele_sfdown=1.    , id_weight_mu_sfup=1.    , id_weight_mu_sfdown=1.;
     float id_weight_ele_ElEnup=1.  , reco_weight_ele_ElEnup=1.  , id_weight_mu_MuEnup=1.  , trk_weight_mu_MuEnup=1.;
     float id_weight_ele_ElEndown=1., reco_weight_ele_ElEndown=1., id_weight_mu_MuEndown=1., trk_weight_mu_MuEndown=1.;
     float btag_sf_JESup=1.  , btag_sf_JERup=1.  , btag_sf_LTagup=1.  , btag_sf_BCTagup=1.;
     float btag_sf_JESdown=1., btag_sf_JERdown=1., btag_sf_LTagdown=1., btag_sf_BCTagdown=1.;
     float fake_weight_FRup=1., fake_weight_FRdown=1.;
     float ZG_SF=1., conv_up=1., conv_down=1.;
     float xsec_up=1., xsec_down=1.;
     float trigger_sf=1., trigger_sf_up=1., trigger_sf_down=1.;
     if(!isData){
       id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
       reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
       id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
       //trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl); - No longer recommendation
       btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

       id_weight_ele_sfup       = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl,  1);
       id_weight_ele_sfdown     = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl, -1);
       id_weight_mu_sfup        = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl,  1);
       id_weight_mu_sfdown      = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl, -1);

       id_weight_ele_ElEnup     = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnUpColl);
       reco_weight_ele_ElEnup   = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
       id_weight_mu_MuEnup      = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnUpColl);
       //trk_weight_mu_MuEnup     = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnUpColl); - No longer recommendation

       id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnDownColl);
       reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);
       id_weight_mu_MuEndown    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnDownColl);
       //trk_weight_mu_MuEndown   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnDownColl); - No longer recommendation

       if(TriMu){
         float trigger_sf1          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v");
         float trigger_sf2          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v");
         float trigger_sf1_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg1SystUp");
         float trigger_sf2_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg1SystUp");
         float trigger_sf1_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg2SystUp");
         float trigger_sf2_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg2SystUp");
         float trigger_sf1_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg1SystDown");
         float trigger_sf2_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg1SystDown");
         float trigger_sf1_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg2SystDown");
         float trigger_sf2_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg2SystDown");
         float trigger_sf2_dzup     = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","DZSystUp");
         float trigger_sf2_dzdown   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","DZSystDown");
         float trigger_sf1_up   = trigger_sf1+sqrt(pow(trigger_sf1_leg1up  -trigger_sf1,2)+pow(trigger_sf1_leg2up  -trigger_sf1,2));
         float trigger_sf1_down = trigger_sf1-sqrt(pow(trigger_sf1_leg1down-trigger_sf1,2)+pow(trigger_sf1_leg2down-trigger_sf1,2));
         float trigger_sf2_up   = trigger_sf2+sqrt(pow(trigger_sf2_leg1up  -trigger_sf2,2)+pow(trigger_sf2_leg2up  -trigger_sf2,2)+pow(trigger_sf2_dzup  -trigger_sf2,2));
         float trigger_sf2_down = trigger_sf2-sqrt(pow(trigger_sf2_leg1down-trigger_sf2,2)+pow(trigger_sf2_leg2down-trigger_sf2,2)+pow(trigger_sf2_dzdown-trigger_sf2,2));

         trigger_sf      = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         trigger_sf_up   = ((Pass_TriggerBG ? trigger_sf1_up:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_up:0.)*LumiH)/LumiBH;
         trigger_sf_down = ((Pass_TriggerBG ? trigger_sf1_down:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_down:0.)*LumiH)/LumiBH;
       }
       else if(EMuMu){
         float trigger_sf1          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
         float trigger_sf2          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
         float trigger_sf1_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg1SystUp");
         float trigger_sf2_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg1SystUp");
         float trigger_sf1_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg2SystUp");
         float trigger_sf2_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg2SystUp");
         float trigger_sf1_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg1SystDown");
         float trigger_sf2_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg1SystDown");
         float trigger_sf1_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg2SystDown");
         float trigger_sf2_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg2SystDown");
         float trigger_sf2_dzup     = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","DZSystUp");
         float trigger_sf2_dzdown   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","DZSystDown");
         float trigger_sf1_up   = trigger_sf1+sqrt(pow(trigger_sf1_leg1up  -trigger_sf1,2)+pow(trigger_sf1_leg2up  -trigger_sf1,2));
         float trigger_sf1_down = trigger_sf1-sqrt(pow(trigger_sf1_leg1down-trigger_sf1,2)+pow(trigger_sf1_leg2down-trigger_sf1,2));
         float trigger_sf2_up   = trigger_sf2+sqrt(pow(trigger_sf2_leg1up  -trigger_sf2,2)+pow(trigger_sf2_leg2up  -trigger_sf2,2)+pow(trigger_sf2_dzup  -trigger_sf2,2));
         float trigger_sf2_down = trigger_sf2-sqrt(pow(trigger_sf2_leg1down-trigger_sf2,2)+pow(trigger_sf2_leg2down-trigger_sf2,2)+pow(trigger_sf2_dzdown-trigger_sf2,2));

         trigger_sf      = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
         trigger_sf_up   = ((Pass_TriggerBG ? trigger_sf1_up:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_up:0.)*LumiH)/LumiBH;
         trigger_sf_down = ((Pass_TriggerBG ? trigger_sf1_down:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_down:0.)*LumiH)/LumiBH;
       }

       btag_sf_LTagup = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpLTag");
       btag_sf_BCTagup= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpBCTag");
       btag_sf_JESup  = BTagScaleFactor_1a(jetJESUpColl, snu::KJet::CSVv2, snu::KJet::Medium);
       btag_sf_JERup  = BTagScaleFactor_1a(jetJERUpColl, snu::KJet::CSVv2, snu::KJet::Medium);

       btag_sf_LTagdown = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownLTag");
       btag_sf_BCTagdown= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownBCTag");
       btag_sf_JESdown  = BTagScaleFactor_1a(jetJESDownColl, snu::KJet::CSVv2, snu::KJet::Medium);
       btag_sf_JERdown  = BTagScaleFactor_1a(jetJERDownColl, snu::KJet::CSVv2, snu::KJet::Medium);

       xsec_up  +=GetXsecUncertainty(k_sample_name);
       xsec_down-=GetXsecUncertainty(k_sample_name);
       if(TriMu){
         if(k_sample_name.Contains("ZGto2LG")){ ZG_SF=0.845; conv_up+=0.15; conv_down-=0.15; }
         else if(k_sample_name.Contains("TTG")){ conv_up+=0.5; conv_down-=0.5; }
       }
       if(EMuMu){
         if(k_sample_name.Contains("ZGto2LG")){ ZG_SF=0.961; conv_up+=0.092; conv_down-=0.092; }
         else if(k_sample_name.Contains("TTG")){ conv_up+=0.5; conv_down-=0.5; }
       }
     }
     else if(k_running_nonprompt){
       fake_weight = GetFakeWeight_Data(muonLooseColl, electronLooseColl, "POGTIsop6IPp2p1sig4" , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "TrkIsoVVLConeSUSY");
     }
     weight*=ZG_SF;
       
     float systweight_central=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
                                                                                                                                                        
     float systweight_Trigup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf_up*fake_weight;
     float systweight_ElIDup =weight*id_weight_ele_sfup  *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
     float systweight_MuIDup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_sfup  *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
     float systweight_ElEnup =weight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup*id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
     float systweight_MuEnup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_MuEnup*trk_weight_mu_MuEnup*btag_sf        *trigger_sf   *fake_weight;
     float systweight_JESup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JESup  *trigger_sf   *fake_weight;
     float systweight_JERup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JERup  *trigger_sf   *fake_weight;
     float systweight_LTagup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_LTagup *trigger_sf   *fake_weight;
     float systweight_BCTagup=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_BCTagup*trigger_sf   *fake_weight;
     float systweight_PUup   =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight*pileup_reweight_systup;
     float systweight_Xsecup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight*xsec_up;
     float systweight_Convup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight*conv_up;

     float systweight_Trigdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf_down*fake_weight;
     float systweight_ElIDdown =weight*id_weight_ele_sfdown  *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight;
     float systweight_MuIDdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_sfdown  *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight;
     float systweight_ElEndown =weight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown*id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight;
     float systweight_MuEndown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_MuEndown*trk_weight_mu_MuEndown*btag_sf          *trigger_sf     *fake_weight;
     float systweight_JESdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JESdown  *trigger_sf     *fake_weight;
     float systweight_JERdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JERdown  *trigger_sf     *fake_weight;
     float systweight_LTagdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_LTagdown *trigger_sf     *fake_weight;
     float systweight_BCTagdown=weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_BCTagdown*trigger_sf     *fake_weight;
     float systweight_PUdown   =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight*pileup_reweight_systdown;
     float systweight_Xsecdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight*xsec_down;
     float systweight_Convdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight*conv_down;


     TString ChannelString = EMuMu? "EMuMu":"TriMu";
     //----------------------------------------------------------------------------------------------------------------------//


     DoSystRun(Cycle, ChannelString+"",
               electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
               systweight_central);
     if(!isData){
       DoSystRun(Cycle, ChannelString+"SystUpElID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_ElIDup);
       DoSystRun(Cycle, ChannelString+"SystUpMuID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_MuIDup);
       DoSystRun(Cycle, ChannelString+"SystUpTrig",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Trigup);
       DoSystRun(Cycle, ChannelString+"SystUpPU",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_PUup);
       DoSystRun(Cycle, ChannelString+"SystUpUncl",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                 systweight_central);
       DoSystRun(Cycle, ChannelString+"SystUpJES",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetJESUpColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                 systweight_JESup);
       DoSystRun(Cycle, ChannelString+"SystUpJER",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetJERUpColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                 systweight_JERup);
       DoSystRun(Cycle, ChannelString+"SystUpBTag_L",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_LTagup);
       DoSystRun(Cycle, ChannelString+"SystUpBTag_BC",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_BCTagup);
       DoSystRun(Cycle, ChannelString+"SystUpElEn",
                 EleTElEnUpColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                 systweight_ElEnup);
       DoSystRun(Cycle, ChannelString+"SystUpMuEn",
                 electronColl, electronLooseColl, MuTMuEnUpColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                 systweight_MuEnup);
       DoSystRun(Cycle, ChannelString+"SystUpXsec",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Xsecup);
       DoSystRun(Cycle, ChannelString+"SystUpConv",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Convup);



       DoSystRun(Cycle, ChannelString+"SystDownElID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_ElIDdown);
       DoSystRun(Cycle, ChannelString+"SystDownMuID",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_MuIDdown);
       DoSystRun(Cycle, ChannelString+"SystDownTrig",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Trigdown);
       DoSystRun(Cycle, ChannelString+"SystDownPU",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_PUdown);
       DoSystRun(Cycle, ChannelString+"SystDownUncl",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                 systweight_central);
       DoSystRun(Cycle, ChannelString+"SystDownJES",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetJESDownColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                 systweight_JESdown);
       DoSystRun(Cycle, ChannelString+"SystDownJER",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetJERDownColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                 systweight_JERdown);
       DoSystRun(Cycle, ChannelString+"SystDownBTag_L",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_LTagdown);
       DoSystRun(Cycle, ChannelString+"SystDownBTag_BC",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_BCTagdown);
       DoSystRun(Cycle, ChannelString+"SystDownElEn",
                 EleTElEnDownColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                 systweight_ElEndown);
       DoSystRun(Cycle, ChannelString+"SystDownMuEn",
                 electronColl, electronLooseColl, MuTMuEnDownColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                 systweight_MuEndown);
       DoSystRun(Cycle, ChannelString+"SystDownXsec",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Xsecdown);
       DoSystRun(Cycle, ChannelString+"SystDownConv",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_Convdown);
     }
   }

   if(SSDilepCR && !SystRun){

//     if(!isData && (k_sample_name.Contains("TT_powheg") || k_sample_name.Contains("DYJets") || k_sample_name.Contains("WJets"))){
//       int NLepCand=0;
//       for(int i=0; i<(int) muonColl.size(); i++){
//         int LepType=GetLeptonType(muonColl.at(i), truthColl);
//         if(LepType>0 || LepType<-4) NLepCand++;
//       }
//       for(int i=0; i<(int) electronColl.size(); i++){
//         int LepType=GetLeptonType(electronColl.at(i), truthColl);
//         if(LepType>0 || LepType<-4) NLepCand++;
//       }
//       if(NLepCand!=2) return;
//     }
   
     CheckSSDilepCRs(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
                     weight, "", "");
   }
   if(TrilepCR && !SystRun){
     if(!SystRun){
       CheckTrilepCRs(muonColl, muonLooseColl, electronColl, electronLooseColl, jetNoVetoColl, bjetNoVetoColl, met, met*cos(metphi), met*sin(metphi),
                      weight*fake_weight, "", "");
     }
     else{
       // Syst Sel. and Syst Corr. ----------------------------------------------------------------------------------------------------//
         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
         eventbase->GetMuonSel()->SetdxySigMax(4.);
         eventbase->GetMuonSel()->SetChiNdof(4.);
       std::vector<snu::KMuon> MuTMuEnUpColl; eventbase->GetMuonSel()->Selection(MuTMuEnUpColl,true);
         SetMuonResCorrection(MuTMuEnUpColl, "SystUp");

         eventbase->GetMuonSel()->SetID(BaseSelection::MUON_POG_TIGHT);
         eventbase->GetMuonSel()->SetPt(10.);                    eventbase->GetMuonSel()->SetEta(2.4);
         eventbase->GetMuonSel()->SetRelIsoType("PFRelIso04");   eventbase->GetMuonSel()->SetRelIso(0.20);
         eventbase->GetMuonSel()->SetBSdxy(0.01);                eventbase->GetMuonSel()->SetBSdz(0.05);
         eventbase->GetMuonSel()->SetdxySigMax(4.);
         eventbase->GetMuonSel()->SetChiNdof(4.);
       std::vector<snu::KMuon> MuTMuEnDownColl; eventbase->GetMuonSel()->Selection(MuTMuEnDownColl,true);
         SetMuonResCorrection(MuTMuEnDownColl, "SystDown");


         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025); eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
         eventbase->GetElectronSel()->SetdxySigMax(4.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnUpColl; eventbase->GetElectronSel()->Selection(EleTElEnUpColl, "SystUpElEn");

         eventbase->GetElectronSel()->SetID(BaseSelection::ELECTRON_POG_MVA_WP90);
         eventbase->GetElectronSel()->SetHLTSafeCut("CaloIdL_TrackIdL_IsoVL");
         eventbase->GetElectronSel()->SetPt(10.);                eventbase->GetElectronSel()->SetEta(2.5);
         eventbase->GetElectronSel()->SetBETrRegIncl(false);
         eventbase->GetElectronSel()->SetRelIsoType("Default");  eventbase->GetElectronSel()->SetRelIsoBEMax(0.06, 0.06);
         eventbase->GetElectronSel()->SetdxyBEMax(0.025, 0.025);   eventbase->GetElectronSel()->SetdzBEMax(0.1, 0.1);
         eventbase->GetElectronSel()->SetdxySigMax(4.);
         eventbase->GetElectronSel()->SetApplyConvVeto(true);
       std::vector<snu::KElectron> EleTElEnDownColl; eventbase->GetElectronSel()->Selection(EleTElEnDownColl, "SystDownElEn");

       LeptonVeto=true;
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESUpColl; eventbase->GetJetSel()->Selection(jetJESUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERUpColl; eventbase->GetJetSel()->Selection(jetJERUpColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystUpJER");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJESDownColl; eventbase->GetJetSel()->Selection(jetJESDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJES");
         eventbase->GetJetSel()->SetID(BaseSelection::PFJET_LOOSE);
         eventbase->GetJetSel()->SetPt(25.);                     eventbase->GetJetSel()->SetEta(2.4);
       std::vector<snu::KJet> jetJERDownColl; eventbase->GetJetSel()->Selection(jetJERDownColl, LeptonVeto, muonLooseColl, electronLooseColl, "SystDownJER");

       std::vector<snu::KJet> bjetJESUpColl   = SelBJets(jetJESUpColl  , "Medium");
       std::vector<snu::KJet> bjetJERUpColl   = SelBJets(jetJERUpColl  , "Medium");
       std::vector<snu::KJet> bjetJESDownColl = SelBJets(jetJESDownColl, "Medium");
       std::vector<snu::KJet> bjetJERDownColl = SelBJets(jetJERDownColl, "Medium");

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

       //Scale Factors--------------------------------------------------------------------------------//      
       float systweight=weight; //Lumi weight+Prescale weight+Period weight+PU weight applied by default;
       float id_weight_ele_ElEnup=1.  , reco_weight_ele_ElEnup=1.  , id_weight_mu_MuEnup=1.  , trk_weight_mu_MuEnup=1.;
       float id_weight_ele_ElEndown=1., reco_weight_ele_ElEndown=1., id_weight_mu_MuEndown=1., trk_weight_mu_MuEndown=1.;
       float btag_sf_JESup=1.  , btag_sf_JERup=1.  , btag_sf_LTagup=1.  , btag_sf_BCTagup=1.;
       float btag_sf_JESdown=1., btag_sf_JERdown=1., btag_sf_LTagdown=1., btag_sf_BCTagdown=1.;
       float fake_weight_FRup=1., fake_weight_FRdown=1.;
       float ZG_SF=1.;
       float trigger_sf=1., trigger_sf_up=1., trigger_sf_down=1.;
       if(!isData){
         id_weight_ele   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", electronColl);
         reco_weight_ele = mcdata_correction->ElectronRecoScaleFactor(electronColl);
         id_weight_mu    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", muonColl);
         trk_weight_mu   = mcdata_correction->MuonTrackingEffScaleFactor(muonColl);
         btag_sf         = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium);

         id_weight_ele_ElEnup   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnUpColl);
         reco_weight_ele_ElEnup = mcdata_correction->ElectronRecoScaleFactor(EleTElEnUpColl);
         id_weight_mu_MuEnup    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnUpColl);
         trk_weight_mu_MuEnup   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnUpColl);

         id_weight_ele_ElEndown   = mcdata_correction->ElectronScaleFactor("ELECTRON_HctoWA_TIGHT", EleTElEnDownColl);
         reco_weight_ele_ElEndown = mcdata_correction->ElectronRecoScaleFactor(EleTElEnDownColl);
         id_weight_mu_MuEndown    = mcdata_correction->MuonScaleFactor("MUON_HctoWA_TIGHT", MuTMuEnDownColl);
         trk_weight_mu_MuEndown   = mcdata_correction->MuonTrackingEffScaleFactor(MuTMuEnDownColl);

         if(TriMu){
           float trigger_sf1          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v");
           float trigger_sf2          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v");
           float trigger_sf1_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg1SystUp");
           float trigger_sf2_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg1SystUp");
           float trigger_sf1_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg2SystUp");
           float trigger_sf2_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg2SystUp");
           float trigger_sf1_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg1SystDown");
           float trigger_sf2_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg1SystDown");
           float trigger_sf1_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_v"   ,"Leg2SystDown");
           float trigger_sf2_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","Leg2SystDown");
           float trigger_sf2_dzup     = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","DZSystUp");
           float trigger_sf2_dzdown   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu17_TrkIsoVVL_Mu8ORTkMu8_TrkIsoVVL_DZ_v","DZSystDown");
           float trigger_sf1_up   = trigger_sf1+sqrt(pow(trigger_sf1_leg1up  -trigger_sf1,2)+pow(trigger_sf1_leg2up  -trigger_sf1,2));
           float trigger_sf1_down = trigger_sf1-sqrt(pow(trigger_sf1_leg1down-trigger_sf1,2)+pow(trigger_sf1_leg2down-trigger_sf1,2));
           float trigger_sf2_up   = trigger_sf2+sqrt(pow(trigger_sf2_leg1up  -trigger_sf2,2)+pow(trigger_sf2_leg2up  -trigger_sf2,2)+pow(trigger_sf2_dzup  -trigger_sf2,2));
           float trigger_sf2_down = trigger_sf2-sqrt(pow(trigger_sf2_leg1down-trigger_sf2,2)+pow(trigger_sf2_leg2down-trigger_sf2,2)+pow(trigger_sf2_dzdown-trigger_sf2,2));

           trigger_sf      = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
           trigger_sf_up   = ((Pass_TriggerBG ? trigger_sf1_up:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_up:0.)*LumiH)/LumiBH;
           trigger_sf_down = ((Pass_TriggerBG ? trigger_sf1_down:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_down:0.)*LumiH)/LumiBH;
         }
         else if(EMuMu){
           float trigger_sf1          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
           float trigger_sf2          = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
           float trigger_sf1_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg1SystUp");
           float trigger_sf2_leg1up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg1SystUp");
           float trigger_sf1_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg2SystUp");
           float trigger_sf2_leg2up   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg2SystUp");
           float trigger_sf1_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg1SystDown");
           float trigger_sf2_leg1down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg1SystDown");
           float trigger_sf1_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"   ,"Leg2SystDown");
           float trigger_sf2_leg2down = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","Leg2SystDown");
           float trigger_sf2_dzup     = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","DZSystUp");
           float trigger_sf2_dzdown   = mcdata_correction->GetTriggerSF(electronColl, muonColl, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","DZSystDown");
           float trigger_sf1_up   = trigger_sf1+sqrt(pow(trigger_sf1_leg1up  -trigger_sf1,2)+pow(trigger_sf1_leg2up  -trigger_sf1,2));
           float trigger_sf1_down = trigger_sf1-sqrt(pow(trigger_sf1_leg1down-trigger_sf1,2)+pow(trigger_sf1_leg2down-trigger_sf1,2));
           float trigger_sf2_up   = trigger_sf2+sqrt(pow(trigger_sf2_leg1up  -trigger_sf2,2)+pow(trigger_sf2_leg2up  -trigger_sf2,2)+pow(trigger_sf2_dzup  -trigger_sf2,2));
           float trigger_sf2_down = trigger_sf2-sqrt(pow(trigger_sf2_leg1down-trigger_sf2,2)+pow(trigger_sf2_leg2down-trigger_sf2,2)+pow(trigger_sf2_dzdown-trigger_sf2,2));

           trigger_sf      = ((Pass_TriggerBG ? trigger_sf1:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2:0.)*LumiH)/LumiBH;
           trigger_sf_up   = ((Pass_TriggerBG ? trigger_sf1_up:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_up:0.)*LumiH)/LumiBH;
           trigger_sf_down = ((Pass_TriggerBG ? trigger_sf1_down:0.)*LumiBG+(Pass_TriggerH ? trigger_sf2_down:0.)*LumiH)/LumiBH;
         }

         btag_sf_LTagup = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpLTag");
         btag_sf_BCTagup= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystUpBCTag");
         btag_sf_JESup  = BTagScaleFactor_1a(jetJESUpColl, snu::KJet::CSVv2, snu::KJet::Medium);
         btag_sf_JERup  = BTagScaleFactor_1a(jetJERUpColl, snu::KJet::CSVv2, snu::KJet::Medium);

         btag_sf_LTagdown = BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownLTag");
         btag_sf_BCTagdown= BTagScaleFactor_1a(jetColl, snu::KJet::CSVv2, snu::KJet::Medium, -1, "SystDownBCTag");
         btag_sf_JESdown  = BTagScaleFactor_1a(jetJESDownColl, snu::KJet::CSVv2, snu::KJet::Medium);
         btag_sf_JERdown  = BTagScaleFactor_1a(jetJERDownColl, snu::KJet::CSVv2, snu::KJet::Medium);

         if(TriMu && k_sample_name.Contains("ZGto2LG")) ZG_SF=0.86191;
         if(EMuMu && k_sample_name.Contains("ZGto2LG")) ZG_SF=0.963169;
       }
       else if(k_running_nonprompt){
         fake_weight = GetFakeWeight_Data(muonLooseColl, electronLooseColl, "POGTIsop6IPp2p1sig4" , "POGTIsop20IPp01p05sig4Chi4", "HctoWAFakeLoose", "POGWP90Isop06IPp025p1sig4", "TrkIsoVVLConeSUSY");
       }
       weight*=ZG_SF;
       
       float systweight_central=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
                                                                                                                                                          
       float systweight_Trigup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf_up*fake_weight;
       float systweight_ElEnup =weight*id_weight_ele_ElEnup*reco_weight_ele_ElEnup*id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight;
       float systweight_MuEnup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu_MuEnup*trk_weight_mu_MuEnup*btag_sf        *trigger_sf   *fake_weight;
       float systweight_JESup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JESup  *trigger_sf   *fake_weight;
       float systweight_JERup  =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_JERup  *trigger_sf   *fake_weight;
       float systweight_LTagup =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_LTagup *trigger_sf   *fake_weight;
       float systweight_BCTagup=weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf_BCTagup*trigger_sf   *fake_weight;
       float systweight_PUup   =weight*id_weight_ele       *reco_weight_ele       *id_weight_mu       *trk_weight_mu       *btag_sf        *trigger_sf   *fake_weight*pileup_reweight_systup;

       float systweight_Trigdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf_down*fake_weight;
       float systweight_ElEndown =weight*id_weight_ele_ElEndown*reco_weight_ele_ElEndown*id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight;
       float systweight_MuEndown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu_MuEndown*trk_weight_mu_MuEndown*btag_sf          *trigger_sf     *fake_weight;
       float systweight_JESdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JESdown  *trigger_sf     *fake_weight;
       float systweight_JERdown  =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_JERdown  *trigger_sf     *fake_weight;
       float systweight_LTagdown =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_LTagdown *trigger_sf     *fake_weight;
       float systweight_BCTagdown=weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf_BCTagdown*trigger_sf     *fake_weight;
       float systweight_PUdown   =weight*id_weight_ele         *reco_weight_ele         *id_weight_mu         *trk_weight_mu         *btag_sf          *trigger_sf     *fake_weight*pileup_reweight_systdown;


       TString ChannelString = EMuMu? "EMuMu":"TriMu";
       //----------------------------------------------------------------------------------------------------------------------//


       DoSystRun("TrilepCR", ChannelString+"",
                 electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                 systweight_central);
       if(!isData){
         DoSystRun("TrilepCR", ChannelString+"SystUpTrig",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_Trigup);
         DoSystRun("TrilepCR", ChannelString+"SystUpPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUup);
         DoSystRun("TrilepCR", ChannelString+"SystUpUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Unclup, met_Unclup*cos(metphi), met_Unclup*sin(metphi),
                   systweight_central);
         DoSystRun("TrilepCR", ChannelString+"SystUpJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESUpColl, bjetJESUpColl, met_JESup, met_JESup*cos(metphi), met_JESup*sin(metphi),
                   systweight_JESup);
         DoSystRun("TrilepCR", ChannelString+"SystUpJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERUpColl, bjetJERUpColl, met_JERup, met_JERup*cos(metphi), met_JERup*sin(metphi),
                   systweight_JERup);
         DoSystRun("TrilepCR", ChannelString+"SystUpBTag_L",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_LTagup);
         DoSystRun("TrilepCR", ChannelString+"SystUpBTag_BC",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_BCTagup);
         DoSystRun("TrilepCR", ChannelString+"SystUpElEn",
                   EleTElEnUpColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEnup, met_ElEnup*cos(metphi), met_ElEnup*sin(metphi),
                   systweight_ElEnup);
         DoSystRun("TrilepCR", ChannelString+"SystUpMuEn",
                   electronColl, electronLooseColl, MuTMuEnUpColl, muonLooseColl, jetColl, bjetColl, met_MuEnup, met_MuEnup*cos(metphi), met_MuEnup*sin(metphi),
                   systweight_MuEnup);


         DoSystRun("TrilepCR", ChannelString+"SystDownTrig",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_Trigdown);
         DoSystRun("TrilepCR", ChannelString+"SystDownPU",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_PUdown);
         DoSystRun("TrilepCR", ChannelString+"SystDownUncl",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_Uncldown, met_Uncldown*cos(metphi), met_Uncldown*sin(metphi),
                   systweight_central);
         DoSystRun("TrilepCR", ChannelString+"SystDownJES",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJESDownColl, bjetJESDownColl, met_JESdown, met_JESdown*cos(metphi), met_JESdown*sin(metphi),
                   systweight_JESdown);
         DoSystRun("TrilepCR", ChannelString+"SystDownJER",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetJERDownColl, bjetJERDownColl, met_JERdown, met_JERdown*cos(metphi), met_JERdown*sin(metphi),
                   systweight_JERdown);
         DoSystRun("TrilepCR", ChannelString+"SystDownBTag_L",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_LTagdown);
         DoSystRun("TrilepCR", ChannelString+"SystDownBTag_BC",
                   electronColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi),
                   systweight_BCTagdown);
         DoSystRun("TrilepCR", ChannelString+"SystDownElEn",
                   EleTElEnDownColl, electronLooseColl, muonColl, muonLooseColl, jetColl, bjetColl, met_ElEndown, met_ElEndown*cos(metphi), met_ElEndown*sin(metphi),
                   systweight_ElEndown);
         DoSystRun("TrilepCR", ChannelString+"SystDownMuEn",
                   electronColl, electronLooseColl, MuTMuEnDownColl, muonLooseColl, jetColl, bjetColl, met_MuEndown, met_MuEndown*cos(metphi), met_MuEndown*sin(metphi),
                   systweight_MuEndown);

       }
       else if( isData && k_running_nonprompt ){
//         DoSystRun("Closure", "", electronLooseColl, electronLooseColl, muonLooseColl, muonLooseColl, jetColl, bjetColl, met, met*cos(metphi), met*sin(metphi), systweight);
       }
     }//End of SystRun
  }//End of Closure


return;
}// End of execute event loop
  


void Mar2018_TrilepCR::CheckSSDilepCRs(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  std::vector<snu::KTruth> truthColl; eventbase->GetTruthSel()->Selection(truthColl);
  if(!isData && (k_sample_name.Contains("TT_powheg") || k_sample_name.Contains("DYJets") || k_sample_name.Contains("WJets") ||
                 k_sample_name.Contains("TTLL") || k_sample_name.Contains("TTLJ"))){
    int NLepCand=0;
    int NPr=0;
    for(int i=0; i<(int) MuTColl.size(); i++){
      int LepType=GetLeptonType(MuTColl.at(i), truthColl);
      if(LepType>0 || LepType<-4) NLepCand++;//For Conv+Flip Est
      //if(LepType<0 && LepType>-5) NLepCand++;//For MCFake Check
      //if(LepType>0 && LepType<4 ) NPr++; //For MCFake Check
    }
    for(int i=0; i<(int) EleTColl.size(); i++){
      int LepType=GetLeptonType(EleTColl.at(i), truthColl);
      if(LepType>0 || LepType<-4) NLepCand++;//For Conv+Flip Est
      //if(LepType<0 && LepType>-5) NLepCand++;//For MCFake Check
      //if(LepType>0 && LepType<4 ) NPr++;//For MCFake Check
    }
    if(NLepCand!=2) return;
    //if(NLepCand==0) return;//For MC FakeCheck
    //if(k_sample_name.Contains("TT_powheg") && NPr!=0) return;//For MC Fake Check
  }

  bool EMu=Option.Contains("EMu")||Option.Contains("EMuMu"), MuMu=!EMu;//By default trimu.
  std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
  std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");

  const int NSignal = 7;
  float Mmumu_CentSeed[NSignal]  = {15., 25., 35., 45., 55., 65., 75.};
  //float Mmumu_WidthSeed[NSignal] = {1., 1., 1.3, 1.5, 1.5, 2.0, 2.5};
  //float Mmumu_Step[NSignal-1]    = {1., 1., 1.2, 1.2, 1.2, 1.8};
  float Mmumu_WidthSeed[NSignal] = {0.5,  0.7, 0.8,  1., 1.2, 1.5, 1.8};
  float Mmumu_Step[NSignal-1]    = {  0.45, 0.55, 0.6, 0.75, 0.9, 1.15};
  vector<float> Mmumu_CentFull, Mmumu_WidthFull;
  for(int it_seed=0; it_seed<NSignal; it_seed++){
    float CurrentCenter = Mmumu_CentSeed[it_seed];
    float CurrentWidth  = Mmumu_WidthSeed[it_seed];
    Mmumu_CentFull.push_back(CurrentCenter);
    Mmumu_WidthFull.push_back(CurrentWidth);
    for(int it_step=0; it_step<100; it_step++){
      if(it_seed==NSignal-1) break;
      CurrentCenter+= Mmumu_Step[it_seed];
      CurrentWidth  = (CurrentCenter-Mmumu_CentSeed[it_seed])/10.*Mmumu_WidthSeed[it_seed+1]
                     +(Mmumu_CentSeed[it_seed+1]-CurrentCenter)/10.*Mmumu_WidthSeed[it_seed];
      if(CurrentCenter<Mmumu_CentSeed[it_seed]+10.){
        Mmumu_CentFull.push_back(CurrentCenter);
        Mmumu_WidthFull.push_back(CurrentWidth);
      }
      else break;
    }
  }

  if(MuMu){

  if( !(MuLColl.size()==2 && EleLColl.size()==0) ) return;
  if( !(MuTColl.size()==2) ) return;
  if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10) ) return;
  if(  BJetColl.size()==0  ) return;
  if( !(JetColl.size()>=3) ) return;
  if( !(MuTColl.at(0).Charge()==MuTColl.at(1).Charge()) ) return;
  float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();

  if( fabs(Mmumu-91.2)<10 || Mmumu<12) return;

  FillHist("PTMu1" +Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
  FillHist("PTMu2" +Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
  FillHist("EtaMu1"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("EtaMu2"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
  FillHist("Mmumu" +Label, Mmumu, weight, 0., 200., 40);
  FillHist("Nj"    +Label, JetColl.size(), weight, 0., 10., 10);
  FillHist("Nb"    +Label, BJetColl.size(), weight, 0., 10., 10);
  FillHist("MET"   +Label, MET, weight, 0., 200., 40);

  if(!isData && (k_sample_name.Contains("TT_powheg") || k_sample_name.Contains("DYJets") || k_sample_name.Contains("WJets"))){
    for(int i=0; i<(int) MuTColl.size(); i++){
      int LepType=GetLeptonType(MuTColl.at(i), truthColl);
      if(!Label.Contains("syst")) FillHist("MuType"+Label, LepType, weight, -10., 10., 20);
    }
    for(int i=0; i<(int) EleTColl.size(); i++){
      int LepType=GetLeptonType(EleTColl.at(i), truthColl);
      if(!Label.Contains("syst")) FillHist("EleType"+Label, LepType, weight, -10., 10., 20);
    }
  }

  //For Mmumu underfluctuation checking
  FillHist("Mmumu_BinW2" +Label, Mmumu, weight, 12., 80., 34);
  for(int it_cent=0; it_cent<34; it_cent++){
    float BinCenter=13.+2.*((float) it_cent);
    float LargeWidth=5.;
    if     ((BinCenter-LargeWidth)<12.)  LargeWidth=BinCenter-12.;
    else if((BinCenter+LargeWidth)>81.2) LargeWidth=81.2-BinCenter;
    if(fabs(Mmumu-BinCenter)<LargeWidth){
      FillHist("Mmumu_BinW2_Width10"+Label, BinCenter, weight, 12., 80., 34);
    }
    if(Label==""){
      float Variation=0.2;
      float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
      if     ((BinCenter-LargeWidthUp)<12. ) LargeWidthUp=BinCenter-12.;
      else if((BinCenter+LargeWidthUp)>81.2) LargeWidthUp=81.2-BinCenter;
      if(fabs(Mmumu-BinCenter)<LargeWidthUp){
        FillHist("Mmumu_BinW2_Width10_systup_Mrange", BinCenter, weight, 12., 80., 34);
      }
      if((BinCenter-LargeWidthDown)<12.) LargeWidthDown=BinCenter-12.;
      if(fabs(Mmumu-BinCenter)<LargeWidthDown){
        FillHist("Mmumu_BinW2_Width10_systdown_Mrange", BinCenter, weight, 12., 80., 34);
      }
    }
  }
  for(int it_cent=0; it_cent<Mmumu_CentFull.size(); it_cent++){
    if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<Mmumu_WidthFull.at(it_cent)){
      FillHist("Mmumu_CntBin"+Label, it_cent+0.01, weight, 0., 95., 95);
    }
    float LargeWidth=5.;
    if((Mmumu_CentFull.at(it_cent)-LargeWidth)<12.) LargeWidth=Mmumu_CentFull.at(it_cent)-12.;
    if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidth){
      FillHist("Mmumu_CntBin_Width10"+Label, it_cent+0.01, weight, 0., 95., 95);
    }
    if(Label==""){
      float Variation=0.2;
      float LargeWidthUp=(1.+Variation)*LargeWidth, LargeWidthDown=(1.-Variation)*LargeWidth;
      if((Mmumu_CentFull.at(it_cent)-LargeWidthUp)<12.) LargeWidthUp=Mmumu_CentFull.at(it_cent)-12.;
      if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthUp){
        FillHist("Mmumu_CntBin_Width10_systup_Mrange", it_cent+0.01, weight, 0., 95., 95);
      }
      if((Mmumu_CentFull.at(it_cent)-LargeWidthDown)<12.) LargeWidthDown=Mmumu_CentFull.at(it_cent)-12.;
      if(fabs(Mmumu-Mmumu_CentFull.at(it_cent))<LargeWidthDown){
        FillHist("Mmumu_CntBin_Width10_systdown_Mrange", it_cent+0.01, weight, 0., 95., 95);
      }
    }
  }


  }//End of SS MuMu
  else if(EMu){

  if( !(MuLColl.size()==1 && EleLColl.size()==1) ) return;
  if( !(MuTColl.size()==1 && EleTColl.size()==1) ) return;
  if( !(MuTColl.at(0).Pt()>10. && EleTColl.at(0).Pt()>25.) ) return;
  if(  BJetColl.size()==0  ) return;
  if( !(JetColl.size()>=3) ) return;
  
  if( !(MuTColl.at(0).Charge()==EleTColl.at(0).Charge()) ) return;

  FillHist("PTMu1" +Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
  FillHist("PTEl1" +Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
  FillHist("EtaMu1"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("EtaEl1"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
  FillHist("Nj"    +Label, JetColl.size(), weight, 0., 10., 10);
  FillHist("Nb"    +Label, BJetColl.size(), weight, 0., 10., 10);
  FillHist("MET"   +Label, MET, weight, 0., 200., 40);

  if(!isData && (k_sample_name.Contains("TT_powheg") || k_sample_name.Contains("DYJets") || k_sample_name.Contains("WJets"))){
    for(int i=0; i<(int) MuTColl.size(); i++){
      int LepType=GetLeptonType(MuTColl.at(i), truthColl);
      if(!Label.Contains("syst")) FillHist("MuType"+Label, LepType, weight, -10., 10., 20);
    }
    for(int i=0; i<(int) EleTColl.size(); i++){
      int LepType=GetLeptonType(EleTColl.at(i), truthColl);
      if(!Label.Contains("syst")) FillHist("EleType"+Label, LepType, weight, -10., 10., 20);
    }
  }

  }//End of SS EMu
}


void Mar2018_TrilepCR::CheckTTZAnomaly(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  bool EMuMu=Option.Contains("EMuMu"), TriMu=!EMuMu;//By default trimu.
   std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
   std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");

  if(TriMu){
    FillHist("CutFlow", 0., weight, 0., 10., 10);
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>25 && MuTColl.at(1).Pt()>20 && MuTColl.at(2).Pt()>20) ) return;
  
    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
  
    if(BJetColl.size()==0) return;
    if(JetColl.size()<2) return;


    int   IdxSSZ    = fabs(MOSSS1-91.2)<fabs(MOSSS2-91.2)? IdxSS1:IdxSS2;
    int   IdxSSNonZ = 3-IdxOS-IdxSSZ;
    float MOSSSZ = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSZ)).M();
    float M3l    = (MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)).M();
    float MTW    = sqrt(2)*sqrt(MET*MuTColl.at(IdxSSNonZ).Pt()-METx*MuTColl.at(IdxSSNonZ).Px()-METy*MuTColl.at(IdxSSNonZ).Py());

    int Idxj1_W=-1, Idxj2_W=-1;
    snu::KParticle v; v.SetPxPyPzE(METx, METy, 0, sqrt(METx*METx+METy*METy));
    float Pzv1 = GetvPz(v, MuTColl.at(0), 1), Pzv2=GetvPz(v, MuTColl.at(0), 2);
    float Pzv  = fabs(Pzv1)<fabs(Pzv2)? Pzv1:Pzv2;
          v.SetXYZM(METx,METy,Pzv,0.);

    float M3lv=(MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)+v).M(), M2l2j=0.; 
    Idxj1_W=GetDaughterCandIdx(JetColl, "W", -1., "Lead");
    Idxj2_W=GetDaughterCandIdx(JetColl, "W", -1., "Subl");
    M2l2j=(MuTColl.at(IdxSSZ)+MuTColl.at(IdxOS)+JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W)).M();
    float Mmumu=MOSSSZ;

  
    bool HasBJet=false, OnZ=false, OffZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false, ZbbSel=false;
    bool ttZSelTight=false;
    if( BJetColl.size()!=0                 ) HasBJet     = true;
    if( fabs(MOSSSZ-91.2)<10               ) OnZ         = true;
    if( fabs(M3l-91.2)<10                  ) OnZG        = true;
    if( OnZ && HasBJet && JetColl.size()>1 ) ttZSel      = true;
    if(ttZSel && MuTColl.at(IdxSSNonZ).Pt()>25.){
      FillHist("PTmu1_ttZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_ttZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_ttZSel"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_ttZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_ttZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_ttZSel"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);

      FillHist("NChannel_ttZSel"+Label, 1., weight, 0., 2., 2);
      FillHist("PTlW_ttZSel"+Label, MuTColl.at(IdxSSNonZ).Pt(), weight, 0., 200., 40);
      int IdxlZ1=min(IdxSSZ,IdxOS), IdxlZ2=max(IdxSSZ,IdxOS);
      FillHist("PTlZ1_ttZSel"+Label, MuTColl.at(IdxlZ1).Pt(), weight, 0., 300., 60);
      FillHist("PTlZ2_ttZSel"+Label, MuTColl.at(IdxlZ2).Pt(), weight, 0., 200., 40);
      FillHist("EtalW_ttZSel"+Label, MuTColl.at(IdxSSNonZ).Eta(), weight, -5., 5., 20);
      FillHist("EtalZ1_ttZSel"+Label, MuTColl.at(IdxlZ1).Eta(), weight, -5., 5., 20);
      FillHist("EtalZ2_ttZSel"+Label, MuTColl.at(IdxlZ2).Eta(), weight, -5., 5., 20);

      FillHist("MOSSSZ_ttZSel"+Label, MOSSSZ, weight, 60., 120., 30);
      FillHist("Mmumu_ttZSel"+Label, MOSSSZ, weight, 60., 120., 30);
      FillHist("M3l_ttZSel"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_ttZSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ttZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_ttZSel"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_ttZSel"+Label, MTW, weight, 0., 200., 20);
      FillHist("M2l2j_ttZSel"+Label, M2l2j, weight, 0., 1000., 100);
      FillHist("M3lv_ttZSel"+Label, M3lv, weight, 0., 1000., 100);

      if(fabs(Mmumu-91.2)<5){
        FillHist("M2l2j_ttZSel_TWin"+Label, M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_ttZSel_TWin"+Label, M3lv, weight, 0., 1000., 100);
      }
      if(BJetColl.size()>1){
        FillHist("M2l2j_ttZSel_TB"+Label, M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_ttZSel_TB"+Label, M3lv, weight, 0., 1000., 100);
      }
    }
  }
  else if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>20) ) return;
  
    float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
    float M3l=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();

    if(Mmumu<12) return;

    if(BJetColl.size()==0) return;
    if(JetColl.size()<2) return;
    

    int Idxj1_W=-1, Idxj2_W=-1;
    snu::KParticle v; v.SetPxPyPzE(METx, METy, 0, sqrt(METx*METx+METy*METy));
    float Pzv1=GetvPz(v, EleTColl.at(0), 1), Pzv2=GetvPz(v, EleTColl.at(0), 2);
    float Pzv=fabs(Pzv1)<fabs(Pzv2)? Pzv1:Pzv2;
    v.SetXYZM(METx,METy,Pzv,0.);

    float M3lv=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)+v).M(), M2l2j=0.; 
    Idxj1_W=GetDaughterCandIdx(JetColl, "W", -1., "Lead");
    Idxj2_W=GetDaughterCandIdx(JetColl, "W", -1., "Subl");
    M2l2j=(MuTColl.at(0)+MuTColl.at(1)+JetColl.at(Idxj1_W)+JetColl.at(Idxj2_W)).M();

    
    bool HasBJet=false, OnZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false, ZbbSel=false;
    if( BJetColl.size()!=0 )                      HasBJet = true;
    if( fabs(Mmumu-91.2)<10     )                 OnZ     = true;
    if( OnZ && HasBJet && JetColl.size()>1)       ttZSel  = true;
    if(ttZSel){
      FillHist("PTe_ttZSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_ttZSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_ttZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_ttZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_ttZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_ttZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("dRemu1_ttZSel"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRemu2_ttZSel"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRmu1mu2_ttZSel"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);

      FillHist("NChannel_ttZSel"+Label, 0., weight, 0., 2., 2);
      FillHist("PTlW_ttZSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("EtalW_ttZSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTlZ1_ttZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTlZ2_ttZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("EtalZ1_ttZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("EtalZ2_ttZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);

      FillHist("Mmumu_ttZSel"+Label, Mmumu, weight, 60., 120., 30);
      FillHist("M3l_ttZSel"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
      FillHist("Nj_ttZSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ttZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_ttZSel"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_ttZSel"+Label, MTW, weight, 0., 200., 20);

      FillHist("M2l2j_ttZSel"+Label, M2l2j, weight, 0., 1000., 100);
      FillHist("M3lv_ttZSel"+Label, M3lv, weight, 0., 1000., 100);

      if(fabs(Mmumu-91.2)<5){
        FillHist("M2l2j_ttZSel_TWin"+Label, M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_ttZSel_TWin"+Label, M3lv, weight, 0., 1000., 100);
      }
      if(BJetColl.size()>1){
        FillHist("M2l2j_ttZSel_TB"+Label, M2l2j, weight, 0., 1000., 100);
        FillHist("M3lv_ttZSel_TB"+Label, M3lv, weight, 0., 1000., 100);
      }
    }
  }

}

void Mar2018_TrilepCR::CheckTrilepCRs(std::vector<snu::KMuon> MuTColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KElectron> EleTColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KJet> JetNoVetoColl, std::vector<snu::KJet> BJetNoVetoColl, float MET, float METx, float METy, float weight, TString Label, TString Option){

  bool EMuMu=Option.Contains("EMuMu"), TriMu=!EMuMu;//By default trimu.
   std::vector<snu::KJet> JetColl  = SkimJetColl(JetNoVetoColl,  EleLColl, MuLColl, "EleMuVeto");
   std::vector<snu::KJet> BJetColl = SelBJets(JetColl, "Medium");

  if(TriMu){
    FillHist("CutFlow", 0., weight, 0., 10., 10);
    if( !(MuLColl.size()==3 && EleLColl.size()==0) ) return;
    if( !(MuTColl.size()==3) ) return;
    if( fabs(SumCharge(MuTColl))!=1 ) return;
    if( !(MuTColl.at(0).Pt()>20 && MuTColl.at(1).Pt()>10 && MuTColl.at(2).Pt()>10) ) return;
  
    int   IdxOS  = TriMuChargeIndex(MuTColl, "OS");
    int   IdxSS1 = TriMuChargeIndex(MuTColl, "SS1");
    int   IdxSS2 = TriMuChargeIndex(MuTColl, "SS2");
    float MOSSS1 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS1)).M();
    float MOSSS2 = (MuTColl.at(IdxOS)+MuTColl.at(IdxSS2)).M();
    if( !(MOSSS1>12 && MOSSS2>12) ) return;
    FillHist("CutFlow", 1., weight, 0., 10., 10);
  
    int   IdxSSZ    = fabs(MOSSS1-91.2)<fabs(MOSSS2-91.2)? IdxSS1:IdxSS2;
    int   IdxSSNonZ = 3-IdxOS-IdxSSZ;
    float MOSSSZ = (MuTColl.at(IdxOS)+MuTColl.at(IdxSSZ)).M();
    float M3l    = (MuTColl.at(0)+MuTColl.at(1)+MuTColl.at(2)).M();
    float MTW    = sqrt(2)*sqrt(MET*MuTColl.at(IdxSSNonZ).Pt()-METx*MuTColl.at(IdxSSNonZ).Px()-METy*MuTColl.at(IdxSSNonZ).Py());
  
    bool HasBJet=false, OnZ=false, OffZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false, ZbbSel=false;
    bool ttZSelTight=false;
    if( BJetColl.size()!=0                 ) HasBJet     = true;
    if( fabs(MOSSSZ-91.2)<10               ) OnZ         = true;
    if( fabs(MOSSS1-91.2)>10 && fabs(MOSSS2-91.2)>10)OffZ= true;
    if( fabs(M3l-91.2)<10                  ) OnZG        = true;
    if( OnZ && M3l>101.2 && MET>50         ) WZSel       = true;
    if( OnZG && OffZ && (MOSSS1<81.2 || MOSSS2<81.2) && MET<50 ) ZGSel = true;
    if( OnZ && HasBJet && JetColl.size()>2 ) ttZSel      = true;
    if( OnZ && BJetColl.size()>1           ) ttZSelTight = true;
  
    //General 3l Selection
    if(!HasBJet){
      FillHist("CutFlow", 1., weight, 0., 10., 10);

      FillHist("PTmu1_3lOS"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_3lOS"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_3lOS"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_3lOS"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_3lOS"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_3lOS"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSS1_3lOS"+Label, MOSSS1, weight, 0., 200., 40);
      FillHist("MOSSS2_3lOS"+Label, MOSSS2, weight, 0., 200., 40);
      FillHist("M3l_3lOS"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_3lOS"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOS"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOS"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_3lOS"+Label, MTW, weight, 0., 200., 20);
    }
    //3l OnZ ; Fake CR
    if(OnZ){
      FillHist("CutFlow", 2., weight, 0., 10., 10);

      FillHist("PTmu1_3lOSOnZ"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_3lOSOnZ"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_3lOSOnZ"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_3lOSOnZ"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_3lOSOnZ"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_3lOSOnZ"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSSZ_3lOSOnZ"+Label, MOSSSZ, weight, 60., 120., 30);
      FillHist("M3l_3lOSOnZ"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_3lOSOnZ"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOSOnZ"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOSOnZ"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_3lOSOnZ"+Label, MTW, weight, 0., 200., 20);
    }
    //3l OnZ & HasBJet
    if(HasBJet && OnZ){
      FillHist("CutFlow", 3., weight, 0., 10., 10);

      FillHist("PTmu1_3lOSOnZHasB"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_3lOSOnZHasB"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_3lOSOnZHasB"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_3lOSOnZHasB"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_3lOSOnZHasB"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_3lOSOnZHasB"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSSZ_3lOSOnZHasB"+Label, MOSSSZ, weight, 60., 120., 30);
      FillHist("M3l_3lOSOnZHasB"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_3lOSOnZHasB"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOSOnZHasB"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOSOnZHasB"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_3lOSOnZHasB"+Label, MTW, weight, 0., 200., 20);
    }
    if(WZSel){
      FillHist("CutFlow", 4., weight, 0., 10., 10);

      FillHist("PTmu1_WZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_WZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_WZSel"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_WZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_WZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_WZSel"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSSZ_WZSel"+Label, MOSSSZ, weight, 60., 120., 30);
      FillHist("M3l_WZSel"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_WZSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_WZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_WZSel"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_WZSel"+Label, MTW, weight, 0., 200., 20);
    }
    if(ZGSel){
      FillHist("CutFlow", 5., weight, 0., 10., 10);

      FillHist("PTmu1_ZGSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_ZGSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_ZGSel"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_ZGSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_ZGSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_ZGSel"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSSZ_ZGSel"+Label, MOSSSZ, weight, 0., 200., 40);
      FillHist("M3l_ZGSel"+Label, M3l, weight, 60., 120., 30);
      FillHist("Nj_ZGSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ZGSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_ZGSel"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_ZGSel"+Label, MTW, weight, 0., 200., 20);
      FillHist("Count_ZGSel"+Label, 0., weight, 0., 1., 1);
    }
    if(ttZSel){
      FillHist("CutFlow", 6., weight, 0., 10., 10);

      FillHist("PTmu1_ttZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_ttZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_ttZSel"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_ttZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_ttZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_ttZSel"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);

      FillHist("NChannel_ttZSel"+Label, 1., weight, 0., 2., 2);
      FillHist("PTlW_ttZSel"+Label, MuTColl.at(IdxSSNonZ).Pt(), weight, 0., 200., 40);
      int IdxlZ1=min(IdxSSZ,IdxOS), IdxlZ2=max(IdxSSZ,IdxOS);
      FillHist("PTlZ1_ttZSel"+Label, MuTColl.at(IdxlZ1).Pt(), weight, 0., 300., 60);
      FillHist("PTlZ2_ttZSel"+Label, MuTColl.at(IdxlZ2).Pt(), weight, 0., 200., 40);
      FillHist("EtalW_ttZSel"+Label, MuTColl.at(IdxSSNonZ).Eta(), weight, -5., 5., 20);
      FillHist("EtalZ1_ttZSel"+Label, MuTColl.at(IdxlZ1).Eta(), weight, -5., 5., 20);
      FillHist("EtalZ2_ttZSel"+Label, MuTColl.at(IdxlZ2).Eta(), weight, -5., 5., 20);

      FillHist("MOSSSZ_ttZSel"+Label, MOSSSZ, weight, 60., 120., 30);
      FillHist("Mmumu_ttZSel"+Label, MOSSSZ, weight, 60., 120., 30);
      FillHist("M3l_ttZSel"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_ttZSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ttZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_ttZSel"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_ttZSel"+Label, MTW, weight, 0., 200., 20);
    }
    if(OffZ && BJetColl.size()==1 && JetColl.size()==1){
      FillHist("PTmu1_TopFakeCR"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_TopFakeCR"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_TopFakeCR"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_TopFakeCR"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_TopFakeCR"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_TopFakeCR"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSS1_TopFakeCR"+Label, MOSSS1, weight, 0., 300., 60);
      FillHist("MOSSS2_TopFakeCR"+Label, MOSSS2, weight, 0., 200., 40);
      FillHist("M3l_TopFakeCR"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_TopFakeCR"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_TopFakeCR"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_TopFakeCR"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_TopFakeCR"+Label, MTW, weight, 0., 200., 20);
    }
    if( ( (!isData) || k_running_nonprompt )
        && OffZ && HasBJet && JetColl.size()>1 ){
      
      FillHist("PTmu1_SR2j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_SR2j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_SR2j"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_SR2j"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_SR2j"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_SR2j"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSS1_SR2j"+Label, MOSSS1, weight, 0., 300., 60);
      FillHist("MOSSS2_SR2j"+Label, MOSSS2, weight, 0., 200., 40);
      FillHist("M3l_SR2j"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_SR2j"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_SR2j"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_SR2j"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_SR2j"+Label, MTW, weight, 0., 200., 20);
    }
    if( ( (!isData) || k_running_nonprompt )
        && OffZ && HasBJet && JetColl.size()>2 ){
      
      FillHist("PTmu1_SR3j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_SR3j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("PTmu3_SR3j"+Label, MuTColl.at(2).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_SR3j"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_SR3j"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("Etamu3_SR3j"+Label, MuTColl.at(2).Eta(), weight, -5., 5., 20);
    
      FillHist("MOSSS1_SR3j"+Label, MOSSS1, weight, 0., 300., 60);
      FillHist("MOSSS2_SR3j"+Label, MOSSS2, weight, 0., 200., 40);
      FillHist("M3l_SR3j"+Label, M3l, weight, 0., 500., 50);
      FillHist("Nj_SR3j"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_SR3j"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_SR3j"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_SR3j"+Label, MTW, weight, 0., 200., 20);
    }
  }
  else if(EMuMu){
    if( !(EleLColl.size()==1 && MuLColl.size()==2) ) return;
    if( !(EleTColl.size()==1 && MuTColl.size()==2) ) return;
    if( !(MuTColl.at(0).Charge()!=MuTColl.at(1).Charge()) ) return;
    if( !(EleTColl.at(0).Pt()>25 && MuTColl.at(0).Pt()>10 && MuTColl.at(1).Pt()>10) ) return;

    float Mmumu=(MuTColl.at(0)+MuTColl.at(1)).M();
    float MTW = sqrt(2)*sqrt(MET*EleTColl.at(0).Pt()-METx*EleTColl.at(0).Px()-METy*EleTColl.at(0).Py());
    float M3l=(EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M();
    if(Mmumu<12) return;
    
    
    bool HasBJet=false, OnZ=false, OnZG=false, WZSel=false, ZGSel=false, ttZSel=false, ZbbSel=false;
    if( BJetColl.size()!=0 )                      HasBJet = true;
    if( fabs(Mmumu-91.2)<10     )                 OnZ     = true;
    if( fabs(M3l-91.2)<10       )                 OnZG    = true;
    if( OnZ && M3l>101.2 && MET>50 )              WZSel   = true;
    if( OnZG && Mmumu<81.2 && MET<50 )            ZGSel   = true;
    if( OnZ && HasBJet && JetColl.size()>2)       ttZSel  = true;
    if( OnZ && HasBJet && JetColl.size()<3)       ZbbSel  = true;
     

    //General 3l Selection
    if(!HasBJet){
      FillHist("PTe_3lOS"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_3lOS"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_3lOS"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_3lOS"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_3lOS"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_3lOS"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("dRemu1_3lOS"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRemu2_3lOS"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRmu1mu2_3lOS"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
    
      FillHist("Mmumu_3lOS"+Label, Mmumu, weight, 0., 200., 40);
      FillHist("M3l_3lOS"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
      FillHist("Nj_3lOS"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOS"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOS"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_3lOS"+Label, MTW, weight, 0., 200., 20);
    }
    //3l OnZ ; Fake CR
    if(OnZ){
      FillHist("PTe_3lOSOnZ"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_3lOSOnZ"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_3lOSOnZ"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_3lOSOnZ"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_3lOSOnZ"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_3lOSOnZ"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    
      FillHist("Mmumu_3lOSOnZ"+Label, Mmumu, weight, 60., 120., 30);
      FillHist("M3l_3lOSOnZ"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
      FillHist("Nj_3lOSOnZ"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOSOnZ"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOSOnZ"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_3lOSOnZ"+Label, MTW, weight, 0., 200., 20);
    }
    //3l OnZ & HasBJet
    if(HasBJet && OnZ){
      FillHist("PTe_3lOSOnZHasB"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_3lOSOnZHasB"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_3lOSOnZHasB"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_3lOSOnZHasB"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_3lOSOnZHasB"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_3lOSOnZHasB"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    
      FillHist("Mmumu_3lOSOnZHasB"+Label, Mmumu, weight, 60., 120., 30);
      FillHist("M3l_3lOSOnZHasB"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
      FillHist("Nj_3lOSOnZHasB"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_3lOSOnZHasB"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_3lOSOnZHasB"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_3lOSOnZHasB"+Label, MTW, weight, 0., 200., 20);
    }
    if(!OnZ && BJetColl.size()==1 && JetColl.size()==1){
      FillHist("PTe_TopFakeCR"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_TopFakeCR"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_TopFakeCR"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_TopFakeCR"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_TopFakeCR"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_TopFakeCR"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    
      FillHist("Mmumu_TopFakeCR"+Label, Mmumu, weight, 0., 200., 40);
      FillHist("M3l_TopFakeCR"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
      FillHist("Nj_TopFakeCR"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_TopFakeCR"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_TopFakeCR"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_TopFakeCR"+Label, MTW, weight, 0., 200., 20);
    }
    if(WZSel){
      FillHist("PTe_WZSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_WZSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_WZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_WZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_WZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_WZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    
      FillHist("Mmumu_WZSel"+Label, Mmumu, weight, 60., 120., 30);
      FillHist("M3l_WZSel"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
      FillHist("Nj_WZSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_WZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_WZSel"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_WZSel"+Label, MTW, weight, 0., 200., 20);
    }
    if(ZGSel){
      FillHist("PTe_ZGSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_ZGSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_ZGSel"+Label, MuTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("PTmu2_ZGSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_ZGSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_ZGSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
    
      FillHist("Mmumu_ZGSel"+Label, Mmumu, weight, 0., 200., 40);
      FillHist("M3l_ZGSel"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 60., 120., 30);
      FillHist("Nj_ZGSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ZGSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_ZGSel"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_ZGSel"+Label, MTW, weight, 0., 200., 20);
      FillHist("Count_ZGSel"+Label, 0., weight, 0., 1., 1);
    }
    if(ttZSel){
      FillHist("PTe_ttZSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_ttZSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_ttZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_ttZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_ttZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_ttZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("dRemu1_ttZSel"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRemu2_ttZSel"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRmu1mu2_ttZSel"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);


      FillHist("NChannel_ttZSel"+Label, 0., weight, 0., 2., 2);
      FillHist("PTlW_ttZSel"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("EtalW_ttZSel"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTlZ1_ttZSel"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTlZ2_ttZSel"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("EtalZ1_ttZSel"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("EtalZ2_ttZSel"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);

      FillHist("Mmumu_ttZSel"+Label, Mmumu, weight, 60., 120., 30);
      FillHist("M3l_ttZSel"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
      FillHist("Nj_ttZSel"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_ttZSel"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_ttZSel"+Label, MET, weight, 0., 300., 30);
      FillHist("MTW_ttZSel"+Label, MTW, weight, 0., 200., 20);
    }
    if(!OnZ){
      FillHist("PTe_OffZ"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
      FillHist("Etae_OffZ"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("PTmu1_OffZ"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
      FillHist("PTmu2_OffZ"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
      FillHist("Etamu1_OffZ"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
      FillHist("Etamu2_OffZ"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
      FillHist("dRemu1_OffZ"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRemu2_OffZ"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
      FillHist("dRmu1mu2_OffZ"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
    
      FillHist("Mmumu_OffZ"+Label, Mmumu, weight, 0., 200., 40);
      FillHist("M3l_OffZ"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
      FillHist("Nj_OffZ"+Label, JetColl.size(), weight, 0., 10., 10);
      FillHist("Nb_OffZ"+Label, BJetColl.size(), weight, 0., 10., 10);
      FillHist("MET_OffZ"+Label, MET, weight, 0., 200., 20);
      FillHist("MTW_OffZ"+Label, MTW, weight, 0., 200., 20);
    }
    if( !isData || k_running_nonprompt){
      if(BJetColl.size()>0 && JetColl.size()>2){
        FillHist("PTe_SRincl"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
        FillHist("Etae_SRincl"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
        FillHist("PTmu1_SRincl"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
        FillHist("PTmu2_SRincl"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
        FillHist("Etamu1_SRincl"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
        FillHist("Etamu2_SRincl"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
        FillHist("dRemu1_SRincl"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
        FillHist("dRemu2_SRincl"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
        FillHist("dRmu1mu2_SRincl"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
      
        FillHist("Mmumu_SRincl"+Label, Mmumu, weight, 0., 200., 40);
        FillHist("M3l_SRincl"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
        FillHist("Nj_SRincl"+Label, JetColl.size(), weight, 0., 10., 10);
        FillHist("Nb_SRincl"+Label, BJetColl.size(), weight, 0., 10., 10);
        FillHist("MET_SRincl"+Label, MET, weight, 0., 300., 30);
        FillHist("MTW_SRincl"+Label, MTW, weight, 0., 200., 20);
      }
      if( !OnZ && BJetColl.size()>0 && JetColl.size()>1){
        FillHist("PTe_SR2j"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
        FillHist("Etae_SR2j"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
        FillHist("PTmu1_SR2j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
        FillHist("PTmu2_SR2j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
        FillHist("Etamu1_SR2j"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
        FillHist("Etamu2_SR2j"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
        FillHist("dRemu1_SR2j"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
        FillHist("dRemu2_SR2j"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
        FillHist("dRmu1mu2_SR2j"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
      
        FillHist("Mmumu_SR2j"+Label, Mmumu, weight, 0., 200., 40);
        FillHist("M3l_SR2j"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
        FillHist("Nj_SR2j"+Label, JetColl.size(), weight, 0., 10., 10);
        FillHist("Nb_SR2j"+Label, BJetColl.size(), weight, 0., 10., 10);
        FillHist("MET_SR2j"+Label, MET, weight, 0., 300., 30);
        FillHist("MTW_SR2j"+Label, MTW, weight, 0., 200., 20);
  
        FillHist("Count_SR2j"+Label, 0., weight, 0., 1., 1);
      }
      if( !OnZ && BJetColl.size()>0 && JetColl.size()>2){
        FillHist("PTe_SR3j"+Label, EleTColl.at(0).Pt(), weight, 0., 200., 40);
        FillHist("Etae_SR3j"+Label, EleTColl.at(0).Eta(), weight, -5., 5., 20);
        FillHist("PTmu1_SR3j"+Label, MuTColl.at(0).Pt(), weight, 0., 300., 60);
        FillHist("PTmu2_SR3j"+Label, MuTColl.at(1).Pt(), weight, 0., 200., 40);
        FillHist("Etamu1_SR3j"+Label, MuTColl.at(0).Eta(), weight, -5., 5., 20);
        FillHist("Etamu2_SR3j"+Label, MuTColl.at(1).Eta(), weight, -5., 5., 20);
        FillHist("dRemu1_SR3j"+Label, MuTColl.at(0).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
        FillHist("dRemu2_SR3j"+Label, MuTColl.at(1).DeltaR(EleTColl.at(0)), weight, 0., 5., 50);
        FillHist("dRmu1mu2_SR3j"+Label, MuTColl.at(0).DeltaR(MuTColl.at(1)), weight, 0., 5., 50);
      
        FillHist("Mmumu_SR3j"+Label, Mmumu, weight, 0., 200., 40);
        FillHist("M3l_SR3j"+Label, (EleTColl.at(0)+MuTColl.at(0)+MuTColl.at(1)).M(), weight, 0., 500., 50);
        FillHist("Nj_SR3j"+Label, JetColl.size(), weight, 0., 10., 10);
        FillHist("Nb_SR3j"+Label, BJetColl.size(), weight, 0., 10., 10);
        FillHist("MET_SR3j"+Label, MET, weight, 0., 300., 30);
        FillHist("MTW_SR3j"+Label, MTW, weight, 0., 200., 20);
      }
    } 
  }

}


void Mar2018_TrilepCR::DoSystRun(TString Cycle, TString Mode, std::vector<snu::KElectron> EleColl, std::vector<snu::KElectron> EleLColl, std::vector<snu::KMuon> MuColl, std::vector<snu::KMuon> MuLColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, float METx, float METy, float weight){

  int SystDir=0; TString SystKindLabel="", SystDirLabel="", ChannelLabel="";
  if(Mode.Contains("Syst")){
    if(isData&& !k_running_nonprompt) return;
    if     (Mode.Contains("Up"))      {SystDir= 1; SystDirLabel="_systup";}
    else if(Mode.Contains("Down"))    {SystDir=-1; SystDirLabel="_systdown";}
    if     (Mode.Contains("Up"))      {SystDir= 1; SystDirLabel="_systup";}
    else if(Mode.Contains("Down"))    {SystDir=-1; SystDirLabel="_systdown";}
    if     (Mode.Contains("PU"))      SystKindLabel="_PU";
    else if(Mode.Contains("Nvtx"))    SystKindLabel="_Nvtx";
    else if(Mode.Contains("JES"))     SystKindLabel="_JES";
    else if(Mode.Contains("JER"))     SystKindLabel="_JER";
    else if(Mode.Contains("Uncl"))    SystKindLabel="_Uncl";
    else if(Mode.Contains("ElID"))    SystKindLabel="_ElID";
    else if(Mode.Contains("MuID"))    SystKindLabel="_MuID";
    else if(Mode.Contains("ElEn"))    SystKindLabel="_ElEn";
    else if(Mode.Contains("MuEn"))    SystKindLabel="_MuEn";
    else if(Mode.Contains("FR"))      SystKindLabel="_FR";
    else if(Mode.Contains("BTag_L"))  SystKindLabel="_BTag_L";
    else if(Mode.Contains("BTag_BC")) SystKindLabel="_BTag_BC";
    else if(Mode.Contains("Trig"))    SystKindLabel="_Trig";
    else if(Mode.Contains("Xsec"))    SystKindLabel="_Xsec";
    else if(Mode.Contains("Conv"))    SystKindLabel="_Conv";
  }
  if     (Mode.Contains("EMuMu"))   ChannelLabel="EMuMu";
  else if(Mode.Contains("TriMu"))   ChannelLabel="TriMu";


  if(Cycle=="TrilepCR"){
    CheckTrilepCRs(MuColl, MuLColl, EleColl, EleLColl, JetColl, BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, ChannelLabel);
  }//End of Closure Cycle
  else if(Cycle=="SSDilepCR"){
    CheckSSDilepCRs(MuColl, MuLColl, EleColl, EleLColl, JetColl, BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, ChannelLabel);
  }//End of Closure Cycle
  else if(Cycle=="TTZAnomaly"){
    CheckTTZAnomaly(MuColl, MuLColl, EleColl, EleLColl, JetColl, BJetColl, MET, METx, METy, weight, SystDirLabel+SystKindLabel, ChannelLabel);
  }//End of Closure Cycle


  return;
}


void Mar2018_TrilepCR::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void Mar2018_TrilepCR::BeginCycle() throw( LQError ){
  
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

Mar2018_TrilepCR::~Mar2018_TrilepCR() {
  
  Message("In Mar2018_TrilepCR Destructor" , INFO);
//  if(!k_isdata)delete reweightPU;
  
}




float Mar2018_TrilepCR::FakeRateData(snu::KMuon Mu, TString Option){

  float FR=0., TightIsoCut=0.;
  if(Option.Contains("ConeSUSY")){
    if     (Option.Contains("TIsop15")) TightIsoCut = 0.15;
    else if(Option.Contains("TIsop20")) TightIsoCut = 0.2;
    else if(Option.Contains("TIsop25")) TightIsoCut = 0.25;
  }
  float PTCorr=ConeCorrectedPT(Mu, TightIsoCut);
    if(Option.Contains("FOPt")) PTCorr=Mu.Pt();
  float fEta=fabs(Mu.Eta());

  if(Option=="POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.103166 ;
      else if(PTCorr<20)  FR=0.0800024;
      else if(PTCorr<25)  FR=0.0847176;
      else if(PTCorr<35)  FR=0.0680578;
      else                FR=0.0658314;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.095482 ;
      else if(PTCorr<20)  FR=0.0713299;
      else if(PTCorr<25)  FR=0.102061 ;
      else if(PTCorr<35)  FR=0.114623 ;
      else                FR=0.0850804;
    }
    else{
      if     (PTCorr<15)  FR=0.0603521;
      else if(PTCorr<20)  FR=0.0918451;
      else if(PTCorr<25)  FR=0.0671922;
      else if(PTCorr<35)  FR=0.128504 ;
      else                FR=0.0964251;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.282211;
      else if(PTCorr<20)  FR=0.2111  ;
      else if(PTCorr<25)  FR=0.187283;
      else if(PTCorr<35)  FR=0.195845;
      else                FR=0.193727;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.321604;
      else if(PTCorr<20)  FR=0.255115;
      else if(PTCorr<25)  FR=0.253692;
      else if(PTCorr<35)  FR=0.197637;
      else                FR=0.238313;
    }
    else{
      if     (PTCorr<15)  FR=0.35875 ;
      else if(PTCorr<20)  FR=0.303843;
      else if(PTCorr<25)  FR=0.298489;
      else if(PTCorr<35)  FR=0.255767;
      else                FR=0.266638;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGLIsop4IPp5p1Chi100_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.276247;
      else if(PTCorr<20)  FR=0.191163;
      else if(PTCorr<25)  FR=0.158219;
      else if(PTCorr<35)  FR=0.155122;
      else                FR=0.153071;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.309303;
      else if(PTCorr<20)  FR=0.219685;
      else if(PTCorr<25)  FR=0.211428;
      else if(PTCorr<35)  FR=0.170437;
      else                FR=0.19083 ;
    }
    else{
      if     (PTCorr<15)  FR=0.339582;
      else if(PTCorr<20)  FR=0.263034;
      else if(PTCorr<25)  FR=0.244504;
      else if(PTCorr<35)  FR=0.222135;
      else                FR=0.219002;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1NoChi_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.104934 ;
      else if(PTCorr<20)  FR=0.0810565;
      else if(PTCorr<25)  FR=0.085597 ;
      else if(PTCorr<35)  FR=0.0684988;
      else                FR=0.0665389;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.0977776;
      else if(PTCorr<20)  FR=0.0729936;
      else if(PTCorr<25)  FR=0.10299  ;
      else if(PTCorr<35)  FR=0.116561 ;
      else                FR=0.0859202;
    }
    else{
      if     (PTCorr<15)  FR=0.0610895;
      else if(PTCorr<20)  FR=0.0923621;
      else if(PTCorr<25)  FR=0.0671893;
      else if(PTCorr<35)  FR=0.128802 ;
      else                FR=0.0967269;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1NoChi_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.308976;
      else if(PTCorr<20)  FR=0.229748;
      else if(PTCorr<25)  FR=0.202024;
      else if(PTCorr<35)  FR=0.211457;
      else                FR=0.215925;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.351577;
      else if(PTCorr<20)  FR=0.275117;
      else if(PTCorr<25)  FR=0.266206;
      else if(PTCorr<35)  FR=0.207133;
      else                FR=0.252437;
    }
    else{
      if     (PTCorr<15)  FR=0.364771;
      else if(PTCorr<20)  FR=0.309155;
      else if(PTCorr<25)  FR=0.302773;
      else if(PTCorr<35)  FR=0.265689;
      else                FR=0.275625;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.30168 ;
      else if(PTCorr<20)  FR=0.205724;
      else if(PTCorr<25)  FR=0.167485;
      else if(PTCorr<35)  FR=0.163753;
      else                FR=0.16684 ;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.336895;
      else if(PTCorr<20)  FR=0.234467;
      else if(PTCorr<25)  FR=0.219452;
      else if(PTCorr<35)  FR=0.176959;
      else                FR=0.200026;
    }
    else{
      if     (PTCorr<15)  FR=0.34518 ;
      else if(PTCorr<20)  FR=0.267018;
      else if(PTCorr<25)  FR=0.247139;
      else if(PTCorr<35)  FR=0.228795;
      else                FR=0.225253;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1NoChi_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.0730267;
      else if(PTCorr<20)  FR=0.0422678;
      else if(PTCorr<25)  FR=0.0394787;
      else if(PTCorr<35)  FR=0.0292311;
      else                FR=0.0247411;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.075509 ;
      else if(PTCorr<20)  FR=0.0421607;
      else if(PTCorr<25)  FR=0.0535583;
      else if(PTCorr<35)  FR=0.0579744;
      else                FR=0.0374653;
    }
    else{
      if     (PTCorr<15)  FR=0.0488572;
      else if(PTCorr<20)  FR=0.0598638;
      else if(PTCorr<25)  FR=0.036773 ;
      else if(PTCorr<35)  FR=0.0710471;
      else                FR=0.0479654;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1NoChi_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.243067 ;
      else if(PTCorr<20)  FR=0.120063 ;
      else if(PTCorr<25)  FR=0.100958 ;
      else if(PTCorr<35)  FR=0.0993144;
      else                FR=0.0935684;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.288401;
      else if(PTCorr<20)  FR=0.159155;
      else if(PTCorr<25)  FR=0.150465;
      else if(PTCorr<35)  FR=0.106834;
      else                FR=0.122566;
    }
    else{
      if     (PTCorr<15)  FR=0.308748;
      else if(PTCorr<20)  FR=0.195025;
      else if(PTCorr<25)  FR=0.182063;
      else if(PTCorr<35)  FR=0.153384;
      else                FR=0.148075;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.23621  ;
      else if(PTCorr<20)  FR=0.1076   ;
      else if(PTCorr<25)  FR=0.0816986;
      else if(PTCorr<35)  FR=0.0747883;
      else                FR=0.068686 ;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.275343 ;
      else if(PTCorr<20)  FR=0.135587 ;
      else if(PTCorr<25)  FR=0.121004 ;
      else if(PTCorr<35)  FR=0.0901332;
      else                FR=0.0936505;
    }
    else{
      if     (PTCorr<15)  FR=0.291057;
      else if(PTCorr<20)  FR=0.169294;
      else if(PTCorr<25)  FR=0.145199;
      else if(PTCorr<35)  FR=0.130402;
      else                FR=0.118149;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.220591 ;
      else if(PTCorr<20)  FR=0.0957646;
      else if(PTCorr<25)  FR=0.0756069;
      else if(PTCorr<35)  FR=0.0664369;
      else                FR=0.0596936;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.251618 ;
      else if(PTCorr<20)  FR=0.115786 ;
      else if(PTCorr<25)  FR=0.102281 ;
      else if(PTCorr<35)  FR=0.0847535;
      else                FR=0.0776838;
    }
    else{
      if     (PTCorr<15)  FR=0.278727;
      else if(PTCorr<20)  FR=0.147944;
      else if(PTCorr<25)  FR=0.125158;
      else if(PTCorr<35)  FR=0.116875;
      else                FR=0.10205 ;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1sig4NoChi_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.333358;
      else if(PTCorr<20)  FR=0.236995;
      else if(PTCorr<25)  FR=0.23349 ;
      else if(PTCorr<35)  FR=0.17589 ;
      else                FR=0.168586;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.323098;
      else if(PTCorr<20)  FR=0.207417;
      else if(PTCorr<25)  FR=0.261432;
      else if(PTCorr<35)  FR=0.288036;
      else                FR=0.21438 ;
    }
    else{
      if     (PTCorr<15)  FR=0.191651;
      else if(PTCorr<20)  FR=0.237084;
      else if(PTCorr<25)  FR=0.170018;
      else if(PTCorr<35)  FR=0.283506;
      else                FR=0.221901;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1sig4NoChi_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.434692;
      else if(PTCorr<20)  FR=0.303991;
      else if(PTCorr<25)  FR=0.246176;
      else if(PTCorr<35)  FR=0.265866;
      else                FR=0.267146;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.472631;
      else if(PTCorr<20)  FR=0.360498;
      else if(PTCorr<25)  FR=0.33543 ;
      else if(PTCorr<35)  FR=0.252456;
      else                FR=0.315846;
    }
    else{
      if     (PTCorr<15)  FR=0.463157;
      else if(PTCorr<20)  FR=0.373889;
      else if(PTCorr<25)  FR=0.363318;
      else if(PTCorr<35)  FR=0.321106;
      else                FR=0.339224;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop4IPp2p1sig4NoChi_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.432992;
      else if(PTCorr<20)  FR=0.2985  ;
      else if(PTCorr<25)  FR=0.244148;
      else if(PTCorr<35)  FR=0.248196;
      else                FR=0.247938;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.468875;
      else if(PTCorr<20)  FR=0.344572;
      else if(PTCorr<25)  FR=0.323177;
      else if(PTCorr<35)  FR=0.259591;
      else                FR=0.29687 ;
    }
    else{
      if     (PTCorr<15)  FR=0.455616;
      else if(PTCorr<20)  FR=0.359923;
      else if(PTCorr<25)  FR=0.338585;
      else if(PTCorr<35)  FR=0.314778;
      else                FR=0.318862;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.235023 ;
      else if(PTCorr<20)  FR=0.126737 ;
      else if(PTCorr<25)  FR=0.10877  ;
      else if(PTCorr<35)  FR=0.0810774;
      else                FR=0.0657765;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.23597  ;
      else if(PTCorr<20)  FR=0.121983 ;
      else if(PTCorr<25)  FR=0.135113 ;
      else if(PTCorr<35)  FR=0.146711 ;
      else                FR=0.0987709;
    }
    else{
      if     (PTCorr<15)  FR=0.150243 ;
      else if(PTCorr<20)  FR=0.151052 ;
      else if(PTCorr<25)  FR=0.0953203;
      else if(PTCorr<35)  FR=0.159848 ;
      else                FR=0.114909 ;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.345149;
      else if(PTCorr<20)  FR=0.164104;
      else if(PTCorr<25)  FR=0.124837;
      else if(PTCorr<35)  FR=0.126013;
      else                FR=0.118599;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.389373;
      else if(PTCorr<20)  FR=0.208879;
      else if(PTCorr<25)  FR=0.189129;
      else if(PTCorr<35)  FR=0.131404;
      else                FR=0.155852;
    }
    else{
      if     (PTCorr<15)  FR=0.39351 ;
      else if(PTCorr<20)  FR=0.237086;
      else if(PTCorr<25)  FR=0.217314;
      else if(PTCorr<35)  FR=0.183934;
      else                FR=0.185743;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.342129;
      else if(PTCorr<20)  FR=0.161717;
      else if(PTCorr<25)  FR=0.124337;
      else if(PTCorr<35)  FR=0.124102;
      else                FR=0.112075;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.3853  ;
      else if(PTCorr<20)  FR=0.200014;
      else if(PTCorr<25)  FR=0.181963;
      else if(PTCorr<35)  FR=0.143144;
      else                FR=0.149275;
    }
    else{
      if     (PTCorr<15)  FR=0.386585;
      else if(PTCorr<20)  FR=0.228739;
      else if(PTCorr<25)  FR=0.203977;
      else if(PTCorr<35)  FR=0.190597;
      else                FR=0.180573;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp2p1sig4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.323899 ;
      else if(PTCorr<20)  FR=0.145267 ;
      else if(PTCorr<25)  FR=0.114264 ;
      else if(PTCorr<35)  FR=0.102888 ;
      else                FR=0.0929075;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.355195;
      else if(PTCorr<20)  FR=0.170462;
      else if(PTCorr<25)  FR=0.15116 ;
      else if(PTCorr<35)  FR=0.127486;
      else                FR=0.118954;
    }
    else{
      if     (PTCorr<15)  FR=0.370076;
      else if(PTCorr<20)  FR=0.199035;
      else if(PTCorr<25)  FR=0.170545;
      else if(PTCorr<35)  FR=0.161002;
      else                FR=0.146339;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_NearB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.276577 ;
      else if(PTCorr<20)  FR=0.142112 ;
      else if(PTCorr<25)  FR=0.117096 ;
      else if(PTCorr<35)  FR=0.0858086;
      else                FR=0.0684839;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.351973;
      else if(PTCorr<20)  FR=0.154119;
      else if(PTCorr<25)  FR=0.154107;
      else if(PTCorr<35)  FR=0.162744;
      else                FR=0.109306;
    }
    else{
      if     (PTCorr<15)  FR=0.325502;
      else if(PTCorr<20)  FR=0.252336;
      else if(PTCorr<25)  FR=0.144032;
      else if(PTCorr<35)  FR=0.212277;
      else                FR=0.151429;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_AwayB_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.363658;
      else if(PTCorr<20)  FR=0.169303;
      else if(PTCorr<25)  FR=0.127615;
      else if(PTCorr<35)  FR=0.128469;
      else                FR=0.121185;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.423617;
      else if(PTCorr<20)  FR=0.219131;
      else if(PTCorr<25)  FR=0.196251;
      else if(PTCorr<35)  FR=0.134896;
      else                FR=0.160779;
    }
    else{
      if     (PTCorr<15)  FR=0.475586;
      else if(PTCorr<20)  FR=0.269348;
      else if(PTCorr<25)  FR=0.242709;
      else if(PTCorr<35)  FR=0.19999 ;
      else                FR=0.201519;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_TrkIsoVVLConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.362155;
      else if(PTCorr<20)  FR=0.167222;
      else if(PTCorr<25)  FR=0.125889;
      else if(PTCorr<35)  FR=0.120206;
      else                FR=0.110192;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.422128;
      else if(PTCorr<20)  FR=0.213486;
      else if(PTCorr<25)  FR=0.189325;
      else if(PTCorr<35)  FR=0.140276;
      else                FR=0.151332;
    }
    else{
      if     (PTCorr<15)  FR=0.473036;
      else if(PTCorr<20)  FR=0.268132;
      else if(PTCorr<25)  FR=0.232456;
      else if(PTCorr<35)  FR=0.201764;
      else                FR=0.193886;
    }
  }
  else if(Option=="POGTIsop20IPp01p05sig4Chi4_POGTIsop6IPp01p05sig4Chi4_NoFilterConeSUSY"){
    if(fEta<0.9){
      if     (PTCorr<15)  FR=0.341034 ;
      else if(PTCorr<20)  FR=0.151259 ;
      else if(PTCorr<25)  FR=0.118491 ;
      else if(PTCorr<35)  FR=0.106629 ;
      else                FR=0.0954289;
    }
    else if(fEta<1.6){
      if     (PTCorr<15)  FR=0.38962 ;
      else if(PTCorr<20)  FR=0.183194;
      else if(PTCorr<25)  FR=0.160408;
      else if(PTCorr<35)  FR=0.133272;
      else                FR=0.124299;
    }
    else{
      if     (PTCorr<15)  FR=0.455303;
      else if(PTCorr<20)  FR=0.232687;
      else if(PTCorr<25)  FR=0.197538;
      else if(PTCorr<35)  FR=0.182767;
      else                FR=0.163027;
    }
  }
  else {cout<<"No Such FR!"<<endl;}



  return FR;
}


float  Mar2018_TrilepCR::GetFakeWeight(std::vector<snu::KMuon>& MuLColl, TString MuLID, TString MuTID, vector<snu::KJet>& BJetNoVetoColl, TString Option){

  float fakeweight=-1.; int NLooseNotTight=0;

  TString FilterInfo="", ConeMethod="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";


  for(int i=0; i<(int) MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID, "Roch")){
      float FR=0.;
      bool IsNearB = IsNearBJet(MuLColl.at(i), BJetNoVetoColl);
      TString FlavOpt = IsNearB? "_NearB":"_AwayB";

      FR=FakeRateData(MuLColl.at(i),MuTID+"_"+MuLID+FlavOpt+"_"+FilterInfo+ConeMethod);
      //cout<<i<<FlavOpt<<" MuFR"<<FR<<endl;
      fakeweight*=-FR/(1.-FR);
      NLooseNotTight++;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;


}

float Mar2018_TrilepCR::GetFakeWeight(std::vector<snu::KMuon>& MuLColl, TString MuLID, TString MuTID, TString Option){

  TString FilterInfo="", ConeMethod="";
  if     (Option.Contains("NoFilter"))  FilterInfo="NoFilter";
  else if(Option.Contains("TrkIsoVVL")) FilterInfo="TrkIsoVVL";
  if     (Option.Contains("ConeE"))     ConeMethod="ConeE";
  else if(Option.Contains("ConeSUSY"))  ConeMethod="ConeSUSY";
  else if(Option.Contains("FOPt"))      ConeMethod="FOPt";

  float fakeweight=-1.; int NLooseNotTight=0;
  for(int i=0; i<MuLColl.size(); i++){
    if(!PassIDCriteria(MuLColl.at(i), MuTID,"Roch")){
      float FR=0.;
      FR=FakeRateData(MuLColl.at(i),MuTID+"_"+MuLID+"_"+FilterInfo+ConeMethod);
      fakeweight*=-FR/(1.-FR);
      NLooseNotTight++;
    }
  }
  if(NLooseNotTight==0) return 0.;

  return fakeweight;
}



int Mar2018_TrilepCR::GetFakeLepSrcIdx(snu::KMuon Mu, std::vector<snu::KTruth> TruthColl){

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



int Mar2018_TrilepCR::GetFakeLepJetSrcType(snu::KMuon Mu, std::vector<snu::KJet> JetColl){
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


float Mar2018_TrilepCR::ConeCorrectedPT(snu::KElectron Ele, float TightIsoCut){

  //float PTCorr=Ele.Pt()*(1+Ele.PFRelIso(0.3));
  float PTCorr=Ele.Pt()*(1+min((float) 0,(float) Ele.PFRelIso(0.3)-TightIsoCut));
  return PTCorr;

}


float Mar2018_TrilepCR::ConeCorrectedPT(snu::KMuon Mu, float TightIsoCut){

  //float PTCorr=Mu.Pt();
  //float PTCorr=Mu.Pt()*(1.+RochIso(Mu,"0.4"));
  float PTCorr=Mu.Pt()*(1.+max((float) 0.,RochIso(Mu,"0.4")-TightIsoCut));
  return PTCorr;

}

bool Mar2018_TrilepCR::IsNearBJet(snu::KMuon Mu, std::vector<snu::KJet> bjetNoVetoColl){

  bool IsNearB=false;
  for(int i=0; i<bjetNoVetoColl.size(); i++){
    if(Mu.DeltaR(bjetNoVetoColl.at(i))<0.4){
      IsNearB=true; break;
    }
  }

  return IsNearB;
}




int Mar2018_TrilepCR::StepPassed(std::vector<snu::KMuon> MuColl, std::vector<snu::KElectron> EleColl, std::vector<snu::KJet> JetColl, std::vector<snu::KJet> BJetColl, float MET, TString Option){

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


int Mar2018_TrilepCR::NPromptFake_Mu(std::vector<snu::KMuon> MuColl, std::vector<snu::KTruth> truthColl, TString Option){

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

int Mar2018_TrilepCR::NPromptFake_Ele(std::vector<snu::KElectron> EleColl, std::vector<snu::KTruth> truthColl, TString Option){

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



bool Mar2018_TrilepCR::IsConvCand(snu::KElectron Ele, std::vector<snu::KTruth> TruthColl, TString Option){

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



void Mar2018_TrilepCR::FillCutFlow(TString cut, float weight){
  
  if(GetHist("cutflow_W") && GetHist("cutflow_N")){
    GetHist("cutflow_W")->Fill(cut,weight);
    GetHist("cutflow_N")->Fill(cut,1);
  }
  else{
    if(!GetHist("cutflow_W")){
      AnalyzerCore::MakeHistograms("cutflow_W", 11, 0., 11.);
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(1,"NoCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(2,"Trigger");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(3,"METFilter");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_W")->GetXaxis()->SetBinLabel(5,"PreSel");
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
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(2,"Trigger");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(3,"METFilter");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(4,"VertexCut");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(5,"PreSel");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(6,"M(#mu#mu)>4");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(7,"#geq1b");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(8,"2j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(9,"3j");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(10,"ZVeto");
      GetHist("cutflow_N")->GetXaxis()->SetBinLabel(11,"Nlj2");
    }
  }
}



void Mar2018_TrilepCR::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}


void Mar2018_TrilepCR::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();


  Message("Made histograms", INFO);
  // **
  // *  Remove//Overide this Mar2018_TrilepCRCore::MakeHistograms() to make new hists for your analysis
  // **
  
}



void Mar2018_TrilepCR::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}
